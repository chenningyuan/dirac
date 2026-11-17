"""Functions to prepare psi4 input and run calculations."""
from __future__ import absolute_import

import os
import re
import subprocess
import warnings
import psi4
from psi4.driver.procrouting import proc_util
import numpy as np
from ._utils import *
from ._AOMObases import *
import sys
np.set_printoptions(threshold=sys.maxsize)

class InputError(Exception):
    pass


def atom_d2h_occupations(Z, charge, multiplicity, p_order=("x", "y", "z")):
    """Explicit D2h occupations (DOCC, SOCC) for a free open-shell atom of the first three rows (s and p shells).

    The three p orbitals of an open-shell p atom are degenerate; psi4 then picks which one is doubly occupied
    (ROHF) in a way that depends on the psi4 version, and the IAO/IBO charges depend on that choice. Fixing the
    occupations in D2h makes the reference deterministic: p electrons beyond a half-filled shell doubly occupy the
    axes in the order p_order, the others are singly occupied in that order.
    Returns None if the atom is not treated (then the C1 calculation is done as before).
    """
    nel = Z - charge
    if not 1 <= Z <= 18 or nel < 1: return None
    # electrons per shell, aufbau 1s 2s 2p 3s 3p
    cap = [("s", 2), ("s", 2), ("p", 6), ("s", 2), ("p", 6)]
    shells = []
    left = nel
    for l, c in cap:
        n = min(left, c); left -= n; shells.append((l, n))
    if left > 0: return None
    docc = [0]*8; socc = [0]*8
    irrep = {"z": 5, "y": 6, "x": 7}                    # D2h: Ag B1g B2g B3g Au B1u B2u B3u
    for l, n in shells:
        if l == "s":
            if n == 2: docc[0] += 1
            elif n == 1: socc[0] += 1
        else:
            na = min(n, 3); nb = n - na
            axes = list(p_order)
            for a in axes[:nb]: docc[irrep[a]] += 1
            for a in axes[nb:na]: socc[irrep[a]] += 1
    if sum(socc) != multiplicity - 1 or 2*sum(docc) + sum(socc) != nel: return None
    return docc, socc

def run_psi4_mocoeff(molecule,
                     basis=None,
                     spherical=False,
                     relativistic=False,
                     restricted=True,
                     openshell=False,
                     fcidump=False,
                     save=False,
                     restart_wfn=False,
                     MAXITER=100,
                     e_convergence=1e-6,
                     d_convergence=1e-6,
                     uncontract=True,
                     just_write=False):
    """This function runs a Psi4 calculation.

    Args:
        molecule: An instance of the MolData class.
        spherical: use spherical coordinates for the Gaussian type orbitals. If False, it uses cartesian coordinates.
        restricted: use of restricted algorithm ; i.e. the coefficient matrix will span half the space compare to the unrestricted calculation.
        openshell: perform an openshell calculation (rohf or uhf depending on the value of the boolean "restricted").
        fcidump: extract the fcidump file.
        restart_wfn: ask to restart the calculation with a given wfn.
        MAXITER: maximum number of iterations
        e_convergence: Convergence criterion for SCF energy
        d_convergence: Convergence criterion for SCF density, which is defined as the RMS value of the orbital gradient
        uncontract: use uncontracted basis set.

    Returns:
        molecule: The MolecularData object.

    Raises:
        psi4 errors: An error from psi4.
    """

    if (basis is None): raise ValueError("basis must be specified.")

    psi4.core.clean()

    # Create Psi4 geometry string.
    molecule_string = create_geometry_string(molecule.geometry)
    xyz_file = molecule.filename + '.xyz'
    with open(xyz_file, 'w') as f:
     f.write(str(molecule.n_atoms)+'\n')
     f.write(molecule.filename + ' # anything can be in this line\n')
     f.write(molecule_string)

    d2h = None
    if restricted and openshell and not relativistic and molecule.n_atoms == 1 and molecule.protons[0] > 0:
        d2h = atom_d2h_occupations(molecule.protons[0], molecule.charge, molecule.get_multiplicity())
    add_input = "\nsymmetry d2h\nno_reorient\nno_com" if d2h else "\nsymmetry c1\nno_reorient\nno_com"
    for _opt in ("DOCC", "SOCC"):                              # do not inherit occupations of an earlier call
        psi4.core.revoke_global_option_changed(_opt)
        psi4.core.revoke_local_option_changed("SCF", _opt)
    if uncontract: basis = basis + "-decon"
    if spherical: 
       puream = "true"
       pureamd = 0
       pureamf = 0
    else:
       puream = "false"
       pureamd = 1
       pureamf = 1
    if restricted and not openshell: scf = "rhf"
    elif restricted and openshell: scf = "rohf"
    elif not restricted: scf = "uhf"
    if not restart_wfn:
       guess = "sad"
    elif restart_wfn is not False:
       guess = "read"
       scf_wfn = restart_wfn

    psi4.core.set_output_file(molecule.filename + ".out", False)

    if relativistic:
       psi4.set_options({'basis':basis,
                      'basis_relativistic': basis,
                      'reference': scf,
                      'relativistic': "X2C",
                      'e_convergence': e_convergence,
                      'd_convergence': d_convergence,
                      'scf_type': 'direct',
                      'guess': guess,
                      'MAXITER': MAXITER,
                      'puream':puream})
    else:
       psi4.set_options({'basis': basis,
                      'reference': scf,
                      'e_convergence': e_convergence,
                      'd_convergence': d_convergence,
                      'scf_type': 'direct',
                      'guess': guess,
                      'MAXITER': MAXITER,
                      'puream':puream})
    
    if d2h: psi4.set_options({'docc': d2h[0], 'socc': d2h[1]})
    if (molecule.multiplicity == None) and (molecule.charge == 0):
        mol = psi4.geometry(molecule_string + add_input)
    else:
        mol = psi4.geometry(str(molecule.charge) + " " + str(molecule.get_multiplicity()) + "\n" + molecule_string + add_input)
    
    if just_write is False:
       # Get the SCF wavefunction & energies
       scf_e, scf_wfn = psi4.energy('scf', return_wfn=True)
       print("Hartree-Fock energy: {}".format(scf_e))
       # Write the wavefunction into file.
       if not restart_wfn:
         # Write a fchk file. This is nice to check if the .psi4 extension is accurate or not...
         fchk_writer = psi4.core.FCHKWriter(scf_wfn)
         fchk_writer.write(molecule.filename + '.fchk')
       if save:
         scf_wfn.to_file(molecule.filename)
    else:
       scf_wfn = just_write["wfn_object"]
       scf_e = scf_wfn.energy()

    # Write hdf5 interface file
    molecule.write_hdf5(file_name=molecule.filename+'.h5')
    basis = scf_wfn.basisset()
    nshells = basis.nshell()
    basis_funcs = []
    for i in range(nshells):
        orb_momentum = basis.shell(i).am+1
        exponent     = [basis.shell(i).exp(j)  for j in range(basis.shell(i).nprimitive)]
        coefficient  = [basis.shell(i).coef(j) for j in range(basis.shell(i).nprimitive)]
        coord=[mol.fx(basis.shell(i).ncenter),
               mol.fy(basis.shell(i).ncenter),
               mol.fz(basis.shell(i).ncenter)]
        basis_func = BasisFunc(orb_momentum=orb_momentum,
                       exponent=exponent,
                       coefficient=coefficient,
                       coord=coord)
        basis_funcs.append(basis_func)

    if not spherical:
        basis_set = AOBasis(1,basis_funcs)
    else:
        basis_set = AOBasis(2,basis_funcs)
    if restricted: # can be passed on as is
       if d2h:   # symmetry-blocked orbitals: back to the AO basis, ordered by energy
          eigenvalues=scf_wfn.epsilon_a_subset("AO","ALL").np
          coeffs=scf_wfn.Ca_subset("AO","ALL").np.T
       else:
          eigenvalues=scf_wfn.epsilon_a().np
          coeffs=scf_wfn.Ca().np.T
    else: # should be concatenated before calling MOBasis
       eigenvalues=np.concatenate((scf_wfn.epsilon_a().np,scf_wfn.epsilon_b().np))
       coeffs=np.concatenate((scf_wfn.Ca().np.T,scf_wfn.Cb().np.T))
        
    mo_basis = MOBasis(aobasis=basis_set,
                       eigenvalues=eigenvalues,
                       algebra='r',
                       restricted=restricted,
                       coeffs=coeffs)
    mo_basis.write_hdf5(file_name=molecule.filename+'.h5')
    # h5 file completed

#   Extract the FCIDUMP if asked.
    if fcidump: psi4.fcidump(scf_wfn)

    return molecule, scf_wfn


#--------------------------------------------------------------------------
#-------------------------- STORE IBO WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

def replace_psi4_MOcoefficients(scf_wfn,basis,return_wfn=False):
   """
   Replace the MO coefficients in a psi4 wavefunction archive with those
   on the hdf5 basis_file and reload the reduced wavefunction.
   Returns number of mo's that were replaced. If return_wfn is True, also
   returns the reduced wavefunction.
   """
   # fetch ibo energies and coefficients from the hdf5 file
   ibo = MOBasisFromFile(basis_file=basis+'.h5')
   nao = ibo.MOs[0].aobasis.nao
   nmo = ibo.nmo
   mostring = ibo.make_mostring(mo_range=(1,nmo))
   scf_wfn.to_file(basis)
   wfn_data = np.load(basis+'.npy', allow_pickle=True).item()
   if ibo.nbar == 1: # restricted
       eigenvalues = ibo.get_eigenvalues(mostring)
       mocoeffs    = ibo.get_mocoeffs(mostring).reshape((nmo,nao)).T
       wfn_data["vector"]["epsilon_a"] = eigenvalues.copy()
       wfn_data["vector"]["epsilon_b"] = eigenvalues.copy()
       wfn_data["matrix"]["Ca"] = mocoeffs.copy()
       wfn_data["matrix"]["Cb"] = mocoeffs.copy()
   else: # unrestricted
       nmo = nmo//2
       wfn_data["vector"]["epsilon_a"] = ibo.get_eigenvalues(mostring)[::2].copy()
       wfn_data["vector"]["epsilon_b"] = ibo.get_eigenvalues(mostring)[1::2].copy()
       wfn_data["matrix"]["Ca"] = ibo.get_mocoeffs(mostring)[::2,:].reshape((nmo,nao)).T.copy()
       wfn_data["matrix"]["Cb"] = ibo.get_mocoeffs(mostring)[1::2,:].reshape((nmo,nao)).T.copy()
   wfn_data["dimension"]["nmopi"] = (nmo,)
   wfn_data["dimension"]["frzvpi"] = (0,)
   wfn_data["int"]["nmo"] = nmo
   np.save(basis+'.npy', wfn_data)
   reduced_wfn = psi4.core.Wavefunction.from_file(basis)

   if return_wfn:
      return nmo, reduced_wfn
   return nmo
