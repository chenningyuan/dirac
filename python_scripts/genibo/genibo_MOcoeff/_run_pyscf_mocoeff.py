"""Functions to prepare pyscf input and run calculations."""
from __future__ import absolute_import

import os
import re
import copy
import subprocess
import warnings
from pyscf import gto
from pyscf import scf
from pyscf import lib
from pyscf.scf.chkfile import dump_scf
import numpy as np
from ._utils import *
from ._AOMObases import *
import sys
np.set_printoptions(threshold=sys.maxsize)

class InputError(Exception):
    pass

def gaussian_norm(exponent,l):
    """
      Taking over from the original Fortran pyscf interface of B. Senjean
      as I could not find the right pyscf function to provide this.

      ! DIRAC norm:
      norm_cons = (2.D0*expo/pi)**(0.75D0)*dsqrt(2.D0**l)*dsqrt(2.D0*expo)**l
      ! Addition for PYSCF.
      ! I don't know why, but one should not apply this additional factor to s and p shells.
      if (l .ge. 2) then
        norm_cons = norm_cons*2.D0*dsqrt(pi/(double_factorial(2*l+1)))
      end if
    """

    # DIRAC norm
    norm_cons =  np.power(2.*exponent/np.pi,0.75)*np.power(4.*exponent,l/2)
    if l>1:
        # constant factor for the higher l-values
        norm_cons *= np.sqrt(4*np.pi) 
        # evaluate the double factorial (2l+1)!! factor
        for ll in range(1,l+1):
            norm_cons /=  np.sqrt(2*ll+1)
    return norm_cons

def run_pyscf_mocoeff(molecule,
                      basis=None,
                      relativistic=False,
                      restricted=True,
                      spherical=False,
                      symmetry=False,
                      uncontract=True,
                      restart_wfn=False,
                      openshell=False,
                      save=False,
                      just_write=False):
                      
    """This function runs a PySCF calculation.

    Args:
        molecule: An instance of the MolData class.

    Returns:
        molecule: The MolecularData object.

    Raises:
        pyscf errors: An error from pyscf.
    """

    # Create PySCF geometry string.
    molecule_string = create_geometry_string(molecule.geometry)
    xyz_file = molecule.filename + '.xyz'
    with open(xyz_file, 'w') as f:
     f.write(str(molecule.n_atoms)+'\n')
     f.write(molecule.filename + ' # anything can be in this line\n')
     f.write(molecule_string)

    if uncontract: basis = "unc-"+basis
    
    if symmetry: sym = 1
    if not symmetry: sym = 0

    if spherical:
      pureamd = 0
      pureamf = 0
      print("**WARNING** ROSE does not work yet with Spherical instead of Cartesian functions. (though I didn't see any difference using the PYSCF interface for now.)")
    else:
      pureamd = 1
      pureamf = 1

    if not spherical: cart = True
    if spherical: cart = False

    mol = gto.M(
          atom = molecule_string, # it can also be molecule.geometry ! PySCF handles both format.
          basis = basis,
          symmetry = sym,
          spin = molecule.get_spin(),
          charge = molecule.charge,
          verbose = 5,
          cart = cart,
          output = molecule.filename + ".out")

    if not relativistic:
       if restricted and not openshell: mf = scf.RHF(mol)
       if restricted and openshell: mf = scf.ROHF(mol)
       if not restricted: mf = scf.UHF(mol)
    else:
       if restricted and not openshell: mf = scf.RHF(mol).x2c()
       if restricted and openshell: mf = scf.ROHF(mol).x2c()
       if not restricted: mf = scf.UHF(mol).x2c()
    
    if save:
      mf.chkfile = molecule.filename + ".chk"

    #When HOMO and LUMO orbitals are degenerated or quasi-degenerated, fill
    #fractional number of electrons in the degenerated orbitals.
    #This is the equivalent of the average of configuration in DIRAC.
    mf = scf.addons.frac_occ(mf)

    # just write a new .pyscf file without running any calculation
    if just_write is False:
      if restart_wfn is not False:
        dm = mf.init_guess_by_chkfile(restart_wfn + ".chk")
        mf.max_cycle = 0
        mf.kernel(dm)
      else:
        mf.kernel()
    else:
      mf = just_write['wfn_object']
    
    # Write hdf5 interface file
    molecule.write_hdf5(file_name=molecule.filename+'.h5')
    nshells = mol.nbas
    basis_funcs = []
    for i in range(nshells):
        coord=[mol.bas_coord(i)[0],
               mol.bas_coord(i)[1],
               mol.bas_coord(i)[2]]
        orb_momentum = mol.bas_angular(i)+1
        exponent     = mol.bas_exp(i)
        contractions = mol.bas_ctr_coeff(i) 
        # contractions is a 2d array as one may use general contraction
        for j in range(contractions.shape[1]): 
            coefficient = contractions[:,j]

            # coefficients need to be multiplied with normalization constant to be compatible
            for k in range(contractions.shape[0]): 
                coefficient[k] *= gaussian_norm(exponent[k],orb_momentum-1)

            basis_func = BasisFunc(orb_momentum=orb_momentum,
                       exponent=exponent,
                       coefficient=coefficient,
                       coord=coord)
            basis_funcs.append(basis_func)

    basis_set = AOBasis(1,basis_funcs)
    eigenvalues = mf.mo_energy.copy()
    coeffs      = mf.mo_coeff.copy().T
    mo_basis = MOBasis(aobasis=basis_set,
                       eigenvalues=eigenvalues,
                       algebra='r',
                       restricted=True,
                       coeffs=coeffs)
    mo_basis.write_hdf5(file_name=molecule.filename+'.h5')
    # h5 file completed

    scf_e = mf.e_tot
    print("Hartree-Fock energy: {}".format(scf_e))

    return molecule, mol, mf

#--------------------------------------------------------------------------
#-------------------------- RETRIEVE  WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

def fetch_MOcoefficients(basis):
   """
   Fetch the MO coefficients from the hdf5 basis file.
   """
   # fetch ibo energies and coefficients from the hdf5 file
   ibo = MOBasisFromFile(basis_file=basis+'.h5')
   nao = ibo.MOs[0].aobasis.nao
   nmo = ibo.nmo
   mostring = ibo.make_mostring(mo_range=(1,nmo))
   if ibo.nbar == 1: # restricted
       eigenvalues = ibo.get_eigenvalues(mostring)
       mocoeffs    = ibo.get_mocoeffs(mostring).reshape((nmo,nao)).T
   else: # unrestricted
       nmo = nmo/2
       eigenvalues = [ibo.get_eigenvalues(mostring)[::2]]
       eigenvalues.append(ibo.get_eigenvalues(mostring)[1::2])
       mocoeffs    = [ibo.get_mocoeffs(mostring)[::2,:].reshape((nmo,nao)).T]
       mocoeffs.append(ibo.get_mocoeffs(mostring)[1::2,:].reshape((nmo,nao)).T)
   return nmo, eigenvalues, mocoeffs

#--------------------------------------------------------------------------
#-------------------------- STORE IBO WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

def replace_pyscf_MOcoefficients(scf_wfn,mol,basis):
   """
   Replace the MO coefficients on the pyscf object with those
   on the hdf5 basis_file and write these to file.
   Returns number of mo's that were replaced.
   """
   ibo_wfn = copy.copy(scf_wfn)
   # fetch ibo energies and coefficients from the hdf5 file
   ibo = MOBasisFromFile(basis_file=basis+'.h5')
   nao = ibo.MOs[0].aobasis.nao
   nmo = ibo.nmo
   mostring = ibo.make_mostring(mo_range=(1,nmo))
   # replace the scf orbitals, keeping the rest of the wfn object the same
   if ibo.nbar == 1: # restricted
       ibo_wfn.mo_energy[:nmo]  = ibo.get_eigenvalues(mostring)
       ibo_wfn.mo_coeff[:,:nmo] = ibo.get_mocoeffs(mostring).reshape((nmo,nao)).T
   else: # unrestricted
       nmo = nmo/2
       ibo_wfn.mo_energy[0][:nmo]  = ibo.get_eigenvalues(mostring)[::2]
       ibo_wfn.mo_energy[1][:nmo]  = ibo.get_eigenvalues(mostring)[::2]
       ibo_wfn.mo_coeff[0][:,:nmo] = ibo.get_mocoeffs(mostring)[::2,:].reshape((nmo,nao)).T
       ibo_wfn.mo_coeff[1][:,:nmo] = ibo.get_mocoeffs(mostring)[1::2,:].reshape((nmo,nao)).T
   # store in pyscf format
   e_tot = 0.0 # just set e_tot to something so this can be written out
   dump_scf(mol, "ribo.chk", e_tot, ibo_wfn.mo_energy, ibo_wfn.mo_coeff,  ibo_wfn.mo_occ)

   return nmo, ibo_wfn
