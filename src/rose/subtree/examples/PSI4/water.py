from genibo_MOcoeff import *
import os
import sys
import subprocess
import numpy as np
import psi4

# LV: WARNING !!!
# LV: this example needs fixing, the psi4 atomic run produces non-symmetric orbitals as fractional occupation is not yet used.
# LV: this asymmetry makes the result depend on the version of psi4 and also yields artefacts in the produced IAOs and IBOs.

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory             = os.getenv('ROSEDIR')
rose_executable            = rose_directory + "/bin/rose.x"
plams_executable           = "plams"
data_directory             = os.getcwd()
xyz_directory              = rose_directory + "/examples/xyz_files/"
version                    = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent                   = [2,3,4][2]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted                 = True      # Change from unrestricted to restricted formalism.
relativistic               = False     # scalar X2C
basis1                     = 'cc-pVDZ' # Basis set for the full problem.
basis2                     = 'cc-pVDZ' # Basis set for the fragments
charge                     = 0         # Charge of the full problem
multiplicity               = 1         # Multiplicity of the full problem
spherical                  = False     # Spherical or Cartesian coordinate GTOs. If True, it will do the calculation using spherical and store information using cartesians.
uncontract                 = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
openshell                  = True      # Openshell or closed-shell system
test                       = True      # Perform tests which will be printed in the terminal.
restart                    = True      # Restart calculation with IAO and IBO.
description                = ""        # Add an optional description to the default name of the files
fcidump                    = False     # Extract FCIDUMP file from PSI4
save                       = True      # Save psi4 wavefunction into .npy file.
include_core               = False     # Include core in the localization procedure, default : False
additional_virtuals_cutoff = 2.0       # Defaut is 2.0 hartrees
frag_threshold             = 2.0       # Energy Cutoff to determine number of reference fragment hard virtuals to use
run_postSCF                = True      # Run post-SCF calculation (CASCI/CASSCF) with original, IBO and AVAS orbitals.
avas_frag                  = []        # Fragment IAO file to be extracted from ROSE (used for AVAS for instance). Keep empty otherwise.
mo_avas_frag               = []        # list of MOs (starting at 1) to consider in AVAS. [[mo_frag1],[mo_frag2],...]
avas_threshold             = "1D-14"   # threshold to truncated the active space in AVAS.

#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

with open(xyz_directory + "water.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   geometry = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

fragments = []
list_atom = []
for token in geometry:
    list_atom.append(token[0]) if token[0] not in list_atom else list_atom
for token in list_atom:
    fragments += [[(str(token),(0.0,0.0,0.0))]] # No need to specify coordinates for atomic fragments.

charge_frag = [0,0]
multiplicity_frag = [1,2]

#--------------------------------------------------------------------------
#-------------------------- END OF MANUAL SET UP --------------------------
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#------------------------------ CREATE INPUT ------------------------------
#--------------------------------------------------------------------------

with open("INPUT_GENIBO","w") as f:
 f.write("**ROSE\n")
 f.write(".VERSION\n")
 f.write(version+"\n")
 f.write(".CHARGE\n")
 f.write(str(charge)+"\n")
 f.write(".EXPONENT\n")
 f.write(str(exponent)+"\n")
 if not restricted:
    f.write(".UNRESTRICTED\n")
 f.write(".FILE_FORMAT\n")
 f.write("h5\n")
 if test == 1:
    f.write(".TEST  \n")
 if include_core:
    f.write(".INCLUDE_CORE\n")
 f.write(".ADDITIONAL_VIRTUALS_CUTOFF\n")
 f.write(str(additional_virtuals_cutoff)+"\n")
 f.write(".FRAG_THRESHOLD\n")
 f.write(str(frag_threshold)+"\n")
 f.write("*END OF INPUT\n")

#--------------------------------------------------------------------------
#---------------------------- SCF CALCULATIONS ----------------------------
#--------------------------------------------------------------------------

molecule = moldata(geometry=geometry,
                   multiplicity=multiplicity,
                   charge=charge,
                   data_directory=data_directory)

molecule, scf_wfn = run_psi4_mocoeff(molecule,
                            basis=basis1,
                            relativistic=relativistic,
                            spherical=spherical,
                            restricted=restricted,
                            openshell=openshell,
                            fcidump=fcidump,
                            save=save,
                            uncontract=uncontract)

subprocess.check_call("mv -f {}.h5 CHECKPOINT.h5".format(molecule.name), shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)
if save:
 subprocess.check_call("mv -f {}.npy MOLECULE.npy".format(molecule.name), shell=True, cwd=data_directory)


frag_size_total = 0
for i in range(len(fragments)):
 frag_geo = fragments[i]
 fragment = moldata(geometry=frag_geo,
                     charge=charge_frag[i],
                     multiplicity=multiplicity_frag[i],
                     data_directory=data_directory)
 fragment, scf_wfn_frag = run_psi4_mocoeff(fragment,
                            basis=basis2,
                            relativistic=relativistic,
                            spherical=spherical,
                            restricted=restricted,
                            openshell=openshell,
                            save=False,
                            uncontract=uncontract)

 subprocess.check_call("mv -f {}.h5 {:03d}.h5".format(fragment.name,fragment.protons[0]), shell=True, cwd=data_directory)
 subprocess.check_call("mv -f {}.xyz {:03d}.xyz".format(fragment.name,fragment.protons[0]), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- STORE IBO WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

if save:
   nmo = replace_psi4_MOcoefficients(scf_wfn,'ribo')

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:
 print("Restart the calculation with the new wavefunction. Check the ouput file.")
 molecule = moldata(geometry=geometry,
                    multiplicity=multiplicity,
                    charge=charge,
                    data_directory=data_directory)

 molecule, wfn = run_psi4_mocoeff(molecule,
                             basis=basis1,
                             relativistic=relativistic,
                             spherical=spherical,
                             restricted=restricted,
                             openshell=openshell,
                             fcidump=fcidump,
                             restart_wfn=scf_wfn,
                             uncontract=uncontract)

#--------------------------------------------------------------------------
#------------------------ POST-SCF CALCULATIONS ---------------------------
#--------------------------------------------------------------------------

if run_postSCF:
   print()
   print("--------------------------------------------------------------------------")
   print("------------------------ POST-SCF CALCULATIONS ---------------------------")
   print("--------------------------------------------------------------------------")
   print()
if run_postSCF and not restricted: print("post SCF tests are for restricted calculations for now.")
if run_postSCF and restricted:

 original_wfn = psi4.core.Wavefunction.from_file("MOLECULE")
 ribo_wfn = psi4.core.Wavefunction.from_file("ribo")

 n_active = 4
 nocc_active = n_active // 2               # CAUTION: this needs adjustment for active spaces that are not constructed as bonding+antibonding
 n_frozen = scf_wfn.nalpha() - nocc_active # compute the number of doubly occupied orbitals by subtracting the active occupied space
 n_virtual = nmo - n_frozen - n_active     # the remaining orbitals are inactive virtuals


 psi4.set_options({'ACTIVE'      : [ n_active ],\
                   'FROZEN_DOCC' : [ n_frozen ],\
                   'FROZEN_UOCC' : [ n_virtual],\
                 })

 PSI4_CASCI_energy = psi4.energy('fci', ref_wfn=original_wfn,  return_wfn=False)
 print("CAS-CI energy (original):",PSI4_CASCI_energy)
 PSI4_CASCI_energy = psi4.energy('fci', ref_wfn=ribo_wfn,  return_wfn=False)
 print("CAS-CI energy (ribo    ):",PSI4_CASCI_energy)
