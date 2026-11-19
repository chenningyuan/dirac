from genibo_MOcoeff import *
import os
import sys
import subprocess
import numpy as np
import psi4

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
basis1                     = 'sto-3g'  # Basis set for the full problem.
basis2                     = 'sto-3g'  # Basis set for the fragments
charge                     = 0         # Charge of the full problem
multiplicity               = 1         # Multiplicity of the full problem
spherical                  = False     # Spherical or Cartesian coordinate GTOs. If True, it will do the calculation using spherical and store information using cartesians.
uncontract                 = False     # Decontract basis sets.
openshell                  = True      # Openshell or closed-shell system
test                       = True      # Perform tests which will be printed in the terminal.
restart                    = True      # Restart calculation with IAO and IBO.
fcidump                    = False     # Extract FCIDUMP file from PSI4
description                = ""        # Add an optional description to the default name of the files
save                       = True      # Save psi4 wavefunction into .npy file.
include_core               = False     # Include core in the localization procedure, default : False
additional_virtuals_cutoff = 2.0       # Defaut is 2.0 hartrees
frag_threshold             = 10.0       # Energy Cutoff to determine number of reference fragment hard virtuals to use
run_postSCF                = True      # Run post-SCF calculation (CASCI/CASSCF) with original, IBO and AVAS orbitals.
avas_frag                  = [0,1]       # Fragment IAO file to be extracted from ROSE (used for AVAS for instance). Keep empty otherwise.
mo_avas_frag               = [[4,5],[4,5]] # list of MOs (starting at 1) to consider in AVAS. [[mo_frag1],[mo_frag2],...]
avas_threshold             = "1D-1"   # threshold to truncated the active space in AVAS.

#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

with open(xyz_directory+"water-ammonia.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   geometry = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

with open(xyz_directory+"water-ammonia_frag0.xyz","r") as f:
   f.readline()
   f.readline()
   frag1 = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

with open(xyz_directory+"water-ammonia_frag1.xyz","r") as f:
   f.readline()
   f.readline()
   frag2 = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

fragments = [frag1,frag2]
nfragments = len(fragments)

charge_frag = [0,0]
multiplicity_frag = [1,1]

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
 f.write(".NFRAGMENTS\n")
 f.write(str(nfragments)+"\n")
 if include_core:
    f.write(".INCLUDE_CORE\n")
 f.write(".ADDITIONAL_VIRTUALS_CUTOFF\n")
 f.write(str(additional_virtuals_cutoff)+"\n")
 f.write(".FRAG_THRESHOLD\n")
 f.write(str(frag_threshold)+"\n")
 for frag, mos in zip(avas_frag,mo_avas_frag):
     f.write(".FRAG_AVAS\n")
     f.write(str(frag+1)+"\n")
     f.write(str(len(mos))+"\n")
     f.writelines("{:3d}\n".format(mo) for mo in mos)
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

subprocess.check_call("mv -f {}.h5 MOLECULE.h5".format(molecule.name), shell=True, cwd=data_directory)
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

 subprocess.check_call("mv -f {}.h5 frag{:d}.h5".format(fragment.name,i), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- STORE IBO WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

if save:
   nmo = replace_psi4_MOcoefficients(scf_wfn,'ibo')

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
#--------------------------- AVAS CONSTRUCTION ----------------------------
#--------------------------------------------------------------------------

if len(avas_frag) != 0:
 nmo = replace_psi4_MOcoefficients(scf_wfn,'avas')

#--------------------------------------------------------------------------
#------------------------ POST-SCF CALCULATIONS ---------------------------
#--------------------------------------------------------------------------

if run_postSCF:
   print()
   print("--------------------------------------------------------------------------")
   print("------------------------ POST-SCF CALCULATIONS ---------------------------")
   print("--------------------------------------------------------------------------")
   print()
if run_postSCF and not restricted: print("Reference UHF for DETCI is not available.")
if run_postSCF and restricted:

 original_wfn = psi4.core.Wavefunction.from_file("MOLECULE")
 ibo_wfn = psi4.core.Wavefunction.from_file("ibo")
 avas_wfn = psi4.core.Wavefunction.from_file("avas")

 n_active = len(sum(mo_avas_frag, []))     # the sum function does just flatten this list of lists
 nocc_active = n_active // 2               # CAUTION: this needs adjustment for active spaces that are not constructed as bonding+antibonding
 n_frozen = scf_wfn.nalpha() - nocc_active # compute the number of doubly occupied orbitals by subtracting the active occupied space
 n_virtual = nmo - n_frozen - n_active     # the remaining orbitals are inactive virtuals

 psi4.set_options({'ACTIVE'      : [ n_active ],\
                   'FROZEN_DOCC' : [ n_frozen ],\
                   'FROZEN_UOCC' : [ n_virtual],\
                 })

 PSI4_CASCI_energy = psi4.energy('fci', ref_wfn=original_wfn,  return_wfn=False)
 print("CAS-CI energy (original):",PSI4_CASCI_energy)
 PSI4_CASCI_energy = psi4.energy('fci', ref_wfn=ibo_wfn,  return_wfn=False)
 print("CAS-CI energy (ibo     ):",PSI4_CASCI_energy)
 PSI4_CASCI_energy = psi4.energy('fci', ref_wfn=avas_wfn,  return_wfn=False)
 print("CAS-CI energy (avas    ):",PSI4_CASCI_energy)
