from genibo_MOcoeff import *
import os
import sys
import subprocess
import numpy as np
import pyscf

#--------------------------------------------------------------------------
#----------------------------  SYSTEM TO WORK ON --------------------------
#--------------------------------------------------------------------------

system = str(sys.argv[1])
try:
  charge_str = sys.argv[2][1:len(sys.argv[2])-1].split(',')
  multiplicity_str = sys.argv[3][1:len(sys.argv[3])-1].split(',')
  charge_list = []
  multiplicity_list = []
  for i in charge_str: charge_list.append(int(i))
  for i in multiplicity_str: multiplicity_list.append(int(i))
except:
  sys.exit("execution should be: python atomfrags.py xyz_name [charge_mol,charge_frag0,charge_frag1,...] [multiplicity_mol,multiplicity_frag0,multiplicity_frag1,...]")

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory    = os.getenv('ROSEDIR')
rose_executable   = rose_directory + "/bin/rose.x"
plams_executable  = "plams"
data_directory    = os.getcwd()
xyz_directory     = rose_directory + "/test/xyz_files/"
xyz_file          = xyz_directory + system + ".xyz"
hamiltonian       = ["NONREL","X2C","X2Cmmf"][0]
version           = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent          = [2,3,4][2]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted        = True      # Change from unrestricted to restricted formalism.
relativistic      = False     # ScalarX2C
basis1            = 'STO-3G'  # Basis set for the full problem.
basis2            = 'STO-3G'  # Basis set for the fragments
charge            = charge_list[0]
multiplicity      = multiplicity_list[0]
spherical         = False     # Spherical or Cartesian coordinate GTOs. If True, it will do the calculation using spherical and store information using cartesians.
uncontract        = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
openshell         = True      # Openshell or closed-shell system
test              = True      # Perform tests which will be printed in the terminal.
restart           = True      # Restart calculation with IAO and IBO.
save              = True      # Save psi4 wavefunction into .npy file.

#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

with open(xyz_file,"r") as f:
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
nfragments = len(geometry)

charge_frag = charge_list[1:]
multiplicity_frag = multiplicity_list[1:]

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
 f.write("*END OF INPUT\n")
#--------------------------------------------------------------------------
#---------------------------- SCF CALCULATIONS ----------------------------
#--------------------------------------------------------------------------

molecule = moldata(geometry=geometry,
                        multiplicity=multiplicity,
                        charge=charge,
                        data_directory=data_directory)

molecule, mol, original_wfn = run_pyscf_mocoeff(molecule,
                            basis=basis1,
                            relativistic=relativistic,
                            spherical=spherical,
                            restricted=restricted,
                            openshell=openshell,
                            save=save,
                            uncontract=uncontract)

subprocess.check_call("mv -f {}.h5 CHECKPOINT.h5".format(molecule.name), shell=True, cwd=data_directory)
if save:
 subprocess.check_call("mv -f {}.chk MOLECULE.chk".format(molecule.name), shell=True, cwd=data_directory)

frag_size_total = 0
for i in range(len(fragments)):
 frag_geo = fragments[i]
 fragment = moldata(geometry=frag_geo,
                         multiplicity=multiplicity_frag[i],
                         charge=charge_frag[i],
                         data_directory=data_directory)

 fragment, frag, scf_wfn_frag = run_pyscf_mocoeff(fragment,
                             basis=basis2,
                             relativistic=relativistic,
                             spherical=spherical,
                             openshell=openshell,
                             restricted=restricted,
                             uncontract=uncontract)

 subprocess.check_call("mv -f {}.h5 {:03d}.h5".format(fragment.name,fragment.protons[0]), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- STORE IBO WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

if save:
   nmo = replace_pyscf_MOcoefficients(original_wfn,mol,'ribo')

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:
 restart_wfn = data_directory + "/ribo" # to get the ribo.chk file
 molecule = moldata(geometry=geometry,
                         multiplicity=multiplicity,
                         charge=charge,
                         data_directory=data_directory)
 
 molecule, new_mol, new_wfn = run_pyscf_mocoeff(molecule,
                             basis=basis1,
                             relativistic=relativistic,
                             spherical=spherical,
                             restricted=restricted,
                             save=save,
                             restart_wfn=restart_wfn,
                             openshell=openshell,
                             uncontract=uncontract)
