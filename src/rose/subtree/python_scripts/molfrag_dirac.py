from genibo_MOcoeff import run_dirac_mocoef, moldata, average_of_configuration_openshell
import os
import sys
import subprocess

#--------------------------------------------------------------------------
#----------------------------  SYSTEM TO WORK ON --------------------------
#--------------------------------------------------------------------------

try:
  systems          = sys.argv[1][1:len(sys.argv[1])-1].split(',')                                               
  charge_str       = sys.argv[2][1:len(sys.argv[2])-1].split(',')
  xyz_list = []
  charge_list = []
  for i in charge_str: charge_list.append(int(i))                                                               
except: 
  sys.exit("execution should be: python molfrags.py [mol,frag0,frag1,...] [charge_mol,charge_frag0,charge_frag1,...]")

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory    = os.getenv('ROSEDIR')
rose_executable   = rose_directory + "/bin/rose.x"
plams_executable  = "plams"
data_directory    = os.getcwd()
xyz_directory     = rose_directory + "/test/xyz_files/"
hamiltonian       = ["NONREL","X2C","X2Cmmf"][0]
version           = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent          = [2,3,4][2]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted        = True      # Change from unrestricted to restricted formalism.
spatial_orbs      = True      # Spatial orbital is the same whatever the spin.
spinfree          = False     # Use Dyall's spin-free Hamiltonian to obtain results without spin-orbit coupling for the four- or two-component Hamiltonian.
basis1            = 'STO-3G'  # Basis set for the full problem.
basis2            = 'STO-3G'  # Basis set for the fragments
special_basis     = None      # Can be used to set another basis for some atoms. An example is ["STO-3G","H cc-pVDZ"]
charge            = charge_list[0]
symmetry          = False     # Use symmetry or not (GENIBO works without symmetry for now)
point_nucleus     = True      # Use the point-nucleus model
uncontract        = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
speed_of_light    = False     # Set the speed of light in a.u. (default is 137 a.u.)
test              = True      # Perform tests which will be printed in the terminal.
restart           = True      # use the --put=DFCOEF option of Dirac with the IAO and IBO coefficients.
get_mrconee       = False     # Extract the MRCONEE file containing the one-electron integrals, essential for the tests.
description       = ""        # Add an optional description to the default name of the files

if spatial_orbs and (hamiltonian != "NONREL"): sys.exit("spatial_orbs = True cannot be set for a relativistic hamiltonian.")
if (not spatial_orbs) and (hamiltonian == "NONREL"): print("*** WARNING *** spatial_orbs = False shouldn't be set for a non relativistic hamiltonian if you want fchk files.")
if test: get_mrconee = True

#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

# Zeroth system is the (super)molecule. Read the xyz file that defines it.
xyz_file = xyz_directory + systems[0] + ".xyz"
with open(xyz_file,"r") as f:
   n_atoms = f.readline()
   f.readline()
   geometry = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

# Loop over the subsystems and also read their geometry
fragments = []
for subsystem in systems[1:]:
    xyz_file = xyz_directory + subsystem + ".xyz"
    with open(xyz_file,"r") as f:
       n_atoms = f.readline()
       f.readline()
       fragments.append([(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()])
       f.close()

nfragments = len(fragments)

charge_frag = charge_list[1:]
openshell = [False] * nfragments
openshell_fullmolecule = False

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
 f.write("dirac\n")
 if test == 1:
    f.write(".TEST  \n")
 f.write(".NFRAGMENTS\n")
 f.write(str(nfragments)+"\n")
 f.write("*END OF INPUT\n")

#--------------------------------------------------------------------------
#---------------------------- SCF CALCULATIONS ----------------------------
#--------------------------------------------------------------------------

molecule = moldata(geometry=geometry,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac_mocoef(molecule,
                     basis=basis1,
                     special_basis=special_basis,
                     symmetry=symmetry,
                     hamiltonian=hamiltonian,
                     spinfree=spinfree,
                     point_nucleus=point_nucleus,
                     uncontract=uncontract,
                     speed_of_light=speed_of_light,
                     get_mrconee=get_mrconee,
                     openshell=openshell_fullmolecule)

subprocess.check_call("mv -f CHECKPOINT.h5 supermolecule.h5", shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.xyz".format(molecule.name), shell=True, cwd=data_directory)

print('Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_DIRAC_scf_energy()[0],molecule.get_DIRAC_scf_energy()[1]))

frag_size_total = 0
for i in range(len(fragments)):
 frag_geo = fragments[i]
 fragment = moldata(geometry=frag_geo,
                          charge=charge_frag[i],
                          description=description,
                          data_directory=data_directory)
 fragment = run_dirac_mocoef(fragment,
                  basis=basis2,
                  special_basis=special_basis,
                  symmetry=symmetry,
                  hamiltonian=hamiltonian,
                  spinfree=spinfree,
                  point_nucleus=point_nucleus,
                  uncontract=uncontract,
                  speed_of_light=speed_of_light,
                  openshell=openshell[i])

 subprocess.check_call("mv -f CHECKPOINT.h5 frag{:d}.h5".format(i+1), shell=True, cwd=data_directory)
 subprocess.check_call("mv -f {}.xyz frag{:d}.xyz".format(fragment.name,i+1), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call("mv -f supermolecule.h5 CHECKPOINT.h5", shell=True, cwd=data_directory)
subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:

 subprocess.check_call("mv ribo.h5 CHECKPOINT.h5", shell=True, cwd=data_directory)

 molecule = moldata(geometry=geometry,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)
 molecule = run_dirac_mocoef(molecule,
                     basis=basis1,
                     special_basis=special_basis,
                     symmetry=symmetry,
                     hamiltonian=hamiltonian,
                     spinfree=spinfree,
                     point_nucleus=point_nucleus,
                     uncontract=uncontract,
                     speed_of_light=speed_of_light,
                     restart=restart,
                     openshell=openshell_fullmolecule)

 print(' (IBO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_DIRAC_scf_energy()[0],molecule.get_DIRAC_scf_energy()[1]))

 print('\n Calculation finished: modified MO coefficient files available in {}'.format(data_directory))
