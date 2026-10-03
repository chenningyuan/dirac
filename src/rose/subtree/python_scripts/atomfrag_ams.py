from genibo_MOcoeff import moldata, generate_adf_input
import os
import sys
import subprocess
import shutil

#--------------------------------------------------------------------------
#----------------------------  SYSTEM TO WORK ON --------------------------
#--------------------------------------------------------------------------

system = str(sys.argv[1])
charge = int(sys.argv[2])

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory    = os.getenv('ROSEDIR')
rose_executable   = rose_directory + "/bin/rose.x"
ams_script        = rose_directory + "/python_scripts/AMS.py"
plams_executable  = "plams"
data_directory    = os.getcwd()
xyz_directory     = rose_directory + "/test/xyz_files/"
version           = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent          = [2,3,4][0]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted        = True      # Change from unrestricted to restricted formalism.
basis1            = "DZP"     # Basis set for the full problem.
basis2            = "DZP"     # Basis set for the fragments
openshell         = False     # Openshell or closed-shell system
spin              = 0         # Multiplicity of the full problem
restart           = True      # Restart calculation with IBO.
test              = False     # Perform tests which will be printed in the terminal.
avas_frag         = []        # Not working with AMS interface, keep empty
nmo_avas          = []        # Not working with AMS interface, keep empty

#--------------------------------------------------------------------------
#------------------------ REPLAY THE STORED INPUTS ------------------------
#--------------------------------------------------------------------------
# By default this test does not run AMS at all. Everything ROSE needs was
# converted once from an ADF calculation with python_scripts/rose_kf_bridge.py
# and stored under test/ams_atomfrags/input/<case>/: CHECKPOINT.h5 for the
# supermolecule, one h5 per atomic fragment named by atomic number, and the
# INPUT_GENIBO that selects them. Keeping the ADF run out of the test suite is
# the whole point of step 3 of the restructuring -- neither the build nor the
# tests may depend on a licensed code.
#
# The stored CHECKPOINT.h5 also carries ADF's own AO overlap, so every run
# re-checks ROSE's Slater overlaps against it; see check_ao_overlap.
#
# Set ROSE_AMS_RUN=1 to take the original path below instead, which drives ADF
# through PLAMS and writes rkf and TAPE15 files into the work directory. It needs
# AMS, and it does NOT refresh what is stored under input/: turning a fresh
# plams_workdir into a fixture is a separate and deliberate step, one
# rose_kf_bridge.py "forward" call for the supermolecule (with --t15, so that the
# AO overlap reference travels with it) and one per fragment.
#
# Nothing is printed here on purpose: the captured stdout has to be ROSE's
# alone, because it is what the stored reference is compared against.

stored_inputs = os.path.join(rose_directory,"test","ams_atomfrags","input",system)

if not os.getenv("ROSE_AMS_RUN"):
   if not os.path.isdir(stored_inputs):
      sys.exit("rose: no stored inputs in " + stored_inputs +
               "\n     set ROSE_AMS_RUN=1 to run ADF instead, then convert its output"
               "\n     with python_scripts/rose_kf_bridge.py")
   for name in sorted(os.listdir(stored_inputs)):
      shutil.copy(os.path.join(stored_inputs,name), data_directory)
   subprocess.check_call(rose_executable, shell=True, cwd=data_directory)
   sys.exit(0)


#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

with open(xyz_directory + system + ".xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   mol_geometry = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

frag_geometry = []
frag_charge   = []
frag_spin     = []
for token in mol_geometry:
    frag_geometry += [[token]]
    frag_charge   += [0]
    frag_spin     += [0]

nfragments = len(frag_geometry)

#--------------------------------------------------------------------------
#-------------------------- END OF MANUAL SET UP --------------------------
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#------------------------------ CREATE INPUT ------------------------------
#--------------------------------------------------------------------------

with open(data_directory + "/INPUT_GENIBO","w") as f:
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
 f.write("adf\n")
 f.write("*END OF INPUT\n")

#--------------------------------------------------------------------------
#---------------------------- SCF CALCULATIONS ----------------------------
#--------------------------------------------------------------------------

molecule = moldata(geometry=mol_geometry,
                       charge=charge,
                       multiplicity=int(2*spin)+1,
                       data_directory=data_directory)

frags_list = []
for i in range(nfragments):
 fragment = moldata(geometry=frag_geometry[i],
                          charge=frag_charge[i],
                          multiplicity=int(2*frag_spin[i])+1,
                          description=str(i),
                          data_directory=data_directory)
 frags_list.append(fragment)


#--------------------------------------------------------------------------
#-------GENERATE XYZ AND AMS INPUT TO RUN AMS.py WITH PLAMS ---------------
#--------------------------------------------------------------------------

generate_adf_input(molecule,basis1,frags_list,basis2,restricted,openshell,False)
subprocess.check_call(plams_executable + " " + ams_script, shell=True, cwd=data_directory) 

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory) 

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:

   # change restart=False to restart=True in INPUT_ADF:
   with open("INPUT_ADF", 'r') as f:
       lines = f.read().splitlines()
   lines[8+nfragments] = "True"
   with open("INPUT_ADF", 'w+') as f:
       f.write('\n'.join(lines))

   # Run plams AMS.py to restart the calculation
   subprocess.check_call(plams_executable + " " + ams_script, shell=True, cwd=data_directory)
