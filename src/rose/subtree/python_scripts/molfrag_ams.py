from genibo_MOcoeff import moldata, generate_adf_input
import os
import sys
import subprocess
import shutil

#--------------------------------------------------------------------------
#----------------------------  SYSTEM TO WORK ON --------------------------
#--------------------------------------------------------------------------

try:
  systems          = sys.argv[1][1:len(sys.argv[1])-1].split(',')
  charge_str       = sys.argv[2][1:len(sys.argv[2])-1].split(',')
  multiplicity_str = sys.argv[3][1:len(sys.argv[3])-1].split(',')
  xyz_list = []
  charge_list = []
  multiplicity_list = []
  for i in charge_str: charge_list.append(int(i))
  for i in multiplicity_str: multiplicity_list.append(int(i))
except:
  sys.exit("execution should be: python molfrags.py [mol,frag0,frag1,...] [charge_mol,charge_frag0,charge_frag1,...]")

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
charge            = charge_list[0]
multiplicity      = multiplicity_list[0]
openshell         = False     # Openshell or closed-shell system
spin              = 0         # Multiplicity of the full problem
restart           = True      # Restart calculation with IBO.
test              = False     # Perform tests which will be printed in the terminal.
avas_frag         = []        # Fragment IAO file to be extracted from ROSE (used for AVAS for instance). Keep empty otherwise.
mo_avas_frag      = [[],[]]   # list of MOs (starting at 1) to consider in AVAS. [[mo_frag1],[mo_frag2],...]
avas_threshold    = "1D-1"    # threshold to truncate the active space in AVAS.

#--------------------------------------------------------------------------
#------------------------ REPLAY THE STORED INPUTS ------------------------
#--------------------------------------------------------------------------
# By default this test does not run AMS at all. Everything ROSE needs was
# converted once from an ADF calculation with python_scripts/rose_kf_bridge.py
# and stored under test/ams_molfrags/input/<case>/: CHECKPOINT.h5 for the
# supermolecule, one h5 per molecular fragment (frag1.h5, frag2.h5, ...), and
# the INPUT_GENIBO that selects them. Keeping the ADF run out of the test suite
# is the whole point of step 3 of the restructuring -- neither the build nor
# the tests may depend on a licensed code.
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

stored_inputs = os.path.join(rose_directory,"test","ams_molfrags","input",systems[0])

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
multiplicity_frag = multiplicity_list[1:]

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
 if test == 1:
    f.write(".TEST  \n")
 f.write(".NFRAGMENTS\n")
 f.write(str(nfragments)+"\n")
 f.write(".AVAS_THRESHOLD \n")
 f.write(str(avas_threshold)+"\n")
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
                       charge=charge,
                       multiplicity=multiplicity,
                       data_directory=data_directory)

frags_list = []
for i in range(nfragments):
 fragment = moldata(geometry=fragments[i],
                          charge=charge_frag[i],
                          multiplicity=multiplicity_frag[i],
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
