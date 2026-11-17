from genibo_MOcoeff import moldata, generate_adf_input
import os
import sys
import subprocess

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory    = os.getenv('ROSEDIR')
rose_executable   = rose_directory + "/bin/rose.x"
adf_script        = rose_directory + "/python_scripts/AMS.py"
plams_executable  = "$AMSBIN/plams"
data_directory    = os.getcwd()
xyz_directory     = data_directory + "/../xyz_files/"
version           = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent          = [2,3,4][0]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted        = True      # Change from unrestricted to restricted formalism.
basis1            = "DZP"     # Basis set for the full problem.
basis2            = "DZP"     # Basis set for the fragments
openshell         = False     # Openshell or closed-shell system
multiplicity      = 1         # Multiplicity of the full problem
charge            = 0
restart           = True      # Restart calculation with IBO.
test              = False     # Perform tests which will be printed in the terminal.

#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

with open(xyz_directory + "Irppy3.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   mol_geometry = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

frag_geometry = []

with open(xyz_directory + "Ir.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   frag_geometry += [[(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]]
   f.close()

with open(xyz_directory + "ppy1.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   frag_geometry += [[(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]]
   f.close()

with open(xyz_directory + "ppy2.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   frag_geometry += [[(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]]
   f.close()

with open(xyz_directory + "ppy3.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   frag_geometry += [[(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]]
   f.close()

nfragments        = 4
frag_charge       = [3,-1,-1,-1]
frag_multiplicity = [1, 1, 1, 1]

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
 f.write("*END OF INPUT\n")

#--------------------------------------------------------------------------
#---------------------------- SCF CALCULATIONS ----------------------------
#--------------------------------------------------------------------------

molecule = moldata(geometry=mol_geometry,
                   charge=charge,
                   multiplicity=multiplicity,
                   data_directory=data_directory)

frags_list = []
for i in range(nfragments):
 fragment = moldata(geometry=frag_geometry[i],
                    charge=frag_charge[i],
                    multiplicity=frag_multiplicity[i],
                    description=str(i),
                    data_directory=data_directory)
 frags_list.append(fragment)

#--------------------------------------------------------------------------
#-------GENERATE XYZ AND ADF INPUT TO RUN ADF.py WITH PLAMS ---------------
#--------------------------------------------------------------------------

generate_adf_input(molecule,basis1,frags_list,basis2,restricted,openshell,False)
subprocess.check_call(plams_executable + " " + adf_script, shell=True, cwd=data_directory)

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

   # Run plams ADF.py to restart the calculation
   subprocess.check_call(plams_executable + " " + adf_script, shell=True, cwd=data_directory)
