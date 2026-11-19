from genibo_MOcoeff import moldata, generate_adf_input
import os
import sys
import subprocess

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory              = os.getenv('ROSEDIR')
rose_executable             = rose_directory + "/bin/rose.x"
adf_script                  = rose_directory + "/python_scripts/AMS.py"
plams_executable            = "$AMSBIN/plams"
data_directory              = os.getcwd()
xyz_directory               = data_directory + "/../xyz_files/"
version                     = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent                    = [2,3,4][0]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted                  = True      # Change from unrestricted to restricted formalism.
basis1                      = "DZP"     # Basis set for the full problem.
basis2                      = "DZP"     # Basis set for the fragments
openshell                   = False     # Openshell or closed-shell system
multiplicity                = 1         # Multiplicity of the full problem
charge                      = 0
restart                     = False      # Restart calculation with IBO.
test                        = False     # Perform tests which will be printed in the terminal.
avas_frag                   = []        # Not working with ADF interface, keep empty
nmo_avas                    = []        # Not working with ADF interface, keep empty
include_core                = True     # Include core in the localization procedure, default : False
#additional_virtuals_cutoff  = 2.0      # Energy Cutoff to determine number of hard virtuals to keep
#frag_threshold              = 10.0      # Energy Cutoff to determine number of reference fragment hard virtuals to use
frag_valence                = []     # Number of Valence orbitals for each fragment, [[fragment, number of orbitals]..], leave empy if want to use defaults 
frag_core                   = []     # Number of Core orbitals for each fragment, [[fragment, number of orbitals]..], leave empy if want to use defaults
#--------------------------------------------------------------------------
#----------------- INPUT GEOMETRY AND FRAGMENTATION -----------------------
#--------------------------------------------------------------------------

with open(xyz_directory + "C2H4_C2F4_d0.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   mol_geometry = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

frag_geometry = []

with open(xyz_directory + "C2H4_C2F4_d0_frag1.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   frag_geometry += [[(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]]
   f.close()

with open(xyz_directory + "C2H4_C2F4_d0_frag2.xyz","r") as f:
   n_atoms = f.readline()
   f.readline()
   frag_geometry += [[(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]]
   f.close()

nfragments        = 2
frag_charge       = [0, 0]
frag_multiplicity = [1, 1]

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
 if include_core:
    f.write(".INCLUDE_CORE\n")
# f.write(".ADDITIONAL_VIRTUALS_CUTOFF\n")
# f.write(str(additional_virtuals_cutoff)+"\n") 
# f.write(".FRAG_THRESHOLD\n")
# f.write(str(frag_threshold)+"\n")
 if (len(frag_valence) >= 1):
    for item in frag_valence:
       f.write(".FRAG_VALENCE\n")
       f.write(str(item[0])+"\n")
       f.write(str(item[1])+"\n")
 if (len(frag_core) >= 1):                                                                                          
    for item in frag_core:                                                                                 
       f.write(".FRAG_CORE\n")
       f.write(str(item[0])+"\n")
       f.write(str(item[1])+"\n")
 if (len(avas_frag) >= 1):      
     f.write(".AVAS  \n")
     f.write(str(len(avas_frag))+"\n")
     f.writelines("{:3d}".format(item) for item in avas_frag)
     f.write("\n")
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
