from gen_dfcoef import run_dirac_dfcoef, MolData_Dirac
import os
import subprocess

#----------------------------  MANUAL SET UP -----------------------------
basis1 = 'cc-pVDZ'
basis2 = 'cc-pVQZ'
charge = 0
geniao_executable="geniao.x"
data_directory="/Users/bsenjean/Documents/Git_repositories/dirac-caar-ibos/test/geniao/work_dir/"
data_directory=os.getcwd()
point_nucleus = False
geometry = [
('N',(-1.39559,      -0.02156,       0.00004    )),
('H',(-1.62981,       0.96110,      -0.10622    )),
('H',(-1.86277,      -0.51254,      -0.75597    )),
('H',(-1.83355,      -0.33077,       0.86231    )),
('O',(1.56850,       0.10589,       0.00001    )),
('H',(0.60674,      -0.03396,      -0.00063    )),
('H',(1.94052,      -0.78000,       0.00022    ))
]
fragment_list = ['H','N','O'] # from the lowest atomic number to the highest
version = ["Stndrd_2013","Simple_2013","Simple_2014"][0]
exponent = [2,3,4][2]
# by default :
#symmetry = False
#unconstract = True
#speed_of_light = False
#hamiltonian = False (means the default of Dirac, i.e. Dirac-Coulomb Hamiltonian using the Simple Coulombic Correction.
#point_nucleus = False
#version = "Stndrd_2013"
#restart = False

#-------------------------- END OF MANUAL SET UP -------------------------


#------------------------- CREATE THE DFCOEF FILES -----------------------
molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               data_directory=data_directory)

# set get_mrconee=True in run_dirac_dfcoef if you want to extract the MRCONEE, used for testing the HF energy.
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     exponent=exponent,
                     point_nucleus=point_nucleus)

subprocess.check_call("mv -f DFCOEF MOLECULE.dfcoef", shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)
subprocess.check_call("mv -f INPUT_IAO INPUT", shell=True, cwd=data_directory)

print('Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

charge0 = 0
j=0
for i in range(len(geometry)):
 if (j <= len(fragment_list)):
   if (molecule.atoms[i] == fragment_list[j]): # avoid to calculate several time the same fragment...
     j = j+1
     frag_geo = [(molecule.atoms[i],(0.,0.,0.))]
     fragment = MolData_Dirac(geometry=frag_geo,
                                    basis=basis2,
                                    charge=charge0,
                                    data_directory=data_directory)
     fragment = run_dirac_dfcoef(fragment,
                          exponent=exponent,
                          version=version,
                          point_nucleus=point_nucleus)
     subprocess.check_call("mv -f DFCOEF {:03d}.dfcoef".format(molecule.protons[i]), shell=True, cwd=data_directory)
     subprocess.check_call("mv -f {}.xyz {}.xyz".format(fragment.name,molecule.atoms[i]), shell=True, cwd=data_directory)

#-------------------------  IAO CONSTRUCTION -----------------------
subprocess.check_call(geniao_executable, shell=True, cwd=data_directory)

#-------------------------  CHECK HF ENERGY -----------------------
subprocess.check_call("cp -f iao.dfcoef DFCOEF", shell=True, cwd=data_directory)
molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               data_directory=data_directory)
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     exponent=exponent,
                     restart=True,
                     point_nucleus=point_nucleus) 

print('(IAO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

# For IBO
subprocess.check_call("cp -f ibo.dfcoef DFCOEF", shell=True, cwd=data_directory)

molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               data_directory=data_directory)
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     exponent=exponent,
                     restart=True,
                     point_nucleus=point_nucleus)

print('(IBO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))
