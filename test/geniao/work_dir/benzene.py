from gen_dfcoef import run_dirac_dfcoef, MolData_Dirac
import os
import subprocess

#----------------------------  MANUAL SET UP -----------------------------
basis1 = 'cc-pVDZ'
basis2 = 'cc-pVQZ'
charge = 0
geniao_executable="geniao.x"
data_directory=os.getcwd()
geometry = [
('C',(      -1.39880212,       0.00000000,       0.00000000)),\
('C',(       1.39880212,       0.00000000,       0.00000000)),\
('C',(      -0.69940106,      -1.21139817,       0.00000000)),\
('C',(      -0.69940106,       1.21139817,       0.00000000)),\
('C',(       0.69940106,       1.21139817,       0.00000000)),\
('C',(       0.69940106,      -1.21139817,       0.00000000)),\
('H',(      -2.49009098,       0.00000000,       0.00000000)),\
('H',(      -1.24504549,       2.15648204,       0.00000000)),\
('H',(       1.24504549,       2.15648204,       0.00000000)),\
('H',(       2.49009098,       0.00000000,       0.00000000)),\
('H',(       1.24504549,      -2.15648204,       0.00000000)),\
('H',(      -1.24504549,      -2.15648204,       0.00000000))
]
fragment_list = ['H','C'] # from the lowest atomic number to the highest
point_nucleus = True
version = ["Stndrd_2013","Simple_2013","Simple_2014"][0]
exponent=[2,3,4][2]
get_mrconee=True
get_dfpcmo=True
test=True
hamiltonian = "X2C"
spinfree=True
speed_of_light = False
restricted = False
# by default :
#symmetry = False
#unconstract = True
#speed_of_light = False
#point_nucleus = False
#version = "Stndrd_2013"
#restart = False

#-------------------------- END OF MANUAL SET UP -------------------------


#------------------------- CREATE THE DFCOEF FILES -----------------------
molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               hamiltonian=hamiltonian,
                               speed_of_light=speed_of_light,
                               data_directory=data_directory)

# set get_mrconee=True in run_dirac_dfcoef if you want to extract the MRCONEE, used for testing the HF energy.
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     exponent=exponent,
                     speed_of_light=speed_of_light,
                     get_mrconee=get_mrconee,
                     get_dfpcmo=get_dfpcmo,
                     test=test,
                     hamiltonian=hamiltonian,
                     spinfree=spinfree,
                     restricted=restricted,
                     point_nucleus=point_nucleus)

subprocess.check_call("mv -f DFCOEF MOLECULE.dfcoef", shell=True, cwd=data_directory)
subprocess.check_call("mv -f DFPCMO MOLECULE.dfpcmo", shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)
subprocess.check_call("mv -f INPUT_IAO INPUT", shell=True, cwd=data_directory)

print('Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

j = 0
charge0=0
for i in range(len(geometry)):
 if (j < len(fragment_list)):
   if (molecule.atoms[i] == fragment_list[j]): # avoid to calculate several time the same fragment...
     j = j+1
     frag_geo = [(molecule.atoms[i],(0.,0.,0.))]
     fragment = MolData_Dirac(geometry=frag_geo,
                                    basis=basis2,
                                    charge=charge0,
                                    speed_of_light=speed_of_light,
                                    hamiltonian=hamiltonian,
                                    data_directory=data_directory)
     fragment = run_dirac_dfcoef(fragment,
                          version=version,
                          hamiltonian=hamiltonian,
                          spinfree=spinfree,
                          exponent=exponent,
                          speed_of_light=speed_of_light,
                          point_nucleus=point_nucleus)
     subprocess.check_call("mv -f DFCOEF {:03d}.dfcoef".format(molecule.protons[i]), shell=True, cwd=data_directory)
     subprocess.check_call("mv -f {}.xyz {}.xyz".format(fragment.name,molecule.atoms[i]), shell=True, cwd=data_directory)

#-------------------------  IAO CONSTRUCTION -----------------------
subprocess.check_call(geniao_executable, shell=True, cwd=data_directory)

#-------------------------  CHECK HF ENERGY -----------------------
# For IAO
subprocess.check_call("cp -f iao.dfcoef DFCOEF", shell=True, cwd=data_directory)
molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               hamiltonian=hamiltonian,
                               speed_of_light=speed_of_light,
                               data_directory=data_directory)
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     hamiltonian=hamiltonian,
                     spinfree=spinfree,
                     speed_of_light=speed_of_light,
                     exponent=exponent,
                     restart=True,
                     point_nucleus=point_nucleus)

print('(IAO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

# For IBO
subprocess.check_call("cp -f ibo.dfcoef DFCOEF", shell=True, cwd=data_directory)

molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               hamiltonian=hamiltonian,
                               charge=charge,
                               speed_of_light=speed_of_light,
                               data_directory=data_directory)
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     hamiltonian=hamiltonian,
                     spinfree=spinfree,
                     speed_of_light=speed_of_light,
                     exponent=exponent,
                     restart=True,
                     point_nucleus=point_nucleus)

print('(IBO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))
