#/Users/bsenjean/psi4conda/bin/python

from gen_dfcoef import run_dirac_dfcoef, MolData_Dirac
import os
import subprocess

#----------------------------  MANUAL SET UP -----------------------------
basis1 = 'STO-3G'
basis2 = 'cc-pVTZ'
bond_length = 1.7
charge = 0
geniao_executable="geniao.x"
data_directory=os.getcwd()
geometry = [('He', (0., 0., 0.)), ('He', (0., 0., bond_length))]
point_nucleus = True
version = ["Stndrd_2013","Simple_2013","Simple_2014"][2]
exponent = [2,3,4][0]
get_mrconee=True
test=True
hamiltonian="X2C"
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
                               data_directory=data_directory)

# set get_mrconee=True in run_dirac_dfcoef if you want to extract the MRCONEE, used for testing the HF energy.
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     exponent=exponent,
                     hamiltonian=hamiltonian,
                     get_mrconee=get_mrconee,
                     test=test,
                     point_nucleus=point_nucleus)

subprocess.check_call("mv -f DFCOEF MOLECULE.dfcoef", shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)
subprocess.check_call("mv -f INPUT_IAO INPUT", shell=True, cwd=data_directory)

print('Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

charge0 = 0
for i in range(len(geometry)):
   frag_geo = [(molecule.atoms[i],(0.,0.,0.))]
   fragment = MolData_Dirac(geometry=frag_geo,
                                  basis=basis2,
                                  charge=charge0,
                                  hamiltonian=hamiltonian,
                                  data_directory=data_directory)
   fragment = run_dirac_dfcoef(fragment,
                        version=version,
                        exponent=exponent,
                        hamiltonian=hamiltonian,
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
                               data_directory=data_directory)
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     hamiltonian=hamiltonian,
                     exponent=exponent,
                     restart=True,
                     point_nucleus=point_nucleus) 

print('(IAO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

# For IBO
subprocess.check_call("cp -f ibo.dfcoef DFCOEF", shell=True, cwd=data_directory)

molecule = MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               hamiltonian=hamiltonian,
                               data_directory=data_directory)
molecule = run_dirac_dfcoef(molecule,
                     version=version,
                     exponent=exponent,
                     restart=True,
                     hamiltonian=hamiltonian,
                     point_nucleus=point_nucleus)

print('(IBO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))
