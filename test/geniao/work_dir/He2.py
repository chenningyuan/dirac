#/Users/bsenjean/psi4conda/bin/python

#from gen_dfcoef import run_dirac_dfcoef, MolData_Dirac
import gen_dfcoef
import os
import subprocess

#----------------------------  MANUAL SET UP -----------------------------
basis1 = 'STO-3G'
basis2 = 'CC-pVTZ'
bond_length = 1.7
charge = 0
data_directory="/Users/bsenjean/Documents/Git_repositories/dirac-caar-ibos/test/geniao/work_dir/"
geometry = [('He', (0., 0., 0.)), ('He', (0., 0., bond_length))]
point_nucleus = True
version = ["Stndrd_2013","Simple_2013","Simple_2014"][2]
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
molecule = gen_dfcoef.MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               charge=charge,
                               data_directory=data_directory)

molecule = gen_dfcoef.run_dirac_dfcoef(molecule,
                     version=version,
                     get_mrconee=True,
                     point_nucleus=point_nucleus)

subprocess.check_call("mv -f DFCOEF MOLECULE.dfcoef", shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)
subprocess.check_call("mv -f INPUT_IAO INPUT", shell=True, cwd=data_directory)

print('Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))

charge = 0
for i in range(len(geometry)):
   frag_geo = [(molecule.atoms[i],(0.,0.,0.))]
   fragment = gen_dfcoef.MolData_Dirac(geometry=frag_geo,
                                  basis=basis2,
                                  charge=charge,
                                  data_directory=data_directory)
   fragment = gen_dfcoef.run_dirac_dfcoef(fragment,
                        version=version,
                        point_nucleus=point_nucleus)
   subprocess.check_call("mv -f DFCOEF {:03d}.dfcoef".format(molecule.protons[i]), shell=True, cwd=data_directory)
   subprocess.check_call("mv -f {}.xyz {}.xyz".format(fragment.name,molecule.atoms[i]), shell=True, cwd=data_directory)

#-------------------------  IAO CONSTRUCTION -----------------------
subprocess.check_call("geniao.x", shell=True, cwd=data_directory)
subprocess.check_call("cp -f iao.dfcoef DFCOEF", shell=True, cwd=data_directory)


#-------------------------  CHECK HF ENERGY -----------------------
molecule = gen_dfcoef.MolData_Dirac(geometry=geometry,
                               basis=basis1,
                               data_directory=data_directory)
molecule = gen_dfcoef.run_dirac_dfcoef(molecule,
                     version=version,
                     restart=True,
                     point_nucleus=point_nucleus) 

print('Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_scf()[0],molecule.get_scf()[1]))
