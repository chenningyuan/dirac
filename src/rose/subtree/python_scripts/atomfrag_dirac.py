from genibo_MOcoeff import run_dirac_mocoef, moldata, average_of_configuration_openshell
import os
import sys
import subprocess

#--------------------------------------------------------------------------
#----------------------------  SYSTEM TO WORK ON --------------------------
#--------------------------------------------------------------------------

system =  str(sys.argv[1])
charge =  int(sys.argv[2])

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
spatial_orbs      = True      # Spatial orbital is the same whatever the spin.
spinfree          = False     # Use Dyall's spin-free Hamiltonian to obtain results without spin-orbit coupling for the four- or two-component Hamiltonian.
basis1            = 'STO-3G'  # Basis set for the full problem.
basis2            = 'STO-3G'  # Basis set for the fragments
special_basis     = None      # Can be used to set another basis for some atoms. An example is ["STO-3G","H cc-pVDZ"]
symmetry          = False     # Use symmetry or not (GENIBO works without symmetry for now)
point_nucleus     = True      # Use the point-nucleus model
uncontract        = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
speed_of_light    = False     # Set the speed of light in a.u. (default is 137 a.u.)
test              = True      # Perform tests which will be printed in the terminal.
restart           = True      # use the --put=DFCOEF option of Dirac with the IAO and IBO coefficients.
get_mrconee       = False     # Extract the MRCONEE file containing the one-electron integrals, essential for the tests.
description       = ""        # Add an optional description to the default name of the files
avas_frag         = []        # Fragment IAO file to be extracted from ROSE (used for AVAS for instance). Keep empty otherwise.
nmo_avas          = []        # list of spatial MOs (or spinors if restricted = False) to consider in AVAS.
run_postSCF       = False     # Run post-SCF calculation (CASCI/CASSCF) with original, IBO and AVAS orbitals.

if spatial_orbs and (hamiltonian != "NONREL"): sys.exit("spatial_orbs = True cannot be set for a relativistic hamiltonian.")
if (not spatial_orbs) and (hamiltonian == "NONREL"): print("*** WARNING *** spatial_orbs = False shouldn't be set for a non relativistic hamiltonian if you want fchk files.")
if test: get_mrconee = True

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
openshell = []
openshell_fullmolecule = False
for token in geometry:
    list_atom.append(token[0]) if token[0] not in list_atom else list_atom
for token in list_atom:
    fragments += [[(str(token),(0.0,0.0,0.0))]] # No need to specify coordinates for atomic fragments.
    openshell += [average_of_configuration_openshell[str(token)]]
nfragments = len(geometry)

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
 if not spatial_orbs:
    f.write(".SPINORS\n")
 if len(avas_frag) > 0:
    f.write(".AVAS  \n")
    f.write(str(len(avas_frag))+"\n")
    f.writelines("{:3d}".format(item) for item in avas_frag)
 f.write("*END OF INPUT\n")

with open("INPUT_AVAS","w") as f:                                                            
    f.write(str(len(geometry)) + " # natoms\n")                                                 
    f.write(str(charge) + " # charge\n")                                                        
    f.write("1 # restricted\n") if restricted else f.write("0 # restricted\n")                  
    f.write("1 # spatial orbs\n") if spatial_orbs else f.write("0 # spatial_orbs\n")            
    f.write("dirac # MO file for the full molecule\n")                                         
    f.write("dirac # MO file for the fragments\n")                                             
    f.write(str(len(nmo_avas)) + " # number of valence MOs in B2\n")                            
    f.writelines("{:3d}".format(item) for item in nmo_avas)                                     

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
                          charge=0,
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

 subprocess.check_call("mv -f CHECKPOINT.h5 {:03d}.h5".format(fragment.protons[0]), shell=True, cwd=data_directory)
 subprocess.check_call("mv -f {}.xyz {:03d}.xyz".format(fragment.name,fragment.protons[0]), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call("mv -f supermolecule.h5 CHECKPOINT.h5", shell=True, cwd=data_directory)
subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:

 if os.path.exists(os.path.join(data_directory, "ribo.h5")):
    subprocess.check_call("mv -f ribo.h5 CHECKPOINT.h5", shell=True, cwd=data_directory)

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

#--------------------------------------------------------------------------
#--------------------------- AVAS CONSTRUCTION ----------------------------
#--------------------------------------------------------------------------

if len(avas_frag) != 0:
 subprocess.check_call("cp -f frag1_iao.h5 B2.h5", shell=True, cwd=data_directory)
 subprocess.check_call(avas_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#------------------------ POST-SCF CALCULATIONS ---------------------------
#--------------------------------------------------------------------------
if run_postSCF:
   print()
   print("--------------------------------------------------------------------------")
   print("------------------------ POST-SCF CALCULATIONS ---------------------------")
   print("--------------------------------------------------------------------------")
   print()
if run_postSCF and (hamiltonian!="NONREL" or not restricted): print("For now works for NONREL and restricted calculations")
if run_postSCF and hamiltonian=="NONREL" and restricted:
 with open("OUTPUT_AVAS","r") as f:
   token = f.read().splitlines()
   N_frozen = int(token[1].split()[-1])
   N_active_occ = int(token[2].split()[-1])
   N_active_vir = int(token[3].split()[-1])
   N_active = N_active_occ + N_active_vir
   N_virtual = int(token[4].split()[-1])

 subprocess.check_call("cp -f avas.h5 CHECKPOINT.h5", shell=True, cwd=data_directory)

 with open("DIRRCI.inp","w") as f:
   f.write("""**DIRAC
.4INDEX
.WAVE FUNCTION
**WAVE FUNCTION
#.SCF
.DIRRCI
*SCF
.MAXITR
500
**HAMILTONIAN
.{0}
**INTEGRALS
.NUCMOD
 1
*READIN
.UNCONTRACT
**MOLTRA
.ACTIVE
{1}..{2}
**MOLECULE
*CHARGE
.CHARGE
 {3}
*SYMMETRY
.NOSYM
*BASIS
.DEFAULT
{4}
*END OF INPUT

 &DIRECT  MAXITER=25 CONVERE=1.0D-8 CONVERR=1.0D-9 &END
 &RASORB  NELEC={5} NRAS1=0,0, NRAS2={6},{6}, MAXH1=0, MAXE3=0  &END
 &CIROOT IREPNA='A  0', NROOTS=1 &END
""".format(hamiltonian,N_frozen+1,N_frozen+N_active,charge,basis1,N_active_occ*2,N_active))

 subprocess.check_call("pam --inp=DIRRCI.inp --mol=" + str(xyz_file) + " --put=CHECKPOINT.h5 --silent --noarch",shell=True,cwd=data_directory)
