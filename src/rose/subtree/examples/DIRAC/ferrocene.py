from genibo_MOcoeff import run_dirac_mocoef, moldata
import os
import sys
import subprocess

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory    = os.getenv('ROSEDIR')
rose_executable   = rose_directory + "/bin/rose.x"
plams_executable  = "plams"
data_directory    = os.getcwd()
xyz_directory     = rose_directory + "/examples/xyz_files/"
xyz_file          = xyz_directory + "Ferrocene.xyz"
hamiltonian       = ["NONREL","X2C","X2Cmmf"][0]
version           = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent          = [2,3,4][2]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted        = True      # Change from unrestricted to restricted formalism.
spatial_orbs      = False     # Spatial orbital is the same whatever the spin.
spinfree          = False     # Use Dyall's spin-free Hamiltonian to obtain results without spin-orbit coupling for the four- or two-component Hamiltonian.
basis1            = 'STO-3G'  # Basis set for the full problem.
basis2            = 'STO-3G'  # Basis set for the fragments
special_basis     = None      # Can be used to set another basis for some atoms. An example is ["STO-3G","H cc-pVDZ"]
charge            = 0         # Charge of the full problem
symmetry          = False     # Use symmetry or not (GENIBO works without symmetry for now)
point_nucleus     = True      # Use the point-nucleus model
uncontract        = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
speed_of_light    = False     # Set the speed of light in a.u. (default is 137 a.u.)
openshell         = True      # If True, please go to the INPUT GEOMETRY AND FRAGMENTATION part and specify openshells by hand ! Otherwise, closed shell are always considered.
test              = True      # Perform tests which will be printed in the terminal.
restart           = True      # use the --put=DFCOEF option of Dirac with the IAO and IBO coefficients.
get_mrconee       = False     # Extract the MRCONEE file containing the one-electron integrals, essential for the tests.
description       = ""        # Add an optional description to the default name of the files
run_postSCF       = True      # Run post-SCF calculation (CASCI) with original, IBO and AVAS orbitals.
avas_frag         = [1]       # Fragment IAO file to be extracted from ROSE (used for AVAS for instance). Keep empty otherwise.
mo_avas_frag      = [[6,7,8]]   # list of spatial MOs (or spinors if restricted = False) to consider in AVAS.
avas_threshold    = "1D-14"


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

with open(xyz_directory+"Ferrocene_frag0.xyz","r") as f:
   f.readline()
   f.readline()
   frag1 = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

with open(xyz_directory+"Ferrocene_frag1.xyz","r") as f:
   f.readline()
   f.readline()
   frag2 = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

with open(xyz_directory+"Ferrocene_frag2.xyz","r") as f:
   f.readline()
   f.readline()
   frag3 = [(line.split()[0],(line.split()[1],line.split()[2],line.split()[3])) for line in f.readlines()]
   f.close()

fragments = [frag1,frag2,frag3]
nfragments = len(fragments)
charge_frag = [2,-1,-1]

openshell = [False]*nfragments
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
 f.write(".AVAS_THRESHOLD \n")
 f.write(str(avas_threshold)+"\n")
 for frag, mos in zip(avas_frag,mo_avas_frag):
     f.write(".FRAG_AVAS\n")
     f.write(str(frag)+"\n")
     f.write(str(len(mos))+"\n")
     f.writelines("{:3d}\n".format(mo) for mo in mos)
 f.write("\n*END OF INPUT\n")
 
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

subprocess.check_call("mv -f CHECKPOINT.h5 MOLECULE.h5", shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)

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

 subprocess.check_call("mv -f CHECKPOINT.h5 frag{:d}.h5".format(i), shell=True, cwd=data_directory)
 subprocess.check_call("mv -f {}.xyz frag{:d}.xyz".format(fragment.name,i), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:
 subprocess.check_call("cp -f ibo.h5 CHECKPOINT.h5", shell=True, cwd=data_directory)

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

 print('(IBO) Hartree-Fock energy of {} Hartree in {} iterations'.format(molecule.get_DIRAC_scf_energy()[0],molecule.get_DIRAC_scf_energy()[1]))

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

 # The next three parameters depend on the outcome of the avas projection, but when selecting 3 occupied fragment orbitals
 # as reference on may expect an active space of 3 occupied and 3 virtual molecular orbitals, giving a total space of 6 orbitals.
 # The number of frozen orbitals is of course system-specific.
 N_frozen = 45
 N_active = 6
 N_active_occ = 3
 
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

 subprocess.check_call("pam --inp=DIRRCI.inp --mol=" + str(xyz_file) + " --incmo --silent --noarch",shell=True,cwd=data_directory)
