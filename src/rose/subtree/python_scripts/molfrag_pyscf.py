from genibo_MOcoeff import *
import os
import sys
import subprocess
import numpy as np
import pyscf

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
  sys.exit("execution should be: python molfrags.py [mol,frag0,frag1,...] [charge_mol,charge_frag0,charge_frag1,...] [multiplicity_mol,multiplicity_frag0,multiplicity_frag1,...]")

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory             = os.getenv('ROSEDIR')
rose_executable            = rose_directory + "/bin/rose.x"
plams_executable           = "plams"
data_directory             = os.getcwd()
xyz_directory              = rose_directory + "/examples/xyz_files/"
version                    = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent                   = [2,3,4][2]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted                 = True      # Change from unrestricted to restricted formalism.
relativistic               = False     # scalar X2C
basis1                     = 'STO-3G'  # Basis set for the full problem.
basis2                     = 'STO-3G'  # Basis set for the fragments
charge                     = charge_list[0]
multiplicity               = multiplicity_list[0]
spherical                  = False     # Spherical or Cartesian coordinate GTOs. If True, it will do the calculation using spherical and store information using cartesians.
uncontract                 = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
openshell                  = True      # Openshell or closed-shell system
test                       = True      # Perform tests which will be printed in the terminal.
restart                    = True      # Restart calculation with IAO and IBO.
save                       = True      # Save pyscf wavefunction into .chk file
include_core               = False     # Include core in the localization procedure, default : False
additional_virtuals_cutoff = 2.0       # Energy Cutoff to determine number of hard virtuals to keep
frag_threshold             = 10.0      # Energy Cutoff to determine number of reference fragment hard virtuals to use
run_postSCF                = True      # Run post-SCF calculation (CASCI/CASSCF) with original, IBO and AVAS orbitals.
avas_frag                  = [0,1]       # Fragment IAO file to be extracted from ROSE (used for AVAS for instance). Keep empty otherwise.
mo_avas_frag               = [[4,5],[4,5]] # list of MOs (starting at 1) to consider in AVAS. [[mo_frag1],[mo_frag2],...]
avas_threshold             = "1D-14"   # threshold to truncated the active space in AVAS.

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
 f.write("h5\n")
 if test == 1:
    f.write(".TEST  \n")
 f.write(".NFRAGMENTS\n")
 f.write(str(nfragments)+"\n")
 if include_core:
    f.write(".INCLUDE_CORE\n")
 f.write(".ADDITIONAL_VIRTUALS_CUTOFF\n")
 f.write(str(additional_virtuals_cutoff)+"\n")
 f.write(".FRAG_THRESHOLD\n")
 f.write(str(frag_threshold)+"\n")
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
                        multiplicity=multiplicity,
                        charge=charge,
                        data_directory=data_directory)

molecule, mol, original_wfn = run_pyscf_mocoeff(molecule,
                            basis=basis1,
                            relativistic=relativistic,
                            spherical=spherical,
                            restricted=restricted,
                            openshell=openshell,
                            save=save,
                            uncontract=uncontract)

subprocess.check_call("mv -f {}.h5 MOLECULE.h5".format(molecule.name), shell=True, cwd=data_directory)
if save:
 subprocess.check_call("mv -f {}.chk MOLECULE.chk".format(molecule.name), shell=True, cwd=data_directory)

frag_size_total = 0
for i in range(len(fragments)):
 frag_geo = fragments[i]
 fragment = moldata(geometry=frag_geo,
                         multiplicity=multiplicity_frag[i],
                         charge=charge_frag[i],
                         data_directory=data_directory)

 fragment, frag, scf_wfn_frag = run_pyscf_mocoeff(fragment,
                             basis=basis2,
                             relativistic=relativistic,
                             spherical=spherical,
                             openshell=openshell,
                             restricted=restricted,
                             uncontract=uncontract)

 subprocess.check_call("mv -f {}.h5 frag{:d}.h5".format(fragment.name,i+1), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

nmo, ribo_energies, ribo_coeff = fetch_MOcoefficients('ribo')
if len(avas_frag) != 0: nmo, avas_energies, avas_coeff = fetch_MOcoefficients('avas')

#--------------------------------------------------------------------------
#------------------------ POST-SCF CALCULATIONS ---------------------------
#--------------------------------------------------------------------------

if run_postSCF:
   from pyscf import mcscf
   print()
   print("--------------------------------------------------------------------------")
   print("------------------------ POST-SCF CALCULATIONS ---------------------------")
   print("--------------------------------------------------------------------------")
   print()
if run_postSCF and not restricted: print("post SCF tests are for restricted calculations for now.")
if run_postSCF and restricted:

 n_active = 2 * len(sum(mo_avas_frag, [])) # the sum function does just flatten this list of lists
 nocc_active = n_active // 2               # CAUTION: this needs adjustment for active spaces that are not constructed as bonding+antibonding
 n_active_elec = 2 * nocc_active

 casci = original_wfn.CASCI(n_active, n_active_elec)
 casci_original = casci.casci()[0]
 casci_ribo = casci.casci(ribo_coeff)[0]
 casci_avas = casci.casci(avas_coeff)[0]

 #casscf = original_wfn.CASSCF(N_active, N_active_elec)
 #casscf_original = casscf.kernel()[0]
 #casscf_ribo = casscf.kernel(ibo_coeff)[0]
 #casscf_avas = casscf.kernel(avas_coeff)[0]

 print("CASCI (original):",casci_original)
 print("CASCI (ribo    ):",casci_ribo)
 print("CASCI (avas    ):",casci_avas)
 #print("CASSCF(original):",casscf_original)
 #print("CASSCF(ibo     ):",casscf_ibo)
 #print("CASSCF(avas    ):",casscf_avas)

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

# At the end because the replace appears to change the original object
# deepcopy is preferable, but this gave an issue.

if restart:
 nmo, ribo_wfn = replace_pyscf_MOcoefficients(original_wfn,mol,'ribo')
 restart_wfn = data_directory + "/ribo" # to get the ribo.chk file
 molecule = moldata(geometry=geometry,
                         multiplicity=multiplicity,
                         charge=charge,
                         data_directory=data_directory)
 
 molecule, new_mol, new_wfn = run_pyscf_mocoeff(molecule,
                             basis=basis1,
                             relativistic=relativistic,
                             spherical=spherical,
                             restricted=restricted,
                             save=save,
                             restart_wfn=restart_wfn,
                             openshell=openshell,
                             uncontract=uncontract)

