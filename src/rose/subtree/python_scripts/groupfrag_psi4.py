# Atomic IAO fragments, grouped into coarser recanonization groups.
#
# Usage:  python groupfrag_psi4.py <system>:<groups_file> [charges] [multiplicities]
#
# <system>      names an xyz file in test/xyz_files, as for atomfrag_psi4.py.
# <groups_file> names a labelled xyz file in the same directory: the same geometry
#               with a fifth field on each atom line giving the name of the group
#               that atom belongs to. Atoms sharing a label are recanonized together,
#               while the IAOs themselves stay atomic.
#
# The labelled file supplies labels only. Its coordinates are used to match lines
# against the centers of the AO basis and are then discarded, so it cannot disagree
# with the orbitals about geometry - it can only fail to match, which ROSE reports.

from genibo_MOcoeff import *
import os
import sys
import subprocess
import numpy as np
import psi4

#--------------------------------------------------------------------------
#----------------------------  SYSTEM TO WORK ON --------------------------
#--------------------------------------------------------------------------

if ':' not in str(sys.argv[1]):
  sys.exit("execution should be: python groupfrag_psi4.py xyz_name:groups_xyz_name [charges] [multiplicities]")
system, groups_name = str(sys.argv[1]).split(':', 1)
try:
  charge_str = sys.argv[2][1:len(sys.argv[2])-1].split(',')
  multiplicity_str = sys.argv[3][1:len(sys.argv[3])-1].split(',')
  charge_list = []
  multiplicity_list = []
  for i in charge_str: charge_list.append(int(i))
  for i in multiplicity_str: multiplicity_list.append(int(i))
except:
  sys.exit("execution should be: python groupfrag_psi4.py xyz_name:groups_xyz_name [charge_mol,charge_frag0,...] [multiplicity_mol,multiplicity_frag0,...]")

#--------------------------------------------------------------------------
#----------------------------  FIXED PARAMETERS ---------------------------
#--------------------------------------------------------------------------

rose_directory    = os.getenv('ROSEDIR')
rose_executable   = rose_directory + "/bin/rose.x"
plams_executable  = "plams"
data_directory    = os.getcwd()
xyz_directory     = rose_directory + "/test/xyz_files/"
xyz_file          = xyz_directory + system + ".xyz"
groups_file       = xyz_directory + groups_name + ".xyz"
groups_local      = "FRAGMENT_GROUPS.XYZ"   # name the labelled file gets in the run directory
hamiltonian       = ["NONREL","X2C","X2Cmmf"][0]
version           = ["Stndrd_2013","Simple_2013","Simple_2014"][0] # Different versions for the localization. They are translated from Knizia's program but in complex algebra.
exponent          = [2,3,4][2]# Exponent used in the localization procedure (exponent = 2 for Pipek-Mezey. Knizia uses exponent = 4)
restricted        = True      # Change from unrestricted to restricted formalism.
relativistic      = False     # ScalarX2C
basis1            = 'STO-3G'  # Basis set for the full problem.
basis2            = 'STO-3G'  # Basis set for the fragments
charge            = charge_list[0]
multiplicity      = multiplicity_list[0]
spherical         = False     # Spherical or Cartesian coordinate GTOs. If True, it will do the calculation using spherical and store information using cartesians.
uncontract        = True      # Decontract basis sets. (Two-component quasirelativistic Hamiltonians (like X2C) work only with decontracted basis set)
openshell         = True      # Openshell or closed-shell system
test              = True      # Perform tests which will be printed in the terminal.
restart           = True      # Restart calculation with IAO and IBO.
fcidump           = False     # Extract FCIDUMP file from PSI4
description       = ""        # Add an optional description to the default name of the files
save              = True      # Save psi4 wavefunction into .npy file.

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
for token in geometry:
    list_atom.append(token[0]) if token[0] not in list_atom else list_atom
for token in list_atom:
    fragments += [[(str(token),(0.0,0.0,0.0))]] # No need to specify coordinates for atomic fragments.
nfragments = len(geometry)

charge_frag = charge_list[1:]
multiplicity_frag = multiplicity_list[1:]

# One charge/multiplicity is needed per DISTINCT element, in order of first
# appearance in the xyz file, after the entry for the molecule itself.
if len(charge_frag) < len(fragments) or len(multiplicity_frag) < len(fragments):
   sys.exit("{} distinct elements ({}) need that many entries after the molecule in "
            "both lists; got {} charges and {} multiplicities."
            .format(len(fragments), ','.join(list_atom),
                    len(charge_frag), len(multiplicity_frag)))

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
 f.write(".FRAGMENT_XYZ\n")
 f.write(groups_local+"\n")
 f.write("*END OF INPUT\n")

# Put the labelled xyz next to the other files ROSE reads from the run directory.
import shutil
if not os.path.isfile(groups_file):
   sys.exit("could not find the labelled xyz file: " + groups_file)
shutil.copyfile(groups_file, os.path.join(data_directory, groups_local))

#--------------------------------------------------------------------------
#---------------------------- SCF CALCULATIONS ----------------------------
#--------------------------------------------------------------------------

molecule = moldata(geometry=geometry,
                        multiplicity=multiplicity,
                        charge=charge,
                        data_directory=data_directory)

molecule, scf_wfn = run_psi4_mocoeff(molecule,
                            basis=basis1,
                            relativistic=relativistic,
                            spherical=spherical,
                            restricted=restricted,
                            openshell=openshell,
                            fcidump=fcidump,
                            save=save,
                            uncontract=uncontract)

subprocess.check_call("mv -f {}.h5 CHECKPOINT.h5".format(molecule.name), shell=True, cwd=data_directory)
subprocess.check_call("mv -f {}.xyz MOLECULE.XYZ".format(molecule.name), shell=True, cwd=data_directory)
if save:
 subprocess.check_call("mv -f {}.npy MOLECULE.npy".format(molecule.name), shell=True, cwd=data_directory)

# The atomic fragments are either taken from a library of pre-made <Z>.h5 files or
# generated here with psi4. Point ROSE_ATOMFRAG_LIBRARY at the directory holding the
# library to use the first route, which is what the regression tests do.
atomfrag_library = os.environ.get("ROSE_ATOMFRAG_LIBRARY")
if atomfrag_library is not None and not os.path.isdir(atomfrag_library):
   print(" WARNING: the atomic fragment library")
   print("             " + atomfrag_library)
   print("          does not exist.")
   atomfrag_library = None
if atomfrag_library is None:
   print(" WARNING: the atomic fragments are generated with psi4. An open-shell atomic SCF")
   print("          without fractional occupations does not yield spherically symmetric")
   print("          orbitals, and which of the degenerate solutions it converges to may")
   print("          differ between machines. The resulting reference orbitals can introduce")
   print("          artifacts in the localized orbitals and make results hard to reproduce.")
   print("          Set ROSE_ATOMFRAG_LIBRARY to a directory of pre-made <Z>.h5 atoms,")
   print("          for instance the DIRAC-generated ones under examples/PSI4, to avoid this.")

frag_size_total = 0
for i in range(len(fragments)):
 frag_geo = fragments[i]
 fragment = moldata(geometry=frag_geo,
                          charge=charge_frag[i],
                          multiplicity=multiplicity_frag[i],
                          data_directory=data_directory)

 if atomfrag_library is not None:

    # Take this element from the library. A missing file is an error rather than a
    # fallback: mixing library atoms with psi4-generated ones in a single run would
    # be worse than either choice on its own.
    source = os.path.join(atomfrag_library, "{:03d}.h5".format(fragment.protons[0]))
    if not os.path.isfile(source):
       sys.exit("could not find the atomic fragment file: " + source)
    target = os.path.join(data_directory, "{:03d}.h5".format(fragment.protons[0]))
    if not (os.path.exists(target) and os.path.samefile(source, target)):
       shutil.copyfile(source, target)

 else:

    # Generate this element with psi4. Catch an impossible charge/multiplicity pairing
    # here, rather than several frames deep inside psi4's molecule parser. n unpaired
    # electrons need multiplicity n+1, so the electron count and multiplicity-1 must
    # share parity.
    nunpaired = fragment.multiplicity - 1
    if nunpaired > fragment.n_electrons or (fragment.n_electrons - nunpaired) % 2 != 0:
       sys.exit("fragment {} ({}) has {} electrons, which cannot form multiplicity {}. "
                "Atomic fragments take their own ground-state multiplicity: H=2, C=3, N=4, O=3."
                .format(i+1, list_atom[i], fragment.n_electrons, fragment.multiplicity))

    fragment, scf_wfn_frag = run_psi4_mocoeff(fragment,
                               basis=basis2,
                               relativistic=relativistic,
                               spherical=spherical,
                               restricted=restricted,
                               openshell=openshell,
                               save=False,
                               uncontract=uncontract)

    subprocess.check_call("mv -f {}.h5 {:03d}.h5".format(fragment.name,fragment.protons[0]), shell=True, cwd=data_directory)
    subprocess.check_call("mv -f {}.xyz {:03d}.xyz".format(fragment.name,fragment.protons[0]), shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- IAO/IBO CONSTRUCTION --------------------------
#--------------------------------------------------------------------------

subprocess.check_call(rose_executable, shell=True, cwd=data_directory)

#--------------------------------------------------------------------------
#-------------------------- STORE IBO WAVEFUNCTION ------------------------
#--------------------------------------------------------------------------

if save:
   nmo = replace_psi4_MOcoefficients(scf_wfn,'ribo')

#--------------------------------------------------------------------------
#------------------------------ CHECK HF ENERGY ---------------------------
#--------------------------------------------------------------------------

if restart:
 print("Restart the calculation with the new wavefunction. Check the ouput file.")
 molecule = moldata(geometry=geometry,
                         multiplicity=multiplicity,
                         charge=charge,
                         data_directory=data_directory)

 molecule, wfn = run_psi4_mocoeff(molecule,
                             basis=basis1,
                             relativistic=relativistic,
                             spherical=spherical,
                             restricted=restricted,
                             openshell=openshell,
                             fcidump=fcidump,
                             restart_wfn=scf_wfn,
                             uncontract=uncontract)
