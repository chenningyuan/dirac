"""Functions to prepare Dirac input and run calculations."""
from __future__ import absolute_import

import os
import re
import subprocess
import warnings
from ._utils import *

class SpecialBasisError(Exception):
    pass
class SpeedOfLightError(Exception):
    pass
class SpinFreeError(Exception):
    pass

average_of_configuration_openshell = {
    'H':["0","1","1/2"], 'He':False,
    'Li':["2","1","1/2"], 'Be':False,
    'B':["4","1","1/6"], 'C':["4","1","2/6"], 'N':["4","1","3/6"], 'O':["4","1","4/6"], 'F':["4","1","5/6"], 'Ne':False,
    'Na':["10","1","1/2"], 'Mg':False,
    'Al':["12","1","1/6"], 'Si':["12","1","2/6"], 'P':["12","1","3/6"], 'S':["12","1","4/6"], 'Cl':["12","1","5/6"], 'Ar':False,
    'K':["18","1","1/2"], 'Ca':False,
    'Sc':["20","1","1/10"], 'Ti':["20","1","2/10"], 'V':["20","1","3/10"], 'Cr':["20","1","4/10"], 'Mn':["20","1","5/10"],
    'Fe':["20","1","6/10"], 'Co':["20","1","7/10"], 'Ni':["20","1","8/10"], 'Cu':["20","1","9/10"], 'Zn':False,
    'Ga':["30","1","1/6"], 'Ge':["30","1","2/6"], 'As':["30","1","3/6"], 'Se':["30","1","4/6"], 'Br':["30","1","5/6"], 'Kr':False,
    'Rb':["36","1","1/2"], 'Sr':False,
    'Y':["38","1","1/10"], 'Zr':["38","1","2/10"], 'Nb':["38","1","3/10"], 'Mo':["36","2","5/10\\1/2"], 'Tc':["36","2","6/10\\1/2"],
    'Ru':["36","2","7/10\\1/2"], 'Rh':["36","2","8/10\\1/2"], 'Pd':False, 'Ag':["46","1","1/2"], 'Cd':False,
    'In':["48","1","1/6"], 'Sn':["48","1","2/6"], 'Sb':["48","1","3/6"], 'Te':["48","1","4/6"], 'I':["48","1","5/6"], 'Xe':False,
    'Cs':["54","1","1/2"], 'Ba':False,
    'La':["56","1","1/10"], 'Ce':["56","2","1/14\\1/10"], 'Pr':["56","1","3/14"], 'Nd':["56","1","4/14"], 'Pm':["56","1","5/14"], 'Sm':["56","1","6/14"], 'Eu':["56","1","7/14"],
    'Gd':["56","1","8/14"],'Tb':["56","1","9/14"], 'Dy':["56","1","10/14"], 'Ho':["56","1","11/14"], 'Er':["56","1","12/14"], 'Tm':["56","1","13/14"], 'Yb':False,
    'Lu':["70","1","1/10"], 'Hf':["70","1","2/10"], 'Ta':["70","1","3/10"], 'W':["70","1","4/10"], 'Re':["70","1","5/10"],
    'Os':["70","1","6/10"], 'Ir':["70","1","7/10"], 'Pt':["70","1","8/10"], 'Au':["70","1","9/10"], 'Hg':False,
    'Tl':["80","1","1/6"], 'Pb':["80","1","2/6"], 'Bi':["80","1","3/6"], 'Po':["80","1","4/6"], 'At':["80","1","5/6"], 'Rn':False,
    'Fr':["86","1","1/2"], 'Ra':False,
    'Ac':["88","1","1/10"], 'Th':["88","1","2/10"], 'Pa':["88","2","2/14\\1/10"], 'U':["88","2","3/14\\1/10"], 'Np':["88","2","4/14\\1/10"], 'Pu':["88","1","6/14"], 'Am':["88","1","7/14"], 'Cm':["88","1","8/14"],
    'Bk':["88","1","9/14"], 'Cf':["88","1","10/14"], 'Es':["88","1","11/14"], 'Fm':["88","1","12/14"], 'Md':["88","1","13/14"], 'No':False,
    'Lr':["102","1","1/10"], 'Rf':["102","1","2/10"], 'Db':["102","1","3/10"], 'Sg':["102","1","4/10"], 'Bh':["102","1","5/10"],
    'Hs':["102","1","6/10"], 'Mt':["102","1","7/10"], 'Ds':["102","1","8/10"], 'Rg':["102","1","9/10"], 'Cn':False,
    'Nh':["112","1","1/6"], 'Fl':["112","1","2/6"], 'Mc':["112","1","3/6"], 'Lv':["112","1","4/6"], 'Ts':["112","1","5/6"], 'Og':False}

def generate_dirac_input(molecule,
                        basis,
                        special_basis,
                        symmetry,
                        hamiltonian,
                        spinfree,
                        point_nucleus,
                        uncontract,
                        speed_of_light,
                        amfich,
                        maxiter,
                        get_mrconee,
                        openshell):
    """This function creates and saves a Dirac input file.

    Args:
        molecule: An instance of the MolData class.
        symmetry: Boolean to specify the use of symmetry
        basis: A string giving the basis set. An example is 'cc-pvtz'.
                Only optional if loading from file.
        special_basis: A list of two strings giving the default 
                and special basis set. An example is ["STO-3G","H cc-PVDZ"] 
                to specify that the hydrogens are in cc-pVDZ basis.
                Only optional if loading from file.
        hamiltonian: string defining the Hamiltonian (NONREL, X2C or X2CMMF)
        spinfree: Get rid of spin-coupling for 4- and 2-component relativistic hamiltonians
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
        uncontract: use uncontracted basis set.
        speed_of_light: Real value for the speed of light (137 a.u.) in atomic unit,
                        to be changed if wanted in order to increase or decrease relativistic
                        effects.
        get_mrconee: Extract the MRCONEE unformatted file containing the one-electron integrals.
        openshell: perform an openshell calculation

    Returns:
        input_file: A string giving the name of the saved input file, and the xyz file.
    """
    # Create Dirac geometry string.
    molecule_string = create_geometry_string(molecule.geometry)
    xyz_file = molecule.filename + '.xyz'
    with open(xyz_file, 'w') as f:
     f.write(str(molecule.n_atoms)+'\n')
     f.write(molecule.filename + ' # Auto-generated by ROSE\n')
     f.write(molecule_string)

    if basis == "special" and special_basis is None:
       raise SpecialBasisError('special_basis should be specified')
    elif basis == "special" and len(special_basis) != 2:
       raise SpecialBasisError('special_basis should be a list of two strings corresponding to the default and special basis, respectively')
    elif (basis != "special") and (special_basis is not None):
       raise SpecialBasisError('special_basis was specified without setting basis = "special"')

    if speed_of_light is not False and hamiltonian=="NONREL":
       raise SpeedOfLightError('A given speed of light has been specified without setting a relativistic Hamiltonian')

    if spinfree and hamiltonian=="NONREL":
       raise SpinFreeError('spinfree keyword should be set to False for nonrelativistic Hamiltonian. This keyword is only used for 4- and 2-component relativistic Hamiltonians.')
       
    # Write input file and return handle.
    input_file = molecule.filename + '.inp'
    with open(input_file, 'w') as f:
      f.write("**DIRAC\n")
      f.write(".4INDEX\n")
      f.write(".WAVE FUNCTION\n")
      f.write("**WAVE FUNCTION\n")
      f.write(".SCF\n")
      f.write("*SCF\n")
      f.write(".MAXITR\n")
      f.write(str(int(maxiter))+"\n")
      if (openshell is not False):
         f.write(".CLOSED SHELL\n"+openshell[0]+"\n")
         f.write(".OPEN SHELL\n"+openshell[1]+"\n"+openshell[2]+"\n")
      f.write("**HAMILTONIAN\n")
      f.write("."+ hamiltonian + "\n")
      if spinfree:
       f.write(".SPINFREE\n")
      if amfich is not False:
       f.write("*AMFI\n")
       f.write(".AMFICH\n")
       f.write(str(amfich) + "\n")
      f.write("**INTEGRALS\n")
      if point_nucleus:
       f.write(".NUCMOD\n")
       f.write(" 1\n")
      if uncontract:
       f.write("*READIN\n")
       f.write(".UNCONTRACT\n")
      f.write("**GENERAL\n")
      f.write(".ACMOUT\n")
      if speed_of_light is not False:
       f.write(".CVALUE\n")
       f.write(" " + str(speed_of_light) + "\n")
      f.write("**MOLTRA\n")
      f.write(".NO4IND\n")
      if get_mrconee:
      # all active has to be set for the MRCONEE and CHECKPOINT files to have the same number of spinors.
      # This is important for the test only, otherwise we don't need the MRCONEE file.
       f.write(".ACTIVE\n") 
       f.write(" all\n")
      f.write("**MOLECULE\n")
      f.write("*CHARGE\n")
      f.write(".CHARGE\n")
      f.write(" " + str(molecule.charge) + "\n")
      if not symmetry:
       f.write("*SYMMETRY\n")
       f.write(".NOSYM\n")
      f.write("*BASIS\n")
      f.write(".DEFAULT\n")
      if basis == "special":
       f.write(special_basis[0] + "\n")
       f.write(".SPECIAL\n")
       f.write(special_basis[1] + "\n")
      else:
       f.write(basis + "\n")
      f.write("*END OF INPUT\n")

    return input_file, xyz_file

def run_dirac_mocoef(molecule,
             basis=None,
             special_basis=None,
             symmetry=False,
             hamiltonian="NONREL",
             spinfree=False,
             point_nucleus=False,
             uncontract=True,
             speed_of_light=False,
             amfich=False,
             maxiter=50,
             get_mrconee=False,
             get_checkpoint=True,
             restart=False,
             memory=False,
             openshell=False):
    """This function runs a Dirac calculation.

    Args:
        molecule: An instance of the MolData class.
        symmetry: Optional boolean to remove symmetry in the calculation
        hamiltonian: string defining the Hamiltonian (NONREL, X2C or X2CMMF)
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
        speed_of_light: Optional real to give another value to the speed of light

    Returns:
        molecule: The updated MolData_Dirac object.
    """

    if basis is None: raise ValueError("basis must be specified.")

    add_to_name = "_" + basis + "_" + hamiltonian
    if spinfree: add_to_name += "_spinfree"
    if speed_of_light is not False: add_to_name += "_c"+str(speed_of_light)
    if not symmetry: add_to_name += "_nosym"
    
    molecule.filename = molecule.filename + add_to_name
    molecule.name = molecule.name + add_to_name

    # Prepare input.
    input_file, xyz_file = generate_dirac_input(molecule,
                        basis,
                        special_basis,
                        symmetry,
                        hamiltonian,
                        spinfree,
                        point_nucleus,
                        uncontract,
                        speed_of_light,
                        amfich,
                        maxiter,
                        get_mrconee,
                        openshell)

    # Run Dirac
    print('Starting Dirac calculation\n')
    get=''
    if get_mrconee: get = " --get MRCONEE"
    mb = ''
    if memory is not False:
     mb = "--mb="+str(memory)
    if not restart:
     subprocess.check_call("pam " + mb + " --mol=" + xyz_file + " --inp=" + input_file + get + " --outcmo --silent --noarch", shell=True, cwd=molecule.data_directory)
    else:
     subprocess.check_call("pam " + mb + " --mol=" + xyz_file + " --inp=" + input_file + " --incmo --silent --noarch", shell=True, cwd=molecule.data_directory)

    return molecule
