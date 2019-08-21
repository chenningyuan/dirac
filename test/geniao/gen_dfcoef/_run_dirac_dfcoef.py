"""Functions to prepare Dirac input and run calculations."""
from __future__ import absolute_import

import os
import re
import subprocess
import warnings

class SpecialBasisError(Exception):
    pass
class SpeedOfLightError(Exception):
    pass

def create_geometry_string(geometry):
    """This function converts MolData_Dirac geometry to Dirac geometry.
       This function is taken from Openfermion-Psi4 interface.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances in
            angstrom. Use atomic symbols to specify atoms.

    Returns:
        molecule_string: A string giving the geometry for each atom on a line, e.g.:
            H 0. 0. 0.
            H 0. 0. 0.7414
    """
    molecule_string = ''
    for item in geometry:
        atom = item[0]
        coordinates = item[1]
        line = '{} {} {} {}'.format(atom,
                                    coordinates[0],
                                    coordinates[1],
                                    coordinates[2])
        if len(molecule_string) > 0:
            molecule_string += '\n'
        molecule_string += line
    return molecule_string

def generate_dirac_input(molecule,
                        symmetry,
                        hamiltonian,
                        point_nucleus,
                        uncontract,
                        speed_of_light,
                        version,
                        get_mrconee,
                        test,
                        exponent):
    """This function creates and saves a Dirac input file.

    Args:
        molecule: An instance of the MolData_Dirac class.
        symmetry: Boolean to specify the use of symmetry
        hamiltonian: string defining the Hamiltonian (NONREL, X2C or X2CMMF)
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
        speed_of_light: Real value for the speed of light (137 a.u.) in atomic unit,
                        to be changed if wanted in order to increase or decrease relativistic
                        effects.

    Returns:
        input_file: A string giving the name of the saved input file, and the xyz file.
    """
    # Create Dirac geometry string.
    molecule_string = create_geometry_string(molecule.geometry)
    xyz_file = molecule.filename + '.xyz'
    with open(xyz_file, 'w') as f:
     f.write(str(molecule.n_atoms)+'\n')
     f.write(molecule.filename + ' # anything can be in this line\n')
     f.write(molecule_string)

    if molecule.basis == "special" and molecule.special_basis is None:
       raise SpecialBasisError('special_basis should be specified')
    elif molecule.basis == "special" and len(molecule.special_basis) != 2:
       raise SpecialBasisError('special_basis should be a list of two strings corresponding to the default and special basis, respectively')
    elif (molecule.basis != "special") and (molecule.special_basis is not None):
       raise SpecialBasisError('special_basis was specified without setting basis = "special"')

    if speed_of_light is not False and hamiltonian=="NONREL":
       raise SpeedOfLightError('A given speed of light has been specified without setting a relativistic Hamiltonian')

    if test and not get_mrconee:
        warnings.warn('MRCONEE is needed for the test. If you have not provided it, please switch on get_mrconee=True')
       

    nelec=int(molecule.n_electrons)
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
      f.write("500\n")
      f.write("**HAMILTONIAN\n")
      f.write("."+ hamiltonian + "\n")
      f.write("**INTEGRALS\n")
      if point_nucleus:
       f.write(".NUCMOD\n")
       f.write(" 1\n")
      if uncontract:
       f.write("*READIN\n")
       f.write(".UNCONTRACT\n")
      if speed_of_light is not False:
       f.write("**GENERAL\n")
       f.write(".CVALUE\n")
       f.write(" " + str(speed_of_light) + "\n")
      f.write("**MOLTRA\n")
      f.write(".NO4IND\n")
      if get_mrconee:
      # all active has to be set for the MRCONEE and DFCOEF files to have the same number of spinors.
      # This is important for the test only, otherwise we don't need the MRCONEE file.
       f.write(".MDCINT\n")
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
      if molecule.basis == "special":
       f.write(molecule.special_basis[0] + "\n")
       f.write(".SPECIAL\n")
       f.write(molecule.special_basis[1] + "\n")
      else:
       f.write(molecule.basis + "\n")
      f.write("*END OF INPUT\n")

    with open(molecule.data_directory+"/INPUT_IAO","w") as f:
      f.write(version+"\n")
      f.write(str(molecule.charge)+"\n")
      f.write(str(exponent)+"\n")
      if test==True:
        f.write("1")
      else:
        f.write("0")
    
    return input_file, xyz_file

def rename(molecule):
    output_file_dirac = molecule.name + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    os.rename(output_file_dirac,output_file)

def run_dirac_dfcoef(molecule,
             symmetry=False,
             hamiltonian="NONREL",
             point_nucleus=False,
             uncontract=True,
             speed_of_light=False,
             version="Stndrd_2013",
             exponent=4,
             get_mrconee=False,
             test=False,
             restart=False):
    """This function runs a Dirac calculation.

    Args:
        molecule: An instance of the MolData_Dirac class.
        symmetry: Optional boolean to remove symmetry in the calculation
        hamiltonian: string defining the Hamiltonian (NONREL, X2C or X2CMMF)
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
        speed_of_light: Optional real to give another value to the speed of light

    Returns:
        molecule: The updated MolData_Dirac object.
    """

    # Prepare input.
    input_file, xyz_file = generate_dirac_input(molecule,
                        symmetry,
                        hamiltonian,
                        point_nucleus,
                        uncontract,
                        speed_of_light,
                        version,
                        get_mrconee,
                        test,
                        exponent)

    # Run Dirac
    print('Starting Dirac calculation\n')
    if not restart and not get_mrconee:
     subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --get='DFCOEF' --silent --noarch", shell=True, cwd=molecule.data_directory)
    elif not restart and get_mrconee:
     subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --get='DFCOEF MRCONEE' --silent --noarch", shell=True, cwd=molecule.data_directory)
    elif restart:
     subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --put='DFCOEF' --silent --noarch", shell=True, cwd=molecule.data_directory)

    try:
     rename(molecule)
    except:
     print("Rename didn't work")

    return molecule
