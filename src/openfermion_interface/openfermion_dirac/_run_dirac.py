# OpenFermion plugin to interface with Psi4.
# Copyright 2017 The OpenFermion Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Functions to prepare dirac input and run calculations."""
from __future__ import absolute_import

import os
import re
import subprocess
import warnings


def create_geometry_string(geometry):
    """This function converts MolecularData geometry to psi4 geometry.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances in
            angstrom. Use atomic symbols to specify atoms.

    Returns:
        geo_string: A string giving the geometry for each atom on a line, e.g.:
            H 0. 0. 0.
            H 0. 0. 0.7414
    """
    geo_string = ''
    for item in geometry:
        atom = item[0]
        coordinates = item[1]
        line = '{} {} {} {}'.format(atom,
                                    coordinates[0],
                                    coordinates[1],
                                    coordinates[2])
        if len(geo_string) > 0:
            geo_string += '\n'
        geo_string += line
    return geo_string


def generate_dirac_input(molecule,
                        run_scf,
                        run_mp2,
                        run_cisd,
                        run_ccsd,
                        run_fci,
                        relativistic,
                        verbose,
                        tolerate_error,
                        memory,
                        template_file):
    """This function creates and saves a psi4 input file.

    Args:
        molecule: An instance of the MolecularData class.
        run_scf: Boolean to run SCF calculation.
        run_mp2: Boolean to run MP2 calculation.
        run_cisd: Boolean to run CISD calculation.
        run_ccsd: Boolean to run CCSD calculation.
        run_fci: Boolean to FCI calculation.
        verbose: Boolean whether to print calculation results to screen.
        tolerate_error: Whether to fail or merely warn when Psi4 fails.
        memory: Int giving amount of memory to allocate in MB.
        template_file(str): Specify the filename of a dirac template

    Returns:
        input_file: A string giving the name of the saved input file, and the xyz file.
    """
    # Create Dirac geometry string.
    geo_string = create_geometry_string(molecule.geometry)
    xyz_file = molecule.filename + '.xyz'
    with open(xyz_file, 'w') as f:
     f.write(str(molecule.n_atoms)+'\n')
     f.write(molecule.filename + ' # anything can be in this line\n')
     f.write(geo_string)
    f.close()

    # Find the dirac_directory.
    dirac_directory = os.path.dirname(os.path.realpath(__file__))
    # Write input file and return handle.
    input_file = molecule.filename + '.inp'
    with open(input_file, 'w') as f:
      f.write("**DIRAC\n")
      f.write(".4INDEX\n")
      f.write(".WAVE FUNCTION\n")
      f.write("**WAVE FUNCTION\n")
      f.write(".SCF\n")
      if run_ccsd:
       f.write(".RELCCSD\n")
      if relativistic is False:
       f.write("**HAMILTONIAN\n")
       f.write(".NONREL\n")
       f.write("**INTEGRALS\n")
       f.write(".NUCMOD\n")
       f.write(" 1\n")
      f.write("**MOLTRA\n")
      f.write(".MDCINT\n")
      f.write("**MOLECULE\n")
      f.write("*BASIS\n")
      f.write(".DEFAULT\n")
      if molecule.basis == "explicit" or molecule.basis == "special":
       f.write(" sto-3g\n")
       f.write(molecule.special_basis+"\n")
      else:
       f.write(" " + molecule.basis + "\n")
      f.write("*END OF INPUT\n")
    f.close()

    return input_file, xyz_file

def clean_up(molecule, delete_input=True, delete_xyz=True, delete_output=False, delete_MRCONEE=True,
             delete_MDCINT=True, delete_FCIDUMP=False):
    os.remove("FCITABLE")
    input_file = molecule.filename + '.inp'
    xyz_file = molecule.filename + '.xyz'
    output_file_dirac = molecule.filename + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    os.rename("FCIDUMP", "FCIDUMP_" + molecule.name)
    os.rename(output_file_dirac,output_file)
    run_directory = os.getcwd()
    for local_file in os.listdir(run_directory):
        if local_file.endswith('.clean'):
            os.remove(run_directory + '/' + local_file)
    try:
        os.remove('timer.dat')
    except:
        pass
    if delete_input:
        os.remove(input_file)
    if delete_output:
        os.remove(output_file)
    if delete_xyz:
        os.remove(xyz_file)
    if delete_MRCONEE:
        os.remove("MRCONEE")
    if delete_MDCINT:
        os.remove("MDCINT")
    if delete_FCIDUMP:
        os.remove("FCIDUMP" + molecule.name)


def run_dirac(molecule,
             run_scf=True,
             run_mp2=False,
             run_cisd=False,
             run_ccsd=False,
             run_fci=False,
             relativistic=True,
             verbose=False,
             tolerate_error=False,
             delete_input=False,
             delete_output=False,
             delete_xyz=False,
             delete_MRCONEE=False,
             delete_MDCINT=False,
             delete_FCIDUMP=False,
             memory=8000,
             template_file=None):
    """This function runs a dirac calculation.

    Args:
        molecule: An instance of the MolecularData class.
        run_scf: Optional boolean to run SCF calculation.
        run_mp2: Optional boolean to run MP2 calculation.
        run_cisd: Optional boolean to run CISD calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        run_fci: Optional boolean to FCI calculation.
        verbose: Boolean whether to print calculation results to screen.
        tolerate_error: Optional boolean to warn or raise when Psi4 fails.
        delete_input: Optional boolean to delete psi4 input file.
        delete_output: Optional boolean to delete psi4 output file.
        memory: Optional int giving amount of memory to allocate in MB.
        template_file(str): Path to Psi4 template file

    Returns:
        molecule: The updated MolecularData object.

    Raises:
        psi4 errors: An error from psi4.
    """
    # Prepare input.
    input_file, xyz_file = generate_dirac_input(molecule,
                        run_scf,
                        run_mp2,
                        run_cisd,
                        run_ccsd,
                        run_fci,
                        relativistic,
                        verbose,
                        tolerate_error,
                        memory,
                        template_file)

    # Run Dirac
    try:
        print('Starting Dirac calculation\n')
        subprocess.call('pam --mol=' + xyz_file + ' --inp=' + input_file + ' --get="MRCONEE MDCINT" --silent --noarch > output_script', shell=True)
    except:
        print('Dirac calculation for {} has failed. Most certainely a problem with pam.'.format(molecule.name))

#    try:
#        subprocess.call('grep "exit           : normal" pam_file | wc -l')

    # Compute the FCIDUMP
    try:
        print('Creation of the FCIDUMP file\n')
        subprocess.call('dirac_openfermion_mointegral_export.x >> output_script', shell=True)
    except:
        print('Creation of FCIDUMP has failed.')
    else:
        clean_up(molecule, delete_input, delete_xyz, delete_output, delete_MRCONEE, delete_MDCINT, delete_FCIDUMP)

    # Return updated molecule instance.
    #try:
    #    molecule.load()
    #except:
    #    warnings.warn('No calculation saved. '
    #                  'Dirac segmentation fault possible.',
    #                  Warning)
    return molecule
