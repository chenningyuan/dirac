#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""Class and functions to store quantum chemistry data from a Dirac calculation"""

#import h5py
import numpy
import os
import re
#import uuid

from openfermion.config import *
from openfermion.ops import InteractionOperator, InteractionRDM


"""NOTE ON PQRS CONVENTION:
  The data structures which hold fermionic operators / integrals /
  coefficients assume a particular convention which depends on how integrals
  are labeled:
  h[p,q]=\int \phi_p(x)* (T + V_{ext}) \phi_q(x) dx
  h[p,q,r,s]=\int \phi_p(x)* \phi_q(y)* V_{elec-elec} \phi_r(y) \phi_s(x) dxdy
  With this convention, the molecular Hamiltonian becomes
  H =\sum_{p,q} h[p,q] a_p^\dagger a_q
    + 0.5 * \sum_{p,q,r,s} h[p,q,r,s] a_p^\dagger a_q^\dagger a_r a_s
"""

# Define a compatible basestring for checking between Python 2 and 3
try:
    basestring
except NameError:  # pragma: no cover
    basestring = str


# Define error objects which inherit from Exception.
class MoleculeNameError(Exception):
    pass

class MissingCalculationError(Exception):
    pass

# Functions to change from Bohr to angstroms and back.
def bohr_to_angstroms(distance):
    # Value defined so it is the inverse to numerical precision of angs to bohr
    return 0.5291772458017723 * distance

def angstroms_to_bohr(distance):
    return 1.889726 * distance

# The Periodic Table as a python list and dictionary.
periodic_table = [
    '?',
    'H', 'He',
    'Li', 'Be',
    'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
    'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
    'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
periodic_hash_table = {}
for atomic_number, atom in enumerate(periodic_table):
    periodic_hash_table[atom] = atomic_number

# Spin polarization of atoms on period table.
periodic_polarization = [-1,
                         1, 0,
                         1, 0, 1, 2, 3, 2, 1, 0,
                         1, 0, 1, 2, 3, 2, 1, 0,
                         1, 0, 1, 2, 3, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0,
                         1, 0, 1, 2, 5, 6, 5, 8, 9, 0, 1, 0, 1, 2, 3, 2, 1, 0]

def name_molecule(geometry,
                  basis,
                  multiplicity,
                  charge,
                  relativistic,
                  symmetry,
                  speed_of_light,
                  description):
    """Function to name molecules.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
            Distances in angstrom. Use atomic symbols to specify atoms.
        basis: A string giving the basis set. An example is 'cc-pvtz'.
        multiplicity: An integer giving the spin multiplicity.
        charge: An integer giving the total molecular charge.
        relativistic: A boolean which tells you if the calculation is
                      relativistic or not.
        symmetry: A boolean to use symmetry or not
        speed_of_light: Real number to specify the speed of light manually.
        description: A string giving a description. As an example,
            for dimers a likely description is the bond length (e.g. 0.7414).

    Returns:
        name: A string giving the name of the instance.

    Raises:
        MoleculeNameError: If spin multiplicity is not valid.
    """
    if not isinstance(geometry, basestring):
        # Get sorted atom vector.
        atoms = [item[0] for item in geometry]
        atom_charge_info = [(atom, atoms.count(atom)) for atom in set(atoms)]
        sorted_info = sorted(atom_charge_info,
                             key=lambda atom: periodic_hash_table[atom[0]])

        # Name molecule.
        name = '{}{}'.format(sorted_info[0][0], sorted_info[0][1])
        for info in sorted_info[1::]:
            name += '-{}{}'.format(info[0], info[1])
    else:
        name = geometry

    # Add basis.
    name += '_{}'.format(basis)

    # Add multiplicity.
    multiplicity_dict = {1: 'singlet',
                         2: 'doublet',
                         3: 'triplet',
                         4: 'quartet',
                         5: 'quintet',
                         6: 'sextet',
                         7: 'septet',
                         8: 'octet',
                         9: 'nonet',
                         10: 'dectet',
                         11: 'undectet',
                         12: 'duodectet'}
    if (multiplicity not in multiplicity_dict):
        raise MoleculeNameError('Invalid spin multiplicity provided.')
    else:
        name += '_{}'.format(multiplicity_dict[multiplicity])

    # Add charge.
    if charge > 0:
        name += '_{}+'.format(charge)
    elif charge < 0:
        name += '_{}-'.format(charge)

    # Optionally add descriptive tag and return.
    if description:
        name += '_{}'.format(description)
    if relativistic:
        name += '_rel'
    if symmetry is False:
        name += '_nosym'
    if speed_of_light is not False:
        name += '_c'+str(speed_of_light)
    return name


def geometry_from_file(file_name):
    """Function to create molecular geometry from text file.

    Args:
        file_name: a string giving the location of the geometry file.
            It is assumed that geometry is given for each atom on line, e.g.:
            H 0. 0. 0.
            H 0. 0. 0.7414

    Returns:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
            Distances in angstrom. Use atomic symbols to specify atoms.
    """
    geometry = []
    with open(file_name, 'r') as stream:
        for line in stream:
            data = line.split()
            if len(data) == 4:
                atom = data[0]
                coordinates = (float(data[1]), float(data[2]), float(data[3]))
                geometry += [(atom, coordinates)]
    return geometry


class MolecularData_Dirac(object):

    """Class for storing molecule data from a fixed basis set at a fixed
    geometry that is obtained from classical electronic structure
    packages. Not every field is filled in every calculation. All data
    that can (for some instance) exceed 10 MB should be saved
    separately. Data saved in HDF5 format.

    Attributes:
        geometry: A list of tuples giving the coordinates of each atom. An
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances
            in angstrom. Use atomic symbols to specify atoms.
        basis: A string giving the basis set. An example is 'cc-pvtz'.
        charge: An integer giving the total molecular charge. Defaults to 0.
        multiplicity: An integer giving the spin multiplicity.
        description: An optional string giving a description. As an example,
            for dimers a likely description is the bond length (e.g. 0.7414).
        name: A string giving a characteristic name for the instance.
        filename: The name of the file where the molecule data is saved.
        n_atoms: Integer giving the number of atoms in the molecule.
        n_electrons: Integer giving the number of electrons in the molecule.
        atoms: List of the atoms in molecule sorted by atomic number.
        protons: List of atomic charges in molecule sorted by atomic number.
        hf_energy: Energy from open or closed shell Hartree-Fock.
        nuclear_repulsion: Energy from nuclei-nuclei interaction.
        canonical_orbitals: numpy array giving canonical orbital coefficients.
        n_orbitals: Integer giving total number of spatial orbitals.
        n_qubits: Integer giving total number of qubits that would be needed.
        orbital_energies: Numpy array giving the canonical orbital energies.
        fock_matrix: Numpy array giving the Fock matrix.
        overlap_integrals: Numpy array of AO overlap integrals
        one_body_integrals: Numpy array of one-electron integrals
        two_body_integrals: Numpy array of two-electron integrals
        mp2_energy: Energy from MP2 perturbation theory.
        cisd_energy: Energy from configuration interaction singles + doubles.
        cisd_one_rdm: Numpy array giving 1-RDM from CISD calculation.
        cisd_two_rdm: Numpy array giving 2-RDM from CISD calculation.
        fci_energy: Exact energy of molecule within given basis.
        fci_one_rdm: Numpy array giving 1-RDM from FCI calculation.
        fci_two_rdm: Numpy array giving 2-RDM from FCI calculation.
        ccsd_energy: Energy from coupled cluster singles + doubles.
        ccsd_single_amps: Numpy array holding single amplitudes
        ccsd_double_amps: Numpy array holding double amplitudes
        general_calculations: A dictionary storing general calculation results
            for this system annotated by the key.
    """
    def __init__(self, geometry=None, basis=None, special_basis=None, multiplicity=None,
                 charge=0, description="", filename="", data_directory=None, relativistic=False,
                 symmetry=True, speed_of_light=False):
        """Initialize molecular metadata which defines class.

        Args:
            geometry: A list of tuples giving the coordinates of each atom.
                An example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
                Distances in angstrom. Use atomic symbols to
                specify atoms. Only optional if loading from file.
            basis: A string giving the basis set. An example is 'cc-pvtz'.
                Only optional if loading from file.
            special_basis: A string giving a special basis set.
                Only optional if loading from file.
            charge: An integer giving the total molecular charge. Defaults
                to 0.  Only optional if loading from file.
            multiplicity: An integer giving the spin multiplicity.  Only
                optional if loading from file.
            description: A optional string giving a description. As an
                example, for dimers a likely description is the bond length
                (e.g. 0.7414).
            filename: An optional string giving name of file.
                If filename is not provided, one is generated automatically.
            data_directory: Optional data directory to change from default
                data directory specified in config file.
            relativistic: A boolean which tells you if the calculation is relativistic
                or not.
            speed_of_light: real number of specify the speed of light manually.
            symmetry: boolean to specify the use of symmetry or not
        """
        # Check appropriate data as been provided and autoload if requested.
        if ((geometry is None) or
                (basis is None) or
                (multiplicity is None)):
            if filename:
                if filename[-5:] == '.hdf5':
                    self.filename = filename[:(len(filename) - 5)]
                else:
                    self.filename = filename
                self.load()
                self.init_lazy_properties()
                return
            else:
                raise ValueError("Geometry, basis, multiplicity must be"
                                 "specified when not loading from file.")

        # Metadata fields which must be provided.
        self.geometry = geometry
        self.basis = basis
        self.multiplicity = multiplicity

        # Metadata fields with default values.
        self.charge = charge
        if (not isinstance(description, basestring)):
            raise TypeError("description must be a string.")
        self.description = description
        self.relativistic = relativistic
        self.symmetry = symmetry
        self.speed_of_light = speed_of_light
        self.special_basis = special_basis

        # Name molecule and get associated filename
        self.name = name_molecule(geometry, basis, multiplicity,
                                  charge, relativistic, symmetry, speed_of_light, description)
        if filename:
            if filename[-5:] == '.hdf5':
                filename = filename[:(len(filename) - 5)]
            self.filename = filename
        else:
            if data_directory is None:
                self.filename = DATA_DIRECTORY + '/' + self.name
            else:
                self.filename = data_directory + '/' + self.name

        # Attributes generated automatically by class.
        if not isinstance(geometry, basestring):
            self.n_atoms = len(geometry)
            self.atoms = sorted([row[0] for row in geometry],
                                key=lambda atom: periodic_hash_table[atom])
            self.protons = [periodic_hash_table[atom] for atom in self.atoms]
            self.n_electrons = sum(self.protons) - charge
        else:
            self.n_atoms = 0
            self.atoms = []
            self.protons = 0
            self.n_electrons = 0

    def get_integrals_FCIDUMP(self):
        if os.path.exists("FCIDUMP_" + self.name) == True:
             E_core = 0
             spinor = {}
             one_body_int = {}
             two_body_int = {}
             num_lines = sum(1 for line in open('FCIDUMP_' + self.name))
             with open("FCIDUMP_" + self.name) as f:
               start_reading=0
               for line in f:
                 start_reading+=1
                 if "&END" in line:
                   break
               listed_values = [[token for token in line.split()] for line in f.readlines()] 
               for row in range(num_lines-start_reading):
                 a_1 = int(listed_values[row][1])
                 a_2 = int(listed_values[row][2])
                 a_3 = int(listed_values[row][3])
                 a_4 = int(listed_values[row][4])
                 if a_4 == 0 and a_3 == 0:
                   if a_2 == 0:
                     if a_1 == 0:
                       E_core = float(listed_values[row][0])
                     else:
                       spinor[a_1] = float(listed_values[row][0])
                   else:
                     one_body_int[a_1,a_2] = float(listed_values[row][0])
                 else:
                   two_body_int[a_1,a_2,a_3,a_4] = float(listed_values[row][0])
             f.close()
        else:
             raise FileNotFoundError('FCIDUMP not found, first make a run_dirac calculation')
        return E_core, spinor, one_body_int, two_body_int

    def get_energies(self):
        hf_energy = None
        mp2_energy = None
        ccsd_energy = None
        if os.path.exists(self.name + '.out') == True:
           with open(self.name + '.out', "r") as f:
             for line in f:
                if re.search("Total energy                             :", line):
                  hf_energy=line.rsplit(None, 1)[-1]
                if re.search("@ Total MP2 energy", line):
                  mp2_energy=line.rsplit(None, 1)[-1]
                if re.search("@ Total CCSD energy", line):
                  ccsd_energy=line.rsplit(None, 1)[-1]
        else:
           raise FileNotFoundError('output not found, check your run_dirac calculation')
        return hf_energy, mp2_energy, ccsd_energy

    def get_molecular_hamiltonian(self, relativistic=False):
        """Output arrays of the second quantized Hamiltonian coefficients.

        Args:
            relativistic: Boolean to specify if the calculation is relativistic or not,
                          which is important to sort the integrals.

        Returns:
            molecular_hamiltonian: An instance of the MolecularOperator class.

        Note:
           OpenFermion requires all integrals, without accounting for permutation symmetry 
           or restricted cases.
           Hence, I modified dirac_mointegral_export.F90 into dirac_openfermion_mointegral_export.F90 
           to generate the FCIDUMP which write explicitly the
           integrals even if they are equivalent by permutation symmetry. 
           The FCIDUMP of a relativistic calculation contains every orbitals needed.
           For the non relativistic case, we have to add the symmetries due to the restricted
           calculation. To do so, I also had to add a subroutine in 
           dirac_openfermion_mointegral_export.F90 to index the orbitals as I wanted.
           Note that the integrals in Dirac and Openfermion are not sorted the same way. In Openfermion,
           they correspond to the h_{p,q,r,s} term, which is equal to < pq | sr > = ( ps | qr ),
           where <.|.> refers to the physicist's notation, and (.|.) to the chemist's one, used in Dirac.
           So p,q,r,s in Openfermion reads p,s,q,r in Dirac, or reversely,
              p,q,r,s in Dirac       reads p,r,s,q in Openfermion.
        """
        # Get active space integrals.
        E_core, spinor, one_body_integrals, two_body_integrals = self.get_integrals_FCIDUMP()
        n_qubits = len(one_body_integrals)
        # Initialize Hamiltonian coefficients.
        one_body_coefficients = numpy.zeros((n_qubits, n_qubits))
        two_body_coefficients = numpy.zeros((n_qubits, n_qubits,
                                             n_qubits, n_qubits))

        if relativistic:
          for p in range(n_qubits):
            for q in range(n_qubits):
                if (p+1,q+1) in one_body_integrals: 
                   one_body_coefficients[p,q] = one_body_integrals[p+1,q+1]
                for r in range(n_qubits):
                    for s in range(n_qubits):
                        if (p+1,q+1,r+1,s+1) in two_body_integrals and relativistic is True:
                           two_body_coefficients[p,r,s,q] = two_body_integrals[p+1,q+1,r+1,s+1] / 2.0
        else:
          for p in range(n_qubits):
            for q in range(n_qubits):
                if (p+1,q+1) in one_body_integrals:
                   one_body_coefficients[p,q] = one_body_integrals[p+1,q+1]
                   one_body_coefficients[q,p] = one_body_integrals[p+1,q+1]

          #permutation symmetry
          for p in range(n_qubits//2):
            for q in range(n_qubits//2):
                for r in range(n_qubits//2):
                    for s in range(n_qubits//2):
                        if (2*p+1,2*q+1,2*r+1,2*s+1) in two_body_integrals:
                           two_body_coefficients[2*p,2*r,2*s,2*q] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*q,2*r,2*s,2*p] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*p,2*s,2*r,2*q] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*q,2*s,2*r,2*p] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*r,2*p,2*q,2*s] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*s,2*p,2*q,2*r] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*r,2*q,2*p,2*s] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
                           two_body_coefficients[2*s,2*q,2*p,2*r] = two_body_integrals[2*p+1,2*q+1,2*r+1,2*s+1] / 2.0
        # restricted calculation
          for p in range(n_qubits//2):
            for q in range(n_qubits//2):
                for r in range(n_qubits//2):
                    for s in range(n_qubits//2):
                        two_body_coefficients[2*p+1,2*q,2*r,2*s+1] = two_body_coefficients[2*p,2*q,2*r,2*s]
                        two_body_coefficients[2*p,2*q+1,2*r+1,2*s] = two_body_coefficients[2*p,2*q,2*r,2*s]
                        two_body_coefficients[2*p+1,2*q+1,2*r+1,2*s+1] = two_body_coefficients[2*p,2*q,2*r,2*s]

        # Truncate.
        one_body_coefficients[
            numpy.absolute(one_body_coefficients) < EQ_TOLERANCE] = 0.
        two_body_coefficients[
            numpy.absolute(two_body_coefficients) < EQ_TOLERANCE] = 0.

        # Cast to InteractionOperator class and return.
        molecular_hamiltonian = InteractionOperator(
            E_core, one_body_coefficients, two_body_coefficients)

        return molecular_hamiltonian
