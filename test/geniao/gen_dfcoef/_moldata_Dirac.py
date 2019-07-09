import numpy
import os
import re
import uuid

# Define a compatible basestring for checking between Python 2 and 3
try:
    basestring
except NameError:  # pragma: no cover
    basestring = str

periodic_table = [
    'foo',
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
                  charge,
                  hamiltonian,
                  symmetry,
                  speed_of_light,
                  description):
    """Function to name molecules. Inspired from _molecule_data.py of OpenFermion.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
            Distances in angstrom. Use atomic symbols to specify atoms.
        basis: A string giving the basis set. An example is 'cc-pvtz'.
        charge: An integer giving the total molecular charge.
        hamiltonian: string defining the Hamiltonian (NONREL, X2C or X2CMMF)
        symmetry: A boolean to use symmetry or not
        speed_of_light: Real number to specify the speed of light manually.
        description: A string giving a description. As an example,
            for dimers a likely description is the bond length (e.g. 0.7414).

    Returns:
        name: A string giving the name of the instance.
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

    # Add charge.
    if charge > 0:
        name += '_{}+'.format(charge)
    elif charge < 0:
        name += '_{}-'.format(charge)

    # Optionally add descriptive tag and return.
    if description:
        name += '_{}'.format(description)
    name += '_{}'.format(hamiltonian)
    if symmetry is False:
        name += '_nosym'
    if speed_of_light is not False:
        name += '_c'+str(speed_of_light)
    return name

def geometry_from_file(file_name):
    """Function to create molecular geometry from text file. This function is the same as in _molecule_data.py of OpenFermion.

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

class MolData_Dirac(object):

    """Attributes:
        geometry: A list of tuples giving the coordinates of each atom. An
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances
            in angstrom. Use atomic symbols to specify atoms.
        basis: A string giving the basis set. An example is 'cc-pvtz'.
        charge: An integer giving the total molecular charge. Defaults to 0.
        description: An optional string giving a description. As an example,
            for dimers a likely description is the bond length (e.g. 0.7414).
        name: A string giving a characteristic name for the instance.
        filename: The name of the file where the molecule data is saved.
        n_atoms: Integer giving the number of atoms in the molecule.
        n_electrons: Integer giving the number of electrons in the molecule.
        atoms: List of the atoms in molecule sorted by atomic number.
        protons: List of atomic charges in molecule sorted by atomic number.
        scf_energy: Energy from open or closed shell Hartree-Fock.
        scf_iterations: Number of iterations of the SCF procedure.
    """
    def __init__(self, geometry=None, basis=None, special_basis=None,
                 charge=0, description="", filename="", data_directory=None, hamiltonian="NONREL",
                 symmetry=True, speed_of_light=False):
        """Initialize molecular metadata which defines class.

        Args:
            geometry: A list of tuples giving the coordinates of each atom.
                An example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
                Distances in angstrom. Use atomic symbols to
                specify atoms. Only optional if loading from file.
            basis: A string giving the basis set. An example is 'cc-pVTZ'.
                Only optional if loading from file.
            special_basis: A list of two strings giving the default 
                and special basis set. An example is ["STO-3G","H cc-PVDZ"] 
                to specify that the hydrogens are in cc-pVDZ basis.
                Only optional if loading from file.
            charge: An integer giving the total molecular charge. Defaults
                to 0.  Only optional if loading from file.
            description: A optional string giving a description. As an
                example, for dimers a likely description is the bond length
                (e.g. 0.7414).
            filename: An optional string giving name of file.
                If filename is not provided, one is generated automatically.
            data_directory: Optional data directory to change from default
                data directory specified in config file.
            hamiltonian: string defining the Hamiltonian (NONREL, X2C or X2CMMF)
            speed_of_light: real number of specify the speed of light manually.
            symmetry: boolean to specify the use of symmetry or not
        """
        # Check appropriate data as been provided and autoload if requested.
        if ((geometry is None) or
                (basis is None)):
            if filename:
                if filename[-5:] == '.hdf5':
                    self.filename = filename[:(len(filename) - 5)]
                else:
                    self.filename = filename
                self.load()
                self.init_lazy_properties()
                return
            else:
                raise ValueError("Geometry, basis must be specified.")

        # Metadata fields which must be provided.
        self.geometry = geometry
        self.basis = basis
        self.data_directory = data_directory

        # Metadata fields with default values.
        self.charge = charge
        if (not isinstance(description, basestring)):
            raise TypeError("description must be a string.")
        self.description = description
        self.hamiltonian = hamiltonian
        self.symmetry = symmetry
        self.speed_of_light = speed_of_light
        self.special_basis = special_basis

        # Name molecule and get associated filename
        self.name = name_molecule(geometry, basis, charge, hamiltonian, symmetry, speed_of_light, description)
        if filename:
            if filename[-5:] == '.hdf5':
                filename = filename[:(len(filename) - 5)]
            self.filename = filename
        else:
            if data_directory is None:
                self.data_directory = os.getcwd()
            self.filename = self.data_directory + '/' + self.name

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

    def get_scf(self):
        self.scf_energy = None
        self.scf_iterations = None
        if os.path.exists(self.filename + '.out'):
           with open(self.filename + '.out', "r") as f:
             for line in f:
                if re.search("Total energy                             :", line):
                  self.scf_energy=line.rsplit(None, 1)[-1]
                elif re.search("Convergence after", line):
                  self.scf_iterations=line.rsplit()[-2]
        else:
           raise FileNotFoundError('output not found, check your run_dirac calculation')
        return self.scf_energy, self.scf_iterations
