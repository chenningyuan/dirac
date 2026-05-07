"""Class and functions to store quantum chemistry data."""
from __future__ import absolute_import

import h5py
import numpy
import os
import uuid
from ._utils import periodic_table, periodic_hash_table

# Define a compatible basestring for checking between Python 2 and 3
try:
    basestring
except NameError:  # pragma: no cover
    basestring = str
class MoleculeNameError(Exception):
    pass

def name_molecule(geometry,
                  multiplicity,
                  charge,
                  description):
    """Function to name molecules.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
            Distances in angstrom. Use atomic symbols to specify atoms.
        multiplicity: An integer giving the spin multiplicity.
        charge: An integer giving the total molecular charge.
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

    # Add multiplicity.
    if (multiplicity != None):
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

class moldata(object):

    """Class for storing molecule data from a fixed basis set at a fixed geometry.

    Attributes:
        geometry: A list of tuples giving the coordinates of each atom. An
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances
            in angstrom. Use atomic symbols to specify atoms.
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
    """
    def __init__(self, geometry=None, multiplicity=None, charge=0, description="", data_directory=None):
        """Initialize molecular metadata which defines class.

        Args:
            geometry: A list of tuples giving the coordinates of each atom.
                An example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
                Distances in angstrom. Use atomic symbols to
                specify atoms. Only optional if loading from file.
            charge: An integer giving the total molecular charge. Defaults
                to 0.
            multiplicity: An integer giving the spin multiplicity.  Defaults
                to None (relativistic calculation).
            description: A optional string giving a description. As an
                example, for dimers a likely description is the bond length
                (e.g. 0.7414).
            data_directory: Optional data directory to change from default
                data directory specified in config file.
        """
        if (geometry is None): raise ValueError("Geometry must be specified.")

        # Metadata fields which must be provided.
        self.geometry = geometry
        self.multiplicity = multiplicity
        self.data_directory = data_directory

        # Metadata fields with default values.
        self.charge = charge
        if (not isinstance(description, basestring)):
            raise TypeError("description must be a string.")
        self.description = description

        # Name molecule and get associated filename
        self.name = name_molecule(geometry, multiplicity, charge, description)
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

    def get_multiplicity(self):
        """Return (guess of) multiplicity"""
        if self.multiplicity is None:
            return self.n_electrons%2 + 1
        else:
            return self.multiplicity

    def get_spin(self):
        """Return (guess of) spin"""
        if self.multiplicity is None:
            return self.n_electrons%2
        else:
            return self.multiplicity - 1

    def get_n_alpha_electrons(self):
        """Return number of alpha electrons."""
        if self.multiplicity is None:
            return 0
        else:
            return int((self.n_electrons + (self.multiplicity - 1)) // 2)

    def get_n_beta_electrons(self):
        """Return number of beta electrons."""
        if self.multiplicity is None:
            return 0
        else:
            return int((self.n_electrons - (self.multiplicity - 1)) // 2)

    def write_hdf5(self,file_name='AOMO.h5',section='/input/molecule/'):
        """
        Write molecule information to hdf5-type file
        """
        import h5py
        h5file = h5py.File(file_name, 'w')
        # Arrays are written as extendable datasets, single strings or scalars are fixed size
        h5file.create_dataset(section+'n_atoms',data=[self.n_atoms])
        # Write geometry in original (unsorted) form
        nuc_charge = [periodic_hash_table[row[0]] for row in self.geometry]
        xyz        = [float(coord) for row in self.geometry for coord in row[1][:]]
        h5file.create_dataset(section+'nuc_charge',data=nuc_charge)
        h5file.create_dataset(section+'geometry',data=xyz)
        h5file.close()


    # This function does not belong here and should get energy from CHECKPOINT, 
    # keeping for now to make tests pass
    def get_DIRAC_scf_energy(self):
        import re
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

