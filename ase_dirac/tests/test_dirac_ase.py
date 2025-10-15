"""Tests for the ASE-DIRAC interface"""
import os
import glob
import shutil
import subprocess
import pytest

from ase import Atoms
from ase.build import molecule
from ase.units import Ha
from ase_dirac import DIRAC


# Check first if the `pam` executable is available
if not shutil.which("pam"):
    pytest.skip("pam executable not found or not in your path. "
                "Skipping tests.",
                allow_module_level=True)


@pytest.fixture(scope="session", autouse=True)
def clean_up_h5_files():
    """Runs at the end of all tests."""
    yield
    # Define the pattern for files to delete
    file_pattern = "*.h5"
    # Get a list of files that match the pattern
    matching_files = glob.glob(file_pattern)
    # Loop through the matching files and delete each one
    for file_path in matching_files:
        os.remove(file_path)


def create_test_atoms():
    """Create a test ASE Atoms object."""
    return Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])


@pytest.fixture
def dirac_calculator():
    """Fixture to set up a test instance of the DIRAC calculator."""
    atoms = create_test_atoms()
    calc = DIRAC(molecule={'*basis': {'.default': '6-31G'}})
    atoms.calc = calc
    return atoms


@pytest.fixture(scope="session", autouse=True)
def clean_up_files():
    """Remove DIRAC calculation outputs at the end of each test."""
    yield
    command = "rm *.xyz* *.inp* *.out* *.h5* *.tgz* MDCINT* " \
        "MRCONEE* FCIDUMP* AOMOMAT* FCI*"
    subprocess.run(command, shell=True, capture_output=True)


def test_DIRAC_energy_rks():
    """Test case # 1 - RKS/B3LYP H2 molecule."""
    # define molecule geometry
    h2_molecule = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.7284]])

    # run calc and convert electronic energy into atomic units
    h2_molecule.calc = DIRAC(dirac={'.wave function': ''},
                             hamiltonian={'.nonrel': '', '.dft': 'b3lyp'},
                             molecule={'*basis': {'.default': 'STO-3G'}}
                             )
    energy_Eh = h2_molecule.get_potential_energy() / Ha

    # compare with the energy obtained using dirac alone
    # => assuming convergence up to ~1e-6
    assert energy_Eh == pytest.approx(-1.1587324754260060, 1e-6)


def test_DIRAC_energy_hf():
    """Test case # 2 - Testing ASE-DIRAC default parameters - HF/STO-3G."""
    h_atom = Atoms('H')

    h_atom.calc = DIRAC()
    energy_Eh = h_atom.get_potential_energy() / Ha

    # compare with the energy obtained using dirac alone
    # => assuming convergence up to ~1e-6
    assert energy_Eh == pytest.approx(-0.466581849557275, 1e-6)


def test_DIRAC_energy_mp2():
    """Test case # 3 - MP2/STO-3G Relativistic H2O.

    Notes:
        Adapted from the original DIRAC test set.
    """
    h2o_molecule = molecule('H2O')

    h2o_molecule.calc = DIRAC(hamiltonian={'.lvcorr': ''},
                              wave_function={'.scf': '', '.mp2': '',
                                             '*scf': {'.itrint': '5 50',
                                                      '.maxitr': '25'},
                                             '*mp2cal': {'.occup': '2..5',
                                                         '.virtual': 'all',
                                                         '.virthr': '2.0D00'}},
                              molecule={'*basis': {'.default': 'STO-3G'}})
    energy_Eh = h2o_molecule.get_potential_energy() / Ha

    # compare with the energy retrieved from DIRAC test set results.
    assert energy_Eh == pytest.approx(-75.043050906542774, 1e-6)


def test_DIRAC_energy_ccsdt():
    """Test case # 4 - CCSD(T)/STO-3G Relativistic H2O.

    Notes:
        Adapted from the original DIRAC test set.
        Uses integral transformation option **MOLTRA.
    """
    h2o_molecule = molecule('H2O')

    h2o_molecule.calc = DIRAC(dirac={'.wave function': '', '.4index': ''},
                              hamiltonian={'.lvcorr': ''},
                              wave_function={'.scf': '', '.relccsd': '',
                                             '*scf': {'.evccnv': '1.0E-8'}},
                              moltra={'.active': 'all'}
                              )
    energy_Eh = h2o_molecule.get_potential_energy() / Ha

    # compare with the energy retrieved from DIRAC test set results.
    assert energy_Eh == pytest.approx(-75.05762412870262, 1e-6)


def test_DIRAC_energy_open_shell():
    """Test case # 5 - open shell HF/STO-3G Relativistic C.

    Notes:
        Adapted from the original DIRAC open-shell test set.
    """
    c_atom = Atoms('C')

    c_atom.calc = DIRAC(hamiltonian={'.x2c': ''},
                        integrals={'*readin': {'.uncontracted': "#"}},
                        molecule={'*basis': {'.default': 'STO-3G'},
                                  '*charge': {'.charge': '0'},
                                  '*symmetry': {'.d2h': '#'}},
                        wave_function={'.scf': '',
                                       '*scf':
                                       {'.closed shell': '4 0',
                                        '.open shell': '2\n1/0,2\n1/0,4',
                                        '.kpsele':
                                        '3\n-1 1 -2\n4 0 0\n0 2 0\n0 0 4'}}
                        )
    energy_Eh = c_atom.get_potential_energy() / Ha

    # compare with the energy obtained using dirac alone
    # => assuming convergence up to ~1e-6
    assert energy_Eh == pytest.approx(-37.253756513429018, 1e-6)
