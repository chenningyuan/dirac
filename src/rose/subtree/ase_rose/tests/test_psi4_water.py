"""
Rose tests for psi4.
"""
import os
import pytest

from ase_rose import ROSE
from ase_rose import ROSETargetMolecule, ROSEFragment

from .rose_test_functions import clean_up, extract_number, read_output


def run_water_rose():
    """Water atom_frag example calculation."""
    H2O = ROSETargetMolecule(
        name='water',
        atoms=[('O', (0.,  0.00000,  0.59372)),
               ('H', (0.,  0.76544, -0.00836)),
               ('H', (0., -0.76544, -0.00836))],
        basis='sto-3g',
        uncontract=True
    )

    oxygen = ROSEFragment(
        name='oxygen',
        atoms=[('O', (0, 0, 0))],
        multiplicity=1, basis='sto-3g',
        uncontract=True
    )

    hydrogen = ROSEFragment(
        name='hydrogen',
        atoms=[('H', (0, 0, 0))],
        multiplicity=2, basis='sto-3g',
        uncontract=True
    )

    H2O_calculator = ROSE(rose_calc_type='atom_frag',
                          exponent=4,
                          rose_target=H2O,
                          rose_frags=[oxygen, hydrogen],
                          test=True,
                          avas_frag=[0],
                          nmo_avas=[2, 3, 4],
                          save_data=False,
                          restricted=True,
                          openshell=True,
                          rose_mo_calculator='psi4')
    # run the calculator
    H2O_calculator.calculate()


EXPECTED_OUTPUT_FILE = 'test_rose_ase_psi4-water.stdout'
ACTUAL_OUTPUT_FILE = 'OUTPUT_ROSE'

AVAS_EXPECTED_OUTPUT_FILE = 'test_rose_avas_ase_psi4-water.stdout'

CHARGE_REGEX = 'Partial charges of fragments.*\n{}'.format(
    ('*\n.' * 6) + '*\n'
)
MO_ENERGIES_REGEX = 'Recanonicalized virtual '
'energies of fragments.*\n{}'.format(('*\n.' * 9) + '*\n')

AVAS_EIGENVALUES = 'Eigenvalues of the Projected overlap.*\n{}'.format(('*\n.' * 7) + '*\n')


def run_calculation():
    """Calculation to run."""
    run_water_rose()


@pytest.fixture(scope="session", autouse=True)
def clean_up_files():
    """Runs always at the end of all tests."""
    yield
    clean_up()


def test_partial_charges():
    """Test partial charges of fragments."""
    # run the calculation
    run_calculation()

    # read the expected and actual outputs
    expected_output = read_output(
        os.path.join(os.path.dirname(__file__), 'data', EXPECTED_OUTPUT_FILE)
    )
    actual_output = read_output(
        os.path.join(os.path.dirname(__file__), ACTUAL_OUTPUT_FILE)
    )

    # check the actual output against the expected output
    expected_charges = extract_number(CHARGE_REGEX, expected_output)
    actual_charges = extract_number(CHARGE_REGEX, actual_output)
    assert actual_charges == pytest.approx(expected_charges, rel=1e-2)


def test_virtual_energies():
    """Test recanonicalized virtual energies of fragments."""
    # read the expected and actual outputs
    expected_output = read_output(
        os.path.join(os.path.dirname(__file__), 'data', EXPECTED_OUTPUT_FILE)
    )
    actual_output = read_output(
        os.path.join(os.path.dirname(__file__), ACTUAL_OUTPUT_FILE)
    )

    # check the actual output against the expected output
    expected_energies = extract_number(MO_ENERGIES_REGEX, expected_output)
    actual_energies = extract_number(MO_ENERGIES_REGEX, actual_output)
    assert actual_energies == pytest.approx(expected_energies, rel=1.0e-4)


def test_avas():
    """Test avas output."""
    # read the expected and actual outputs
    expected_output = read_output(
        os.path.join(os.path.dirname(__file__), 'data',
                     AVAS_EXPECTED_OUTPUT_FILE)
    )
    actual_output = read_output(
        os.path.join(os.path.dirname(__file__), ACTUAL_OUTPUT_FILE)
    )
    # check the actual output against the expected output
    expected_eigenvalues = extract_number(AVAS_EIGENVALUES, expected_output)
    actual_eigenvalue = extract_number(AVAS_EIGENVALUES, actual_output)
    assert actual_eigenvalue == pytest.approx(expected_eigenvalues, rel=1.0e-4)
