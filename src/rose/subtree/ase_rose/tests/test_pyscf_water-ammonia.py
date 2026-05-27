import os
import pytest

from ase_rose import ROSE
from ase_rose import ROSETargetMolecule, ROSEFragment

from .rose_test_functions import clean_up, extract_number, read_output


def run_water_ammonia():
    """Water-Ammonia mol_frag example calculation."""
    NH3OH2 = ROSETargetMolecule(name='water-ammonia',
                                atoms=[('N', (-1.39559, -0.02156,  0.00004)),
                                       ('H', (-1.62981,  0.96110, -0.10622)),
                                       ('H', (-1.86277, -0.51254, -0.75597)),
                                       ('H', (-1.83355, -0.33077,  0.86231)),
                                       ('O', (1.56850,  0.10589,  0.00001)),
                                       ('H', (0.60674, -0.03396, -0.00063)),
                                       ('H', (1.94052, -0.78000,  0.000))],
                                basis='sto-3g')

    NH3 = ROSEFragment(name='ammonia',
                       atoms=[('N', (-1.39559, -0.02156,  0.00004)),
                              ('H', (-1.62981,  0.96110, -0.10622)),
                              ('H', (-1.86277, -0.51254, -0.75597)),
                              ('H', (-1.83355, -0.33077,  0.86231))],
                       basis='sto-3g')

    H2O = ROSEFragment(name='water',
                       atoms=[('O', (1.56850,  0.10589,  0.00001)),
                              ('H', (0.60674, -0.03396, -0.00063)),
                              ('H', (1.94052, -0.78000,  0.000))],
                       basis='sto-3g')

    NH3OH2_calculator = ROSE(rose_calc_type='mol_frag',
                             exponent=2,
                             rose_target=NH3OH2,
                             rose_frags=[NH3, H2O],
                             test=True,
                             additional_virtuals_cutoff=2.0,
                             frag_threshold=10.0,
                             frag_valence=[[1, 7], [2, 6]],
                             frag_core=[[1, 1], [2, 1]],
                             avas_frag=[1,2],
                             nmo_avas=[[4, 5],[4, 5]])

    NH3OH2_calculator.calculate()


EXPECTED_OUTPUT_FILE = 'test_rose_ase_pyscf-water-ammonia.stdout'
ACTUAL_OUTPUT_FILE = 'OUTPUT_ROSE'

CHARGE_REGEX = 'Partial charges of fragments.*\n{}'.format(
    ('*\n.' * 6) + '*\n')
MO_ENERGIES_REGEX = 'Fragment          MO              Energy.*\n{}'.format(
    ('*\n.' * 6) + '*\n')
AVAS_EIGENVALUES = 'Eigenvalues of the Projected overlap.*\n{}'.format(('*\n.' * 10) + '*\n')


def run_calculation():
    """Calculation to run."""
    run_water_ammonia()


@pytest.fixture(scope="session", autouse=True)
def clean_up_files():
    """Runs always at the end of each test."""
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
    assert actual_charges == pytest.approx(expected_charges, rel=1e-3)


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
    assert actual_energies == pytest.approx(expected_energies, rel=1.0e-6)


def test_avas():
    """Test avas output."""
    # read the expected and actual outputs
    expected_output = read_output(
        os.path.join(os.path.dirname(__file__), 'data',
                     EXPECTED_OUTPUT_FILE)
    )
    actual_output = read_output(
        os.path.join(os.path.dirname(__file__), ACTUAL_OUTPUT_FILE)
    )
    # check the actual output against the expected output
    expected_eigenvalues = extract_number(AVAS_EIGENVALUES, expected_output)
    actual_eigenvalue = extract_number(AVAS_EIGENVALUES, actual_output)
    assert actual_eigenvalue == pytest.approx(expected_eigenvalues, rel=1.0e-6)
