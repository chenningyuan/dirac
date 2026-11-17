from ase_rose import ROSE
from ase_rose import ROSETargetMolecule, ROSEFragment

NH3OH2 = ROSETargetMolecule(
    name='water-ammonia',
    atoms=[('N', (-1.39559, -0.02156,  0.00004)),
           ('H', (-1.62981,  0.96110, -0.10622)),
           ('H', (-1.86277, -0.51254, -0.75597)),
           ('H', (-1.83355, -0.33077,  0.86231)),
           ('O', (1.56850,  0.10589,  0.00001)),
           ('H', (0.60674, -0.03396, -0.00063)),
           ('H', (1.94052, -0.78000,  0.000))],
    basis='sto-3g'
)

NH3 = ROSEFragment(
    name='ammonia',
    atoms=[('N', (-1.39559, -0.02156,  0.00004)),
           ('H', (-1.62981,  0.96110, -0.10622)),
           ('H', (-1.86277, -0.51254, -0.75597)),
           ('H', (-1.83355, -0.33077,  0.86231))],
    basis='sto-3g'
)

H2O = ROSEFragment(
    name='water',
    atoms=[('O', (1.56850,  0.10589,  0.00001)),
           ('H', (0.60674, -0.03396, -0.00063)),
           ('H', (1.94052, -0.78000,  0.000))],
    basis='sto-3g'
)

H2O_calculator = ROSE(
    rose_calc_type='mol_frag',
    exponent=2,
    rose_target=NH3OH2,
    rose_frags=[NH3, H2O],
    test=True,
    additional_virtuals_cutoff=2.0,
    frag_threshold=10.0,
    frag_valence=[[1, 7], [2, 6]],
    frag_core=[[1, 1], [2, 1]],
    avas_frag=[1],
    nmo_avas=[[2, 3, 4]],
    rose_mo_calculator='psi4'
)

H2O_calculator.calculate()
