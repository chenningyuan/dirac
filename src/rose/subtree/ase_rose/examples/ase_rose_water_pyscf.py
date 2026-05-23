from ase_rose import ROSE
from ase_rose import ROSETargetMolecule, ROSEFragment

H2O = ROSETargetMolecule(
    name='water',
    atoms=[('O', (0.,  0.00000,  0.59372)),
           ('H', (0.,  0.76544, -0.00836)),
           ('H', (0., -0.76544, -0.00836))],
    basis='sto-3g'
)

oxygen = ROSEFragment(
    name='oxygen',
    atoms=[('O', (0, 0, 0))],
    multiplicity=3, basis='sto-3g'
)

hydrogen = ROSEFragment(
    name='hydrogen',
    atoms=[('H', (0, 0, 0))],
    multiplicity=2, basis='sto-3g'
)

H2O_calculator = ROSE(
    rose_calc_type='atom_frag',
    exponent=4,
    rose_target=H2O,
    rose_frags=[oxygen, hydrogen],
    test=True,
    avas_frag=[0],
    nmo_avas=[2, 3, 4],
    rose_mo_calculator='pyscf'
)

H2O_calculator.calculate()
