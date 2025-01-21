from ase import Atoms
from ase.build import molecule
from ase.units import Ha
from ase_dirac import DIRAC

results = []
for mol in ['H2O','Cl2','C2H4']:
    small_molecule = molecule(mol)
    try:
        small_molecule.calc = DIRAC(hamiltonian={'.x2cmmf': ''},
                          wave_function={'.scf': '', '.mp2': '',
                                         '*scf': {'.ergcnv': '1.E-8 1.E-6',
                                                  '.maxitr': '35'},
                                         '*mp2cal': {'.occup': 'all',
                                                     '.virtual': 'all',
                                                      }},
                          molecule={'*basis': {'.default': 'cc-pVDZ'}})
        energy_Eh = small_molecule.get_potential_energy() / Ha
        results.append([mol,energy_Eh])
    except:
        print (f'Molecule {mol} could not be run correctly, please check')
    
print('Relativistic MP2 energies using the molecular mean field approach')
for result in results:
    print (f'Molecule: {result[0]} Energy: {result[1]:12.6f} Hartree.')

