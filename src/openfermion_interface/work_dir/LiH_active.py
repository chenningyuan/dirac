# Test set up for generating Hamiltonian for H2.

from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.hamiltonians import MolecularData
from openfermion.transforms import jordan_wigner, project_onto_sector, bravyi_kitaev
from openfermion.utils import count_qubits,eigenspectrum
from openfermionpsi4 import run_psi4
import os
import subprocess

# Set molecule parameters.
basis = 'sto-3g'
bond_length = 2.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()
delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True 
delete_FCIDUMP = False
geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*60)
print('NONREL Dirac calculation with active orbitals')
print('#'*60)
print()
run_scf = 1
if run_scf==1:
 description = 'R' + str(bond_length) + '_scf_dirac'
run_ccsd = 1
if run_ccsd==1:
 description = 'R' + str(bond_length) + '_ccsd_dirac'

active = [-2.0,0.3,0.1]
description += '_active'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_scf=run_scf,
                    run_ccsd=run_ccsd,
                    active=active)

fermion_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(min(evs_dirac)))

subprocess.call("mkdir OUTPUT_PSI4 INTEGRALS OUTPUT_DIRAC/", shell=True)
subprocess.call("mv -f *dirac*out OUTPUT_DIRAC/", shell=True)
subprocess.call("mv -f *psi4*out *psi4*hdf5 OUTPUT_PSI4/", shell=True)
subprocess.call("mv -f FCIDUMP_* INTEGRALS/", shell=True)
