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
multiplicity = 1
charge = 0
data_directory=os.getcwd()

run_fci = 1 
if run_fci==1:
 description = 'fci_psi4'
delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
delete_FCIDUMP = False
geometry = [('Be', (0., 0., 0.))]

molecule = MolecularData(geometry=geometry,
                          basis=basis,
                          multiplicity=multiplicity,
                          charge=charge,
                          description=description,
                          data_directory = data_directory)


print('#'*40)
print('Psi4 calculation')
print('#'*40)
print()

molecule = run_psi4(molecule,  
                    delete_input=delete_input,
                    delete_output=delete_output,
                    run_fci=run_fci)

fermion_hamiltonian = molecule.get_molecular_hamiltonian()
qubit_hamiltonian_psi4 = jordan_wigner(fermion_hamiltonian)
evs_psi4 = eigenspectrum(qubit_hamiltonian_psi4)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.hf_energy))
print('MP2 energy of {} Hartree.'.format(molecule.mp2_energy))
print('CCSD energy of {} Hartree.'.format(molecule.ccsd_energy))
print('FCI energy of {} Hartree.'.format(molecule.fci_energy))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_psi4))



print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
run_ccsd = 1
if run_ccsd==1:
 description = 'ccsd_dirac'

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
                    run_ccsd=run_ccsd)

fermion_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))
#print('Solving the Qubit Hamiltonian (Bravyi-Kitaev): {}'.format(min(evs_dirac_BK)))

print()
print('#'*40)
print('REL Dirac calculation')
print('#'*40)
print()
relativistic=True

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               relativistic=relativistic,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_ccsd=run_ccsd,
                    relativistic=relativistic)

fermion_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
#evs_dirac_BK = eigenspectrum(qbit_hamiltonian_dirac_BK)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))


speed_of_light = 3.0
print()
print('#'*60)
print('REL Dirac calculation with speed of light = ', speed_of_light, ' a.u.')
print('#'*60)
print()

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               speed_of_light=speed_of_light,
                               description=description,
                               relativistic=relativistic,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_ccsd=run_ccsd,
                    relativistic=relativistic,
                    speed_of_light=speed_of_light)

fermion_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))

subprocess.call("mkdir OUTPUT_PSI4 INTEGRALS OUTPUT_DIRAC/", shell=True)
subprocess.call("mv -f *dirac*out OUTPUT_DIRAC/", shell=True)
subprocess.call("mv -f *psi4*out *psi4*hdf5 OUTPUT_PSI4/", shell=True)
subprocess.call("mv -f FCIDUMP_* INTEGRALS/", shell=True)
