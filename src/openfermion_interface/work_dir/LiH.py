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

run_fci = 1 
if run_fci==1:
 description = 'R' + str(bond_length) + '_fci_psi4'
delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
delete_FCIDUMP = False
geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., bond_length))]

molecule1 = MolecularData(geometry=geometry,
                          basis=basis,
                          multiplicity=multiplicity,
                          charge=charge,
                          description=description,
                          data_directory = data_directory)


print('#'*40)
print('Psi4 calculation')
print('#'*40)
print()

molecule_psi4 = run_psi4(molecule1,  
                    delete_input=delete_input,
                    delete_output=delete_output,
                    run_fci=run_fci)

fermion_hamiltonian = molecule_psi4.get_molecular_hamiltonian()
qubit_hamiltonian_psi4 = jordan_wigner(fermion_hamiltonian)
evs_psi4 = eigenspectrum(qubit_hamiltonian_psi4)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_psi4.hf_energy))
print('MP2 energy of {} Hartree.'.format(molecule_psi4.mp2_energy))
print('CCSD energy of {} Hartree.'.format(molecule_psi4.ccsd_energy))
print('FCI energy of {} Hartree.'.format(molecule_psi4.fci_energy))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_psi4))



print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
run_scf = 1
if run_scf==1:
 description = 'R' + str(bond_length) + '_scf_dirac'
run_ccsd = 1
if run_ccsd==1:
 description = 'R' + str(bond_length) + '_ccsd_dirac'

molecule2 = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule_dirac = run_dirac(molecule2,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_scf=run_scf,
                    run_ccsd=run_ccsd)

fermion_hamiltonian = molecule_dirac.get_molecular_hamiltonian()
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_dirac.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule_dirac.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule_dirac.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))
#print('Solving the Qubit Hamiltonian (Bravyi-Kitaev): {}'.format(min(evs_dirac_BK)))

print()
print('#'*40)
print('REL Dirac calculation')
print('#'*40)
print()
relativistic=True

molecule2 = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               relativistic=relativistic,
                               description=description,
                               data_directory=data_directory)

molecule_dirac = run_dirac(molecule2,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_scf=run_scf,
                    run_ccsd=run_ccsd,
                    relativistic=relativistic)

# Don't forget to specify the relativistic argument in get_molecular_hamiltonian if relavistic is True.
fermion_hamiltonian = molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic)
qubit_hamiltonian_dirac = jordan_wigner(molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic))
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
#evs_dirac_BK = eigenspectrum(qbit_hamiltonian_dirac_BK)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_dirac.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule_dirac.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule_dirac.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))


speed_of_light = 3.0
print()
print('#'*60)
print('REL Dirac calculation with speed of light = ', speed_of_light, ' a.u.')
print('#'*60)
print()

molecule2 = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               speed_of_light=speed_of_light,
                               description=description,
                               data_directory=data_directory)

molecule_dirac = run_dirac(molecule2,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_scf=run_scf,
                    run_ccsd=run_ccsd,
                    relativistic=relativistic,
                    speed_of_light=speed_of_light)

fermion_hamiltonian = molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic)
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_dirac.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule_dirac.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule_dirac.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))


print()
print('#'*40)
print('NONREL Dirac calculation, cc-pVDZ basis for H')
print('#'*40)
print()
relativistic=False
basis="special"
special_basis=".SPECIAL\nH BASIS cc-pVDZ\n"
molecule2 = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               special_basis=special_basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)


molecule_dirac = run_dirac(molecule2,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_scf=run_scf,
                    run_ccsd=run_ccsd,
                    relativistic=relativistic,
                    speed_of_light=speed_of_light)

fermion_hamiltonian = molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic)
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_dirac.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule_dirac.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule_dirac.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))

subprocess.call("mkdir OUTPUT_PSI4 INTEGRALS OUTPUT_DIRAC/", shell=True)
subprocess.call("mv -f *dirac*out OUTPUT_DIRAC/", shell=True)
subprocess.call("mv -f *psi4*out *.hdf5 OUTPUT_PSI4/", shell=True)
subprocess.call("mv -f FCIDUMP_* INTEGRALS/", shell=True)
