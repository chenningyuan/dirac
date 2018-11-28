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
bond_length = 1.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()

# COMMENTS : For the moment, I don't manage to get the exact same molecular Hamiltonian in psi4 and Dirac...
# In consequence, the result is wrong with the Dirac interface. It does not seems like a very big issue,
# but it needs to be fixed of course. Otherwise, everything runs smoothly...

# Set calculation parameters. For the moment, I only interfaced the SCF calculation.
run_scf = 1
if run_scf==1:
 description = 'R' + str(bond_length) + '_scf'
run_mp2 = 0
if run_mp2==1:
 description = 'R' + str(bond_length) + '_mp2'
run_cisd = 0
if run_cisd==1:
 description = 'R' + str(bond_length) + '_cisd'
run_ccsd = 1
if run_ccsd==1:
 description = 'R' + str(bond_length) + '_ccsd'
run_fci = 1
delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
delete_FCIDUMP = False
verbose = False
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
                    verbose=verbose,
                    run_scf=run_scf,
                    run_mp2=run_mp2,
                    run_cisd=run_cisd,
                    run_ccsd=run_ccsd,
                    run_fci=run_fci)

with open("integrals_LiH_psi4","w") as f:
 f.write(str(molecule_psi4.get_molecular_hamiltonian()))
f.close
qbit_hamiltonian_psi4 = jordan_wigner(molecule_psi4.get_molecular_hamiltonian())
#print('Qubit Hamiltonian:\n')
#print(qbit_hamiltonian_psi4, '\n')
#new_qbit_hamiltonian_psi4 = project_onto_sector(qbit_hamiltonian_psi4, qubits=[1,3], sectors=[0,0])
#print('After projection onto sector 0:\n')
#print(new_qbit_hamiltonian_psi4, '\n')
evs_psi4 = eigenspectrum(qbit_hamiltonian_psi4)
#new_evs_psi4 = eigenspectrum(new_qbit_hamiltonian_psi4)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_psi4.hf_energy))
print('MP2 energy of {} Hartree.'.format(molecule_psi4.mp2_energy))
print('CCSD energy of {} Hartree.'.format(molecule_psi4.ccsd_energy))
print('FCI energy of {} Hartree.'.format(molecule_psi4.fci_energy))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): {}'.format(min(evs_psi4)))
#print('minimum eigenvalue with psi4 is, after projection: {}'.format(min(new_evs_psi4)))
#print('\n\n')
#print('One-body integrals from psi4:\n')
#print(molecule1.get_from_file("one_body_integrals"), '\n')
#print('two-body integrals from psi4:\n')
#print(molecule1.get_from_file("two_body_integrals"), '\n')
#print()

print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
relativistic=False
print(basis)

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
                    verbose=verbose,
                    run_scf=run_scf,
                    run_mp2=run_mp2,
                    run_cisd=run_cisd,
                    run_ccsd=run_ccsd,
                    run_fci=run_fci,
                    relativistic=relativistic)

with open("integrals_LiH_dirac","w") as f:
 f.write(str(molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic)))
f.close
qbit_hamiltonian_dirac = jordan_wigner(molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic))
#print('Qubit Hamiltonian:\n')
#print(qbit_hamiltonian_dirac,'\n')
evs_dirac = eigenspectrum(qbit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_dirac.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule_dirac.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule_dirac.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): {}'.format(min(evs_dirac)))

print()
print('#'*40)
print('REL Dirac calculation')
print('#'*40)
print()
relativistic=True
description += '_rel'

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
                    verbose=verbose,
                    run_scf=run_scf,
                    run_mp2=run_mp2,
                    run_cisd=run_cisd,
                    run_ccsd=run_ccsd,
                    run_fci=run_fci,
                    relativistic=relativistic)

with open("integrals_LiH_dirac_rel","w") as f:
 f.write(str(molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic)))
f.close
qbit_hamiltonian_dirac = jordan_wigner(molecule_dirac.get_molecular_hamiltonian(relativistic=relativistic))
#print('Qubit Hamiltonian:\n')
#print(qbit_hamiltonian_dirac,'\n')
evs_dirac = eigenspectrum(qbit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule_dirac.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule_dirac.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule_dirac.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): {}'.format(min(evs_dirac)))
