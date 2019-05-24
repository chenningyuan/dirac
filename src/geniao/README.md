Program to construct the Intrinsic Atomic Orbitals (IAOs).

Executable:
./geniao.x

Input files:
MOLECULE.XYZ (geometry file)
MOLECULE.dfcoef (DFCOEF file for the full molecule)
{x}.dfcoef ({x} being the atomic number of the fragments, i.e. 001.dfcoef for hydrogen, 008.dfcoef for oxygen, etc...)

Output files:
DFCOEF_iao and DFPCMO_iao which contain the coefficient matrix of the IAOs in term of the AOs.

Remark:
The fragment.dfcoef should be obtained with a large basis set.
The MOLECULE.dfcoef can be obtained from any standard basis set.
