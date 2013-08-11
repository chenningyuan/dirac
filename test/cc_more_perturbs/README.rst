Perturbation operators with the KU-CCSD(T) method
==================================================

With the keyword ADD_FINITE_FIELD and having operator
previously transformed in PRPTRA the user can add 
multiple operators to the Kramers unrestricted Fock matrix.

If operator is antisymmetric (what means imaginary), you have to include it
with the imaginary field strength to keep the Fock matrix Hermitian.

Also don't forget to check the proper point group symmetry of your perturbation operators.

Calculations are on your own risk.

In this test example, however, chosen perturbations are hurting the CCSD iteration process,
giving different energies for C1 and C2 symmetries (C2v is here only for comparison).

( C2v:  -192.357293399350084 )
   C2:  -192.357359545648364
   C1:   -192.357413401593305
