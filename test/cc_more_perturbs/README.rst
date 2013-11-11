Perturbation operators with the KU-CCSD(T) method
==================================================

With the RELCCSD keyword ADD_FINITE_FIELD and having operator
previously transformed in the PRPTRA module the user can add 
multiple one-electron operators to the (Kramers unrestricted) Fock matrix and
do the perturbed Kramers unrestricted Coupled Cluster procedure with. 

If operator is antisymmetric (what means imaginary), you have to include it
with the imaginary field strength to keep the Fock matrix Hermitian.

Also don't forget to check the proper point group symmetry of your perturbation operators.

Calculations are on your own risk.

In this test example, however, chosen perturbations are hurting the CCSD iteration process,
giving different energies for C1 and C2 symmetries (C2v is here only for comparison).

( C2v:  -192.357293399350084 )
   C2:  -192.357359545648364
   C1:  -192.357413401593305


ZDIPVEL (as real paerturb):
           CCSD(T)                 CCSD
   C2v: -192.357277468035534   -192.357277468035534
    C2: -192.357277507880099   -192.357277507880099
    C1: -192.357277538249690   -192.357277538249690

Z betasigma operator (as imaginary perturb):
              CCSD(T)
  C2v:   -192.357277468035591
   C2:   -192.357277468035477
   C1:   -192.357277468209219

