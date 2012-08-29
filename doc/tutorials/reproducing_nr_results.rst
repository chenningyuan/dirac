

How to reproduce nonrelativistic results
========================================

DIRAC nonrelativistic Hamiltonians
-----------------------------------
To approach the nonrelativistic limit you can
either switch to the 4-component Lévy-Leblond Hamiltonian::

  **DIRAC
  [...]
  **HAMILTONIAN
  .LEVY-LEBLOND

or the 2-component nonrelativistic (Schrodinger) Hamiltonian::

  **DIRAC
  [...]
  **HAMILTONIAN
  .NONREL

Changing speed of light
-----------------------
Another option to reproduce nonrelativistic data with the Dirac relativistic Hamiltonian
is to increase the speed of light (here to 2000 a.u.)::

  **DIRAC
  [...]
  **GENERAL
  .CVALUE
  2000.0

Adjusting the nucleus model
---------------------------
Most nonrelativistic programs employ a point nucleus model.
You can force point nuclei also in DIRAC::

  **DIRAC
  [...]
  **INTEGRALS
  .NUCMOD
   1

Basis sets
----------
With the above mentioned nonrelativistic Hamiltonians the user can also employ variety of standard
nonrelativistic basis sets present in the "basis_dalton" directory.

Altogether, to reproduce results of most nonrelativistic packages the user should set nonrelativistic
Hamiltonian, switch to point nucleus if necessary for the compatibility, 
and choose some nonrelativistic basis set.

Recommended method of choice is the closed shell Hartree-Fock SCF, which can be followed
by the CCSD(T) correlation method. In this case you are to check the occupied shells, and
active space, respectively.

