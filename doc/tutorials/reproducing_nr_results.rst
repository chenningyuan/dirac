

How to reproduce nonrelativistic results
========================================

To approach the nonrelativistic limit you can
either switch to the 4-component Lévy-Leblond Hamiltonian::

  **DIRAC
  [...]
  **HAMILTONIAN
  .LEVY-LEBLOND

or the 2-component nonrelativistic Hamiltonian::

  **DIRAC
  [...]
  **HAMILTONIAN
  .NONREL

or increase the speed of light (here to 2000 a.u.)::

  **DIRAC
  [...]
  **GENERAL
  .CVALUE
  2000.0

Most nonrelativistic programs employ a point nucleus model.
You can force point nuclei also in DIRAC::

  **DIRAC
  [...]
  **INTEGRALS
  .NUCMOD
   1
