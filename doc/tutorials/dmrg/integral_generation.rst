

Generate a FCIDUMP file for external DMRG calculations
======================================================

This tutorial briefly demonstrates how to generate a FCIDUMP integral file which can be used to perform relativistic DMRG
calculations with an external DMRG code. The interface and the first relativistic DMRG implementation and 
calculation are described in [Knecht2014a]_. 

In the spinfree case one should use::

  **HAMILTONIAN
  .X2C
  .SPINFREE


This input should be used with the following pam line::

  $ ./pam ... --put=X2CMAT.* 

