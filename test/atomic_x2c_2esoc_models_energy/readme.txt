=================== description of test =====================

This test calculates the ground state energy of (LiHe)+ at
the HF level using the exact two-component (X2C) hamiltonian 
in the implementation interface to InteRest by M. Repisky.
the calculations are done using various X2C "flavors" listed in the following 
and described in more detail in the paper XXX:

--> molecular X2C approach (using the bare-nucleus h1 hamiltonian to derive the transformation matrix U) 
  + AMFI 2e-spin-orbit corrections                            (molecular-x2c-amfisoc.inp)
  or
  + atomic spin-orbit corrections from 4c-DC SCF calculations (molecular-x2c-aoosoc.inp)

--> molecular X2C approach using the transformation matrix U from a molecular 4c-DC SCF (h1 == fock operator)
  + atomic spin-orbit corrections from 4c-DC SCF calculations (molecular-fragment-x2c-aoosoc.inp)

--> atomic X2C approach using the transformation matrices U from atomic calculations 
  + atomic spin-orbit corrections from 4c-DC SCF calculations (atomic-x2c-aoosoc.inp)

bug reports to dirac-users@googlegroups.com or stefan.knecht@phys.chem.ethz.ch
=================== =================== =====================
