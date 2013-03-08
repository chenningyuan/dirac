=================== description of test =====================

This test calculates the ground state energy of (LiHe)+ at
the HF level using the exact two-component (X2C) hamiltonian 
in the implementation interface to InteRest by M. Repisky.
the calculation is done using a 
- molecular X2C approach                                                         
  +
	- AMFI 2e-spin-orbit corrections                            (molecular-x2c-amfisoc.inp)
	- atomic spin-orbit corrections from 4c-DC SCF calculations (molecular-x2c-aoosoc.inp)
- atomic X2C approach using the transformation matrices from atomic calculations 
  + 
	- atomic spin-orbit corrections from 4c-DC SCF calculations (atomic-x2c-aoosoc.inp)

bug reports to dirac-users@googlegroups.com or stefan.knecht@phys.chem.ethz.ch
=================== =================== =====================
