

13.4 (this is the next patch)
-----------------------------

- Nicer error messages if the mol file is wrongly formatted.
- Bugfix for dipole moment calculations with the KRCI property module.
  The core-electronic part was not added to the electronic contribution 
  of the total dipole moment.
- Bugfix for X2Cmmf calculation (thanks to Kirk Peterson for pointing it out)


13.3 (2014-06-02)
-----------------

- Fix for more than nine buffer files in integral trnasformation
  (problem reported on dirac-users mailing list).
- Bugfix for correlated spinfree X2C calculations with low symmetry molecules (NZ > 1).
- Adapted testing framework to older versions of python.


13.2 (2014-04-14)
-----------------

- Simplified restart option for X2Cmmf calculations.
- Adapted testing framework to older versions of python.


13.1 (2013-12-16)
-----------------

- Magnetizabilities at the Hartree-Fock level with London atomic orbitals
  (LAOs) were accidentally deactivated in 13.0.
- Enabled numerical gradient geometry optimization for a selected Fock-space Coupled-Cluster
  state (to be picked from the last active sector).


13.0 (2013-12-12)
-----------------

- DIRAC13 release, see doc/release/release-statement.rst
