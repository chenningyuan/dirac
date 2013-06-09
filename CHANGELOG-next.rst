

12.7
----

- DSO integrals can be limited by the .SELECT keyword. Previously DSO integral code
  ignored .SELECT.
- Extended addressing of DSO integrals up to 333 atom centers.
- Compilation fix for 64bit integer Intel MKL.
- DIRAC compiles using Portland group compilers.
- Introduced block for any structure optimization runs using the xyz format.
- Removed --force option from setup script.
- Pople basis sets containing "*" have been renamed to "-star" file names. However you do
  not need to change your input files. DIRAC will understand "*" characters in Pople basis sets
  and internally translate them.
- Bugfix for permanent dipole moment calculations with KRCI. Missing nuclear contribution term 
  have been added and the correct sign for the dipole integrals (length gauge) has been introduced.
- Bugfix for MP2-NO calculations. correct writing of 2-electron file MDCINT.
- Workaround for memory allocator segmentation fault with gfortran can now be triggered with:
  ./setup -D CPP="-DGFORTRAN_ALLOCATOR_WORKAROUND"
- improved final state analysis for MCSCF wave functions: print natural orbital occupation numbers by default
- BUGfix for parallel CC/MP2 calculations using the RELCCSD module and the new input style (the old Namelist input was
  fine); thanks to Kirk Peterson for pointing it out.
