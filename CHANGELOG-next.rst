

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
