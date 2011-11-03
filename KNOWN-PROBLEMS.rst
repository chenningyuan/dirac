

g95 does not support cray pointers (MPI-2)
------------------------------------------

MPI-2 compatibility test fails with g95 and DIRAC is compiled without
-DVAR_MPI2.  There is not much to be done about it, apart from changing
compiler.


Modules not running in parallel
-------------------------------

* RELADC is not parallelized


Modules not running with 32 bit integers
----------------------------------------

* LUCITA
* KR-MCSCF (using LUCIAREL)
* KR-CI modules
* large Coupled-Cluster calculations


gfortran: Can't open module file
--------------------------------

Perhaps you are using gfortran older than 4.3? Please update your gfortran compiler.


ifort and "dummy reads"
-----------------------

Certain versions of the Intel Fortran compiler implement "dummy reads" in a
very stupid way. This makes seeking very slow in large files, which can affect
calculations on large molecules. The issue is reported to Intel, but affect
ifort at least up to version 10. As far as we know there is no solution except
to use a different compiler.
