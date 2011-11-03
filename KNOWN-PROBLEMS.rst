

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
