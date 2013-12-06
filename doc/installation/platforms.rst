

Supported platforms and compilers
=================================

DIRAC requires Fortran 90, C, and C++ compilers.  The program is designed to
run on a unix-like operating system, and uses MPI for parallel calculations.
Optionally DIRAC can (and should) be linked to platform specific BLAS and
LAPACK libraries.


Tested on the following platforms
---------------------------------

* CrunchBang 11, GNU 4.7.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Ubuntu 12.04, GNU 4.6.3, MKL 11.0, OpenMPI 1.6.5
* Ubuntu 12.04, Intel 14.0.0, MKL 11.0, OpenMPI 1.6.5


Compiles but does not pass all tests
------------------------------------

* Windows 7/8, system native BLAS/LAPACK, GNU MinGW32/64, only tiny portion of tests fails


Compiles but untested on the following environments
---------------------------------------------------

* Ubuntu 12.04, Intel 14.0.0, MKL 11.0
* SGI Altix ICE 8200, Intel 11.1
* Cray XE6, Intel 12.1.5


Currently does not compile
--------------------------

* IBM AIX, XL 7.1
