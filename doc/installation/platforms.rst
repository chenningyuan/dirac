

Supported platforms and compilers
=================================

DIRAC requires Fortran 90, C, and C++ compilers.  The program is designed to
run on a unix-like operating system, and uses MPI for parallel calculations.
Optionally DIRAC can (and should) be linked to platform specific BLAS and
LAPACK libraries.


Tested on the following platforms
---------------------------------

* ClustervisionOS 4.0, GNU 4.3.6, system native BLAS/LAPACK
* ClustervisionOS 4.0, Intel 11.1, MKL 11.1, system native BLAS/LAPACK
* CrunchBang 11, GNU 4.7.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Ubuntu 12.04, GNU 4.6.3, MKL 11.0, OpenMPI 1.6.5
* Ubuntu 12.04, Intel 14.0.0, MKL 11.0, OpenMPI 1.6.5
* Mac OS X Lion 10.7.5, GNU 4.6.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Mac OS X Mavericks 10.9, GNU 4.8.2, system native BLAS/LAPACK
* OpenSUSE 12.2, Intel 12.1.0, MKL 10.3.6, OpenMPI 1.6.4
* Linux, Intel 13.0.1, MKL 11.0, Intel MPI 4.1
* Linux, Intel 12.1.9, MKL 10.3, Bull X MPI 1.1


Compiles but does not pass all tests
------------------------------------

* Windows 7/8, system native BLAS/LAPACK, GNU MinGW32/64, only tiny portion of tests fails


Compiles but untested on the following environments
---------------------------------------------------

* SGI Altix ICE 8200, Intel 11.1
* Cray XE6, Intel 12.1.5
* IBM AIX, XL Fortran V13.1, XL C V11.1, g++, builtin math

Additional dependencies due to PCMSolver
----------------------------------------

It is optionally possible to enable the polarizable continuum model functionality. 
The functionality is provided using the external module 
`PCMSolver <http://pcmsolver.github.io/pcmsolver-doc/>`_. The module requires some additional 
dependencies:

* the Python interpreter, version 2.6 or higher;
* the Python development libraries and headers, version 2.6 or higher;
* the zlib compression library, version 1.2 or higher.  

By default, compilation of the module is enabled. CMake will check that the additional
dependencies are met. In case they are not, compilation of PCMSolver will be disabled.
CMake will print out which dependencies were not satisfied.
Passing the option -DENABLE_PCMSOLVER=OFF to the setup script will disable compilation
of the module and skip detection of additional dependencies.
