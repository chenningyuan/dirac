:orphan:
 

Supported platforms and compilers
=================================

DIRAC requires Fortran 90, C, and C++ compilers.  The program is designed to
run on a unix-like operating system, and uses MPI for parallel calculations.
Optionally DIRAC can (and should) be linked to platform specific BLAS and
LAPACK libraries.


Additional dependencies due to PCMSolver
----------------------------------------

It is optionally possible to enable the polarizable continuum model functionality. 
The functionality is provided using the external module 
`PCMSolver <http://pcmsolver.github.io/pcmsolver-doc/>`_. The module requires additional 
dependencies:

* the zlib compression library, version 1.2 or higher;
* the Boost set of libraries, version 1.54 or higher.

The module performs checks on the system to find an installation of Boost that has
the suitable version. The search is directed towards conventional system directories.
If the search fails, the module will build a copy of Boost on the spot.
If your system has Boost installed in a nonstandard location, e.g. /boost/location, you have to direct the
search by adding the following flags to the ``setup`` script::

  $ ./setup -DBOOST_INCLUDEDIR=/boost/location/include -DBOOST_LIBRARYDIR=/boost/location/lib

By default, compilation of the module is enabled. CMake will check that the additional
dependencies are met. In case they are not, compilation of PCMSolver will be disabled.
CMake will print out which dependencies were not satisfied.
Passing the option ``-DENABLE_PCMSOLVER=OFF`` to the ``setup`` script will disable compilation
of the module and skip detection of additional dependencies.
