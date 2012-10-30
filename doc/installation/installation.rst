

Installation on Linux, Unix, and Mac
====================================

DIRAC is configured using CMake, typically via the setup script,
and subsequently compiled using make or gmake.


Dependencies
------------

To compile and run DIRAC you need CMake 2.8 or higher
and python 2.4 or higher.


What to do if CMake is not available or too old?
------------------------------------------------

It is your machine and you have Ubuntu or Debian::

  $ sudo apt-get install cmake cmake-curses-gui cmake-qt-gui

Similar mechanisms exist for other distributions or 
operating systems. Please consult Google.

If it is a cluster, please ask the Administrator to install/upgrade CMake.

If it is a cluster, but you prefer to install it yourself (it's easy):

1. Download the latest precompiled tarball from http://www.cmake.org/cmake/resources/software.html
2. Extract the tarball to, say, ~/cmake-2.8.5-Linux-i386
3. Set correct PATH variable


Basic installation using the setup script
-----------------------------------------

DIRAC is configured using CMake.
The setup script is a useful front-end to CMake.
It does nothing else than create the directory "build" and call
CMake with appropriate environment variables and flags::

  $ ./setup [--flags]
  $ cd build
  $ make

To see all available options::

  $ ./setup --help

You can see the CMake command using::

  $ ./setup [--flags] --show

and use it directly to call CMake without setup.


Most typical examples
---------------------

In order to get familiar with the configuration setup, let us demonstrate
some of the most typical configuration scenarios.

Configure for parallel compilation using MPI (make sure to properly export MPI
paths)::

  $ ./setup --fc=mpif90 --cc=mpicc

Configure for parallel compilation using MPI and 64-bit integers::

  $ ./setup --fc=mpif90 --cc=mpicc --int64

Configure for sequential compilation using ifort/icc::

  $ ./setup --fc=ifort --cc=icc

Configure for sequential compilation using gfortran/gcc::

  $ ./setup --fc=gfortran --cc=gcc

You get the idea. The configuration is normally pretty good at detecting math libraries
automatically, provided you export proper paths, see :ref:`linking_math`.


Out of source compilation
-------------------------

By default CMake builds out of source.  This means that all object files and
the final binary are generated outside of the src/ directory.  Typically the
build directory is called "build", but you can modify this default behavior::

  $ .setup [--flags] alternative-build-path

The out of source compilation is in contrast to previous DIRAC
releases.  This strategy offers several advantages. One obvious advantage is
that you can now build several binaries with the same source::

  $ cd /sourcepath
  $ ./setup --fc=gfortran --cc=gcc /gfortran-buildpath
  $ cd /gfortran-buildpath
  $ make
  $ cd /sourcepath
  $ ./setup --fc=ifort --cc=icc /ifort-buildpath
  $ cd /ifort-buildpath
  $ make


Basic installation without the setup script
-------------------------------------------

If you are familiar with CMake you don't have to use the setup script.
The setup script does nothing else than calling CMake with appropriate
environment variables and flags, it is a convenient front-end.

Minimal example::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make

The two following strategies are completely
equivalent:

Using CMake directly::

  $ mkdir build
  $ cd build
  $ FC=mpif90 CC=mpicc cmake -DENABLE_MPI=1 -DCMAKE_BUILD_TYPE=Release ..
  $ make

Using setup::

  $ ./setup --fc=mpif90 --cc=mpicc
  $ cd build
  $ make

If the compiler contains "mpi", then you can omit the flag --mpi, setup will set
it in this case automatically.

Please note that the defaults for performance optimization are different for
setup and direct CMake: by default setup configures for optimization, whereas
direct CMake commands configure code without optimization. Both defaults can be
changed.

There is nothing special about the directory "build".
You can do this instead::

  $ mkdir /buildpath
  $ cd /buildpath
  $ cmake /sourcepath
  $ make


.. _linking_math:


Linking to external math libraries
----------------------------------

Typically you will want to link to external math (BLAS and LAPACK) libraries,
for instance provided by MKL or Atlas.

The CMake configuration script will automatically find them if you define MATH_ROOT::

  $ export MATH_ROOT='/opt/intel/mkl'

Do not use full path MATH_ROOT='/opt/intel/mkl/lib/ia32'. CMake will append the
correct paths depending on the processor and the default integer type.  If the
MKL libraries that you want to use reside in
/opt/intel/mkl/10.0.3.020/lib/em64t, then MATH_ROOT is defined as::

  $ export MATH_ROOT='/opt/intel/mkl/10.0.3.020'

Then::

  $ ./setup [--flags]                 # do not need to define --math
  $ cd build
  $ make

Alternatively::

  $ cd build
  $ [FC=gfortran CC=gcc] MATH_ROOT='/opt/intel/mkl' cmake ..
  $ make

If automatic detection of math libraries fails for whatever reason, you can
always call the libraries explicitly like here::

  $ ./setup --explicit-libs="-L/path -lfoo -lbar"


Linking to Atlas libraries on SUSE or Fedora
--------------------------------------------

DIRAC may not be able to detect Atlas BLAS/Lapack
libraries on SUSE or Fedora because of nonstandard suffix::

  $ cd /usr/lib64/atlas
  $ ls -l

  lrwxrwxrwx. 1 root root   17 Oct 30 14:41 libatlas.so.3 -> ./libatlas.so.3.0
  -rwxr-xr-x. 1 root root 5.0M Sep  2  2011 libatlas.so.3.0
  lrwxrwxrwx. 1 root root   17 Oct 30 14:41 libcblas.so.3 -> ./libcblas.so.3.0
  -rwxr-xr-x. 1 root root 128K Sep  2  2011 libcblas.so.3.0
  lrwxrwxrwx. 1 root root   19 Oct 30 14:41 libclapack.so.3 -> ./libclapack.so.3.0
  -rwxr-xr-x. 1 root root  96K Sep  2  2011 libclapack.so.3.0
  lrwxrwxrwx. 1 root root   19 Oct 30 14:41 libf77blas.so.3 -> ./libf77blas.so.3.0
  -rwxr-xr-x. 1 root root 117K Sep  2  2011 libf77blas.so.3.0
  lrwxrwxrwx. 1 root root   18 Oct 30 14:41 liblapack.so.3 -> ./liblapack.so.3.0
  -rwxr-xr-x. 1 root root 5.4M Sep  2  2011 liblapack.so.3.0
  lrwxrwxrwx. 1 root root   19 Oct 30 14:41 libptcblas.so.3 -> ./libptcblas.so.3.0
  -rwxr-xr-x. 1 root root 128K Sep  2  2011 libptcblas.so.3.0
  lrwxrwxrwx. 1 root root   21 Oct 30 14:41 libptf77blas.so.3 -> ./libptf77blas.so.3.0
  -rwxr-xr-x. 1 root root 118K Sep  2  2011 libptf77blas.so.3.0

This can be fixed by creating soft links with .so suffix::

  $ su
  $ for file in *3; do ln -s $file ${file%.3}; done
  $ ls

  lrwxrwxrwx. 1 root root      13 Oct 30 17:45 libatlas.so -> libatlas.so.3
  lrwxrwxrwx. 1 root root      17 Oct 30 14:41 libatlas.so.3 -> ./libatlas.so.3.0
  -rwxr-xr-x. 1 root root 5225848 Sep  2  2011 libatlas.so.3.0
  lrwxrwxrwx. 1 root root      13 Oct 30 17:45 libcblas.so -> libcblas.so.3
  lrwxrwxrwx. 1 root root      17 Oct 30 14:41 libcblas.so.3 -> ./libcblas.so.3.0
  -rwxr-xr-x. 1 root root  130824 Sep  2  2011 libcblas.so.3.0
  lrwxrwxrwx. 1 root root      15 Oct 30 17:45 libclapack.so -> libclapack.so.3
  lrwxrwxrwx. 1 root root      19 Oct 30 14:41 libclapack.so.3 -> ./libclapack.so.3.0
  -rwxr-xr-x. 1 root root   97752 Sep  2  2011 libclapack.so.3.0
  lrwxrwxrwx. 1 root root      15 Oct 30 17:45 libf77blas.so -> libf77blas.so.3
  lrwxrwxrwx. 1 root root      19 Oct 30 14:41 libf77blas.so.3 -> ./libf77blas.so.3.0
  -rwxr-xr-x. 1 root root  119760 Sep  2  2011 libf77blas.so.3.0
  lrwxrwxrwx. 1 root root      14 Oct 30 17:45 liblapack.so -> liblapack.so.3
  lrwxrwxrwx. 1 root root      18 Oct 30 14:41 liblapack.so.3 -> ./liblapack.so.3.0
  -rwxr-xr-x. 1 root root 5566480 Sep  2  2011 liblapack.so.3.0
  lrwxrwxrwx. 1 root root      15 Oct 30 17:45 libptcblas.so -> libptcblas.so.3
  lrwxrwxrwx. 1 root root      19 Oct 30 14:41 libptcblas.so.3 -> ./libptcblas.so.3.0
  -rwxr-xr-x. 1 root root  130824 Sep  2  2011 libptcblas.so.3.0
  lrwxrwxrwx. 1 root root      17 Oct 30 17:45 libptf77blas.so -> libptf77blas.so.3
  lrwxrwxrwx. 1 root root      21 Oct 30 14:41 libptf77blas.so.3 -> ./libptf77blas.so.3.0
  -rwxr-xr-x. 1 root root  119824 Sep  2  2011 libptf77blas.so.3.0

Now DIRAC will correctly detect them.


Running CMake using GUI
-----------------------

You prefer GUI? No problem. You can configure with GUI::

  $ cd build
  $ cmake ..
  $ cmake-gui ..

You may have to install cmake-gui for it, on debian/ubuntu::

  $ sudo apt-get install cmake cmake-curses-gui cmake-qt-gui


Running tests
-------------

You can (and should) run the test suite with::

  $ make test

It is HIGHLY recommended to run the test set after you have compiled
DIRAC to make sure that your binary correctly reproduces reference results.


Make install
------------

Make install is very useful to make DIRAC available to other users on the same
machine::

  $ ./setup [--flags] --prefix=/path
  $ cd build
  $ make
  $ make install


Where should $PATH point to? Source directory or build directory?
-----------------------------------------------------------------

We recommend to let $PATH point to the install directory::

  $ ./setup [--flags] --prefix=/install/path
  $ cd build
  $ make
  $ make install

This way everything (binary, scripts, basis sets) will be at the right place
under /install/path and $PATH should contain /install/path.


Compiling in verbose mode
-------------------------

Sometimes you want to see the actual compiler flags and definitions::

  $ make VERBOSE=1


Compiling on many cores
-----------------------

Yes, it works. Try::

  $ make -j4

We have successfully compiled DIRAC on 64 cores. With the new configuration
based on CMake, compilation race condition errors do not appear.


How can I change optimization flags?
------------------------------------

If you want to turn optimization off (debug mode), there are several ways to do that.

Either use setup::

  $ ./setup --debug [other flags]
  $ cd build
  $ make

Or use CMake directly (default here is debug mode)::

  $ mkdir build
  $ cd build
  $ [FC=ifort CC=icc] cmake ..
  $ make
