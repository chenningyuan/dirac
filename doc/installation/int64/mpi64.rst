:orphan:

How to build MPI libraries for 64-bit integers
----------------------------------------------

If you want to use 64-bit integers in parallel DIRAC, you also need 64-bit integers in your MPI installation.

Please check the integer type of your currently available Open MPI installation (perhaps it can do 64-bit integers already).
For this, type:

::

 ompi_info -a | grep 'Fort integer size'

if the output contains 8 for 64-bit integers

::

  Fort integer size: 8

then you have a suitable 64-bit Open MPI installation which uses 64-bit integers by default. If the outpus shows 4, you have 32-bit integers.

If you decided that you want to build one yourself (not difficult),
then download the latest stable source from http://www.open-mpi.org/ and extract it.
Then enter the directory and configure Open MPI (edit prefix path).

For Intel compilers use:

::

 ./configure CXX=icpc CC=icc F77=ifort FC=ifort FFLAGS=-i8  FCFLAGS=-i8  CFLAGS=-m64  CXXFLAGS=-m64 --prefix=/path/to/your_openmpi

For GNU compilers type:

::

 ./configure  CXX=g++ CC=gcc F77=gfortran FC=gfortran FFLAGS="-m64 -fdefault-integer-8" FCFLAGS="-m64 -fdefault-integer-8" CFLAGS=-m64 CXXFLAGS=-m64 --prefix=/path/to/your_openmpi

Once the configure step is finished (lot of output, at the end there should be no error), build the library:

::

 make -jN
 make all install


Finally export the path locations in .bashrc (or in the corresponding file if you don't use bash-shell):

::

  export PATH=$PATH:'/path/to/your_openmpi/bin'
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'/path/to/your_openmpi/lib'
  export MANPATH=$MANPATH:'/path/to/your_openmpi/share/man'

