

How to build MPI libraries for 64-bit integers
==============================================

Before you continue please verify whether you really need 64-bit integers.
It is extra work and perhaps not needed.

We are currently aware of two MPI distributions that can either be
installed as 64-bit integer MPI version or offer a pre-installed 64-bit
integer configuration:

-  Intel MPI (4.0.1 and higher)
-  Open MPI (please check below for more details)


Open MPI with 64-bit integers
-----------------------------

Please check the integer type of your currently available Open MPI
installation (perhaps it can do 64-bit integers already).

For this, type::

  ompi_info -a | grep 'Fort integer size'

if the output is::

  Fort integer size: 8 

then you have a suitable 64-bit Open MPI installation which uses 64-bit
integers by default.

If not and you decided that you want to build one yourself (not
difficult), then download the latest stable source from
`http://www.open-mpi.org/ <http://www.open-mpi.org/>`_ and extract it.

Then enter the directory and configure Open MPI (edit prefix path).

For Intel compilers use::

  ./configure --prefix=/path/to/openmpi \
              CXX=icpc                  \
              CC=icc                    \
              F77=ifort                 \
              FC=ifort                  \
              FFLAGS=-i8                \
              FCFLAGS=-i8 

For GNU compilers use::

  ./configure --prefix=/path/to/openmpi          \
              CXX=g++                            \
              CC=gcc                             \
              F77=gfortran                       \
              FC=gfortran                        \
              FFLAGS="-m64 -fdefault-integer-8"  \
              FCFLAGS="-m64 -fdefault-integer-8" \
              CFLAGS=-m64                        \
              CXXFLAGS=-m64

Once the configure step is finished (lot of output, at the end there
should be no error), build the library::

  make all install

Finally export the path locations (adapt paths) in .bashrc (or in the
corresponding file if you don't use bash)::

  export PATH=$PATH':/path/to/openmpi/intel/bin'
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/path/to/openmpi/intel/lib'

and verify the integer type::

  ompi_info -a | grep 'Fort integer size'

if the output is::

  Fort integer size: 8

Then you are set. If it is different, then something went wrong. Please
double-check and consult Open MPI documentation.
