How to build MPI libraries for 64-bit integers
==============================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Before you continue please read the
`FAQ </dirac/index.php/FAQ#I_am_about_to_install_DIRAC_on_a_64-bit_machine.2C_should_I_compile_using_64-bit_integers_or_32-bit_integers.3F>`_
and verify whether you really need 64-bit integers.

It is extra work and perhaps not needed.

We are currently aware of two MPI distributions that can either be
installed as 64-bit integer MPI version or offer a pre-installed 64-bit
integer configuration:

-  Intel MPI (v. 4.0.1 and higher)
-  Open MPI (please check below for more details)

Open MPI with 64-bit integers
-----------------------------

Please check the integer type of your currently available Open MPI
installation (perhaps it can do 64-bit integers already).

For this, type

::

    ompi_info -a | grep 'Fort integer size'

if the output is

::

    Fort integer size: 8 

then you have a suitable 64-bit Open MPI installation which uses 64-bit
integers by default.

If not and you decided that you want to build one yourself (not
difficult), then download the latest stable source from
`http://www.open-mpi.org/ <http://www.open-mpi.org/>`_ and extract it.

Then enter the directory and configure Open MPI (edit prefix path).

For Intel compilers use:

::

    ./configure --prefix=/path/to/openmpi CXX=icpc CC=icc F77=ifort FC=ifort FFLAGS=-i8 FCFLAGS=-i8 

For GNU compilers use:

::

    ./configure --prefix=/path/to/openmpi CXX=g++ CC=gcc F77=gfortran FC=gfortran \
                                          FFLAGS="-m64 -fdefault-integer-8"       \
                                          FCFLAGS="-m64 -fdefault-integer-8"      \
                                          CFLAGS=-m64                             \
                                          CXXFLAGS=-m64

If you want to add possibility for debugging with **valgrind**, please
add these extra flags to configure:

::

     --enable-debug --enable-memchecker --with-valgrind(=/path/to/valgrind)

You get, however, slower performace, for more visit
`http://www.open-mpi.org/faq/?category=debugging#memchecker\_how <http://www.open-mpi.org/faq/?category=debugging#memchecker_how>`_.

Once the configure step is finished (lot of output, at the end there
should be no error), build the library:

::

    make all install

Finally export the path locations (adapt paths) in .bashrc (or in the
corresponding file if you don't use bash):

::

    export PATH=$PATH':/path/to/openmpi/intel/bin'
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/path/to/openmpi/intel/lib'

and verify the integer type:

::

    ompi_info -a | grep 'Fort integer size'

if the output is

::

    Fort integer size: 8

Then you are set. If it is different, then something went wrong. Please
double-check and consult Open MPI documentation.

