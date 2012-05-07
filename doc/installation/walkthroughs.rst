

Installation walkthroughs on some example machines
==================================================


Hyperion (2012-02-13)
---------------------

First set the correct environment. For this unload all possibly
conflicting modules:

::

    module unload default-intel-cc default-intel-fc default-intel-mkl default-mpt-sgi

Then load:

::

    module load intel-cc-12/12.0.4 intel-fc-12/12.0.4 intel-mkl/12.0.4

Also hook up to own 64bit integer OpenMPI build (ask Radovan for correct
path):

::

    export PATH=$PATH':/path/to/openmpi/bin'
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/path/to/openmpi/lib'

::

    export MATH_ROOT='/divers/intel/Compiler/12.0.4/mkl'
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/divers/intel/Compiler/12.0.4/mkl/lib/intel64'

Now setup DIRAC:

::

    ./setup --fc=mpif90 --cc=mpicc --cxx=icpc --int64 -D ENABLE_LARGE_MOLECULES=1

Finally build the code:

::

    cd build
    make -j16

