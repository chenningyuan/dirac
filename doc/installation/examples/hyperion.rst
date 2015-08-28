:orphan:
 

Installation on `Hyperion <http://www.calmip.cict.fr/spip/spip.php?rubrique90>`_ (Toulouse)
===========================================================================================

We will use SGI Message Passing Toolkit (MPT) parallelization.


Configuration
-------------

When loggin in on the machine (Sep 2013) the available compilers and libraries are::

   $ module list
   Currently Loaded Modulefiles:
     1) default-intel-cc    3) default-intel-mkl
     2) default-intel-fc    4) default-mpt-sgi

Using the command ``module show`` we find that this corresponds to version 11.1 of the Intel compilers and MKL library 
as well as version 1.26 of SGI MPT. We go for some more recent versions and set (to see available modules do ``module avail``) ::

   $ module purge
   $ module load intel-cc-12/12.0 intel-fc-12/12.0  intel-mkl/12.0.4 sgi-mpt/2.0.4

In order for cmake to find the MKL library we also have to set::

   $ export MATH_ROOT='/divers/intel/Compiler/12.0.4/mkl'

For maximum performance we will add hardware-specific optimization flags::

  diff --git a/cmake/compilers/CFlags.cmake b/cmake/compilers/CFlags.cmake
  index f9ea370..26b74ea 100644
  --- a/cmake/compilers/CFlags.cmake
  +++ b/cmake/compilers/CFlags.cmake
  @@ -14,7 +14,7 @@ endif()
   if(CMAKE_C_COMPILER_ID MATCHES Intel)
       set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
       set(CMAKE_C_FLAGS_DEBUG   "-O0")
  -    set(CMAKE_C_FLAGS_RELEASE "-O2")
  +    set(CMAKE_C_FLAGS_RELEASE "-O2 -xSSE4.2")
       set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
   endif()

  diff --git a/cmake/compilers/FortranFlags.cmake b/cmake/compilers/FortranFlags.cmake
  index a2d2436..08a5721 100644
  --- a/cmake/compilers/FortranFlags.cmake
  +++ b/cmake/compilers/FortranFlags.cmake
  @@ -40,7 +40,7 @@ endif()
   if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
       set(CMAKE_Fortran_FLAGS         "-w -assume byterecl -DVAR_IFORT -g -traceback")
       set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
  -    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
  +    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xSSE4.2")
       if(ENABLE_STATIC_LINKING)
           set(CMAKE_Fortran_FLAGS
               "${CMAKE_Fortran_FLAGS} -static-libgcc -static-intel"

Finally we configure and build the binary::

  $ ./setup.py --fc=ifort --cc=icc --cxx=icpc --int64 --sgi-mpt
  $ cd build
  $ make


Setting up pam
--------------

In your top directory you may configure the pam run script by adding a file ``.diracrc`` with contents::

  --scratch=/tmpdir
  --machfile=$PBS_NODEFILE
  --noarch
  --debugger=/usr/bin/gdb
  --basis="$HOME/Dirac/basis:$HOME/Dirac/basis_dalton"
  --profiler=/usr/bin/gprof


Example run script
------------------

Here is a simple example of a run script::

  #PBS -S /bin/bash 
  #PBS -N fullHe
  #PBS -l select=1:ncpus=8:mpiprocs=4
  #PBS -l walltime=01:00:00
  #PBS -M trond.saue@irsamc.ups-tlse.fr -m ae 

  cd $PBS_O_WORKDIR
  nproc=`cat $PBS_NODEFILE | wc -l`
  inp=HF
  mol=H2O
  pam --mpi=$nproc --inp=${inp} --mol=${mol}
  exit 0

The meta-launcher paul
----------------------

Use the meta-launcher ``paul`` as a good starting point::

  $ /users/p1224/bast/bin/paul --help
