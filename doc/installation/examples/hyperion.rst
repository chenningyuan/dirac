

Installation on Hyperion (Toulouse)
===================================

We will use SGI MPT parallelization.


Configuration
-------------

First we load some environment variables::

  $ . /etc/profile.d/modules.sh

  $ module load python

  $ export MATH_ROOT='/divers/intel/Compiler/12.0.4/mkl'
  $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/divers/intel/Compiler/12.0.4/mkl/lib/intel64'

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

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --int64 --sgi-mpt -D MPI_INCLUDE_PATH=/divers/sgi_mpt/opt/sgi/mpt/mpt-1.26/include
  $ make


Example run script
------------------

Use the meta-launcher ``paul`` as a good starting point::

  $ /users/p1224/bast/bin/paul --help
