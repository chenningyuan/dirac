Intel MKL library
==================

The Intel MKL library can be linked in a sequential or in a multithreaded (default) mode; the latter is
advantegous for multi-CPU architectures.

In order to profir from the paralelization of MKL, the user has to provide environment variables.

Examples:                                       

- Safe-play default for processor unit with N cores 

::

 export MKL_NUM_THREADS="N"
 export MKL_DYNAMIC="FALSE"
 export OMP_NUM_THREADS="1" 

- Allow for 4 threads per core in MKL/BLAS routines

::

  export MKL_NUM_THREADS="4"
  export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=4"
  export OMP_NUM_THREADS="1"
  export MKL_DYNAMIC="FALSE"
  export OMP_DYNAMIC="FALSE" 

Recommended settings for calling Intel MKL routines from multi-threaded applications
is on Intel web page, 
http://software.intel.com/en-us/articles/recommended-settings-for-calling-intelr-mkl-routines-from-multi-threaded-applications/


