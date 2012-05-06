How to set good environment variables for the Intel MKL library
===============================================================

The Intel MKL library can be linked in a sequential or in a multithreaded (default) mode; 
the latter is advantegous for multi-CPU architectures.

In order to profit from the paralelization of MKL, 
the user should provide appropriate environment variables.

Examples of setting:
--------------------

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
is on the `Intel web page <http://software.intel.com/en-us/articles/recommended-settings-for-calling-intelr-mkl-routines-from-multi-threaded-applications/>`_.

Notes
-----

Sometimes you could get this error:

::
 OMP: Error #18: Setting environment variable "__KMP_REGISTERED_LIB_12973" failed:
 OMP: Hint: Seems application required too much memory.

This happens with statically linked multithreaded MKL library, 
see this thread `<http://software.intel.com/en-us/forums/showthread.php?t=104947>`_.

