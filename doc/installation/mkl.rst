
How to set good environment variables for the Intel MKL library
===============================================================

The Intel MKL library can be linked in a sequential or in a multi-threaded mode
(default); the latter is advantageous for multi-CPU architectures.

In order to benefit from the parallelization of MKL,
the user should provide appropriate environment variables.


Examples
--------

Safe-play default for processor unit with N cores::

  export MKL_NUM_THREADS="N"
  export MKL_DYNAMIC="FALSE"
  export OMP_NUM_THREADS=1

Allow for 4 threads per core in MKL/BLAS routines::

  export MKL_NUM_THREADS=4
  export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=4"
  export OMP_NUM_THREADS=1
  export MKL_DYNAMIC="FALSE"
  export OMP_DYNAMIC="FALSE"

Note that is important to properly balance MPI and MKL threads. 
If all cores are already taken by the MPI, the the "unwanted" feature of the threaded-MKL  
is (sometimes) significant slow down of the code's execution run.

As an example of proper splitting of CPUs between MPI and MKL, 
assume a 16 core node where 8 cores are assigned to MPI and the rest to MKL threads::

  export MKL_NUM_THREADS="2"
  export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=2"
  export OMP_NUM_THREADS="1"
  export MKL_DYNAMIC="FALSE"
  export OMP_DYNAMIC="FALSE"

If you would like to run a sequential job on the same node it would then read::

  export MKL_NUM_THREADS="16"
  export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=16"
  export OMP_NUM_THREADS="1"
  export MKL_DYNAMIC="FALSE"
  export OMP_DYNAMIC="FALSE"

You can find recommended settings for calling Intel MKL routines from multi-threaded applications
on the `Intel web page <http://software.intel.com/en-us/articles/recommended-settings-for-calling-intelr-mkl-routines-from-multi-threaded-applications/>`_.


Troubleshooting
---------------

Sometimes you could get this error::

  OMP: Error #18: Setting environment variable "__KMP_REGISTERED_LIB_12973" failed:
  OMP: Hint: Seems application required too much memory.

This can happen with statically linked multithreaded MKL library,
see `here <http://software.intel.com/en-us/forums/showthread.php?t=104947>`_.
