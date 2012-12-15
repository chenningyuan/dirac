

.. _linking_math:
Linking math libraries
======================

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

Finally set::

  $ export MATH_ROOT=/usr/lib64/atlas

Now DIRAC will correctly detect them.


Linking to Atlas libraries on Linux Mint
----------------------------------------

DIRAC may not be able to detect Atlas BLAS/Lapack
libraries on Linux Mint because of nonstandard (?) suffix::

After having installed them::

  $ sudo apt-get install liblapack3 libatlas3-base

the libraries were not detected by DIRAC. I had to do this::

  $ cd /usr/lib/atlas-base
  $ su
  $ for file in *3; do ln -s $file ${file%.3}; done

and this::

  $ cd /usr/lib
  $ sudo ln -s liblapack.so.3 liblapack.so

With this the automatic detection of math libraries worked.
