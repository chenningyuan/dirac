Installing static dirac.x executable
====================================

Theoretically, it's easy to install the static dirac.x executable with the CMake buildup system.
However, not all static executables are performing correctly. Therefore we show some configurations giving
correct and wrong executables.

Correct static executables
----------------------------
ilias@194.160.135.47:
Linux-2.6.30-1-amd64, ownmath, i32lp64, GNU Fortran/gcc (Debian 4.4.5-10) 4.4.5; -static -fpic, nooptim
Linux-2.6.30-1-amd64, ownmath, i32lp64, GNU Fortran/gcc (Debian 4.6.2-9); -static -fpic, nooptim

ilias@frpd2.utcpd.sk:
Linux-2.6.32-5-amd64; ownmath; i32lp64, GNU Fortran/gcc (Debian 4.6.1-4); -static -fpic, nooptim

Based on these results we may say that GNU compilers produce proper static executable
NOT LINKED against GNU external static libraries (libblas-dev,liblapack-dev,libatlas).

Partially wrong executable
----------------------------
- crashing lot of tests

milias@gbb4.grid.umb.sk:
export MATH_ROOT='/home/chemia/ucitelia/milias/bin/intel_stuff/intel/mkl'
./setup --static -D VERBOSE_OUTPUT=ON  --int64
Linux-2.6.26-2-xen-amd64;ilp64,MKL ilp64 static math, ifort (IFORT) 12.0.2 20110112, gcc (Debian 4.3.2-1.1) 4.3.2
(failures: DFTINPUT, CI/RELCCSD input readings)

Wrong executables
------------------
- crashing already from beginning 

ilias@frpd2.utcpd.sk:
Linux-2.6.32-5-amd64, ATLAS BLAS+LAPACK; i32lp64, GNU Fortran/gcc (Debian 4.6.1-4) 4.6.1; opt

ilias@194.160.135.47:
Linux-2.6.30-1-amd64, i32lp64, GNU Fortran/gcc (Debian 4.4.5-10) 4.4.5; -static -fpic, nooptim
- does not compile with static external /usr/lib/libblas.a+liblapack.a 
/usr/lib/liblapack.a(xerbla.o):function xerbla_: error: undefined reference to '_gfortran_transfer_character_write'
/usr/lib/liblapack.a(xerbla.o):function xerbla_: error: undefined reference to '_gfortran_transfer_integer_write'
Likewise:
/usr/bin/ld: error: /usr/lib/libblas.a(xerbla.o): multiple definition of 'xerbla_'
/usr/bin/ld: /usr/lib/liblapack.a(xerbla.o): previous definition here

The upgrade is compilable, but does not work with static ATLAS (as frpd2):
Linux-2.6.30-1-amd64, ATLAS BLAS+LAPACK; i32lp64, GNU Fortran/gcc (Debian 4.6.2-9); -static -fpic, nooptim

One can summarize that linking to external static libraries from GNU does not give proper static executable.

