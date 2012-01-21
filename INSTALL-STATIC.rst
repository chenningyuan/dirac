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

ilias@grid3.ui.savba.sk:
Linux-2.6.18-274.3.1.el5; ownath; x86_64, ilp64, GNU Fortran (GCC) 4.4.4 20100726 (Red Hat 4.4.4-13)';  -g -static -fpic; opt
Note: this gfortran shows error in test program for NAMELIST reading, however, DIRAC runs OK...
NB: gfortran44 gives crashes for these OpenRSP tests: openrsp_pv,openrsp_jones,openrsp_cme.
For "pure" runs you have to havehigher version of GNU compilers.

ilias@miro_ilias_desktop:
Linux-3.0.0-14-generic;x86_64; i32lp64; GNU Fortran/gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1; 
/usr/lib/libblas.a+/usr/lib/liblapack.a; -g -static -fpic -O0
(NB: as superuser remove xerbla.o from /usr/lib/liblapack.a -mixing with /usr/lib/liblas.a)

ilias@miro_ilias_desktop:
Linux-3.0.0-14-generic;x86_64; i32lp64; GNU Fortran/gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1; 
ATLAS-static; -g -static -fpic -O0 (fixed bug below)

Based on these results we may say that GNU compilers produce proper static executable
mostly not linked against GNU external static libraries (libblas-dev,liblapack-dev,libatlas).

Partially wrong executable
----------------------------
- crashing substantial portion of tests

(failures: DFTINPUT, CI/RELCCSD input readings)

milias@gbb4.grid.umb.sk:
export MATH_ROOT='/home/chemia/ucitelia/milias/bin/intel_stuff/intel/mkl'
./setup --static -D VERBOSE_OUTPUT=ON  --int64
Linux-2.6.26-2-xen-amd64;ilp64,MKL ilp64 static math, ifort (IFORT) 12.0.2 20110112, gcc (Debian 4.3.2-1.1) 4.3.2

ilias@pd.uniza.sk:
./setup --static -D VERBOSE_OUTPUT=ON  --internal-math (--debug)
Linux-2.6.32-5-amd64, i32lp64, ownmath; GNU Fortran/gcc (Debian 4.4.5-8) 4.4.5

milias@grafix.fpv.umb.sk:
./setup --static -D VERBOSE_OUTPUT=ON  --debug
Linux-2.6.32-28-generic; x86_64; i32lp64; GNU Fortran/Gcc (Ubuntu 4.4.3-4ubuntu5) 4.4.3; ATLAS-static
./setup -D VERBOSE_OUTPUT=ON --debug --static --int64 --internal-math
Linux-2.6.32-28-generic; x86_64; ilp64; GNU Fortran/Gcc (Ubuntu 4.4.3-4ubuntu5) 4.4.3; ownmath

Note: Gfortran error in NAMELIST reading is decribed here:  http://comments.gmane.org/gmane.comp.gcc.bugs/289633 .
Such 'ill' Gcc nstallation is probably causing mentioned failures of tests.

Wrong executables
------------------
- crashing already from the beginning 

ilias@miro_ilias_desktop:
Linux-3.0.0-14-generic;x86_64; i32lp64; GNU Fortran/gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1; 
ATLAS-static; -g -static -fpic -O0
Program received signal 11 (SIGSEGV): Segmentation fault.

Backtrace for this error:
  + function __restore_rt (0x15C4220)
    from file sigaction.c
Segmentation fault
(see http://gcc.gnu.org/bugzilla/show_bug.cgi?format=multiple&id=42477)


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


* Static parallel dirac.x - not working at all

 Linking Fortran executable dirac.x gives:
/home/ilias/bin/openmpi_ilp64/lib/libopen-pal.a(dlopen.o): In function `vm_open':
dlopen.c:(.text+0x148): warning: Using 'dlopen' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
/home/ilias/bin/openmpi_ilp64/lib/libopen-rte.a(plm_rsh_module.o): In function `setup_launch':
plm_rsh_module.c:(.text+0x821): warning: Using 'getpwuid' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
/home/ilias/bin/openmpi_ilp64/lib/libmpi.a(btl_tcp_component.o): In function `mca_btl_tcp_component_create_listen':
btl_tcp_component.c:(.text+0x104): warning: Using 'getaddrinfo' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
[100%] Built target dirac.x
