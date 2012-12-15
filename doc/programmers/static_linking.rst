
Static linking
==============
Theoretically, it's easy to tell the "setup" script to prepare the dirac.x executable statically linked (with the "--static" flag, 
thanks to the intervined python scripts and CMake buildup system.

However, neither all buildups are ending up happily with static dirac.x, nor - in cases of succesfull compilation - 
resulting static executables are performing correctly. 

Therefore we show some example configurations giving correct, partially correct, totally wrong and non-compilable executables.
We hope that this info shall help developers to identify possible cracks in the generation of their own static executables.

Developers and users are invited to extend this collection of experiences of static buildups on various platforms. 

Correct static executables
--------------------------
* ilias, 194.160.135.47:
  Linux-2.6.30-1-amd64, ownmath, i32lp64, GNU Fortran/gcc (Debian 4.4.5-10) 4.4.5; -static -fpic, nooptim
  Linux-2.6.30-1-amd64, ownmath, i32lp64, GNU Fortran/gcc (Debian 4.6.2-9); -static -fpic, nooptim

* ilias, frpd2.utcpd.sk:
  Linux-2.6.32-5-amd64; ownmath; i32lp64, GNU Fortran/gcc (Debian 4.6.1-4); -static -fpic, nooptim

* ilias, grid3.ui.savba.sk:
  Linux-2.6.18-274.3.1.el5; ownath; x86_64, ilp64, GNU Fortran (GCC) 4.4.4 20100726 (Red Hat 4.4.4-13)';  -g -static -fpic; opt
  Note: this gfortran shows error in test program for NAMELIST reading, however, DIRAC runs OK...
  NB: gfortran44 gives crashes for these OpenRSP tests: openrsp_pv,openrsp_jones,openrsp_cme.
  For "pure" runs you have to have higher version of GNU compilers.

* ilias, miro_ilias_desktop (home PC Desktop, with most recent Ubuntu) 
  a) Linux-3.0.0-14-generic; x86_64; i32lp64; GNU Fortran/gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1; 
     /usr/lib/libblas.a+/usr/lib/liblapack.a; -g -static -fpic -O0
     (NB: as superuser remove xerbla.o from /usr/lib/liblapack.a -mixing with /usr/lib/liblas.a)

  b) Linux-3.0.0-14-generic;x86_64; i32lp64; GNU Fortran/gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1; 
     ATLAS-static; -g -static -fpic -O0 (fixed bug below)

* ilias, opteron.tau.ac.il:
  ./setup  --debug  -D VERBOSE_OUTPUT=ON --build=build_static  --static --internal-math  --int64
  /opt/intel/fce/10.1.015/bin/ifort (ifort (IFORT) 10.1 20080312); /usr/bin/gcc - gcc (GCC) 4.1.2 20080704 (Red Hat 4.1.2-51)

* ilias, 194.160.34.207, Miro's old PC Desktop with updated Ubuntu
   Linux-3.0.0-16-generic; GNU Fortran (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1 -   both serial and parallel OpenMPI dirac.x are working 

Based on these achievements we may say that GNU compilers produce proper static executable
mostly not linked against GNU external static libraries (libblas-dev,liblapack-dev,libatlas).
Likewise the Intel compiler - ifort - gives correct static executable.

Partially wrong executables
---------------------------
One can state that static executables are partially wrong when large portion of tests is crashing. 
For instance, failures are happening in input readings of DFTINPUT, CI/RELCCSD.

* ilias, gbb4.grid.umb.sk:
  export MATH_ROOT='/home/chemia/ucitelia/milias/bin/intel_stuff/intel/mkl'
  setup --static -D VERBOSE_OUTPUT=ON --int64
  Linux-2.6.26-2-xen-amd64;ilp64,MKL ilp64 static math, ifort (IFORT) 12.0.2 20110112, gcc (Debian 4.3.2-1.1) 4.3.2

* ilias, pd.uniza.sk:
  setup --static -D VERBOSE_OUTPUT=ON  --internal-math (--debug)
  Linux-2.6.32-5-amd64, i32lp64, ownmath; GNU Fortran/gcc (Debian 4.4.5-8) 4.4.5

* ilias, grafix.fpv.umb.sk:
  a) setup --static -D VERBOSE_OUTPUT=ON  --debug
     Linux-2.6.32-28-generic; x86_64; i32lp64; GNU Fortran/Gcc (Ubuntu 4.4.3-4ubuntu5) 4.4.3; ATLAS-static
  b) setup -D VERBOSE_OUTPUT=ON --debug --static --int64 --internal-math
     Linux-2.6.32-28-generic; x86_64; ilp64; GNU Fortran/Gcc (Ubuntu 4.4.3-4ubuntu5) 4.4.3; ownmath

Note that the Gfortran error in NAMELIST reading is described  http://comments.gmane.org/gmane.comp.gcc.bugs/289633 .
Most probably this ill Gcc installation is probably causing mentioned failures of tests.

Wrong executables
------------------
Here we describe cases of static dirac.x executables crashing already at the start of the run.

* ilias, miro_ilias_desktop:
  Linux-3.0.0-14-generic;x86_64; i32lp64; GNU Fortran/gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1; 
  ATLAS-static; -g -static -fpic -O0, giving ::
   Program received signal 11 (SIGSEGV): Segmentation fault.
   Backtrace for this error:
    + function __restore_rt (0x15C4220)
      from file sigaction.c
   Segmentation fault
   (see also http://gcc.gnu.org/bugzilla/show_bug.cgi?format=multiple&id=42477)

* ilias, frpd2.utcpd.sk:
  Linux-2.6.32-5-amd64, ATLAS BLAS+LAPACK; i32lp64, GNU Fortran/gcc (Debian 4.6.1-4) 4.6.1; opt

* ilias, 194.160.135.47: 
  Linux-2.6.30-1-amd64, i32lp64, GNU Fortran/gcc (Debian 4.4.5-10) 4.4.5; -static -fpic, nooptim
  does not compile with static external /usr/lib/libblas.a+/usr/lib/liblapack.a, giving :: 
     /usr/lib/liblapack.a(xerbla.o):function xerbla_: error: undefined reference to '_gfortran_transfer_character_write'
     /usr/lib/liblapack.a(xerbla.o):function xerbla_: error: undefined reference to '_gfortran_transfer_integer_write'
and :: 
     /usr/bin/ld: error: /usr/lib/libblas.a(xerbla.o): multiple definition of 'xerbla\_'
     /usr/bin/ld: /usr/lib/liblapack.a(xerbla.o): previous definition here

The upgrade is compilable, but does not work with static ATLAS (as frpd2):
Linux-2.6.30-1-amd64, ATLAS BLAS+LAPACK; i32lp64, GNU Fortran/gcc (Debian 4.6.2-9); -static -fpic, nooptim

Executables not compilable
----------------------------
We list cases when one can not get static dirac.x executable based on some linking error.  

Linking Fortran executable dirac.x gives for parallel, OpenMPI, for example ::
  /home/ilias/bin/openmpi_ilp64/lib/libopen-pal.a(dlopen.o): In function `vm_open':
  dlopen.c:(.text+0x148): warning: Using 'dlopen' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
  /home/ilias/bin/openmpi_ilp64/lib/libopen-rte.a(plm_rsh_module.o): In function `setup_launch':
  plm_rsh_module.c:(.text+0x821): warning: Using 'getpwuid' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
  /home/ilias/bin/openmpi_ilp64/lib/libmpi.a(btl_tcp_component.o): In function `mca_btl_tcp_component_create_listen':
  btl_tcp_component.c:(.text+0x104): warning: Using 'getaddrinfo' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
  [100%] Built target dirac.x

* ilias, opteron.tau.ac.il
   export MATH_ROOT='/opt/intel/mkl/10.0.3.020/lib/em64t' , static buildup does not work with MKL libs, giving ::
     /opt/intel/fce/10.1.015/bin/ifort -static -Wl,-E -w -assume byterecl -DVART -g -traceback -static-libgcc -static-intel -i8 -O0 CMakeFiles/dirac.x.dir/main/main.F90.o  -o dirac.x lib/libdirac.a lib/libxcfun.a  -Wl,--start-g   /opt/intel/mkl/10.0.3.020/lib/em64t/libmkl_lapack.a    /opt/intel/mkl/10.0.3.020/lib/em64t/libmkl_core.a   /opt/intel/mkl/10.0.3.020/lib/em64t/libmkls_ilp64.a  /opt/intel/mkl/10.0.3.020/lib/em64t/libmkl_sequential.a /opt/intel/mkl/10.0.3.020/lib/em64t/libmkl_em64t.a   /opt/intel/mkl/10.0.3.020/lib/emibguide.a  /usr/lib64/libpthread.a /usr/lib64/libm.a -Wl,--end-group
     ld: cannot find libmkl_intel_lp64.a

