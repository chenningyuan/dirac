   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:16:28 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_fdekern.inp 
   PID         : 12744
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av2z.12744


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040b8e28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c21a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 628 bytes at 0000000103c21788 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.12744 on main (master) node fpcmtaspg:
total 71928
-rw-r--r--  1 gomes  wheel    271436 Mar 26 20:16 AOMOMAT
-rw-r--r--  1 gomes  wheel    993456 Mar 26 20:16 AOPROPER
-rw-r--r--  1 gomes  wheel     23336 Mar 26 20:16 DFCMOS
-rw-r--r--  1 gomes  wheel    136102 Mar 26 20:16 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 20:16 DFCYCL
-rw-r--r--  1 gomes  wheel    788858 Mar 26 20:16 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:16 DFDIIS
-rw-r--r--  1 gomes  wheel    163296 Mar 26 20:16 DFEVEC
-rw-r--r--  1 gomes  wheel    788874 Mar 26 20:16 DFFCK1
-rw-r--r--  1 gomes  wheel    788874 Mar 26 20:16 DFFCK2
-rw-r--r--  1 gomes  wheel    209952 Mar 26 20:16 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 20:16 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:16 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:16 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:16 FRZDN1
-rw-r--r--  1 gomes  wheel     23456 Mar 26 20:16 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 20:16 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:16 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    756864 Mar 26 20:18 PAMBOE
-rw-r--r--  1 gomes  wheel     23336 Mar 26 20:16 PAMFCK
-rw-r--r--  1 gomes  wheel       228 Mar 26 20:18 PAMRSP
-rw-r--r--  1 gomes  wheel     94452 Mar 26 20:18 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:16 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:16 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 35M	/tmp/gomes/gomes/ar.av2z.12744
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.12744 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:18:59 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:18:59 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_nofde.inp 
   PID         : 12854
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av2z.12854


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 000000010409a628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104078228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c214f8 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 628 bytes at 0000000103c21748 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21628 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c215b8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000103c21568 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c216b8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 628 bytes at 0000000103c21c88 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c219e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010407a028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 0000000104079628 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104096828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 1256 bytes at 000000010407be28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104094a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104098628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104099628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 000000012abdb028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010409e628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010409d628 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 000000010409c628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010409b628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010409f628 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.12854 on main (master) node fpcmtaspg:
total 71552
-rw-r--r--  1 gomes  wheel    271436 Mar 26 20:19 AOMOMAT
-rw-r--r--  1 gomes  wheel    894184 Mar 26 20:19 AOPROPER
-rw-r--r--  1 gomes  wheel     23336 Mar 26 20:19 DFCMOS
-rw-r--r--  1 gomes  wheel    136102 Mar 26 20:19 DFCOEF
-rw-r--r--  1 gomes  wheel       920 Mar 26 20:19 DFCYCL
-rw-r--r--  1 gomes  wheel    788858 Mar 26 20:19 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:19 DFDIIS
-rw-r--r--  1 gomes  wheel    139968 Mar 26 20:19 DFEVEC
-rw-r--r--  1 gomes  wheel    788874 Mar 26 20:19 DFFCK1
-rw-r--r--  1 gomes  wheel    788874 Mar 26 20:19 DFFCK2
-rw-r--r--  1 gomes  wheel    186624 Mar 26 20:19 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 20:18 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:18 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:18 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:18 FRZDN1
-rw-r--r--  1 gomes  wheel     23456 Mar 26 20:19 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 20:19 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:18 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    705024 Mar 26 20:20 PAMBOE
-rw-r--r--  1 gomes  wheel     23336 Mar 26 20:19 PAMFCK
-rw-r--r--  1 gomes  wheel       156 Mar 26 20:20 PAMRSP
-rw-r--r--  1 gomes  wheel     94452 Mar 26 20:20 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:19 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:19 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 35M	/tmp/gomes/gomes/ar.av2z.12854
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.12854 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:20:19 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:20:19 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_nokern.inp 
   PID         : 12955
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av2z.12955


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c21a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 628 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c21788 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 00000001288b2028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.12955 on main (master) node fpcmtaspg:
total 72048
-rw-r--r--  1 gomes  wheel    271436 Mar 26 20:20 AOMOMAT
-rw-r--r--  1 gomes  wheel    993456 Mar 26 20:20 AOPROPER
-rw-r--r--  1 gomes  wheel     23336 Mar 26 20:20 DFCMOS
-rw-r--r--  1 gomes  wheel    136102 Mar 26 20:20 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 20:20 DFCYCL
-rw-r--r--  1 gomes  wheel    788858 Mar 26 20:20 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:20 DFDIIS
-rw-r--r--  1 gomes  wheel    163296 Mar 26 20:20 DFEVEC
-rw-r--r--  1 gomes  wheel    788874 Mar 26 20:20 DFFCK1
-rw-r--r--  1 gomes  wheel    788874 Mar 26 20:20 DFFCK2
-rw-r--r--  1 gomes  wheel    209952 Mar 26 20:20 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 20:20 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:20 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:20 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:20 FRZDN1
-rw-r--r--  1 gomes  wheel     23456 Mar 26 20:20 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 20:20 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:20 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    819072 Mar 26 20:22 PAMBOE
-rw-r--r--  1 gomes  wheel     23336 Mar 26 20:20 PAMFCK
-rw-r--r--  1 gomes  wheel       252 Mar 26 20:22 PAMRSP
-rw-r--r--  1 gomes  wheel     94452 Mar 26 20:22 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:20 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:20 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 35M	/tmp/gomes/gomes/ar.av2z.12955
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.12955 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:22:32 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:22:32 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_fdekern.inp 
   PID         : 13061
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v2z.13061


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0e28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104095a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c21708 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c21298 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 496 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 0000000104090a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.13061 on main (master) node fpcmtaspg:
total 67720
-rw-r--r--  1 gomes  wheel    142988 Mar 26 20:22 AOMOMAT
-rw-r--r--  1 gomes  wheel    621080 Mar 26 20:22 AOPROPER
-rw-r--r--  1 gomes  wheel     10376 Mar 26 20:22 DFCMOS
-rw-r--r--  1 gomes  wheel     71770 Mar 26 20:22 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 20:22 DFCYCL
-rw-r--r--  1 gomes  wheel    492122 Mar 26 20:22 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:22 DFDIIS
-rw-r--r--  1 gomes  wheel     82944 Mar 26 20:22 DFEVEC
-rw-r--r--  1 gomes  wheel    492138 Mar 26 20:22 DFFCK1
-rw-r--r--  1 gomes  wheel    492138 Mar 26 20:22 DFFCK2
-rw-r--r--  1 gomes  wheel     93312 Mar 26 20:22 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 20:22 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:22 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:22 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:22 FRZDN1
-rw-r--r--  1 gomes  wheel     10496 Mar 26 20:22 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       529 Mar 26 20:22 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 20:22 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    336960 Mar 26 20:24 PAMBOE
-rw-r--r--  1 gomes  wheel     10376 Mar 26 20:22 PAMFCK
-rw-r--r--  1 gomes  wheel       228 Mar 26 20:24 PAMRSP
-rw-r--r--  1 gomes  wheel     47796 Mar 26 20:24 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:22 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:22 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 33M	/tmp/gomes/gomes/ar.v2z.13061
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.13061 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:24:33 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:24:33 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_nofde.inp 
   PID         : 13161
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v2z.13161


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000104083028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104079c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c214f8 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 496 bytes at 0000000103c21748 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21628 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c215b8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000103c21568 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c216b8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 496 bytes at 0000000103c21b88 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c21968 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 0000000104079228 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010407f228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 992 bytes at 000000010407ac28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104081028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104082028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 000000012384c028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 00000001040ece28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104087028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104086028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104085028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104084028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000104088028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040b9e28 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.13161 on main (master) node fpcmtaspg:
total 67744
-rw-r--r--  1 gomes  wheel    142988 Mar 26 20:24 AOMOMAT
-rw-r--r--  1 gomes  wheel    559032 Mar 26 20:24 AOPROPER
-rw-r--r--  1 gomes  wheel     10376 Mar 26 20:24 DFCMOS
-rw-r--r--  1 gomes  wheel     71770 Mar 26 20:24 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 20:24 DFCYCL
-rw-r--r--  1 gomes  wheel    492122 Mar 26 20:24 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:24 DFDIIS
-rw-r--r--  1 gomes  wheel     93312 Mar 26 20:24 DFEVEC
-rw-r--r--  1 gomes  wheel    492138 Mar 26 20:24 DFFCK1
-rw-r--r--  1 gomes  wheel    492138 Mar 26 20:24 DFFCK2
-rw-r--r--  1 gomes  wheel     93312 Mar 26 20:24 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 20:24 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:24 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:24 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:24 FRZDN1
-rw-r--r--  1 gomes  wheel     10496 Mar 26 20:24 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       529 Mar 26 20:24 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 20:24 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    404352 Mar 26 20:25 PAMBOE
-rw-r--r--  1 gomes  wheel     10376 Mar 26 20:24 PAMFCK
-rw-r--r--  1 gomes  wheel       132 Mar 26 20:25 PAMRSP
-rw-r--r--  1 gomes  wheel     47796 Mar 26 20:25 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:24 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:24 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 33M	/tmp/gomes/gomes/ar.v2z.13161
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.13161 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.

   Exit date: Fri Mar 26 20:25:31 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:25:31 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_nokern.inp 
   PID         : 13264
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v2z.13264


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c21708 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 496 bytes at 0000000103c21298 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 00000001288b2028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104095a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 992 bytes at 000000010406ea28 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.13264 on main (master) node fpcmtaspg:
total 67728
-rw-r--r--  1 gomes  wheel    142988 Mar 26 20:25 AOMOMAT
-rw-r--r--  1 gomes  wheel    621080 Mar 26 20:25 AOPROPER
-rw-r--r--  1 gomes  wheel     10376 Mar 26 20:25 DFCMOS
-rw-r--r--  1 gomes  wheel     71770 Mar 26 20:25 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 20:25 DFCYCL
-rw-r--r--  1 gomes  wheel    492122 Mar 26 20:25 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:25 DFDIIS
-rw-r--r--  1 gomes  wheel     82944 Mar 26 20:25 DFEVEC
-rw-r--r--  1 gomes  wheel    492138 Mar 26 20:25 DFFCK1
-rw-r--r--  1 gomes  wheel    492138 Mar 26 20:25 DFFCK2
-rw-r--r--  1 gomes  wheel     93312 Mar 26 20:25 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 20:25 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:25 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:25 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:25 FRZDN1
-rw-r--r--  1 gomes  wheel     10496 Mar 26 20:25 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       529 Mar 26 20:25 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 20:25 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    342144 Mar 26 20:27 PAMBOE
-rw-r--r--  1 gomes  wheel     10376 Mar 26 20:25 PAMFCK
-rw-r--r--  1 gomes  wheel       228 Mar 26 20:27 PAMRSP
-rw-r--r--  1 gomes  wheel     47796 Mar 26 20:27 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:25 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:25 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 33M	/tmp/gomes/gomes/ar.v2z.13264
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.13264 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:27:09 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:27:09 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_fdekern.inp 
   PID         : 13364
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av3z.13364


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1032 bytes at 000000010407fa28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1032 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.13364 on main (master) node fpcmtaspg:
total 89176
-rw-r--r--  1 gomes  wheel    825740 Mar 26 20:27 AOMOMAT
-rw-r--r--  1 gomes  wheel   2674600 Mar 26 20:27 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 20:27 DFCMOS
-rw-r--r--  1 gomes  wheel    413530 Mar 26 20:27 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 20:27 DFCYCL
-rw-r--r--  1 gomes  wheel   2130138 Mar 26 20:27 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:27 DFDIIS
-rw-r--r--  1 gomes  wheel    640000 Mar 26 20:27 DFEVEC
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 20:27 DFFCK1
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 20:27 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 20:27 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 20:27 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:27 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:27 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:27 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 20:27 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       745 Mar 26 20:27 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:27 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1771200 Mar 26 20:32 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 20:27 PAMFCK
-rw-r--r--  1 gomes  wheel       220 Mar 26 20:32 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 20:32 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:27 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:27 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 44M	/tmp/gomes/gomes/ar.av3z.13364
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.13364 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:32:26 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:32:26 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_nofde.inp 
   PID         : 13468
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av3z.13468


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000129083228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129080228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000128900108 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 640000 bytes at 000000012384c028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000012907a828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128900808 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000128900798 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000128900748 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128900898 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 2200 bytes at 0000000129075428 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128f42028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000012907e428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000012907c628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129081228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129082228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1032 bytes at 00000001290a9028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129087228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129086228 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000129085228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129084228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000129088228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 0000000129099628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012908d228 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.13468 on main (master) node fpcmtaspg:
total 88312
-rw-r--r--  1 gomes  wheel    825740 Mar 26 20:32 AOMOMAT
-rw-r--r--  1 gomes  wheel   2407264 Mar 26 20:32 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 20:32 DFCMOS
-rw-r--r--  1 gomes  wheel    413530 Mar 26 20:32 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 20:32 DFCYCL
-rw-r--r--  1 gomes  wheel   2130138 Mar 26 20:32 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:32 DFDIIS
-rw-r--r--  1 gomes  wheel    560000 Mar 26 20:32 DFEVEC
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 20:32 DFFCK1
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 20:32 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 20:32 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 20:32 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:32 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:32 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:32 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 20:32 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       745 Mar 26 20:32 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:32 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1676736 Mar 26 20:36 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 20:32 PAMFCK
-rw-r--r--  1 gomes  wheel       132 Mar 26 20:36 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 20:36 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:32 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:32 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 43M	/tmp/gomes/gomes/ar.av3z.13468
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.13468 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:36:15 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:36:15 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_nokern.inp 
   PID         : 13568
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av3z.13568


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128996028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2064 bytes at 0000000104080028 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.13568 on main (master) node fpcmtaspg:
total 89640
-rw-r--r--  1 gomes  wheel    825740 Mar 26 20:36 AOMOMAT
-rw-r--r--  1 gomes  wheel   2674600 Mar 26 20:36 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 20:36 DFCMOS
-rw-r--r--  1 gomes  wheel    413530 Mar 26 20:36 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 20:36 DFCYCL
-rw-r--r--  1 gomes  wheel   2130138 Mar 26 20:36 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:36 DFDIIS
-rw-r--r--  1 gomes  wheel    640000 Mar 26 20:36 DFEVEC
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 20:36 DFFCK1
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 20:36 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 20:36 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 20:36 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:36 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:36 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:36 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 20:36 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       745 Mar 26 20:36 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:36 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2007360 Mar 26 20:41 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 20:36 PAMFCK
-rw-r--r--  1 gomes  wheel       220 Mar 26 20:41 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 20:41 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:36 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:36 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 44M	/tmp/gomes/gomes/ar.av3z.13568
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.13568 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:41:22 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:41:22 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_fdekern.inp 
   PID         : 13671
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v3z.13671


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c21b48 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 776 bytes at 0000000103c21818 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1552 bytes at 000000010407ee28 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.13671 on main (master) node fpcmtaspg:
total 77240
-rw-r--r--  1 gomes  wheel    422284 Mar 26 20:41 AOMOMAT
-rw-r--r--  1 gomes  wheel   1514920 Mar 26 20:41 AOPROPER
-rw-r--r--  1 gomes  wheel     37000 Mar 26 20:41 DFCMOS
-rw-r--r--  1 gomes  wheel    211610 Mar 26 20:41 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 20:41 DFCYCL
-rw-r--r--  1 gomes  wheel   1204442 Mar 26 20:41 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:41 DFDIIS
-rw-r--r--  1 gomes  wheel    295936 Mar 26 20:41 DFEVEC
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 20:41 DFFCK1
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 20:41 DFFCK2
-rw-r--r--  1 gomes  wheel    332928 Mar 26 20:41 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 20:41 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:41 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:41 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:41 FRZDN1
-rw-r--r--  1 gomes  wheel     37120 Mar 26 20:41 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       660 Mar 26 20:41 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 20:41 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1137600 Mar 26 20:44 PAMBOE
-rw-r--r--  1 gomes  wheel     37000 Mar 26 20:41 PAMFCK
-rw-r--r--  1 gomes  wheel       236 Mar 26 20:44 PAMRSP
-rw-r--r--  1 gomes  wheel    130740 Mar 26 20:44 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:41 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:41 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 38M	/tmp/gomes/gomes/ar.v3z.13671
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.13671 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:44:11 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:44:11 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_nofde.inp 
   PID         : 13774
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v3z.13774


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000104064628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104061628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21088 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c000e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c21288 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21168 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c210f8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211f8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 776 bytes at 0000000103c218e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c215b8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405bc28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010400f628 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405da28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104062628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104063628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 000000012384c028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067628 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000104069628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 000000010407aa28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406e628 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.13774 on main (master) node fpcmtaspg:
total 76816
-rw-r--r--  1 gomes  wheel    422284 Mar 26 20:44 AOMOMAT
-rw-r--r--  1 gomes  wheel   1363552 Mar 26 20:44 AOPROPER
-rw-r--r--  1 gomes  wheel     37000 Mar 26 20:44 DFCMOS
-rw-r--r--  1 gomes  wheel    211610 Mar 26 20:44 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 20:44 DFCYCL
-rw-r--r--  1 gomes  wheel   1204442 Mar 26 20:44 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:44 DFDIIS
-rw-r--r--  1 gomes  wheel    258944 Mar 26 20:44 DFEVEC
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 20:44 DFFCK1
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 20:44 DFFCK2
-rw-r--r--  1 gomes  wheel    332928 Mar 26 20:44 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 20:44 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:44 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:44 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:44 FRZDN1
-rw-r--r--  1 gomes  wheel     37120 Mar 26 20:44 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       660 Mar 26 20:44 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 20:44 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1108800 Mar 26 20:45 PAMBOE
-rw-r--r--  1 gomes  wheel     37000 Mar 26 20:44 PAMFCK
-rw-r--r--  1 gomes  wheel       132 Mar 26 20:45 PAMRSP
-rw-r--r--  1 gomes  wheel    130740 Mar 26 20:45 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:44 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:44 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 38M	/tmp/gomes/gomes/ar.v3z.13774
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.13774 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:45:48 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:45:48 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_nokern.inp 
   PID         : 13876
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v3z.13876


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c21b48 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 776 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c21818 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128935028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.13876 on main (master) node fpcmtaspg:
total 77296
-rw-r--r--  1 gomes  wheel    422284 Mar 26 20:45 AOMOMAT
-rw-r--r--  1 gomes  wheel   1514920 Mar 26 20:45 AOPROPER
-rw-r--r--  1 gomes  wheel     37000 Mar 26 20:46 DFCMOS
-rw-r--r--  1 gomes  wheel    211610 Mar 26 20:46 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 20:46 DFCYCL
-rw-r--r--  1 gomes  wheel   1204442 Mar 26 20:46 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:46 DFDIIS
-rw-r--r--  1 gomes  wheel    295936 Mar 26 20:46 DFEVEC
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 20:45 DFFCK1
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 20:46 DFFCK2
-rw-r--r--  1 gomes  wheel    332928 Mar 26 20:46 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 20:45 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:45 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:45 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:45 FRZDN1
-rw-r--r--  1 gomes  wheel     37120 Mar 26 20:46 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       660 Mar 26 20:45 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 20:45 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1166400 Mar 26 20:48 PAMBOE
-rw-r--r--  1 gomes  wheel     37000 Mar 26 20:46 PAMFCK
-rw-r--r--  1 gomes  wheel       236 Mar 26 20:48 PAMRSP
-rw-r--r--  1 gomes  wheel    130740 Mar 26 20:48 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:45 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:45 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 38M	/tmp/gomes/gomes/ar.v3z.13876
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.13876 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 20:48:15 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 20:48:15 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av4z.mol lda_ll_fdekern.inp 
   PID         : 13978
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av4z.13978


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av4z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1684 bytes at 000000010407fe28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1684 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.av4z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av4z.13978 on main (master) node fpcmtaspg:
total 131984
-rw-r--r--  1 gomes  wheel   2263436 Mar 26 20:48 AOMOMAT
-rw-r--r--  1 gomes  wheel   7109096 Mar 26 20:48 AOPROPER
-rw-r--r--  1 gomes  wheel    225800 Mar 26 20:48 DFCMOS
-rw-r--r--  1 gomes  wheel   1132786 Mar 26 20:48 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 20:48 DFCYCL
-rw-r--r--  1 gomes  wheel   5671802 Mar 26 20:48 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 20:48 DFDIIS
-rw-r--r--  1 gomes  wheel   2032128 Mar 26 20:48 DFEVEC
-rw-r--r--  1 gomes  wheel   5671818 Mar 26 20:48 DFFCK1
-rw-r--r--  1 gomes  wheel   5671818 Mar 26 20:48 DFFCK2
-rw-r--r--  1 gomes  wheel   2032128 Mar 26 20:48 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 20:48 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 20:48 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 20:48 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 20:48 FRZDN1
-rw-r--r--  1 gomes  wheel    225920 Mar 26 20:48 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       897 Mar 26 20:48 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 20:48 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   3153600 Mar 26 21:00 PAMBOE
-rw-r--r--  1 gomes  wheel    225800 Mar 26 20:48 PAMFCK
-rw-r--r--  1 gomes  wheel       252 Mar 26 21:00 PAMRSP
-rw-r--r--  1 gomes  wheel    389940 Mar 26 21:00 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 20:48 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 20:48 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 64M	/tmp/gomes/gomes/ar.av4z.13978
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av4z.13978 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 21:00:59 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 21:01:00 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av4z.mol lda_ll_nofde.inp 
   PID         : 14089
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av4z.14089


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av4z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 000000010405b428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104002628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21058 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 7600 bytes at 0000000104000828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21188 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21118 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000103c210c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21218 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 0000000104003628 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 1684 bytes at 000000010409fe28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1684 bytes at 0000000104081a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010405f428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010405e428 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 000000010405d428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010405c428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000104060428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 0000000104071828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104066428 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 1684 bytes at 0000000104081228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010409c028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3368 bytes at 0000000104082228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010409a228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010409de28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.av4z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av4z.14089 on main (master) node fpcmtaspg:
total 129968
-rw-r--r--  1 gomes  wheel   2263436 Mar 26 21:01 AOMOMAT
-rw-r--r--  1 gomes  wheel   6398400 Mar 26 21:01 AOPROPER
-rw-r--r--  1 gomes  wheel    225800 Mar 26 21:01 DFCMOS
-rw-r--r--  1 gomes  wheel   1132786 Mar 26 21:01 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 21:01 DFCYCL
-rw-r--r--  1 gomes  wheel   5671802 Mar 26 21:01 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 21:01 DFDIIS
-rw-r--r--  1 gomes  wheel   1580544 Mar 26 21:01 DFEVEC
-rw-r--r--  1 gomes  wheel   5671818 Mar 26 21:01 DFFCK1
-rw-r--r--  1 gomes  wheel   5671818 Mar 26 21:01 DFFCK2
-rw-r--r--  1 gomes  wheel   2032128 Mar 26 21:01 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 21:01 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 21:01 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 21:01 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 21:01 FRZDN1
-rw-r--r--  1 gomes  wheel    225920 Mar 26 21:01 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       897 Mar 26 21:01 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 21:01 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   3283200 Mar 26 21:11 PAMBOE
-rw-r--r--  1 gomes  wheel    225800 Mar 26 21:01 PAMFCK
-rw-r--r--  1 gomes  wheel       180 Mar 26 21:11 PAMRSP
-rw-r--r--  1 gomes  wheel    389940 Mar 26 21:11 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 21:01 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 21:01 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 63M	/tmp/gomes/gomes/ar.av4z.14089
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av4z.14089 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 21:11:56 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 21:11:56 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av4z.mol lda_ll_nokern.inp 
   PID         : 14193
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av4z.14193


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av4z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128a4b028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3368 bytes at 0000000104080628 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.av4z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av4z.14193 on main (master) node fpcmtaspg:
total 132160
-rw-r--r--  1 gomes  wheel   2263436 Mar 26 21:12 AOMOMAT
-rw-r--r--  1 gomes  wheel   7109096 Mar 26 21:12 AOPROPER
-rw-r--r--  1 gomes  wheel    225800 Mar 26 21:12 DFCMOS
-rw-r--r--  1 gomes  wheel   1132786 Mar 26 21:12 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 21:12 DFCYCL
-rw-r--r--  1 gomes  wheel   5671802 Mar 26 21:12 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 21:12 DFDIIS
-rw-r--r--  1 gomes  wheel   2032128 Mar 26 21:12 DFEVEC
-rw-r--r--  1 gomes  wheel   5671818 Mar 26 21:12 DFFCK1
-rw-r--r--  1 gomes  wheel   5671818 Mar 26 21:12 DFFCK2
-rw-r--r--  1 gomes  wheel   2032128 Mar 26 21:12 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 21:11 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 21:11 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 21:11 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 21:11 FRZDN1
-rw-r--r--  1 gomes  wheel    225920 Mar 26 21:12 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       897 Mar 26 21:11 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 21:11 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   3240000 Mar 26 21:24 PAMBOE
-rw-r--r--  1 gomes  wheel    225800 Mar 26 21:12 PAMFCK
-rw-r--r--  1 gomes  wheel       252 Mar 26 21:24 PAMRSP
-rw-r--r--  1 gomes  wheel    389940 Mar 26 21:24 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 21:11 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 21:11 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 65M	/tmp/gomes/gomes/ar.av4z.14193
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av4z.14193 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 21:24:14 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 21:24:14 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v4z.mol lda_ll_fdekern.inp 
   PID         : 14301
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v4z.14301


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v4z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1244 bytes at 000000010407fa28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1244 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.v4z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v4z.14301 on main (master) node fpcmtaspg:
total 100840
-rw-r--r--  1 gomes  wheel   1174476 Mar 26 21:24 AOMOMAT
-rw-r--r--  1 gomes  wheel   3883648 Mar 26 21:24 AOPROPER
-rw-r--r--  1 gomes  wheel    111400 Mar 26 21:24 DFCMOS
-rw-r--r--  1 gomes  wheel    588006 Mar 26 21:24 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 21:24 DFCYCL
-rw-r--r--  1 gomes  wheel   3095162 Mar 26 21:24 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 21:24 DFDIIS
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 21:24 DFEVEC
-rw-r--r--  1 gomes  wheel   3095178 Mar 26 21:24 DFFCK1
-rw-r--r--  1 gomes  wheel   3095178 Mar 26 21:24 DFFCK2
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 21:24 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 21:24 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 21:24 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 21:24 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 21:24 FRZDN1
-rw-r--r--  1 gomes  wheel    111520 Mar 26 21:24 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       792 Mar 26 21:24 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 21:24 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2332800 Mar 26 21:33 PAMBOE
-rw-r--r--  1 gomes  wheel    111400 Mar 26 21:24 PAMFCK
-rw-r--r--  1 gomes  wheel       284 Mar 26 21:33 PAMRSP
-rw-r--r--  1 gomes  wheel    260340 Mar 26 21:33 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 21:24 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 21:24 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 49M	/tmp/gomes/gomes/ar.v4z.14301
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v4z.14301 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 21:33:10 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 21:33:10 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v4z.mol lda_ll_nofde.inp 
   PID         : 14433
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v4z.14433


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v4z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000129029c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129003a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000128a00278 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 640000 bytes at 0000000128901028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 1244 bytes at 0000000129001e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129000028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128a003a8 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000128a00338 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000128a002e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128a00438 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 2200 bytes at 0000000129002828 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129025e28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129024028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129027c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129028c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1244 bytes at 0000000129050028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012902dc28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012902cc28 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 000000012902bc28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012902ac28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000012902ec28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 0000000129040028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129033c28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.v4z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v4z.14433 on main (master) node fpcmtaspg:
total 99592
-rw-r--r--  1 gomes  wheel   1174476 Mar 26 21:33 AOMOMAT
-rw-r--r--  1 gomes  wheel   3495472 Mar 26 21:33 AOPROPER
-rw-r--r--  1 gomes  wheel    111400 Mar 26 21:33 DFCMOS
-rw-r--r--  1 gomes  wheel    588006 Mar 26 21:33 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 21:33 DFCYCL
-rw-r--r--  1 gomes  wheel   3095162 Mar 26 21:33 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 21:33 DFDIIS
-rw-r--r--  1 gomes  wheel    779744 Mar 26 21:33 DFEVEC
-rw-r--r--  1 gomes  wheel   3095178 Mar 26 21:33 DFFCK1
-rw-r--r--  1 gomes  wheel   3095178 Mar 26 21:33 DFFCK2
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 21:33 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 21:33 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 21:33 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 21:33 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 21:33 FRZDN1
-rw-r--r--  1 gomes  wheel    111520 Mar 26 21:33 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       792 Mar 26 21:33 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 21:33 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2304000 Mar 26 21:39 PAMBOE
-rw-r--r--  1 gomes  wheel    111400 Mar 26 21:33 PAMFCK
-rw-r--r--  1 gomes  wheel       156 Mar 26 21:39 PAMRSP
-rw-r--r--  1 gomes  wheel    260340 Mar 26 21:39 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 21:33 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 21:33 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 49M	/tmp/gomes/gomes/ar.v4z.14433
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v4z.14433 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 21:39:39 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 21:39:39 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v4z.mol lda_ll_nokern.inp 
   PID         : 14540
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v4z.14540


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v4z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 00000001289b3028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2488 bytes at 0000000104080028 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.v4z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v4z.14540 on main (master) node fpcmtaspg:
total 101008
-rw-r--r--  1 gomes  wheel   1174476 Mar 26 21:39 AOMOMAT
-rw-r--r--  1 gomes  wheel   3883648 Mar 26 21:39 AOPROPER
-rw-r--r--  1 gomes  wheel    111400 Mar 26 21:40 DFCMOS
-rw-r--r--  1 gomes  wheel    588006 Mar 26 21:40 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 21:40 DFCYCL
-rw-r--r--  1 gomes  wheel   3095162 Mar 26 21:40 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 21:40 DFDIIS
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 21:40 DFEVEC
-rw-r--r--  1 gomes  wheel   3095178 Mar 26 21:39 DFFCK1
-rw-r--r--  1 gomes  wheel   3095178 Mar 26 21:40 DFFCK2
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 21:40 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 21:39 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 21:39 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 21:39 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 21:39 FRZDN1
-rw-r--r--  1 gomes  wheel    111520 Mar 26 21:40 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       792 Mar 26 21:39 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 21:39 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2419200 Mar 26 21:48 PAMBOE
-rw-r--r--  1 gomes  wheel    111400 Mar 26 21:40 PAMFCK
-rw-r--r--  1 gomes  wheel       292 Mar 26 21:48 PAMRSP
-rw-r--r--  1 gomes  wheel    260340 Mar 26 21:48 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 21:39 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 21:39 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 49M	/tmp/gomes/gomes/ar.v4z.14540
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v4z.14540 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 21:48:10 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 21:48:10 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av5z.mol lda_ll_fdekern.inp 
   PID         : 14649
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av5z.14649


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av5z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2672 bytes at 0000000104080628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2672 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.av5z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av5z.14649 on main (master) node fpcmtaspg:
total 230712
-rw-r--r--  1 gomes  wheel   5600652 Mar 26 21:48 AOMOMAT
-rw-r--r--  1 gomes  wheel  17880120 Mar 26 21:48 AOPROPER
-rw-r--r--  1 gomes  wheel    549160 Mar 26 21:50 DFCMOS
-rw-r--r--  1 gomes  wheel   2801958 Mar 26 21:50 DFCOEF
-rw-r--r--  1 gomes  wheel      1380 Mar 26 21:50 DFCYCL
-rw-r--r--  1 gomes  wheel  14279258 Mar 26 21:50 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 21:50 DFDIIS
-rw-r--r--  1 gomes  wheel   4942368 Mar 26 21:50 DFEVEC
-rw-r--r--  1 gomes  wheel  14279274 Mar 26 21:48 DFFCK1
-rw-r--r--  1 gomes  wheel  14279274 Mar 26 21:50 DFFCK2
-rw-r--r--  1 gomes  wheel   4942368 Mar 26 21:50 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 21:48 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 21:48 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 21:48 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 21:48 FRZDN1
-rw-r--r--  1 gomes  wheel    549280 Mar 26 21:50 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel      1112 Mar 26 21:48 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 21:48 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   5059584 Mar 26 22:20 PAMBOE
-rw-r--r--  1 gomes  wheel    549160 Mar 26 21:50 PAMFCK
-rw-r--r--  1 gomes  wheel       276 Mar 26 22:21 PAMRSP
-rw-r--r--  1 gomes  wheel    633588 Mar 26 22:21 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 21:48 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 21:48 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory:113M	/tmp/gomes/gomes/ar.av5z.14649
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  780Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av5z.14649 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 22:21:18 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 22:21:18 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av5z.mol lda_ll_nofde.inp 
   PID         : 14772
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av5z.14772


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av5z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 000000012900e428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900b428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000128e00238 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 640000 bytes at 0000000128cd2028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129000028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128e00368 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000128e002f8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000128e002a8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128e00478 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 7600 bytes at 0000000129009628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129007828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900c428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900d428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2672 bytes at 0000000129034c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129012428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129011428 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000129010428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900f428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000129013e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 0000000129013428 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 0000000129025228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129018e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000129019e28 allocated at line 428 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_ar.av5z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av5z.14772 on main (master) node fpcmtaspg:
total 226416
-rw-r--r--  1 gomes  wheel   5600652 Mar 26 22:21 AOMOMAT
-rw-r--r--  1 gomes  wheel  16092504 Mar 26 22:21 AOPROPER
-rw-r--r--  1 gomes  wheel    549160 Mar 26 22:22 DFCMOS
-rw-r--r--  1 gomes  wheel   2801958 Mar 26 22:22 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 22:22 DFCYCL
-rw-r--r--  1 gomes  wheel  14279258 Mar 26 22:22 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 22:22 DFDIIS
-rw-r--r--  1 gomes  wheel   4393216 Mar 26 22:22 DFEVEC
-rw-r--r--  1 gomes  wheel  14279274 Mar 26 22:21 DFFCK1
-rw-r--r--  1 gomes  wheel  14279274 Mar 26 22:22 DFFCK2
-rw-r--r--  1 gomes  wheel   4942368 Mar 26 22:22 DFFOCK
-rw-r--r--  1 gomes  wheel       413 Mar 26 22:21 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 22:21 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 22:21 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 22:21 FRZDN1
-rw-r--r--  1 gomes  wheel    549280 Mar 26 22:22 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel      1112 Mar 26 22:21 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 22:21 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   5200128 Mar 26 23:22 PAMBOE
-rw-r--r--  1 gomes  wheel    549160 Mar 26 22:22 PAMFCK
-rw-r--r--  1 gomes  wheel       876 Mar 26 23:22 PAMRSP
-rw-r--r--  1 gomes  wheel    633588 Mar 26 23:22 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 22:21 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 22:21 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory:111M	/tmp/gomes/gomes/ar.av5z.14772
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  780Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av5z.14772 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
 *** WARNING-Maximum number of microiterations,  50, reached
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
 *** WARNING-Maximum number of microiterations,  50, reached
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
 *** WARNING-Maximum number of microiterations,  50, reached

   Exit date: Fri Mar 26 23:22:43 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 23:22:43 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av5z.mol lda_ll_nokern.inp 
   PID         : 14919
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.av5z.14919


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.av5z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128ab5028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2672 bytes at 0000000104080628 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_ar.av5z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av5z.14919 on main (master) node fpcmtaspg:
total 232360
-rw-r--r--  1 gomes  wheel   5600652 Mar 26 23:23 AOMOMAT
-rw-r--r--  1 gomes  wheel  17880120 Mar 26 23:23 AOPROPER
-rw-r--r--  1 gomes  wheel    549160 Mar 26 23:24 DFCMOS
-rw-r--r--  1 gomes  wheel   2801958 Mar 26 23:25 DFCOEF
-rw-r--r--  1 gomes  wheel      1380 Mar 26 23:25 DFCYCL
-rw-r--r--  1 gomes  wheel  14279258 Mar 26 23:24 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 23:24 DFDIIS
-rw-r--r--  1 gomes  wheel   4942368 Mar 26 23:24 DFEVEC
-rw-r--r--  1 gomes  wheel  14279274 Mar 26 23:23 DFFCK1
-rw-r--r--  1 gomes  wheel  14279274 Mar 26 23:24 DFFCK2
-rw-r--r--  1 gomes  wheel   4942368 Mar 26 23:25 DFFOCK
-rw-r--r--  1 gomes  wheel       677 Mar 26 23:22 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 23:22 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 23:22 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 23:22 FRZDN1
-rw-r--r--  1 gomes  wheel    549280 Mar 26 23:25 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel      1112 Mar 26 23:22 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 23:22 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   5902848 Mar 26 23:58 PAMBOE
-rw-r--r--  1 gomes  wheel    549160 Mar 26 23:25 PAMFCK
-rw-r--r--  1 gomes  wheel       292 Mar 26 23:58 PAMRSP
-rw-r--r--  1 gomes  wheel    633588 Mar 26 23:58 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 23:22 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 23:22 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory:113M	/tmp/gomes/gomes/ar.av5z.14919
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  780Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av5z.14919 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Fri Mar 26 23:58:59 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 23:58:59 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v5z.mol lda_ll_fdekern.inp 
   PID         : 15049
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v5z.15049


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v5z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_fdekern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f9a028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1976 bytes at 000000010407fe28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1976 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_ar.v5z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v5z.15049 on main (master) node fpcmtaspg:
total 156264
-rw-r--r--  1 gomes  wheel   3003660 Mar 26 23:59 AOMOMAT
-rw-r--r--  1 gomes  wheel   9784656 Mar 26 23:59 AOPROPER
-rw-r--r--  1 gomes  wheel    288808 Mar 27 00:00 DFCMOS
-rw-r--r--  1 gomes  wheel   1503030 Mar 27 00:00 DFCOEF
-rw-r--r--  1 gomes  wheel      1610 Mar 27 00:00 DFCYCL
-rw-r--r--  1 gomes  wheel   7809242 Mar 27 00:00 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 27 00:00 DFDIIS
-rw-r--r--  1 gomes  wheel   2599200 Mar 27 00:00 DFEVEC
-rw-r--r--  1 gomes  wheel   7809258 Mar 26 23:59 DFFCK1
-rw-r--r--  1 gomes  wheel   7809258 Mar 27 00:00 DFFCK2
-rw-r--r--  1 gomes  wheel   2599200 Mar 27 00:00 DFFOCK
-rw-r--r--  1 gomes  wheel       685 Mar 26 23:58 DIRAC.INP
-rw-r--r--  1 gomes  wheel   4807596 Mar 26 23:59 EMBGRD
-rw-r--r--  1 gomes  wheel   5701041 Mar 26 23:59 EMBPOT
-rw-r--r--  1 gomes  wheel   6849763 Mar 26 23:59 FRZDN1
-rw-r--r--  1 gomes  wheel    288928 Mar 27 00:00 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       986 Mar 26 23:59 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 23:58 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   4012416 Mar 27 00:16 PAMBOE
-rw-r--r--  1 gomes  wheel    288808 Mar 27 00:00 PAMFCK
-rw-r--r--  1 gomes  wheel       252 Mar 27 00:17 PAMRSP
-rw-r--r--  1 gomes  wheel    446964 Mar 27 00:17 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 23:59 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 23:59 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 76M	/tmp/gomes/gomes/ar.v5z.15049
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  151Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v5z.15049 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !
*** WARNING - Possible instability detected in reduced system !

   Exit date: Sat Mar 27 00:17:08 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Sat Mar 27 00:17:08 CET 2010
   Command     : ../../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v5z.mol lda_ll_nofde.inp 
   PID         : 15161
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot
   Work dir    : /tmp/gomes/gomes/ar.v5z.15161


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/ar.v5z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/lda_ll_nofde.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/qzpot/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

