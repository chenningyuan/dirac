   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:29:56 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_fdekern_unc.inp 
   PID         : 9163
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.9163


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 720 bytes at 0000000103c21ac8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21078 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cd8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 720 bytes at 0000000103c214c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211e8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c210e8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21158 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 720 bytes at 0000000103c217c8 allocated at line 538 of gp/memory_allocator.F90
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
Remaining memory: 1440 bytes at 000000010407ee28 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_unc_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.9163 on main (master) node fpcmtaspg:
total 70312
-rw-r--r--  1 gomes  wheel    576140 Mar 26 18:30 AOMOMAT
-rw-r--r--  1 gomes  wheel   1304736 Mar 26 18:30 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 18:30 DFCMOS
-rw-r--r--  1 gomes  wheel    288730 Mar 26 18:30 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 18:30 DFCYCL
-rw-r--r--  1 gomes  wheel   1036890 Mar 26 18:30 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:30 DFDIIS
-rw-r--r--  1 gomes  wheel    640000 Mar 26 18:30 DFEVEC
-rw-r--r--  1 gomes  wheel   1036906 Mar 26 18:30 DFFCK1
-rw-r--r--  1 gomes  wheel   1036906 Mar 26 18:30 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 18:30 DFFOCK
-rw-r--r--  1 gomes  wheel       707 Mar 26 18:29 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:29 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:29 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:29 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 18:30 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       593 Mar 26 18:30 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:29 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1771200 Mar 26 18:34 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 18:30 PAMFCK
-rw-r--r--  1 gomes  wheel       252 Mar 26 18:34 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 18:34 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:29 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:29 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 34M	/tmp/gomes/gomes/ar.av2z.9163
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.9163 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 18:34:59 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:35:00 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_nofde_unc.inp 
   PID         : 9271
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.9271


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nofde_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000104064628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104061628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21078 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cd8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 720 bytes at 0000000103c21278 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21158 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c210e8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211e8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 720 bytes at 0000000103c21878 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 720 bytes at 0000000103c21578 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405bc28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010400f628 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405da28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104062628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104063628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128b7b028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067628 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000104069628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 000000010407aa28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406e628 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_unc_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.9271 on main (master) node fpcmtaspg:
total 69872
-rw-r--r--  1 gomes  wheel    576140 Mar 26 18:35 AOMOMAT
-rw-r--r--  1 gomes  wheel   1174368 Mar 26 18:35 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 18:35 DFCMOS
-rw-r--r--  1 gomes  wheel    288730 Mar 26 18:35 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 18:35 DFCYCL
-rw-r--r--  1 gomes  wheel   1036890 Mar 26 18:35 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:35 DFDIIS
-rw-r--r--  1 gomes  wheel    640000 Mar 26 18:35 DFEVEC
-rw-r--r--  1 gomes  wheel   1036906 Mar 26 18:35 DFFCK1
-rw-r--r--  1 gomes  wheel   1036906 Mar 26 18:35 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 18:35 DFFOCK
-rw-r--r--  1 gomes  wheel       435 Mar 26 18:35 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:35 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:35 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:35 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 18:35 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       593 Mar 26 18:35 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:35 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1676736 Mar 26 18:38 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 18:35 PAMFCK
-rw-r--r--  1 gomes  wheel       156 Mar 26 18:38 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 18:38 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:35 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:35 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 34M	/tmp/gomes/gomes/ar.av2z.9271
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.9271 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 18:38:31 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:38:31 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_nokern_unc.inp 
   PID         : 9376
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.9376


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nokern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21158 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c210e8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21078 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cd8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 720 bytes at 0000000103c21ac8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211e8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 720 bytes at 0000000103c214c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 720 bytes at 0000000103c217c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 000000012886f028 allocated at line 374 of gp/memory_allocator.F90
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
   Creating archive file lda_ll_nokern_unc_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.9376 on main (master) node fpcmtaspg:
total 70496
-rw-r--r--  1 gomes  wheel    576140 Mar 26 18:38 AOMOMAT
-rw-r--r--  1 gomes  wheel   1304736 Mar 26 18:38 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 18:38 DFCMOS
-rw-r--r--  1 gomes  wheel    288730 Mar 26 18:38 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 18:38 DFCYCL
-rw-r--r--  1 gomes  wheel   1036890 Mar 26 18:38 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:38 DFDIIS
-rw-r--r--  1 gomes  wheel    640000 Mar 26 18:38 DFEVEC
-rw-r--r--  1 gomes  wheel   1036906 Mar 26 18:38 DFFCK1
-rw-r--r--  1 gomes  wheel   1036906 Mar 26 18:38 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 18:38 DFFOCK
-rw-r--r--  1 gomes  wheel       699 Mar 26 18:38 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:38 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:38 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:38 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 18:38 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       593 Mar 26 18:38 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:38 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1865664 Mar 26 18:42 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 18:38 PAMFCK
-rw-r--r--  1 gomes  wheel       260 Mar 26 18:42 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 18:42 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:38 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:38 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 34M	/tmp/gomes/gomes/ar.av2z.9376
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.9376 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 18:42:28 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:42:28 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_ll_fdekern_noauto.inp 
   PID         : 9476
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.9476


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_noauto.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 40 bytes at 0000000103c210e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21078 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 63536 bytes at 00000001040b8e28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c21768 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c214c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21138 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211a8 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 96 bytes at 0000000103c21238 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 628 bytes at 0000000103c21a08 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 1256 bytes at 000000010407ee28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_noauto_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.9476 on main (master) node fpcmtaspg:
total 62688
-rw-r--r--  1 gomes  wheel    271436 Mar 26 18:42 AOMOMAT
-rw-r--r--  1 gomes  wheel    993456 Mar 26 18:42 AOPROPER
-rw-r--r--  1 gomes  wheel     23336 Mar 26 18:42 DFCMOS
-rw-r--r--  1 gomes  wheel    136102 Mar 26 18:42 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 18:42 DFCYCL
-rw-r--r--  1 gomes  wheel    788858 Mar 26 18:42 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:42 DFDIIS
-rw-r--r--  1 gomes  wheel    163296 Mar 26 18:42 DFEVEC
-rw-r--r--  1 gomes  wheel    788874 Mar 26 18:42 DFFCK1
-rw-r--r--  1 gomes  wheel    788874 Mar 26 18:42 DFFCK2
-rw-r--r--  1 gomes  wheel    209952 Mar 26 18:42 DFFOCK
-rw-r--r--  1 gomes  wheel       681 Mar 26 18:42 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:42 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:42 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:42 FRZDN1
-rw-r--r--  1 gomes  wheel     23456 Mar 26 18:42 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 18:42 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:42 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    663552 Mar 26 18:43 PAMBOE
-rw-r--r--  1 gomes  wheel     23336 Mar 26 18:42 PAMFCK
-rw-r--r--  1 gomes  wheel       268 Mar 26 18:43 PAMRSP
-rw-r--r--  1 gomes  wheel     94452 Mar 26 18:43 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:42 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:42 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 31M	/tmp/gomes/gomes/ar.av2z.9476
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.9476 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 18:43:45 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:43:45 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_fdekern_unc.inp 
   PID         : 9579
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.9579


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21078 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cd8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211e8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c210e8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21158 allocated at line 706 of dft/dft_ao.F90
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
Remaining memory: 1132 bytes at 000000010407fa28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1132 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1132 bytes at 000000010407f428 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_unc_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.9579 on main (master) node fpcmtaspg:
total 91576
-rw-r--r--  1 gomes  wheel   1358540 Mar 26 18:43 AOMOMAT
-rw-r--r--  1 gomes  wheel   3217496 Mar 26 18:43 AOPROPER
-rw-r--r--  1 gomes  wheel    180008 Mar 26 18:44 DFCMOS
-rw-r--r--  1 gomes  wheel    680230 Mar 26 18:44 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 18:44 DFCYCL
-rw-r--r--  1 gomes  wheel   2562938 Mar 26 18:44 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:44 DFDIIS
-rw-r--r--  1 gomes  wheel   1620000 Mar 26 18:44 DFEVEC
-rw-r--r--  1 gomes  wheel   2562954 Mar 26 18:43 DFFCK1
-rw-r--r--  1 gomes  wheel   2562954 Mar 26 18:44 DFFCK2
-rw-r--r--  1 gomes  wheel   1620000 Mar 26 18:44 DFFOCK
-rw-r--r--  1 gomes  wheel       707 Mar 26 18:43 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:43 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:43 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:43 FRZDN1
-rw-r--r--  1 gomes  wheel    180128 Mar 26 18:44 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       747 Mar 26 18:43 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:43 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2699136 Mar 26 18:50 PAMBOE
-rw-r--r--  1 gomes  wheel    180008 Mar 26 18:44 PAMFCK
-rw-r--r--  1 gomes  wheel       228 Mar 26 18:50 PAMRSP
-rw-r--r--  1 gomes  wheel    343284 Mar 26 18:50 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:43 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:43 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 45M	/tmp/gomes/gomes/ar.av3z.9579
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.9579 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 18:50:59 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:50:59 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_nofde_unc.inp 
   PID         : 9694
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.9694


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nofde_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000104084e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104081e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c214d8 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 7600 bytes at 000000010407c428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21608 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21598 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000103c21548 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21698 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 2200 bytes at 000000010400f628 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104080028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010407e228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104082e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104083e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1132 bytes at 00000001040aac28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104088e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104087e28 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104086e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104085e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000104089e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 000000010409b228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010408ee28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 000000010408fe28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 1132 bytes at 00000001040ab828 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1132 bytes at 00000001040ab228 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128b00028 allocated at line 563 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_unc_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.9694 on main (master) node fpcmtaspg:
total 91624
-rw-r--r--  1 gomes  wheel   1358540 Mar 26 18:51 AOMOMAT
-rw-r--r--  1 gomes  wheel   2895960 Mar 26 18:51 AOPROPER
-rw-r--r--  1 gomes  wheel    180008 Mar 26 18:51 DFCMOS
-rw-r--r--  1 gomes  wheel    680230 Mar 26 18:51 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 18:51 DFCYCL
-rw-r--r--  1 gomes  wheel   2562938 Mar 26 18:51 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:51 DFDIIS
-rw-r--r--  1 gomes  wheel   1620000 Mar 26 18:51 DFEVEC
-rw-r--r--  1 gomes  wheel   2562954 Mar 26 18:51 DFFCK1
-rw-r--r--  1 gomes  wheel   2562954 Mar 26 18:51 DFFCK2
-rw-r--r--  1 gomes  wheel   1620000 Mar 26 18:51 DFFOCK
-rw-r--r--  1 gomes  wheel       435 Mar 26 18:50 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:50 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:50 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:50 FRZDN1
-rw-r--r--  1 gomes  wheel    180128 Mar 26 18:51 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       747 Mar 26 18:51 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:50 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   3041280 Mar 26 18:57 PAMBOE
-rw-r--r--  1 gomes  wheel    180008 Mar 26 18:51 PAMFCK
-rw-r--r--  1 gomes  wheel       156 Mar 26 18:57 PAMRSP
-rw-r--r--  1 gomes  wheel    343284 Mar 26 18:57 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:50 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:51 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 45M	/tmp/gomes/gomes/ar.av3z.9694
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.9694 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 18:57:15 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 18:57:15 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_nokern_unc.inp 
   PID         : 10101
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.10101


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nokern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21158 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c210e8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21078 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cd8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211e8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 00000001288ce028 allocated at line 374 of gp/memory_allocator.F90
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
Remaining memory: 2264 bytes at 0000000104080028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 1132 bytes at 000000010407fa28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_unc_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.10101 on main (master) node fpcmtaspg:
total 92696
-rw-r--r--  1 gomes  wheel   1358540 Mar 26 18:57 AOMOMAT
-rw-r--r--  1 gomes  wheel   3217496 Mar 26 18:57 AOPROPER
-rw-r--r--  1 gomes  wheel    180008 Mar 26 18:57 DFCMOS
-rw-r--r--  1 gomes  wheel    680230 Mar 26 18:57 DFCOEF
-rw-r--r--  1 gomes  wheel      1265 Mar 26 18:57 DFCYCL
-rw-r--r--  1 gomes  wheel   2562938 Mar 26 18:57 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 18:57 DFDIIS
-rw-r--r--  1 gomes  wheel   1620000 Mar 26 18:57 DFEVEC
-rw-r--r--  1 gomes  wheel   2562954 Mar 26 18:57 DFFCK1
-rw-r--r--  1 gomes  wheel   2562954 Mar 26 18:57 DFFCK2
-rw-r--r--  1 gomes  wheel   1620000 Mar 26 18:57 DFFOCK
-rw-r--r--  1 gomes  wheel       699 Mar 26 18:57 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 18:57 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 18:57 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 18:57 FRZDN1
-rw-r--r--  1 gomes  wheel    180128 Mar 26 18:57 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       747 Mar 26 18:57 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 18:57 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   3269376 Mar 26 19:04 PAMBOE
-rw-r--r--  1 gomes  wheel    180008 Mar 26 18:57 PAMFCK
-rw-r--r--  1 gomes  wheel       204 Mar 26 19:04 PAMRSP
-rw-r--r--  1 gomes  wheel    343284 Mar 26 19:04 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 18:57 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 18:57 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 45M	/tmp/gomes/gomes/ar.av3z.10101
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.10101 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:04:26 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:04:26 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_ll_fdekern_noauto.inp 
   PID         : 10508
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.10508


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_noauto.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 40 bytes at 0000000103c21108 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21158 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211c8 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 96 bytes at 0000000103c21258 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2064 bytes at 0000000104080028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 1032 bytes at 000000010407f428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 1032 bytes at 000000010407ee28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_noauto_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.10508 on main (master) node fpcmtaspg:
total 80024
-rw-r--r--  1 gomes  wheel    825740 Mar 26 19:04 AOMOMAT
-rw-r--r--  1 gomes  wheel   2674600 Mar 26 19:04 AOPROPER
-rw-r--r--  1 gomes  wheel     80008 Mar 26 19:04 DFCMOS
-rw-r--r--  1 gomes  wheel    413530 Mar 26 19:04 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 19:04 DFCYCL
-rw-r--r--  1 gomes  wheel   2130138 Mar 26 19:04 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:04 DFDIIS
-rw-r--r--  1 gomes  wheel    640000 Mar 26 19:04 DFEVEC
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 19:04 DFFCK1
-rw-r--r--  1 gomes  wheel   2130154 Mar 26 19:04 DFFCK2
-rw-r--r--  1 gomes  wheel    720000 Mar 26 19:04 DFFOCK
-rw-r--r--  1 gomes  wheel       681 Mar 26 19:04 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:04 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:04 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:04 FRZDN1
-rw-r--r--  1 gomes  wheel     80128 Mar 26 19:04 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       745 Mar 26 19:04 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 19:04 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1723968 Mar 26 19:08 PAMBOE
-rw-r--r--  1 gomes  wheel     80008 Mar 26 19:04 PAMFCK
-rw-r--r--  1 gomes  wheel       228 Mar 26 19:08 PAMRSP
-rw-r--r--  1 gomes  wheel    213684 Mar 26 19:08 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:04 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:04 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 39M	/tmp/gomes/gomes/ar.av3z.10508
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.10508 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:08:34 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:08:34 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_fdekern_unc.inp 
   PID         : 10713
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v2z.10713


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040c8e28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 00000001040bda28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 588 bytes at 0000000103c219e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 588 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 588 bytes at 0000000103c21768 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 00000001040b6028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407f428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104076828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104074a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407d428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407b428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104078628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407a428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407c428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407e428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104080428 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 20000 bytes at 00000001040b7a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 00000001040b7028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 00000001040bca28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_unc_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.10713 on main (master) node fpcmtaspg:
total 64824
-rw-r--r--  1 gomes  wheel    385868 Mar 26 19:08 AOMOMAT
-rw-r--r--  1 gomes  wheel    871776 Mar 26 19:08 AOPROPER
-rw-r--r--  1 gomes  wheel     53800 Mar 26 19:08 DFCMOS
-rw-r--r--  1 gomes  wheel    193486 Mar 26 19:08 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 19:08 DFCYCL
-rw-r--r--  1 gomes  wheel    691578 Mar 26 19:08 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:08 DFDIIS
-rw-r--r--  1 gomes  wheel    430336 Mar 26 19:08 DFEVEC
-rw-r--r--  1 gomes  wheel    691594 Mar 26 19:08 DFFCK1
-rw-r--r--  1 gomes  wheel    691594 Mar 26 19:08 DFFCK2
-rw-r--r--  1 gomes  wheel    484128 Mar 26 19:08 DFFOCK
-rw-r--r--  1 gomes  wheel       707 Mar 26 19:08 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:08 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:08 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:08 FRZDN1
-rw-r--r--  1 gomes  wheel     53920 Mar 26 19:08 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       530 Mar 26 19:08 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:08 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1290240 Mar 26 19:11 PAMBOE
-rw-r--r--  1 gomes  wheel     53800 Mar 26 19:08 PAMFCK
-rw-r--r--  1 gomes  wheel       236 Mar 26 19:11 PAMRSP
-rw-r--r--  1 gomes  wheel    167028 Mar 26 19:11 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:08 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:08 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 32M	/tmp/gomes/gomes/ar.v2z.10713
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.10713 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:11:14 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:11:14 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_nofde_unc.inp 
   PID         : 10819
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v2z.10819


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nofde_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 0000000104094628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104091628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c214f8 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 588 bytes at 0000000103c21748 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21628 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c215b8 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000103c21568 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c216b8 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 588 bytes at 0000000103c21c48 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 588 bytes at 0000000103c219c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010407a028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 0000000104079628 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010408f828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 1176 bytes at 000000010407be28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010408da28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104092628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104093628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 000000012abd0028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 00000001040b5428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 00000001040b4428 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 00000001040b3428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104095628 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 00000001040b6428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 00000001040c7828 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_unc_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.10819 on main (master) node fpcmtaspg:
total 64904
-rw-r--r--  1 gomes  wheel    385868 Mar 26 19:11 AOMOMAT
-rw-r--r--  1 gomes  wheel    784704 Mar 26 19:11 AOPROPER
-rw-r--r--  1 gomes  wheel     53800 Mar 26 19:11 DFCMOS
-rw-r--r--  1 gomes  wheel    193486 Mar 26 19:11 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 19:11 DFCYCL
-rw-r--r--  1 gomes  wheel    691578 Mar 26 19:11 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:11 DFDIIS
-rw-r--r--  1 gomes  wheel    376544 Mar 26 19:11 DFEVEC
-rw-r--r--  1 gomes  wheel    691594 Mar 26 19:11 DFFCK1
-rw-r--r--  1 gomes  wheel    691594 Mar 26 19:11 DFFCK2
-rw-r--r--  1 gomes  wheel    484128 Mar 26 19:11 DFFOCK
-rw-r--r--  1 gomes  wheel       435 Mar 26 19:11 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:11 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:11 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:11 FRZDN1
-rw-r--r--  1 gomes  wheel     53920 Mar 26 19:11 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       530 Mar 26 19:11 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:11 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1474560 Mar 26 19:12 PAMBOE
-rw-r--r--  1 gomes  wheel     53800 Mar 26 19:11 PAMFCK
-rw-r--r--  1 gomes  wheel       188 Mar 26 19:12 PAMRSP
-rw-r--r--  1 gomes  wheel    167028 Mar 26 19:12 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:11 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:11 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 32M	/tmp/gomes/gomes/ar.v2z.10819
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.10819 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.

   Exit date: Fri Mar 26 19:12:55 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:12:55 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_nokern_unc.inp 
   PID         : 10928
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v2z.10928


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nokern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 0000000104074a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 588 bytes at 0000000103c219e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 588 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 588 bytes at 0000000103c21768 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407b428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104076828 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407a428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104078628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128802028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 00000001040bda28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 00000001040bca28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 00000001040b7028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104080428 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 000000010407e428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407c428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407d428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010407f428 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 00000001040b6028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 00000001040b7a28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nokern_unc_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.10928 on main (master) node fpcmtaspg:
total 65040
-rw-r--r--  1 gomes  wheel    385868 Mar 26 19:12 AOMOMAT
-rw-r--r--  1 gomes  wheel    871776 Mar 26 19:12 AOPROPER
-rw-r--r--  1 gomes  wheel     53800 Mar 26 19:13 DFCMOS
-rw-r--r--  1 gomes  wheel    193486 Mar 26 19:13 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 19:13 DFCYCL
-rw-r--r--  1 gomes  wheel    691578 Mar 26 19:13 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:13 DFDIIS
-rw-r--r--  1 gomes  wheel    430336 Mar 26 19:13 DFEVEC
-rw-r--r--  1 gomes  wheel    691594 Mar 26 19:12 DFFCK1
-rw-r--r--  1 gomes  wheel    691594 Mar 26 19:13 DFFCK2
-rw-r--r--  1 gomes  wheel    484128 Mar 26 19:13 DFFOCK
-rw-r--r--  1 gomes  wheel       699 Mar 26 19:12 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:12 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:12 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:12 FRZDN1
-rw-r--r--  1 gomes  wheel     53920 Mar 26 19:13 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       530 Mar 26 19:12 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:12 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1400832 Mar 26 19:15 PAMBOE
-rw-r--r--  1 gomes  wheel     53800 Mar 26 19:13 PAMFCK
-rw-r--r--  1 gomes  wheel       236 Mar 26 19:15 PAMRSP
-rw-r--r--  1 gomes  wheel    167028 Mar 26 19:15 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:12 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:12 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 32M	/tmp/gomes/gomes/ar.v2z.10928
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.10928 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:15:19 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:15:19 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v2z.mol lda_ll_fdekern_noauto.inp 
   PID         : 11136
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v2z.11136


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_noauto.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 40 bytes at 0000000103c21108 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 63536 bytes at 00000001040a0e28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c21708 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21158 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211c8 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 96 bytes at 0000000103c21258 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 496 bytes at 0000000103c21928 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 992 bytes at 000000010406ea28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000104090a28 allocated at line 538 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_noauto_ar.v2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v2z.11136 on main (master) node fpcmtaspg:
total 58608
-rw-r--r--  1 gomes  wheel    142988 Mar 26 19:15 AOMOMAT
-rw-r--r--  1 gomes  wheel    621080 Mar 26 19:15 AOPROPER
-rw-r--r--  1 gomes  wheel     10376 Mar 26 19:15 DFCMOS
-rw-r--r--  1 gomes  wheel     71770 Mar 26 19:15 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 19:15 DFCYCL
-rw-r--r--  1 gomes  wheel    492122 Mar 26 19:15 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:15 DFDIIS
-rw-r--r--  1 gomes  wheel     82944 Mar 26 19:15 DFEVEC
-rw-r--r--  1 gomes  wheel    492138 Mar 26 19:15 DFFCK1
-rw-r--r--  1 gomes  wheel    492138 Mar 26 19:15 DFFCK2
-rw-r--r--  1 gomes  wheel     93312 Mar 26 19:15 DFFOCK
-rw-r--r--  1 gomes  wheel       681 Mar 26 19:15 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:15 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:15 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:15 FRZDN1
-rw-r--r--  1 gomes  wheel     10496 Mar 26 19:15 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       529 Mar 26 19:15 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:15 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    311040 Mar 26 19:16 PAMBOE
-rw-r--r--  1 gomes  wheel     10376 Mar 26 19:15 PAMFCK
-rw-r--r--  1 gomes  wheel       244 Mar 26 19:16 PAMRSP
-rw-r--r--  1 gomes  wheel     47796 Mar 26 19:16 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:15 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:15 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 29M	/tmp/gomes/gomes/ar.v2z.11136
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v2z.11136 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:16:10 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:16:10 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_fdekern_unc.inp 
   PID         : 11236
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v3z.11236


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
Remaining memory: 876 bytes at 0000000103c21c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 876 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 876 bytes at 0000000103c21888 allocated at line 538 of gp/memory_allocator.F90
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
   Creating archive file lda_ll_fdekern_unc_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.11236 on main (master) node fpcmtaspg:
total 77640
-rw-r--r--  1 gomes  wheel    827084 Mar 26 19:16 AOMOMAT
-rw-r--r--  1 gomes  wheel   1929336 Mar 26 19:16 AOPROPER
-rw-r--r--  1 gomes  wheel    111400 Mar 26 19:16 DFCMOS
-rw-r--r--  1 gomes  wheel    414310 Mar 26 19:16 DFCOEF
-rw-r--r--  1 gomes  wheel      1495 Mar 26 19:16 DFCYCL
-rw-r--r--  1 gomes  wheel   1534842 Mar 26 19:16 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:16 DFDIIS
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 19:16 DFEVEC
-rw-r--r--  1 gomes  wheel   1534858 Mar 26 19:16 DFFCK1
-rw-r--r--  1 gomes  wheel   1534858 Mar 26 19:16 DFFCK2
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 19:16 DFFOCK
-rw-r--r--  1 gomes  wheel       707 Mar 26 19:16 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:16 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:16 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:16 FRZDN1
-rw-r--r--  1 gomes  wheel    111520 Mar 26 19:16 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       661 Mar 26 19:16 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:16 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2246400 Mar 26 19:21 PAMBOE
-rw-r--r--  1 gomes  wheel    111400 Mar 26 19:16 PAMFCK
-rw-r--r--  1 gomes  wheel       284 Mar 26 19:21 PAMRSP
-rw-r--r--  1 gomes  wheel    260340 Mar 26 19:21 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:16 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:16 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 38M	/tmp/gomes/gomes/ar.v3z.11236
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.11236 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:21:42 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:21:42 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_nofde_unc.inp 
   PID         : 11443
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v3z.11443


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nofde_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 3800 bytes at 000000012900de28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900ae28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000128a00278 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 640000 bytes at 000000012384c028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000128736028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 876 bytes at 0000000128a004c8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128a003a8 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000128a00338 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 40 bytes at 0000000128a002e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000128a00438 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 876 bytes at 0000000128a00c08 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 876 bytes at 0000000128a00868 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129005428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 0000000129000028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129009028 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000129007228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900be28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900ce28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129011e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000129010e28 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 000000012900fe28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000012900ee28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 0000000129012e28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 63536 bytes at 0000000129024228 allocated at line 374 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_nofde_unc_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.11443 on main (master) node fpcmtaspg:
total 77536
-rw-r--r--  1 gomes  wheel    827084 Mar 26 19:21 AOMOMAT
-rw-r--r--  1 gomes  wheel   1736568 Mar 26 19:21 AOPROPER
-rw-r--r--  1 gomes  wheel    111400 Mar 26 19:21 DFCMOS
-rw-r--r--  1 gomes  wheel    414310 Mar 26 19:21 DFCOEF
-rw-r--r--  1 gomes  wheel      1495 Mar 26 19:21 DFCYCL
-rw-r--r--  1 gomes  wheel   1534842 Mar 26 19:21 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:21 DFDIIS
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 19:21 DFEVEC
-rw-r--r--  1 gomes  wheel   1534858 Mar 26 19:21 DFFCK1
-rw-r--r--  1 gomes  wheel   1534858 Mar 26 19:21 DFFCK2
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 19:21 DFFOCK
-rw-r--r--  1 gomes  wheel       435 Mar 26 19:21 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:21 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:21 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:21 FRZDN1
-rw-r--r--  1 gomes  wheel    111520 Mar 26 19:22 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       661 Mar 26 19:21 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:21 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2390400 Mar 26 19:26 PAMBOE
-rw-r--r--  1 gomes  wheel    111400 Mar 26 19:22 PAMFCK
-rw-r--r--  1 gomes  wheel       188 Mar 26 19:26 PAMRSP
-rw-r--r--  1 gomes  wheel    260340 Mar 26 19:26 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:21 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:21 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 38M	/tmp/gomes/gomes/ar.v3z.11443
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.11443 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:26:32 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:26:33 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_nokern_unc.inp 
   PID         : 11551
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v3z.11551


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_nokern_unc.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21178 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 64 bytes at 0000000103c21108 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 40 bytes at 0000000103c20cf8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 876 bytes at 0000000103c21c28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c21208 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 876 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 876 bytes at 0000000103c21888 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 190000 bytes at 0000000103f6b028 allocated at line 400 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 640000 bytes at 0000000128870028 allocated at line 374 of gp/memory_allocator.F90
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
   Creating archive file lda_ll_nokern_unc_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.11551 on main (master) node fpcmtaspg:
total 77640
-rw-r--r--  1 gomes  wheel    827084 Mar 26 19:26 AOMOMAT
-rw-r--r--  1 gomes  wheel   1929336 Mar 26 19:26 AOPROPER
-rw-r--r--  1 gomes  wheel    111400 Mar 26 19:26 DFCMOS
-rw-r--r--  1 gomes  wheel    414310 Mar 26 19:26 DFCOEF
-rw-r--r--  1 gomes  wheel      1495 Mar 26 19:26 DFCYCL
-rw-r--r--  1 gomes  wheel   1534842 Mar 26 19:26 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:26 DFDIIS
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 19:26 DFEVEC
-rw-r--r--  1 gomes  wheel   1534858 Mar 26 19:26 DFFCK1
-rw-r--r--  1 gomes  wheel   1534858 Mar 26 19:26 DFFCK2
-rw-r--r--  1 gomes  wheel   1002528 Mar 26 19:26 DFFOCK
-rw-r--r--  1 gomes  wheel       699 Mar 26 19:26 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:26 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:26 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:26 FRZDN1
-rw-r--r--  1 gomes  wheel    111520 Mar 26 19:26 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       661 Mar 26 19:26 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:26 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   2246400 Mar 26 19:31 PAMBOE
-rw-r--r--  1 gomes  wheel    111400 Mar 26 19:26 PAMFCK
-rw-r--r--  1 gomes  wheel       284 Mar 26 19:31 PAMRSP
-rw-r--r--  1 gomes  wheel    260340 Mar 26 19:31 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:26 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:26 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 38M	/tmp/gomes/gomes/ar.v3z.11551
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.11551 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:31:36 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:31:36 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.v3z.mol lda_ll_fdekern_noauto.inp 
   PID         : 11755
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.v3z.11755


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.v3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_ll_fdekern_noauto.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : ok
   Allocated nodes for parallel 1 cpu run (content of dirac_machines file):
fpcmtaspg.univ-lille1.fr
   All cpus are on one machine - fpcmtaspg - therefore no distributing files among other nodes.
   Executing ILP32, parallel, 1 cpu run, only on fpcmtaspg in working directory, command:
         " /opt/local/lib/openmpi/bin/mpiexec  -np 1 ./dirac.x" with BASDIR(S) :
      /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/
      /Users/gomes/devel/svn/Dirac_embrsp/basis/
      /Users/gomes/devel/svn/Dirac_embrsp/basis_dalton/ 

Remaining memory: 40 bytes at 0000000103c21108 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21098 allocated at line 682 of dft/dft_ao.F90
Remaining memory: 63536 bytes at 00000001040a0a28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104063228 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 000000010405f628 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c21818 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 776 bytes at 0000000103c214e8 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 64 bytes at 0000000103c21158 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 96 bytes at 0000000103c211c8 allocated at line 706 of dft/dft_ao.F90
Remaining memory: 96 bytes at 0000000103c21258 allocated at line 707 of dft/dft_ao.F90
Remaining memory: 776 bytes at 0000000103c21b48 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 7600 bytes at 0000000104061428 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 1552 bytes at 000000010407ee28 allocated at line 374 of gp/memory_allocator.F90
Remaining memory: 2200 bytes at 000000010406d028 allocated at line 563 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406b028 allocated at line 678 of dft/dft_ao.F90
Remaining memory: 3800 bytes at 0000000104069028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104067028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104065028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104066028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104068028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406a028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 000000010406c028 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 20000 bytes at 000000010406da28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 3800 bytes at 0000000104072a28 allocated at line 538 of gp/memory_allocator.F90
Remaining memory: 45600 bytes at 0000000104073a28 allocated at line 428 of gp/memory_allocator.F90
 ... More segments remain
   Creating archive file lda_ll_fdekern_noauto_ar.v3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.v3z.11755 on main (master) node fpcmtaspg:
total 67936
-rw-r--r--  1 gomes  wheel    422284 Mar 26 19:31 AOMOMAT
-rw-r--r--  1 gomes  wheel   1514920 Mar 26 19:31 AOPROPER
-rw-r--r--  1 gomes  wheel     37000 Mar 26 19:31 DFCMOS
-rw-r--r--  1 gomes  wheel    211610 Mar 26 19:31 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 19:31 DFCYCL
-rw-r--r--  1 gomes  wheel   1204442 Mar 26 19:31 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:31 DFDIIS
-rw-r--r--  1 gomes  wheel    295936 Mar 26 19:31 DFEVEC
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 19:31 DFFCK1
-rw-r--r--  1 gomes  wheel   1204458 Mar 26 19:31 DFFCK2
-rw-r--r--  1 gomes  wheel    332928 Mar 26 19:31 DFFOCK
-rw-r--r--  1 gomes  wheel       681 Mar 26 19:31 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:31 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:31 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:31 FRZDN1
-rw-r--r--  1 gomes  wheel     37120 Mar 26 19:31 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       660 Mar 26 19:31 MNF.INP
-rw-r--r--  1 gomes  wheel       159 Mar 26 19:31 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1008000 Mar 26 19:33 PAMBOE
-rw-r--r--  1 gomes  wheel     37000 Mar 26 19:31 PAMFCK
-rw-r--r--  1 gomes  wheel       252 Mar 26 19:33 PAMRSP
-rw-r--r--  1 gomes  wheel    130740 Mar 26 19:33 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:31 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:31 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 33M	/tmp/gomes/gomes/ar.v3z.11755
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.v3z.11755 on main (master) node fpcmtaspg deleted.
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

   Exit date: Fri Mar 26 19:33:40 CET 2010
   *** NORMAL EXIT
