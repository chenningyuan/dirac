   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:44:56 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_nonrel_fdekern.inp 
   PID         : 11971
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.11971


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_nonrel_fdekern.inp to working directory as DIRAC.INP : ok
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

[fpcmtaspg:12055] *** Process received signal ***
[fpcmtaspg:12055] Signal: Segmentation fault (11)
[fpcmtaspg:12055] Signal code: Address not mapped (1)
[fpcmtaspg:12055] Failing at address: 0x0
[fpcmtaspg:12055] [ 0] 2   libSystem.B.dylib                   0x00007fff813e8eaa _sigtramp + 26
[fpcmtaspg:12055] [ 1] 3   ???                                 0x00000000ffffffff 0x0 + 4294967295
[fpcmtaspg:12055] [ 2] 4   dirac.x                             0x00000001005c591d quit_ + 749
[fpcmtaspg:12055] [ 3] 5   dirac.x                             0x00000001003ec5fa numls1_ + 858
[fpcmtaspg:12055] [ 4] 6   dirac.x                             0x00000001003ed85c numlst_ + 364
[fpcmtaspg:12055] [ 5] 7   dirac.x                             0x000000010095553f ostring_ + 3343
[fpcmtaspg:12055] [ 6] 8   dirac.x                             0x0000000100974948 setxlr_ + 88
[fpcmtaspg:12055] [ 7] 9   dirac.x                             0x000000010097a212 prpxl1_ + 226
[fpcmtaspg:12055] [ 8] 10  dirac.x                             0x000000010097b4de prpxlr_ + 382
[fpcmtaspg:12055] [ 9] 11  dirac.x                             0x000000010092d9a5 pampr1_ + 693
[fpcmtaspg:12055] [10] 12  dirac.x                             0x000000010092ded9 pamprp_ + 473
[fpcmtaspg:12055] [11] 13  dirac.x                             0x000000010082d6cb pamdrv_ + 395
[fpcmtaspg:12055] [12] 14  dirac.x                             0x000000010082d7bb dirctl_ + 187
[fpcmtaspg:12055] [13] 15  dirac.x                             0x000000010082d85d dirac_ + 157
[fpcmtaspg:12055] [14] 16  dirac.x                             0x00000001008307cc MAIN_ + 652
[fpcmtaspg:12055] [15] 17  dirac.x                             0x0000000100b27fbe main + 14
[fpcmtaspg:12055] [16] 18  dirac.x                             0x0000000100000f26 _start + 224
[fpcmtaspg:12055] [17] 19  dirac.x                             0x0000000100000e45 start + 33
[fpcmtaspg:12055] [18] 20  ???                                 0x0000000000000001 0x0 + 1
[fpcmtaspg:12055] *** End of error message ***
--------------------------------------------------------------------------
mpiexec noticed that process rank 0 with PID 12055 on node fpcmtaspg.univ-lille1.fr exited on signal 11 (Segmentation fault).
--------------------------------------------------------------------------

   dirac.x returned non-zero exit code
   Creating archive file lda_nonrel_fdekern_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.11971 on main (master) node fpcmtaspg:
total 54480
-rw-r--r--  1 gomes  wheel     12740 Mar 26 19:44 AOMOMAT
-rw-r--r--  1 gomes  wheel     24856 Mar 26 19:44 AOPROPER
-rw-r--r--  1 gomes  wheel      5840 Mar 26 19:45 DFCMOS
-rw-r--r--  1 gomes  wheel      6718 Mar 26 19:45 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 19:45 DFCYCL
-rw-r--r--  1 gomes  wheel     27002 Mar 26 19:45 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:45 DFDIIS
-rw-r--r--  1 gomes  wheel     40824 Mar 26 19:45 DFEVEC
-rw-r--r--  1 gomes  wheel     27018 Mar 26 19:44 DFFCK1
-rw-r--r--  1 gomes  wheel     27018 Mar 26 19:45 DFFCK2
-rw-r--r--  1 gomes  wheel     52488 Mar 26 19:45 DFFOCK
-rw-r--r--  1 gomes  wheel       671 Mar 26 19:44 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:44 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:44 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:44 FRZDN1
-rw-r--r--  1 gomes  wheel      5960 Mar 26 19:45 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 19:44 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 19:44 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    419904 Mar 26 19:46 PAMBOE
-rw-r--r--  1 gomes  wheel      5840 Mar 26 19:45 PAMFCK
-rw-r--r--  1 gomes  wheel      2644 Mar 26 19:46 PAMRSP
-rw-r--r--  1 gomes  wheel    108664 Mar 26 19:46 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:44 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:44 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 27M	/tmp/gomes/gomes/ar.av2z.11971
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.11971 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
*** WARNING - Microiterations stopped after 20 iterations

   Exit date: Fri Mar 26 19:46:58 CET 2010
   *** ABNORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:46:58 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_nonrel_nofde.inp 
   PID         : 12076
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.12076


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_nonrel_nofde.inp to working directory as DIRAC.INP : ok
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

[fpcmtaspg:12159] *** Process received signal ***
[fpcmtaspg:12159] Signal: Segmentation fault (11)
[fpcmtaspg:12159] Signal code: Address not mapped (1)
[fpcmtaspg:12159] Failing at address: 0x0
[fpcmtaspg:12159] [ 0] 2   libSystem.B.dylib                   0x00007fff813e8eaa _sigtramp + 26
[fpcmtaspg:12159] [ 1] 3   ???                                 0x00000000ffffffff 0x0 + 4294967295
[fpcmtaspg:12159] [ 2] 4   dirac.x                             0x00000001005c591d quit_ + 749
[fpcmtaspg:12159] [ 3] 5   dirac.x                             0x00000001003ec5fa numls1_ + 858
[fpcmtaspg:12159] [ 4] 6   dirac.x                             0x00000001003ed85c numlst_ + 364
[fpcmtaspg:12159] [ 5] 7   dirac.x                             0x000000010095553f ostring_ + 3343
[fpcmtaspg:12159] [ 6] 8   dirac.x                             0x0000000100974948 setxlr_ + 88
[fpcmtaspg:12159] [ 7] 9   dirac.x                             0x000000010097a212 prpxl1_ + 226
[fpcmtaspg:12159] [ 8] 10  dirac.x                             0x000000010097b4de prpxlr_ + 382
[fpcmtaspg:12159] [ 9] 11  dirac.x                             0x000000010092d9a5 pampr1_ + 693
[fpcmtaspg:12159] [10] 12  dirac.x                             0x000000010092ded9 pamprp_ + 473
[fpcmtaspg:12159] [11] 13  dirac.x                             0x000000010082d6cb pamdrv_ + 395
[fpcmtaspg:12159] [12] 14  dirac.x                             0x000000010082d7bb dirctl_ + 187
[fpcmtaspg:12159] [13] 15  dirac.x                             0x000000010082d85d dirac_ + 157
[fpcmtaspg:12159] [14] 16  dirac.x                             0x00000001008307cc MAIN_ + 652
[fpcmtaspg:12159] [15] 17  dirac.x                             0x0000000100b27fbe main + 14
[fpcmtaspg:12159] [16] 18  dirac.x                             0x0000000100000f26 _start + 224
[fpcmtaspg:12159] [17] 19  dirac.x                             0x0000000100000e45 start + 33
[fpcmtaspg:12159] [18] 20  ???                                 0x0000000000000001 0x0 + 1
[fpcmtaspg:12159] *** End of error message ***
--------------------------------------------------------------------------
mpiexec noticed that process rank 0 with PID 12159 on node fpcmtaspg.univ-lille1.fr exited on signal 11 (Segmentation fault).
--------------------------------------------------------------------------

   dirac.x returned non-zero exit code
   Creating archive file lda_nonrel_nofde_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.12076 on main (master) node fpcmtaspg:
total 54472
-rw-r--r--  1 gomes  wheel     12740 Mar 26 19:47 AOMOMAT
-rw-r--r--  1 gomes  wheel     21328 Mar 26 19:47 AOPROPER
-rw-r--r--  1 gomes  wheel      5840 Mar 26 19:47 DFCMOS
-rw-r--r--  1 gomes  wheel      6718 Mar 26 19:47 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 19:47 DFCYCL
-rw-r--r--  1 gomes  wheel     27002 Mar 26 19:47 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:47 DFDIIS
-rw-r--r--  1 gomes  wheel     40824 Mar 26 19:47 DFEVEC
-rw-r--r--  1 gomes  wheel     27018 Mar 26 19:47 DFFCK1
-rw-r--r--  1 gomes  wheel     27018 Mar 26 19:47 DFFCK2
-rw-r--r--  1 gomes  wheel     52488 Mar 26 19:47 DFFOCK
-rw-r--r--  1 gomes  wheel       399 Mar 26 19:46 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:46 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:46 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:46 FRZDN1
-rw-r--r--  1 gomes  wheel      5960 Mar 26 19:47 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 19:47 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 19:46 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    419904 Mar 26 19:48 PAMBOE
-rw-r--r--  1 gomes  wheel      5840 Mar 26 19:47 PAMFCK
-rw-r--r--  1 gomes  wheel      2644 Mar 26 19:48 PAMRSP
-rw-r--r--  1 gomes  wheel    108664 Mar 26 19:48 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:46 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:46 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 27M	/tmp/gomes/gomes/ar.av2z.12076
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.12076 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
*** WARNING - Microiterations stopped after 20 iterations

   Exit date: Fri Mar 26 19:48:32 CET 2010
   *** ABNORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:48:32 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av2z.mol lda_nonrel_nokern.inp 
   PID         : 12181
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av2z.12181


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av2z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_nonrel_nokern.inp to working directory as DIRAC.INP : ok
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

[fpcmtaspg:12264] *** Process received signal ***
[fpcmtaspg:12264] Signal: Segmentation fault (11)
[fpcmtaspg:12264] Signal code: Address not mapped (1)
[fpcmtaspg:12264] Failing at address: 0x0
[fpcmtaspg:12264] [ 0] 2   libSystem.B.dylib                   0x00007fff813e8eaa _sigtramp + 26
[fpcmtaspg:12264] [ 1] 3   ???                                 0x00000000ffffffff 0x0 + 4294967295
[fpcmtaspg:12264] [ 2] 4   dirac.x                             0x00000001005c591d quit_ + 749
[fpcmtaspg:12264] [ 3] 5   dirac.x                             0x00000001003ec5fa numls1_ + 858
[fpcmtaspg:12264] [ 4] 6   dirac.x                             0x00000001003ed85c numlst_ + 364
[fpcmtaspg:12264] [ 5] 7   dirac.x                             0x000000010095553f ostring_ + 3343
[fpcmtaspg:12264] [ 6] 8   dirac.x                             0x0000000100974948 setxlr_ + 88
[fpcmtaspg:12264] [ 7] 9   dirac.x                             0x000000010097a212 prpxl1_ + 226
[fpcmtaspg:12264] [ 8] 10  dirac.x                             0x000000010097b4de prpxlr_ + 382
[fpcmtaspg:12264] [ 9] 11  dirac.x                             0x000000010092d9a5 pampr1_ + 693
[fpcmtaspg:12264] [10] 12  dirac.x                             0x000000010092ded9 pamprp_ + 473
[fpcmtaspg:12264] [11] 13  dirac.x                             0x000000010082d6cb pamdrv_ + 395
[fpcmtaspg:12264] [12] 14  dirac.x                             0x000000010082d7bb dirctl_ + 187
[fpcmtaspg:12264] [13] 15  dirac.x                             0x000000010082d85d dirac_ + 157
[fpcmtaspg:12264] [14] 16  dirac.x                             0x00000001008307cc MAIN_ + 652
[fpcmtaspg:12264] [15] 17  dirac.x                             0x0000000100b27fbe main + 14
[fpcmtaspg:12264] [16] 18  dirac.x                             0x0000000100000f26 _start + 224
[fpcmtaspg:12264] [17] 19  dirac.x                             0x0000000100000e45 start + 33
[fpcmtaspg:12264] [18] 20  ???                                 0x0000000000000001 0x0 + 1
[fpcmtaspg:12264] *** End of error message ***
--------------------------------------------------------------------------
mpiexec noticed that process rank 0 with PID 12264 on node fpcmtaspg.univ-lille1.fr exited on signal 11 (Segmentation fault).
--------------------------------------------------------------------------

   dirac.x returned non-zero exit code
   Creating archive file lda_nonrel_nokern_ar.av2z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av2z.12181 on main (master) node fpcmtaspg:
total 54464
-rw-r--r--  1 gomes  wheel     12740 Mar 26 19:48 AOMOMAT
-rw-r--r--  1 gomes  wheel     24856 Mar 26 19:48 AOPROPER
-rw-r--r--  1 gomes  wheel      5840 Mar 26 19:48 DFCMOS
-rw-r--r--  1 gomes  wheel      6718 Mar 26 19:48 DFCOEF
-rw-r--r--  1 gomes  wheel      1035 Mar 26 19:48 DFCYCL
-rw-r--r--  1 gomes  wheel     27002 Mar 26 19:48 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:48 DFDIIS
-rw-r--r--  1 gomes  wheel     40824 Mar 26 19:48 DFEVEC
-rw-r--r--  1 gomes  wheel     27018 Mar 26 19:48 DFFCK1
-rw-r--r--  1 gomes  wheel     27018 Mar 26 19:48 DFFCK2
-rw-r--r--  1 gomes  wheel     52488 Mar 26 19:48 DFFOCK
-rw-r--r--  1 gomes  wheel       663 Mar 26 19:48 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:48 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:48 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:48 FRZDN1
-rw-r--r--  1 gomes  wheel      5960 Mar 26 19:48 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       592 Mar 26 19:48 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 19:48 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel    412128 Mar 26 19:50 PAMBOE
-rw-r--r--  1 gomes  wheel      5840 Mar 26 19:48 PAMFCK
-rw-r--r--  1 gomes  wheel      2596 Mar 26 19:50 PAMRSP
-rw-r--r--  1 gomes  wheel    108664 Mar 26 19:50 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:48 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:48 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 27M	/tmp/gomes/gomes/ar.av2z.12181
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av2z.12181 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
*** WARNING - Microiterations stopped after 19 iterations

   Exit date: Fri Mar 26 19:50:06 CET 2010
   *** ABNORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:50:06 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_nonrel_fdekern.inp 
   PID         : 12286
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.12286


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_nonrel_fdekern.inp to working directory as DIRAC.INP : ok
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

mpiexec: killing job...

--------------------------------------------------------------------------
mpiexec noticed that process rank 0 with PID 12369 on node fpcmtaspg.univ-lille1.fr exited on signal 0 (Signal 0).
--------------------------------------------------------------------------
mpiexec: clean termination accomplished


   dirac.x returned non-zero exit code
   Creating archive file lda_nonrel_fdekern_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.12286 on main (master) node fpcmtaspg:
total 61776
-rw-r--r--  1 gomes  wheel     47412 Mar 26 19:50 AOMOMAT
-rw-r--r--  1 gomes  wheel     99672 Mar 26 19:50 AOPROPER
-rw-r--r--  1 gomes  wheel     20008 Mar 26 19:50 DFCMOS
-rw-r--r--  1 gomes  wheel     24330 Mar 26 19:50 DFCOEF
-rw-r--r--  1 gomes  wheel      1150 Mar 26 19:50 DFCYCL
-rw-r--r--  1 gomes  wheel    111482 Mar 26 19:50 DFDENS
-rw-r--r--  1 gomes  wheel       820 Mar 26 19:50 DFDIIS
-rw-r--r--  1 gomes  wheel    160000 Mar 26 19:50 DFEVEC
-rw-r--r--  1 gomes  wheel    111498 Mar 26 19:50 DFFCK1
-rw-r--r--  1 gomes  wheel    111498 Mar 26 19:50 DFFCK2
-rw-r--r--  1 gomes  wheel    180000 Mar 26 19:50 DFFOCK
-rw-r--r--  1 gomes  wheel       671 Mar 26 19:50 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:50 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:50 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:50 FRZDN1
-rw-r--r--  1 gomes  wheel     20128 Mar 26 19:50 KRMC_MATRICES
-rw-r--r--  1 gomes  wheel       745 Mar 26 19:50 MNF.INP
-rw-r--r--  1 gomes  wheel       163 Mar 26 19:50 MOLECULE.MOL
-rw-r--r--  1 gomes  wheel   1800720 Mar 26 19:56 PAMBOE
-rw-r--r--  1 gomes  wheel        31 Mar 26 19:50 PAMCYC
-rw-r--r--  1 gomes  wheel     20008 Mar 26 19:50 PAMFCK
-rw-r--r--  1 gomes  wheel      4796 Mar 26 19:56 PAMRSP
-rw-r--r--  1 gomes  wheel   1783008 Mar 26 19:56 PAMSOE
-rw-r--r--  1 gomes  wheel        24 Mar 26 19:50 PAMXVC
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:50 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:50 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 30M	/tmp/gomes/gomes/ar.av3z.12286
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.12286 on main (master) node fpcmtaspg deleted.
WARNING!!! Cartesian Gaussians must be specified
 **** WARNING **** Spinfree properties not tested.
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)
 **WARNING** XPPORD: ZERO METRIC EIGENVALUE(S)

   Exit date: Fri Mar 26 19:56:13 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:56:13 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_nonrel_nofde.inp 
   PID         : 12390
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.12390


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_nonrel_nofde.inp to working directory as DIRAC.INP : ok
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

mpiexec: killing job...

--------------------------------------------------------------------------
mpiexec noticed that process rank 0 with PID 12473 on node fpcmtaspg.univ-lille1.fr exited on signal 0 (Signal 0).
--------------------------------------------------------------------------
mpiexec: clean termination accomplished


   dirac.x returned non-zero exit code
   Creating archive file lda_nonrel_nofde_ar.av3z.tgz


  Contents of work directory /tmp/gomes/gomes/ar.av3z.12390 on main (master) node fpcmtaspg:
total 52904
-rw-r--r--  1 gomes  wheel       399 Mar 26 19:56 DIRAC.INP
-rw-r--r--  1 gomes  wheel   3558155 Mar 26 19:56 EMBGRD
-rw-r--r--  1 gomes  wheel   4219403 Mar 26 19:56 EMBPOT
-rw-r--r--  1 gomes  wheel   4943635 Mar 26 19:56 FRZDN1
-rw-r--r--  1 gomes  wheel       163 Mar 26 19:56 MOLECULE.MOL
-rwxr-xr-x  1 gomes  wheel  14346488 Mar 26 19:56 dirac.x
-rw-r--r--  1 gomes  wheel        25 Mar 26 19:56 dirac_machines
  Total size of the main (master) node fpcmtaspg work directory: 26M	/tmp/gomes/gomes/ar.av3z.12390
  Info on the disk space usage on the master node fpcmtaspg (df -h .):
Filesystem     Size   Used  Avail Capacity  Mounted on
/dev/disk0s2  931Gi  779Gi  152Gi    84%    /
Work directory /tmp/gomes/gomes/ar.av3z.12390 on main (master) node fpcmtaspg deleted.

   Exit date: Fri Mar 26 19:56:14 CET 2010
   *** NORMAL EXIT
   ****************************************
   **** OUTPUT FROM DIRAC SHELL SCRIPT ****
   ****************************************

   $Revision: 10284 $
   User@Host   : gomes@fpcmtaspg.univ-lille1.fr
   Date        : Fri Mar 26 19:56:14 CET 2010
   Command     : ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean ar.av3z.mol lda_nonrel_nokern.inp 
   PID         : 12490
   Input dir   : /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar
   Work dir    : /tmp/gomes/gomes/ar.av3z.12490


   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/ar.av3z.mol to working directory as MOLECULE.MOL : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/lda_nonrel_nokern.inp to working directory as DIRAC.INP : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBPOT to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/FRZDN1 to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/test/48.fde_response_autodiff/ar-ar/EMBGRD to working directory : ok
   Copying /Users/gomes/devel/svn/Dirac_embrsp/dirac.x to working directory : 