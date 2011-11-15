!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      PROGRAM RELCCSD
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description-----------------------------------------------
C
C    Double group symmetry adapted unrestricted coupled cluster code
C    CRAY version uses dynamic memory allocation and variable file buffers
C
C---------------Routines called-------------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      PARAMETER (MAXCOR = 60 000 000)
#include "files.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC(MAXCOR)
      REAL*8 SEC,SEC0
      LOGICAL CA
      CHARACTER*10 INTFCE
C
C---------------Executable code--------------------------------------
C
      CALL CPUUSED(SEC0)
      IN = 5
      IW = 6
      INTFCE = 'MOLFDIR   '
      CA = .FALSE.
#if defined (VAR_MPI)
C
C     MPI initialization
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYPROC,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NMPROC,IERR)
      MASTER = 0
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
#if defined (CC_ALONE)
      CALL CCMAIN(MAXCOR,IN,IW,CA,INTFCE)
#else
      CALL CCMAIN(CC,MAXCOR,IN,IW,CA,INTFCE)
#endif
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
#if defined (VAR_MPI)
C
C     MPI finalization
C
      CALL MPI_FINALIZE(IERR)
#endif
C
 1000 FORMAT (//' CPU time (seconds) used in RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used before RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of RELCCSD Run ---')
C
      END
