C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMCCM (CC,MAXCOR)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description-----------------------------------------------
C
C    Double group symmetry adapted unrestricted coupled cluster code
C    This version is called from within DIRAC
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
      INCLUDE "files.inc"
      INCLUDE "inpt.inc"
      INCLUDE "complex.inc"
#include <priunit.h>
#include <dgroup.h>
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC(MAXCOR)
      REAL*8 SEC,SEC0
C
C---------------Executable code--------------------------------------
C
      CALL CPUUSED(SEC0)
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      IN = LUCMD
      IW = LUPRI
      INTERFACE = 'DIRAC     '
      CARITH = NZ.GE.2
      WRITE (IW,'(//80A1)') ('#',I=1,80)
C
      CALL CCHEADER
C
      CALL CCMSTI
      CALL CCMAIN(CC,MAXCOR)
      CALL CCMSTE
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
      CLOSE(LUCMD)
C
 1000 FORMAT (//' CPU time (seconds) used in RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used before RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of RELCCSD Run ---')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCMSTI */
      SUBROUTINE CCMSTI
C
C     Written by Luuk Visscher, August 1997
C     Call the slaves if we run in parallel
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <mxorb.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#else
#include <fpvm3h.h>
#endif
#include <gnrinf.h>
#include <infpar.h>
C
      IF (.NOT.PARCAL) RETURN
C
      ITASK = 3
#if defined (VAR_MPI)
      CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
      CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCMSTE */
      SUBROUTINE CCMSTE
C
C     Written by Luuk Visscher, August 1997
C     Release the slaves if we run in parallel
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <mxorb.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#else
#include <fpvm3h.h>
#endif
#include <gnrinf.h>
#include <infpar.h>
C
      IF (.NOT.PARCAL) RETURN
C
C     Find the slaves and release them
C
      NTEST = -1
      DO I = 1,NODES
#if defined (VAR_MPI)
         CALL MPI_IRECV(NODE,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,IREQ,IERR)
         CALL MPI_WAIT(IREQ,ISTAT,IERR)
         CALL MPI_SEND(NTEST,1,MPI_INTEGER,NODE,30,
     &                  MPI_COMM_WORLD,IERR)
#else
#endif
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCNOD */
      SUBROUTINE CCNOD(WORK,LWORK)
C
C     Written by Luuk Visscher, August 1997
C     The slave RELCCSD routine
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <mxorb.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#else
#include <fpvm3h.h>
#endif
#include <gnrinf.h>
#include <infpar.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL CCEND
C
      CALL QENTER('CCNOD')
#include <memint.h>
      KFRSAV = KFREE
#if defined (VAR_PVM)
      CALL PVMFMYTID(MYTID)
      IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C
      CALL PVMFPARENT(MPARID)
      IF (MPARID .LT. 0) CALL PFAIL(MPARID)
C
#endif
C
C     At the moment this is a dummy routine, letting the slaves
C     sleep so that they do not consume CPU-time idling for tasks.
C
      CCEND = .FALSE.
#if defined (VAR_MPI)
      CALL MPI_ISEND(MYTID,1,MPI_INTEGER,MPARID,20,
     &              MPI_COMM_WORLD,IREQ,IERR)
C
    1    CONTINUE
         CALL SLEEP(5)
         CALL MPI_IPROBE(MPARID,30,MPI_COMM_WORLD,CCEND,ISTAT,IERR)
         IF (.NOT.CCEND) GOTO 1
C
      CALL MPI_RECV(NTEST,1,MPI_INTEGER,MPARID,30,
     &               MPI_COMM_WORLD,ISTAT,IERR)
#else
#endif
C
      CALL MEMREL('CCNOD',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('CCNOD')
C
      RETURN
C
      END
