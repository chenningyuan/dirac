C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccmain.F
C SOURCE  : $Source$
C
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
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      CALL ASNUNIT(LUCMD,'-Yon',IERR)
      IF (IERR.NE.0)  WRITE (IW,'(/A,I4/)') ' Assign error ',IERR
#endif
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
