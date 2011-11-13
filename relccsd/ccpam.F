C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
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
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C
C
C VERSION : $Revision: 1.24.4.2 $
C DATE    : $Date: 2004/09/01 08:17:54 $
C FILE    : ccmain.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAMCCM (CC,MAXCOR)
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
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbham.h"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "ccpar.inc"
#include "files.inc"
#include "complex.inc"
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
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      CALL ASNUNIT(LUCMD,'-Yon',IERR)
      IF (IERR.NE.0)  WRITE (IW,'(/A,I4/)') ' Assign error ',IERR
#endif
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      IN = LUCMD
      IW = LUPRI
      INTFCE = 'DIRAC     '
      CA = NZ.GE.2.AND..NOT.(SPINFR.OR.LEVYLE)
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
C
C     Summon the slaves, who are waiting in the general menu routine.
C
      CALL CCMSTI
C
C     Tell them the correct arithmetic (NZ is not initialized on the slaves
C     if you do not run HF !)
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         CALL MPI_BCAST(CA,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
      END IF
#endif
C
C     Enter the generic CCMAIN routine (master becomes now primus inter parens)
C
      CALL CCMAIN(CC,MAXCOR,IN,IW,CA,INTFCE)
C
C     Return the slaves to the general menu routine.
C
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
#include "implicit.h"
#include "maxorb.h"
#include "dcbgen.h"
C
C     ITASK = 3 for parallel RELCCSD
C
      IF (PARCAL) CALL DIRAC_PARCTL( 3 ) 
C
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
#include "implicit.h"
#include "maxorb.h"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "dcbgen.h"
#include "infpar.h"
C
      IF (.NOT.PARCAL) RETURN
C
C     Find the slaves and release them
C
      NTEST = -1
#if defined (VAR_MPI)
      DO I = 1,NUMNOD
         CALL MPI_IRECV(NODE,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,IREQ,IERR)
         CALL MPI_WAIT(IREQ,ISTAT,IERR)
         CALL MPI_SEND(NTEST,1,MPI_INTEGER,NODE,30,
     &                  MPI_COMM_WORLD,IERR)
      ENDDO
#endif
C
C     Let the slaves sleep (if relevant for this operating system)
C
      CALL DIRAC_PARCTL( 0 )
      RETURN
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CCNOD */
      SUBROUTINE CCNOD(CC,MAXCOR)
C
C     Written by Luuk Visscher, August 1997
C     The slave RELCCSD routine (DIRAC interface)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "maxorb.h"
#include "dcbham.h"
#include "infpar.h"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "priunit.h"
#include "dgroup.h"
#include "ccpar.inc"
#include "files.inc"
      CHARACTER*10 INTFCE
      LOGICAL CA,CCEND
C
      DIMENSION CC(MAXCOR)
C
C     arrange for the MPI stuff and correct node number
C     to the total number of running invocations.
C
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif

      INTFCE = 'DIRAC     '
C
C     Get the correct arithmetic (NZ is not initialized on slaves if you do not run HF !)
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         CALL MPI_BCAST(CA,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
      END IF
#endif
C
C     create a node-unique filename as output file. Important on
C     shared file systems. Otherwise all the output gets mingled in one
C     file
C
      CCBASF="RELCCSD.OUT"
      
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (CCFILN,'(A11,A1,I1)') CCBASF,'.',MYPROC
         LCCFIL=13
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (CCFILN,'(A11,A1,I2)') CCBASF,'.',MYPROC
         LCCFIL=14
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (CCFILN,'(A11,A1,I3)') CCBASF,'.',MYPROC
         LCCFIL=15
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND CCNOD.F MODULE")
      ENDIF
#else
      CCFILN=CCBASF
      LCCFIL=11
#endif
C
C     open the local input file and the node specific output file.
C     The associated file handle MSLVOUT is transferred to CCMAIN.
C     Every access to the local stdout handle then automatically writes
C     to the corresponding output file.
C
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      OPEN(MSLVOUT,FILE = CCFILN(1:LCCFIL))
      CALL CCMAIN(CC,MAXCOR,LUCMD,MSLVOUT,CA,INTFCE)
      CLOSE (MSLVOUT,STATUS='KEEP')
      CLOSE (LUCMD,STATUS='KEEP')
C
C     The slave is finished, wait for the master
C
      CCEND = .FALSE.
#if defined (VAR_MPI)
      CALL MPI_ISEND(MYTID,1,MPI_INTEGER,MPARID,20,
     &              MPI_COMM_WORLD,IREQ,IERR)
C
    1    CONTINUE
         CALL SLEEP(1)
         CALL MPI_IPROBE(MPARID,30,MPI_COMM_WORLD,CCEND,ISTAT,IERR)
         IF (.NOT.CCEND) GOTO 1
C
      CALL MPI_RECV(NTEST,1,MPI_INTEGER,MPARID,30,
     &               MPI_COMM_WORLD,ISTAT,IERR)
#endif
C
      RETURN
      END
#endif    /* ifdef VAR_MPI */
