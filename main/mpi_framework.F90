!
!...   Copyright (c) 2004 by the authors of Dirac (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   Dirac, a relativistic ab initio electronic structure program,
!...   release DIRAC04.0 (2004),
!...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
!...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
!...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
!...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
!...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
!...   (http://dirac.chem.sdu.dk).
!...
!...   For a suitable BibTEX entry, see:
!...      http://dirac.chem.sdu.dk/doc/reference.shtml
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the
!...   research group of the licence holder may not take any part of Dirac,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For questions concerning this copyright write to:
!...      dirac-admin@dirac.chem.sdu.dk
!...
!...   For information on how to get a licence see:
!...      http://dirac.chem.sdu.dk
!
! VERSION : $Revision$
! DATE    : $Date$
! FILE    : mpi_framework.F90
!
!      /* Deck MPIXINIT */
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (VAR_MPI) 
SUBROUTINE MPIXINIT()
!**********************************************************************************
!
!          called from DIRAC main routine - initializes MPI environment
!
!          written by Stefan Knecht, Duesseldorf 2008
!
!**********************************************************************************
!
use os_utils
#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
#include "infpar.h"
#include "dcbgen.h"
INTEGER             :: IPAREX, IPAREX_SCR, IERR, NUMPROC
CHARACTER (LEN=6)   :: DIRPAR
!
! initialize
!
IERR       = 0
IPAREX_SCR = 6
IPAREX     = 1
MYTID      = 0
NUMNOD     = 0
NUMPROC    = 0
MPARID     = 0
!
!
! check for parallel execution
! default: parallel if DIRPAR not defined
!
IPAREX = GET_ENVIRONMENT_INTEGER('DIRPAR',1)
!
! initialize MPI_COMM_WORLD
!
IF( IPAREX .ne. 0 )THEN
  CALL MPI_INIT( IERR )
  IF( IERR .ne. 0 ) THEN
     WRITE(LUPRI,'(/A)') '  ** ERROR in MPIXINIT ** MPI initialization failed. '
     CALL Abend2( 'Error calling MPI_INIT' )
  END IF
!
! TODO: MPI errorhandler 
! 
! SK: MPI librarys provide much better MPI errorhandlers than we can do.
!     Therefore, I did not implement it by now.
!     12-03-2008.
!     any suggestions?
!
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYTID, IERR)
  IF( IERR .ne. 0 ) THEN
     WRITE(LUPRI,'(/A)') '  ** ERROR in MPIXINIT ** call to MPI_COMM_RANK failed. '
     CALL Abend2( 'Error calling MPI_COMM_RANK' )
  END IF
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC, IERR)
  IF( IERR .ne. 0 ) THEN
     WRITE(LUPRI,'(/A)') '  ** ERROR in MPIXINIT ** call to MPI_COMM_SIZE failed. '
     CALL Abend2( 'Error calling MPI_COMM_SIZE' )
  END IF
!
! we count the number of slave nodes, subtract the master
!
! at present this subtraction is still done in the subroutine DIRAC (dirac.F)
! NUMNOD = NUMPROC - 1
  NUMNOD = NUMPROC
!
  PARCAL = .TRUE.
!
! synchronize loosely, so the nodes may print to output before the master
! goes on printing header ...
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
#if defined (MPE)
!
! initialize MPE if precompiler option MPE is set
!
  CALL MPE_INIT_LOG()
  IF ( MYTID .eq. MPARID ) THEN
     CALL MPE_DESCRIBE_STATE( 1, 2, "Init"          ,  "black")
     CALL MPE_DESCRIBE_STATE( 3, 4, "Wait"          ,    "red")
     CALL MPE_DESCRIBE_STATE( 5, 6, "Send/recv task", "yellow")
     CALL MPE_DESCRIBE_STATE( 7, 8, "Send/recv data", "green" )
     CALL MPE_DESCRIBE_STATE( 9,10, "Calc"          , "blue"  )
  END IF
  CALL MPE_START_LOG()
#endif  /* MPE */
!
! distribute environment variables
!
! ITASK = 0: "BASDIR"
! ITASK = 1: "DIRWRK"
! ITASK = 2: "DIRNOD"
!
  CALL DIST_ENV(1)
  CALL DIST_ENV(2)
ELSE
	PARCAL = .FALSE.
END IF ! IPAREX
END SUBROUTINE MPIXINIT

!
!      /* Deck MPIXFINALIZE */
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE MPIXFINALIZE()
!**********************************************************************************
!
!          called from DIRAC main routine - finalize MPI
!
!          written by Stefan Knecht, Duesseldorf 2008
!
!**********************************************************************************
!
#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
#include "infpar.h"
INTEGER :: IERR
!
IERR = 0
#if defined (MPE)
! stop logging
CALL MPE_STOP_LOG()
CALL MPE_FINISH_LOG("mpilog")
#endif
!
! finalize MPI
!
CALL MPI_Finalize( IERR )
IF( IERR .ne. 0 ) THEN
   WRITE(LUPRI,'(/A)') '  ** ERROR in MPIXFINALIZE ** call to MPI_FINALIZE failed. '
   CALL Abend2( 'Error calling MPI_FINALIZE' )
END IF
END SUBROUTINE MPIXFINALIZE
!
!      /* Deck DIST_ENV */
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE DIST_ENV(ITASK)
!**********************************************************************************
!
!          called from DIRAC main routine - synchronizes environment variables 
!                                           between MASTER and SLAVES.
!
!          ITASK = 0: "BASDIR"
!          ITASK = 1: "DIRWRK"
!          ITASK = 2: "DIRNOD"
!
!          OUTPUT: character strings on COMMON/GNCPAM/ (dcbgen.h)
!
!          written by Stefan Knecht, Duesseldorf 2008
!
!**********************************************************************************
!
#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
#include "infpar.h"
#include "dcbgen.h"
INTEGER                 :: LEN_LOCAL, IERR, LDIRWRK,LDIRNOD
PARAMETER (LDIRWRK=255,LDIRNOD=255)
! local scratch
CHARACTER (LEN=LBASDIR) :: BASDIR_SCR
CHARACTER (LEN=LDIRWRK) :: CDIRWRK_SCR
CHARACTER (LEN=LDIRNOD) :: CDIRNOD_SCR
!
LEN_LOCAL = 0
IERR      = 0
!
IF( ITASK .eq. 0 )THEN
   IF( MYTID .eq. MPARID )THEN
!tsaue: Adding initialization until I have checked what happend for SYS_T3D and SYS_T90....
      LEN_LOCAL = LBASDIR
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV ('BASDIR',6,BASDIR,LBASDIR,IERR)
#else
      CALL GETENV('BASDIR',BASDIR_SCR)
#endif
      LEN_LOCAL = LEN_TRIM(BASDIR_SCR)
!
!     If string is empty add current directory
!
      IF ( LEN_LOCAL .eq. 0 ) THEN
          BASDIR_SCR(1:1) = '.'
          LEN_LOCAL   = 1
      END IF
! add a trailing blank to avoid problems
      BASDIR_SCR(LEN_LOCAL+1:LEN_LOCAL+1) = ' '
   END IF ! MYTID .eq. MPARID
! 
!  update slaves
!
   IF( NUMNOD .gt. 0 )THEN
      CALL MPI_BCAST(LEN_LOCAL,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(BASDIR_SCR,LEN_LOCAL,MPI_CHARACTER,MPARID,MPI_COMM_WORLD,IERR)
   END IF
   BASDIR(1:LEN_LOCAL) = BASDIR_SCR(1:LEN_LOCAL)
!  IF (MYTID .eq. MPARID) WRITE(LUPRI,'(2A)') ' Basis set dir.    : ',BASDIR(1:LEN_LOCAL)
!
ELSE IF ( ITASK .eq. 1 ) THEN ! ITASK == 1: DIRWRK
!
   IF( MYTID .eq. MPARID )THEN
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV ('DIRWRK',6,CDIRWRK_SCR,LDIRWRK,IERR)
#else
      CALL GETENV("DIRWRK",CDIRWRK_SCR)
#endif
      LEN_LOCAL = LDIRWRK
20    CONTINUE
      IF ( LEN_LOCAL .gt. 0 ) THEN
          IF ( CDIRWRK_SCR(LEN_LOCAL:LEN_LOCAL) .eq. ' ' ) THEN
              LEN_LOCAL = LEN_LOCAL - 1
              GOTO 20
          END IF
      END IF
!
!     if string is empty add a blank value
!
      IF ( LEN_LOCAL .eq. 0 ) THEN
          CDIRWRK_SCR(1:1) = ' '
          LEN_LOCAL   = 1
      END IF
   END IF ! MYTID .eq. MPARID
! 
!  update slaves
!
   IF( NUMNOD .gt. 0 )THEN
      CALL MPI_BCAST(LEN_LOCAL,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CDIRWRK_SCR,LEN_LOCAL,MPI_CHARACTER,MPARID,MPI_COMM_WORLD,IERR)
   END IF
! ulfek: CDIRWRK not used below?
! ulfek: comment out   CDIRWRK(1:LEN_LOCAL) = CDIRWRK_SCR(1:LEN_LOCAL)
!
!  WRITE(LUPRI,'(2A,2X,I4)') ' string DIRWRK (MYTID) : ',CDIRWRK(1:LEN_LOCAL),MYTID
!
ELSE IF ( ITASK .eq. 2 ) THEN ! ITASK == 2: DIRNOD
!
   IF( MYTID .eq. MPARID )THEN
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV ('DIRNOD',6,CDIRNOD_SCR,LDIRNOD,IERR)
#else
      CALL GETENV("DIRNOD",CDIRNOD_SCR)
#endif
      LEN_LOCAL = LDIRNOD
30    CONTINUE
      IF ( LEN_LOCAL .gt. 0 ) THEN
          IF ( CDIRNOD_SCR(LEN_LOCAL:LEN_LOCAL) .eq. ' ' ) THEN
              LEN_LOCAL = LEN_LOCAL - 1
              GOTO 30
          END IF
      END IF
!
!     if string is empty add a blank value
!
      IF ( LEN_LOCAL .eq. 0 ) THEN
          CDIRNOD_SCR(1:1) = ' '
          LEN_LOCAL   = 1
      END IF
   END IF ! MYTID .eq. MPARID
! 
!  update slaves
!
   IF( NUMNOD .gt. 0 )THEN
      CALL MPI_BCAST(LEN_LOCAL,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CDIRNOD_SCR,LEN_LOCAL,MPI_CHARACTER,MPARID,MPI_COMM_WORLD,IERR)
   END IF
! ulfek:   CDIRNOD(1:LEN_LOCAL) = CDIRNOD_SCR(1:LEN_LOCAL)
!  WRITE(LUPRI,'(2A,2X,I4)') ' string DIRNOD (MYTID) : ',CDIRNOD(1:LEN_LOCAL),MYTID
!
END IF ! ITASK
!
END SUBROUTINE DIST_ENV

!
!      /* Deck SYNC_NODES */
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE SYNC_NODES(IVARIABLE,ITASK)
!**********************************************************************************
!
!          called from DIRAC main routine - write DIRAC information in ordered
!          way - smarter solutions are welcome.
!
!          written by Stefan Knecht, Duesseldorf 2008
!
!          ITASK == 0: print memory allocation; IVARIABLE = memsize (from
!                      calling routine)
!          ITASK == 1: call timer and print timing results
!
!**********************************************************************************
!
#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
#include "infpar.h"
  INTEGER :: IERR, IJ, LHOST, ILEN
  PARAMETER (ILEN=72)
  CHARACTER (LEN=ILEN) HOST
!
  IERR  = 0
  LHOST = 0
! ordered start-up info print
  IF( ITASK .eq. 0 )THEN
    CALL MPI_GET_PROCESSOR_NAME(HOST,LHOST,IERR)
    IF( MYTID .eq. MPARID ) THEN
       print '(a,a,a,i0,a,i0,a)', 'DIRAC master    (',HOST(1:LHOST), &
             ') allocating ',IVARIABLE,' words (',IVARIABLE/2**17,   &
             ' MB) of memory'
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
    DO IJ = 1, (NUMNOD - 1)
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       IF( MYTID .eq. IJ )  THEN
          print '(a,i4,1x,a,a,a,i0,a,i0,a)', 'DIRAC node ',mytid, &
          '(',HOST(1:LHOST),') allocating ',IVARIABLE,' words (', &
          IVARIABLE/2**17,' MB) of memory'
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
    END DO
!   ordered last printing
  ELSE IF( ITASK .eq. 1 )THEN
    DO IJ = 1, NUMNOD
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IF( MYTID .eq. IJ ) THEN
         CALL STOPWATCH(0,MYTID)
         CALL FLSHFO(LUPRI)
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    END DO
  ELSE IF( ITASK .eq. 2 )THEN
    CALL MPI_GET_PROCESSOR_NAME(HOST,LHOST,IERR)
    IF( MYTID .eq. MPARID ) THEN
       print '(a,a,a,i0,a,i0,a)', 'DIRAC master    (',HOST(1:LHOST), &
             ') will attempt to use at most ',IVARIABLE,' words (',IVARIABLE/2**17,&
             ' MB) of memory in dynamic allocations'
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
    DO IJ = 1, (NUMNOD - 1)
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       IF( MYTID .eq. IJ )  THEN
          print '(a,i4,1x,a,a,a,i0,a,i0,a)', 'DIRAC node ',mytid, &
          '(',HOST(1:LHOST),') will attempt to use at most ',IVARIABLE,' words (', &
          IVARIABLE/2**17,' MB) of memory in dynamic allocations'
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
    END DO
!
  END IF
!
END SUBROUTINE SYNC_NODES
!
#else
! Empty files are not allowed, at least not by ifort /ulfek
SUBROUTINE MPI_STUB
END SUBROUTINE 

#endif /* VAR_MPI */
