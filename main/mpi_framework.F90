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
PARCAL     = .FALSE.
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
   print *,'init_init' 
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
   print *,'mpi_comm_rank' 
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYTID, IERR)
  IF( IERR .ne. 0 ) THEN
     WRITE(LUPRI,'(/A)') '  ** ERROR in MPIXINIT ** call to MPI_COMM_RANK failed. '
     CALL Abend2( 'Error calling MPI_COMM_RANK' )
  END IF
  print *,'mpi_comm_size' 
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC, IERR)
  IF( IERR .ne. 0 ) THEN
     WRITE(LUPRI,'(/A)') '  ** ERROR in MPIXINIT ** call to MPI_COMM_SIZE failed. '
     CALL Abend2( 'Error calling MPI_COMM_SIZE' )
  END IF
!
! we count the number of slave nodes, subtract the master
!
  NUMNOD = NUMPROC - 1
!
  PARCAL = .TRUE.
!
! synchronize loosely, so the nodes may print to output before the master
! goes on printing header ...
!
   print *,'mpi_barrier' 
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
END IF ! IPAREX
!
! distribute environment variables
!
! ITASK = 0: "BASDIR"
! ITASK = 1: "DIRWRK"
! ITASK = 2: "DIRNOD"
!
  CALL DIST_ENV(0) !later in the DIRAC run???
  CALL DIST_ENV(1)
  CALL DIST_ENV(2)
!
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
#include "dcbgen.h"
INTEGER :: IERR
!
IERR = 0
IF( PARCAL )THEN
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
END IF ! PARCAL == .TRUE.
END SUBROUTINE MPIXFINALIZE
!
!      /* Deck MPIXINIT */
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
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV ('BASDIR',6,BASDIR,LBASDIR,IERR)
#else
      CALL GETENV('BASDIR',BASDIR_SCR)
#endif
      LEN_LOCAL = LBASDIR
10    CONTINUE
      IF ( LEN_LOCAL .gt. 0 ) THEN
          IF ( BASDIR_SCR(LEN_LOCAL:LEN_LOCAL) .EQ. ' ' ) THEN
              LEN_LOCAL = LEN_LOCAL - 1
              GOTO 10
          END IF
      END IF
!
!     If string is empty add current directory
!
      IF ( LEN_LOCAL .eq. 0 ) THEN
          BASDIR_SCR(1:1) = '.'
          LEN_LOCAL   = 1
      END IF
   END IF ! MYTID .eq. MPARID
! 
!  update slaves
!
   IF( NUMNOD .gt. 0 )THEN
      CALL MPI_BCAST(LEN_LOCAL,1,my_MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(BASDIR_SCR,LEN_LOCAL,MPI_CHARACTER,MPARID,MPI_COMM_WORLD,IERR)
   END IF
   BASDIR(1:LEN_LOCAL) = BASDIR_SCR(1:LEN_LOCAL)
   IF (MYTID .eq. MPARID) WRITE(LUPRI,'(2A)') ' Basis set dir.    : ',BASDIR(1:LEN_LOCAL)
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
#else
! Empty files are not allowed, at least no by ifort /ulfek
SUBROUTINE MPI_STUB
END SUBROUTINE 

#endif /* VAR_MPI */
