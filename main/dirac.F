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
C FILE    : dirac.F
C SOURCE  : $Source$
C
C*****************************************************************************
      SUBROUTINE DIRAC(IPARCAL,IMYTID,IMPARID,INUMNOD,WORK,LWORK)
#include <implicit.h>
C
#include <dcbgen.h>
C
#if defined (VAR_MPI)
#include <mxorb.h>
#include <infpar.h>
#include <infmpi.h>
#endif
C
      DIMENSION WORK(LWORK)
C
C
      PARCAL = (IPARCAL .NE. 0)
#if defined (VAR_MPI)
C
C     Transfer variables to common blocks
C
      MPARID = IMPARID
C We count the number of slave nodes, subtract the master
      NUMNOD = INUMNOD - 1
      MYTID  = IMYTID
C
      PARCAL = PARCAL.AND.NUMNOD.ne.0
C
C     Let the slaves sleep for now (where implemented) ...
C        ( ITASK = 0 )
C
      IF (PARCAL) CALL DIRAC_PARCTL( 0 )
#endif
C
C     Set basis set directory
C
      CALL GETBASDIR()
C
C     Execute Dirac
C
      CALL TIMDAT
      CALL DIRCTL(WORK,LWORK)
C
C     Send finish signal to slaves
C
      IF (PARCAL) CALL DIRAC_PARCTL( -2 )
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirctl */
      SUBROUTINE DIRCTL(WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
      DIMENSION WORK(LWORK)
C
      CALL QENTER('DIRCTL')
C
C***********************************************************************
C*****   I  N  P  U  T    S  E  C  T  I  O  N   ************************
C***********************************************************************
C
      CALL PAMINP(WORK,LWORK)
C
      IF ( INPTES ) RETURN
C
C*****************************************************************************
C*****  S E T U P  -  M O D U L E  *******************************************
C*****************************************************************************
C
      CALL PAMSET(WORK,LWORK)
C
C*****************************************************************************
C**** M A I N    S E C T I O N  **********************************************
C*****************************************************************************
C
      IF(.NOT.DOHRM) THEN
         IF (MINIMI) THEN
            CALL PAMGEO(WORK,LWORK)
         ELSE
            CALL PAMDRV(WORK,LWORK)
         ENDIF
      END IF
C
C******************************************************************************
C*****  T  E  R  M  I  N  A  T  I  O  N  **************************************
C******************************************************************************
C
      CALL QEXIT('DIRCTL')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamdrv */
      SUBROUTINE PAMDRV(WORK,LWORK)
C*****************************************************************************
C
C     Driver for DIRAC - main part
C
C     Written by T. Saue - November 1994
C     Last revision : Sep 11 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PAMDRV')
C
C******************************************************************************
C*****  G E T    W A V E    F U N C T I O N   *********************************
C******************************************************************************
C
      IF(DOPIN) CALL PCMOIN(WORK,LWORK)
      IF(DOPSI .OR. DOTRA) CALL PAMPSI(WORK,LWORK)
      IF(DOPUT) CALL PCMOUT(WORK,LWORK)
C
C*****************************************************************************
C*****  A N A L Y S I S    M O D U L E  **************************************
C*****************************************************************************
C
      IF(DOANA) CALL PAMANA(WORK,LWORK)
C
C
C*****************************************************************************
C*****  R E S P O N S E    M O D U L E  **************************************
C*****************************************************************************
C
      IF(DOPRP) CALL PAMPRP(WORK,LWORK)
C
C*****************************************************************************
C***** H E I N Z *************************************************************
C*****************************************************************************
C
      IF(DOHEINZ) CALL HEINZ(WORK,LWORK)
C
 666  CONTINUE
      CALL QEXIT('PAMDRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampsi */
      SUBROUTINE PAMPSI(WORK,LWORK)
C*****************************************************************************
C
C     Wave function module
C
C     Written by T.Saue April 1996
C     Last revision: April 27 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
C
#include <dcbgen.h>
#include <dcbpsi.h>
C
      CALL QENTER('PAMPSI')
C
C
C*****************************************************************************
C***** D I R A C - H A R T R E E - F O C K ***********************************
C*****************************************************************************
C
      IF(DODHF .AND. DOPSI) CALL PSIDHF(WORK,LWORK)
      IF(DOMVO) CALL MVOGEN(WORK,LWORK)
C
C*****************************************************************************
C***** R E S O L V E    O P E N - S H E L L  S T A T E S *********************
C*****************************************************************************
C
      IF(DORES) CALL RESOLV(WORK,LWORK)
C
C*****************************************************************************
C***** M O E L L E R - P L E S S E T  ****************************************
C*****************************************************************************
C
      IF(DOMP2) CALL PSIMP2(WORK,LWORK)
C
C*****************************************************************************
C***** I N D E X  T R A N S F O R M A T I O N*********************************
C*****************************************************************************
C
      IF(DOTRA) CALL PAMTRA(WORK,LWORK)
C
C*****************************************************************************
C***** R E L C C S D  M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCCM) CALL PAMCCM(WORK,LWORK)
C
C*****************************************************************************
C***** D I R R C I   M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCIM) CALL PAMCIM(WORK,LWORK)
C
C
      CALL QEXIT('PAMPSI')
      RETURN
C
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnod */
#if defined (VAR_MPI) || defined (VAR_MPI)
      SUBROUTINE DIRNOD(IPARCAL,IMYTID,IMPARID,INUMNOD,WORK,MWORK)
#endif
#if defined (VAR_PVM)
      PROGRAM DIRNOD
#endif
C*****************************************************************************
C
C     Menu file for slave nodes. Check whether the program needs a slave
C     and determine the type of work to be done.
C
C     Written by L.Visscher July 1997
C     Last revision: July 29 1996
C
C*****************************************************************************
#include <implicit.h>
#include <dcbdhf.h>
#include <dcbgen.h>
#include <priunit.h>
C
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
      INTEGER ISTATUS(MPI_STATUS_SIZE)
      LOGICAL NEWTASK
#else
#include <fpvm3h.h>
#endif
C
#if defined (VAR_PVM)
      PARAMETER (MWORK = 6 000 000)
#endif
#include  <mxorb.h>
#include <infpar.h>
C
      DIMENSION WORK(MWORK)
C
      CALL QENTER('DIRNOD')
      LWORK = MWORK
#include <memint.h>
C
#if defined (VAR_MPI)
C
C     Transfer variables to common blocks
C
      MPARID = IMPARID
C We count the number of slave nodes, subtract the master
      NUMNOD = INUMNOD - 1
      MYTID  = IMYTID
      PARCAL = (IPARCAL .NE. 0)
#endif
#if defined (VAR_PVM)
      CALL PVMFMYTID(MYTID)
      IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C
      CALL PVMFPARENT(MPARID)
      IF (MPARID .LT. 0) CALL PFAIL(MPARID)
#endif
C
C     Set basis set directory for this slave
C     (it may be needed for reading the .mol file)
C
      CALL GETBASDIR()
C
chj   write (0,*) 'Calling mpi_iprobe for node',mytid
      call flush(0)
 100  CONTINUE
#if defined (VAR_MPI)
      CALL MPI_IPROBE(MPARID,27,MPI_COMM_WORLD,
     &                NEWTASK,STATUS_MPI,IERR)
      IF (NEWTASK) THEN
chj      write (0,*) 'mpi_iprobe true for node',mytid
chj      write (0,*) 'istatus(source,tag,error):',
chj  &     ( istatus(i),i=1,3)
         CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPARID,27,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
chj      write (0,*) ' itask, node',itask,mytid
chj      write (0,*) 'istatus(source,tag,error):',
chj  &     ( istatus(i),i=1,3)
      ELSE
         CALL SLEEP(1)
         GO TO 100
      END IF
#endif
C
C     Start timing.
C
      CALL GETTIM(CPU1,WALL1)
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
#if defined (VAR_MPI)
chj   CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
      CALL PVMFUNPACK(INTEGER4,ITASK,1,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
#endif
      IF (ITASK .LT. 0) THEN
C
C        Global end of program : go up for MPI_FINALIZE
C
         GOTO 999
      ELSE IF (ITASK .EQ. 0) THEN
C
C        Nothing to be done at the moment. Go to sleep.
C
         CONTINUE
      ELSE IF (ITASK .EQ. 1) THEN
C
C        Start calculating integrals for Fock matrices.
C
         CALL HERNOD(WORK,LWORK)
      ELSE IF (ITASK .EQ. 2) THEN
C
C        Start calculating integrals for MOLTRA.
C
          CALL TRANOD(WORK,LWORK)
      ELSE IF (ITASK .EQ. 3) THEN
C
C        Start coupled cluster calculation.
C
         CALL CCNOD(WORK,LWORK)
      ELSE
C
C        Unknown task. Abort execution.
C
         WRITE(LUPRI,'(A,I5)')
     &     '*** ERROR *** DIRNOD: Unknown task ITASK = ',ITASK
         CALL QUIT(' Unknown task in DIRNOD')
      END IF
C
C     We want the slave program to go until there is nothing
C     more to calculate, K.Ruud, May-95
C     But they may rest for a while to start refreshed for their new job.
C     L. Visscher. Jul-1997
C
      GOTO 100
C
#if defined (VAR_PVM)
      CALL PVMFEXIT(INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
#endif
C
 999  CALL QEXIT('DIRNOD')
      END
#endif
C
#ifdef SYS_SX
C
C     Small dummy program that calls the main c-routine, because
C     some linkers are worthless (eg. NEC SX).
C     If you call the linker through the fortran compiler the 'main'
C     program has one name, and if you call the linker through the C
C     compiler it has another name. Marvellous work!
C
      PROGRAM SXMAIN
      CALL C_MAIN
      STOP
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirctl */
      SUBROUTINE GETBASDIR()
C
#include <implicit.h>
#include <priunit.h>
#include <infmpi.h>
#include <dcbgen.h>
C
C     Check environment variable BASDIR
C
#if defined (SYS_T3D) || defined (SYS_T90)
      CALL PXFGETENV ('BASDIR',6,BASDIR,LBASDR,IERR)
#else
      CALL GETENV ('BASDIR',BASDIR)
#endif
      LENBAS = LBASDIR
 10   CONTINUE
      IF ( LENBAS .GT. 0 ) THEN
         IF ( BASDIR(LENBAS:LENBAS) .EQ. ' ' ) THEN
            LENBAS = LENBAS - 1
            GOTO 10
         END IF
      END IF
C
C     If string is empty add current directory 
C
      IF ( LENBAS .EQ. 0 ) THEN
C
         BASDIR(1:1) = '.'
         LENBAS = 1
C
      END IF
C
      IF (MYTID .EQ. 0) WRITE(LUPRI,'(2A)')
     $     ' Basis set dir.    : ',BASDIR(1:LENBAS)
C 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirac_parctl */
      SUBROUTINE DIRAC_PARCTL(ITASK)
C
C 5-July-2000 Hans Joergen Aa. Jensen
C
C Called by master to control parallel calculation.
C Can be used to send slaves to sleep and to wake them up.
C
C ITASK .gt. 0: wake slaves and get hold of them for task ITASK
C ITASK .eq. 0: send slaves to sleep ("0" makes sure they are in menu)
C ITASK .lt. 0: release slaves from current task and send them to sleep
C
#include <implicit.h>
#include <priunit.h>
C

#if defined (VAR_MPI)

#include <infmpi.h>
#include <infpar.h>
      INCLUDE 'mpif.h'
C
      IF (ITASK .EQ. -1) THEN
C        go back from submenu to main menu via MPI_BCAST
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      ELSE
C        slave should be in main menu, ready to receive task
         DO NODE = 1,NUMNOD
chj         write (0,*) 'sending itask',itask,' to node',node
            CALL MPI_SEND(ITASK,1,MPI_INTEGER,NODE,27,
     &                 MPI_COMM_WORLD,IERR)
         END DO
C        should call wait or barrier here
      END IF

#else /* not VAR_MPI */

      CALL QUIT('DIRAC_PARCTL called '//
     &   'but this is a sequential Dirac version !')

#endif /* VAR_MPI */

      RETURN
      END
