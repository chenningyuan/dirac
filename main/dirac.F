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
C     
C     Direct Iterative Relativistic All-electron Calculations
C     
C     Original version written by:
C     Trond Saue, University of Oslo, Norway, 1992
C
C     Further contributions by:
C
C     NOTE
C
C     This is an experimental code for the relativistic calculation of
C     atomic and molecular electronic wavefunctions.
C     The author accept no responsibility for the performance of the code
C     or for the correctness of the results.
C
C     The code (in whole or part) is not to be reproduced for
C     further distribution without the written permission of the author.
C
C     If results obtained with this code are published, an appropriate
C     citation would be:
C     T.Saue:"DIRAC(Direct Iterative Relativistic All-electron
C     Calculations), a program for relativistic molecular calculations".
C
C     This note is not to be removed or altered.
C
C*****************************************************************************
C
C     R E V I S I O N S:
C
C     3.0 - Property module:
C            - Expectation values
C
C     2.0 - September 1995: T.Saue
C           Symmetry implemented: D2H and subgroups + time reversal
C           Dynamic allocation implemented using environmental variable DIRWRK
C             (modelled after the ABACUS code)
C           Convergence acceleration DIIS implemented (J.Laerdahl & T.Saue)
C
C
C*****************************************************************************
C
C     PAM - Paul Adrien Maurice
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
C
C     Transfer variables to common blocks
C
#if defined (VAR_MPI)
      MPARID = IMPARID
      NUMNOD = INUMNOD
      MYTID  = IMYTID
#endif
      PARCAL = (IPARCAL .NE. 0)
C
C
      CALL TIMDAT
      CALL DIRCTL(WORK,LWORK)
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
#ifdef MOD_GEO
         IF (MINIMI) THEN
            CALL PAMGEO(WORK,LWORK)
         ELSE
            CALL PAMDRV(WORK,LWORK)
         ENDIF
#else
         CALL PAMDRV(WORK,LWORK)
#endif
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
      IF(DODHF) CALL PSIDHF(WORK,LWORK)
#if defined (MOD_TRA)
C
C*****************************************************************************
C***** M O E L L E R - P L E S S E T  ****************************************
C*****************************************************************************
C
      IF(DOMP2) CALL PSIMP2(WORK,LWORK)
#endif
#if defined (MOD_TRA)
C
C*****************************************************************************
C***** I N D E X  T R A N S F O R M A T I O N*********************************
C*****************************************************************************
C
      IF(DOTRA) CALL PAMTRA(WORK,LWORK)
#endif
#if defined (MOD_CC)
C
C*****************************************************************************
C***** R E L C C S D  M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCCM) CALL PAMCCM(WORK,LWORK)
#endif
#if defined (MOD_CI)
C
C*****************************************************************************
C***** D I R R C I   M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCIM) CALL PAMCIM(WORK,LWORK)
C
#endif

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
      NUMNOD = INUMNOD
      MYTID  = IMYTID
      PARCAL = (IPARCAL .NE. 0)
 100  CONTINUE
#endif
#if defined (VAR_PVM)
      CALL PVMFMYTID(MYTID)
      IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C
      CALL PVMFPARENT(MPARID)
      IF (MPARID .LT. 0) CALL PFAIL(MPARID)
C
 100  CONTINUE
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
      CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
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
C        Start calculating integrals.
C
         CALL HERNOD(WORK,LWORK)
      ELSE IF (ITASK .EQ. 2) THEN
C
C        Start calculating integrals.
C
#ifdef MOD_TRA
          CALL TRANOD(WORK,LWORK)
#else
          CALL QUIT('*** MOLTRA is not in this version')
#endif
      ELSE IF (ITASK .EQ. 3) THEN
C
C        Start coupled cluster calculation.
C
#ifdef MOD_CC
         CALL CCNOD(WORK,LWORK)
#else
          CALL QUIT('*** RELCCSD is not in this version')
#endif
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

