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
C VERSION : $Revision: 1.43.2.4 $
C DATE    : $Date: 2004/09/01 08:13:50 $
C FILE    : dirac.F
C
C***********************************************************************
      SUBROUTINE DIRAC(IPARCAL,IMYTID,IMPARID,INUMNOD,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
#include "infpar.h"
C
      DIMENSION WORK(LWORK)
C
C
      PARCAL = (IPARCAL .NE. 0)
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
      CALL TIMDAT
      CALL FLSHFO(LUPRI)
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirctl */
      SUBROUTINE DIRCTL(WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
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
C***********************************************************************
C*****  S E T U P  -  M O D U L E  *************************************
C***********************************************************************
C
      CALL PAMSET(WORK,LWORK)
C
C***********************************************************************
C**** M A I N    S E C T I O N  ****************************************
C***********************************************************************
C
      IF(.NOT.DOHRM) THEN
         IF (OPTIMI) THEN
            CALL PAMGEO(WORK,LWORK)
         ELSE
            CALL PAMDRV(WORK,LWORK)
         ENDIF
      END IF
C
C***********************************************************************
C*****  T  E  R  M  I  N  A  T  I  O  N  *******************************
C***********************************************************************
C
      CALL QEXIT('DIRCTL')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamdrv */
      SUBROUTINE PAMDRV(WORK,LWORK)
C***********************************************************************
C
C     Driver for DIRAC - main part
C
C     Written by T. Saue - November 1994
C     Last revision : Sep 11 1995 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbpsi.h"
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PAMDRV')
C
C***********************************************************************
C*****  G E T    W A V E    F U N C T I O N   **************************
C***********************************************************************
C
      IF(DOPIN) CALL PCMOIN(WORK,LWORK)
      IF(DOPSI .OR. DOTRA) CALL PAMPSI(WORK,LWORK)
      IF(DOACUT) CALL ACMOUT(WORK,LWORK)
      IF(DOPUT) CALL PCMOUT(WORK,LWORK)
C
C***********************************************************************
C*****  A N A L Y S I S    M O D U L E  ********************************
C***********************************************************************
C
      IF(DOANA .AND. .NOT.((DOSCF.AND.DOPSI).OR.(DOKRMC.AND.DOPSI))) 
     &   CALL PAMANA(WORK,LWORK)
C
C
C***********************************************************************
C*****  R E S P O N S E    M O D U L E  ********************************
C***********************************************************************
C
      IF(DOPRP) CALL PAMPRP(WORK,LWORK)
C
 666  CONTINUE
      CALL QEXIT('PAMDRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampsi */
      SUBROUTINE PAMPSI(WORK,LWORK)
C***********************************************************************
C
C     Wave function module
C
C     Written by T.Saue April 1996
C     Revision for MCSCF J. Thyssen - Oct 1998
C     Restructured 2002 by Hans Joergen Aa. Jensen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL   WFCONV
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
C
      CALL QENTER('PAMPSI')
C
C HJAaJ Nov. 2002 TODO: PAMANA only works for closed shell SCF now
C (DHF or DFT), not for other wave functions.
C It ought to be implemented for other wave function types as well.
C TF, Nov. 2003:  Done for KRMC. 
C
C***********************************************************************
C***** D I R A C - H A R T R E E - F O C K *****************************
C***** o r   D I R A C - K O H N - S H A M *****************************
C***********************************************************************
C
      IF(DOSCF .AND. DOPSI) THEN
C
C        Transfer DHF info to common block DCBORB
C        and recalculate DCBORB.
C
         CALL TDHFORB
         CALL SETDC2(0)
c        call setdc2(iprint)
C
C
C        *************************************************
C        *** Call PSIDHF for Roothaan/DIIS iterations  ***
C        *************************************************
C
         IF (MAXITR.GT.0) THEN
            CALL PSIDHF(WFCONV,WORK,LWORK)
         END IF
C
C        *************************************************
C        *** Call PSIOPT for second order optimization ***
C        *************************************************
C
         IF ( .NOT. NOQCDHF .AND. .NOT. DODFT .AND.
     &       (.NOT. WFCONV .AND. MAXITR .GT.0 )) THEN
#ifdef MOD_KRMC
            CALL PSIOPT('SCF  ',WFCONV,WORK,LWORK)
#else
c           CALL QUIT(
c    &   'Second order SCF optimization not included in this version')
c        ... let it stop below ... /hjaaj Feb 2004
#endif
         END IF
         IF(DOANA) CALL PAMANA(WORK,LWORK)
         IF(.NOT.WFCONV) THEN
C        ... DHF not converged;
            WRITE (LUPRI,'(//A)')
     &     ' ERROR, SCF not converged, Dirac stops!'
            CALL QUIT('SCF not converged')
         END IF
      END IF
C
C***********************************************************************
C***** R E S O L V E    O P E N - S H E L L  S T A T E S ***************
C***********************************************************************
C
      IF(DORES.AND.DOPSI) THEN
         IF (.NOT.DOSCF) THEN
C           ... Transfer DHF info to common block DCBORB
C           ... and recalculate DCBORB.
            CALL TDHFORB
            CALL SETDC2(0)
         END IF
         CALL RESOLV(WORK,LWORK)
      END IF
C
C***********************************************************************
C***** M O E L L E R - P L E S S E T ***********************************
C***********************************************************************
C
      IF(DOMP2.AND.DOPSI) THEN
C
C        *************************************************
C        *** Call PSIMP2 for MP2 Energy/wave function  ***
C        *************************************************
C
         CALL PSIMP2(WORK,LWORK)
      END IF
C
C***********************************************************************
C***** M O D I F I E D   V I R T U A L   O R B I T A L S ***************
C***********************************************************************
C
      IF(DOMVO.AND.DOPSI) CALL MVOGEN(WORK,LWORK)
C
C***********************************************************************
C***** K R - M C S C F *************************************************
C***********************************************************************
C      Do MCSCF before CC and CI because user may want
C      to do an MR-CC or MR-CI /hjaaj aug 03
C
      IF(DOKRMC.AND.DOPSI) THEN
#ifdef MOD_KRMC
C
C        Transfer DHF info to common block DCBORB
C        and recalculate DCBORB.
C
         CALL TKRMCORB
         CALL SETDC2(0)
C
C        *************************************************
C        *** Call PSIOPT for KR-MCSCF wave function ******
C        *************************************************
C
         CALL PSIOPT('MCSCF',WFCONV,WORK,LWORK)
         IF(DOANA) CALL PAMANA(WORK,LWORK)
#else
            CALL QUIT(
     &  'Second order KMCSCF optimization not included in this version')
#endif
      END IF
C
C***********************************************************************
C***** I N D E X  T R A N S F O R M A T I O N **************************
C***********************************************************************
C
      IF(DOTRA) CALL PAMTRA(WORK,LWORK)
C
C***********************************************************************
C***** R E L C C S D  M O D U L E **************************************
C***********************************************************************
C
      IF(DOCCM.AND.DOPSI) CALL PAMCCM(WORK,LWORK)
C
C***********************************************************************
C***** D I R R C I   M O D U L E ***************************************
C***********************************************************************
C
      IF(DOCIM.AND.DOPSI) CALL PAMCIM(WORK,LWORK)

C***********************************************************************
C***** L U C I T A  M O D U L E ****************************************
C***********************************************************************
C
      IF(DOLUCT.AND.DOPSI) CALL PAMLUCITA(WORK,LWORK)
C
C
#if defined (UNDEF)
C
C***********************************************************************
C***** C I *************************************************************
C***********************************************************************
C
      IF(DOCI.AND.DOPSI) THEN
C
C        *************************************************
C        *** Call PSICI for CI wave function *************
C        *************************************************
C
         CALL PSICI(WORK,LWORK)
      END IF
#endif
C
      CALL QEXIT('PAMPSI')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tdhforb */
      SUBROUTINE TDHFORB
C***********************************************************************
C
C     Transfer information from common block DCODHF to 
C     common block DCBORB and set common block DCxOPT
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      CALL QENTER('TDHFORB')
C
C
      NASHT = 0
      DO I = 1,NFSYM
         NISH(I) = NISH_DHF(I)
         NASH(I) = NASH_DHF(I)
         NASHT   = NASHT + NASH(I)
         NOCC(I) = NISH(I) + NASH(I)
         NSSH(I) = NSSH_DHF(I)
c         NESH(I) = NESH_DHF(I)
c         NPSH(I) = NPSH_DHF(I)
         NFRO(I) = NFRO_DHF(I)
         NELEC(I) = NELEC_DHF(I)
      END DO
      NELECT = NELECT_DHF
C
C     Transfer open shell average-of-configuration to GAS space notation
C     (for e.g. ESR GOSCI) /July 2002 hjaaj
C
      NAELEC = NAELEC_DHF
      NGAS   = NOPEN
      MAXMK2 = 2*MIN(NAELEC,NASHT) - NAELEC
      MINMK2 = NAELEC-2*MIN(NAELEC,NASHT)
      MK2REF = 0
      MK2DEL = MAXMK2

      CALL IZERO(NGSHT,MXMCGAS)
      DO J = 1, NGAS
         DO I = 1,NFSYM
            NGSH(I,J) = NACSH(I,J)
            NGSHT(J)  = NGSHT(J) + NGSH(I,J)
         END DO
      END DO
C
C     initialize min and max accumlated no. of electrons 
C     in each gas space
C
      XELEC = D0
      DO J = 1, NGAS
         XELEC = XELEC + D2*DF(J)*NGSHT(J)
         NGSP(1,J) = NINT(XELEC)
         NGSP(2,J) = NINT(XELEC)
      END DO
C
C     Set common block DCIOPT
C
      IPROPT = IPRSCF
      IOPT_ITRINT(1) = ITRINT(1)
      IOPT_ITRINT(2) = ITRINT(2)
      IOPT_INTDEF    = INTDEF
      IOPT_INTFLG    = INTFLG
      IOPT_INTBUF    = INTBUF
      MOPT_MXMICRO   = MXMICRO
      MOPT_MXMACRO   = MXMACRO
      MOPT_MAXBCK    = MAXBCK
      IOPT_ISTATE    = 1
C
C     Set common block DCROPT
C
      OPT_CNVINT(1)  = CNVINT(1)
      OPT_CNVINT(2)  = CNVINT(2)
      IF (EVCCNV) THEN
         OPT_THRGRD     = SCFCNV(1)
      ELSE
         OPT_THRGRD     = SQRT(SCFCNV(1))
      END IF
      THRPCI = D0
C
C     Set common block DCLOPT
C
      OPT_SKIPEE     = DHF_SKIPEE
      OPT_SKIPEP     = DHF_SKIPEP
      OPT_NOCI       = .FALSE.
      OPT_UCIBOS     = .FALSE.
      OPT_MKCNST     = .FALSE.
C
C     Set common block DCCOPT
C
      OPT_CIPROGRAM  = 'noone HF'
      DO I = 1, NFSYM
         OPT_FRZSTR(I) = ' '
      END DO
C
C
      CALL QEXIT('TDHFORB')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tkrmcorb */
      SUBROUTINE TKRMCORB
C***********************************************************************
C
C     Transfer information from common block DCOMC to 
C     common block DCBORB and set common block DCxOPT
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbkrmc.h"
#include "dcbopt.h"
C
      CALL QENTER('TKRMCORB')
C
C
      NGAS   = NKRMCGAS
      MK2REF = NKRMC_MK2REF
      MK2DEL = NKRMC_MK2DEL
      MINMK2 = NKRMC_MINMK2
      MAXMK2 = NKRMC_MAXMK2
      CALL IZERO(NGSHT,MXMCGAS)
      DO I = 1,NFSYM
         NISH(I) = NKRMCISH(I)
         NASH(I) = NKRMCASH(I)
         DO J = 1, NGAS
            NGSH(I,J) = NKRMCGSH(I,J)
            NGSHT(J) = NGSHT(J) + NGSH(I,J)
         END DO
         NOCC(I) = NISH(I) + NASH(I)
         NSSH(I) = NKRMCSSH(I)
c         NESH(I) = NKRMCESH(I)
c         NPSH(I) = NKRMCPSH(I)
         NFRO(I) = NKRMCFRO(I)
      END DO
      NAELEC = NKRMCAELEC
C     min and max electrons in each gas space:
      DO J = 1, NGAS
         NGSP(1,J) = NKRMCGSP(1,J)
         NGSP(2,J) = NKRMCGSP(2,J)
      END DO
C
C
C     Set common block DCIOPT
C
      IPROPT = IPRKRMC
      IOPT_ITRINT(1) = IKRMC_ITRINT(1)
      IOPT_ITRINT(2) = IKRMC_ITRINT(2)
      IOPT_INTDEF    = IKRMC_INTDEF
      IOPT_INTFLG    = IKRMC_INTFLG
      IOPT_INTBUF    = IKRMC_INTBUF
      MOPT_MXMICRO   = MKRMC_MXMICRO
      MOPT_MXMACRO   = MKRMC_MXMACRO
      MOPT_MAXBCK    = MKRMC_MAXBCK
      IOPT_ISTATE    = IKRMC_ISTATE
      IOPT_SYMMETRY  = IKRMC_SYMMETRY
C
C     Set common block DCROPT
C
      OPT_CNVINT(1)  = DKRMC_CNVINT(1)
      OPT_CNVINT(2)  = DKRMC_CNVINT(2)
      OPT_THRGRD     = DKRMC_THRGRD
      THRPCI         = DKRMC_THRPCI
C
C     Set common block DCLOPT
C
      OPT_SKIPEE     = KRMC_SKIPEE
      OPT_SKIPEP     = KRMC_SKIPEP
      OPT_NOCI       = KRMC_NOCI
      OPT_UCIBOS     = KRMC_UCIBOS
      OPT_MKCNST     = KRMC_MKCNST
C
      DO I = 1, NOPTFLAGS
         FLAG(I) = MCFLAG(I)
      END DO
C
C     Set common block DCCOPT
C
      OPT_CIPROGRAM  = KRMC_CIPROGRAM
      DO I = 1, NFSYM
         OPT_FRZSTR(I) = KRMC_FRZSTR(I)
      END DO
C
C
      CALL QEXIT('TKRMCORB')
      RETURN
C
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnod */
      SUBROUTINE DIRNOD(IPARCAL,IMYTID,IMPARID,INUMNOD,WORK,MWORK)
C***********************************************************************
C
C     Menu file for slave nodes. Check whether the program needs a slave
C     and determine the type of work to be done.
C
C     Written by L.Visscher July 1997
C     Last revision: July 29 1996
C
C***********************************************************************
#include "implicit.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "priunit.h"
C
      INCLUDE 'mpif.h'
      INTEGER ISTAT(MPI_STATUS_SIZE)
      LOGICAL NEWTASK
C
#include "infpar.h"
C
      DIMENSION WORK(MWORK)
C
      CALL QENTER('DIRNOD')
      LWORK = MWORK
#include "memint.h"
C
C
C     Transfer variables to common blocks
C
      MPARID = IMPARID
C We count the number of slave nodes, subtract the master
      NUMNOD = INUMNOD - 1
      MYTID  = IMYTID
      PARCAL = (IPARCAL .NE. 0)
C
C     Set basis set directory for this slave
C     (it may be needed for reading the .mol file)
C
      CALL GETBASDIR()
C
 100  CONTINUE
C
      CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPARID,27,
     &              MPI_COMM_WORLD,ISTAT,IERR)
C
C     The code below may be useful for some architectures but 
C     fixes a MPI implementation problem and has as side effect
C     a slowdown om massively parallel calculations on small
C     systems with many tasks.  I therefore disabled it.
C     LV 6-8-2004
C
C     Check whether there is a new task waiting : if not
C     the slave may sleep for one second. This is useful on
C     shared memory architectures where MPI_RECV is implemented
C     using an idle loop that consumes unnecessary CPU when
C     waiting in a serial part of the code.
C
C     CALL MPI_IPROBE(MPARID,27,MPI_COMM_WORLD,
C    &                NEWTASK,ISTAT,IERR)
C     IF (NEWTASK) THEN
C        CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPARID,27,
C    &                 MPI_COMM_WORLD,ISTAT,IERR)
C     ELSE
C        CALL SLEEP(1)
C        GO TO 100
C     END IF
C
C     Start timing.
C
      CALL GETTIM(CPU1,WALL1)
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
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
      ELSE IF (ITASK .EQ. 4) THEN
C
C        Enter the DFT module
C
         CALL DFTNOD(WORK,LWORK)
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
C
      GOTO 100
C
 999  CALL QEXIT('DIRNOD')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnod2 */
      SUBROUTINE DIRNOD2(WORK,MWORK)
C***********************************************************************
C
C     Another menu file for slave nodes. This can be used to change the type
C     of work during a calculation (i.e. : if one parallel module calls the
C     other like in RELCCSD when doing the AO Lagrangian).
C
C***********************************************************************
#include "implicit.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "priunit.h"
C
      INCLUDE 'mpif.h'
      INTEGER ISTAT(MPI_STATUS_SIZE)
      LOGICAL NEWTASK
C
#include "infpar.h"
C
      DIMENSION WORK(MWORK)
C
      CALL QENTER('DIRNOD2')
      LWORK = MWORK
#include "memint.h"
C
 100  CONTINUE
C
      CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPARID,27,
     &              MPI_COMM_WORLD,ISTAT,IERR)
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
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
      ELSE IF (ITASK .EQ. 4) THEN
C
C        Enter the DFT module
C
         CALL DFTNOD(WORK,LWORK)
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
C
      GOTO 100
C
 999  CALL QEXIT('DIRNOD2')
      END
#endif /* VAR_MPI */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirctl */
      SUBROUTINE GETBASDIR()
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
#include "infpar.h"
#include "dcbgen.h"
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
      IF (MYTID .EQ. MPARID) WRITE(LUPRI,'(2A)')
     $     ' Basis set dir.    : ',BASDIR(1:LENBAS)
C 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
C
C
#if defined (VAR_MPI)
#include "infpar.h"
      INCLUDE 'mpif.h'
C
      IF (ITASK .EQ. -1) THEN
C        go back from submenu to main menu via MPI_BCAST
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      ELSE
C        slave should be in main menu, ready to receive task
         DO NODE = 1,NUMNOD
            CALL MPI_SEND(ITASK,1,MPI_INTEGER,NODE,27,
     &                 MPI_COMM_WORLD,IERR)
         END DO
      END IF

#else /* not VAR_MPI */

      CALL QUIT('DIRAC_PARCTL called '//
     &   'but this is a sequential Dirac version !')

#endif /* VAR_MPI */
      RETURN
      END
