C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampsi */
      SUBROUTINE PAMPSI(WORK,LWORK)
C*****************************************************************************
C
C     Wave function module
C
C     Written by T.Saue April 1996
C     Last revision: J. Thyssen - Oct 22 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbdhf.h>
C
      CALL QENTER('PAMPSI')
C
C
C*****************************************************************************
C***** D I R A C - H A R T R E E - F O C K ***********************************
C***** o r   D I R A C - K O H N - S H A M ***********************************
C*****************************************************************************
C
      IF(DOSCF .AND. DOPSI) THEN
C
C        Transfer DHF info to common block DCBORB
C        and recalculate DCBORB.
C
         CALL TDHFORB
         CALL SETDC2(5)
c        CALL SETDC2(0)
c        call setdc2(iprint)
C
C
C        *************************************************
C        *** Call PSIDHF for Roothaan/DIIS iterations  ***
C        *************************************************
C
         IF (MAXITR.GT.0) THEN
            CALL PSIDHF(WORK,LWORK)
         END IF
C
C        *************************************************
C        *** Call PSIOPT for second order optimization ***
C        *************************************************
C
         IF ( .NOT. NOQCDHF .AND. .NOT. DODFT .AND.
     &       (.NOT. WFCONV .AND. MAXITR .GT.0 )) THEN
            CALL PSIOPT(WORK,LWORK)
         END IF
      END IF
C
C*****************************************************************************
C***** R E S O L V E    O P E N - S H E L L  S T A T E S *********************
C*****************************************************************************
C
      IF(DORES.AND.DOPSI) CALL RESOLV(WORK,LWORK)
C
C*****************************************************************************
C***** M O E L L E R - P L E S S E T  ****************************************
C*****************************************************************************
C
      IF(DOMP2.AND.DOPSI) THEN
C
C        *************************************************
C        *** Call PSIMP2 for MP2 Energy/wave function  ***
C        *************************************************
C
         IF (DODHF .AND. .NOT. WFCONV) THEN
             CALL QUIT('PAMPSI: calling MP2 without converged DHF'//
     &                 ' wave function.')
         END IF
         CALL PSIMP2(WORK,LWORK)
      END IF
C
C*****************************************************************************
C***** M O D I F I E D   V I R T U A L   O R B I T A L S *********************
C*****************************************************************************
C
      IF(DOMVO.AND.DOPSI) CALL MVOGEN(WORK,LWORK)
C
C*****************************************************************************
C***** I N D E X  T R A N S F O R M A T I O N ********************************
C*****************************************************************************
C
      IF(DOTRA) CALL PAMTRA(WORK,LWORK)
C
C*****************************************************************************
C***** R E L C C S D  M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCCM.AND.DOPSI) CALL PAMCCM(WORK,LWORK)
C
C*****************************************************************************
C***** D I R R C I   M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCIM.AND.DOPSI) CALL PAMCIM(WORK,LWORK)
C
C*****************************************************************************
C***** L U C I T A  M O D U L E **********************************************
C*****************************************************************************
C
      IF(DOLUCT.AND.DOPSI) CALL PAMLUCITA(WORK,LWORK)
Chj TODO: move PAMLUCITA after MCSCF ????
C
C
#if defined (UNDEF)
C
C*****************************************************************************
C***** C I *******************************************************************
C*****************************************************************************
C
      IF(DOCI) THEN
C
C        *************************************************
C        *** Call PSICI for CI wave function *************
C        *************************************************
C
         IF (DODHF .AND. .NOT. WFCONV) THEN
             CALL QUIT('PAMPSI: calling CI without converged DHF'//
     &                 ' wavefunction.')
         END IF
         CALL PSICI(WORK,LWORK)
      END IF
#endif
C
C*****************************************************************************
C***** K R - M C S C F *******************************************************
C*****************************************************************************
C
      IF(DOKRMC.AND.DOPSI) THEN
C
C        Transfer DHF info to common block DCBORB
C        and recalculate DCBORB.
C
         CALL TKRMCORB
         CALL SETDC2(0)
C
C
C        *************************************************
C        *** Call PSIOPT for KR-MCSCF wave function ******
C        *************************************************
C
c         IF (DODHF .AND. .NOT. WFCONV) THEN
c             CALL QUIT('PAMPSI: calling KF-MCSCF without converged '//
c     &                 ' DHF wavefunction.')
c         END IF
         CALL PSIOPT(WORK,LWORK)
      END IF
C
      CALL QEXIT('PAMPSI')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tdhforb */
      SUBROUTINE TDHFORB
C*****************************************************************************
C
C     Transfer information from common block DCODHF to common block DCBORB
C     and set common block DCxOPT
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision: 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <consts.h>
C
#include <dgroup.h>
#include <dcbdhf.h>
#include <dcborb.h>
#include <dcbopt.h>
C
      CALL QENTER('TDHFORB')
C
C
      NASHT = 0
      DO I = 1,NFSYM
         NISH(I) = NISH_DHF(I)
         NASH(I) = NASH_DHF(I)
         NASHT   = NASHT + NASH(I)
         NOCC(I) = NOCC_DHF(I)
         NSSH(I) = NSSH_DHF(I)
c         NESH(I) = NESH_DHF(I)
c         NPSH(I) = NPSH_DHF(I)
         NFRO(I) = NFRO_DHF(I)
         NELEC(I) = NELEC_DHF(I)
      END DO
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
C     initialize min and max accumlated no. of electrons in each gas space
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
         OPT_THRGRD     = SCFCNV
      ELSE
         OPT_THRGRD     = SQRT(SCFCNV)
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
C
      CALL QEXIT('TDHFORB')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tkrmcorb */
      SUBROUTINE TKRMCORB
C*****************************************************************************
C
C     Transfer information from common block DCOMC to common block DCBORB
C     and set common block DCxOPT
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision: 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbkrmc.h>
#include <dcbopt.h>
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
         NOCC(I) = NKRMCOCC(I)
         NSSH(I) = NKRMCSSH(I)
c         NESH(I) = NKRMCESH(I)
c         NPSH(I) = NKRMCPSH(I)
         NFRO(I) = NKRMCFRO(I)
         NELEC(I) = NKRMCELEC(I)
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
