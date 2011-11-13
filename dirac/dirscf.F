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
C VERSION : $Revision: 1.112.2.3 $
C DATE    : $Date: 2004/09/01 08:08:53 $
C FILE    : dirscf.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psidhf */
      SUBROUTINE PSIDHF(WFCONV,WORK,LWORK)
C***********************************************************************
C
C     Solve DHF - equations; print results
C
C     Written by T.Saue November 1994
C     Last revision : Sep 11 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbpsi.h"
#include "blocks.h"
C
      LOGICAL   WFCONV
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PSIDHF')
#include "memint.h"
C
      WFCONV = .FALSE.
      IF(DODFT) THEN
        CALL TITLER('Kohn - Sham calculation','*',125)
      ELSE
        CALL TITLER('Dirac-Hartree-Fock calculation','*',125)
      ENDIF
      CALL FLSHFO(LUPRI)
C
C     Initialize
C     ==========
C
      CALL SETDHF(IPRSCF)
C
C     Memory allocation
C
C     KCMO is N2BBASXQ because it is used for Fock matrix in PREDHF.
C
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,      WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBRP,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFOCK,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
C
C     Solve DHF-equations
C     ===================
C
      CALL DHFSCF(WORK(KCMO),WORK(KEIG),WORK(KIBRP),WORK(KFOCK),
     &            WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('PSIDHF.DHFSCF',WORK,KWORK,KFOCK,KFREE,LFREE)
C
C     Output section
C
      CALL DHFOUT(WORK(KCMO),WORK(KEIG),WORK(KIBRP),WORK(KFREE),
     &            LFREE,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('PSIDHF.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Check if DHF/DFT is converged (to hard or soft threshold):
C
      WFCONV = DHFCONV(1) .OR. DHFCONV(2)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('PSIDHF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck predhf */
      SUBROUTINE PREDHF(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,NPOS,WORK,LWORK)
C*****************************************************************************
C
C     Obtain orthonormal trial vectors for Dirac-Fock
C     either by:
C       i) bare nucleus approximation
C      ii) starting from trial vectors
C     iii) starting from two-electron Fock matrix
C     Also serves as solver module for the Dirac-equation in the
C     finite basis approximation (ONESYS)
C
C     Numerous possible schemes for restart has been implemtented.
C     They are based on the program looking for necessary files.
C     Three modes should be discerned:
C
C     1) One-electron system (ONESYS)
C     2) Many-electron system
C         a) DFCYCL present:
C            The program will restart on as much as possible.
C            (coefficients, convergence acceleration,
C             differential density etc.)
C         b) DFCYCL abscent:
C            Restart is only possible on coefficients/Fock matrix.
C
C     Written by T.Saue November 1994
C     Last revision : Jan 27 1997 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0,D0 = 0.0D0, D4 = 4.00D00,
     &          D2 = 2.0D0, DP5=0.5D0)
C
#include "mxcent.h"
#include "nuclei.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbdhf.h"
#if defined (MOD_DFT)
#include "dftcom.h"
#endif
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbreo.h"
      LOGICAL TOBE,DOSEL,DSCFON,BARNUC,TEST,DIISAO_SAVE
      CHARACTER CPUTID*12,DAYTID*24,SECTID*12
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*),
     &          BMAT(MXDIIS,MXDIIS),DMAT(*),
     &          NPOS(*),WORK(LWORK),
     &          DVOVLP(MXOPEN)
C
      CALL QENTER('PREDHF')
C
#include "memint.h"
C
      TEST = DIISAO.OR.NASHT.GT.0.OR.NAELEC.GT.0
C
C     ****************************************************
C     ***  I N I T I A L I Z E  a n d   R E S T A R T ****
C     ****************************************************
C
      CALL DHFTST(BARNUC,DSCFON,ICHKCNV,BMAT,IPRINT)
      IF(.NOT.ONESYS) THEN
C
C       Initiate SCF cycle
C
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +  '########## START ITERATION NO.',NITER,' ##########',
     +  DAYTID
        CALL FLSHFO(LUPRI)
C
C       Open files necessary for Dirac-Fock calculation
C
        OPEN (LUEVEC, FILE='DFEVEC', FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*N2EVECTQ,STATUS='UNKNOWN')
        CALL OPNFIL(LUFCK2,'DFFCK2','UNKNOWN','PREDHF')
        CALL OPNFIL(LUDENS,'DFDENS','UNKNOWN','PREDHF')
        CALL OPNFIL(LUCMOS,'DFCMOS','UNKNOWN','PREDHF')
      ENDIF
      OPEN (LUFOCK, FILE='DFFOCK', FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*N2TMOTQ,STATUS='UNKNOWN')
      CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','PREDHF')
      CALL OPNFIL(LUFCK1,'DFFCK1','OLD','PREDHF')
C
C     *****************************************************************
C     *****  G E N E R A T E   F U L L    F O C K    M A T R I X  *****
C     *****************************************************************
C
      CALL FLSHFO(LUPRI)
      CALL GETTIM(CTIDSTR,WTIDSTR)
C
C     Get 2-electron Fock matrix (save in CMO)
C     ========================================
C
C     1) - using trial vectors
C     -------------------------
C
      DOSEL = .FALSE.
      TOBE  = .TRUE.
      IF(TRIVEC) THEN
C
C       We calculate a regular Fock matrix so we can check convergence ...
C
        IF (ICHKCNV.EQ.0) ICHKCNV = -1
        IF(IPRSCF.GE.1) THEN
          CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,15,
     &              WORK,KFREE,LFREE)
        ELSE
          CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14,
     &              WORK,KFREE,LFREE)
        ENDIF
C
C       Reorder orbitals if requested
C
        IF (LMOORD) THEN
          CALL REORD(CMO,EIG,IBEIG,IMOORD,IREORD,MXREORD,
     &         WORK(KFREE),LFREE)
C
C         Make sure that reordered coefficients are written back
C         to file.
C
          CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
        END IF
C
C       Determine occupation if this was not given in input
C
        IF (INIOCC) THEN
           CALL DETOCC (EIG,WORK,KFREE,LFREE)
           INIOCC = .FALSE.
        ENDIF
C
C       Prepare for overlap selection
C
        IF(OVLSEL) THEN
          DOSEL = .TRUE.
          INQUIRE(FILE='DFSMOS',EXIST=TOBE)
          IF(.NOT.TOBE .OR. LMOORD) THEN
C            ... new order for overlap selection if reordered (LMOORD)
             CALL OPNFIL(LUTMAT,'DFTMAT','OLD','PREDHF')
             CALL READT(LUTMAT,N2TMT,DMAT)
             CLOSE(LUTMAT,STATUS='KEEP')
             CALL SMOGEN(CMO,DMAT,WORK(KFREE),LFREE)
          ENDIF
        ENDIF
C
C        Generate density matrix/matrices
C        Be aware: All density matrices are normalized to 1,
C                  occupation is in the DF(*) factor
C
         CALL GENDEN(DMAT,CMO,1,IPRSCF)
         DO IOPEN = 1,NOPEN
            CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPRSCF)
         END DO
C
C       Differential density matrix approach
C
        TDF2 = SECOND()
        IF(DSCFON) THEN
          WRITE(LUPRI,'(I5,A)') NITER,
     &       ' *** Differential density matrix ***'
          CALL READNS(LUDENS,FOCK)
          CALL WRIDNS(LUDENS,DMAT)
          IF (NISHT .GT. 0) THEN
            IF(FIXDIF) THEN
              DCOVLP = D1
            ELSE
              DCOVLP = DDOT(N2BBASXQ,FOCK,1,FOCK,1)
              DCOVLP = DDOT(N2BBASXQ,DMAT,1,FOCK,1) / DCOVLP
            ENDIF
            CALL DAXPY(N2BBASXQ,(-DCOVLP),FOCK,1,DMAT,1)
          END IF
          IF (NASHT .GT. 0) THEN
             DO IOPEN = 1,NOPEN
                IF (FIXDIF) THEN
                   DVOVLP(IOPEN) = D1
                ELSE
                   DVOVLP(IOPEN) =
     &                DDOT(N2BBASXQ,FOCK(1+N2BBASXQ*IOPEN),1,
     &                              FOCK(1+N2BBASXQ*IOPEN),1)
                   DVOVLP(IOPEN) =
     &                DDOT(N2BBASXQ,DMAT(1+N2BBASXQ*IOPEN),1,
     &                FOCK(1+N2BBASXQ*IOPEN),1) / DVOVLP(IOPEN)
                END IF
                WRITE(LUPRI,'(I5,A,I2,A,F6.4)') NITER,
     &             ' *** Differential density matrix. DVOVLP(',
     &             IOPEN,') =',DVOVLP(IOPEN)
                CALL DAXPY(N2BBASXQ,(-DVOVLP(IOPEN)),
     &               FOCK(1+N2BBASXQ*IOPEN),1,DMAT(1+N2BBASXQ*IOPEN),1)
             END DO
          END IF
C
C           Note that all Fock matrices are multiplied with 2
C              in TWOFCK.
C
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
C           read old Fock matrix in DMAT to finish new Fock matrix/hj
          CALL REAFCK(LUFCK2,DMAT,.FALSE.,NFMAT)
          IF (NISHT.GT.0)
     &         CALL DAXPY(N2BBASXQ,DCOVLP,DMAT,1,FOCK,1)
          IF (NASHT.GT.0) THEN
            DO IOPEN = 1,NOPEN
              CALL DAXPY(N2BBASXQ,DVOVLP(IOPEN),
     &                   DMAT(1+N2BBASXQ*IOPEN),1,
     &                   FOCK(1+N2BBASXQ*IOPEN),1)
            END DO
          END IF
C         restore density matrix in DMAT/tsaue
          CALL READNS(LUDENS,DMAT)
        ELSE
          CALL WRIDNS(LUDENS,DMAT)
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
        ENDIF
C
C
        IF(DOTSC) CALL TSCORR(FOCK,WORK(KFREE),LFREE,IPRSCF)
C
C
C       Write current two-electron Fock matrix to file
C       ==============================================
C
        CALL WRIFCK(LUFCK2,FOCK,NFMAT)
C
        TDXC = SECOND()
        TDF2 = TDXC - TDF2
#if defined (MOD_DFT)
C       Construct Kohn-Sham contribution
C       ================================
        IF(DFTADD) THEN
C          tsaue: turn off test on numerical integration,
C                 When starting from coefficients from a different
C                 geometry; they are not orthonormal at the new
C                 geometry and thus the correct number of electrons
C                 is generally not integrated
           TMP = DFTELS
           DFTELS = D1
C          NFMAT=1,ISYMOP=1,IHRMOP=1,DOERG=T,DOGRD=F,DOLND=F,DOATR=F
           CALL DFTEXC(DMAT,DUM,FOCK,1,1,1,.TRUE.,.FALSE.,
     &                 .FALSE.,.FALSE.,WORK(KFREE),LFREE,IPRSCF)
           DFTELS = TMP
           TDXC = SECOND() - TDXC
        ENDIF
#endif
C
C       Solvent contribution
C
        IF(SOLVEN) THEN
           CALL SOLFCK(FOCK,DMAT,1,ESOLVE,ESOLVN,
     &                 WORK(KFREE),LFREE,IPRSOL)
        ENDIF
      ENDIF
C
C     2) - using old 2-electron Fock matrix
C     --------------------------------------
C
      IF(TRIFCK) THEN
         CALL REAFCK(LUFCK2,FOCK,.TRUE.,NFMAT)
      END IF
C
C     Get one-electron Fock matrix (save in FOCK)
C     ===========================================
C
      IF(BARNUC) THEN
        CALL REAFCK(LUFCK1,FOCK,.FALSE.,1)        
        IF(BNCRON .AND. .NOT. ONESYS) THEN
           CALL BNCORR(FOCK,WORK(KFREE),LFREE,IPRSCF)
        ENDIF
      ELSE
        CALL REAFCK(LUFCK1,CMO,.FALSE.,1)
C        
C     Calculate correction according to the One-center app. models 
C     ============================================================
C
        IF (TRIVEC) THEN
           IF (ONECAP.AND. .NOT.SMLV1C) THEN
              CALL CATCORR(CMO,FOCK,DMAT,WORK(KFREE),LFREE,IPRSCF)
           END IF 
           CALL ERGCAL(CMO,FOCK,DMAT,WORK(KFREE),LFREE)
        END IF
        CALL DAXPY(N2BBASXQ,D1,CMO,1,FOCK,1)
C
C       Add active two-electron Fock matrices
C
        DO IOPEN = 1,NOPEN
C
C         The factor in only DF and not D2*DF because of the
C         multiplication with D2 in TWOFCK.
C
          CALL DAXPY(N2BBASXQ,DF(IOPEN),
     &               FOCK(1+N2BBASXQ*IOPEN),1,FOCK,1)
        END DO
C
        IF(IPRSCF.GE.5) THEN
          IF (NISHT .GT. 0) THEN
            CALL HEADER('Total FD matrix',-1)
            CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          END IF
          DO IOPEN = 1,NOPEN
            CALL HEADER('Total FV matrix',-1)
            WRITE(LUPRI,'(A,I3)') ' Matrix number',IOPEN
            CALL PRQMAT(FOCK(1+N2BBASXQ*IOPEN),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          END DO
         END IF
      ENDIF
C
C     *************************************************
C     *****   SOLVE  DHF  EIGENVALUE  PROBLEM    ******
C     *************************************************
C
      DIISAO_SAVE = DIISAO
      DIISAO = .TRUE.
C     ... only DIISAO works in PREDHF because no MOs on LUCMOS yet ...
      CALL DFSOLV(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,ICHKCNV,
     &            DOSEL,WORK(KFREE),LFREE)
      DIISAO = DIISAO_SAVE
C
      IF(BARNUC.OR.TRIFCK) CALL ONEERG(EIG,IBEIG)
C
      CALL GETTIM(CTIDEND,WTIDEND)
      CTIDTOT = CTIDEND - CTIDSTR
      WTIDTOT = WTIDEND - WTIDSTR
      IF(IPRSCF.GE.1) THEN
        CPUTID = SECTID(TDDG)
        IF(NZ.EQ.1) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Real diagonalization      (CPU):',CPUTID
        IF(NZ.EQ.2) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Complex diagonalization   (CPU):',CPUTID
        IF(NZ.EQ.4) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Quaternion diagonalization(CPU):',CPUTID
      ENDIF
      IF(ONESYS) THEN
        IF(PHCOEF) THEN
          DO I = 1,NFSYM
            CALL PHATRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ)
          ENDDO
          CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG)
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
        CLOSE(LUFCK1,STATUS='KEEP')
        CLOSE(LUFOCK,STATUS='KEEP')
        CPUTID = SECTID(WTIDTOT)
        WRITE(LUPRI,'(A,A12)')
     &    '* Total time for one-electron system:',CPUTID
      ELSE
        ITRSCF(INTFLG) = ITRSCF(INTFLG) + 1
        SCFTID(INTFLG) = SCFTID(INTFLG) + WTIDTOT
        IF(IPRSCF.GE.1) THEN
          CPUTID = SECTID(TDF2)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* Two-electron Fock matrix  (CPU):',CPUTID
          IF(DFTADD) THEN
          CPUTID = SECTID(TDXC)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* DFT xc matrix             (CPU):',CPUTID
          ENDIF
        ENDIF
        CPUTID = SECTID(CTIDTOT)
        CALL GTINFO(DAYTID)
        WRITE(LUCYCL,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID
        CALL FLSHFO(LUCYCL)
        WRITE(LUPRI, 1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID
      ENDIF
C
      IF (TRIFCK) THEN
C
C       Get density matrices (normed to one)
C       and write them to DFDENS.
C
        CALL GENDEN(DMAT,CMO,1,IPRSCF)
        DO IOPEN = 1,NOPEN
          CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPRSCF)
        END DO
        CALL WRIDNS(LUDENS,DMAT)
      END IF
C
      CALL FLSHFO(LUPRI)
      CALL MEMREL('PREDHF',WORK,KWORK,KWORK,KFREE,LFREE)
C
 666  CONTINUE
      CALL QEXIT('PREDHF')
      RETURN
 1010 FORMAT('It. ',I3,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A9,3X,A10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dfsolv */
      SUBROUTINE DFSOLV(FMO,EIG,IBEIG,FAO,DMAT,BMAT,ICHKCNV,DOSEL,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Solve Dirac-Fock eigenvalue equation
C     On input:
C       FMO - MO-coefficients from previous iteration
C       FAO - Fock-matrix in AO-basis
C     On output:
C       FMO - MO-coefficients from current iteration
C       EIG - eigenvalues
C
C     Input/output:
C
C     Written by T.Saue Mar 09 1995
C     Last revision: Mar 09 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00)
      LOGICAL DOSEL,TEST
      DIMENSION FMO(*),EIG(*),IBEIG(*),FAO(*),DMAT(*),BMAT(*),
     &          WORK(LWORK)
C
#include "memint.h"
C
      TEST = NASHT.GT.0.OR.DIISAO.OR.NAELEC.GT.0
      IF (TEST) THEN
        CALL MEMGET('REAL',KDTEMP,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KDEVEC,N2BBASXQ      ,WORK,KFREE,LFREE)
      ELSE
        KDTEMP = KFREE
        KDEVEC = KFREE
      ENDIF
C
C     Transform AO Fock matrix to reference MO basis
C     ============================================== 
C
      CALL MEMGET('REAL',KTMAT,N2TMT,WORK,KFREE,LFREE)
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','PREDHF')
      CALL READT(LUTMAT,N2TMT,WORK(KTMAT))
      IF (NASHT.GT.0.AND.TRIVEC) THEN
        CALL MEMGET('REAL',KCMOSAV,NCMOTQ,WORK,KFREE,LFREE)
        CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMOSAV),EIG,IBEIG,TOTERG,6,
     &              WORK,KFREE,LFREE)
      ELSE
        KCMOSAV = KFREE
      END IF
      IF (TEST) THEN
C        read LUTMAT rec#2 with SAOMO
C        NZT is the NZ value for the transformation matrix
        ISIZE = 0
        DO I = 1,NFSYM
           ISIZE = ISIZE + NZT*NTMO(I)*NFBAS(I,0)
        END DO
        CALL MEMGET('REAL',KSAOMO,ISIZE,WORK,KFREE,LFREE)
        CALL READT(LUTMAT,ISIZE,WORK(KSAOMO))
        CALL MKMOFK(FMO,FAO,DMAT,WORK(KTMAT),WORK(KCMOSAV),
     &              WORK(KSAOMO),WORK(KDTEMP),WORK(KDEVEC),
     &              WORK(KFREE),LFREE)
      ELSE
        KSAOMO = KFREE
        CALL MKMOFK(FMO,FAO,DUM,WORK(KTMAT),DUM,
     &              DUM,DUM,DUM,WORK(KFREE),LFREE)
      END IF
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Convergence control/acceleration
C     ================================
C
      IF(TEST) THEN
        CALL DHFCNV(FMO,FAO,DMAT,WORK(KDEVEC),WORK(KDTEMP),WORK(KTMAT),
     &              BMAT,ICHKCNV,WORK(KSAOMO),WORK(KFREE),LFREE)
      ELSE
        CALL DHFCNV(FMO,DUM,DUM,FAO,DMAT,DUM,BMAT,ICHKCNV,
     &              DUM,WORK(KFREE),LFREE)
      ENDIF
      IF (DHFEXIT) GO TO 9999
C
C     Diagonalize Fock matrix
C     =======================
C
      TDDG = SECOND()
      CALL DFDIAG (FMO,EIG,IBEIG,FAO,WORK,KFREE,LFREE)
      TDDG = SECOND() - TDDG
C
C
C     Automatic occupation
C     ====================
C
      IF(AUTOCC.OR.INIOCC) CALL DETOCC(EIG,WORK,KFREE,LFREE)
C
C     Overlap selection
C     =================
C
      IF(DOSEL) CALL SMOSEL(EIG,IBEIG,FAO,WORK(KFREE),LFREE)
C
C
C     Backtransform all spinors (temporary storage in FMO)
C     ====================================================
C
      IF(NFROT.NE.0) CALL OPNFIL(LUTMAT,'DFFROZ','OLD','DFSOLV')
      DO 40 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 40
        IF(NFRO(I).EQ.0) THEN
        CALL BCKTR1(FMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              FAO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &              NTMO(I),NZ,NTMO(I),1,NFBAS(I,0),
     &              WORK(KTMAT+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,
     &              IPRSCF)
        ELSE
          CALL MEMGET('INTE',KIND,NFRO(I),WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KBUF,NFBAS(I,0)*NZ,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KEBF,NTMO(I),WORK,KFREE,LFREE)
          CALL BCKTR2(I,LUTMAT,FMO(ICMOQ(I)+1),FAO(I2TMOT(I)+1),
     &                EIG(IORB(I)+1),WORK(KTMAT+I2TMT(I)),
     &                WORK(KIND),WORK(KBUF),WORK(KEBF),IPRSCF)
          CALL MEMREL('DFSOLV.froz',WORK,KIND,KIND,KFREE,LFREE)
        ENDIF
 40   CONTINUE
      IF(NFROT.GT.0) CLOSE(LUTMAT,STATUS='KEEP')
C
C     Print eigenvalues 
C     (after backtransformation in case of frozen orbitals)
C
      IF (IPRSCF.GE.3)
     &   CALL PREIGN(EIG,IBEIG,(IPRSCF .GE. 4),WORK(KFREE),LFREE)
C
C     Write coefficients to file
C     ==========================
C
C     Write untransformed coefficients
      REWIND LUCMOS
      CALL WRITT(LUCMOS,N2TMOTQ,FAO)
C     Write transformed coefficients
      CALL WRICMO(LUCOEF,FMO,EIG,IBEIG,DHFERG)
C
C
C     Dynamic update of overlap selection matrix
C
      IF (DYNSEL) THEN
        CALL SMOGEN(FMO,WORK(KTMAT),WORK(KFREE),LFREE)
      ENDIF
C
 9999 CALL MEMREL('DFSOLV',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfscf */
      SUBROUTINE DHFSCF(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "cbihr1.h"
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*),
     &          WORK(LWORK)
C
C     Generate one-electron Fock matrix and write to file
C     ===================================================
C
      CALL ONEFCK(FOCK,IPRONE,WORK,LWORK)
      CALL OPNFIL(LUFCK1,'DFFCK1','UNKNOWN','ONEFC1')
      CALL WRIFCK(LUFCK1,FOCK,1)
      CLOSE(LUFCK1,STATUS='KEEP')
      IF(ONESYS) THEN
C
C     Solve one-electron system ...
C     =============================
C
        CALL PREDHF(CMO,EIG,IBEIG,FOCK,DUM,DUM,IDUM,WORK,LWORK)
        DHFCONV(1) = .TRUE.
        DHFCONV(2) = .TRUE.
      ELSE
C
C     ... or perform SCF-cycle
C     ========================
C
        CALL SCFCYC(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densty */
      SUBROUTINE DENSTY(DMAT,CMO,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),CMO(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),
     &              CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              NPSH(I)+1,NISH(I),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.10) THEN
        CALL TITLER('Output from DENSTY','*',103)
        DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck denst1 */
      SUBROUTINE DENST1(DMAT,LRD,LCD,CMO,LRC,LCC,JVEC,NVEC,NBAS)
C*****************************************************************************
C
C     Construct density matrix for fermion corep IFRP
C
C     Written by T.Saue Mar 3 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,D0=0.0D0)
#include "dgroup.h"
C
      DIMENSION DMAT(LRD,LCD,NZ),CMO(LRC,LCC,NZ)
      CALL QGEMM(NBAS,NBAS,NVEC,D1,
     &     'N','N',IPQTOQ(1,0),CMO(1,JVEC,1),LRC,LCC,NZ,
     &     'H','N',IPQTOQ(1,0),CMO(1,JVEC,1),LRC,LCC,NZ,
     &          D0,IPQTOQ(1,0),DMAT,LRD,LCD,NZ)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck scfcyc */
      SUBROUTINE SCFCYC(CMO,EIG,IBEIG,FOCK,WORK,LWORK)
C*****************************************************************************
C
C     Perform Dirac-Fock SCF-cycle
C
C     Written by T.Saue Mar 3 1995
C     Last revision: tsaue Mar 3 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*), WORK(LWORK)
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
      CALL QENTER('SCFCYC')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KDMAT,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT,MXDIIS*MXDIIS,WORK,KFREE,LFREE)
      NPOS = 0
      IF(PARCAL) THEN
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
      ENDIF
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
C     Initiate/restart SCF-cycle
C
      CALL SCFCY1(CMO,EIG,IBEIG,FOCK,WORK(KDMAT),WORK(KBMAT),
     &            WORK(KPOS),WORK(KFREE),LFREE)
      CALL MEMREL('SCFCYC.SCFCY1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('SCFCYC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck scfcy1 */
      SUBROUTINE SCFCY1(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,NPOS,WORK,LWORK)
C*****************************************************************************
C
C     Perform Dirac-Fock SCF-cycle
C
C     Written by T.Saue Mar 3 1995
C     Last revision: tsaue Mar 3 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DM1=-1.0D0,D1=1.0D0, DP5 = 0.50D00,
     &          D2 = 2.00D00)
C
#include "mxcent.h"
#include "nuclei.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbdhf.h"
#if defined (MOD_DFT)
#include "dftcom.h"
#endif
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbreo.h"
      CHARACTER CPUTID*12,SECTID*12,DAYTID*24
      DIMENSION FOCK(*),EIG(*),IBEIG(*),CMO(*),DMAT(*),
     &          BMAT(*),NPOS(*),WORK(LWORK)
      DIMENSION DVOVLP(MXOPEN)
      CALL QENTER('SCFCY1')
#include "memint.h"
C
C
C     ***********************************
C     ***** Perform SCF-iterations ******
C     ***********************************
C
C     Initiate/restart SCF cycle
C
      DHFEXIT = .FALSE.
      CALL PREDHF(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,NPOS,WORK,LWORK)
C     ... this is a restarted calculation we might already be converged
      IF(DHFEXIT) GOTO 30
C
C     Start cycle
C
      CALL FLSHFO(LUPRI)
      DO 10 ITER = 2,MAXITR
         NITER  = NITER + 1
         KITER  = KITER + 1
         CALL GTINFO(DAYTID)
         CALL GETTIM(CTIDSTR,WTIDSTR)
         WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +     '########## START ITERATION NO.',NITER,'##########',
     &     DAYTID
         CALL FLSHFO(LUPRI)
C
C        Check what integral types contribute
C        ====================================
C
         CALL INTCON(INTFLG,INTBUF,INTDEF,
     &               ERGVAL,CNVINT,NITER,ITRINT,DHF_INTTYP)
C
C        Construct density matrix
C        ========================
C
C        Get density matrices (normed to one)
C
         CALL GENDEN(DMAT,CMO,1,IPRSCF)
         DO IOPEN = 1,NOPEN
            CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPRSCF)
         END DO
C
C        CMO not needed any more, may be reused between here and DFSOLV/hj
C
C        Construct two-electron Fock matrix
C        (stored in FOCK)
C        ==================================
C
         TDF2 = SECOND()
         IF(DODSCF.AND.INTFLG.EQ.INTBUF) THEN
C           *** Differential density matrix code ***
C           read old density matrix in FOCK before saving new density matrix/hj
            CALL READNS(LUDENS,FOCK)
            CALL WRIDNS(LUDENS,DMAT)
            IF (NISHT .GT. 0) THEN
               IF(FIXDIF) THEN
                  DCOVLP = D1
               ELSE
                  DCOVLP = DDOT(N2BBASXQ,FOCK,1,FOCK,1)
                  DCOVLP = DDOT(N2BBASXQ,DMAT,1,FOCK,1) / DCOVLP
               END IF
               WRITE(LUPRI,'(I5,A,F7.4)') NITER,
     &            ' *** Differential density matrix. DCOVLP     =',
     &            DCOVLP
               IF (DCOVLP .LT. D0 .OR. DCOVLP .GT. D2) THEN
                  WRITE(LUPRI,'(A)') ' INFO: Useless value for DCOVLP,'
     &               //' value reset to: 1.0D0'
                  DCOVLP = D1
               END IF
               CALL DAXPY(N2BBASXQ,(-DCOVLP),FOCK,1,DMAT,1)
            END IF
            IF (NASHT.GT.0) THEN
               DO IOPEN = 1,NOPEN
                  IF (FIXDIF) THEN
                     DVOVLP(IOPEN) = D1
                  ELSE
                     DVOVLP(IOPEN) =
     &                  DDOT(N2BBASXQ,FOCK(1+N2BBASXQ*IOPEN),1,
     &                                FOCK(1+N2BBASXQ*IOPEN),1)
                     DVOVLP(IOPEN) =
     &                  DDOT(N2BBASXQ,DMAT(1+N2BBASXQ*IOPEN),1,
     &                  FOCK(1+N2BBASXQ*IOPEN),1) / DVOVLP(IOPEN)
                  END IF
                  WRITE(LUPRI,'(I5,A,I2,A,F7.4)') NITER,
     &               ' *** Differential density matrix. DVOVLP(',
     &               IOPEN,') =',DVOVLP(IOPEN)
                  IF (DVOVLP(IOPEN) .LT. D0
     &                  .OR. DVOVLP(IOPEN) .GT. D2) THEN
                     WRITE(LUPRI,'(A,I2,A)')
     &                  ' INFO: Useless value for DVOVLP(',IOPEN,'),'
     &                  //' value reset to: 1.0D0'
                     DVOVLP(IOPEN) = D1
                  END IF
                  CALL DAXPY(N2BBASXQ,(-DVOVLP(IOPEN)),
     &               FOCK(1+N2BBASXQ*IOPEN),1,DMAT(1+N2BBASXQ*IOPEN),1)
               END DO
            END IF
C
C           Note that all Fock matrices are multiplied with 2
C              in TWOFCK.
C
            CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                  NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
C           read old Fock matrix in DMAT to finish new Fock matrix/hj
            CALL REAFCK(LUFCK2,DMAT,.FALSE.,NFMAT)
C
            IF (NISHT .GT. 0)
     &         CALL DAXPY(N2BBASXQ,DCOVLP,DMAT,1,FOCK,1)
            IF (NASHT .GT. 0) THEN
               DO IOPEN = 1,NOPEN
                  CALL DAXPY(N2BBASXQ,DVOVLP(IOPEN),
     &               DMAT(1+N2BBASXQ*IOPEN),1,
     &               FOCK(1+N2BBASXQ*IOPEN),1)
               END DO
            END IF
C           restore density matrix in DMAT/hj
            CALL READNS(LUDENS,DMAT)
         ELSE
            CALL WRIDNS(LUDENS,DMAT)
            CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                  NPOS,INTFLG,IPRSCF,WORK(KFREE),LFREE)
         ENDIF
C
         IF(DOTSC) CALL TSCORR(FOCK,WORK(KFREE),LFREE,IPRSCF)
         TDF2 = SECOND() - TDF2
C
C        Write current two-electron Fock matrix to file
C        ==============================================
C
         CALL WRIFCK(LUFCK2,FOCK,NFMAT)
#if defined (MOD_DFT)
C        Construct Kohn-Sham contribution
C        ================================
         IF(DFTADD) THEN
C        NFMAT=1,ISYMOP=1,IHRMOP=1,DOERG=T,DOGRD=F,DOLND=F,DOATR=F
             CALL DFTEXC(DMAT,DUM,FOCK,1,1,1,.TRUE.,.FALSE.,
     &                 .FALSE.,.FALSE.,WORK(KFREE),LFREE,IPRSCF)
         ENDIF
#endif
C
C        Solvent contribution
C        ====================
C
         IF(SOLVEN) THEN
           WRITE(LUPRI,'(A)') '* Calling SOLFCK'
           NF = 1
           CALL SOLFCK(FOCK,DMAT,NF,ESOLVE,ESOLVN,
     &                 WORK(KFREE),LFREE,IPRSOL)
         ENDIF
C
C        Read one-electron Fock matrix from file into CMO
C        ================================================
C
         CALL REAFCK(LUFCK1,CMO,.FALSE.,1)
C        
C        Calculate correction according to the One-center app. models 
C        ============================================================
C
         IF (ONECAP.AND. .NOT.SMLV1C) THEN
              CALL CATCORR(CMO,FOCK,DMAT,WORK(KFREE),LFREE,IPRSCF)
         END IF
C
C
C        Calculate total energy
C        ======================
C
         CALL ERGCAL(CMO,FOCK,DMAT,WORK(KFREE),LFREE)
C
C        Add one-electron Fock matrix to two-electron matrix
C        ===================================================
C
         CALL DAXPY(N2BBASXQ,D1,CMO,1,FOCK,1)
C
C        Add active two-electron Fock matrices
C
         DO IOPEN = 1,NOPEN
C
C           The factor in only DF and not D2*DF because of the
C              multiplication with D2 in TWOFCK.
C
            CALL DAXPY(N2BBASXQ,DF(IOPEN),
     &                 FOCK(1+N2BBASXQ*IOPEN),1,FOCK,1)
         END DO
C
         IF(IPRSCF.GE.5) THEN
            IF (NISHT .GT. 0) THEN
               CALL HEADER('Total FD matrix',-1)
               CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                     NZ,IPQTOQ(1,0),LUPRI)
            END IF
            DO IOPEN = 1,NOPEN
               CALL HEADER('Total FV matrix',-1)
#if defined (SYS_SX)
               WRITE(LUPRI,'(1X,A,I2)') 'Matrix number ',IOPEN
#else
               WRITE(LUPRI,'(X,A,I2)') 'Matrix number ',IOPEN
#endif
               CALL PRQMAT(FOCK(1+N2BBASXQ*IOPEN),NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            END DO
         END IF
C
C        Solve DF eigenvalue problem
C        ===========================
C        Input:  FOCK - Total Fock matrix on AO-basis
C                DMAT - Density matrix
C                BMAT - DIIS matrix
C        Output: CMO  - MO-coefficients for current iteration
C
         ICHKCNV = 1
         CALL DFSOLV(CMO,EIG,IBEIG,FOCK,DMAT,BMAT,ICHKCNV,
     &               OVLSEL,WORK(KFREE),LFREE)
        IF(IPRSCF.GE.1) THEN
          CPUTID = SECTID(TDF2)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* Two-electron Fock matrix  (CPU):',CPUTID
          CPUTID = SECTID(TDDG)
          IF(NZ.EQ.1) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Real diagonalization      (CPU):',CPUTID
          IF(NZ.EQ.2) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Complex diagonalization   (CPU):',CPUTID
          IF(NZ.EQ.4) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Quaternion diagonalization(CPU):',CPUTID
        ENDIF
        CALL GETTIM(CTIDEND,WTIDEND)
        CTIDTOT = CTIDEND - CTIDSTR
        WTIDTOT = WTIDEND - WTIDSTR
        ITRSCF(INTFLG) = ITRSCF(INTFLG) + 1
        SCFTID(INTFLG) = SCFTID(INTFLG) + WTIDTOT
        CPUTID = SECTID(WTIDTOT)
        WRITE(LUPRI,'(A,A12)')
     &    '>>> Total wall time:',CPUTID
        CPUTID = SECTID(CTIDTOT)
        WRITE(LUPRI,'(A,A12)')
     &    '>>> Total CPU-time :',CPUTID
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +    '########## END ITERATION NO.',NITER,' ##########',
     &    DAYTID
        WRITE(LUCYCL,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID(1:10)
        CALL FLSHFO(LUCYCL)
        WRITE(LUPRI, 1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,DHF_INTTYP,DAYTID(1:10)
        CALL FLSHFO(LUPRI)
C
C       switch off onecap approximation ?
        IF (ONECOFF) THEN
           ONECOFF = .FALSE.
C          do continue even if it was converged with ONECAP :
           DHFEXIT = .FALSE.
C          reset intbuf so no old Fock matrices are used :
           INTBUF = 0
           WRITE(LUPRI,'(/A,1P,D15.6/2A/)') 
     &  ' ** INFO : ONECAP. Wave function converged beyond : ONECNV = ',
     &     ONECNV,'          ONECAP will be turned off, and integrals',
     &     ' evaluated on the fly !'
        ELSE
           INTBUF = INTFLG
        END IF
C
        IF(DHFEXIT) GOTO 30
 10   CONTINUE
      WRITE(LUPRI,'(/A)')
     &' ** Exit SCF because maximum number of iterations reached.'
C
C     **********************************
C     ***** End of SCF-iterations ******
C     **********************************
C
 30   CONTINUE
C
C     Reorder converged orbitals if requested
C
      IF (LMOORF) THEN
C
C        It looks like CMO is trashed at this point, so we re-read
C        the orbitals from DFCOEF:
C
         CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,EXXX,14,
     &               WORK,KFREE,LFREE)
C
         CALL REORD(CMO,EIG,IBEIG,IMOORF,IREORF,MXREORD,
     &        WORK(KFREE),LFREE)
C
C        Make sure that reordered coefficients are written back
C        to file.
C
          CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG)
      END IF
C
C     Phase-adjustment of coefficients, if requested...
C
      IF(PHCOEF) THEN
        DO I = 1,NFSYM
          CALL PHATRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ)
        ENDDO
        CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG)
      ENDIF
C
C     Close necessary files
C
      CLOSE(LUCOEF,STATUS='KEEP')
      CLOSE(LUCMOS,STATUS='KEEP')
      CLOSE(LUDENS,STATUS='KEEP')
      CLOSE(LUFOCK,STATUS='KEEP')
      CLOSE(LUFCK1,STATUS='KEEP')
      CLOSE(LUFCK2,STATUS='KEEP')
      CLOSE(LUCYCL,STATUS='KEEP')
      CLOSE(LUEVEC,STATUS='KEEP')
C
      CALL QEXIT('SCFCY1')
      RETURN
 1010 FORMAT('It. ',I3,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A9,3X,A10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ergcal */
      SUBROUTINE ERGCAL(F1,F2,DMAT,WORK,LWORK)
C****************************************************************************
C     PURPOSE:
C        Calculate Dirac-Hartree-Fock energy
C
C      E = D(IJ)(*)(F1(IJ)+0.5F2(IJ))
C
C     Last revision: 971022-jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2=2.0D0, DP5 = 0.5D00,
     &          DP25 = 0.25D00)
#include "mxcent.h"
C
#include "nuclei.h"
#include "frame.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbpsi.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#if defined (MOD_DFT)
#include "dfterg.h"
#endif
      LOGICAL LVCON
      CHARACTER M(4)*2,FMT*6,MXFORM*6
      DIMENSION F1(N2BBASX,NZ),F2(N2BBASX,NZ,NFMAT),
     &          DMAT(N2BBASX,NZ,NFMAT),
     &          E(4,2,0:MXOPEN),WORK(LWORK)
      DATA M/'AR','AI','BR','BI'/
      DATA INUM /0/
      SAVE INUM
C
      CALL QENTER('ERGCAL')
      INUM = INUM + 1
C
      LVCON = DOLVC.AND.
     &        ( (INTFLG.EQ.3) .OR. (ONECAP.AND.INTV1C.EQ.2) )
      CALL DZERO(E,8*(MXOPEN+1))
      ERGBUF  = DHFERG
      E1PART  = D0
      E2PART  = D0
C
C
      DO ISHELL = 0,NOPEN
C
C        Factor D2 comes from Time-Reversal Symmetry.
C
         DFRAC1 = D2*DF(ISHELL)
C
C        Factor DP25 = DP5*DP5 comes from the formula and the
C           fact that all Fock matrices were multiplied with 2
C           in TWOFCK.
C
         DFRAC2 = DP25*DFRAC1*DFRAC1*DA(ISHELL)
C...     Intra-shell contributions
         DO IZ = 1,NZ
            E(IZ,1,ISHELL) =
     &         DFRAC1*DDOT(N2BBASX,DMAT(1,IZ,ISHELL+1),1,F1(1,IZ),1)
            E(IZ,2,ISHELL) =
     &         DFRAC2*DDOT(N2BBASX,DMAT(1,IZ,ISHELL+1),1,
     &                               F2(1,IZ,ISHELL+1),1)
         END DO
C...     Inter-shell contribution
         DO JSHELL = 0,ISHELL-1
C
C           The factor is 0.5 * (2*f_{o1})*(2*f_{o2})
C           Divide by 2 (from TWOFCK)
C           Total factor is DF(ISHELL)*DF(JSHELL)
C
            DFRAC1 = DF(ISHELL)*DF(JSHELL)
            DO IZ = 1,NZ
               ETEMP =
     &            DFRAC1*DDOT(N2BBASX,DMAT(1,IZ,JSHELL+1),1,
     &                                  F2(1,IZ,ISHELL+1),1)
               E(IZ,2,ISHELL) = E(IZ,2,ISHELL) + ETEMP
               E(IZ,2,JSHELL) = E(IZ,2,JSHELL) + ETEMP
            END DO
         END DO
      END DO
C...  If DFT correct for contribution from Kohn-Sham potential
C...  to two-electron energy
      IF(DODFT) THEN
        E(1,2,0) = E(1,2,0) - ERRKS
      ENDIF
      DO IZ = 1,NZ
        DO IOPEN = 0,NOPEN
           E1PART = E1PART + E(IZ,1,IOPEN)
           E2PART = E2PART + E(IZ,2,IOPEN)
        END DO
      END DO
C
      DHFERG = POTNUC
      IF(SOLVEN) THEN
        E2PART = E2PART - ESOLVE
        ERGSOL = ESOLVE + ESOLVN
        DHFERG = DHFERG + ERGSOL
      ENDIF
      IF(LVCON) THEN
         IF (ONECAP .AND. (INTV1C.EQ.2)) THEN
            JPRINT=0
         ELSE
            JPRINT=IPRSCF
         END IF
         CALL LVCORR(DMAT,WORK,LWORK,JPRINT)
         DHFERG = DHFERG + CORRLV
      ENDIF
      ELERGY = E1PART + E2PART
      IF(DODFT) THEN
        ELERGY = ELERGY + EDFTY
      ENDIF
      DHFERG = DHFERG + ELERGY
C
      IF(IPRSCF.GE.1) THEN
         CALL HEADER('Output from ERGCAL',-1)
         ERGMAX = MAX(ABS(E1PART),ABS(E2PART))*10
         FMT = MXFORM(ERGMAX,20)
            WRITE(LUPRI,'(A,'//FMT//')')
     +       'Total electronic energy            :  ',ELERGY,
     +       'Nuclear potential energy           :  ',POTNUC
         IF(LVCON) WRITE(LUPRI,'(A,'//FMT//')')
     +       'SS Coulombic correction            :  ',CORRLV
         IF(SOLVEN) WRITE(LUPRI,'(A,'//FMT//')')
     +       'Solvent energy                     :  ',ERGSOL
         WRITE(LUPRI,'(A,'//FMT//')')
     +       'Total energy                       :  ',DHFERG
         IF(IPRSCF.GE.2) THEN
            WRITE(LUPRI,'(6X,A,'//FMT//')')
     &         '* 1-electron energy          :  ',E1PART
            IF(NZ.GT.1) THEN
               DO 20 IZ = 1,NZ
               DO IOPEN = 0,NOPEN
                  WRITE(LUPRI,'(6X,A,A2,A,I2,A,'//FMT//')')
     +               '            from ',M(IZ),'(Shell ',
     +               IOPEN,'):  ',E(IZ,1,IOPEN)
               END DO
   20          CONTINUE
            ENDIF
            WRITE(LUPRI,'(6X,A,'//FMT//')')
     +         '* 2-electron energy          :  ',E2PART
            IF(NZ.GT.1) THEN
               DO 30 IZ = 1,NZ
               DO IOPEN = 0,NOPEN
                  WRITE(LUPRI,'(6X,A,A2,A,I2,A,'//FMT//')')
     +               '            from ',M(IZ),'(Shell ',
     +               IOPEN,'):  ',E(IZ,2,IOPEN)
               END DO
   30          CONTINUE
            ENDIF
            IF (DODFT) WRITE(LUPRI,'(6X,A,'//FMT//')')
     +         '* DFT xc energy              :  ',EDFTY
         END IF
      ENDIF
C
      CALL QEXIT('ERGCAL')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intcon */
      SUBROUTINE INTCON(INTFLG,INTOLD,INTDEF,
     &                  CNV,CNVINT,ITR,ITRINT,INTTYP)
C***********************************************************************
C
C     Check what integrals contribute in a given SCF-iteration
C       CNVINT has precedence over ITRINT !!!
C       The information is bitpacked in INTDEF:
C         001   LL-contribution
C         010   SL-contribution
C         100   SS-contribution
C       and all combinations.....
C
C     Written by T.Saue Mar 09 1995
C     Last revision: Aug 24 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      LOGICAL LBIT
      CHARACTER INTTYP*9
      DIMENSION CNVINT(2),ITRINT(2)
      INTTYP = '         '
C
C     LL-integrals
C     ============
C
      IF(LBIT(INTDEF,1)) THEN
        ILLINT = 1
      ELSE
        ILLINT = 0
      ENDIF
C
C     SL-integrals
C     ============
C
      IF(LBIT(INTDEF,2)) THEN
        ISLINT = 1
        IF(.NOT.LBIT(INTOLD,2)) THEN
          IF(ABS(CNV).GT.CNVINT(1)) ISLINT = 0
          IF(ITR     .LT.ITRINT(1)) ISLINT = 0
        ENDIF
      ELSE
        ISLINT = 0
      ENDIF
C
C     SS-integrals
C     ============
C
      IF(LBIT(INTDEF,3)) THEN
        ISSINT = 1
        IF(.NOT.LBIT(INTOLD,3)) THEN
          IF(ABS(CNV).GT.CNVINT(2)) ISSINT = 0
          IF(ITR     .LT.ITRINT(2)) ISSINT = 0
        ENDIF
      ELSE
        ISSINT = 0
      ENDIF
C
      IF(ILLINT.EQ.1) INTTYP(1:3) = 'LL '
      IF(ISLINT.EQ.1) INTTYP(4:6) = 'SL '
      IF(ISSINT.EQ.1) INTTYP(7:9) = 'SS '
      INTFLG= ILLINT + 2*ISLINT + 4*ISSINT
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smogen */
      SUBROUTINE SMOGEN(CMO,TMAT,WORK,LWORK)
C***********************************************************************
C
C     Multiply overlap matrix onto occupied vectors and write to file;
C     this forms the basis for overlap selection
C
C     Written by T.Saue Nov 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C  DCBBAS: N2BBASX,NFBAS
C  DCBORB: NOCC
C  DHFCAL: NSMOTQ
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION CMO(*),TMAT(*),WORK(LWORK)
C
      CALL QENTER('SMOGEN')
#include "memint.h"
      NBFDIM = 0
      DO I = 1,NFSYM
        NBFDIM = MAX(NBFDIM,NOCC(I)*NFBAS(I,0))
      ENDDO
      NBFDIM = NBFDIM*NZ
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMO,NSMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,NBFDIM,WORK,KFREE,LFREE)
      CALL SMOGE1(WORK(KSMO),CMO,TMAT,WORK(KSMAT),WORK(KBUF),
     &            WORK(KFREE),LFREE)
      CALL MEMREL('SMOGEN.SMOGE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SMOGEN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smoge1 */
      SUBROUTINE SMOGE1(SMO,CMO,TMAT,SMAT,BUF,WORK,LWORK)
C***********************************************************************
C
C     Multiply overlap matrix onto occupied vectors and write to file;
C     this forms the basis for overlap selection
C
C     Overlap selection is based on a previous set CPRE of vectors in AO-basis
C     We form the matrix
C          SMO  =  T^(+) * S * CPRE
C     where S is the overlap matrix and T the MO-transformation matrix.
C     In each SCF-iteration we form the matrix
C          SEL  = CMO^(+) * SMO
C     using the current CMO-coefficients in MO-basis.
C
C     Written by T.Saue Nov 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
C Used from COMMON blocks:
C  PAMIOU: LUOVLP
C  DCBBAS:
C  DCBORB:
C  DGROUP: NFSYM
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbdhf.h"
      DIMENSION SMO(*),CMO(*),TMAT(*),SMAT(*),BUF(*),WORK(LWORK)
C
C     Get overlap matrix
C
      CALL GTOVLX(SMAT,SSMTRC,WORK,LWORK)
C
      DO I = 1,NFSYM
      IF(NOCC(I).NE.0) THEN
        IOFF = ICMOQ(I) + NPSH(I)*NFBAS(I,0)
        CALL QGEMM(NFBAS(I,0),NOCC(I),NFBAS(I,0),D1,
     &             'N','N',IPQTOQ(1,0),SMAT(I2BASX(I,I)+1),
     &                     NTBAS(0),NTBAS(0),1,
     &             'N','N',IPQTOQ(1,0),CMO(IOFF+1),
     &                     NFBAS(I,0),NORB(I),NZ,
     &        D0,IPQTOQ(1,0),BUF,NFBAS(I,0),NOCC(I),NZ)
        CALL QGEMM(NORB(I),NOCC(I),NFBAS(I,0),D1,
     &             'H','N',IPQTOQ(1,0),TMAT(I2TMT(I)+1),
     &                     NFBAS(I,0),NORB(I),NZT,
     &             'N','N',IPQTOQ(1,0),BUF,NFBAS(I,0),NOCC(I),NZ,
     &        D0,IPQTOQ(1,0),SMO(ISMOQ(I)+1),NORB(I),NOCC(I),NZ)
        IF(IPRSCF.GE.6) THEN
          CALL HEADER('SMOGE1: SMO-matrix for overlap selection',-1)
          IF (NFSYM.GT.1) WRITE(LUPRI,'(A,I2)') ' Fermion symmetry',I
          CALL PRQMAT(SMO(ISMOQ(I)+1),NORB(I),NOCC(I),NORB(I),NOCC(I),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
      ENDIF
      ENDDO
C
C     Write to file
C
      CALL OPNFIL(LUSMOS,'DFSMOS','UNKNOWN','SMOGE1')
      REWIND (LUSMOS)
      CALL WRITT(LUSMOS,NSMOTQ,SMO)
      CLOSE(LUSMOS,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smosel */
      SUBROUTINE SMOSEL(EIG,IBEIG,CFU,WORK,LWORK)
C***********************************************************************
C
C     Select vectors based on overlap.
C     Input is coefficients in MO-basis.
C
C     Written by T.Saue Nov 11 1995
C     Last revision May 2002 hjaaj: include EIG and IBEIG
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG(*),IBEIG(*),CFU(*),WORK(LWORK)
C Used from COMMON blocks:
C  DHFCAL: NSMOTQ
C
#include "dcbdhf.h"
      CALL QENTER('SMOSEL')
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KSMO,NSMOTQ,WORK,KFREE,LFREE)
      CALL SMOSE1(EIG,IBEIG,CFU,WORK(KSMO),WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('SMOSEL.SMOSE1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('SMOSEL')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smose1 */
      SUBROUTINE SMOSE1(EIG,IBEIG,CFU,SMO,WORK,LWORK)
C***********************************************************************
C
C     Select vectors based on overlap
C
C     Written by T.Saue Nov 11 1995
C     Last revision May 2002 hjaaj: include EIG and IBEIG
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG(*),IBEIG(*),CFU(*),SMO(*),WORK(LWORK)
C
C Used from COMMON blocks:
C  DCBORB: NSMOTQ,NOCC
C  DGROUP: NFSYM
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
      CALL QENTER('SMOSE1')
#include "memint.h"
C
C     Read vectors to select from
C
      CALL OPNFIL(LUSMOS,'DFSMOS','OLD','SMOSE1')
      CALL READT(LUSMOS,NSMOTQ,SMO)
      CLOSE(LUSMOS,STATUS='KEEP')
C
      DO I = 1,NFSYM
      IF(NOCC(I).NE.0) THEN
C       Memory allocation
        CALL MEMGET('REAL',KEBUF,NESH(I),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEB,NESH(I),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KCBUF,NESH(I)*NORB(I)*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSBUF,NESH(I)*NOCC(I)*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRSEL,NESH(I),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KCSEL,NESH(I),WORK,KFREE,LFREE)
        CALL SMOSE2(I,EIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &              CFU(I2ORBT(I)+1),SMO(ISMOQ(I)+1),
     &              WORK(KEBUF),WORK(KIBEB),WORK(KCBUF),WORK(KSBUF),
     &              WORK(KRSEL),WORK(KCSEL),ISELM,
     &              NORB(I),NOCC(I),NESH(I),WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('SMOSE1.SMOSE2',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
      ENDDO
      CALL QEXIT('SMOSE1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smose2 */
      SUBROUTINE SMOSE2(IFRP,EIG,IBEIG,CFU,SMO,EBUF,IBEBUF,CBUF,SBUF,
     &                  RCSEL,ICSEL,ISELM,KNORB,KNOCC,KNESH,WORK,LWORK)
C***********************************************************************
C
C     Select vectors based on overlap
C
C     Written by T.Saue Nov 11 1995
C     Major revisions: Jan 7 1998 - jth NORB(IFRP) etc. is now parameter
C                      May 2002 hjaaj: EIG and IBEIG added
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
C Used from COMMON blocks:
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbdhf.h"
C
C     KNORB = NORB(IFRP)
C     KNOCC = NOCC(IFRP)
C     KNESH = NESH(IFRP)
C
      CHARACTER FSYM*2
      DIMENSION EIG(KNORB),IBEIG(KNORB),
     &          CFU(KNORB,KNORB,NZ),
     &          SMO(KNORB,KNOCC,NZ),
     &          SBUF(KNESH,KNOCC,NZ),
     &          EBUF(KNESH),IBEBUF(KNESH),
     &          CBUF(KNORB,KNESH,NZ),
     &          RCSEL(KNESH),ICSEL(KNESH),
     &          WORK(LWORK)
C
C     Form product CFU * (SMAT * COCC)
C
      CALL QGEMM(NESH(IFRP),NOCC(IFRP),NORB(IFRP),D1,
     &       'H','N',IPQTOQ(1,0),CFU(1,NPSH(IFRP)+1,1),
     &                NORB(IFRP),NORB(IFRP),NZ,
     &       'N','N',IPQTOQ(1,0),SMO,NORB(IFRP),NOCC(IFRP),NZ,
     &       D0,IPQTOQ(1,0),SBUF,NESH(IFRP),NOCC(IFRP),NZ)
      IF(IPRSCF.GE.10) THEN
        CALL HEADER('SMOSE1: SBUF-matrix',-1)
        CALL PRQMAT(SBUF,NESH(IFRP),NOCC(IFRP),NESH(IFRP),NOCC(IFRP),
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      IF(NZ.EQ.1) THEN
        DO J = 1,NOCC(IFRP)
          DO I = 1,NESH(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
          ENDDO
        ENDDO
      ELSEIF(NZ.EQ.2) THEN
        DO J = 1,NOCC(IFRP)
          DO I = 1,NESH(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
     &                  + SBUF(I,J,2)*SBUF(I,J,2)
          ENDDO
        ENDDO
      ELSEIF(NZ.EQ.4) THEN
        DO J = 1,NOCC(IFRP)
          DO I = 1,NESH(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
     &                  + SBUF(I,J,2)*SBUF(I,J,2)
     &                  + SBUF(I,J,3)*SBUF(I,J,3)
     &                  + SBUF(I,J,4)*SBUF(I,J,4)
          ENDDO
        ENDDO
      ENDIF
      IF(IPRSCF.GE.6) THEN
        CALL HEADER('SMOSE1: Squared-matrix',-1)
        CALL PRQMAT(SBUF,NESH(IFRP),NOCC(IFRP),NESH(IFRP),NOCC(IFRP),
     &              1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Perform selection
C
      N2DIM = NORB(IFRP)*NESH(IFRP)
      DO IZ = 1,NZ
        CALL DCOPY(N2DIM,CFU(1,NPSH(IFRP)+1,IZ),1,CBUF(1,1,IZ),1)
      ENDDO
      CALL DCOPY(KNESH,EIG(NPSH(IFRP)+1),1,EBUF,1)
      CALL ICOPY(KNESH,IBEIG(NPSH(IFRP)+1),1,IBEBUF,1)
      ISELM = 0
      DO J = 1,NOCC(IFRP)
        JJ = NPSH(IFRP) + J
        ISEL = IDAMAX(NESH(IFRP),SBUF(1,J,1),1)
        ICSEL(J) = ISEL
        RCSEL(J) = SBUF(ISEL,J,1)
        ISELM = MAX(ISELM,ISEL)
C       Zero out the selected overlap
        CALL DCOPY(NOCC(IFRP),D0,0,SBUF(ISEL,1,1),NESH(IFRP))
C       Copy buffer into vector
        DO IZ = 1,NZ
          CALL DCOPY(NORB(IFRP),CBUF(1,ISEL,IZ),1,CFU(1,JJ,IZ),1)
        ENDDO
C       Copy eigenvalue and (approx.) boson symmetry to new place:
        EIG(JJ) = EBUF(ISEL)
        IBEIG(JJ) = IBEBUF(ISEL)
      ENDDO
C
C     Copy in missing vectors
C
      IF(ISELM.GT.NOCC(IFRP)) THEN
        JSTR = NOCC(IFRP)+1
        DO J = JSTR,ISELM
          JJ = NPSH(IFRP) + J
          ISEL = IDAMAX(ISELM,SBUF(1,J,1),1)
          ICSEL(J) = ISEL
          RCSEL(J) = SBUF(ISEL,J,1)
C         Zero out the selected overlap
          CALL DCOPY(NOCC(IFRP),D0,0,SBUF(ISEL,1,1),NESH(IFRP))
C         Copy buffer into vector
          DO IZ = 1,NZ
            CALL DCOPY(NORB(IFRP),CBUF(1,ISEL,IZ),1,CFU(1,JJ,IZ),1)
          ENDDO
C         Copy eigenvalue and (approx.) boson symmetry to new place:
          EIG(JJ) = EBUF(ISEL)
          IBEIG(JJ) = IBEBUF(ISEL)
        ENDDO
      ENDIF
C
C     Print section
C
      WRITE(FSYM,'(I2)') IFRP
      CALL HEADER('SMOSEL: Reordered vectors in fermion ircop '
     &            //FSYM,-1)
      NVEC = 0
      DO J = 1,NOCC(IFRP)
      IF(ICSEL(J).NE.J) THEN
        NVEC = NVEC + 1
        WRITE(LUPRI,'(I5,A,I5,A,F20.10)') ICSEL(J),' --> ',J,
     &                                    ' OVERLAP ',RCSEL(J)
      ENDIF
      ENDDO
      IF(NVEC.EQ.0) WRITE(LUPRI,'(A)') 'No reordering.'
      IF(IPRSCF.GE.2) THEN
        WRITE(LUPRI,'(20I6)') (ICSEL(I),I=1,NOCC(IFRP))
        WRITE(LUPRI,'(10F12.3)') (RCSEL(I),I=1,NOCC(IFRP))
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck denorb */
      SUBROUTINE DENORB(DMAT,JVEC,NVEC,IFRP,CMO,IPRINT)
C***********************************************************************
C
C     Construct density matrix for NVEC orbitals of fermion corep IFRP
C     starting with orbital JVEC (with respect to NPSH(IFRP))
C
C     Written by T.Saue Oct 22 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,D0=0.0D0)
      DIMENSION DMAT(*),CMO(*)
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      JOFF = NPSH(IFRP) + JVEC
      CALL DZERO(DMAT,N2BBASXQ)
      CALL DENST1(DMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),
     &              CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &              JOFF,NVEC,NFBAS(IFRP,0))
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL TITLER('Output from DENORB','*',103)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',IFRP
        CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1),NFBAS(IFRP,0),
     &              NFBAS(IFRP,0),NTBAS(0),NTBAS(0),
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck oneerg */
      SUBROUTINE ONEERG(EIG,IBEIG)
C**********************************************************************
C     PURPOSE:
C        Calculate one-electron energy as sum of eigenvalues
C
C        Written by T.Saue Jan 28 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2=2.0D0)
C
#include "mxcent.h"
#include "frame.h"
#include "nuclei.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbdhf.h"
      CHARACTER FMT*6,MXFORM*6
      DIMENSION EIG(*)
C
      CALL QENTER('ONEERG')
C
      ERGBUF = DHFERG
      E1PART = D0
      EA1PART = D0
      DO IFRP = 1,NFSYM
        IF(NOCC(IFRP).GT.0) THEN
           IF (NISHT .GT. 0) THEN
              E1PART = E1PART + D2*DSUM(NISH(IFRP),
     &                         EIG(IORB(IFRP)+NPSH(IFRP)+1),1)
           END IF
           IF (NASHT .GT. 0) THEN
             IDX = IORB(IFRP)+NPSH(IFRP)+NISH(IFRP)+1
             DO IOPEN = 1,NOPEN
                E1PART = E1PART + D2*DF(IOPEN)*DSUM(NACSH(IFRP,IOPEN),
     &                    EIG(IDX),1)
                IDX = IDX + NACSH(IFRP,IOPEN)
             END DO
           END IF
        ENDIF
      ENDDO
      ELERGY = E1PART
      DHFERG = ELERGY + POTNUC
C
      IF(IPRSCF.GE.1) THEN
        Call HEADER('Output from ONEERG',-1)
        ERGMAX = ABS(E1PART+EA1PART)*10
        FMT = MXFORM(ERGMAX,20)
        IF (IPRSCF .GE. 2) THEN
           WRITE(LUPRI,'(A,'//FMT//')')
     &     'Inactive energy           :  ',E1PART
        END IF
        WRITE(LUPRI,'(A,'//FMT//')')
     +     'Electronic energy         :  ',ELERGY,
     +     'Nuclear potential energy  :  ',POTNUC,
     +     'Total energy              :  ',DHFERG
      ENDIF
C
      CALL QEXIT('ONEERG')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reord */
      SUBROUTINE REORD(CMO,EIG,IBEIG,IMOORD,IREORD,MXREORD,WORK,LWORK)
C**********************************************************************
C     PURPOSE:
C        Reorder orbitals
C
C        Written by J. Thyssen Nov 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION WORK(LWORK),CMO(*),EIG(*),IBEIG(*)
      DIMENSION IMOORD(MXREORD,2),IREORD(2)
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
C
#include "memint.h"
C
      CALL QENTER('REORD')
C

      DO 10 I = 1,NFSYM
         IF (IREORD(I).EQ.0) GOTO 10
         IF (IPRSCF.GE.2) THEN
            CALL HEADER('REORD: before reordering',-1)
            IF (NFSYM.GT.1) WRITE(LUPRI,'(5X,A,I2/)')
     &         'for fermion ircop',I
            WRITE(LUPRI,'(3X,A)')
     &         '* Electron eigenvalues and boson classification ...'
            WRITE(LUPRI,'(3X,F15.8,I5)')
     &         (EIG(J),IBEIG(J),J=1+IORB(I)+NPSH(I),IORB(I)+NORB(I))
         END IF
         IF (IPRSCF.GE.10) THEN
            WRITE(LUPRI,'(3X,A)') '* Eigenvectors...'
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Save coefficients
C
         CALL MEMGET('REAL',KTCMO,NCMOQ(I),WORK,KFREE,LFREE)
         CALL DCOPY(NCMOQ(I),CMO(1+ICMOQ(I)),1,WORK(KTCMO),1)
         CALL MEMGET('REAL',KTEIG,NORB(I),WORK,KFREE,LFREE)
         CALL DCOPY(NORB(I),EIG(1+IORB(I)),1,WORK(KTEIG),1)
         CALL MEMGET('INTE',KTIBE,NORB(I),WORK,KFREE,LFREE)
         CALL ICOPY(NORB(I),IBEIG(1+IORB(I)),1,WORK(KTIBE),1)
C
C        Reorder coefficients
C
         CALL REOR1(CMO(1+ICMOQ(I)),WORK(KTCMO),
     &              EIG(1+IORB(I)),WORK(KTEIG),
     &              IBEIG(1+IORB(I)),WORK(KTIBE),
     &              NFBAS(I,0),NORB(I),NZ,
     &              IMOORD(1,I),IREORD(I),NPSH(I),NESH(I))
C
C        Release memory
C
         CALL MEMREL('REORD',WORK,KWORK,KWORK,KFREE,LFREE)
C
         IF (IPRSCF.GE.2) THEN
            CALL HEADER('REORD: after reordering',-1)
            IF (NFSYM.GT.1) WRITE(LUPRI,'(5X,A,I2/)')
     &         'for fermion ircop',I
            WRITE(LUPRI,'(3X,A)')
     &         '* Electron eigenvalues and boson classification ...'
            WRITE(LUPRI,'(3X,F15.8,I5)')
     &         (EIG(J),IBEIG(J),J=1+IORB(I)+NPSH(I),IORB(I)+NORB(I))
         END IF
         IF (IPRSCF.GE.10) THEN
            WRITE(LUPRI,'(3X,A)') '* Eigenvectors...'
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
         END IF
 10   CONTINUE
      CALL FLSHFO(LUPRI)
      CALL QEXIT('REORD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reor1 */
      SUBROUTINE REOR1(CMO,CMOBCK,EIG,EIGBCK,IBEIG,IBEIGBCK,
     &                 NBASI,NORBI,NZ,IORDER,NORDER,
     &                 NPOSITRON,NELECTRON)
C**********************************************************************
C     PURPOSE:
C        Reorder orbitals
C
C        Written by J. Thyssen Nov 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO(NBASI,NORBI,NZ), CMOBCK(NBASI,NORBI,NZ),
     &          EIG(*),EIGBCK(*), IBEIG(*), IBEIGBCK(*),
     &          IORDER(NORDER)
C
C
      CALL QENTER('REOR1')
C
      ICHK_S = 0
      ICHK_D = 0
      DO I = 1,NORDER
         IORB_S = IORDER(I)+NPOSITRON
         IORB_D = I+NPOSITRON
         IF (IORB_D .GT. NELECTRON+NPOSITRON) THEN
            CALL QUIT('REOR1 reordering to non-existent orbital')
         ELSE
            ICHK_S = ICHK_S + IORB_S
            ICHK_D = ICHK_D + IORB_D
C
C           Copy eigenvalue and boson classification
C
            EIG(IORB_D)   = EIGBCK(IORB_S)
            IBEIG(IORB_D) = IBEIGBCK(IORB_S)
C
C           Copy eigenvector
C
            DO IZ = 1,NZ
               CALL DCOPY(NBASI,CMOBCK(1,IORB_S,IZ),1,
     &                          CMO(1,IORB_D,IZ),1)
            END DO
         END IF
      END DO
C
C     Make sure we are not copying more than one orbital to the same 
C     destination,
C     or that IORDER(I) does not belong to 1..NORDER
C     (not 100% test: if user specified e.g. IORDER=2,2,2 this would 
C     not catch it,
C     but that is quite unlikely) /hjaaj May 2002
C
      IF (ICHK_S .NE. ICHK_D) THEN
         CALL QUIT('REOR1: destination orbs. .ne. source orbitals')
      END IF
C
      CALL QEXIT('REOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKMOFK */
      SUBROUTINE MKMOFK2(FMO,FAO,TMAT,WORK,LWORK)
C**********************************************************************
C     PURPOSE:
C        Construct MO Fock matrix from Fock matrices in AO-basis
C     WARNING:
C        This is a small part of the routine MKMOFK distributed
C           with the open shell code.
C
C        Written by J. Thyssen Dec 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00)
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION FMO(*),FAO(N2BBASXQ,NFMAT),TMAT(*),WORK(*)
      CHARACTER CTMP*3
C
#include "memint.h"
C
      CALL QENTER('MKMOFK2')
C
C     Transform FD to MO-basis
C
      DO I = 1,NFSYM
      IF(NTMO(I).GT.0) THEN
        IF(SUB_BL) THEN
          CALL TSUBBL(I,FAO,FMO,TMAT,TMAT,WORK(KFREE),LFREE,IPRSCF)
        ELSE
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NTMO(I),NTMO(I),
     &            FAO(I2BASX(I,I)+1,1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRSCF)
        ENDIF
        IF(IPRSCF.GE.10 ) THEN
          WRITE(CTMP,'(I1)') I
          CALL HEADER('DFSOLV: D Fock matrix in MO-basis'//
     &                ', corep '//CTMP,-1)
          CALL PRQMAT(FMO(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &                NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
      END IF
      END DO
      CALL QEXIT('MKMOFK2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck EXTRSB */
      SUBROUTINE EXTRSB (IMODE,I,NORB_SUB_L,M_SUB_BL_L,N_SUB_BL_L,
     &                   NORB_SUB_R,M_SUB_BL_R,N_SUB_BL_R,
     &                   ISUB,JSUB,ICOMP,JCOMP,IZS,IZE,
     &                   FSUB,MSUB,NSUB,NZS,
     &                   FFUL,MFUL,NFUL,NZ)
C**********************************************************************
C     PURPOSE:
C        Extract or insert blocks out of (into) a full matrix 
C
C        General parameters :
C        IMODE                 - 1 : extract, 2 : insert
C        I                     Parity (fermion irrep)
C        IZS                   First quaternion unit
C        IZE                   Last quaternion unit
C        NZ                    Number of quaternion units
C        NZS                   Number of active quaternion units
C
C        Left index :
C        NORB_SUB_L            Block dimensions left index
C        M_SUB_BL_L            First dimension of NORB_SUB_L
C        N_SUB_BL_L            Number of blocks left index
C        ISUB                  Active block left (zero : all blocks)
C        ICOMP                 Active component left (zero : all comp.)
C
C        Right index :
C        NORB_SUB_R            Block dimensions right index
C        M_SUB_BL_R            First dimension of NORB_SUB_R
C        N_SUB_BL_R            Number of blocks right index
C        JSUB                  Active block right (zero : all blocks)
C        JCOMP                 Active component right (zero : all comp.)
C
C        Input and output :
C        FSUB                  Sub matrix
C        FFUL                  Full matrix
C
C        Written for for boson symmetry (august 2000
C        Generalized to subblocks (november 2003)
C        L. Visscher
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FSUB(MSUB,NSUB,NZS),FFUL(MFUL,NFUL,NZ)
      DIMENSION NORB_SUB_L(M_SUB_BL_L,2,0:2)
      DIMENSION NORB_SUB_R(M_SUB_BL_R,2,0:2)
      LOGICAL SUBI, SUBJ
C
C     The full matrix is ordered as (function index, component type, sub block),
C     the submatrix should be ordered as (function index, sub blocK)
C
C     if isub is zero we extract blocks from all sub blocks with parity I
C     if isub has a non-zero value we just extract the block of this type.
C
      IF (ISUB.EQ.0) THEN
         ISS = 1
         ISE = N_SUB_BL_L
      ELSE
         ISS = ISUB
         ISE = ISUB
      ENDIF
C
      IF (JSUB.EQ.0) THEN
         JSS = 1
         JSE =  N_SUB_BL_R
      ELSE
         JSS = JSUB
         JSE = JSUB
      ENDIF
C
C     Extract the subblock by looping over both indices in the right order
C     while using different counting variables for the full and the sub matrices
C     Whenever we encounter a block that fulfills the criteria 
C     we copy it into or out from the subbblock.
C
      DO IZ = IZS, IZE
        JXF = 0
        JXS = 0
C
C       Loop over all blocks for the right index.
C
        DO JS = 1, N_SUB_BL_R
C
C         The small component functions come first, loop from S to L
C
          DO JC = 2, 1, -1
            SUBJ = JC.EQ.JCOMP.AND.JS.GE.JSS.AND.JS.LE.JSE
            DO 100 JX = 1, NORB_SUB_R(JS,I,JC)
              JXF = JXF + 1
              IF (SUBJ) THEN
C                Only when we are in a selected block we update this pointer
                 JXS = JXS + 1
              ELSE
C
C             We are not going to run the loops for the left index 
C             and have already updated JXF so we might as well skip 
C             the inner loops
C
                 GOTO 100
              ENDIF
              IXF = 0
              IXS = 0
C
C             Loop over all blocks for the left index.
C
              DO IS = 1, N_SUB_BL_L
C
C               Loop from S to L for the left index.
C
                DO IC = 2, 1, -1
                  SUBI = IC.EQ.ICOMP.AND.IS.GE.ISS.AND.IS.LE.ISE
                  DO IX = 1, NORB_SUB_L(IS,I,IC)
                    IXF = IXF + 1
C
C                   Check the criteria and fill in
C
                    IF (SUBI.AND.SUBJ) THEN
                       IXS = IXS + 1
                       IF (IMODE.EQ.1) THEN
                          FSUB(IXS,JXS,IZ) = FFUL(IXF,JXF,IZ)
                       ELSEIF (IMODE.EQ.2) THEN
                          FFUL(IXF,JXF,IZ) = FSUB(IXS,JXS,IZ)
                       ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
 100        CONTINUE
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck EXTRSBR */
      SUBROUTINE EXTRSBR (IMODE,I,ISYM,JSYM,ICOMP,JCOMP,IZS,IZE,
     &                    FSUB,MSUB,NSUB,NZS,
     &                    FFUL,MFUL,NFUL)
C**********************************************************************
C     PURPOSE:
C        Extract or insert a sub block out of the full matrix 
C        that is ordered on boson irreps for the left and on small, 
C        large or positron, electron on the right index.
C
C        Written by L. Visscher
C        September 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION FSUB(MSUB,NSUB,NZS),FFUL(MFUL,NFUL,NZ)
      LOGICAL SUBI, SUBJ
C
C     Calculate number of boson symmetries : the full matrix is ordered
C     as (function index, component type, boson irrep), the submatrix
C     should be ordered as (function index, boson irrep)
C
      NBRP = 4 / NZ
C
C     When isym is zero we extract blocks from all boson irreps that
C     have the right parity, if isym has a non-zero value we just extract
C     the block of this symmetry.
C
      IF (ISYM.EQ.0) THEN
         ISS = 1
         ISE = 4 / NZ
      ELSE
         ISS = ISYM
         ISE = ISYM
      ENDIF
C
      IF (JSYM.EQ.0) THEN
         JSS = 1
         JSE = 4 / NZ
      ELSE
         JSS = JSYM
         JSE = JSYM
      ENDIF
C
C     Extract the subblock by looping over both indices in the 
C     right order while using different counting variables for the 
C     full and the sub matrices. Whenever we encounter a block that 
C     fulfills the criteria we copy it into or out from the subbblock.
C
      DO IZ = IZS, IZE
        JXF = 0
        JXS = 0
C
C       The small component functions come first, loop from S to L
C
        DO JC = 2, 1, -1
C
C       Loop over the boson irreps for the right index.
C
          DO JS = 1, NBRP
            SUBJ = JC.EQ.JCOMP.AND.JS.GE.JSS.AND.JS.LE.JSE
            DO 100 JX = 1, NBORB(JS,I,JC)
              JXF = JXF + 1
              IF (SUBJ) THEN
C                Only when we are in a subblock we update this pointer
                 JXS = JXS + 1
              ELSE
C
C             We are not going to run the loops for the left index 
C             and have already updated JXF so we might as well skip 
C             the inner loops
C
                 GOTO 100
              ENDIF
              IXF = 0
              IXS = 0
C
C             Loop over the boson irreps for the left index.
C
              DO IS = 1, NBRP
C
C               Loop from S to L for the left index.
C
                DO IC = 2, 1, -1
                  SUBI = IC.EQ.ICOMP.AND.IS.GE.ISS.AND.IS.LE.ISE
                  DO IX = 1, NBORB(IS,I,IC)
                    IXF = IXF + 1
C
C                   Check the criteria and fill in
C
                    IF (SUBI.AND.SUBJ) THEN
                       IXS = IXS + 1
                       IF (IMODE.EQ.1) THEN
                          FSUB(IXS,JXS,IZ) = FFUL(IXF,JXF,IZ)
                       ELSEIF (IMODE.EQ.2) THEN
                          FFUL(IXF,JXF,IZ) = FSUB(IXS,JXS,IZ)
                       ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
 100        CONTINUE
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck selbos */
      SUBROUTINE SELBOS (NO,INDX,EIG,IBEIG,EVEC,RDUM,IDUM)
C***********************************************************************
C     PURPOSE:
C        Sort eigenvectors from block ordering to energy ordering
C        On entry one should provide an array IBEIG with identification
C        of the spinors.
C
C        Written by L. Visscher
C        April 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG(NO),INDX(NO),IBEIG(NO),IDUM(NO)
      DIMENSION EVEC(NO,NO),RDUM(NO,NO)
C
C     Establish the order
C
      CALL INDEXX (NO,EIG,INDX)
C
C     Sort the eigenvalues and their IDs
C
      DO J = 1, NO
         RDUM(J,1) = EIG(INDX(J))
         IDUM(J) = IBEIG(INDX(J))
      ENDDO
      CALL DCOPY(NO,RDUM,1,EIG,1)
      CALL ICOPY(NO,IDUM,1,IBEIG,1)
C
C     Sort the eigenvectors
C
      DO J = 1, NO
         CALL DCOPY(NO,EVEC(1,INDX(J)),1,RDUM(1,J),1)
      ENDDO
      CALL DCOPY(NO*NO,RDUM,1,EVEC,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck appbos */
      SUBROUTINE APPBOS (IFRP,NO,INDX,IBEIG,EVEC,RDUM)
C***********************************************************************
C     PURPOSE:
C        Identify vectors by their largest boson character
C        This approximate labeling can be used in the CI code.
C
C        Works only for the modified Dirac equation
C        (restricted kinetic balance)
C
C        Written by L. Visscher
C        October 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION IBEIG(NO),INDX(NO)
      DIMENSION EVEC(NO,NO,NZ),RDUM(NO,4)
C
      IF (URKBAL) RETURN
C
      NBRP = 4 / NZ
C
C     Calculate the weight of each boson irrep in each vector.
C     Because we work in the orthonormal basis the total weights
C     add up to one.
C
      CALL DZERO(RDUM,4*NO)
      DO I = 1, NO
         IBEIG(I) = -1
      ENDDO
C
      IX = 1
      DO ISYM = 1, NBRP
        NBO = NBORB(ISYM,IFRP,0)
        IF (NBO.GT.0) THEN
          DO IZ = 1, NZ
             DO I = 1, NO
                RDUM(I,ISYM) = RDUM(I,ISYM) +
     &                         DDOT(NBO,EVEC(IX,I,IZ),1,EVEC(IX,I,IZ),1)
             ENDDO
          ENDDO
          IX = IX + NBO
        ENDIF
      ENDDO
C
      DO ISYM = 1, NBRP
         IREP = JFSYM(ISYM,IFRP) - 1
         NBO = NBORB(ISYM,IFRP,0)
         NSKIP = 0
C
C        Establish the order, indexx will give the vectors with the
C        highest occupation for this boson irreps last
C
         CALL INDEXX (NO,RDUM(1,ISYM),INDX)
C
C        Label the vectors, skip one if already claimed by another
C        irrep. This is not fullproof but will at least satisfy the
C        constraint that every boson irrep gets the correct number
C        of vectors.
C
         DO J = NO, NO-NBO+1, -1
            IBOLD = IBEIG(INDX(J))
            IF (IBOLD.EQ.-1) THEN
               IBEIG(INDX(J)) = IREP
            ELSE
               NSKIP = NSKIP + 1
            ENDIF
         ENDDO
C
         ISKIP = 0
         DO J = NO-NBO,1,-1
            IBOLD = IBEIG(INDX(J))
            IF (IBOLD.EQ.-1.AND.ISKIP.LT.NSKIP) THEN
               ISKIP = ISKIP + 1
               IBEIG(INDX(J)) = IREP
            ENDIF
         ENDDO
      ENDDO
C
c     do i = no/2+1, no
c     write (LUPRI,'(3(A,I4),A,4F8.2)')
c    & 'fermion irrep ',ifrp,' electron function ',i-no/2,' labeled ',
c    &  ibeig(i),' percentages ',(rdum(i,ix),ix=1,nbrp)
c     enddo
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mvogen */
      SUBROUTINE MVOGEN(WORK,LWORK)
C***********************************************************************
C
C     Generate modified virtual orbitals
C
C     Written by T.Saue Sep 18 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "aovec.h"
C
#include "dcbbas.h"
#include "dcborb.h"
C
#include "dcbfir.h"
#include "dcbgen.h"
#include "dcbmvo.h"
#include "blocks.h"
C
      LOGICAL BF1,BF2,BF3
      DIMENSION WORK(*)
C
#include "memint.h"
      CALL TITLER('Modified virtual orbitals','*',125)
C
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBRP,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFOCK,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      NPOS = 0
      IF (PARCAL) THEN
        BF1  = FIRST1
        BF2  = FIRST2
        BF3  = FIRST3
        FIRST1 = .TRUE.
        FIRST2 = .TRUE.
        FIRST3 = .TRUE.
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
      ENDIF
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
      CALL MVOGE1(WORK(KCMO),WORK(KEIG),WORK(KIBRP),
     &            WORK(KFOCK),WORK(KDMAT),
     &            WORK(KPOS),WORK(KFREE),LFREE)
      CALL MEMREL('MVOGEN',WORK,KWORK,KWORK,KFREE,LFREE)
C
      IF (PARCAL) THEN
         FIRST1 = BF1
         FIRST2 = BF2
         FIRST3 = BF3
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mvoge1 */
      SUBROUTINE MVOGE1(CMO,EIG,IBEIG,FOCK,DMAT,NPOS,WORK,LWORK)
C***********************************************************************
C     Generate modified virtual orbitals
C
C     Written by T.Saue June 18 1998
C     Last revision : June 18 1998 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DM1=-1.0D0,DM2=-2.0D0,D2=2.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbmvo.h"
#include "dcbham.h"
      LOGICAL TOBE
      DIMENSION CMO(*),EIG(*),IBEIG(*),FOCK(*),DMAT(*),NPOS(*),WORK(*)
      DIMENSION IC(2),IE(2),ID(2)
      CHARACTER CTMP*1
      CALL QENTER('MVOGEN')
#include "memint.h"
C
C     Initialization
C
      LUBUF = 22
C     full Fock matrix
      IFCMVO = 1
C     only direct (Coulomb) contributions 
      IF(IF2MVO.EQ.1) IFCMVO = 2
C     only exchange contributions
      IF(IF2MVO.EQ.2) IFCMVO = 3
      ISYMVO = 1
      IHRMVO = 1
C
C     Get all coefficients
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('MVOGEN: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','MVOGEN')
         CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14,
     &               WORK,KFREE,LFREE)
      ENDIF
      IF(IONMVO) THEN
C
C       Generate MVOs by the Bauschlicher method (ionic Fock)
C
        CALL DENMVO(DMAT,CMO,FOCK,WORK(KFREE),LFREE)
        CALL TWOFCK(ISYMVO,IHRMVO,IFCMVO,FOCK,DMAT,1,NPOS,
     &              INTMVO,IPRMVO,WORK(KFREE),LFREE)
        INQUIRE(FILE='DFFCK1',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFFCK1','OLD','MVOGE1')
          CALL REAFCK(LUBUF,DMAT,.TRUE.,1)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
          CALL ONEFCK(DMAT,IPRINT,WORK(KFREE),LFREE)
        ENDIF
        CALL DAXPY(N2BBASXQ,D1,DMAT,1,FOCK,1)
      ELSE
C
C       Construct density matrix for chosen vectors
C
        CALL DENMVO(DMAT,CMO,FOCK,WORK(KFREE),LFREE)
        IF(ADDREP) THEN
          CALL DENMVO(FOCK,CMO,DMAT,WORK(KFREE),LFREE)
          CALL GENDEN(DMAT,CMO,1,IPRMVO)
          CALL DAXPY(N2BBASXQ,DM2,FOCK,1,DMAT,1)
        ELSE
          CALL DENMVO(DMAT,CMO,FOCK,WORK(KFREE),LFREE)
          CALL DSCAL(N2BBASXQ,DM1,DMAT,1)
        ENDIF
        IF(IFCMVO.EQ.3) CALL DSCAL(N2BBASXQ,DM1,DMAT,1)
C
C       Construct two-electron Fock matrix
C
        CALL TWOFCK(ISYMVO,IHRMVO,IFCMVO,FOCK,DMAT,1,NPOS,
     &              INTMVO,IPRMVO,WORK(KFREE),LFREE)
      ENDIF
C
C       First transform to virtual MO - basis;
C       this liberates FOCK. Transformed matrix in DMAT
C
      ND = 1
      DO I = 1,NFSYM
      IF (NSSH(I).GT.0 ) THEN
        IC(I) = ICMOQ(I)+NFBAS(I,0)*(NOCC(I)+NPSH(I))+1
        IE(I) = IORB(I) + NOCC(I) + NPSH(I) + 1
        ID(I) = ND
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &           NSSH(I),NSSH(I),
     &           FOCK(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DMAT(ID(I)),NSSH(I),NSSH(I),NZ,IPQTOQ(1,0),
     &           CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRMVO)
        IF (IPRMVO.GE.10 ) THEN
          WRITE(CTMP,'(I1)') I
          CALL HEADER('MVOGEN: Fock matrix in virtual MO-basis'//
     &                 ', corep '//CTMP,-1)
          CALL PRQMAT(DMAT(ID(I)),NSSH(I),NSSH(I),
     &                NSSH(I),NSSH(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
        ND = ND + NSSH(I)*NFBAS(I,0)*NZ
      ENDIF
      ENDDO
C
C     Diagonalize matrix and backtransform vectors
C
      DO I = 1,NFSYM
      IF (NSSH(I).GT.0 ) THEN
        CALL QDIAG(NZ,NSSH(I),DMAT(ID(I)),NSSH(I),NSSH(I),
     &             EIG(IE(I)),1,FOCK(ID(I)),NSSH(I),
     &             NSSH(I),WORK(KFREE),LFREE,IERR)
        CALL BCKTR1(DMAT(ID(I)),NFBAS(I,0),NSSH(I),
     &              FOCK(ID(I)),NSSH(I),NSSH(I),
     &              NSSH(I),NZ,NSSH(I),1,NFBAS(I,0),
     &              CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,
     &              IPRMVO)
        IF(IPRMVO.GE.2) THEN
          WRITE(LUPRI,'(A,A3)')
     &     '* Virtual eigenvalues of ircop ',FREP(I)
          WRITE(LUPRI,'(5F15.8)') (EIG(IE(I)+J),J=0,(NSSH(I)-1))
        ENDIF
C
C       Insert new coefficients
C
        NDIM = NSSH(I)*NFBAS(I,0)
        ICC = IC(I)
        IDD = ID(I)
        DO IZ = 1,NZ
          CALL DCOPY(NDIM,DMAT(IDD),1,CMO(ICC),1)
          ICC = ICC + NCMO(I)
          IDD = IDD + NFBAS(I,0)*NSSH(I)
        ENDDO
C
C       Insert dummy eigenvalues...
C
        CALL DCOPY(NSSH(I),DUMMY,0,EIG(IE(I)),1)
      END IF
      END DO
      CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL QEXIT('MVOGEN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck denmvo */
      SUBROUTINE DENMVO(DMAT,CMO,BUF,WORK,LWORK)
C***********************************************************************
C
C     Generate density matrix for chosen MVO vectors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmvo.h"
      DIMENSION DMAT(*),CMO(*),BUF(*),WORK(*)
      DIMENSION NR(2),IR(2,2),KVEC(2,2),NVEC(2,0:2)
      DIMENSION ICMV(2)
#include "memint.h"
      KFRSAV = KFREE
C
C     Set up index arrays for coefficients contributing
C     =================================================
C
      IC   = 1
      DO I = 1,NFSYM
        NR(I)   = NESH(I)
        IR(1,I) = 0
        IR(2,I) = NESH(I)
      ENDDO
      DO IFRP = 1,NFSYM
        CALL MEMGET('INTE',KVEC(IFRP,1),NR(IFRP),WORK,KFREE,LFREE)
        NVEC(IFRP,0) = 1
        CALL  NUMLST(VECMVO(IFRP),WORK(KVEC(IFRP,1)),
     &               NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &               IFRP,NVEC(IFRP,0),WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KVEC(IFRP,1)),NVEC(IFRP,0),IFRP,
     &              NVEC(IFRP,2),NVEC(IFRP,1))
        NVEC(IFRP,0) = NVEC(IFRP,1) + NVEC(IFRP,2)
        CALL MEMREL('DENMVO',WORK,KFRSAV,KVEC(IFRP,1),KFREE,LFREE)
        CALL MEMGET('INTE',KVEC(IFRP,1),NVEC(IFRP,0),
     &              WORK,KFREE,LFREE)
      ENDDO
C
C     Calculate dimensions of the coefficient array
C
      NCMV = 0
      DO IFRP = 1,NFSYM
        ICMV(IFRP) = NCMV
        NCMV       = NCMV + NFBAS(IFRP,0)*NVEC(IFRP,0)*NZ
      ENDDO
C
C     Select the set that we need
C
      DO IFRP = 1, NFSYM
         CALL SELCFS (CMO(ICMOQ(IFRP)+1),IFRP,BUF(ICMV(IFRP)+1),
     &                NVEC(IFRP,0),
     &                WORK(KVEC(IFRP,1)),NVEC(IFRP,2),NVEC(IFRP,1),
     &                NFBAS(IFRP,0),NORB(IFRP))
      ENDDO
C
C     Get the density matrix
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NVEC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),
     &              BUF(ICMV(I)+1),NFBAS(I,0),NVEC(I,0),
     &              1,NVEC(I,0),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRMVO.GE.4) THEN
        CALL TITLER('Output from DENMVO','*',103)
        DO 20 I = 1,NFSYM
        IF(NVEC(I,0).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Dhftst */
      SUBROUTINE DHFTST(BARNUC,DSCFON,ICHKCNV,BMAT,IPRINT)
C***********************************************************************
C
C     Restart checks etc. for DHF-cycle
C     Based on old code.
C     Written by T.Saue Oct 1 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dftcom.h"
C
      DIMENSION BMAT(MXDIIS*MXDIIS)
      LOGICAL TOBE,OLDFMO,DSCFON,BARNUC
C
      BARNUC  = .TRUE.
      CACC    = '        '
      DSCFON  = .FALSE.
      OLDFMO  = .FALSE.
      ICHKCNV = 0
      KITER   = 1

C
C     1. One-electron system
C     ======================
C
      IF(ONESYS) THEN
        TRIVEC = .FALSE.
        TRIFCK = .FALSE.
        NITER  = 1
        RETURN
      ENDIF
C
C
C     2. Many-electron system
C     =======================
C
C
C     ***************************************
C     **** Looking for old wave function ****
C     ***************************************
C
C     * TRIFCK: Looking for old 2-electron Fock matrix
C     ------------------------------------------------
C
      IF(TRIFCK.AND.NASHT.GT.0) THEN
C         jth April 98
         WRITE(LUPRI,'(/A/A/)') '*** WARNING *** '//
     &      'Trial 2-electron Fock matrix not possible for open shell.',
     &      'Looking for trial vectors instead'
         TRIFCK = .FALSE.
         TRIVEC = .TRUE.
      END IF
      IF(TRIFCK) THEN
        INQUIRE(FILE='DFFCK2',EXIST=TOBE)
        IF(TOBE) THEN
          DHF_INTTYP = 'Old FOCK2'
          BARNUC = .FALSE.
        ELSE
          WRITE(LUPRI,'(/A/)')
     &    '*** WARNING *** No trial 2-electron Fock matrix found.'//
     &    'Looking for trial vectors instead'
          TRIFCK = .FALSE.
          TRIVEC = .TRUE.
        ENDIF
      ENDIF
C
C     * TRIVEC: Looking for old coefficients
C     --------------------------------------
C
      IF(TRIVEC) THEN
        INQUIRE(FILE='DFCOEF',EXIST=TOBE)
        IF(TOBE) THEN
          BARNUC = .FALSE.
        ELSE
          WRITE(LUPRI,'(/A)')
     &    '*** WARNING *** No trial vectors found.'//
     &    'Using bare nucleus approximation instead.'
          IF (BNCRON) WRITE(LUPRI,'(A/)')
     &    '                Improved by a crude '//
     &    'estimate of the two-electron repulsion.'
          TRIVEC = .FALSE.
C     Do not use density-based atomic size adjustments in this case:
          ATSIZE = .FALSE.
        ENDIF
      ENDIF
C
C     ***********************************
C     **** Looking for cycle history ****
C     ***********************************
C
      INQUIRE(FILE='DFCYCL',EXIST=TOBE)
      TOBE = TOBE.AND.(TRIVEC.OR.TRIFCK)
      IF(TOBE) THEN
        OPEN(LUCYCL,FILE ='DFCYCL',STATUS='OLD',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
        REWIND LUCYCL
 10     CONTINUE
        READ(LUCYCL,1020,END=20) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &              CACC,CPUTID,DHF_INTTYP,DAYTID
        GOTO 10
 20     CONTINUE
#if defined (SYS_IBM) || defined (SYS_CRAY) || defined (VAR_MFDS)
C        backspace on multifile systems for correct positioning
C        at the end of the file
        BACKSPACE LUCYCL
#endif
        WRITE(LUPRI,'(A)') '*** WARNING *** : SCFCYC restart info:'
        WRITE(LUPRI,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &            CACC,CPUTID,DHF_INTTYP,DAYTID
        ILLINT = 0
        IF(DHF_INTTYP(1:2).EQ.'LL') ILLINT = 1
        ISLINT = 0
        IF(DHF_INTTYP(4:5).EQ.'SL') ISLINT = 1
        ISSINT = 0
        IF(DHF_INTTYP(7:8).EQ.'SS') ISSINT = 1
        INTBUF = ILLINT + 2*ISLINT + 4*ISSINT
        NITER = NITER + 1
        CALL INTCON(INTFLG,INTBUF,INTDEF,
     &            ERGVAL,CNVINT,NITER,ITRINT,DHF_INTTYP)
C
C       * Check on convergence acceleration
C       -----------------------------------
C
        INQUIRE(FILE='DFFOCK',EXIST=OLDFMO)
        OLDFMO = OLDFMO.AND.(DODIIS.OR.DODAMP)
        IF(DODIIS) THEN
          DIISON = OLDFMO.AND.INTFLG.EQ.INTBUF.AND.NITER.GT.1
          INQUIRE(FILE='DFDIIS',EXIST=TOBE)
          DIISON = DIISON.AND.TOBE
          INQUIRE(FILE='DFCMOS',EXIST=TOBE)
          DIISON = DIISON.AND.TOBE
          INQUIRE(FILE='DFEVEC',EXIST=TOBE)
          DIISON = DIISON.AND.TOBE
          IF(DIISON) THEN
            CALL OPNFIL(LUDIIS,'DFDIIS','OLD','SCFINT')
            READ(LUDIIS) ITDIIS,MDIIS2,NELMBM,BMAT
            IF(MDIIS2.EQ.MXDIIS) THEN
              OLDFMO = .TRUE.
              WRITE(LUPRI,'(A)') '*** WARNING *** : DIIS RESTART'
              WRITE(LUPRI,'(3X,A,I5)')
     &        '* Previous DIIS iteration        :', ITDIIS,
     &        '* Number of elements in B-matrix :',NELMBM
            ELSE
              ITDIIS = 0
            ENDIF
          ELSE
            ITDIIS = 0
          ENDIF
        ENDIF
C
C       Check for differential density matrix
C
        IF(DODSCF.AND.(INTBUF.EQ.INTFLG)) THEN
          INQUIRE(FILE='DFFCK2',EXIST=DSCFON)
          INQUIRE(FILE='DFDENS',EXIST=TOBE)
          DSCFON = DSCFON.AND.TOBE
          IF(DSCFON) THEN
            WRITE(LUPRI,'(A)')
     &      '*** WARNING *** : Differential density restart'
          ENDIF
        ENDIF
      ELSE
        NITER  = 1
        OPEN(LUCYCL,FILE ='DFCYCL',STATUS='UNKNOWN',
     &     ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
        IF(TRIVEC) CALL INTCON(INTFLG,INTBUF,INTDEF,
     &                  DUMMY,CNVINT,NITER,ITRINT,DHF_INTTYP)
      ENDIF
C
C     If (oldfmo) we can check all convergence criteria
C
      IF (OLDFMO) ICHKCNV = 1
      RETURN
 1010 FORMAT('It. ',I3,1P,G23.13,3D10.2,3X,A8,3X,A12,3X,A9,3X,A10)
 1020 FORMAT(4X,    I3,   G23.13,3D10.2,3X,A8,3X,A12,3X,A9,3X,A10)
C     ... 1020 is used for reading, thus no '1P' !!! /hjaaj
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck fndocc */
      SUBROUTINE FNDOCC(EIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Find DHF occupation based on eigenvalues
C     using the build-up principle.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
C
      DIMENSION EIG(*),WORK(*)
      DIMENSION IEIGMIN(2)
C
      DO I = 1,NFSYM
         IEIGMIN(I) = 1+IORB(I) + NPSH(I)
         NISH(I) = 0
         NASH(I) = 0
         NACSH(I,1) = 0
      END DO
C
      NELEC2GO = NELECT_DHF
      IF ((NELECT_DHF+1)/2 .GT. NESHT) THEN
         WRITE(LUPRI,'(I3,A,I3,A)')
     &      NELECT_DHF,' electrons but only ',NESHT,' orbitals.'
         CALL QUIT('ERROR: Not enough eigenvalues for FNDOCC.')
      END IF
C
C     100 loop is over electron pairs (closed shell orbitals)
C
 100  CONTINUE
C        Find the irrep that has the lowest eigenvalue
         IEIGSYM = 0
 110     IEIGSYM = IEIGSYM + 1
         IF (NISH(IEIGSYM) .GE. NESH(IEIGSYM)) GOTO 110
         J = IEIGSYM
         DO I = J+1,NFSYM
            IF (EIG(IEIGMIN(I)) .LT. EIG(IEIGMIN(IEIGSYM)))
     &         IEIGSYM = I
         END DO
C        We found the lowest yet unoccupied MO, occupy it
         NISH(IEIGSYM) = NISH(IEIGSYM) + 1
         IEIGMIN(IEIGSYM) = IEIGMIN(IEIGSYM) + 1
         NELEC2GO = NELEC2GO - 2
C        More orbitals to be filled ?
      IF (NELEC2GO .GT. 0) GOTO 100
C
C     Same procedure for the open shell (works only for one e in one open shell)
C
      IF (NAELEC_DHF .GT. 0) THEN
         IEIGSYM = 0
 200     IEIGSYM = IEIGSYM + 1
            IF (NISH(IEIGSYM) .GE. NESH(IEIGSYM)) GOTO 200
         J = IEIGSYM
         DO I = J+1,NFSYM
            IF (EIG(IEIGMIN(I)) .LT. EIG(IEIGMIN(IEIGSYM)))
     &         IEIGSYM = I
         END DO
         NASH(IEIGSYM) = 1
         NACSH(IEIGSYM,1) = 1
      END IF
C
      NASHMFT = 0
      DO I = 1, NFSYM
         NOCC(I) = NISH(I) + NASH(I)
         NISHMF(I) = NISH(I)
C
C        LV : CAUTION, POSSIBLE BUG
C
C        After the MCSCF merge we now have two arrays that seem to be
C        used interchangeably at various places, I just copy the info
C        to both so that the auto-occupation feature works for DHF but
C        this needs to be checked by the MCSCF bunch !!!
C
         NISH_DHF(I) = NISH(I)
         NASH_DHF(I) = NASH(I)
      END DO
C     Subtract frozen orbitals (if any) from orbital classes
      CALL SUBFRZ(WORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dfdiag */
      SUBROUTINE DFDIAG(FMO,EIG,IBEIG,FAO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Diagonalize Dirac-Fock matrix
C     On input:
C       FMO - MO-coefficients from previous iteration
C     On output:
C       FAO - MO-coefficients from current iteration
C       EIG - eigenvalues
C
C     Taken out of DFSOLV and modified by L.Visscher Jan 10 2000
C     Last revision: Jan 10 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      PARAMETER (D0 = 0.00D00, D1 = 1.00D00)
      LOGICAL SPECIAL
      DIMENSION FMO(*),EIG(*),IBEIG(*),FAO(*),WORK(*)
C
C     Check for special cases :
C     - Spinfree or Levy-Leblond; diagonalize boson irreps (subblocks)
C     - Linear symmetry; diagonalize Omega irreps (subblocks)
C     - ZORA; first diagonalize SS-block, than form and diagonalize
C       modified large component block
C
      SPECIAL = (ZORA.OR.SUB_BL)
      CALL DZERO(FAO,N2BBASXQ)
      NBRP = 4 / NZ
C
      DO 30 I = 1,NFSYM
C
        IF(NTMO(I).EQ.0) GOTO 30
C
C       The next array is used in subblock calcs, initialize to zero.
C
        DO J = 1, NORB(I)
          IBEIG(IORB(I)+J) = 0
        ENDDO
C
C       Diagonalize Fock-matrix in MO-basis
C         - store coefficients in FAO
C       ===================================
C
C       1 : Normal diagonalization, use matrix as supplied by dfsolv
C
        IF (.NOT.SPECIAL) THEN
          IF(DOJACO) THEN
            CALL RSJACO(NTMO(I),NTMO(I),NTMO(I),FMO(I2TMOT(I)+1),
     &                  EIG(IORB(I)+1),1,1,0,FAO(I2TMOT(I)+1),
     &                  WORK(KFREE),LFREE)
          ELSE
            CALL QDIAG(NZ,NTMO(I),FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),
     &                 EIG(IORB(I)+1),1,FAO(I2TMOT(I)+1),NTMO(I),
     &                 NTMO(I),WORK(KFREE),LFREE,IERR)
          ENDIF
C
C
C       2 : Subblock diagonalization (e.g. Levy-Leblond, spinfree, linear symmetry)
C           use matrix as supplied by dfsolv but diagonalize blockwise
C
        ELSEIF (.NOT.ZORA) THEN
C
C         Spinfree calculations can use the boson blocking of the
C         Fock matrix. For quaternion groups we use the fact that
C         the matrix becomes real. For linear symmetry we also have
C         real matrices, but a blocking given by another array.
C         Initialize FAO to zero.
C         --------------------------------------------------------
C
          NORB2 = NORB(I)*NORB(I)
          CALL DZERO(FAO(1+I2ORBT(I)),NORB2*NZ)
C
          IBO = 0
          DO ISUB = 1, N_SUB_BL(I)
           NBO = NORB_SUB(ISUB,I,0)
           IF(NBO.GT.0) THEN
            IMAT = I2ORBT(I) + (NORB(I)+1)*IBO + 1
            IEIG = IORB(I) + IBO + 1
C
C           Diagonalize the isub'th block of parity i
C
            IF (DOJACO) THEN
              CALL RSJACO(NORB(I),NBO,NBO,FMO(IMAT),
     &                    EIG(IEIG),1,1,0,FAO(IMAT),
     &                    WORK(KFREE),LFREE)
            ELSE
              CALL QDIAG(1,NBO,FMO(IMAT),NORB(I),NORB(I),
     &                   EIG(IEIG),1,FAO(IMAT),NORB(I),NORB(I),
     &                   WORK(KFREE),LFREE,IERR)
            ENDIF
            IBO = IBO + NORB_SUB(ISUB,I,0)
           ENDIF
          ENDDO
C
C         Sort the eigenvalues according to energy : keep track of
C         the block identification. (Use FMO as scratch space for the
C         vector sort)
C
          CALL INIBOS(IBEIG(IORB(I)+1),I)
          CALL MEMGET('INTE',KINDX,NORB(I),WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KIDUM,NORB(I),WORK,KFREE,LFREE)
          CALL SELBOS (NORB(I),WORK(KINDX),EIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1),FAO(I2ORBT(I)+1),
     &                 FMO(I2ORBT(I)+1),WORK(KIDUM))
          CALL MEMREL('DFDIAG',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       3 : ZORA; diagonalize modified LL block
C
        ELSEIF (ZORA.AND.(.NOT.SUB_BL)) THEN
C
C         Set dimensions
C
          NL = NESH(I)
          NS = NPSH(I)
          NT = NORB(I)
          NLL = NL * NL
          NLS = NL * NS
          NSS = NS * NS
          NTT = NT * NT
          IEIGP = IORB(I) + 1
          IEIGE = IEIGP + NS
C
C         Allocate memory for the ZORA diagonalization. We need to
C         have the LL, LS, SL and SS blocks in separate blocks.
C
          CALL MEMGET('REAL',KFLL,NZ*NLL,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KFLS,NZ*NLS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KFSS,NZ*NSS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KVLL,NZ*NLL,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KVSL,NZ*NLS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KVSS,NZ*NSS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KWLL,NZ*NLL,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KWLS,NZ*NLS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KWSL,NZ*NLS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KWSS,NZ*NSS,WORK,KFREE,LFREE)
C
C         Extract LL block of Fock matrix
          CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                 0,0,1,1,1,NZ,
     &                 WORK(KFLL),NL,NL,NZ,
     &                 FMO(I2ORBT(I)+1),NT,NT,NZ)
C         Extract LS block of Fock matrix
          CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                 0,0,1,2,1,NZ,
     &                 WORK(KFLS),NL,NS,NZ,
     &                 FMO(I2ORBT(I)+1),NT,NT,NZ)
C         Extract SS block of Fock matrix
          CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                 0,0,2,2,1,NZ,
     &                 WORK(KFSS),NS,NS,NZ,
     &                 FMO(I2ORBT(I)+1),NT,NT,NZ)
C
C         Do the diagonalization using the ZORA approximation
C
          CALL ZORADG(NZ,NL,NS,WORK(KFLL),WORK(KFLS),WORK(KFSS),
     &         WORK(KVLL),WORK(KVSL),WORK(KVSS),
     &         WORK(KWLL),WORK(KWLS),WORK(KWSL),WORK(KWSS),
     &         EIG(IEIGE),EIG(IEIGP),
     &         WORK,KFREE,LFREE)
C
C         Put the MO vectors in the FAO matrix, actually only the
C         electron solutions are calculated (and stored in FAO(LL),
C         FAO(SL)) but we put the eigenvectors of the equation
C         (VSS-2C^2) CS = E CS in the SS block and zero the LS block.
C
C         Note that we cannot use EXTRSB as this routine will sort the
C         eigenctors assuming a boson block structure which does not
C         apply here. The point is that we now have the solutions
C         on the right index (ordered as positron, electron) while
C         the basis is on the left ordered via boson irreps first.
C         In the spinfree ZORA case we do use EXTRSB because the
C         boson subblocking is now also (enforced) present for the
C         solutions.
C
          CALL EXTRSBR (2,I,0,0,1,1,1,NZ,WORK(KVLL),NL,NL,NZ,
     &                 FAO(I2ORBT(I)+1),NT,NT)
          CALL EXTRSBR (2,I,0,0,2,1,1,NZ,WORK(KVSL),NS,NL,NZ,
     &                 FAO(I2ORBT(I)+1),NT,NT)
          CALL EXTRSBR (2,I,0,0,2,2,1,NZ,WORK(KVSS),NS,NS,NZ,
     &                 FAO(I2ORBT(I)+1),NT,NT)
          CALL DZERO (WORK(KWLS),NZ*NLS)
          CALL EXTRSBR (2,I,0,0,1,2,1,NZ,WORK(KWLS),NL,NS,NZ,
     &                 FAO(I2ORBT(I)+1),NT,NT)
C
C         Done; release the extra memory used in ZORA
C
          CALL MEMREL('DFDIAG',WORK,KFLL,KFLL,KFREE,LFREE)
C
C       4 : Spinfree ZORA : Make modified matrix but do it block by block
C
        ELSEIF (ZORA.AND.SPINFR) THEN
C
C           The total dimension of the fermion irrep can be set now
C
            NT = NORB(I)
            NTT = NT * NT
C
            IEIGP = IORB(I) + 1
            DO 21 ISUB = 1, N_SUB_BL(I)
C
C             Set dimensions within the boson subblock
C
              NBO = NORB_SUB(ISUB,I,0)
              IF (NBO.EQ.0) GOTO 21
              NBO2 = NBO * NBO
              NL = NORB_SUB(ISUB,I,1)
              NS = NORB_SUB(ISUB,I,2)
              NLL = NL * NL
              NLS = NL * NS
              NSS = NS * NS
              NZS = 1
              IEIGE = IEIGP + NS
C
C             Allocate memory for the ZORA diagonalization. We need to
C             have the LL, LS, SL and SS blocks in separate blocks.
C
              CALL MEMGET('REAL',KFLL,NZ*NLL,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KFLS,NZ*NLS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KFSS,NZ*NSS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KVLL,NZ*NLL,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KVSL,NZ*NLS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KVSS,NZ*NSS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KWLL,NZ*NLL,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KWLS,NZ*NLS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KWSL,NZ*NLS,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KWSS,NZ*NSS,WORK,KFREE,LFREE)
C
C             Extract LL block of Fock matrix
              CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,1,1,1,NZS,
     &                     WORK(KFLL),NL,NL,NZS,
     &                     FMO(I2ORBT(I)+1),NT,NT,NZ)
C             Extract LS block of Fock matrix
              CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,1,2,1,NZS,
     &                     WORK(KFLS),NL,NS,NZS,
     &                     FMO(I2ORBT(I)+1),NT,NT,NZ)
C             Extract SS block of Fock matrix
              CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,2,2,1,NZS,
     &                     WORK(KFSS),NS,NS,NZS,
     &                     FMO(I2ORBT(I)+1),NT,NT,NZ)
C
C             Do the diagonalization using the ZORA approximation
C
              CALL ZORADG(NZ,NL,NS,WORK(KFLL),WORK(KFLS),WORK(KFSS),
     &             WORK(KVLL),WORK(KVSL),WORK(KVSS),
     &             WORK(KWLL),WORK(KWLS),WORK(KWSL),WORK(KWSS),
     &             EIG(IEIGE),EIG(IEIGP),
     &             WORK,KFREE,LFREE)
C
C             Put the MO vectors in the FAO matrix, actually only the
C             electron solutions are calculated (and stored in FAO(LL),
C             FAO(SL)) but we put the eigenvectors of the equation
C             (VSS-2C^2) CS = E CS in the SS block and zero the LS block.
C
C             Insert LL block of Fock matrix
              CALL EXTRSB (2,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,1,1,1,NZS,
     &                     WORK(KVLL),NL,NL,NZS,
     &                     FAO(I2ORBT(I)+1),NT,NT,NZ)
C             Insert SL block of Fock matrix
              CALL EXTRSB (2,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,2,1,1,NZS,
     &                     WORK(KVSL),NS,NL,NZS,
     &                     FAO(I2ORBT(I)+1),NT,NT,NZ)
C             Insert SS block of Fock matrix
              CALL EXTRSB (2,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,2,2,1,NZS,
     &                     WORK(KVSS),NS,NS,NZS,
     &                     FAO(I2ORBT(I)+1),NT,NT,NZ)
C             Insert LS block of Fock matrix
              CALL DZERO (WORK(KWLS),NZ*NLS)
              CALL EXTRSB (2,I,NBORB,4,NBRP,NBORB,4,NBRP,
     &                     ISUB,ISUB,1,2,1,NZS,
     &                     WORK(KWLS),NL,NS,NZS,
     &                     FAO(I2ORBT(I)+1),NT,NT,NZ)
C
C             Done; release the extra memory used in ZORA
C
              CALL MEMREL('DFDIAG',WORK,KFLL,KFLL,KFREE,LFREE)
C
C             Update pointers for the eigenvalues
C
              IEIGP = IEIGP + NBO
  21        CONTINUE
C
C           Sort the eigenvalues according to energy : keep track of
C           the boson symmetry. (Use FMO as scratch space for the
C           vector sort)
C
            CALL INIBOS(IBEIG(IORB(I)+1),I)
            CALL MEMGET('INTE',KINDX,NORB(I),WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KIDUM,NORB(I),WORK,KFREE,LFREE)
            CALL SELBOS (NORB(I),WORK(KINDX),EIG(IORB(I)+1),
     &                   IBEIG(IORB(I)+1),FAO(I2ORBT(I)+1),
     &                   FMO(I2ORBT(I)+1),WORK(KIDUM))
            CALL MEMREL('DFDIAG',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       End of different diagonalization cases
C
        ENDIF
C
C       Check if any positronic eigenvalues are below -2mc^2
C
        IF(IPRSCF.GE.2) CALL EIGCHK(EIG(IORB(I)+1),NPSH(I))
C
C       Do an approximate identification of the spinors by their
C       boson character for use in the CI.
C
        IF (.NOT.SPINFR .AND. .NOT.LEVYLE .AND. .NOT.SUB_BL) THEN
           CALL MEMGET('INTE',KINDX,NORB(I),WORK,KFREE,LFREE)
           CALL MEMGET('REAL',KRDUM,4*NORB(I),WORK,KFREE,LFREE)
           CALL APPBOS (I,NTMO(I),WORK(KINDX),
     &                  IBEIG(IORB(I)+1),FAO(I2TMOT(I)+1),
     &                  WORK(KRDUM))
           CALL MEMREL('DFDIAG',WORK,KINDX,KINDX,KFREE,LFREE)
        END IF
C
 30   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck zoradg */
      SUBROUTINE ZORADG(NZS,NL,NS,FLL,FLS,FSS,VLL,VSL,VSS,
     &                  WLL,WLS,WSL,WSS,EIGE,EIGP,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     Diagonalize the Dirac matrix using the ZORA approximation :
C     first fold in the small components and build an effective large
C     component equation.
C
C     Taken out of DFDIAG and modified by L.Visscher Aug 28 2000
C     Last revision: Aug 28 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      PARAMETER (D0 = 0.00D00, D1 = 1.00D00)
      DIMENSION EIGP(NS),EIGE(NL),WORK(*)
      DIMENSION FLL(NL,NL,NZS),FSS(NS,NS,NZS),FLS(NL,NS,NZS)
      DIMENSION VLL(NL,NL,NZS),VSS(NS,NS,NZS),VSL(NS,NL,NZS)
      DIMENSION WLL(NL,NL,NZS),WSS(NS,NS,NZS)
      DIMENSION WLS(NL,NS,NZS),WSL(NS,NL,NZS)
C
C     Diagonalize SS block of the matrix
C
      CALL QDIAG(NZS,NS,FSS,NS,NS,
     &           EIGP,1,VSS,NS,NS,
     &           WORK(KFREE),LFREE,IERR)
C
C     Form U * (2c^2-V)**-1/2
C
      CALL DZERO (WSS,NZS*NS*NS)
      DO IZ = 1, NZS
         DO J = 1, NS
            EIGINV = D1 / SQRT(-EIGP(J))
            CALL DAXPY(NS,EIGINV,VSS(1,J,IZ),1,WSS(1,J,IZ),1)
         ENDDO
      ENDDO
C
C     Form c(sigma.p) * U * (2c^2-V)**-1/2
C
      CALL QGEMM(NL,NS,NS,D1,
     &     'N','N',IPQTOQ(1,0),FLS,NL,NS,NZS,
     &     'N','N',IPQTOQ(1,0),WSS,NS,NS,NZS,
     &      D0,IPQTOQ(1,0),WLS,NL,NS,NZS)
C
C     Form c(sigma.p) * U * (2c^2-V)**-1 * U+ * c(sigma.p)
C     This is the ZORA modified K.E. matrix
C
      CALL QGEMM(NL,NL,NS,D1,
     &     'N','N',IPQTOQ(1,0),WLS,NL,NS,NZS,
     &     'H','N',IPQTOQ(1,0),WLS,NL,NS,NZS,
     &      D0,IPQTOQ(1,0),WLL,NL,NL,NZS)
C
C     Add the ZORA modified K.E. matrix to the LL-block
C
      CALL DAXPY (NZS*NL*NL,D1,WLL,1,FLL,1)
C
C     Diagonalize the ZORA matrix
C
      CALL QDIAG(NZS,NL,FLL,NL,NL,
     &           EIGE,1,VLL,NL,NL,
     &           WORK(KFREE),LFREE,IERR)
C
C
C     Form C[S] = U * (2c^2-V)**-1 * U+ * c(sigma.p) * C[L]
C               = {U * (2c^2-V)**-1/2} *
C                 {(2c^2-V)**-1/2 * U+ * c(sigma.p)} * C[L]
C
C     The coefficients for the small components.
C
      CALL QGEMM(NS,NL,NS,D1,
     &     'N','N',IPQTOQ(1,0),WSS,NS,NS,NZS,
     &     'H','N',IPQTOQ(1,0),WLS,NL,NS,NZS,
     &      D0,IPQTOQ(1,0),WSL,NS,NL,NZS)
C
      CALL QGEMM(NS,NL,NL,D1,
     &     'N','N',IPQTOQ(1,0),WSL,NS,NL,NZS,
     &     'N','N',IPQTOQ(1,0),VLL,NL,NL,NZS,
     &      D0,IPQTOQ(1,0),VSL,NS,NL,NZS)
C
C     Determine scaling factors
C
      CALL QGEMM(NL,NL,NS,D1,
     &     'H','N',IPQTOQ(1,0),VSL,NS,NL,NZS,
     &     'N','N',IPQTOQ(1,0),VSL,NS,NL,NZS,
     &      D0,IPQTOQ(1,0),WLL,NL,NL,NZS)
C
C     Normalize the vectors for ZORA-4
C
      IF (ZORA4)  THEN
         DO J = 1, NL
            FACTOR = D1 / SQRT(D1+WLL(J,J,1))
            DO IZ = 1, NZS
               CALL DSCAL (NL,FACTOR,VLL(1,J,IZ),1)
               CALL DSCAL (NS,FACTOR,VSL(1,J,IZ),1)
            ENDDO
         ENDDO
      ENDIF
C
C     Calculate the scaled ZORA eigenvalues
C
      IF (ZORASC)  THEN
         DO J = 1, NL
            FACTOR = D1 / (D1 + WLL(J,J,1))
            EIGE(J) = EIGE(J) * FACTOR
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inibos */
      SUBROUTINE INIBOS(IBOS,IFRP)
C***********************************************************************
C
C     Fill an array with boson irrep identification
C
C     A Saue cut and paste job, thanks to Luuk ...April 2003
C     Generalized to general subblocks
C
C***********************************************************************
#include "dcborb.h"
      DIMENSION IBOS(*)
C
      NSUB = N_SUB_BL(IFRP)
      IX = 1
      DO ISUB = 1, NSUB
        ID = ID_SUB_BL(ISUB,IFRP)
        DO IC = 2, 1, -1
          NBO = NORB_SUB(ISUB,IFRP,IC)
          IF(NBO.GT.0) THEN
            CALL ICOPY(NBO,ID,0,IBOS(IX),1)
            IX = IX + NBO
          ENDIF
        ENDDO
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck detocc */
      SUBROUTINE DETOCC(EIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Handle automatic occupation
C     A cut-and-paste job by T. Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
      LOGICAL TEST
      DIMENSION EIG(*),WORK(*)
      DIMENSION NISHOLD(2),NASHOLD(2)
C
C     Save old occupation
C
      CALL ICOPY(2,NISH,1,NISHOLD,1)
      CALL ICOPY(2,NASH,1,NASHOLD,1)
C
C     Determine new occupation
C
      CALL FNDOCC(EIG,WORK,KFREE,LFREE)
C
      TEST = .FALSE.
      DO IFRP = 1,NFSYM
        TEST = TEST.OR.(NISH(IFRP).NE.NISHOLD(IFRP))
        TEST = TEST.OR.(NASH(IFRP).NE.NASHOLD(IFRP))
      END DO
      IF (TEST.OR.INIOCC) THEN
C
C        AUTOCC has changed occupation...
C        Update common block info about N2ISHT,...
C
         CALL SETDC2(0)
         CALL SETDHF(0)
         IF(INIOCC) THEN
           INIOCC = .FALSE.
           WRITE(LUPRI,'(/A/)') '* AUTOCC( 0) : Initial occupation:'
         ELSE
           NEWOCC = NEWOCC + 1
           WRITE(LUPRI,'(/A,I2,A/)')
     &            '* AUTOCC(',NEWOCC,') : New occupation:'
         END IF
C
         CALL PROCC()
C
         IF(NEWOCC .GT. 10) THEN
           WRITE(LUPRI,'(/A,/A,I2,A)') 'Problems with AUTOCC:',
     &            'DHF occupation has changed ',NEWOCC, ' times.'
           CALL QUIT('DHF occupation has changed too many times.')
         END IF
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Subfrz */
      SUBROUTINE SUBFRZ(WORK,KFREE,LFREE)
C***********************************************************************
C
C     For each orbital class determine the number of frozen orbitals
C
C     Written by T. Saue Apr 16 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION WORK(*)
C
      IF(NFROT.GT.0) CALL OPNFIL(LUTMAT,'DFFROZ','OLD','SUBFRZ')
      DO IFRP = 1,NFSYM
        NISHMF(IFRP) = NISH(IFRP)
        DO IOPEN = 1,NOPEN
          NACSHMF(IFRP,IOPEN) = NACSH(IFRP,IOPEN)
        END DO
        IF(NFRO(IFRP).GT.0) THEN
          CALL MEMGET('INTE',KIND,NFRO(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KIBF,NORB(IFRP),WORK,KFREE,LFREE)
          CALL ORBARR(IFRP,WORK(KIBF))
          CALL READI(LUTMAT,NFRO(IFRP),WORK(KIND))
          READ(LUTMAT)
          CALL SUBFR1(IFRP,WORK(KIND),WORK(KIBF))
          CALL MEMREL('SUBFRZ.froz',WORK,KIND,KIND,KFREE,LFREE)
        ENDIF
      ENDDO
      IF(NFROT.GT.0) CLOSE(LUTMAT,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Subfr1 */
      SUBROUTINE SUBFR1(IFRP,IND,IBUF)
C***********************************************************************
C
C     For each orbital class determine the number of frozen orbitals
C
C     Written by T. Saue Apr 16 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION IND(*),IBUF(*)
C
      DO I = 1,NFRO(IFRP)
        IP = IBUF(IND(I))
        IF(IP.EQ.0) THEN
          NISHMF(IFRP) = NISHMF(IFRP) - 1
        ELSEIF(IP.GT.0) THEN
          NACSHMF(IFRP,IP) = NACSHMF(IFRP,IP) - 1
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck orbstr */
      SUBROUTINE ORBARR(IFRP,IBUF)
C***********************************************************************
C
C     Make an array to distinguish orbital classes
C     
C     Written by T. Saue Apr 16 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION IBUF(*)
C.....Negative-energy orbitals:
      II = 1
      CALL ICOPY(NPSH(IFRP),-2,0,IBUF(II),1)
C.....Inactive orbitals:
      II = II + NPSH(IFRP)
      CALL ICOPY(NISH(IFRP),0,0,IBUF(II),1)
C.....Active orbitals:
      II = II + NISH(IFRP)
      DO IOPEN = 1,NOPEN
        CALL ICOPY(NACSH(IFRP,IOPEN),IOPEN,0,IBUF(II),1)
        II = II + NACSH(IFRP,IOPEN)
      ENDDO
C.....Secondary positive-energy orbitals:
      CALL ICOPY(NSSH(IFRP),-1,0,IBUF(II),1)
C
      RETURN
      END



