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
C VERSION : $Revision: 1.41.2.2 $
C DATE    : $Date: 2004/09/01 08:16:13 $
C FILE    : pamset.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxop */
      SUBROUTINE SETXOP(INDSTR,SKIPEE,SKIPEP,GAS,IOPSY,NOROT,
     &                  NXOPE,NHOPE,KXOPE, NXOPP,KXOPP,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Set orbital rotation index vectors
C
C     Orbitals are classified as:
C       inactive  : doubly occupied in all determinants
C       active    : singly occupied in some determinants
C       secondary : unoccupied in all determinants
C       positronic: positronic orbitals
C
C     Ordering is:
C       e-e rotations: (skipped if SKIPEE = T)
C          IA : inactive-active   
C          IS : inactive-secondary
C          AS : active-secondary
C          AB : active-active
C       e-p rotations: (skipped if SKIPEP = T)
C          IP : inactive-positronic
C          AP : active-positronic
C         (SP : secondary-positronic   -  only for QED=T)
C
C     NXOPE is the number of non-redundant e-e rotations.
C     NHOPE is the total number of rotations.
C
C     Written by  Hans Joergen Aa. Jensen and T.Saue 27-Jan-1995
C     Last revision: July 21 1997 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from common blocks:
C  DGROUP: NFSYM
C  DCBHAM: QED
C  DCBORB: IORB,NPSH,NFRO,NISH,NASH
C
#include "dgroup.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      LOGICAL SKIPEP,SKIPEE,GAS
      CHARACTER INDSTR(3,2)*72
      DIMENSION WORK(*),NSTR(0:2,3,2)
      DIMENSION NOROT(*)
C
      CALL QENTER('SETXOP')
      KXOPE = KFREE
      KXOPP = KFREE
      NXOPE = 0
      NXOPP = 0
      NHOPE = 0
C
C     Allocate buffer areas
C     =====================
C
      NBUF = N2BBASX
C
C     ***********************************
C     *** Electron-electron rotations ***
C     ***********************************
C      
 10   CONTINUE
      IF(SKIPEE) GOTO 20
      CALL MEMGET('INTE',KBUF ,NBUF          ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,MXFBAS*3*NFSYM,WORK,KFREE,LFREE)
C     ... Unpack strings
      CALL UNPSTR(INDSTR,3,WORK(KIBUF),MXFBAS,NSTR,
     &            WORK,KFREE,LFREE,IPRINT)
      CALL SETXEE(WORK(KBUF),IOPSY,NXOPE,NHOPE,WORK(KIBUF),
     &     NSTR,GAS,IPRINT,NOROT)
      IF(NHOPE.GT.0) THEN
        CALL MEMREL('SETXOP.e-e',WORK,1,KBUF,KFREE,LFREE)
        CALL MEMGET('INTE',KXOPE,2*NHOPE,WORK,KFREE,LFREE)
        NBUF = N2BBASX - 2*(NHOPE+1)
      ENDIF
C
C     ***********************************
C     *** Electron-positron rotations ***
C     ***********************************
C      
 20   CONTINUE
      IF(SKIPEP) GOTO 30
      CALL MEMGET('INTE',KBUF ,NBUF          ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,MXFBAS*3*NFSYM,WORK,KFREE,LFREE)
C     ... Unpack strings
      CALL UNPSTR(INDSTR,3,WORK(KIBUF),MXFBAS,NSTR,
     &            WORK,KFREE,LFREE,0)
      CALL SETXEP(WORK(KBUF),IOPSY,NXOPP,QED,WORK(KIBUF),NSTR,IPRINT,
     &     NOROT)
      CALL MEMREL('SETXOP.e-p',WORK,1,KBUF,KFREE,LFREE)
      CALL MEMGET('INTE',KXOPP,2*NXOPP,WORK,KFREE,LFREE)
C
 30   CONTINUE
C
      CALL QEXIT('SETXOP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setxee */
      SUBROUTINE SETXEE(JXOPE,IOPSY,NXOPE,NHOPE,IND,NIND,GAS,IPRINT,
     &                  NOROT)
C   
C     Generate electron-electron rotations
C
C     Written by T.Saue July 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dcbibn.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      DIMENSION JXOPE(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      DIMENSION NOROT(*)
      LOGICAL GAS
C
      DIMENSION JACAC(MAXASH,MAXASH)
C
      CALL QENTER('SETXEE')
C
C
C     ************************************
C     *** Find non-redundant rotations ***
C     ************************************
C
C     No active-active rotations are included here,
C     if not CAS (told with GAS variable) these are added
C     in RASGWOP below.
C     
      NXOPE = 0
      DO 10 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
C        
         ISTI  = NIND(1,1,ISYM) + 1
         IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
C        
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
         JSTS  = NIND(1,3,JSYM) + 1
         JENDS = NIND(1,3,JSYM) + NIND(2,3,JSYM)
C        
         IF(NIND(1,1,ISYM).NE.0) THEN
            WRITE(LUPRI,'(A,I5)')
     &           '*** WARNING(SETXEE): Inactive positrons detected:',
     &           NIND(1,1,ISYM)
         ENDIF
         IF(NIND(2,1,ISYM).EQ.0) GOTO 300
C        
C        Inactive - active elements
C        ==========================
C        
 100     CONTINUE
         IF(NIND(1,2,JSYM).NE.0) THEN
            WRITE(LUPRI,'(A,I5)')
     &           '*** WARNING(SETXEE): Active positrons detected:',
     &           NIND(1,2,JSYM)
         ENDIF
         IF (NIND(2,1,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 200
         DO IN = ISTI,IENDI
            INDI = IND(IN,1,ISYM)
            DO JA = JSTA,JENDA
               INDA = IND(JA,2,JSYM)
               IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDA) .EQ. 0) THEN
                  NXOPE          = NXOPE + 1
                  JXOPE(1,NXOPE) = INDI
                  JXOPE(2,NXOPE) = INDA
               END IF
            ENDDO
         ENDDO
C        
C        Inactive-secondary elements
C        ===========================
C        
 200     CONTINUE
         IF (NIND(2,1,ISYM).EQ.0.OR.NIND(2,3,JSYM).EQ.0) GO TO 300
         DO IN = ISTI,IENDI
            INDI = IND(IN,1,ISYM)
            DO JS = JSTS,JENDS
               INDS = IND(JS,3,JSYM)
               IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDS) .EQ. 0) THEN
                  NXOPE          = NXOPE + 1
                  JXOPE(1,NXOPE) = INDI
                  JXOPE(2,NXOPE) = INDS
               END IF
            ENDDO
         ENDDO
C        
C        Active-secondary elements
C        =========================
C        
 300     CONTINUE
         IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,3,JSYM).EQ.0) GO TO 10
         DO IA = ISTA,IENDA
            INDA = IND(IA,2,ISYM)
            DO JS = JSTS,JENDS
               INDS = IND(JS,3,JSYM)
               IF (NOROT(INDA) .EQ. 0 .AND. NOROT(INDS) .EQ. 0) THEN
                  NXOPE          = NXOPE+1
                  JXOPE(1,NXOPE) = INDA
                  JXOPE(2,NXOPE) = INDS
               END IF
            ENDDO
         ENDDO
 10   CONTINUE
C     
C     
C     **********************************************
C     *** Find redundant active-active rotations ***
C     **********************************************
C     
C     
      NHOPE = NXOPE
      CALL IZERO(JACAC,MAXASH**2)
C     
      IF (GAS) THEN
C        
C        Find non-redundant active-active rotations.
C        
c        CALL RGASWOP(JXOPE,IOPSY,NXOPE,IND,NIND,JACAC,100,NOROT)
         CALL RGASWOP(JXOPE,IOPSY,NXOPE,IND,NIND,JACAC,IPRINT,NOROT)
C        
      END IF
C     
      NHOPE = NXOPE
C     
      DO 20 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
         IF (JSYM .LT. ISYM) GOTO 20
C        
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
C        
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
C        
C        Active-active elements
C        ----------------------
C        
         IF (ISYM .EQ. JSYM) THEN
            IF (NIND(2,2,ISYM).LE.1) GO TO 20
            DO IA = ISTA, IENDA
               IAA = IND(IA,2,ISYM)
               DO JA = IA + 1, IENDA
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.
     &                 NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NHOPE = NHOPE + 1
                     JXOPE(1,NHOPE) = IAA
                     JXOPE(2,NHOPE) = JAA
                  END IF
               END DO
            END DO
         ELSE
            IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 20
            DO IA = ISTA,IENDA
               IAA = IND(IA,2,ISYM)
               DO JA = JSTA,JENDA
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.
     &                 NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NHOPE = NHOPE + 1
                     JXOPE(1,NHOPE) = IAA
                     JXOPE(2,NHOPE) = JAA
                  END IF
               ENDDO
            ENDDO
         END IF
 20   CONTINUE
C     
C     Print section
C     =============
C     
      IF(IPRINT.GE.4) THEN
         CALL HEADER('Output from SETXEE',-1)        
         WRITE(LUPRI,'(A,I10)') 'e-e rotations (non-redundant): ',NXOPE
         WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,NXOPE)
         WRITE(LUPRI,'(A,I10)') 'extra e-e rotations (redundant): ',
     &        NHOPE - NXOPE
         WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=NXOPE+1,NHOPE)
      ENDIF
C     
      CALL QEXIT('SETXEE')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setxep */
      SUBROUTINE SETXEP(JXOPP,IOPSY,NXOPP,QED,IND,NIND,IPRINT,NOROT)
C***********************************************************************
C   
C     Generate electron-positron rotations
C
C     Written by T.Saue July 21 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      LOGICAL QED
      DIMENSION JXOPP(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      DIMENSION NOROT(*)
C
      CALL QENTER('SETXEP')
C
      NXOPP = 0
      DO 10 ISYM = 1,NFSYM
        JSYM  = MOD(ISYM+IOPSY,2) + 1
C
        ISTI  = NIND(1,1,ISYM) + 1
        IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
        ISTA  = NIND(1,2,ISYM) + 1
        IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
C       
C       Inactive-positronic elements
C       ===========================
C       
 100    CONTINUE 
        IF (NIND(2,1,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 200
        DO IN = ISTI,IENDI
           INDI =  IND(IN,1,ISYM)
           DO JP = 1,NIND(1,3,JSYM)
              INDP = IND(JP,3,JSYM)
              IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDP) .EQ. 0) THEN
                 NXOPP          = NXOPP+1
                 JXOPP(1,NXOPP) = INDI
                 JXOPP(2,NXOPP) = INDP
              END IF
           ENDDO
        ENDDO
C       
C       Active-positronic elements
C       ==========================
C       
 200    CONTINUE
        IF (NIND(2,2,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 300
        DO IA = ISTA,IENDA
           INDA = IND(IA,2,ISYM)
           DO JP = 1,NIND(1,3,JSYM)
              INDP = IND(JP,3,JSYM)
              IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDP) .EQ. 0) THEN
                 NXOPP          = NXOPP+1
                 JXOPP(1,NXOPP) = INDA
                 JXOPP(2,NXOPP) = INDP
              END IF
           ENDDO
        ENDDO
C       
C       Secondary-positronic elements
C       =============================
C       
 300    CONTINUE
        IF (QED) THEN
           IF (NIND(1,3,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 10
           ISTS  = NIND(1,3,ISYM) + 1
           IENDS = NIND(1,3,ISYM) + NIND(2,3,ISYM)
           DO IS = ISTS,IENDS
              INDS = IND(IS,1,ISYM)
              DO JP = 1,NIND(1,3,JSYM)
                 INDP = IND(JP,3,JSYM)
                 IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDP) .EQ. 0) THEN
                    NXOPP          = NXOPP+1
                    JXOPP(1,NXOPP) = INDS
                    JXOPP(2,NXOPP) = INDP
                 END IF
              ENDDO
           ENDDO
        ENDIF
C       
 10   CONTINUE
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from SETXEP',-1)        
        WRITE(LUPRI,'(A,I10)') 'e-p rotations: ',NXOPP
        WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,NXOPP)
      ENDIF
C
      CALL QEXIT('SETXEP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck getfck */
      SUBROUTINE GETFCK(FMO,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Get Fock-matrix in MO-basis; use coefficients for MO-transformation
C
C     Written by T.Saue, Sep 16 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      LOGICAL TOBE
      DIMENSION FMO(*),WORK(*)
C
      CALL QENTER('GETFCK')
      KFRSAV = KFREE
C
C     First check if Fock-matrix is on file PAMFCK
C
      LUBUF = 22
      INQUIRE(FILE='PAMFCK',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'PAMFCK','OLD','GETFC1')
C
C       Read FCmo and FVmo
C
        CALL READT(LUBUF,N2ORBXQ,FMO)
        IF (NASHT .GT. 0)
     &     CALL READT(LUBUF,N2ORBXQ,FMO(1+N2ORBXQ))
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
C
         IF (NASHT .EQ. 0) THEN
            CALL MEMGET('REAL',KFAO,N2BBASXQ,WORK,KFREE,LFREE)      
         ELSE
            CALL MEMGET('REAL',KFAO,2*N2BBASXQ,WORK,KFREE,LFREE)      
         END IF
         CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
         IF(SPINFR.OR.LEVYLE) THEN
           CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
         ELSE
           KBEIG = KFREE
         ENDIF
C
         CALL GETFC1(FMO,WORK(KFAO),WORK(KCMO),WORK(KBEIG),
     &             IPRINT,WORK,KFREE,LFREE)
C
         CALL OPNFIL(LUBUF,'PAMFCK','NEW','GETFC1')
         CALL WRITT(LUBUF,N2ORBXQ,FMO)
         IF (NASHT .GT. 0) 
     &      CALL WRITT(LUBUF,N2ORBXQ,FMO(1+N2ORBXQ))
         CLOSE(LUBUF,STATUS='KEEP')
         CALL MEMREL('GETFCK',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
          CALL HEADER('GETFCK:Fock (FCmo) matrix in MO-basis',-1)
          CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),LUPRI)
          IF (NASHT .GT. 0) THEN
             CALL HEADER('GETFCK:Fock (FVmo) matrix in MO-basis',-1)
             CALL PRQMAT(FMO(1+N2ORBXQ),NORBT,NORBT,NORBT,NORBT,
     &                   NZ,IPQTOQ(1,0),LUPRI)
          END IF
      ENDIF
C
      CALL QEXIT('GETFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETFC1 */
      SUBROUTINE GETFC1(FMO,FAO,CMO,IBEIG,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Get Fock-matrix in MO-basis
C
C     Written by T.Saue, Sep 16 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0, DP5 = 0.50D00,D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbfir.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "blocks.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "cbihr2.h"
CTOBS: Note the use of dhfcal ! Should be generalized !!!
      LOGICAL TOBE,BF1,BF2,BF3
      DIMENSION FMO(*),FAO(*),CMO(*),IBEIG(*),WORK(*)
C
#include "ibtfun.h"
C
      CALL QENTER('GETFC1')
      LUBUF = 22
C
C     1. Get coefficients
C     ===================
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','GETFC1')
      IF(SPINFR.OR.LEVYLE) THEN
        IOPT = 10
      ELSE
        IOPT = 2
      ENDIF
      IF (IPRINT.GE.10) IOPT = IOPT + 1
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUM,IBEIG,TOTERG,IOPT,
     &            WORK,KFREE,LFREE)
      CLOSE(LUBUF,STATUS='KEEP')
      NFMAT = 1
      IF (NASHT .GT. 0) NFMAT = 2
C
C     1. Get two-electron Fock matrix in AO-basis
C     ===========================================
C
      INQUIRE(FILE='DFFCK2',EXIST=TOBE)   
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK2','OLD','GETFC1')
        WRITE(LUPRI,'(/A)')
     &    '* GETFCK: reading old two-electron Fock matrix'
        IF (NASHT .EQ. 0) THEN
           CALL REAFCK(LUBUF,FMO,.TRUE.,1)
         ELSE
           CALL REAFCK(LUBUF,FMO,.TRUE.,2)
         END IF
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        INTFLG = INTDEF
Chjaaj MAERKE 01-Feb-2001: this is the best we can do now; 
C      this assumes the user uses the same input for restart
C      with *DHFCAL specified as in the .DHF calculation such that
C      any user-specified change of INTFLG (e.g. no SSSS) is used.
C >>>  THIS MUST BE CHANGED FOR MCSCF !!!!
C
        ILLINT = IBTAND(INTFLG,1)
        ISLINT = IBTAND(INTFLG/2,1)
        ISSINT = IBTAND(INTFLG/4,1)
        WRITE(LUPRI,'(/A,3I2)')
     &    '*** WARNING *** GETFCK: No 2-Fock found.'//
     &    ' Regenerating with .INTFLG =',ILLINT,ISLINT,ISSINT
        NPOS = 0
        IF(PARCAL) THEN
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
C       Closed shell Hartree - Fock or open shell Hartree-Fock
C       ======================================================
CTROND: This must be generalized !!!
C       Consider it done...
        DO IMAT = 1,NFMAT
C          Totally symmetric operator
           ISYMOP(IMAT) = 1
C          Fock matrix type
           IFCKOP(IMAT) = 1
C          Hermitian operator
           IHRMOP(IMAT) = 1
        END DO
        INQUIRE(FILE='DFDENS',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFDENS','OLD','GETFC1')
          CALL READNS(LUBUF,FAO)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
c         CALL DENSTY(FAO,CMO,IPRINT)        
          CALL GENDEN(FAO,CMO,1,IPRINT)
          IF (NASHT .GT. 0) 
     &       CALL GENDEN(FAO(1+N2BBASXQ),CMO,2,IPRINT)
        ENDIF
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMO,FAO,NFMAT,
     &       WORK(KPOS),INTFLG,IPRINT,WORK(KFREE),LFREE)
        IF(PARCAL) THEN
          FIRST1 = BF1
          FIRST2 = BF2
          FIRST3 = BF3
        ENDIF
      ENDIF
C
C     3a. Solvent: get solvent contribution
C     =====================================
C
      IF(SOLVEN) THEN
C        CALL DENSTY(FAO,CMO,IPRINT)        
CJTH check this
        NF = 1
        CALL SOLFCK(FMO,FAO,NF,ESOLVE,ESOLVN,WORK(KFREE),LFREE,IPRSOL)
      ENDIF           
C
C     3. Get one-electron Fock matrix
C     ===============================
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK1','OLD','GETFC1')
        WRITE(LUPRI,'(/A)')
     &    '* GETFCK: reading old one-electron Fock matrix'
        CALL REAFCK(LUBUF,FAO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(/A)')
     &    '*** WARNING *** GETFCK: No 1-Fock found. Regenerating.'
        CALL ONEFCK(FAO,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
C     4. Add 1- and 2-electron matrix and transform to MO-basis
C     =========================================================
C
      CALL DAXPY(N2BBASXQ,D1,FMO,1,FAO,1) 
      IF (NASHT .GT. 0) THEN
         CALL DCOPY(N2BBASXQ,FMO(1+N2BBASXQ),1,FAO(1+N2BBASXQ),1)
C
C        Scale active Fock matrix with 0.5 (all Fock matrices are
C        (multiplied with 2 in TWOFCK)
C
         CALL DSCAL(N2BBASXQ,DP5,FAO(1+N2BBASXQ),1)
      END IF
      IOFF_MO = 0
      IOFF_AO = 0
      DO IMAT = 1,NFMAT
        CALL DZERO(FMO(1+IOFF_MO),N2ORBXQ)
        DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        IF(SPINFR) THEN
          CALL QBTRANS(0,'AOMO','S',D0,
     &      NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &      FAO(I2BASX(I,I)+1+IOFF_AO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &      FMO(I2ORBX(I,I)+1+IOFF_MO),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      IBEIG(IORB(I)+1),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      IBEIG(IORB(I)+1),
     &      WORK(KFREE),LFREE,IPRINT)
        ELSE
          CALL QTRANS('AOMO','S',D0,
     &      NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &      FAO(I2BASX(I,I)+1+IOFF_AO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &      FMO(I2ORBX(I,I)+1+IOFF_MO),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,IPRINT)
        ENDIF
 10     CONTINUE
        IOFF_MO = IOFF_MO + N2ORBXQ
        IOFF_AO = IOFF_AO + N2BBASXQ
      END DO
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GETFC1:Fock (FCmo) matrix in MO-basis',-1)
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &    IPQTOQ(1,0),LUPRI)
        IF (NASHT .GT. 0) THEN
           CALL HEADER('GETFC1:Fock (FVmo) matrix in MO-basis',-1)
           CALL PRQMAT(FMO(1+N2ORBXQ),NORBT,NORBT,NORBT,NORBT,NZ,
     &       IPQTOQ(1,0),LUPRI)
        END IF
      ENDIF
C
      CALL QEXIT('GETFC1')
      RETURN
C
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GETFC1: Coefficient file not found !'
      CALL QUIT('GETFC1: Coefficients not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck unpstr */
      SUBROUTINE UNPSTR(INDSTR,NTYP,IND,LIND,NIND,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     General routine for unpacking of orbital strings
C
C     Written by T.Saue July 19 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"      
      CHARACTER INDSTR(NTYP,NFSYM)*72
      DIMENSION IND(LIND,NTYP,NFSYM),NIND(2,NTYP,NFSYM),
     &          WORK(*)
C
      DO IFRP = 1,NFSYM
        DO ITYP = 1,NTYP
          NVEC = 1
          CALL  NUMLST(INDSTR(ITYP,IFRP),IND(1,ITYP,IFRP),
     &                 NFBAS(IFRP,0),-NFBAS(IFRP,2),NFBAS(IFRP,1),
     &                 IFRP,NVEC,WORK(KFREE),LFREE)
          CALL ORBCNT(IND(1,ITYP,IFRP),NVEC,IFRP,
     &                NIND(1,ITYP,IFRP),NIND(2,ITYP,IFRP))
C
C         Absolute indices of positronic orbitals
C
          IOFF = NPSH(IFRP) + 1 + IORB(IFRP)
          DO I = 1,NIND(1,ITYP,IFRP)
            IND(I,ITYP,IFRP) = IOFF + IND(I,ITYP,IFRP)
          ENDDO
C
C         Absolute indices of electronic orbitals
C
          IOFF = NPSH(IFRP) + IORB(IFRP)
          JOFF = NIND(1,ITYP,IFRP)
          DO I = 1,NIND(2,ITYP,IFRP)
            IND(I+JOFF,ITYP,IFRP) = IOFF + IND(I+JOFF,ITYP,IFRP) 
          ENDDO
        ENDDO
      ENDDO
C
C     Print section
C     =============      
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Output from UNPSTR',-1)
        DO ITYP = 1,NTYP
          WRITE(LUPRI,'(A,I5)') '* Orbital string no.',ITYP
          DO IFRP = 1,NFSYM
            WRITE(LUPRI,'(3X,A,A3)') '* Fermion ircop ',FREP(IFRP)
            IF(NIND(1,ITYP,IFRP).EQ.0) THEN
              WRITE(LUPRI,'(5X,A)') 'No positronic orbitals.'
            ELSE
              WRITE(LUPRI,'(5X,A)') 'Positronic orbitals:'
              WRITE(LUPRI,'(8X,12I5)')
     &             (IND(I,ITYP,IFRP),I=1,NIND(1,ITYP,IFRP))
            ENDIF
            IF(NIND(2,ITYP,IFRP).EQ.0) THEN
              WRITE(LUPRI,'(5X,A)') 'No electronic orbitals.'
            ELSE
              JOFF = NIND(1,ITYP,IFRP)
              WRITE(LUPRI,'(5X,A)') 'Electronic orbitals:'
              WRITE(LUPRI,'(8X,12I5)')
     &             (IND(I+JOFF,ITYP,IFRP),I=1,NIND(2,ITYP,IFRP))
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck REDPAR */
      SUBROUTINE REDPAR(GPCI,GPOE,GPOP,JCONF,JXOPE,JXOPP,
     &                  THRCOM,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Reduce number of parameters by deleting parameters
C     for which the corresponding gradient element is
C     below threshold THRCOM
C
C     Written by T.Saue - July 22 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DC = 100.0D0)
C
#include "dcbxrs.h"
#include "dgroup.h"
      DIMENSION GPCI(*),GPOE(*),GPOP(*),WORK(*),
     &          JCONF(*),JXOPE(2,NZXOPE),JXOPP(2,NZXOPP)
C
      KFRSAV = KFREE
      IF(IPRINT.GE.4) CALL HEADER('Output from REDPAR',-1)        

C
C     Configurational part
C     ====================
C
      IF(NZCONF.GT.0) THEN
        WRITE(LUPRI,'(A)') '* REDPAR: No configurational code.'
        CALL QUIT('REDPAR: No conf.code')
      ENDIF
C
C     e-e orbital rotations
C     =====================
C
      IF(NZXOPE.GT.0) THEN
        CALL MEMGET('REAL',KBUF,NZXOPEQ,WORK,KFREE,LFREE)
        II = 0
        DO I = 1,NZXOPE
          DTST = D0
          IZOFF = 0
          DO IZ = 1,NZ
            DTST = DTST + GPOE(I+IZOFF)*GPOE(I+IZOFF)  
            IZOFF = IZOFF + NZXOPE
          ENDDO
          IF(SQRT(DTST).GE.THRCOM) THEN
            II = II + 1
            JXOPE(1,II) = JXOPE(1,I)
            JXOPE(2,II) = JXOPE(2,I)
            IOFF = KBUF-1
            JOFF = 0
            DO IZ = 1,NZ
              WORK(IOFF+II) = GPOE(JOFF+I)
              IOFF = IOFF + NZXOPE
              JOFF = JOFF + NZXOPE
            ENDDO
          ENDIF
        ENDDO
        IF(II.LT.NZXOPE) THEN
          IOFF = KBUF
          JOFF = 1
          DO IZ = 1,NZ
            CALL DCOPY(II,WORK(IOFF),1,GPOE(JOFF),1)
            IOFF = IOFF + NZXOPE
            JOFF = JOFF + II
          ENDDO
          IF(IPRINT.GE.4) THEN
            DTST = DC*(NZXOPE-II)/NZXOPE
            WRITE(LUPRI,'(A,I10,A,F5.1,A)') 
     &            'e-e rotations: ',II,' (reduction: ',DTST,'%)'
            WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,II)
            CALL PRQMAT(GPOE,II,1,II,1,NZ,1,LUPRI)
          ENDIF
          NZXOPE  = II
          NZXOPEQ = NZXOPE*NZ
        ENDIF
        CALL MEMREL('REDPAR.ee',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF      
C
C     e-p orbital rotations
C     =====================
C
      IF(NZXOPP.GT.0) THEN
        CALL MEMGET('REAL',KBUF,NZXOPPQ,WORK,KFREE,LFREE)
        II = 0
        DO I = 1,NZXOPP
          DTST = D0
          IZOFF = 0
          DO IZ = 1,NZ
            DTST = DTST + GPOP(I+IZOFF)*GPOP(I+IZOFF)  
            IZOFF = IZOFF + NZXOPP
          ENDDO
          IF(SQRT(DTST).GE.THRCOM) THEN
            II = II + 1
            JXOPP(1,II) = JXOPP(1,I)
            JXOPP(2,II) = JXOPP(2,I)
            IOFF = KBUF-1
            JOFF = 0
            DO IZ = 1,NZ
              WORK(IOFF+II) = GPOP(JOFF+I)
              IOFF = IOFF + NZXOPP
              JOFF = JOFF + NZXOPP
            ENDDO
          ENDIF
        ENDDO
        IF(II.LT.NZXOPP) THEN
          IOFF = KBUF
          JOFF = 1
          DO IZ = 1,NZ
            CALL DCOPY(II,WORK(IOFF),1,GPOP(JOFF),1)
            IOFF = IOFF + NZXOPP
            JOFF = JOFF + II
          ENDDO
          IF(IPRINT.GE.4) THEN
            DTST = DC*(NZXOPP-II)/NZXOPP
            WRITE(LUPRI,'(A,I10,A,F5.1,A)') 
     &            'e-p rotations: ',II,'(reduction: ',DTST,'%)'
            WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,II)
          ENDIF
          NZXOPP  = II
          NZXOPPQ = NZXOPP*NZ
        ENDIF
        CALL MEMREL('REDPAR.ep',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF      
C
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT *NZ
C
C     Total number of parameters
C     ==========================
C
      NZVAR    = NZCONF + NZXOPT
      NZVARQ   = NZVAR*NZ
      RETURN
      END
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DEFORB */
      SUBROUTINE DEFORB(IORBCL,NIND,NTYP,INDSTR,WORK,KFREE,LFREE,
     &                  IPRINT)
C***********************************************************************
C
C     This routine reads orbital strings and classifies
C     and indexes orbitals
C
C     Written by T.Saue Oct 3 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER INDSTR(NTYP,2)*72
      DIMENSION IORBCL(4,NORBT),NIND(2,NTYP,2),WORK(*)
C
      CALL QENTER('DEFORB')
      KFRSAV = KFREE
C
C     Allocate buffer areas
C     =====================
C
      CALL MEMGET('INTE',KIND ,MXFBAS,WORK,KFREE,LFREE)
      CALL DEFOR1(IORBCL,NIND,NTYP,INDSTR,WORK(KIND),
     &            WORK,KFREE,LFREE,IPRINT)
      CALL MEMREL('DEFORB',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DEFORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DEFOR1 */
      SUBROUTINE DEFOR1(IORBCL,NIND,NTYP,INDSTR,IND,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     This routine reads orbital strings and classifies
C     and indexes orbitals
C
C     Three types are defined:
C     ITYP = 1   inactive
C     ITYP = 2   active
C     ITYP = 3   secondary
C
C     NIND(1,ITYP,IFRP) - number of positronic orbitals 
C     NIND(2,ITYP,IFRP) - number of electronic orbitals
C       ...of type ITYP and fermion ircop IFRP  
C
C     IORB is the absolute index of the orbital
C     IORBCL(1,IORB) - ITYP
C     IORBCL(2,IORB) - index within fermion irrep
C     IORBCL(3,IORB) - fermion irrep IFRP
C     IORBCL(4,IORB) - absolute index for given type ITYP
C                      and electron/positron
C
C     Written by T.Saue Oct 3 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER INDSTR(NTYP,2)*72
      DIMENSION IORBCL(4,NORBT),IND(MXFBAS),
     &          NIND(2,NTYP,2),WORK(KFREE)
C
      NDIM = NORBT*3
      CALL IZERO(IORBCL,NDIM)
C
C     Unpack strings
C     ==============
C
      DO ITYP = 1,NTYP
        KOFE = 0
        KOFP = 0
        DO IFRP = 1,NFSYM
          NVEC = 1
          CALL NUMLST(INDSTR(ITYP,IFRP),IND,
     &                NFBAS(IFRP,0),-NFBAS(IFRP,2),NFBAS(IFRP,1),
     &                IFRP,NVEC,WORK(KFREE),LFREE)
          CALL ORBCNT(IND,NVEC,IFRP,NIND(1,ITYP,IFRP),NIND(2,ITYP,IFRP))
C
C         Absolute indices of positronic orbitals
C
          IOFF = NPSH(IFRP)+IORB(IFRP)+1
          DO I = 1,NIND(1,ITYP,IFRP)
            J           = IOFF+IND(I)
            IORBCL(1,J) = ITYP
            IORBCL(2,J) = NPSH(IFRP)+IND(I)+1
            IORBCL(3,J) = IFRP
            KOFP = KOFP + 1
            IORBCL(4,J) = KOFP
            
          ENDDO
C
C         Absolute indices of electronic orbitals
C
          IOFF = NPSH(IFRP)+IORB(IFRP)
          JOFF = NIND(1,ITYP,IFRP)
          DO I = 1,NIND(2,ITYP,IFRP)
            II = I + JOFF
            J = IOFF + IND(II) 
            IORBCL(1,J) = ITYP
            IORBCL(2,J) = IND(II)
            IORBCL(3,J) = IFRP
            KOFE = KOFE + 1
            IORBCL(4,J) = KOFE
          ENDDO
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from DEFORB',-1)
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A3)') '* Fermion ircop:',FREP(IFRP)
          WRITE(LUPRI,'(5X,3I5)')
     &    (I,IORBCL(1,IORB(IFRP)+I),IORBCL(2,IORB(IFRP)+I),
     &     I = 1,NORB(IFRP))
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rgaswop */
      SUBROUTINE RGASWOP(JXOPE,IOPSY,NXOPE,IND,NIND,JACAC,IPRINT,NOROT)
C***********************************************************************
C   
C     Generate non-redundant active-active rotations.
C
C     Written by J. Thyssen - Jan 18 2001
C
C     Revised GAS setup, Timo Fleig, May 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dcbibn.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
C
      DIMENSION JXOPE(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      DIMENSION JACAC(MAXASH,MAXASH)
      DIMENSION NOROT(*)
C
      DIMENSION IGAS(MAXASH)
      LOGICAL GASCAS(MXMCGAS)
C
C
c     iprint = 20
      NXOPES = NXOPE
C
C     Set up help arrays
C     ------------------
C
C     IGAS point from active orbital to which GAS space it belongs.
C
      II = 0
      DO ISYM = 1, NFSYM
         DO I = 1, NGAS
            DO J = 1, NGSH(ISYM,I)
               II = II + 1
               IGAS(II) = I
            END DO
         END DO
      END DO
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Output from RGASWOP',-1)
         WRITE(LUPRI,'(1X,A)')
     &        'Pointer from active orbitals to GAS spaces:'
         DO I = 1, NASHT
            WRITE(LUPRI,'(5X,A,I3,A,I3)')
     &           'Orbital ',I,' is in GAS space ',IGAS(I)
         END DO
      END IF
C
C     Determine whether a GAS space is a CAS space:
C
      DO I = 1, NGAS
C
C        Number of orbitals in this GAS shell
C
         NGST = 0
         DO ISYM = 1, NFSYM
            NGST = NGST + NGSH(ISYM,I)
         END DO
C
C        Minimum number of electrons in this GAS space:
C
         MINCAS = 2*NGST - (NGSP(2,I) - NGSP(1,I)) 
         IF (I.EQ.NGAS.AND.NGAS.GT.1) THEN
           MINCAS = NAELEC - NGSP(2,I-1)
         ELSE IF (I.EQ.NGAS.AND.NGAS.EQ.1) then
           MINCAS = NAELEC
         END IF
c        print*,'GAS ',I,' MINCAS,MAXCAS ',MINCAS,MAXCAS
C
C        Maximum number of holes in this GAS space:
C
         MAXCAS = 2*NGST - MINCAS
C
C        This is a GAS space if the minimum number of electrons 
C        after this space is less than MINCAS - and -
C        if the maximum number of holes in this GAS space is larger
C        than MAXCAS:
C
         IF ( (NGSP(1,I) .LE. MINCAS .AND.
     &        (NGSP(2,I) - NGSP(1,I)) .GE. MAXCAS) .OR.
     &        I .EQ. NGAS ) THEN
            GASCAS(I) = .TRUE.
         ELSE
            GASCAS(I) = .FALSE.
         END IF
      END DO
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(1X,A)')
     &        'GAS space is a CAS space?'
         DO I = 1, NGAS
            WRITE(LUPRI,'(5X,A,I3,A,L1)')
     &           'GAS space ',I,' : ',GASCAS(I)
         END DO
      END IF
C
C     Algorithm:
C     ----------
C
C     No inter-GAS-space active-active rotations
C     No active-active rotations into CAS spaces
C     No active-active rotations from a CAS space into the NEXT GAS space.
C
C     Active-active elements
C     ----------------------
C
      DO 20 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
         IF (JSYM .LT. ISYM) GOTO 20
C        
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
C        
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
         IF (ISYM .EQ. JSYM) THEN
            IF (NIND(2,2,ISYM).LE.1) GO TO 20
            DO IA = ISTA, IENDA
               DO JA = IA + 1, IENDA
                  IAA = IND(IA,2,ISYM)
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.
     &                 IGAS(IDXG2U(IAA)) .NE. IGAS(IDXG2U(JAA)) .AND.
     &                 .NOT. GASCAS(IGAS(IDXG2U(JAA))) .AND.
     &                 .NOT. (
     &                    GASCAS(IGAS(IDXG2U(IAA))) .AND.
     &                 IGAS(IDXG2U(JAA)) .LE. IGAS(IDXG2U(IAA)+1)) .AND.
     &                 NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NXOPE = NXOPE + 1
                     JXOPE(1,NXOPE) = IAA
                     JXOPE(2,NXOPE) = JAA
                     JACAC(IDXG2U(IAA),IDXG2U(JAA)) = 1
                     JACAC(IDXG2U(JAA),IDXG2U(IAA)) = 1
                  END IF
               END DO
            END DO
         ELSE
            IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 20
            DO IA = ISTA,IENDA
               DO JA = JSTA,JENDA
                  IAA = IND(IA,2,ISYM)
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.
     &                 IGAS(IDXG2U(IAA)) .NE. IGAS(IDXG2U(JAA)) .AND.
     &                 .NOT. GASCAS(IDXG2U(JAA)) .AND.
     &                 .NOT. (
     &                    GASCAS(IDXG2U(IAA)) .AND.
     &                 IGAS(IDXG2U(JAA)) .LE. IGAS(IDXG2U(IAA)+1))) THEN
                     NXOPE = NXOPE + 1
                     JXOPE(1,NXOPE) = IAA
                     JXOPE(2,NXOPE) = JAA
                     JACAC(IDXG2U(IAA),IDXG2U(JAA)) = 1
                     JACAC(IDXG2U(JAA),IDXG2U(IAA)) = 1
                  END IF
               ENDDO
            ENDDO
         END IF
 20   CONTINUE
C
C     
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(1X,A,I3,A)')
     &        'The following ',NXOPE - NXOPES,
     &        ' active-active rotations were considered non-redundant: '
         DO I = NXOPES + 1, NXOPE
            WRITE(LUPRI,'(5X,2I4)')
     &           JXOPE(1,I), JXOPE(2,I)
         END DO
      END IF
C     stop 'after GASCAS determination.'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Ostring */
      SUBROUTINE OSTRING(INDSTR,SKIPEE,SKIPEP,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Set orbital strings         
C     Stolen by T.Saue Aug 24 2000
C
C      Revisions: Aug 2002 hjaaj
C                 May 2003 MI
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
C
      CHARACTER*72 INDSTR(3,2)
      LOGICAL SKIPEE,SKIPEP,ORBDEF(3)
      DIMENSION WORK(*),IR(2,2),NR(2)

      CALL QENTER('OSTRING')
C
      DO I = 1,3
         ORBDEF(I) = .FALSE.
         DO J = 1, NFSYM
            ORBDEF(I) = ORBDEF(I) .OR. (INDSTR(I,J) .NE. ' ')
         END DO
      END DO
C
C     Inactive orbitals...
C
      IF (.NOT. ORBDEF(1)) THEN
         DO I = 1,NFSYM
            INDSTR(1,I) = ' '
            IF (NISH(I).GT.0) THEN
CMI            WRITE(INDSTR(1,I),'(I4,A2,I4)') 1,'..',NISH(I)
               WRITE(INDSTR(1,I),'(I6,A2,I6)') 1,'..',NISH(I)
            ENDIF
         ENDDO
      END IF
C
C     Active orbitals...
C
      IF (.NOT. ORBDEF(2)) THEN
         DO I = 1,NFSYM
            INDSTR(2,I) = ' '
            IF (NASH(I).GT.0) THEN
CMI            WRITE(INDSTR(2,I),'(I4,A2,I4)')
               WRITE(INDSTR(2,I),'(I6,A2,I6)')
     &               NISH(I)+1,'..',NISH(I)+NASH(I)
            ENDIF
         ENDDO
      END IF
C
C     Secondary orbitals...
C
      IF (.NOT. ORBDEF(3)) THEN
         IF (SKIPEP.AND.SKIPEE) THEN
            DO I = 1,NFSYM
               INDSTR(3,I) = ' '
            ENDDO
         ELSEIF (SKIPEP) THEN
            DO I = 1,NFSYM
               INDSTR(3,I) = ' '
               IF (NSSH(I).GT.0) THEN
CMI               WRITE(INDSTR(3,I),'(I4,A2,I4)')
                  WRITE(INDSTR(3,I),'(I6,A2,I6)')
     &                  NOCC(I)+1,'..',NESH(I)
               ENDIF
            ENDDO
         ELSEIF (SKIPEE) THEN
            DO I = 1,NFSYM
               INDSTR(3,I) = ' '
CMI            WRITE(INDSTR(3,I),'(I4,A2,I4)')
               WRITE(INDSTR(3,I),'(I6,A2,I6)')
     &             -NPSH(I),'..',-1
            ENDDO
         ELSE
            DO I = 1,NFSYM
               INDSTR(3,I) = ' '
               IF (NSSH(I).GT.0) THEN
                  IF (NPSH(I) .GT. 0) THEN
CMI                  WRITE(INDSTR(3,I),'(I4,A5,I4,A2,I4)')
                     WRITE(INDSTR(3,I),'(I6,A5,I6,A2,I6)')
     &                   -NPSH(I),'..-1,',NOCC(I)+1,'..',NESH(I)
                  ELSE
CMI                  WRITE(INDSTR(3,I),'(I4,A2,I4)')MI
                     WRITE(INDSTR(3,I),'(I6,A2,I6)')
     &                     NOCC(I)+1,'..',NESH(I)
                  END IF
               ELSE
                  IF (NPSH(I) .GT. 0) THEN
CMI                  WRITE(INDSTR(3,I),'(I4,A2)') -NPSH(I),'..-1'
                     WRITE(INDSTR(3,I),'(I6,A2)') -NPSH(I),'..-1'
                  END IF
               END IF
            ENDDO
         ENDIF
      END IF
C
C
C
C      Definition of orbital classes
C      -----------------------------
C
      WRITE(LUPRI,'(/4X,A)') '* Orbital classes:'
      DO I = 1,NFSYM
        IF(SKIPEP) THEN
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ELSE
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDIF
        WRITE(LUPRI,'(4X,A,A3)') 'Fermion ircop:',FREP(I)
C
C       Inactive orbitals
C       -----------------
C
        NVEC = 0
        CALL  NUMLST(INDSTR(1,I),IDUMMY,NR(I),
     &             IR(1,I),IR(2,I),I,NVEC,WORK(KFREE),LFREE)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A)') 
     &       '- Inactive orbitals  : No electrons'
        ELSE
          WRITE(LUPRI,'(4X,A,A56)') 
     &       '- Inactive orbitals  : ', INDSTR(1,I)
        ENDIF
C
C       Active orbitals
C       ---------------
C
        NVEC = 0
        CALL  NUMLST(INDSTR(2,I),IDUMMY,NR(I),
     &             IR(1,I),IR(2,I),I,NVEC,WORK(KFREE),LFREE)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A)') 
     &       '- Active orbitals    : No electrons'
        ELSE
          WRITE(LUPRI,'(4X,A,A56)') 
     &       '- Active orbitals    : ', INDSTR(2,I)
        ENDIF
C
C       Secondary orbitals
C       ------------------
C
        NVEC = 0
        CALL  NUMLST(INDSTR(3,I),IDUMMY,NR(I),
     &            IR(1,I),IR(2,I),I,NVEC,WORK(KFREE),LFREE)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A)') 
     &       '- Secondary orbitals : No electrons'
        ELSE
          WRITE(LUPRI,'(4X,A,A56)') 
     &       '- Secondary orbitals : ', INDSTR(3,I)
        ENDIF
      ENDDO
C
      CALL QEXIT('OSTRING')
      RETURN
      END
