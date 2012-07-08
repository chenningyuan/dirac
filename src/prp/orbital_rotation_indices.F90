!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck setxop */
      SUBROUTINE SETXOP(INDSTR,SKIPEE,SKIPEP,GAS,IOPSY,NOROT,           &
                   NXOPE,NHOPE,KXOPE, NXOPP,KXOPP,                      &
                   WORK,KFREE,LFREE,IPRINT)
!***********************************************************************
!
!     Set orbital rotation index vectors
!
!     Orbitals are classified as:
!       inactive  : doubly occupied in all determinants
!       active    : singly occupied in some determinants
!       secondary : unoccupied in all determinants
!       positronic: positronic orbitals
!
!     Ordering is:
!       e-e rotations: (skipped if SKIPEE = T)
!          IA : inactive-active   
!          IS : inactive-secondary
!          AS : active-secondary
!          AB : active-active
!       e-p rotations: (skipped if SKIPEP = T)
!          IP : inactive-positronic
!          AP : active-positronic
!         (SP : secondary-positronic   -  only for QED=T)
!
!     NXOPE is the number of non-redundant e-e rotations.
!     NHOPE is the total number of rotations.
!
!     Written by  Hans Joergen Aa. Jensen and T.Saue 27-Jan-1995
!     Last revision: July 21 1997 - tsaue
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
! Used from common blocks:
!  DGROUP: NFSYM
!  DCBHAM: QED
!  DCBORB: IORB,NPSH,NFRO,NISH,NASH
!
#include "dgroup.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbbas.h"
!
      LOGICAL SKIPEP,SKIPEE,GAS
      CHARACTER INDSTR(3,2)*72
      DIMENSION WORK(*),NSTR(0:2,3,2)
      DIMENSION NOROT(*)
!
      CALL QENTER('SETXOP')
      KXOPE = KFREE
      KXOPP = KFREE
      NXOPE = 0
      NXOPP = 0
      NHOPE = 0
!
!     Allocate buffer areas
!     =====================
!
      NBUF = N2BBASX
!
!     ***********************************
!     *** Electron-electron rotations ***
!     ***********************************
!      
 10   CONTINUE
      IF(SKIPEE) GOTO 20
      CALL MEMGET('INTE',KBUF ,NBUF          ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,MXFBAS*3*NFSYM,WORK,KFREE,LFREE)
!     ... Unpack strings
      CALL UNPSTR(INDSTR,3,WORK(KIBUF),MXFBAS,NSTR,                     &
             WORK,KFREE,LFREE,IPRINT)
      CALL SETXEE(WORK(KBUF),IOPSY,NXOPE,NHOPE,WORK(KIBUF),             &
      NSTR,GAS,IPRINT,NOROT)
      IF(NHOPE.GT.0) THEN
        CALL MEMREL('SETXOP.e-e',WORK,1,KBUF,KFREE,LFREE)
        CALL MEMGET('INTE',KXOPE,2*NHOPE,WORK,KFREE,LFREE)
        NBUF = N2BBASX - 2*(NHOPE+1)
      ENDIF
!
!     ***********************************
!     *** Electron-positron rotations ***
!     ***********************************
!      
 20   CONTINUE
      IF(SKIPEP) GOTO 30
      CALL MEMGET('INTE',KBUF ,NBUF          ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,MXFBAS*3*NFSYM,WORK,KFREE,LFREE)
!     ... Unpack strings
      CALL UNPSTR(INDSTR,3,WORK(KIBUF),MXFBAS,NSTR,                     &
             WORK,KFREE,LFREE,0)
      CALL SETXEP(WORK(KBUF),IOPSY,NXOPP,QED,WORK(KIBUF),NSTR,IPRINT,   &
      NOROT)
      CALL MEMREL('SETXOP.e-p',WORK,1,KBUF,KFREE,LFREE)
      CALL MEMGET('INTE',KXOPP,2*NXOPP,WORK,KFREE,LFREE)
!
 30   CONTINUE
!
      CALL QEXIT('SETXOP')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck setxee */
      SUBROUTINE SETXEE(JXOPE,IOPSY,NXOPE,NHOPE,IND,NIND,GAS,IPRINT,    &
                   NOROT)
!   
!     Generate electron-electron rotations
!
!     Written by T.Saue July 18 1997
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
!
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbibn.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbxrs.h"
!
      DIMENSION JXOPE(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      DIMENSION NOROT(*)
      LOGICAL GAS
!
      DIMENSION JACAC(MAXASH,MAXASH)
!
      CALL QENTER('SETXEE')
!
!
!     ************************************
!     *** Find non-redundant rotations ***
!     ************************************
!
!     No active-active rotations are included here,
!     if not CAS (told with GAS variable) these are added
!     in RASGWOP below.
!     
      NXOPE = 0
      DO 10 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
!        
         ISTI  = NIND(1,1,ISYM) + 1
         IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!        
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
         JSTS  = NIND(1,3,JSYM) + 1
         JENDS = NIND(1,3,JSYM) + NIND(2,3,JSYM)
!        
         IF(NIND(1,1,ISYM).NE.0) THEN
            WRITE(LUPRI,'(A,I5)')                                       &
            '*** WARNING(SETXEE): Inactive positrons detected:',        &
            NIND(1,1,ISYM)
         ENDIF
         IF(NIND(2,1,ISYM).EQ.0) GOTO 300
!
!        For Sternheimer complement add inactive-inactive elements
!        =========================================================
         IF(STERNC) THEN
           JSTI  = NIND(1,1,JSYM) + 1
           JENDI = NIND(1,1,JSYM) + NIND(2,1,JSYM)
           DO IN = ISTI,IENDI
             INDI = IND(IN,1,ISYM)
             DO JA = JSTI,JENDI
               INDA = IND(JA,1,JSYM)
               IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDA) .EQ. 0) THEN
                  NXOPE          = NXOPE + 1
                  JXOPE(1,NXOPE) = INDI
                  JXOPE(2,NXOPE) = INDA
               END IF
             ENDDO
           ENDDO
         ENDIF
!        
!        Inactive - active elements
!        ==========================
!        
 100     CONTINUE
         IF(NIND(1,2,JSYM).NE.0) THEN
            WRITE(LUPRI,'(A,I5)')                                       &
            '*** WARNING(SETXEE): Active positrons detected:',          &
            NIND(1,2,JSYM)
         ENDIF
!                 inactive               active
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
!        
!        Inactive-secondary elements
!        ===========================
!        
 200     CONTINUE
!                 inactive               secondary
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
!        
!        Active-secondary elements
!        =========================
!        
 300     CONTINUE
!                 active                 secondary
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
!     
!     
!     **********************************************
!     *** Find redundant active-active rotations ***
!     **********************************************
!     
!     
      NHOPE = NXOPE
      CALL IZERO(JACAC,MAXASH**2)
!     
      IF (GAS) THEN
!        
!        Find non-redundant active-active rotations.
!        
         CALL RGASWOP(JXOPE,IOPSY,NXOPE,IND,NIND,JACAC,IPRINT,NOROT)
!        
      END IF
!     
      NHOPE = NXOPE
!     
      DO 20 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
         IF (JSYM .LT. ISYM) GOTO 20
!        
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!        
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
!        
!        Active-active elements
!        ----------------------
!        
         IF (ISYM .EQ. JSYM) THEN
            IF (NIND(2,2,ISYM).LE.1) GO TO 20
            DO IA = ISTA, IENDA
               IAA = IND(IA,2,ISYM)
               DO JA = IA + 1, IENDA
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.       &
                  NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
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
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.       &
                  NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NHOPE = NHOPE + 1
                     JXOPE(1,NHOPE) = IAA
                     JXOPE(2,NHOPE) = JAA
                  END IF
               ENDDO
            ENDDO
         END IF
 20   CONTINUE
!     
!     Print section
!     =============
!     
      IF(IPRINT.GE.2) THEN
        CALL HEADER('Output from SETXEE',-1)        
        WRITE(LUPRI,'(A,I10)')                                          &
        '# of e-e rotations (non-redundant): ',NXOPE
        IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,NXOPE)
        ENDIF
        WRITE(LUPRI,'(A,I10)') '# of extra e-e rotations (redundant): ',&
         NHOPE - NXOPE
        IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=NXOPE+1,NHOPE)
        ENDIF
      ENDIF
!     
      CALL QEXIT('SETXEE')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck setxep */
      SUBROUTINE SETXEP(JXOPP,IOPSY,NXOPP,QED,IND,NIND,IPRINT,NOROT)
!***********************************************************************
!   
!     Generate electron-positron rotations
!
!     Written by T.Saue July 21 1997
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
!
      LOGICAL QED
      DIMENSION JXOPP(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      DIMENSION NOROT(*)
!
      CALL QENTER('SETXEP')
!
      NXOPP = 0
      INDI  = 0 
      DO 10 ISYM = 1,NFSYM
        JSYM  = MOD(ISYM+IOPSY,2) + 1
!
        ISTI  = NIND(1,1,ISYM) + 1
        IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
        ISTA  = NIND(1,2,ISYM) + 1
        IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!       
!       Inactive-positronic elements
!       ===========================
!       
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
!       
!       Active-positronic elements
!       ==========================
!       
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
!       
!       Secondary-positronic elements
!       =============================
!       
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
!       
 10   CONTINUE
!
!     Print section
!     =============
!
      IF(IPRINT.GE.2) THEN
        CALL HEADER('Output from SETXEP',-1)        
        WRITE(LUPRI,'(A,I10)') '# of e-p rotations: ',NXOPP
        IF(IPRINT.GE.3) THEN
          WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,NXOPP)
        ENDIF
      ENDIF
!
      CALL QEXIT('SETXEP')
      RETURN
      END
