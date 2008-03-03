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
C VERSION : $Revision$
C DATE    : $Date$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetno */
      SUBROUTINE RGETNO(CREF,OCCNO,UNO,CMO,NATORB,MOISNO,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate natural orbital occupation numbers and, 
C     if NATORB is true, transform orbital to natural orbitals.
C
C     Input:
C        CREF - CI reference vector
C        CMO  - MO orbitals
C        NATORB - transform to natural orbitals
C
C     Output:
C        OCCNO - natural orbital occupation numbers
C        UNO   - transformation matrix
C        MOISNO - MO orbitals was already NO orbitals.
C        CMO   - if NATORB, then NO orbitals, else unchanged
C
C     Written by J. Thyssen - Dec 21 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION CREF(*)
      DIMENSION OCCNO(*), UNO(NASHT,NASHT,*)
      DIMENSION CMO(*)
      DIMENSION WORK(*)
C
      LOGICAL NATORB, MOISNO
C
      CALL QENTER('RGETNO')
      KFRSAV = KFREE
C
C     Calculate DV density matrix
C     ---------------------------
C
      CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NASHT,WORK,KFREE,LFREE)
      CALL RMAKDV(WORK(KDV),CREF,WORK,KFREE,LFREE)
C
C     Calculate max off-diagonal elements of density matrix
C
      DVMAX = D0
      DO I = 2, NASHT
         DO J = 1, I - 1
            DVMAX = MAX(DVMAX, ABS(WORK(KDV + (J-1) + (I-1)*NASHT)))
         END DO
      END DO
      DO IZ = 2, NZ
         DO I = 1, NASHT
            DO J = 1, I
               DVMAX = MAX(DVMAX,
     &              ABS(WORK(KDV + (J-1) + (I-1)*NASHT+(IZ-1)*N2ASHX)))
            END DO
         END DO
      END DO
      IF (DVMAX .LE. THRZER) THEN
         MOISNO = .TRUE.
         IF (IPRINT .GE. 2) WRITE(LUPRI,'(1X,A)')
     &      '(RGETNO) The orbitals are already natural orbitals.'
         CALL DZERO(UNO,N2ASHXQ)
         CALL DUNIT(UNO,NASHT)
         DO I = 1, NASHT
            WORK(KEIG+I-1) = WORK(KDV+(I-1)+(I-1)*NASHT)
         END DO
      ELSE
         MOISNO = .FALSE.
      END IF
C
C     Diagonalize density matrix
C     --------------------------
C
      IF (.NOT. MOISNO) THEN
         IF (NZ .EQ. 1) THEN
            IJOB = 1
            IORDER = 0
            IPACK = 0
            CALL RSJACO(NASHT,NASHT,NASHT,WORK(KDV),WORK(KEIG),
     &           IJOB,IORDER,IPACK,UNO,WORK(KFREE),LFREE)
         ELSE
            MATZ = 1
            IF( NATOLCR ) MATZ = 0
            CALL QDIAG(NZ,NASHT,WORK(KDV),NASHT,NASHT,
     &           WORK(KEIG),MATZ,UNO,NASHT,NASHT,WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/,1X,2A,I4)')
     &              '*** ERROR in RGETNO ***: ',
     &              'QDIAG failed with error code ',IERR
               CALL QUIT('*** ERROR in RGETNO ***')
            END IF
         END IF
C
C        FIXME: select positive phase of transformation vectors
C
      END IF
C
      JOCC = 0
      JACT = 0
      JACT_TOT = 0
*
*     FIXME: right order of active orbitals if NZ == 2 and 
*            NO occ. numbers are calculated with LUCIAREL
*
      IF( NATOLCR .and. NZ .eq. 2 ) THEN
        DO I = 1, NFSYM
           DO J = 1, NISH(I)
              JOCC = JOCC + 1
              OCCNO(JOCC) = D2
           END DO
           JACT_TOT = JACT_TOT + NASH(I)
           DO J = 1, NASH(I)
              JOCC = JOCC + 1
              JACT = JACT + 1
              OCCNO(JOCC) = WORK(KEIG + JACT_TOT - JACT )
           END DO
        END DO
      ELSE
        DO I = 1, NFSYM
           DO J = 1, NISH(I)
              JOCC = JOCC + 1
              OCCNO(JOCC) = D2
           END DO
           DO J = 1, NASH(I)
              JOCC = JOCC + 1
              JACT = JACT + 1
              OCCNO(JOCC) = WORK(KEIG + JACT - 1)
           END DO
        END DO
      END IF
C
      IF ( IPRINT .GE. 1 .or. NATOLCR ) THEN
         IF( NATOLCR )THEN
            WRITE(LUPRI,'(//A/A/A)')
     &' **************************************************************',
     &' *************** CI natural orbital calculation ***************',
     &' **************************************************************'
         END IF
         ICI_RT = IDENSLR_STATE
         WRITE(LUPRI,'(//A,I3,/A)')
     &        '   Natural orbital occupation numbers for state:',ICI_RT,
     &        '   ------------------------------------------------'
         DO I = 1, NFSYM
            NOCCI = NOCC(I)
            IF (NOCCI .EQ. 0) THEN
               WRITE(LUPRI,9001) I
            ELSE
               IOCCI = IOCC(I)
               OCCSUM = DSUM(NOCCI,OCCNO(IOCCI+1),1)
               WRITE(LUPRI,9002) I, (OCCNO(IOCCI+J),J=1,NOCCI)
               WRITE(LUPRI,9003) OCCSUM
            END IF
         END DO
      END IF
      IPRINT = 0
C
 9001 FORMAT(/'   Symmetry',I3,//,'   No occupied orbitals')
 9002 FORMAT(/'   Symmetry',I3,//,(5F14.9))
 9003 FORMAT(/'   Sum =',T15,F14.9)
C
      IF ( IPRINT .GE. 5 .AND. .NOT. MOISNO ) THEN
         DO I = 1, NFSYM
            WRITE(LUPRI,'(/1X,A,I3)')
     &           '(RGETNO) Transformation to NO '//
     &           '(active orb.) for symmetry ',I
            CALL PRQMAT(UNO(1+IASH(I),1+IASH(I),1),NASH(I),NASH(I),
     &           NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
      CALL MEMREL('RGETNO.diag',WORK,1,KFRSAV,KFREE,LFREE)
C
C     Transform to natural orbitals
C     -----------------------------
C
      IF (.NOT. MOISNO) THEN
         IF (.NOT. NATORB) THEN
            IF ( IPRINT .GE. 1 )
     &           WRITE(LUPRI,'(/1X,A)')
     &           '(RGETNO) active orbitals are NOT transformed '//
     &           'to natural orbitals.'
         ELSE
            IF ( IPRINT .GE. 1 )
     &           WRITE(LUPRI,'(/1X,A)')
     &           '(RGETNO) active orbitals are transformed '//
     &           'to natural orbitals.'
C
            CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
            CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
            DO 10 I = 1, NFSYM
               IF (NASH(I) .EQ. 0) GO TO 10
C
               IF ( IPRINT .GE. 20 ) THEN
                  WRITE(LUPRI,'(/1X,A,I3)')
     &                 '(RGETNO) Before transformation: '//
     &                 'MO coefficients, symmetry ',I
                  CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                 NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Transform active orbitals:
C
C              CMO(NO) = CMO * UNO
C
               NBASI = NFBAS(I,0)
               NASHI = NASH(I)
               NORBI = NORB(I)
               JUNO = 1 + IASH(I)
               ICMOA = 1 + ICMOQ(I) + (NPSH(I)+NISH(I)) * NBASI
               CALL QGEMM(NBASI,NASHI,NASHI,D1,
     &              'N','N',IPQTOQ(1,0),
     &              WORK(KCMO+ICMOA-1),NBASI,NORBI,NZ,
     &              'N','N',IPQTOQ(1,0),
     &              UNO(JUNO,JUNO,1),NASHT,NASHT,NZ,
     &              D0,IPQTOQ(1,0),
     &              CMO(ICMOA),NBASI,NORBI,NZ)
C
               IF ( IPRINT .GE. 20 ) THEN
                  WRITE(LUPRI,'(/1X,A,I3)')
     &                 '(RGETNO) After transformation: '//
     &                 'MO coefficients, symmetry ',I
                  CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                 NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
 10         CONTINUE
C
            MOISNO = .TRUE.
            CALL MEMREL('RGETNO.transform',WORK,1,KFRSAV,KFREE,LFREE)
         ENDIF
      ENDIF
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RGETNO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmccno */
      SUBROUTINE KRMCCNO(KEYCNO,RTRACI,CREF,CMO,IPRINT,EIG,
     &                   WORK,KFREE,LFREE)
C***********************************************************************
C
C     Transform orbitals to natural and/or Fock type orbitals.
C
C     Input:
C        CREF - CI reference vector
C        CMO  - MO orbitals
C
C     Output:
C        CREF, CMO - transformed quantities.
C        EIG,  eigenvalues of blocked FD matrix (if TRTOFD true)
C              primarily for use in QCDHF
C
C     Written by J. Thyssen - Mar 4 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "mx2orb.h"
#include "dcbidx.h"
C
      DIMENSION CREF(*), CMO(*), EIG(*), WORK(*)
      CHARACTER KEYCNO*6
      LOGICAL RTRACI
C
      LOGICAL TRTOFD, TRTOCO, TRTONO
      DIMENSION J1MO(2), N1MO(2), J2MO(2), N2MO(2)
C
      CALL QENTER('KRMCCNO')
      KFRSAV = KFREE
C
C
      IF (KEYCNO .EQ. 'ONLYNO' .OR. KEYCNO .EQ. 'NATORB') THEN
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 4) WRITE(LUPRI,'(/1X,2A)')
     &           '(KRMCCNO) No transformation.'
            GOTO 2000
         END IF
         IF (IPROPT .GE. 4) WRITE(LUPRI,'(/1X,2A)')
     &        '(KRMCCNO) Transformation of active ',
     &        'orbitals to natural orbitals.'
         TRTOFD = .FALSE.
         TRTOCO = .FALSE.
         TRTONO = .TRUE.
      ELSE IF (KEYCNO .EQ. 'ONLYFD' .OR. KEYCNO .EQ. 'ONLYFC' .OR.
     &         KEYCNO .EQ. 'ONLYCO' .OR. KEYCNO .EQ. 'CANORB') THEN
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 4) WRITE(LUPRI,'(/1X,3A)')
     &           '(KRMCCNO) Transformation of ',
     &           'inactive/virtual/positronic ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .FALSE.
         ELSE
            IF (IPROPT .GE. 4) WRITE(LUPRI,'(/1X,2A)')
     &           '(KRMCCNO) Transformation of all ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .TRUE.
            TRTONO = .FALSE.
         END IF
      ELSE
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 4) WRITE(LUPRI,'(/1X,3A)')
     &           '(KRMCCNO) Transformation of ',
     &           'inactive/virtual/positronic ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .FALSE.
         ELSE
            IF (IPROPT .GE. 4) WRITE(LUPRI,'(/1X,4A)')
     &           '(KRMCCNO) Transformation of active ',
     &           'orbitals to natural orbitals and ',
     &           'inactive/virtual/positronic orbitals ',
     &           'to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .TRUE.
         END IF
      END IF
C
C
C     ***************************************
C     *** Transform to Fock type orbitals ***
C     ***************************************
C
C     Calculate one-electron density matrix
C     -------------------------------------
C
C     Needed for natural orbitals or FV.
C
      CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
      IF (NASHT .GE. 1) THEN
         CALL RMAKDV(WORK(KDV),CREF,WORK,KFREE,LFREE)
      END IF
C
      IF ( .NOT. (TRTOFD .OR. TRTOCO) ) GOTO 1000
C
      CALL MEMGET('REAL',KEIGVC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFD,N2ORBXQ,WORK,KFREE,LFREE)
      IF (NASHT .GT. 0) THEN
         LFV = N2ORBXQ
      ELSE
         LFV = 0
      END IF
      CALL MEMGET('REAL',KFV,LFV,WORK,KFREE,LFREE)
C
C
C     Calculate Fock matrices
C     -----------------------
C
      KFC = KFD
      CALL RFCKMAT(WORK(KDV),CMO,EIX,EAX,WORK(KFC),WORK(KFV),
     &     WORK(KFREE),LFREE)
C     CALL RFCKMAT(DV,CMO,EMCMY,EMCACT,FC,FV,WORK,LWORK)
C
C     Calculate FD
C
      IF (NASHT .GT. 0) CALL DAXPY(N2ORBXQ,D1,WORK(KFV),1,WORK(KFD),1)
C
      IF (TRTOCO) THEN
C
C        Extract active-active part of FD for RTRACI.
C
         CALL RGETAC(WORK(KFD),WORK(KDV),IPROPT)
C        CALL RGETAC(FD,FDAC,IPRINT)
C
         IF (.NOT. TRTOFD) THEN
            CALL MEMREL('RKRMCCNO.after RGETAC for FCAC',
     &                  WORK,1,KEIGVC,KFREE,LFREE)
            GO TO 1000
         END IF
      END IF
      CALL QUIT('KRMCCNO: TRTOFD not implemented correctly yet')
C     hjaaj Feb 2004 TODO: because QDIAG sorts eigenvalues code
C     must be modified to only call QDIAG with a block at a time
C     (p-p block, i-i block, s-s block) instead of zeroing inter-shell
C     couplings and active-active block.
C
C     Zero inter-shell coupling-blocks.
C     ---------------------------------
C
C     e-p and p-e:
C
      DO I = 1, NFSYM
         J1MO(I) = IORB(I) + NPSH(I) + 1
         N1MO(I) = NESH(I)
         J2MO(I) = IORB(I) + 1
         N2MO(I) = NPSH(I)
      END DO
C
      CALL BLKZERO(J1MO,N1MO,J2MO,N2MO,WORK(KFD),NORBT,NORBT,
     &     NZ,NFSYM,1)
      CALL BLKZERO(J2MO,N2MO,J1MO,N1MO,WORK(KFD),NORBT,NORBT,
     &     NZ,NFSYM,1)
C
C     i-{u,a}, {u,a}-i:
C
      DO I = 1, NFSYM
         J1MO(I) = IORB(I) + NPSH(I) + NISH(I) + 1
         N1MO(I) = NASH(I) + NSSH(I)
         J2MO(I) = IORB(I) + NPSH(I) + 1
         N2MO(I) = NISH(I)
      END DO
C
      CALL BLKZERO(J1MO,N1MO,J2MO,N2MO,WORK(KFD),NORBT,NORBT,
     &     NZ,NFSYM,1)
      CALL BLKZERO(J2MO,N2MO,J1MO,N1MO,WORK(KFD),NORBT,NORBT,
     &     NZ,NFSYM,1)
C
C     u-a, a-u:
C
      DO I = 1, NFSYM
         J1MO(I) = IORB(I) + NPSH(I) + NOCC(I) + 1
         N1MO(I) = NSSH(I)
         J2MO(I) = IORB(I) + NPSH(I) + NISH(I) + 1
         N2MO(I) = NASH(I)
      END DO
C
      CALL BLKZERO(J1MO,N1MO,J2MO,N2MO,WORK(KFD),NORBT,NORBT,
     &     NZ,NFSYM,1)
      CALL BLKZERO(J2MO,N2MO,J1MO,N1MO,WORK(KFD),NORBT,NORBT,
     &     NZ,NFSYM,1)
C
      IF (NASHT .GT. 0) THEN
C        Zero active-active block, this one is treated below in
C        RTRACI section /fix feb 2004 hjaaj
C        hjaaj: QUESTION / TODO? can we be sure
C               RSJACO/QDIAG returns a unit matrix as eigenvectors
C               for a zero matrix ???
         CALL BLKZERO(J2MO,N2MO,J2MO,N2MO,WORK(KFD),NORBT,NORBT,
     &        NZ,NFSYM,1)
      END IF
C
C
C     Diagonalize FD
C     --------------
C
C
      IF (NZ .EQ. 1) THEN
         IJOB = 1
         IORDER = 0
         IPACK = 0
         CALL RSJACO(NORBT,NORBT,NORBT,WORK(KFD),EIG,
     &        IJOB,IORDER,IPACK,WORK(KEIGVC),WORK(KFREE),LFREE)
      ELSE
         MATZ = 1
         CALL QDIAG(NZ,NORBT,WORK(KFD),NORBT,NORBT,
     &        EIG,MATZ,WORK(KEIGVC),NORBT,NORBT,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/,1X,2A,I4)')
     &           '*** ERROR in RGETNO ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in RGETNO ***')
         END IF
      END IF
C
      CALL MEMREL('RKRMCCNO.after FD diag',WORK,1,KFD,KFREE,LFREE)
C
C     Transform orbitals
C     ------------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
      DO I = 1, NFSYM
C
C        CMO(new) = CMO(old) * EIGVC
C
         NBASI = NFBAS(I,0)
         NORBI = NORB(I)
C
C        Transform all orbitals
C
         CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KCMO+ICMOQ(I)),NBASI,NORBI,NZ,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KEIGVC+I2ORBX(I,I)),NORBT,NORBT,NZ,
     &        D0,IPQTOQ(1,0),
     &        CMO(1+ICMOQ(I)),NBASI,NORBI,NZ)
C
C        Copy back active orbitals
C
         ICMOA = ICMOQ(I) + (NPSH(I)+NISH(I)) * NBASI
         CALL DCOPY(NBASI*NASH(I),WORK(KCMO+ICMOA),1,
     &        CMO(1+ICMOA),1)
C
      END DO
C
C     Deallocate everything but the density matrix.
C
      CALL MEMREL('RKRMCCNO.after FD rot',WORK,1,KEIGVC,KFREE,LFREE)
C
 1000 CONTINUE
C
C
C     **********************************************************
C     *** Transform to natural or canonical active  orbitals ***
C     **********************************************************
C
      IF (.NOT. (TRTONO .OR. TRTOCO) .OR. NZCONF.EQ.0) GOTO 2000
C
C
      call quit('*** ERROR in KRMCCNO ***' //
     &     'TRTOCO/TRTONO not implemeted yet!')
C
C
      CALL MEMGET('REAL',KEIG,NASHT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,N2ASHXQ,WORK,KFREE,LFREE)
C
C     Find transformation matrix
C     --------------------------
C
      IF (MCTYPE .EQ. JGAS) THEN
C
C        Calculate pseudo-natural or ``pseudo''-Fock type
C        transformation matrix.
C
C        CALL RGASNO(WORK(KDV),WORK(KEIGVC),WORK(KEIG))
C
         call quit('*** ERROR in KRMCCNO *** ' //
     &        'NO not implemented for GAS')
C
      ELSE
C
C        CAS wave function.
C
C        Find transformation matrix.
C        (WORK(KDV) is either active part of FD or DV).
C
         call quit('*** ERROR in KRMCCNO *** ' //
     &        'NO not implemented for CAS yet')
C        hjaaj TODO: at least NORBT must be NASHT in RSJACO/QDIAG calls
C        check for other errors ... /feb 2004
         IF (NZ .EQ. 1) THEN
            IJOB = 1
            IORDER = 0
            IPACK = 0
            CALL RSJACO(NORBT,NORBT,NORBT,WORK(KDV),WORK(KEIG),
     &           IJOB,IORDER,IPACK,WORK(KEIGVC),WORK(KFREE),LFREE)
         ELSE
            MATZ = 1
            CALL QDIAG(NZ,NORBT,WORK(KDV),NORBT,NORBT,
     &           WORK(KEIG),MATZ,WORK(KEIGVC),NORBT,NORBT,
     &           WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/,1X,2A,I4)')
     &              '*** ERROR in RGETNO ***: ',
     &              'QDIAG failed with error code ',IERR
               CALL QUIT('*** ERROR in RGETNO ***')
            END IF
         END IF
C
      END IF
C
C     Transform orbitals
C     ------------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
      DO I = 1, NFSYM
         IF (NASH(I) .NE. 0) THEN
C
C           CMO(NO) = CMO * UNO
C
            NBASI = NFBAS(I,0)
            NASHI = NASH(I)
            NORBI = NORB(I)
            ICMOA = 1 + ICMOQ(I) + (NPSH(I)+NISH(I)) * NBASI
            CALL QGEMM(NBASI,NASHI,NASHI,D1,
     &           'N','N',IPQTOQ(1,0),
     &           WORK(KCMO+ICMOA-1),NBASI,NORBI,NZ,
     &           'N','N',IPQTOQ(1,0),
     &           WORK(KEIGVC+I2ASHX(I,I)),NASHT,NASHT,NZ,
     &           D0,IPQTOQ(1,0),
     &           CMO(ICMOA),NBASI,NORBI,NZ)
C
         END IF
      END DO
C
C     Counter-rotate CI vector
C     ------------------------
C
C     CALL RTRACI(CREF,WORK(KEIGVC))
C
 2000 CONTINUE
      CALL QEXIT('KRMCCNO')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck canorb */
      SUBROUTINE CANORB(C,FD,U,WORK,LWORK)
C***********************************************************************
C
C     Diagonalize FD and return transformation matrix U.
C     C controls which part of FD we diagonalize.
C
C     Written by J. Thyssen - Mon Mar 19 16:37:10 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      CHARACTER C*1
      DIMENSION FD(NORBT,NORBT,NZ)
      DIMENSION U(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('CANORB')
C
C
      IF (C .EQ. 'E') THEN
         CALL MEMGET('REAL',KFD,NESHT*NESHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NESHT*NESHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NESHT,NESHT,IDXE2G,NZ)
         NDIM = NESHT
      ELSE IF (C .EQ. 'P') THEN
         CALL MEMGET('REAL',KFD,NPSHT*NPSHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NPSHT*NPSHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NPSHT,NPSHT,IDXP2G,NZ)
         NDIM = NPSHT
      ELSE IF (C .EQ. 'I') THEN
         CALL MEMGET('REAL',KFD,NISHT*NISHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NISHT*NISHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NISHT,NISHT,IDXI2G,NZ)
         NDIM = NISHT
      ELSE IF (C .EQ. 'A') THEN
         CALL MEMGET('REAL',KFD,NASHT*NASHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NASHT*NASHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NASHT,NASHT,IDXU2G,NZ)
         NDIM = NASHT
      ELSE IF (C .EQ. 'V') THEN
         CALL MEMGET('REAL',KFD,NSSHT*NSSHT*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFD),NSSHT*NSSHT*NZ)
         CALL MATGAT(FD,NORBT,NORBT,WORK(KFD),NSSHT,NSSHT,IDXV2G,NZ)
         NDIM = NSSHT
      ELSE
         WRITE(LUPRI,'(2A)')
     &        '*** ERROR in MATGAT *** Unknown control character: ',C
         CALL QUIT('*** ERROR in MATGAT *** Unknown control character')
      END IF
C
C
C     Diagonalize matrix
C     ------------------
C
      CALL MEMGET('REAL',KEIG,NDIM,WORK,KFREE,LFREE)
C
      IF (NZ .EQ. 1) THEN
         IJOB = 1
         IORDER = 0
         IPACK = 0
         CALL RSJACO(NDIM,NDIM,NDIM,WORK(KFD),WORK(KEIG),
     &        IJOB,IORDER,IPACK,U,WORK(KFREE),LFREE)
      ELSE
         MATZ = 1
         CALL QDIAG(NZ,NDIM,WORK(KFD),NDIM,NDIM,
     &        WORK(KEIG),MATZ,U,NDIM,NDIM,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/,1X,2A,I4)')
     &           '*** ERROR in CANORB ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in CANORB ***')
         END IF
      END IF
C
C     write(6,*) 'eigenvalues'
C     CALL OUTPUT(WORK(KEIG),1,NDIM,1,1,NDIM,1,1,6)
C
c     write(6,*) 'eigenvectors'
c     CALL PRQMAT(U,NDIM,NDIM,NDIM,NDIM,NZ,IPQTOQ(1,0),LUPRI)
C
      CALL MEMREL('CANORB',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('CANORB')
      RETURN
      END
