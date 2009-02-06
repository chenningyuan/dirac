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
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION CREF(*)
      DIMENSION OCCNO(*), UNO(NASHT,NASHT,*)
      DIMENSION CMO(*)
      DIMENSION WORK(*)
C
      LOGICAL NATORB, MOISNO, EX, FILEOPEN
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
      CALL MEMGET('REAL',KOCCNO2,NORBT  ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KOCCNO2),NORBT)
*
*     FIXME: right order of active orbitals if NZ == 2 and 
*            NO occ. numbers are calculated with LUCIAREL
*
      IF( NATOLCR .and. NZ .eq. 2 ) THEN
        WRITE(LUPRI,'(/A/)') 
     & ' *** WARNING *** QDIAG sorts eigenvalues - NO occupation'//
     & ' might appear in wrong order.'
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
C             OCCNO(JOCC) = WORK(KEIG + JACT - 1)
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
C        
C        save occupancies for later use 
         DO I = 1, NFSYM
           IOFF1 = IOCC(I) + 1
           IOFF2 = IORB(I) + 1
           IF(NOCC(I).gt.0)THEN
             CALL DCOPY(NOCC(I),OCCNO(IOFF1),1,WORK(KOCCNO2+IOFF2-1),1)
           END IF
         END DO
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

               CALL ORDER3(WORK(KFREE),OCCNO(IOCCI+1),NOCCI,NOCCI,0,-1)
               WRITE(LUPRI,9004) I, (OCCNO(IOCCI+J),J=1,NOCCI)
            END IF
         END DO
      END IF
C
 9001 FORMAT(/'   Symmetry',I3,//,'   No occupied orbitals')
 9002 FORMAT(/'   Symmetry',I3,//,(5F14.9))
 9003 FORMAT(/'   Sum =',T15,F14.9)
 9004 FORMAT(/'   Symmetry',I3,' ** REORDERED DESCENDING OCCUPATION **',
     &            //,(5F14.9)/)
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
C           write back to file with label MCCINATO
            LUNI = -1
            INQUIRE(FILE='KRMCSCF',EXIST=EX,OPENED=FILEOPEN,NUMBER=LUNI)
            IF(EX.and.FILEOPEN)THEN
              CALL WRTKRMC(LUNI,'MCCINATO',CMO,NCMOTQ)
              CALL WRTKRMC(LUNI,'MCNATOCC',WORK(KOCCNO2),NORBT)
            ELSE IF( EX.and.(.not. FILEOPEN))THEN
              LUNI = 99
              CALL OPNFIL(LUNI,'KRMCSCF','OLD    ','RGETNO')
              CALL WRTKRMC(LUNI,'MCCINATO',CMO,NCMOTQ)
              CALL WRTKRMC(LUNI,'MCNATOCC',WORK(KOCCNO2),NORBT)
              CLOSE(LUNI,STATUS='KEEP')
            ELSE IF(.NOT. EX)THEN
              LUNI = 99
              CALL OPNFIL(LUNI,'KRMCSCF','UNKNOWN','RGETNO')
              CALL NEWLAB('*RGETNO ',LUNI,LUPRI)
              CALL WRTKRMC(LUNI,'MCCINATO',CMO,NCMOTQ)
              CALL WRTKRMC(LUNI,'MCNATOCC',WORK(KOCCNO2),NORBT)
              CLOSE(LUNI,STATUS='KEEP')
            END IF
            WRITE(LUPRI,'(/A)')
     &   '  Natural orbitals saved with label MCCINATO on KRMCSCF'
C
C           copy back MO coefficients
            CALL DZERO(CMO,NCMOTQ)
            CALL DCOPY(NCMOTQ,WORK(KCMO),1,CMO,1)
            MOISNO = .TRUE.
         ENDIF
      ENDIF
      CALL MEMREL('RGETNO.transform',WORK,1,KFRSAV,KFREE,LFREE)
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
C                S. Knecht  - Jul 16 2008
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
      IF (KEYCNO .EQ. 'ONLYNO' .OR. KEYCNO .EQ. 'NATORB') THEN
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,2A)')
     &           '(KRMCCNO) No transformation.'
            GOTO 2000
         END IF
         IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,2A)')
     &        '(KRMCCNO) Transformation of active ',
     &        'orbitals to natural orbitals.'
         TRTOFD = .FALSE.
         TRTOCO = .FALSE.
         TRTONO = .TRUE.
      ELSE IF (KEYCNO .EQ. 'ONLYFD' .OR. KEYCNO .EQ. 'ONLYFC' .OR.
     &         KEYCNO .EQ. 'ONLYCO' .OR. KEYCNO .EQ. 'CANORB') THEN
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,3A)')
     &           '(KRMCCNO) Transformation of ',
     &           'positronic/inactive/virtual ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .FALSE.
         ELSE
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,2A)')
     &           '(KRMCCNO) Transformation of all ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .TRUE.
            TRTONO = .FALSE.
         END IF
      ELSE
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,3A)')
     &           '(KRMCCNO) Transformation of ',
     &           'positronic/inactive/virtual ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .FALSE.
         ELSE
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,4A)')
     &           '(KRMCCNO) Transformation of active ',
     &           'orbitals to natural orbitals and ',
     &           'positronic/inactive/virtual orbitals',
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
C     1. Calculate one-electron density matrix
C     ----------------------------------------
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
C     ... initialize
      CALL DZERO(WORK(KEIGVC),N2ORBXQ)
      CALL DZERO(WORK(KFD),   N2ORBXQ)
C
      IF (NASHT .GT. 0) THEN
         LFV = N2ORBXQ
      ELSE
         LFV = 0
      END IF
      CALL MEMGET('REAL',KFV,LFV,WORK,KFREE,LFREE)
C
C
C     2. Calculate Fock matrices
C     --------------------------
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
C
CSK   CALL QUIT('KRMCCNO: TRTOFD not implemented correctly yet')
C     hjaaj Feb 2004 TODO: because QDIAG sorts eigenvalues code
C     must be modified to only call QDIAG with a block at a time
C     (p-p block, i-i block, s-s block) instead of zeroing inter-shell
C     couplings and active-active block.
C
CSK   DONE! July 2008
C
      IF (IPROPT .GE. 5) THEN
         DO I = 1, NFSYM
            WRITE(LUPRI,*) 'RKRMCCNO: FOCK matrix for symmetry',I
            CALL PRQMAT(WORK(KFD+I2ORBX(I,I)),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         ENDDO
      END IF
C
C     3. Diagonalize Fock matrix 
C     --------------------------
C
C     3.1 p-p block
C     -------------
C
      IF ( NPSHT .gt. 0 ) THEN
C
C        ... allocate for p-p block diagonalization of Fock-matrix
         CALL MEMGET('REAL',KFDP,   N2PSHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVCP,N2PSHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGP,  NPSHT    ,WORK,KFREE,LFREE)
C        ... initialize
         CALL DZERO(WORK(KFDP),   N2PSHX*NZ)
         CALL DZERO(WORK(KEIGVCP),N2PSHX*NZ)
         CALL DZERO(WORK(KEIGP),  NPSHT)
C
C        extract p-p part of WORK(KFD)
C
         CALL MATGAT(WORK(KFD),NORBT,NORBT,WORK(KFDP),NPSHT,NPSHT,
     &               IDXP2G,NZ)
C
C        diagonalize
C
         CALL DIAGC('P',WORK(KFDP),WORK(KEIGP),WORK(KEIGVCP),
     &              NPSHT,WORK(KFREE),LFREE)
C
C        insert new p-p matrix block into complete transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCP),
     &               NPSHT,NPSHT,IDXP2G,NZ)
C
         CALL MEMREL('RKRMCCNO.after p-p',WORK,1,KFDP,KFREE,LFREE)
      END IF
C
C     3.2 i-i block
C     -------------
C
      IF ( NISHT .gt. 0 ) THEN
C
C        ... allocate for i-i block diagonalization of Fock-matrix
         CALL MEMGET('REAL',KFDI,   N2ISHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVCI,N2ISHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGI,  NISHT    ,WORK,KFREE,LFREE)
C        ... initialize
         CALL DZERO(WORK(KFDI),   N2ISHX*NZ)
         CALL DZERO(WORK(KEIGVCI),N2ISHX*NZ)
         CALL DZERO(WORK(KEIGI),  NISHT)
C
C        extract i-i part of WORK(KFD)
C
         CALL MATGAT(WORK(KFD),NORBT,NORBT,WORK(KFDI),NISHT,NISHT,
     &               IDXI2G,NZ)
C
C        diagonalize
C
         CALL DIAGC('I',WORK(KFDI),WORK(KEIGI),WORK(KEIGVCI),
     &              NISHT,WORK(KFREE),LFREE)
C
C        insert new i-i matrix block into complete transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCI),
     &               NISHT,NISHT,IDXI2G,NZ)
C
         CALL MEMREL('RKRMCCNO.after i-i',WORK,1,KFDI,KFREE,LFREE)
      END IF
C
C     3.3 'active' secondary orbitals: v-v block
C     -------------
C
      IF ( NVSHT .gt. 0 ) THEN
C
C        ... allocate for v-v block diagonalization of Fock-matrix
         CALL MEMGET('REAL',KFDVI,   N2VSHXQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVCVI,N2VSHXQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVI,  NVSHT  ,WORK,KFREE,LFREE)
C        ... initialize
         CALL DZERO(WORK(KFDVI),   N2VSHXQ)
         CALL DZERO(WORK(KEIGVCVI),N2VSHXQ)
         CALL DZERO(WORK(KEIGVI),  NVSHT  )
C
C        extract s-s part of WORK(KFD)
C
         CALL MATGAT(WORK(KFD),NORBT,NORBT,WORK(KFDVI),NVSHT,
     &               NVSHT,IDXV2G,NZ)
C
C        diagonalize
C
         CALL DIAGC('V',WORK(KFDVI),WORK(KEIGVI),WORK(KEIGVCVI),
     &              NVSHT,WORK(KFREE),LFREE)
C
C
C        insert new s-s matrix block into transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCVI),
     &               NVSHT,NVSHT,IDXV2G,NZ)
C
         CALL MEMREL('RKRMCCNO.after v-v',WORK,1,KFDVI,KFREE,LFREE)
      END IF
C
C     3.4 deleted secondary orbitals: d-d block 
C     -------------
C
      IF( NDSHT .gt. 0 )THEN
C        ... allocate
         CALL MEMGET('REAL',KEIGVCDL,N2DSHXQ,WORK,KFREE,LFREE)
C
C        ... initialize and make unit matrix
C
         CALL DUNIT2(WORK(KEIGVCDL),NDSHT,NDSHT,NDSHT,NZ)
C
C        insert d-d matrix block into complete transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCDL),
     &               NDSHT,NDSHT,IDXD2G,NZ)
C
         CALL MEMREL('RKRMCCNO.after d-d',WORK,1,KEIGVCDL,KFREE,LFREE)
      END IF
C
C     3.5 unit matrix for a-a block
C     -----------------------------
C
      IF (NASHT .GT. 0) THEN
C        ... allocate
         CALL MEMGET('REAL',KEIGVCA,N2ASHXQ,WORK,KFREE,LFREE)
C
C        ... initialize and make unit matrix
C
         CALL DUNIT2(WORK(KEIGVCA),NASHT,NASHT,NASHT,NZ)
C
C        insert a-a matrix block into complete transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCA),
     &               NASHT,NASHT,IDXU2G,NZ)
C
         CALL MEMREL('RKRMCCNO.after a-a',WORK,1,KEIGVCA,KFREE,LFREE)
      END IF
C
      CALL MEMREL('RKRMCCNO.after FD diag',WORK,1,KFD,KFREE,LFREE)
C
C     4. Transform orbitals
C     ---------------------
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
         IF( IPROPT .ge. 2 )THEN
           WRITE(LUPRI,'(A,I3)') ' RKRMCCNO: '//
     &     'Coefficients before transformation for symmetry',I
           CALL PRQMAT(CMO(1+ICMOQ(I)),NBASI,NORBI,
     &          NBASI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
           WRITE(LUPRI,'(A,I3)')
     &     ' RKRMCCNO: transformation matrix for symmetry',I
           CALL PRQMAT(WORK(KEIGVC+I2ORBX(I,I)),NORBI,NORBI,
     &          NORBI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
         END IF
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
         IF( NDSH(I) .gt. 0 )THEN
C
C          Copy back deleted orbitals from s-s block
C
           ICMOD = ICMOQ(I) + 
     &             (NPSH(I)+NISH(I)+NASH(I)+NVSH(I)) * NBASI
           CALL DCOPY(NBASI*NDSH(I),WORK(KCMO+ICMOD),1,
     &          CMO(1+ICMOD),1)
         END IF
C
         IF( IPROPT .ge. 2 )THEN
           WRITE(LUPRI,'(A,I3)') ' RKRMCCNO: '//
     &     'MO coefficients after FD transformation for symmetry',I
           CALL PRQMAT(CMO(1+ICMOQ(I)),NBASI,NORBI,
     &          NBASI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END DO
C
C     Deallocate everything but the density matrix.
C
      CALL MEMREL('RKRMCCNO.after FD rot',WORK,1,KEIGVC,KFREE,LFREE)
C
 1000 CONTINUE
C
      IF (.NOT. (TRTONO .OR. TRTOCO) .OR. NZCONF.EQ.0) GOTO 2000
C
C     **********************************************************
C     *** Transform to natural or canonical active  orbitals ***
C     **********************************************************
C
C
      CALL Abend2('*** ERROR in KRMCCNO ***' //
     &            'TRTOCO/TRTONO not implemeted yet!')
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
C
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
C/* Deck RPROP */
      SUBROUTINE KRCI_PROP(CMO,WORK,LWORK)
C***********************************************************************
C
C     1. calculate property MO representation for list of (symmetric) 
C        one-electron (and later also two-electron) operators.
C
C     2. calculate properties , e.g., transition dipole moments. 
C
C     3. obtain final quantities , e.g., oscillator strength f. 
C
C     Written by S. Knecht - Aug 2008
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
#include "dcbopt.h"
#include "krciprop.h"
#include "dcbkrci.h"
C
      DIMENSION CMO(*), WORK(*)
C
#include "memint.h"
      MZ = MIN(2,NZ)
C
C     introduce property task
C
      CALL TITLER('KR-CI property calculation ','*',199)
C
      WRITE(LUPRI,9000)
 9000 FORMAT(
     &     3X,'This is output from the DIRAC KR-CI property module',/,
     &     //,
     &     3X,'Written by: ',/,
     &     3X,'  Stefan Knecht and Hans Joergen Aa. Jensen ',//,
     &     1X,79('*'))
C
      WRITE (LUPRI,'(/A/A)')
     &   '  Symmetry   # of electronic states',
     &   '  ________   ______________________'
C               1              4
      NPROP_ROOTS_KRCI = 0
      DO I = 1, NKRCI_MAX_SYM
         NPROP_ROOTS_KRCI = NPROP_ROOTS_KRCI + NKRCI_CIROOTS(I)
         WRITE(LUPRI,'(I7,I15)') I, NKRCI_CIROOTS(I)
      END DO
         WRITE(LUPRI,'(A/A,I15/)')
     &   '  _________________________________',
     &   '  Total', NPROP_ROOTS_KRCI
C
C     eigenstate <--> symmetry relation
      IXROOTS     = 0
      IXROOTS_SYM = 0
      CALL IZERO(ISYMEIG_KRCI,NPROP_ROOTS_KRCI)
      DO I = 1, NKRCI_MAX_SYM
         IXROOTS_SYM = I
         DO J = 1, NKRCI_CIROOTS(I)
            IXROOTS = IXROOTS + 1
            ISYMEIG_KRCI(IXROOTS) = IXROOTS_SYM
         END DO
      END DO
C
C     KRCI_CVECS.x vector file extensions
      CALL DETXFLAB(XSYMFLAB,NKRCI_MAX_SYM)
C
C     ***************************************************
C     *** Memory allocation *****************************
C     ***************************************************
C
      LXPRPKRCI = MZ * NPROP_KRCI * NPROP_ROOTS_KRCI**2
      CALL MEMGET('REAL',KXPRPKRCI,LXPRPKRCI,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPRPCORE, NPROP_KRCI,WORK,KFREE,LFREE)
C
C     ... initialize
      CALL DZERO(WORK(KXPRPKRCI),LXPRPKRCI)
      CALL DZERO(WORK(KPRPCORE),NPROP_KRCI)
C
C     ***************************************************
C     *** operator(s) in MO representation **************
C     ***************************************************
C
C     active-active part in MOLFDIR format saved on 
C     file KRMC_FOCK
C
      CALL KRCIPRPMO(CMO,WORK(KPRPCORE),WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** calculate properties **************************
C     ***************************************************
C
      CALL KRCIGETPRP(WORK(KXPRPKRCI),WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** get spectroscopic quantities ******************
C     ***************************************************
C
C     array of eigenvalues 
      CALL MEMGET('REAL',KXROOTS,NPROP_ROOTS_KRCI,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KXROOTS),NPROP_ROOTS_KRCI)
C
      CALL KRCIPRPANA(WORK(KXPRPKRCI),WORK(KXROOTS),WORK(KPRPCORE),
     &                WORK(KFREE),LFREE)
C
      CALL MEMREL('RPROP_KRCI',WORK,KWORK,KWORK,KFREE,LFREE)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck KRCIPRPMO */
      SUBROUTINE KRCIPRPMO(CMO,PRPCOREE,WORK,LWORK)
C***********************************************************************
C
C     Transform property integrals to MO basis.
C     Integrals are saved on file KRMC_FOCK.
C
C     Input:
C        CMO     - orbital coefficients.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C write-protection for DUMMY.
#include "dummy.h"
#include "consts.h"
#include "maxorb.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "krciprop.h"
#include "dcbxpr.h"
      DIMENSION CMO(*), PRPCOREE(*), WORK(*)
      CHARACTER*8 FILELAB
      LOGICAL LOPEN
C
#include "memint.h"
C
      MZ = MIN(2,NZ)
C
      CALL QENTER('KRCIPRPMO')
C
      INQUIRE(FILE='AOPROPER',OPENED=LOPEN)
      IF (.NOT. LOPEN) THEN
         OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &        FILE='AOPROPER')
      ENDIF
C
C     ***************************************************
C     *** Get property matrices *************************
C     ***************************************************
C
      DO I = 1, NPROP_KRCI
C
         INDXPR  = LPROP_KRCI(I)
         ISYM    = IPRPSYM(INDXPR)
         ITRSYM  = IPRPTIM(INDXPR)
         IOPSY   = JBTOF(ISYM-1,1)
         FILELAB = PRPNAM(INDXPR)(1:8)
C
         CALL MEMGET('REAL',KRPRPMO,N2ORBXQ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KRPRPMO),N2ORBXQ)
         CALL PRPMAT(INDXPR,IOPSY,WORK(KRPRPMO),CMO,
     &               DUMMY,ICMOQ,NORB,WORK,KFREE,LFREE,IPROPT)
C
C        print full property matrix
C
         IF( IPROPT .ge. 5 ) THEN
            WRITE(LUPRI,'(/2A)')'   KRCIPRPMO: MO matrix of property: ',
     &                              FILELAB
            CALL PRQMAT(WORK(KRPRPMO),NORBT,NORBT,NORBT,
     &                  NORBT,NZ,IPQTOQ(1,ISYM-1),LUPRI )
         ENDIF
C
C        calculate PRPCORE(I) = 2.0D0 * SUM_(all inactive) P_ii
C
         IF( NISHT .gt. 0)THEN
            TEMP_PCORE  = SUM_DIAG_ELM('I',WORK(KRPRPMO),NORBT,NORBT)
            PRPCOREE(I) = D2 * TEMP_PCORE
         ELSE
            PRPCOREE(I) = D0
         END IF
         IF (IPROPT .GE. 0)
     &      WRITE(LUPRI,9000) FILELAB,PRPCOREE(I)
 9000 FORMAT(' (KRCIPRPMO): ',A,' property core value: ',F25.15/)
C
C        get active-active part
C
         CALL MEMGET('REAL',KRPRPMOAC,   N2ASHXQ,WORK,KFREE,LFREE)
CSK      IPROPT = 50
         CALL RGETAC(WORK(KRPRPMO),WORK(KRPRPMOAC),IPROPT)
C
C        transform to MOLFDIR format
         CALL MEMGET('REAL',KRPRPMOACM,4*N2ASHX*MZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KRPRPMOACM),  4*N2ASHX*MZ)
C        WRITE(LUPRI,'(/2A)')'   KRCIPRPMO: MO matrix of property: ',
C    &                           FILELAB
C        WRITE(LUPRI,'(A,3I5)')
C    &   ' - KRCIPRPMO: Boson, Fermion, and Time reversal symmetry :',
C    &   ISYM,IOPSY,ITRSYM
         CALL QFC2MFC(WORK(KRPRPMOAC),WORK(KRPRPMOACM),ISYM,
     &                ITRSYM,IPROPT)
C             QFC2MFC(DQFC,DMFC,ISYM,ITIM,IPRINT)
C
C        save on file KRMC_FOCK
         CALL KRCI_PRPFILE(LUKRM3,FILELAB,WORK(KRPRPMOACM),
     &                     4*N2ASHX*MZ,0)
         CALL MEMREL('KRCIPRPMO',WORK,KWORK,KWORK,KFREE,LFREE)
      END DO
C
      CLOSE(LU1INT,STATUS = 'KEEP')
      CALL QEXIT('KRCIPRPMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck KRCIPRPANA */
      SUBROUTINE KRCIPRPANA(XPROPMAT,EVECS,EPROP_CORE,WORK,LWORK)
C***********************************************************************
C
C     Calculate requested spectroscopic quantities.
C
C     implemented: 1. oscillator strength for electronic transitions
C                  2. omega expectation values (for j_z = s_z + l_z)
C
C     Input: property matrix XPROPMAT
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( XTHRESH = 1.0D-7 )
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "krciprop.h"
#include "dcbkrci.h"
#include "dcbxpr.h"
#include "debye.h"
      DOUBLE COMPLEX XMOM
      REAL*8 OSCL
      DIMENSION XPROPMAT(NPROP_ROOTS_KRCI,NPROP_ROOTS_KRCI,NPROP_KRCI,*)
      DIMENSION WORK(*), OSCL(3), EVECS(*), EPROP_CORE(*)
      INTEGER IND_DIP(3)
      CHARACTER XLABEL*16
C
#include "memint.h"
      CALL QENTER('XKRCIPRPANA')
C
      MZ = MIN(NZ,2)
C
C     1. (transition dipole moments) and oscillator strengths
C     =======================================================
C
      IF( DOOSCILLST )THEN
C
C       get energies from file KRCI_CVECS.INFO
        CALL GET_KRCICVECS_E(EVECS,NKRCI_CIROOTS,XSYMFLAB,NKRCI_MAX_SYM)
C
C       calculate and print oscillator strength - length gauge
C       ------------------------------------------------------
C
C       determine index for (x,y,z)
        XLABEL     = 'X dipole length'
        IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Y dipole length'
        IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Z dipole length'
        IND_DIP(3) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XREAL = 0.0D0
        XCOMP = 0.0D0
C
C       Check sum rule for oscillator strength
        SUML = 0.0D0
        WRITE(LUPRI,'(/A/A)')
     & '    initial state             final state   '//
     & '      oscillator strength (length gauge)',
     & '  __________________       _________________'//
     & '      __________________________________'
        DO 200 IVEC=1,NPROP_ROOTS_KRCI
          DO 100 JVEC= IVEC+1, NPROP_ROOTS_KRCI
            IF( EVECS(IVEC) .lt. EVECS(JVEC) )THEN
              INIST = IVEC
              JFLST = JVEC
            ELSE IF( EVECS(IVEC) .eq. EVECS(JVEC) )THEN
              GOTO 100
            ELSE IF( EVECS(IVEC) .gt. EVECS(JVEC) )THEN
              INIST = JVEC
              JFLST = IVEC
            END IF
            DEIG  = EVECS(JFLST) - EVECS(INIST)
C           WRITE(LUPRI,*) ' INIST, JFLST, DEIG',INIST, JFLST, DEIG
            DO IPRP=1,3
C             length gauge
              XREAL = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),1)
              IF( MZ .eq. 2 ) 
     &        XCOMP = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),2)
              XMOM = DCMPLX(XREAL,XCOMP)
C             WRITE(LUPRI,*) ' XREAL,XCOMP, ABS(XMOM)**2 ',
C    &                         XREAL,XCOMP, ABS(XMOM)**2
C             f_(x,y,z)  = 2.0D0*|<inist|(x,y,z)|jflst>|**2 * (E_fl - E_ini)
              OSCL(IPRP) = 2.0D0*ABS(XMOM)**2*(DEIG)
            END DO
            SUML = SUML + (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0
            ISYMST_IN = ISYMEIG_KRCI(INIST)
            ISYMST_FL = ISYMEIG_KRCI(JFLST)
C
C           print only values above a given threshold - 1.0D-7
            IF( (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0 .lt. XTHRESH ) GOTO 100
            WRITE(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F10.6)')
     &  '   * state',INIST,'(',XREPEIG(ISYMST_IN),')        * state',
     &      JFLST,'(',XREPEIG(ISYMST_FL),')                   ',
     &      (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0 
            WRITE(LUPRI,'(A,1F16.8,A,1F16.8)')
     &  '  E(au)',EVECS(INIST),'    E(au)',EVECS(JFLST)
 100      CONTINUE
 200    CONTINUE
        WRITE(LUPRI,*) 
     & ' -----------------------------------------------------------'//
     & '----------'
        WRITE(LUPRI,'(1A,1F10.6)')
     $     '   Sum of oscillator strength : ',SUML
C
#ifdef CODE_NOT_WORKING
C
C       determine index for (dx,dy,dz)
        XLABEL     = 'X dip vel'
        IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Y dip vel'
        IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Z dip vel'
        IND_DIP(3) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XREAL = 0.0D0
        XCOMP = 0.0D0
C
C       Check sum rule for oscillator strength
        SUMV = 0.0D0
        WRITE(LUPRI,'(/A/A)')
     & '    initial state             final state   '//
     & '      oscillator strength (velocity gauge)',
     & '  __________________       _________________'//
     & '      ____________________________________'
        DO 400 IVEC=1,NPROP_ROOTS_KRCI
          DO 300 JVEC= IVEC+1, NPROP_ROOTS_KRCI
            IF( EVECS(IVEC) .lt. EVECS(JVEC) )THEN
              INIST = IVEC
              JFLST = JVEC
            ELSE IF( EVECS(IVEC) .eq. EVECS(JVEC) )THEN
              GOTO 300
            ELSE IF( EVECS(IVEC) .gt. EVECS(JVEC) )THEN
              INIST = JVEC
              JFLST = IVEC
            END IF
            DEIG  = EVECS(JFLST) - EVECS(INIST)
            DO IPRP=1,3
C             velocity gauge
              XREAL = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),1)
              IF( MZ .eq. 2 ) 
     &        XCOMP = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),2)
              XMOM = DCMPLX(XREAL,XCOMP)
C
C             f_(dx,dy,dz)  = 2.0D0*|<inist|(x,y,z)|jflst>|**2 * (E_fl - E_ini)
              OSCL(IPRP) = 2.0D0*ABS(XMOM)**2*(DEIG)
            END DO
            SUMV = SUMV + (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0
            ISYMST_IN = ISYMEIG_KRCI(INIST)
            ISYMST_FL = ISYMEIG_KRCI(JFLST)
        WRITE(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F10.6)')
     & '   * state',INIST,'(',XREPEIG(ISYMST_IN),')        * state',
     &       JFLST,'(',XREPEIG(ISYMST_FL),')                     ',
     &       (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0 
        WRITE(LUPRI,'(A,1F16.8,A,1F16.8)')
     & '  E(au)',EVECS(INIST),'    E(au)',EVECS(JFLST)
 300      CONTINUE
 400    CONTINUE
        WRITE(LUPRI,*) 
     & ' -----------------------------------------------------------'//
     & '------------'
        WRITE(LUPRI,'(1A,1F10.6/)')
     &     '   Sum of oscillator strength : ',SUMV
        WRITE(LUPRI,'(1A,1F10.6/)')
     &     '   Total sum of oscillator strengths (length + velocity'//
     &     ' gauge)  : ',SUMV+SUML
#endif
      END IF
C
C     2. omega expectation values
C     ===========================
C
      IF( DOJZEXP )THEN
C       ... allocate omega matrix
        CALL MEMGET('REAL',KOMEGAEXP,NZ*NPROP_ROOTS_KRCI**2,WORK,
     &              KFREE,LFREE)
C       ... initialize
        CALL DZERO(WORK(KOMEGAEXP),NZ*NPROP_ROOTS_KRCI**2)
C
C       get and diagonalize omega expectation value matrix
        CALL GET_OMEGA_EXP(WORK(KOMEGAEXP),NPROP_ROOTS_KRCI,NZ,
     &                     XPROPMAT,NPROP_KRCI,XREPEIG,
     &                     MAX_NKRCI_MAX_SYM,NKRCI_CIROOTS,
     &                     WORK(KFREE),LFREE)
C       release memory        
        CALL MEMREL('XKRCIPRPANA',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
C 
      CALL QEXIT('XKRCIPRPANA')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck GET_OMEGA_EXP */
      SUBROUTINE GET_OMEGA_EXP(OMEGAMAT,NOMEGAST,NZ,XPROPMAT,NPROPX,
     &                         XREPOM,IALL_SYM,NKRCI_CIROOTS,
     &                         WORK,LWORK)
C***********************************************************************
C
C     Diagonalize matrix containing j_z = s_z + l_z expectation values.
C
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "krciprop.h"
#include "dcbxpr.h"
      CHARACTER*4 XREPOM(*)
      DIMENSION OMEGAMAT(NOMEGAST,NOMEGAST,*), NKRCI_CIROOTS(*)
      DIMENSION XPROPMAT(NOMEGAST,NOMEGAST,NPROPX,*), WORK(*)
      INTEGER IND_DIP(2)
      CHARACTER XLABEL*16
#include "memint.h"
      CALL QENTER('GET_OMEGA_EXP')
      NXTEST = 0
C
C     determine index for sigma_z and l_z
      XLABEL     = 'z-spin'
      IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROPX,XLABEL)
      XLABEL     = 'z-angular mom'
      IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROPX,XLABEL)
C
C     copy sz-matrix part and lz-matrix part to OMEGAMAT
      DO IZ = 1, NZ
        DO I = 1, NOMEGAST
C         ... extract sz and copy to OMEGAMAT
          CALL DCOPY(NOMEGAST,XPROPMAT(1,I,IND_DIP(1),IZ),1,
     &               OMEGAMAT(1,I,IZ),1)
C         ... now extract and add lz 
          CALL DAXPY(NOMEGAST,1.0D0,XPROPMAT(1,I,IND_DIP(2),IZ),1,
     &               OMEGAMAT(1,I,IZ),1)
        END DO
      END DO
C
C     debug print
      NXTEST = 0
      IF( NXTEST .gt. 0 )THEN
C       complete XPROPMAT
        WRITE(LUPRI,*) ' complete XPROPMAT'
        CALL WRTMATMN(XPROPMAT,1,LXPRPKRCI,1,LXPRPKRCI,LUPRI)
C       complete OMEGAMAT
        WRITE(LUPRI,*) ' complete OMEGAMAT'
        CALL WRTMATMN(OMEGAMAT,1,NZ*NOMEGAST**2,1,NZ*NOMEGAST**2,LUPRI)
      END IF
C
C     print omega expectation values
      WRITE(LUPRI,'(/A/,A)')
     & '    omega expectation values for calculated eigenstates',
     & '  ________________________________________________________'
      WRITE(LUPRI,'(/A/,A/)')
     & '    state #        symmetry        <omega>',
     & '  ----------     -----------   ---------------'
      IX_AC = 0
C
C     loop over all symmetries and block diagonalize
      DO 100 IX = 1, IALL_SYM
C
C       only symmetries with non-zero number of eigenstates
        IF( NKRCI_CIROOTS(IX) .le. 0 ) GOTO 100
C
        NXSTATES = NKRCI_CIROOTS(IX)
C       ... allocate for x-x block diagonalization of <OMEGA>-matrix
        CALL MEMGET('REAL',KOMX,   NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGVCX,NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGX,  NXSTATES      ,WORK,KFREE,LFREE)
C       ... initialize
        CALL DZERO(WORK(KOMX),   NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGVCX),NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGX),  NXSTATES)
C
C       extract x-x part of <OMEGA>-matrix
        CALL MATGATX(OMEGAMAT,NOMEGAST,WORK(KOMX),NXSTATES,IX_AC,NZ)
C
#ifdef MOD_DEBUG
        WRITE(LUPRI,*) ' OMEGAMAT (block) before diagonalization'
        CALL WRTMATMN(WORK(KOMX),1,NZ*NXSTATES**2,1,
     &                NZ*NXSTATES**2,LUPRI)
#endif
C
C       diagonalize
        CALL DIAGC('X',WORK(KOMX),WORK(KEIGX),WORK(KEIGVCX),
     &             NXSTATES,WORK(KFREE),LFREE)
C
        NXTEST = 0
        IF( NXTEST .gt. 0 )THEN
          WRITE(LUPRI,*) ' OMEGAMAT after diagonalization'
          CALL WRTMATMN(WORK(KOMX),1,NZ*NXSTATES**2,1,
     &                  NZ*NXSTATES**2,LUPRI)
          WRITE(LUPRI,*) ' OMEGAEIGV after diagonalization'
          CALL WRTMATMN(WORK(KEIGX),1,NXSTATES,1,NXSTATES,LUPRI)
          WRITE(LUPRI,*) ' OMEGAEIGVEC after diagonalization'
          CALL WRTMATMN(WORK(KEIGVCX),1,NZ*NXSTATES**2,1,
     &                  NZ*NXSTATES**2,LUPRI)
        END IF
        NXTEST = 0
C
C       in the complex case "NZ==2" the eigenvalues are sorted in
C       ascending order - find absolute largest element in each
C       eigenvector to identify the corresponding CI eigenstate.
        IF( NZ .ge. 2 )THEN
C
#ifdef CODE_NOT_WORKING
C         FIXME: this "reordering" code does not work properly ...
C
C         list of CI eigenstates
          CALL MEMGET('REAL',KEIGTX,  NXSTATES      ,WORK,KFREE,LFREE)
          CALL DZERO(WORK(KEIGTX),NXSTATES)
          CALL GET_XOMEGALIST(WORK(KEIGTX),WORK(KEIGVCX),NXSTATES,
     &                        XREPOM(IX),WORK(KEIGX))
#endif
          WRITE(LUPRI,'(/A/)') 
     &    ' ***** OMEGA VALUES ARE SORTED in ascending order *****'
          WRITE(LUPRI,'(A)') 
     &    ' To assign the states properly, rerun the property part'//
     &    ' with less number of roots'
          WRITE(LUPRI,'(A/)') 
     &    ' and do thus a stepwise identification... ;)'
C
C         print section
          DO I = 1, NXSTATES
            WRITE(LUPRI,'(I11,8X,A4,6X,1F16.10)')
     &            I,XREPOM(IX),WORK(KEIGX+I-1)
          END DO
        ELSE
C
C         print section
          DO I = 1, NXSTATES
          WRITE(LUPRI,'(I11,8X,A4,6X,1F16.10)')
     &                I,XREPOM(IX),WORK(KEIGX+I-1)
          END DO
        END IF
C
C       release memory
        CALL MEMREL('RKRMCCNO.after x-x',WORK,1,KOMX,KFREE,LFREE)
C
C       keep track of active states
        IX_AC = IX_AC + NXSTATES
 100  CONTINUE
C
C     ... done
      CALL QEXIT('GET_OMEGA_EXP')
      END
C----------------------------------------------------------------------
C     end of file krmccan.F
C----------------------------------------------------------------------
