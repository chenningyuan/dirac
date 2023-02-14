#define ESR_DEBUG -1
!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

C VERSION : $Revision: 11757 $
C DATE    : $Date: 2010-10-10 09:06:37 +0200 (Sun, 10 Oct 2010) $
C FILE    : pamesr.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ESRINP */
      SUBROUTINE ESRINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for ESR properties
C
C     Written by A.Noerager - Mar 2000
C     Last revision: hjaaj Oct 2022
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from .h files:
C  dcbprp.h : ALLATM, NPATOM, IATOM(:), ...
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcborb.h"
#include "dgroup.h"

      PARAMETER (NTABLE = 20)
      LOGICAL SET, NEWDEF, RESET, GORG_BY_USER, g_10ppm
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      REAL*8    WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'xXXXXXX','.PRINT ','.OPERAT','.MULTIP','.SPINPO',
     &            '.SDCI  ','.MAX CI','.CI ROO','.THRCVE','.THRCI ',
     &            '.SDTCI ','.SDT*CI','.ESR CI','.MY CI ','.THRPCI',
     &            '.KR-CON','.MAXMK2','.G10PPM','.STATE ','.SELECT'/
      DATA SET/.FALSE./
C
      NEWDEF = (WORD .EQ. '*ESR   ')
      IF (SET) THEN
         IF (NEWDEF)
     &      CALL QUIT('Only one "*ESR" input section allowed')
C        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C
C     Initialize /CBIESR/
C     ===================
C
      IPRESR    = 0
      NESRP     = 0
      MULTIPESR = 0
      N_CIESR   = -1
      MAX_CIESR_IT = 20
      THR_CVEC_ESR = 1.0D-3  ! trial vector threshold
      THRPCI_ESR   = 1.0D-3  ! print of CI vector threshold
      THR_CIESR    = -1.0D0  ! default depends on g_10ppm true or false
      STATE_ESR    = 1       ! first CI root belonging to multiplet

C     variables for spin polarization and correlation in CI
      NGAS_CIESR(1)  = 0  ! RAS1_GAS
      LVL_CIESR(:,1) = 0
      NGAS_CIESR(2)  = 1  ! one RAS2-GAS subspace in this version (NOPEN = 1)
      LVL_CIESR(:,2) = -1   ! no occupation restriction in RAS2-GAS
      NGAS_CIESR(3)  = 0  ! RAS3_GAS
      LVL_CIESR(:,3) = 0
      MAXMK2_ESR     = 99999

      NGSH_CIESR(:,:,:) = 0

      USE_KRAMERS_CONJ = .FALSE.
C
C     Process input
C     =========================
C
      NEWDEF = (WORD .EQ. '*ESR   ')
      g_10ppm = .FALSE. ! 1 ppt accuracy by default
      ICHANG = 0
      INPERR = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                      11,12,13,14,15,16,17,18,19,20), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized for *ESR'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword for *ESR')
    1          CONTINUE !
               GO TO 100
    2          CONTINUE ! .PRINT
C&&&& PRINT:  Print level in ESR module
                  READ(LUCMD,*) IPRESR
               GO TO 100
C&&&& OPERATOR: Define additional operator for results section
    3          CONTINUE ! .OPERATOR
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        IPRESR)
                  CALL OP1IND('ESRINP',INDESR,LESRP,NESRP,INDXPR,
     &                        MAXESR)
                  ICHANG = ICHANG + 1
               GO TO 100
C Read multiplicity
    4          CONTINUE ! .MULTIPLICITY
                  READ(LUCMD,*) MULTIPESR
               GO TO 100
C&&&& SPINPOLARIZATION (= SCI)
    5          CONTINUE ! .SPINPOL
                  NGAS_CIESR(1)  = 1
                  LVL_CIESR(1,1) = 1
                  READ(LUCMD,*) (NGSH_CIESR(I,1,1),I=1,NFSYM) ! valence hole orbitals
                  NGAS_CIESR(3)  = 1
                  LVL_CIESR(1,3) = 1
                  READ(LUCMD,*) (NGSH_CIESR(I,1,3),I=1,NFSYM) ! electron orbitals
               GO TO 100
C&&&& Singles and doubles CI
    6          CONTINUE ! .SDCI
                  NGAS_CIESR(1)  = 1
                  LVL_CIESR(1,1) = 2
                  READ(LUCMD,*) (NGSH_CIESR(I,1,1),I=1,NFSYM) ! valence hole orbitals
                  NGAS_CIESR(3)  = 1
                  LVL_CIESR(1,3) = 2
                  READ(LUCMD,*) (NGSH_CIESR(I,1,3),I=1,NFSYM) ! electron orbitals
               GO TO 100
C&&&& MAX CI ITERATIONS
    7          CONTINUE ! .MAX CI ITERATIONS
                  ! negative value is code for doing full CI
                  ! zero value is code for generate list of determinants but no CI
                  READ(LUCMD,*) MAX_CIESR_IT
               GO TO 100
C&&&& CI ROOTS to calculate
    8          CONTINUE ! .CI ROOTS
                  READ(LUCMD,*) N_CIESR
               GO TO 100
C&&&& THRESHOLD FOR CVECS IN GASCIP_RCISTD
    9          CONTINUE  ! .THRCVEC
                  READ(LUCMD,*) THR_CVEC_ESR
               GO TO 100
C&&&& Convergence threshold for CI
   10          CONTINUE  ! .THRCI
                  READ(LUCMD,*) THR_CIESR
               GO TO 100
C&&&& Singles doubles triples CI
   11          CONTINUE ! .SDTCI
                  NGAS_CIESR(1)  = 1
                  LVL_CIESR(1,1) = 3
                  READ(LUCMD,*) (NGSH_CIESR(I,1,1),I=1,NFSYM) ! valence hole orbitals
                  NGAS_CIESR(3)  = 1
                  LVL_CIESR(1,3) = 3
                  READ(LUCMD,*) (NGSH_CIESR(I,1,3),I=1,NFSYM) ! electron orbitals
               GO TO 100
C&&&& MR-CI-SDT/SD
   12          CONTINUE ! .SDT*CI
                  NGAS_CIESR(1)  = 1
                  LVL_CIESR(1,1) = 3
                  READ(LUCMD,*) (NGSH_CIESR(I,1,1),I=1,NFSYM) ! valence hole orbitals
                  NGAS_CIESR(3)  = 1
                  LVL_CIESR(1,3) = 2
                  READ(LUCMD,*) (NGSH_CIESR(I,1,3),I=1,NFSYM) ! electron orbitals
               GO TO 100
C&&&& general specification of ESR CI space
   13          CONTINUE ! .ESR CI
                  READ(LUCMD,*) NGAS_CIESR(1), NGAS_CIESR(3)
                  NERR = 0
                  IF (NGAS_CIESR(1) .LT. 0 .OR.
     &                NGAS_CIESR(1) .GT. maxgas_esr) NERR = NERR + 1
                  IF (NGAS_CIESR(3) .LT. 0 .OR.
     &                NGAS_CIESR(3) .GT. maxgas_esr) NERR = NERR + 1
                  IF (NERR .GT. 0) THEN
                     WRITE(LUPRI,'(//A,I0//A,I0/A,I0)')
     &               ' ".ESR CI" input error, number of GAS spaces'//
     &               ' not in interval 0:',maxgas_esr,
     &         ' Number of RAS1-GAS spaces from input: ',NGAS_CIESR(1),
     &         ' Number of RAS3-GAS spaces from input: ',NGAS_CIESR(3)
                     CALL QUIT('".ESR CI" input error')
                  END IF
                  ! RAS1-GAS subspaces:
                  DO J = NGAS_CIESR(1),1,-1
                     READ(LUCMD,*) LVL_CIESR(J,1) ! max number of holes in this RAS1-GAS space
                     READ(LUCMD,*) (NGSH_CIESR(I,J,1),I=1,NFSYM) ! hole orbitals
                  END DO
                  ! RAS3-GAS subspaces:
                  DO J = 1,NGAS_CIESR(3)
                     READ(LUCMD,*) LVL_CIESR(J,3) ! max number of electrons in this RAS3-GAS space
                     READ(LUCMD,*) (NGSH_CIESR(I,J,3),I=1,NFSYM) ! electron orbitals
                  END DO
               GO TO 100
C&&&& flexible input of the CI you want
   14          CONTINUE ! .MY CI
                  NGAS_CIESR(1) = 2
                  NGAS_CIESR(3) = 1
                  READ(LUCMD,*) LVL_CIESR(2,1),LVL_CIESR(1,1),
     &                          LVL_CIESR(1,3)
                  READ(LUCMD,*) (NGSH_CIESR(I,2,1),I=1,NFSYM) ! core hole orbitals
                  READ(LUCMD,*) (NGSH_CIESR(I,1,1),I=1,NFSYM) ! valence hole orbitals
                  READ(LUCMD,*) (NGSH_CIESR(I,1,3),I=1,NFSYM) ! electron orbitals
               GO TO 100
C&&&& limit for printing CI coefficients in output
   15          CONTINUE ! .THRPCI
                  READ(LUCMD,*) THRPCI_ESR
               GO TO 100
C&&&& use Kramers conjugation on CI vector
   16          CONTINUE ! .KR-CON
                  USE_KRAMERS_CONJ = .TRUE.
               GO TO 100
C&&&& limit MK2 value
   17          CONTINUE ! .MAXMK2
                  READ(LUCMD,*) MAXMK2_ESR
C&&&& calculate g tensor with 10 ppm accuracy
   18          CONTINUE ! .G10PPM
                  g_10ppm = .TRUE.
               GO TO 100
C&&&& first CI root belonging to multiplet
   19          CONTINUE ! .STATE
                  READ(LUCMD,*) STATE_ESR
               GO TO 100
C&&&& select atoms for HFC calculations
   20          CONTINUE ! .SELECT
                  IF (.NOT. ALLATM) THEN
                     WRITE(LUPRI,'(//A/)') 'INPUT ERROR,'//
     &               '.SELECT specified both under **INTEGRALS and *ESR'
                     INPERR = INPERR + 1
                  END IF
                  READ (LUCMD, *) NPATOM
                  IF (NPATOM .GT. MXCENT) THEN
                     WRITE (LUPRI,'(/A,/A,I6,/A,I6)')
     &             ' Too many atoms selected with .SELECT.',
     &             ' Number of atoms selected: ',NPATOM,
     &             ' Number of atoms allowed:  ',MXCENT
                     CALL QUIT('Error in *ESR input')
                  END IF
                  READ (LUCMD, *) (IPATOM(I),I=1,NPATOM)
                  ALLATM = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized for *ESR'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt for *ESR')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(/A,I3)')
     & '* Solving EPR/ESR problem for an effective spin multiplicity of'
     & ,MULTIPESR
      IF (MULTIPESR .LT. 1) THEN
         WRITE(LUPRI,*) 'INPUT ERROR: .MULTIP not (correctly) specified'
         INPERR = INPERR + 1
      END IF
      IF (STATE_ESR .NE. 1) THEN
      WRITE(LUPRI,'(/A,I0)')
     & '* Lowest CI root belonging to multiplet is CI root no.  '
     & ,STATE_ESR
      END IF
      NGAS_CIESR_T = NGAS_CIESR(1) + NGAS_CIESR(3)
      IF (NGAS_CIESR_T .gt. 0 .OR. MULTIPESR .ne. 2) THEN
      !  NGAS_CIESR_T test because o.k. with no CI if just doublet one open-shell DHF or DFT
         IF (N_CIESR < 0) N_CIESR = MULTIPESR + (STATE_ESR-1)
         WRITE(LUPRI,'(/A,I8)')
     & '* Number of roots to converge in CI:', N_CIESR
         IF (N_CIESR .lt. MULTIPESR) THEN
            WRITE(LUPRI,*) 'INPUT ERROR: .MULTIP > .CI ROOTS',
     &         MULTIPESR, N_CIESR
            INPERR = INPERR + 1
         END IF
         IF (THR_CIESR .LT. 0.0D0) THEN
            IF (g_10ppm) THEN
               THR_CIESR = 1.D-6
               WRITE(LUPRI,'(A)')
     & '* g-tensor is solved for 10 ppm accuracy.'
            ELSE
               THR_CIESR = 1.D-4
               WRITE(LUPRI,'(A)')
     & '* g-tensor is solved for 1 ppt accuracy (default).'
            END IF
         ELSE IF (g_10ppm) THEN
            WRITE(LUPRI,'(A)')
     & '* INFO: both .THRCI and .G10PPM specified, smallest value used.'
            THR_CIESR = MIN(1.D-6, THR_CIESR)
         END IF
         IF (MAX_CIESR_IT .gt. 0) THEN
            WRITE(LUPRI,'(A,I8/A,1P,D8.1)')
     & '* Max number of CI iterations      :', MAX_CIESR_IT,
     & '* CI convergence threshold         :', THR_CIESR
         ELSE IF (MAX_CIESR_IT .eq. 0) THEN
            WRITE(LUPRI,'(A)')
     & '* Generate determinant list but no CI iterations.'
         ELSE
            WRITE(LUPRI,'(A)')
     & '* Diagonalize CI Hamiltonian matrix completely.'
         END IF
      END IF
      IF (NGAS_CIESR_T .GT. 0) THEN
         WRITE (LUPRI,'(/A/)')
     &      ' * Spin polarization and correlation CI specification:'

         DO J = NGAS_CIESR(1),1,-1
             WRITE(LUPRI,'(A,I2,A,T50,2I5)')
     &       ' Max',LVL_CIESR(J,1),' holes in RAS1-GAS orbitals: ',
     &       NGSH_CIESR(1:NFSYM,J,1)
         END DO

         WRITE (LUPRI,'(/A,T50,2I5/)')
     &      '  RAS2-GAS active/open shell orbitals:', NASH(1:NFSYM)

         DO J = 1,NGAS_CIESR(3)
             WRITE(LUPRI,'(A,I2,A,T50,2I5)')
     &       ' Max',LVL_CIESR(J,3),' electrons in RAS3-GAS orbitals: ',
     &       NGSH_CIESR(1:NFSYM,J,3)
         END DO

         IF (MAXMK2_ESR .NE. 99999) THEN
             WRITE(LUPRI,'(/A,I0)')
     &       ' Max MK2 value for included Slater determinants: ',
     &       MAXMK2_ESR
         END IF

         NERR = 0
         LVL_RAS1 = 10000
         LVL_RAS3 = 10000
         DO J = 1, maxgas_esr
            IF (LVL_CIESR(J,1) > LVL_RAS1) NERR = NERR + 100
            LVL_RAS1 = LVL_CIESR(J,1)
            IF (LVL_CIESR(J,3) > LVL_RAS3) NERR = NERR + 1000
            LVL_RAS3 = LVL_CIESR(J,3)

            DO I = 1,NFSYM
               IF (NGSH_CIESR(I,J,1) < 0) NERR = NERR + 1
               IF (NGSH_CIESR(I,J,3) < 0) NERR = NERR + 10
            END DO
         END DO
         IF (MAXMK2_ESR .LT. MULTIPESR-1) NERR = NERR + 10000

         IF (NERR > 0) THEN
           WRITE(LUPRI,'(//A,I5/)')
     &        '*ESR: ERROR in specification of CI. Error code',NERR
           WRITE(LUPRI,*) 'NGAS_CIESR    :',NGAS_CIESR(1:3)
           WRITE(LUPRI,*) 'RAS1 LVL_CIESR:',LVL_CIESR(1:maxgas_esr,1)
           WRITE(LUPRI,*) 'RAS3 LVL_CIESR:',LVL_CIESR(1:maxgas_esr,3)
           WRITE(LUPRI,*) 'RAS1 orbitals :',NGSH_CIESR(1,1:maxgas_esr,1)
           WRITE(LUPRI,*) 'RAS3 orbitals :',NGSH_CIESR(1,1:maxgas_esr,3)
           WRITE(LUPRI,*) 'MAXMK2 < MULTIPESR-1?',MAXMK2_ESR,MULTIPESR-1
         END IF
         IF (NFSYM .GT. 1) THEN
            CALL QUIT('Sorry, inversion symmetry is not implemented'//
     &         ' in ESR CI module.')
         END IF
         IF (NERR > 0) CALL QUIT('*ESR: ERROR in specification of CI')
      END IF

      CALL PRSYMB(LUPRI,'=',75,0)

  999 CONTINUE
      IF (INPERR.GT.0) CALL QUIT('Input error(s) in *ESR')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesr */
      SUBROUTINE PRPESR(WF,ESRVAL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     First version by Anette Noerager - Mar 2000
C     Last revision: March 2022, Hans Joergen Aa. Jensen
C
C     Calculate expectation values defined in dcbesr.h(/CBIESR/)
C
C     Input:
C        WF : wave function type (text)
C
C     Output:
C        ESRVAL : the complex quasi-degenerate expectation values
C                 for requested NESRP properties
C
C     Note that the ESR property matrix ESRVAL is complex,
C     also when NZ=1.
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "mxcent.h"
#include "consts.h"
C     for PRQMAT print of complex matrices:
      integer, parameter ::NZCMPLX = 2, IQINDX(2) = (/ 1,2 /)
      integer, parameter ::max_MULTIPESR = 9 ! corresponding to max S_eff = 4
C
      REAL*8    ESRVAL(MULTIPESR,MULTIPESR,NZCMPLX,*),WORK(*)
      CHARACTER WF*4, PPM_OR_PPT*3
      LOGICAL   FILEX, SWITCH_XY
      REAL*8    GTENSOR(3,3),GCORR(3,3),GMATRIX(3,3),GEIG(3),GSHIFT(3),
     &          GEVEC(3,3), GCOR(3)
      INTEGER   JMS(-4:4)
      REAL*8    PHASE_FAC(MULTIPESR,2) ! automatic array
C
C Used from common blocks:
C  dcbgen.h : LUCOEF
C  dcbesr.h : MULTIPESR, ...
C
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbesr.h"
#include "dcbprp.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbopt.h"
C
      CALL QENTER('PRPESR')
!d    write (lupri,*) 'DEBUG 1 NZ, NZ_in_CI',NZ,NZ_in_CI
C
      IPRESR = MAX(IPRESR, ESR_DEBUG )
      KFRSAV = KFREE
!d    print *,'PRPESR  KFREE,LFREE',KFREE,LFREE
!d    call memchk('prpesr 1',work,1)
C
      IF (WF(1:3) .NE. 'DHF' .AND. WF(1:3) .NE. 'DFT') THEN
         WRITE (LUPRI,'(/2A/A)')
     &   ' WARNING: ESR properties are only'//
     &   ' implemented for open shell DHF; not for ',WF,
     &   ' WARNING: .ESR request is ignored !'
         GO TO 9999
      END IF

!d    write (lupri,*) 'hj1 N_CIESR MAX_CIESR_IT MULTIPESR NZCONF',
!d   &                     N_CIESR,MAX_CIESR_IT,MULTIPESR,NZCONF
C
C ** Some parameter checking
C
      IF (N_CIESR .GE. 0) THEN ! not DFT nor HF
         !
         IF (MULTIPESR .gt. NZCONF) THEN
            CALL QUIT('.MULTIPESR specification > # determinants')
         END IF
         IF (MAX_CIESR_IT .EQ. 0) THEN
            CALL HEADER(
     & 'INFO: Calculation stops now as requested with ".MAX CI" of 0',1)
            GO TO 9000
         ELSE IF (MAX_CIESR_IT .GT. 0) THEN
            ! Davidson iterations
            N_CIESR = MAX(N_CIESR,MULTIPESR)
            N_CIESR = MIN(N_CIESR,NZCONF)
!        write(lupri,*) 'N_CIESR,MULTIPESR is', N_CIESR,MULTIPESR
         ELSE
            ! full CI (now primarily for debugging Davidson iterations code)
            N_CIESR = NZCONF
         END IF
      END IF
!d    write (lupri,*) 'hj2 N_CIESR MAX_CIESR_IT MULTIPESR NZCONF',
!d   &                     N_CIESR,MAX_CIESR_ITMULTIPESR,,NZCONF
C
C **  Retrieve molecular orbitals
C
      CALL MEMGET2('REAL','CMO',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL ESR_READ_CMO(WORK(KCMO),IPRESR,WORK,KFREE,LFREE)
C
C **  Solve ESR CI problem
C     on exit: WORK(KEIGVL) - N_CIESR CI eigenvalues
C              WORK(KEIGVC) - N_CIESR CI eigenvectors (NZCONF,N_CIESR,NZ_in_CI)
C
      CALL ESR_RESOLVE(WORK(KCMO),IPRESR,N_CIESR,MAX_CIESR_IT,
     &                 THR_CIESR,WORK,KEIGVL,KEIGVC,KFREE,LFREE)

      call memchk('after esr_resolve',work,1)
C     CALL ESR_RESOLVE(CMO,IPRINT,N_CIROOTS,MAX_CI_IT,THR_CI,
C     &           WORK,K_EIGVAL,K_EIGVEC,KFREE,LFREE)

C
C     Make CI Hamiltonian matrix (diagonal now, but it may be made
C     non-diagonal in ESRPHASE)
C
      CALL MEMGET2('REAL','H_CI',KHCI,
     &   NZ_in_CI*MULTIPESR*MULTIPESR,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHCI),NZ_in_CI*MULTIPESR*MULTIPESR)
      DO I = STATE_ESR, STATE_ESR - 1 + MULTIPESR
         II = KHCI + (I-1)*MULTIPESR + (I-1)
         WORK(II) = WORK(KEIGVL-1 + I)
      END DO

      IF (IPRESR .GT. 1) THEN
         WRITE(lupri,*) 'H_CI after ESR_RESOLVE:'
         CALL OUTPUT(WORK(KHCI),1,MULTIPESR,1,MULTIPESR,
     &               MULTIPESR,MULTIPESR,-1,LUPRI)
      END IF
C
C **  Find state vectors spin Hamiltonian ordered after effective S_z,
C     and fix relative phases of these state vectors.
C
C     Memory allocation for ESRPHASE and PRPESR1
C
      CALL MEMGET2('REAL','PMAT',KPMAT,NZ*N2ASHX ,WORK,KFREE,LFREE)
C
#if ESR_DEBUG > 20
       WRITE(LUPRI,*)
     &   'PRPESR: Real part of the N_CIESR CI vectors'
       CALL OUTPUT(WORK(KEIGVC),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       IF (NZ_in_CI .EQ. 2) THEN
          WRITE(LUPRI,*)
     &      'PRPESR: Imag part of the N_CIESR CI vectors'
          KEIGVC_I = KEIGVC + NZCONF*N_CIESR
          CALL OUTPUT(WORK(KEIGVC_I),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       END IF
#endif
      CALL ESRPHASE(ESRVAL,PHASE_FAC,SWITCH_XY,WORK(KCMO), WORK(KPMAT),
     &              WORK(KHCI),WORK(KEIGVC),WORK,KFREE,LFREE)
C     CALL ESRPHASE(ESRVAL,PHASE_factors,SWITCH_XY,CMO,PMAT,
C    &              HCI,CIEVEC,WORK,KFREE,LFREE)
C
C     Memory allocation for PRPESR1
C
      CALL MEMGET2('LOGI','ESRST',KESRST ,2*NESRP   ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','PMAT ',KPMAT  ,NZ*N2ASHX ,WORK,KFREE,LFREE)

      IF (IPRESR .GT. 1) THEN
        WRITE(LUPRI,'(/A)')
     &    '* PRPESR called with the following operators:'
        DO I = 1,NESRP
          INDXPR = LESRP(I)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A)')
     &       'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      END IF
C
      CALL PRPESR1(ESRVAL,PHASE_FAC,WORK(KESRST),WORK(KCMO),
     &             WORK(KPMAT),WORK(KEIGVC),WORK,KFREE,LFREE)
C     CALL PRPESR1(ESRVAL,PHASE_factors,LESRST,CMO,PMAT,
C    &             CIEVEC,WORK,KFREE,LFREE)
C
      IF (IPRESR .GT. 1) THEN
        CALL HEADER('PRPESR results in CI eigenvector basis'//
     &           ' (test output)',1)
        DO I = 1,NESRP
          INDXPR = LESRP(I)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A)')
     &       'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
          CALL PRQMAT(ESRVAL(1,1,1,I),N_CIESR,N_CIESR,
     &                N_CIESR,N_CIESR,NZCMPLX,IQINDX,LUPRI)
        ENDDO
      END IF
C
      IF (ESRGTENS .OR. ESR_HFCC) THEN
         CALL TITLER('@  DIRAC  - ESR PROPERTIES','*',116)
         WRITE(LUPRI,'(A,I5)')
     &   '@  *  Effective spin multiplicity:',MULTIPESR
      END IF
      IF (ESRGTENS) THEN
C
C     ESRVAL corresponds now to B * g in state vector basis
C     The factors in front of ESRVAL below will cancel
C     the corresponding factors in the Seff_p matrices
C     and multiply with 2 (because H_Zeeman = 1/2 * B * g * Seff )
C
      IF (MULTIPESR .EQ. 2) THEN
C     ....doublet molecules (S_eff = 1/2)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
        jmhlf = 2
        jphlf = 1
        DO K = 1,3
C       ... for B field in direction k (x,y,z)
          I = IPESR(MXCOOR+K,1)
          GTENSOR(K,1) =  4.0D0*ESRVAL(jphlf,jmhlf,1,I)
          GTENSOR(K,2) = -4.0D0*ESRVAL(jphlf,jmhlf,2,I)
          GTENSOR(K,3) =  4.0D0*ESRVAL(jphlf,jphlf,1,I)
          J = IPESR(MXCOOR+K,2)
          GCORR(K,1)   =  4.0D0*ESRVAL(jphlf,jmhlf,1,J)
          GCORR(K,2)   = -4.0D0*ESRVAL(jphlf,jmhlf,2,J)
          GCORR(K,3)   =  4.0D0*ESRVAL(jphlf,jphlf,1,J)
        ENDDO
C
      ELSE IF (MULTIPESR .EQ. 3) THEN
C     .... triplet molecules (S_eff = 1)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
C QUESTION: would it not be better to average jp1-j0 and j0-jm1 values ???
C
        jm1 = 3
        j0  = 2
        jp1 = 1
        D2SQRT2 = 2.0D0*SQRT(2.0D0)
        DO K = 1,3
C       ... for B field in direction k (x,y,z)
          I = IPESR(MXCOOR+K,1)
          GTENSOR(K,1) =  D2SQRT2*ESRVAL(j0 ,jm1,1,I)
          GTENSOR(K,2) = -D2SQRT2*ESRVAL(j0 ,jm1,2,I)
          GTENSOR(K,3) =    2.0D0*ESRVAL(jp1,jp1,1,I)
c??       GTENSOR(K,3) = ESRVAL(jp1,jp1,1,I) - ESRVAL(jm1,jm1,1,I)
c         WRITE(lupri,*)'K =', K
c         WRITE(lupri,*)'GTENSOR(K,1)', GTENSOR(K,1),
c    &                  'GTENSOR(K,2)', GTENSOR(K,2),
c    &                  'GTENSOR(K,3)', GTENSOR(K,3)
          J = IPESR(MXCOOR+K,2)
          GCORR(K,1)   =  D2SQRT2*ESRVAL(j0 ,jm1,1,J)
          GCORR(K,2)   = -D2SQRT2*ESRVAL(j0 ,jm1,2,J)
          GCORR(K,3)   =    2.0D0*ESRVAL(jp1,jp1,1,J)
c         WRITE(lupri,*)'K =', K
c         WRITE(lupri,*)'GCORR(K,1)', GCORR(K,1),
c    &      'GCORR(K,2)', GCORR(K,2),
c    &      'GCORR(K,3)', GCORR(K,3)
        ENDDO
C
C     ***************************************
C     Zero-field-splitting parameters D and E
C     ***************************************
C
        EP1 = WORK(KHCI)
        E0  = WORK(KHCI+  MULTIPESR+1)
        EM1 = WORK(KHCI+2*MULTIPESR+2)
        ZFS_D = -(E0-(EP1+EM1)/2.0D0)
        ZFS_E = (EP1-EM1)/2.0D0
        ZFS_DCM = ZFS_D*XTKAYS
        ZFS_ECM = ZFS_E*XTKAYS

        CALL HEADER(' Zero field splitting parameters : ',1)
        WRITE (LUPRI,'(2(A,F25.8,A/))')
     &   '  ZFS parameter D:      ',ZFS_DCM, ' cm-1 ',
     &   '  ZFS parameter E:      ',ZFS_ECM, ' cm-1 '

      ELSE IF (MULTIPESR .EQ. 4) THEN
C     .... quartet molecules (S_eff = 3/2)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
C       Take mean value of M_Seff=3/2 and 1/2 values
        fac3hlf = 2.0d0/sqrt(3.0d0) ! fac1hlf is 1.0d0
C
        jm3hlf = 4
        jm1hlf = 3
        jp1hlf = 2
        jp3hlf = 1
        DO K = 1,3
C       ... for B field in direction k (x,y,z)
          I = IPESR(MXCOOR+K,1)
!         GTENSOR(K,1) =  2.0D0*ESRVAL(jp1hlf,jm1hlf,1,I)
!         GTENSOR(K,2) = -2.0D0*ESRVAL(jp1hlf,jm1hlf,2,I)
!         GTENSOR(K,3) =  4.0D0*ESRVAL(jp1hlf,jp1hlf,1,I)
          GTENSOR(K,1) =  ESRVAL(jp3hlf,jp1hlf,1,I) * fac3hlf
     &                 +  ESRVAL(jp1hlf,jm1hlf,1,I)
          GTENSOR(K,2) = -ESRVAL(jp3hlf,jp1hlf,2,I) * fac3hlf
     &                 -  ESRVAL(jp1hlf,jm1hlf,2,I)
          GTENSOR(K,3) =  2.0D0*(ESRVAL(jp3hlf,jp3hlf,1,I)/3.0D0
     &                 +         ESRVAL(jp1hlf,jp1hlf,1,I) )

          J = IPESR(MXCOOR+K,2)
!         GCORR(K,1)   =  2.0D0*ESRVAL(jp1hlf,jm1hlf,1,J)
!         GCORR(K,2)   = -2.0D0*ESRVAL(jp1hlf,jm1hlf,2,J)
!         GCORR(K,3)   =  4.0D0*ESRVAL(jp1hlf,jp1hlf,1,J)
          GCORR(K,1)   =  ESRVAL(jp3hlf,jp1hlf,1,J) * fac3hlf
     &                 +  ESRVAL(jp1hlf,jm1hlf,1,J)
          GCORR(K,2)   = -ESRVAL(jp3hlf,jp1hlf,2,J) * fac3hlf
     &                 -  ESRVAL(jp1hlf,jm1hlf,2,J)
          GCORR(K,3)   =  2.0D0*(ESRVAL(jp3hlf,jp3hlf,1,J)/3.0D0
     &                 +         ESRVAL(jp1hlf,jp1hlf,1,J) )
        ENDDO
C
C     ***************************************
C     Zero-field-splitting
C     ***************************************
C
        Ep3 = WORK(KHCI)
        Ep1 = WORK(KHCI+  MULTIPESR+1)
        Em1 = WORK(KHCI+2*MULTIPESR+2)
        Em3 = WORK(KHCI+3*MULTIPESR+3)
        ZFS_D = 0.5D0*(Ep3+Em3-Ep1-Em1)
        ZFS_DCM = ZFS_D*XTKAYS

        CALL HEADER(' Zero field splitting parameter : ',1)
        WRITE (LUPRI,'(A,F25.8,A/)')
     &   '  E(Sz_eff=3/2) - E(Sz_eff=1/2):  ',ZFS_DCM, ' cm-1 '
C
      ELSE
         call quit('Only doublet, triplet, and quartet implemented yet')
      END IF

      IF (SWITCH_XY) THEN
         DO I = 1,3
            TMP = GTENSOR(1,I)
            GTENSOR(1,I) = GTENSOR(2,I)
            GTENSOR(2,I) = TMP
            TMP = GCORR(1,I)
            GCORR(1,I) = GCORR(2,I)
            GCORR(2,I) = TMP
         END DO
      END IF
C
      CALL HEADER(' ESR g shifts : ',1)
C
C     Tests indicated that convergence noise in g-tensors is approximately +/-10*THR_CIESR
C
      IF (THR_CIESR .LT. 9.D-5) THEN
         PPM_OR_PPT = 'ppm'
         FAC_GSHIFT = 1.D6
         WRITE(LUPRI,'(/A,F5.1,A)')
     &   '* g-shifts converged to approximately +/-',
     &   (THR_CIESR/1.D-7),' ppm' !
      ELSE
         PPM_OR_PPT = 'ppt'
         FAC_GSHIFT = 1.D3
         WRITE(LUPRI,'(/A,F5.1,A)')
     &   '* g-shifts converged to approximately +/-',
     &   (THR_CIESR/1.D-4),' ppt' !
      END IF
C
      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   '* The g-tensor for H(Zeeman) = 1/2c B g Seff :',
     &   ((GTENSOR(I,K),K=1,3),I=1,3)

C     Get G-matrix = gtensor*gtensor(transposed)
C
      DO I = 1,3
        DO K = 1,3
          GMATRIX(I,K) = GTENSOR(I,1)*GTENSOR(K,1)
     &                 + GTENSOR(I,2)*GTENSOR(K,2)
     &                 + GTENSOR(I,3)*GTENSOR(K,3)
        ENDDO
      ENDDO


      IF (IPRESR .GE. 1) WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   '  The G-matrix = g * g^T :', ((GMATRIX(I,K),K=1,3),I=1,3)

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - 2.0D0)*FAC_GSHIFT
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / 3.0D0
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / 3.0D0
      WRITE (LUPRI,
     &   '(/A//A//5X,4F18.8//3A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The Dirac equation g values corresponding to g_e = 2',
     &   '  The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   '  The g shifts in ',PPM_OR_PPT,' (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO,
     &   ' The corresponding eigenvectors:',
     &   ((GEVEC(I,K),K=1,3),I=1,3)
C
C     (2) Get the corrected g values by scaling
C
      DO I = 1,3
        GCOR(I) = (GEIG(I)/2)*GFAC
        GSHIFT(I) = (GCOR(I) - GFAC)*FAC_GSHIFT
      ENDDO
      GISO = (GCOR(1) + GCOR(2) + GCOR(3)) / 3.0D0
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / 3.0D0
      WRITE (LUPRI,'(/A,F13.10,A//A//5X,4F18.8//3A//5X,4F18.1))')
     &   '* The corrected g values scaled with',GFAC,' / 2',
     &   '  The corrected g values (x, y, z, iso):',
     &   (GCOR(I),I =1,3), GISO,
     &   '  The corrected g shifts in ',PPM_OR_PPT,' (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO
C
C     (3) Get the corrected g values by adding correction term
C
      DO I = 1,3
      DO K = 1,3
C        ... for B field in direction k (x,y,z)
         GTENSOR(K,I) = GTENSOR(K,I) + GCORR(K,I)
      ENDDO
      ENDDO

      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   '* The g_1 corrected g-tensor for H(Zeeman) = 1/2c B g Seff :',
     &   ((GTENSOR(I,K),K=1,3),I=1,3)

C     Get G-matrix = gtensor*gtensor(transposed)
C
      DO I = 1,3
        DO K = 1,3
          GMATRIX(I,K) = GTENSOR(I,1)*GTENSOR(K,1)
     &                 + GTENSOR(I,2)*GTENSOR(K,2)
     &                 + GTENSOR(I,3)*GTENSOR(K,3)
        ENDDO
      ENDDO


      IF (IPRESR .GE. 1) WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   '  The G-matrix = g * g^T :', ((GMATRIX(I,K),K=1,3),I=1,3)

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - GFAC)*FAC_GSHIFT
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / 3.0D0
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / 3.0D0
      WRITE (LUPRI,
     &   '(/A,F13.10/A//5X,4F18.8//3A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The g_1 corrected Dirac eq. g values corresponding to g_e ='
     &   ,GFAC,
     &   '  The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   '  The g shifts in ',PPM_OR_PPT,' (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO,
     &   '  The corresponding eigenvectors:',
     &   ((GEVEC(I,K),K=1,3),I=1,3)

      call memchk('prpesr 4',work,1)
      ENDIF
C     .. end if (ESRGTENS) then
C
C     Print ESR hyperfine coupling constants from hyperfine Hamiltonian
C
      IF (ESR_HFCC)
     &   CALL PRP_ESR_HFC(ESRVAL,SWITCH_XY,WORK,KFREE,LFREE)
C
C     Memory deallocation
C
 9000 CALL MEMREL('PRPESR',WORK,1,KFRSAV,KFREE,LFREE)
 9999 CALL QEXIT('PRPESR')
      RETURN
C
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPESR: Coefficient file not found !'
      CALL QUIT('PRPESR: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ESRPHASE(ESRVAL,PHASE_factors,SWITCH_XY,CMO,PMAT,
     &                    HCI,CIEVEC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     --- 1) Find state vectors diagonalizing S_z (for a "standard" S_eff basis)
C         2) Modify phases acc. to convention for <M_S | S_x | M_S - 1 >
C         and  <M_S | S_y | M_S - 1 >.
C
C     Input :
C        CMO, HCI, CIEVEC
C
C     Output:
C        PHASE_factors : phase factors for each CI vector in S_eff basis, i.e. in S_z
C                        order, such that S_x and S_y have the conventional structure
C                       (if this is not ensured, some g values and HFCC values
C                        might have incorrect sign or be off-diagonal compared to the standard
C                        spin Hamiltonians based on relativistic perturbations to
C                        to non-relativistic eigenfunctions with exact M_S values.
C                        Wrong signs also affect the calculation of isotropic values.)
C                        NOTE: we do not multiply the CI eigenvectors with this factor,
C                        because this makes it possible to still use real CI eigenvectors
C                        when NZ = 4.
C        CIEVEC : CI eigenvectors transformed to S_eff basis
C        ESRVAL(:,:,2,1:4) : S_x, S_y, S_z, S^2 CI property matrices
C
C      Last revision March 2022 by Hans Joergen Aa. Jensen,
C      based on derivations by Jacob Kaldau.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C     for PRQMAT print of complex matrices:
      integer, parameter ::NZCMPLX = 2, IQINDX(2) = (/ 1,2 /)
C
      PARAMETER( THR_ZERO=1.D-12 )
C
      LOGICAL   SWITCH_XY, REAL_PMLF
      REAL*8    ESRVAL(MULTIPESR,MULTIPESR,2,4),
     &          PHASE_factors(MULTIPESR,2),
     &          CMO(*),PMAT(*),HCI(MULTIPESR,MULTIPESR,NZ_in_CI),
     &          CIEVEC(NZCONF,N_CIESR,NZ_in_CI),
     &          WORK(*)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   IPCMO(2)
C
C Used from common blocks:
C  dcbgen : LU1INT,?
C  dcbopt : KZCONF,NZCONF
C  dgroup : NZ,NZ_in_CI
C  dcbesr : ...,STATE_ESR,...
C
#include "dcbgen.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbesr.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
C
      CALL QENTER('ESRPHASE')
      KFRSAV = KFREE
C
C     Initialize
C     ==========
C     (using first NPSIGMA elements of ESRVAL for IPSIGMA properties,
C     later reused for the NESRP requested properties)
C
      IF (N_CIESR .LT. MULTIPESR+(STATE_ESR-1)) THEN
         WRITE(LUPRI,*) 'ERROR: N_CIESR.lt.MULTIPESR',N_CIESR,MULTIPESR
         CALL QUIT('N_CIESR .lt. MULTIPESR')
      END IF
      IF (NPSIGMA .NE. 3)      CALL QUIT('NPSIGMA .ne. 3')
      N2MULT  = MULTIPESR ** 2
      CALL DZERO(ESRVAL,N2MULT*2*NPSIGMA)
C
C     Step 1: Diagonalize S_z operator and transform to M_S(eff) : -M_S(eff) CI basis
C     ======
C
      I = IPSIGMA(3) ! S_z operator index
      INDXPR = LESRP(I)

#if ESR_DEBUG > 20
       WRITE(LUPRI,*)
     &   'ESRPHASE: Real part of the N_CIESR CI vectors'
       CALL OUTPUT(CIEVEC(1,1,1),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       IF (NZ_in_CI .EQ. 2) THEN
          WRITE(LUPRI,*)
     &      'ESRPHASE: Imag part of the N_CIESR CI vectors'
          CALL OUTPUT(CIEVEC(1,1,2),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       END IF
#endif

      PHASE_factors(:,1) = 1.0D0 ! no phase changes in ESRTRANS yet
      PHASE_factors(:,2) = 0.0D0
      CALL ESRTRANS(INDXPR,ESRVAL,PHASE_factors,
     &              CMO,PMAT,CIEVEC,WORK,KFREE,LFREE)
C       <i|S_z|j> is now in ESRVAL(...,1)
C       temp. use ESRVAL(...,2:3) for eigenvectors and eigenvalues
C       note: using NZ_in_CI in QDIAG and not 2 because we know Sigma_z
C             is a real operator
      CALL QDIAG(NZ_in_CI,MULTIPESR,ESRVAL,MULTIPESR,MULTIPESR,
     &           ESRVAL(1,1,1,3),1,ESRVAL(1,1,1,2),
     &           MULTIPESR,MULTIPESR,WORK(KFREE),LFREE,IERR)
      ! Change order in S_Z to decreasing instead of increasing value
      ! to follow usual convention
      DO I = 1, MULTIPESR/2
         J = MULTIPESR + 1 -I
         TEMP = ESRVAL(I,1,1,3)
         ESRVAL(I,1,1,3) = ESRVAL(J,1,1,3)
         ESRVAL(J,1,1,3) = TEMP
         DO K = 1,NZ_in_CI
            CALL DSWAP(MULTIPESR,ESRVAL(1,J,K,2),1,ESRVAL(1,I,K,2),1)
         END DO
      END DO

      IF (IPRESR .GT. 4) THEN
        DO I = 2,3
          write (lupri,'(A,I2/A)')
     &     ' test output after qdiag for S_z; ESRVAL matrix',I,
     &     ' (2 ~ eigenvectors, 3 ~ eigenvalues)'
          CALL PRQMAT(ESRVAL(1,1,1,I),MULTIPESR,MULTIPESR,
     &                MULTIPESR,MULTIPESR,NZ_in_CI,IQINDX,LUPRI)
        END DO
      END IF

C
C     Transform the HCI energy matrix to the new CI basis
C
      CALL DCOPY(N2MULT*NZ_in_CI,HCI,1,ESRVAL,1) ! using ESRVAL as work space
      CALL QTRANS('AOMO','S',0.0D0,
     &        MULTIPESR,MULTIPESR,MULTIPESR,MULTIPESR,
     &        ESRVAL,         MULTIPESR,MULTIPESR,NZ_in_CI,IQINDX,
     &        HCI,            MULTIPESR,MULTIPESR,NZ_in_CI,IQINDX,
     &        ESRVAL(1,1,1,2),MULTIPESR,MULTIPESR,NZ_in_CI,IQINDX,
     &        ESRVAL(1,1,1,2),MULTIPESR,MULTIPESR,NZ_in_CI,IQINDX,
     &        WORK(KFREE),LFREE,-1)

!     IF (IPRESR .GT. 3) THEN
        CALL HEADER(' H_CI matrix after transformation to S_z basis ',1)
        CALL PRQMAT(HCI,MULTIPESR,MULTIPESR,
     &              MULTIPESR,MULTIPESR,NZ_in_CI,IQINDX,LUPRI)
!     END IF
C
C     Transform the first MULTIPESR CI eigenvectors to S_z basis
C     (S_z is a real operator, so no problem in keeping NZ_in_CI)
C
      L_CIVEC2_R = NZCONF * MULTIPESR
      L_CIVEC2   = L_CIVEC2_R * NZ_in_CI
      CALL MEMGET2('REAL','CIVEC2',KCIVEC2 ,L_CIVEC2,WORK,KFREE,LFREE)
      CALL DCOPY(L_CIVEC2_R,CIEVEC(1,STATE_ESR,1),1,WORK(KCIVEC2),1)
      IF (NZ_in_CI .eq. 2) THEN
         KCIVEC2_I = KCIVEC2 + L_CIVEC2_R
         CALL DCOPY(L_CIVEC2_R,CIEVEC(1,STATE_ESR,2),1,
     &      WORK(KCIVEC2_I),1)
      END IF
      CALL QGEMM(NZCONF,MULTIPESR,MULTIPESR,1.0D0,
     &           'N','N',IQINDX,
     &           WORK(KCIVEC2),NZCONF,MULTIPESR,NZ_in_CI,
     &           'N','N',IQINDX,
     &           ESRVAL(1,1,1,2),MULTIPESR,MULTIPESR,NZ_in_CI,
     &           0.0D0,IQINDX,
     &           CIEVEC(1,STATE_ESR,1),NZCONF,N_CIESR,NZ_in_CI)
      CALL MEMREL('ESRPHASE-1',WORK,1,KCIVEC2,KFREE,LFREE)
      IF (IPRESR .GT. 10) THEN
        write (lupri,*) ' CIEVEC after transform. to S_z basis'
        CALL PRQMAT(CIEVEC(1,STATE_ESR,1),NZCONF,MULTIPESR,
     &              NZCONF,N_CIESR,NZ_in_CI,IQINDX,LUPRI)
      END IF
C
C     Step 2: Calculate S_x, S_y matrices in CIEVEC basis, using ESRVAL(...,1:3)
C     ======
C     (for print levels > 3, also calculate S_z matrix again, only relevant for testing
C     if code above is correct.)
C
C --- current code only ok for I_end = 3, see S_eff test below (S_z from
C above has been destroyed !
!     IF (IPRESR .GT. 3) THEN
         I_end = 3 ! S_x, S_y, S_z - S_z only to verify transformation in Step 1 was correct
!     ELSE
!        I_end = 2 ! only S_x and S_y
         ! S_z matrix is in ESRVAL(:,:,2,3) from Step 1
!     END IF
#if ESR_DEBUG > 3
      WRITE(LUPRI,*) 'Initial phase factors, should be (1,0):'
      do i = 1, multipesr
         write(lupri,*) phase_factors(i,1:2)
      end do
#endif
      DO I = 1,I_end
          I1 = IPSIGMA(I)
          INDXPR = LESRP(I1)
C
C         =============================
C
          IF (IPRESR.GT.1) THEN
            CALL PRSYMB(LUPRI,'=',75,0)
            WRITE(LUPRI,'(A,I3,A,3X,A)')
     &        '* ESRPHASE: calculating property no.',I,' :',
     &        PRPNAM(INDXPR)
            CALL PRSYMB(LUPRI,'.',75,0)
            CALL WRIXPR(I,INDXPR)
          END IF
C
          CALL ESRTRANS(INDXPR,ESRVAL(1,1,1,I),PHASE_factors,
     &                  CMO,PMAT,CIEVEC,WORK,KFREE,LFREE)
      ENDDO
C
C
C     SEFF        = S_eff (one of 0.0, 0.5, 1.0, ...)
C     MULTIPESR   = 2 * S_eff + 1
C     2 for doublet, 3 for triplet, etc.
C
      SEFF = (MULTIPESR - 1)/2.0D0
C
C     Eigenvalues and eigenvectors were sorted after increasing M_Seff value
C     above using QDIAG. Now we test if values are OK
C     (i.e. sufficiently close to -Seff, -Seff+1, ..., +Seff)
C
      WRITE (LUPRI,'(//A/)')
     &  '*ESRPHASE: checking S_Z(eff) and relative phases ...'
      IERR = 0
      DO J = 1,MULTIPESR
         SZJ = ESRVAL(J,J,1,3)
         K   = NINT(-SZJ+SEFF+1.0D0)
         WRITE (LUPRI,'(A,I5,F15.10,I5)')
     &     ' State no., <Sigma_Z>, index S_Z(eff) :',J,SZJ,K
         IF (K .NE. J) THEN
            IERR = IERR + 1
            WRITE(LUPRI,*) 'ERROR, State no. .ne. index S_Z(eff) !!'
         END IF
      END DO
      IF (IERR .GT. 0) THEN
         CALL QUIT('error: S_Z(eff) outside allowable bounds')
      END IF
C
C     ... Check relative phases using Sigma_X
C
C       Use Sigma_X and Sigma_Y to determine phase of M_S .gt. -J_eff
C       using the requirement that Re(<J+1|S_x|J>) .gt. 0 and Im(<J+1|S_y|J>) .gt. 0
C       (as for non-rel <2|S_x|1> = 0.5 and  <2|S_y|1> = 0.5 i)
C
      PHASE_J_R = 1.0D0 ! keep eigenvector for M_S = +S (J = 1) as it is
      PHASE_J_I = 0.0D0
      PHASE_factors(1,1) = PHASE_J_R
      PHASE_factors(1,2) = PHASE_J_I

      DO J = 1,MULTIPESR-1  ! M_S = S,S-1, ...,-S+1
      ! Find phase factor for eigenvector J+1
      ! Start with left transformation using established phase factors for j
      ! i.e. multiply with complex conjugate: PHASE_J_R - i PHASE_J_I
        JP1 = J+1 ! M_S - 1
        SX_R = ESRVAL(J,JP1,1,1)*PHASE_J_R + ESRVAL(J,JP1,2,1)*PHASE_J_I
        SX_I =-ESRVAL(J,JP1,1,1)*PHASE_J_I + ESRVAL(J,JP1,2,1)*PHASE_J_R
        SY_R = ESRVAL(J,JP1,1,2)*PHASE_J_R + ESRVAL(J,JP1,2,2)*PHASE_J_I
        SY_I =-ESRVAL(J,JP1,1,2)*PHASE_J_I + ESRVAL(J,JP1,2,2)*PHASE_J_R

!       IF (IPRESR .GT. 0) THEN
           WRITE (LUPRI,'(/A,F4.1,A)')
     &     "  Initial values before phase transformation"//
     &     "  for M_Seff =",SEFF-DFLOAT(J-1)
           WRITE (LUPRI,'(/A,F15.10,A,F15.10)')
     &     " <M_Seff|S_x|M_Seff - 1>  :",SX_R,' + i',SX_I
           WRITE (LUPRI,'( A,F15.10,A,F15.10)')
     &     " <M_Seff|S_y|M_Seff - 1>  :",SY_R,' + i',SY_I
!        END IF

      ! Now find phase factors for j+1 from the right
      ! so we get the desired structures for S_x and S_y

        ! Using (a+ib)(c+is) = (ac-bs) + i(as+bc)
        ! make <j|S_x|jp1> value purely real
        ! and  <j|S_y|jp1> value purely imaginary
        IF (ABS(SX_I) .gt. THR_ZERO) THEN
           PHI = -ATAN(SX_I,SX_R)
           C_PHI = COS(PHI)
           S_PHI = SIN(PHI)
           SX_R_tmp = C_PHI*SX_R - S_PHI*SX_I
           SX_I_tmp = C_PHI*SX_I + S_PHI*SX_R
           SY_R_tmp = C_PHI*SY_R - S_PHI*SY_I
           SY_I_tmp = C_PHI*SY_I + S_PHI*SY_R
           IF (ABS(SX_I_tmp) .gt. THR_ZERO) THEN
              write(lupri,*) 'ESRPHASE: error in Sx phase',
     &        SX_R_tmp,SX_I_tmp
              CALL QUIT('ESRPHASE: error in Sx phase')
           END IF
           SX_R = SX_R_tmp
           SX_I = 0.0D0
           SY_R = SY_R_tmp
           SY_I = SY_I_tmp

!          PHASE_R =  C_PHI*PHASE_J_R - S_PHI*PHASE_J_I
!          PHASE_I =  C_PHI*PHASE_J_I + S_PHI*PHASE_J_R
           PHASE_R =  C_PHI
           PHASE_I =  S_PHI
        ELSE
           SX_I = 0.0D0 ! remove numerical noise

!          PHASE_R =  PHASE_J_R
!          PHASE_I =  PHASE_J_I
           PHASE_R =  1.0D0
           PHASE_I =  0.0D0
        END IF

!       IF (IPRESR .GT. 1) THEN
           WRITE (LUPRI,'(/A)')
     &     "  After phase rotations; without sign fix:"
           WRITE (LUPRI,'(/A,F15.10,A,F15.10)')
     &     " <M_Seff|S_x|M_Seff - 1>  :",SX_R,' + i',SX_I
           WRITE (LUPRI,'( A,F15.10,A,F15.10)')
     &     " <M_Seff|S_y|M_Seff - 1>  :",SY_R,' + i',SY_I
!       END IF

        ! Next fix sign so S_x and S_y values are according to convention
        IF      (SX_R .gt. 0.0D0 .and. SY_I .lt. 0.0D0) then
           ! OK, do nothing, (a - ib)
            PHASE_JP1_R = +PHASE_R
            PHASE_JP1_I = +PHASE_I
            SWITCH_XY   = .FALSE.
        ELSE IF (SX_R .gt. 0.0D0 .and. SY_I .gt. 0.0D0) then
           ! (a - ib) = (b + ia)*(-i)
            PHASE_JP1_R = +PHASE_I
            PHASE_JP1_I = -PHASE_R
            SWITCH_XY   = .TRUE.
        ELSE IF (SX_R .lt. 0.0D0 .and. SY_I .lt. 0.0D0) then
           ! (a - ib) = (-b - ia)*(i)
            PHASE_JP1_R = -PHASE_I
            PHASE_JP1_I = +PHASE_R
            SWITCH_XY   = .TRUE.
        ELSE IF (SX_R .lt. 0.0D0 .and. SY_I .gt. 0.0D0) then
           ! (a -ib) = (-a + ib)*(-1)
            PHASE_JP1_R = -PHASE_R
            PHASE_JP1_I = -PHASE_I
            SWITCH_XY   = .FALSE.
        ELSE
           CALL QUIT('One of Sx and Sy is zero.'//
     &        ' Probably pseudo-singlet state included.')
        END IF

        PHASE_factors(JP1,1) = PHASE_JP1_R
        PHASE_factors(JP1,2) = PHASE_JP1_I

!       IF (ABS(SY_R) .gt. THR_ZERO) THEN
        IF (ABS(SY_R) .gt. 1.0D-4) THEN
           write(lupri,*)'ESRPHASE error: S_y not purely imaginary'//
     &        ' although S_x is purely real'
           write(lupri,*) 'S_x',SX_R,SX_I
           write(lupri,*) 'S_y',SY_R,SY_I
           CALL QUIT('ESRPHASE error: S_y not purely imaginary'//
     &        ' although S_x is purely real')
        ELSE IF (ABS(SY_R) .gt. 1.0D-6) THEN
           write(lupri,*)'ESRPHASE info: S_y not purely imaginary'//
     &        ' although S_x is purely real'
           write(lupri,*) 'S_x',SX_R,SX_I
           write(lupri,*) 'S_y',SY_R,SY_I
           write(lupri,*)'ESRPHASE info: Calculation continues,'//
     &        ' assuming deviation is numerical noise.'
        END IF

C       re-evaluate S_x and S_y values for JP1 for print
C       1) multiply from left with phase factors for J, complex conjugated:
        SXR = ESRVAL(J,JP1,1,1)*PHASE_J_R + ESRVAL(J,JP1,2,1)*PHASE_J_I
        SXI =-ESRVAL(J,JP1,1,1)*PHASE_J_I + ESRVAL(J,JP1,2,1)*PHASE_J_R
        SYR = ESRVAL(J,JP1,1,2)*PHASE_J_R + ESRVAL(J,JP1,2,2)*PHASE_J_I
        SYI =-ESRVAL(J,JP1,1,2)*PHASE_J_I + ESRVAL(J,JP1,2,2)*PHASE_J_R
C       1) multiply from left with phase factors for J, complex conjugated:
        SX_R = SXR*PHASE_JP1_R - SXI*PHASE_JP1_I
        SX_I = SXR*PHASE_JP1_I + SXI*PHASE_JP1_R
        SY_R = SYR*PHASE_JP1_R - SYI*PHASE_JP1_I
        SY_I = SYR*PHASE_JP1_I + SYI*PHASE_JP1_R

        WRITE (LUPRI,'(/A,F4.1)')
     &     "  Final values after phase transformation"//
     &     " for M_Seff =",SEFF-DFLOAT(J-1)

        IF (SWITCH_XY) THEN
           WRITE(LUPRI,'(/A)') ' Note: S_x and S_y operators have'//
     &     ' been swapped to get correct phases.'
           SY_R_tmp  = SY_R
           SY_R = SX_R
           SX_R = SY_R_tmp
           SY_I_tmp  = SY_I
           SY_I = SX_I
           SX_I = SY_I_tmp
        END IF
        WRITE (LUPRI,'(/A,F15.10,A,F15.10)')
     &     " <M_Seff|S_x|M_Seff - 1>  :",SX_R,' + i',SX_I
        WRITE (LUPRI,'( A,F15.10,A,F15.10)')
     &     " <M_Seff|S_y|M_Seff - 1>  :",SY_R,' + i',SY_I

C       for next J in do loop
        PHASE_J_R = PHASE_JP1_R
        PHASE_J_I = PHASE_JP1_I
C       ... here PHASE is new phase of vector MJ-1

      END DO ! J = 1, MULTIPESR-1

!     IF (IPRESR .GT. 1) THEN
         WRITE(LUPRI,'(/A)') '  M_Seff : Final phase factor on |M_Seff>'
         do j = 1, multipesr
            write(lupri,'(F6.1,A,2F15.10)')
     &         SEFF-DFLOAT(J-1),'   :',phase_factors(j,1:2)
         end do

         CALL ESR_PHASE_TRANS(PHASE_factors, ESRVAL(1,1,1,1), MULTIPESR) ! S_x
         CALL ESR_PHASE_TRANS(PHASE_factors, ESRVAL(1,1,1,2), MULTIPESR) ! S_y
         IF (SWITCH_XY) THEN
            WRITE(LUPRI,'(/A)')
     &      ' Note: S_x and S_y operators are swapped in order'//
     &      ' to get correct phases of final results.'
            L_ESRVAL = MULTIPESR*MULTIPESR*2
            CALL DSWAP(L_ESRVAL,ESRVAL(1,1,1,1),1,ESRVAL(1,1,1,2),1)
         END IF

         CALL HEADER('Phase updated property matrix for S_x',1)
         CALL PRQMAT(ESRVAL(1,1,1,1),MULTIPESR,MULTIPESR,
     &            MULTIPESR,MULTIPESR,2,IQINDX,LUPRI)

         CALL HEADER('Phase updated property matrix for S_y',1)
         CALL PRQMAT(ESRVAL(1,1,1,2),MULTIPESR,MULTIPESR,
     &            MULTIPESR,MULTIPESR,2,IQINDX,LUPRI)
!     END IF
C
C **  Form S**2 matrix in ESRVAL(...,4) in CIEVEC basis
C
      FAC  = 0.0D0
      DO I = 1,3
        CALL QGEMM(MULTIPESR,MULTIPESR,MULTIPESR,1.0D0,
     &           'N','N',IQINDX, ESRVAL(1,1,1,I),MULTIPESR,MULTIPESR,2,
     &           'N','N',IQINDX, ESRVAL(1,1,1,I),MULTIPESR,MULTIPESR,2,
     &           FAC,    IQINDX, ESRVAL(1,1,1,4),MULTIPESR,MULTIPESR,2)
        FAC = 1.0D0
      END DO
C
      CALL HEADER(' S**2 matrix in CI basis ',1)
      CALL PRQMAT(ESRVAL(1,1,1,4),MULTIPESR,MULTIPESR,
     &            MULTIPESR,MULTIPESR,2,IQINDX,LUPRI)

C
      CALL QEXIT('ESRPHASE')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesr1 */
      SUBROUTINE PRPESR1(ESRVAL,PHASE_factors,LESRST,CMO,PMAT,
     &                   CIEVEC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIESR/
C
C     A. Noerager - Mar 2000
C     Based on prpex1 by Trond Saue May 27 1996
C     Revised August 2002 Hans Joergen Aa. Jensen
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      LOGICAL   LESRST(2,*), REAL_PMLF
      REAL*8    ESRVAL(MULTIPESR,MULTIPESR,2,*),
     &          PHASE_factors(MULTIPESR,2),
     &          CMO(*),PMAT(*),
     &          CIEVEC(NZCONF,N_CIESR,*), WORK(*)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   IPCMO(2)
C
C Used from common blocks:
C  dcbgen : LU1INT,?
C  dcbopt : KZCONF,NZCONF
C
#include "dgroup.h"

#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbesr.h"
#include "dcbxpr.h"

#include "dcbprp.h"
#include "dcbdhf.h"
C
      KFRSAV = KFREE
      call memchk('prpesr1 1',work,1)
C
C     Initialize
C     ==========
C
      IF (N_CIESR .lt. MULTIPESR) CALL QUIT('N_CIESR .lt. MULTIPESR')
      CALL DZERO(ESRVAL,(MULTIPESR**2)*2*NESRP)
      N2ESR = 2*NESRP
      CALL LSET(N2ESR,.FALSE.,LESRST)

C
C     Generate density matrix
C     =======================
C
C      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6)

C     Expectation values normal procedure
C    ************************************
C
C     TEMPLATE for later, when calculation of transition density matrices have been implemented
C
c       CALL MEMGET('REAL',KDMAT,NZ*N2BBASX*(NOPEN+1),WORK,KFREE,LFREE)
c       CALL GENDEN(work(kDMAT),CMO,1,IPRESR)
c       DO IOPEN = 1,NOPEN
c          CALL ACTDEN(work(kDMAT+N2BBASXQ*IOPEN),CMO,IOPEN,IPRESR)
c       END DO
c       CALL TITLER('Expectation values','*',127)
c       WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
c       CALL PRSYMB(LUPRI,'-',76,0)
C
C       Allocate memory for property matrix in molfdir format:
C        NZ_in_CI = MIN(2,NZ)
C        CALL MEMGET('REAL',KPMLF  ,NZ_in_CI*(2*NASHT)**2,WORK,KFREE,LFREE)
C        CALL MEMGET('REAL',KMUUUU ,0              ,WORK,KFREE,LFREE)
C        LPCI = N2CIESR
C        IF (NZ_in_CI .EQ. 2) LPCI = 2*LPCI
C        CALL MEMGET('REAL',KPR    ,LPCI           ,WORK,KFREE,LFREE)
C        KPI = KPR + N2CIESR
C
C
        DO I = 1,NESRP
          INDXPR = LESRP(I)
C
C         Extract symmetry information
C         =============================
C
          ISYM = IPRPSYM(INDXPR)
          IREP = ISYM - 1
          ITIM = IPRPTIM(INDXPR)
          IOPSY = JBTOF(ISYM-1,1)
chj       IF (IPRESR.GT.1 .OR. ITIM.GE.0) THEN
          IF (IPRESR.GT.1) THEN
             CALL PRSYMB(LUPRI,'=',75,0)
             WRITE(LUPRI,'(A,I3,A,A/A,3I5)')
     &       '* PRPESR1: calculating property no.',I,' : ',
     &         PRPNAM(INDXPR),
     &       '  Boson symmetry, g/u symmetry, time reversal :',
     &         ISYM,IOPSY,ITIM
          END IF
C
          CALL ESRTRANS(INDXPR,ESRVAL(1,1,1,I),PHASE_factors,
     &                  CMO,PMAT,CIEVEC,WORK,KFREE,LFREE)
C
c
c      print *, ' ..... test output of property matrix in MO basis ...'
c      print *, 'isym,iopsy:',isym,iopsy
c      print *, 'icmoq  :',icmoq
c      print *, 'nfbas  :',nfbas(1,0),nfbas(2,0)
c      print *, 'npsh   :',npsh
c      print *, 'nish   :',nish
c      print *, 'IPCMO,nash:',IPCMO,nash
C
c     call memchk('prpesr1 2',work,1)
C          CALL PRPMAT(INDXPR,IOPSY,PMAT,.TRUE.,WORK,CMO,IBEIG,IPCMO,NASH,
C     &                WORK,KFREE,LFREE,IPRESR)
c    &                WORK,KFREE,LFREE,10)
C      call memchk('prpesr1 3',work,1)
C          CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,40)
chjdbg    CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,IPRESR)
C
C         Calculate CI property matrix over determinants
C
C          IF (NZ.GT.1) THEN
C             IHRM = 1
C             REAL_PMLF = .TRUE.
c          ELSE
C             JQ = IPQTOQ(1,IREP)
C             REAL_PMLF = JQ.EQ.1 .OR. JQ.EQ.3
C             IF (ITIM.EQ.-1) REAL_PMLF = .NOT. REAL_PMLF
C             IF (REAL_PMLF) THEN
C                IHRM = 1
C             ELSE
C                IHRM = -1
C             END IF
C          END IF
C          CALL GASCIP_MAKEH(.FALSE.,IHRM,NZCONF,WORK(KZCONF),
C     &         WORK(KPR),WORK(KPI),WORK(KPMLF),WORK(KMUUUU))
C          IF (ITIM.EQ.1) THEN
C             WRITE (LUPRI,*)' WARNING: Inactive contribution not added!'
C TODO
c            DO I = 1,NZCONF
c               JPR = KPR + (I-1)*N_CIESR + (I-1)
c               WORK(JPR) = WORK(JPR) + PCORE
c            END DO
C          END IF
c         IF (IPRESR .GT. 5) THEN
C             IF (NZ_in_CI .EQ. 2) THEN
C                WRITE (LUPRI,'(/A)')
C     &          ' Real part of det-CI property matrix'
C             ELSE IF (REAL_PMLF) THEN
C                WRITE (LUPRI,'(/A)')
C     &          ' Real det-CI property matrix'
C             ELSE
C                WRITE (LUPRI,'(/A)')
C     &          ' Imaginary det-CI property matrix'
C             END IF
C             CALL OUTPUT(WORK(KPR),1,N_CIESR,1,N_CIESR,
C     &           N_CIESR,N_CIESR,1,LUPRI)
C             IF (NZ_in_CI .EQ. 2) THEN
C                WRITE (LUPRI,'(/A)')
C     &          ' Imag part of det-CI property matrix'
C                CALL OUTPUT(WORK(KPI),1,N_CIESR,1,N_CIESR,
C     &              N_CIESR,N_CIESR,1,LUPRI)
C             END IF
c         END IF
C
C         Transform to eigenvector basis
C
C          IF (NZ_in_CI .EQ. 2) THEN
C             JZ = 1
C          ELSE IF (REAL_PMLF) THEN
C             JZ = 1
Cunoedvendig CALL DZERO(ESRVAL(1,1,2,I),N2CIESR)
C          ELSE
C             JZ = 2
Cunoedvendig CALL DZERO(ESRVAL(1,1,1,I),N2CIESR)
C          END IF
C
C          CALL QTRANS('AOMO','S',0.0D0,N_CIESR,N_CIESR,N_CIESR,N_CIESR,
C     &        WORK(KPR),N_CIESR,N_CIESR,NZ_in_CI,IQINDX(JZ),
C     &        ESRVAL(1,1,1,I),N_CIESR,N_CIESR,2,IQINDX,
C     &        CIEVEC,N_CIESR,N_CIESR,NZ_in_CI,IQINDX,
C     &        CIEVEC,N_CIESR,N_CIESR,NZ_in_CI,IQINDX,
C     &        WORK(KFREE),LFREE,-1)
C
ccc       IF ( MOD(NAELEC,2) .EQ. 1) THEN
C         ... odd number of electrons
ccc          CALL PRPESRODD(ESRVAL(1,1,1,I),IREP,ITIM,WORK(KPMLF))
ccc       ELSE
C         ... even number of electrons
ccc          CALL PRPESREVEN(ESRVAL(1,1,1,I),IREP,ITIM,WORK(KPMLF))
ccc       END IF
c      print *, ' ..... end of test output of MO prp. matrix ...'
Chj-end test output

c         CALL MEMGET('REAL',KPTRI  ,NZ*NNBBASX,WORK,KFREE,LFREE)
c         CALL MEMGET('REAL',KOP1INT,   NNBBASX,WORK,KFREE,LFREE)
c         CALL PRPESR2(LU1INT,INDXPR,ESRVAL(I),LESRST(1,I),
c    &        work(kDMAT),PMAT,work(kPTRI),work(kOP1INT),
c    &        WORK(KFREE),LFREE,IPRESR)
c         FMT = MXFORM(ESRVAL(I),15)
c         WRITE(LUPRI,'(4X,A,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
c    &       PRPNAM(INDXPR),' : ',ESRVAL(I),'a.u.',
c    &       LESRST(1,I),LESRST(2,I)
c         CALL FLSHFO(LUPRI)
        ENDDO
c       CALL PRSYMB(LUPRI,'-',76,0)
c       WRITE(LUPRI,'(4X,A)')
c    &   's0 = T : Expectation value zero by point group symmetry.'
c       WRITE(LUPRI,'(4X,A)')
c    &   't0 = T : Expectation value zero by time reversal symmetry.'
C
C
      call memchk('prpesr1 4',work,1)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ESRTRANS(INDXPR,ESRVAL,PHASE_factors,
     &                    CMO,PMAT,CIEVEC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate one-electron property INDXPR in CI basis (in CIEVEC)
C     and return result in ESRVAL.
C
C     Input:
C       INDXPR : index of property operator
C       PHASE_factors : phase factor for each CI vector
C       CMO    : MO coefficients
C       CIEVEC : CI vectors
C
C     Output:
C       ESRVAL : Property matrix elements in CI basis
C       PMAT   : Property matrix elements in MO basis
C
C     Scratch:
C       WORK(KFREE:KFREE+LFREE-1)
C
C     Based on PRPESR1
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C     for PRQMAT print of complex matrices:
      integer, parameter ::NZCMPLX = 2, IQINDX(2) = (/ 1,2 /)
C
      LOGICAL   REAL_PMLF
      REAL*8    ESRVAL(MULTIPESR,MULTIPESR,2),
     &          PHASE_factors(MULTIPESR,2),CMO(*),PMAT(*),
     &          CIEVEC(NZCONF,N_CIESR,*), WORK(*)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   IPCMO(2)
C
C Used from common blocks:
C  dcbgen : LU1INT,?
C  dcbopt : KZCONF,NZCONF
C  dgroup : NZ, NZ_in_CI (NZ_in_CI = 1 if CI vectors real, =2 if complex)
C
#include "dcbgen.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbesr.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
C
      CALL QENTER('ESRTRANS')
      KFRSAV = KFREE
C
C     Initialize
C     ==========
C
      IF (N_CIESR .lt. MULTIPESR) CALL QUIT('N_CIESR .Lt. MULTIPESR')
C
      CALL DZERO(ESRVAL,(MULTIPESR**2)*2)
C     Allocate memory for property matrix in molfdir format:
      CALL MEMGET2('REAL','PMLF ',
     &   KPMLF  ,NZ_in_CI*(2*NASHT)**2,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','MUUUU',
     &   KMUUUU ,0,WORK,KFREE,LFREE)

!     IF (INDEX(PRPNAM(INDXPR),'magnetic field') .gt. 0 .or.
!    &    INDEX(PRPNAM(INDXPR),'g_1 MAGFLD') .gt. 0 ) THEN
!        ACCURACY = THR_CIESR ! does not make sense to request higher accuracy than CI convergence
!        ! NB! this defines g shifts in 10 ppm if THR_CIESR = 1.D-6, ppt if THR_CIESR = 1.D-4 (default)
!     ELSE
!        ACCURACY = MAX(THR_CIESR,1.D-4) ! 3 digits enough for HFI etc. - other errors are bigger
!     END IF
! 04-Oct-2022 hjaaj: an accuracy of 1.d-4 cause an HFC for I in NI to
! change from 58.xxx to 79.xxx !!! This was unexpected for me. Therefore
! now always THR_CIESR
      ACCURACY = THR_CIESR * 0.1D0
! if converged to e.g. 1.d-4 we will
! have coefficients less than 1.d-4 which are important !!!
C
C     Extract symmetry information
C     =============================
C
      ISYM = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      ITIM = IPRPTIM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
chj   IF (IPRESR.GT.1 .OR. ITIM.GE.0) THEN
      IF (IPRESR.GT.2) THEN
         CALL PRSYMB(LUPRI,'-',75,0)
         WRITE(LUPRI,'(2A/A,3I5/A,1P,D10.2)')
     &       '* ESRTRANS: calculating CI property matrix    : ',
     &         PRPNAM(INDXPR),
     &       '  Boson symmetry, g/u symmetry, time reversal :',
     &         ISYM,IOPSY,ITIM,
     &       '  Requested relative accuracy in CI property matrix :',
     &         ACCURACY
      END IF
C
chj       IF (ITIM .GE. 0) CALL QUIT(
chj  &    'ESRTRANS only implemented for time-antisymmetric operators')
C

      IPCMO(1) = icmoq(1) + nfbas(1,0)*(npsh(1)+nish(1))
      IPCMO(2) = icmoq(2) + nfbas(2,0)*(npsh(2)+nish(2))
Chj-start test output
c         iprtst = 10
c
!d     WRITE(LUPRI,*)
!d   & ' ..... test output of property matrix in MO basis ...'
c      print *, 'isym,iopsy:',isym,iopsy
c      print *, 'isym,iopsy:',isym,iopsy
c      print *, 'icmoq  :',icmoq
c      print *, 'nfbas  :',nfbas(1,0),nfbas(2,0)
c      print *, 'npsh   :',npsh
c      print *, 'nish   :',nish
c      print *, 'IPCMO,nash:',IPCMO,nash
C
      CALL MEMGET2('REAL','IBEIG',KIBEIG ,NORBT,WORK,KFREE,LFREE)
      CALL REACMO(LUCOEF,'DFCOEF',DUMMY,DUMMY,WORK(KIBEIG),
     &                    TOTERG,8)
      CALL PRPMAT(INDXPR,IOPSY,PMAT,.TRUE.,WORK,CMO,WORK(KIBEIG),
     &            IPCMO,NASH,WORK,KFREE,LFREE,IPRESR)
      CALL MEMREL('after PRPMAT',WORK,1,KIBEIG,KFREE,LFREE)
c    &            WORK,KFREE,LFREE,10)
      CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,IPRESR)
C
C     Calculate CI property matrix over determinants
C
      IF (NZ.GT.1) THEN
         REAL_PMLF = .TRUE.
         IHRM = 1
      ELSE
         JQ = IPQTOQ(1,IREP)
         REAL_PMLF = JQ.EQ.1 .OR. JQ.EQ.3
         IF (ITIM.EQ.-1) REAL_PMLF = .NOT. REAL_PMLF
         IF (REAL_PMLF) THEN
            IHRM = 1
         ELSE
            IHRM = -1
         END IF
      END IF

C     In which component of ESRVAL (real or imaginary) should the
C     property matrix be ?

      IF (NZ_in_CI .EQ. 2) THEN
         JZ = 1
      ELSE IF (REAL_PMLF) THEN
         JZ = 1
      ELSE
         JZ = 2
      END IF

      IF (MAX_CIESR_IT .GT. 0) THEN
C
C         Make property matrix in ESRVAL directly from the MULTIPESR first CI vectors
C           CIEVEC(NZCONF, N_CIESR, NZ_in_CI)
C           ESRVAL(MULTIPESR, MULTIPESR, 2)
C

#if ESR_DEBUG > 20
         write(lupri,*)
     &   'Check hermiticity before GASCIP_MAKE_P1MAT ', IHRM
         write(lupri,*)
     &    'N_CIESR,NZCONF,KZCONF,MULTIPESR,IHRM,KPMLF',
     &    N_CIESR,NZCONF,KZCONF,MULTIPESR,IHRM,KPMLF
         WRITE(LUPRI,*)
     &   'GASCIP_MAKE_P1MAT: Real part of the N_CIESR CI vectors'//
     &   ' before MAKE_P1MAT'
         CALL OUTPUT(CIEVEC,1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
         IF (NZ_in_CI .EQ. 2) THEN
          WRITE(LUPRI,*)
     &      'GASCIP_MAKE_P1MAT: Imag part of the N_CIESR CI vectors'//
     &      '  BEFORE MAKE_P1MAT'
          CALL OUTPUT(CIEVEC(1,1,2),1,NZCONF,1,N_CIESR,
     &                NZCONF,N_CIESR,-1,LUPRI)
         END IF
         write(lupri,*) 'N_CIESR, NZCONF', N_CIESR, NZCONF
#endif
         CALL GASCIP_MAKE_P1MAT(CIEVEC(1,STATE_ESR,1),N_CIESR,NZCONF,
     &                         WORK(KZCONF),
     &                         ESRVAL(1,1,JZ),MULTIPESR,IHRM,
     &                         WORK(KPMLF),ACCURACY,
     &                         WORK(KFREE),LFREE,IPRESR)
C       CALL GASCIP_MAKE_P1MAT(CIVEC,N_CIVEC,NDET,IDET,
C    &                         P1MAT_CI,N_P1MAT_CI,IHRM_P1,
C    &                         P1MAT_MO,  ACCURACY, WORK,LWORK, IPRINT)
      ELSE  ! here if MAX_CIESR_IT .le. 0
         ACCURACY = 0.0D0 ! no screening
C
C         Our first non-iterative code with full, explicit CI matrices.
C         Construct complete property matrix over determinants
C         and transform to the desired property matrix in ESRVAL from
C         the MULTIPESR first CI vectors with QTRANS afterwards.
C
         LPCI = NZ_in_CI*NZCONF*NZCONF
         CALL MEMGET2('REAL','PR',KPR,LPCI,WORK,KFREE,LFREE)
         KPI  = KPR + NZCONF*NZCONF
         CALL GASCIP_MAKEH(.FALSE.,IHRM,NZCONF,WORK(KZCONF),
     &         WORK(KPR),WORK(KPI),WORK(KPMLF),WORK(KMUUUU))
         IF (IPRESR .GT. 10) THEN
            IF (NZ_in_CI .EQ. 2) THEN
               WRITE (LUPRI,'(/A)')
     &          ' Real part of det-CI property matrix'
            ELSE IF (REAL_PMLF) THEN
               WRITE (LUPRI,'(/A)')
     &          ' Real det-CI property matrix'
            ELSE
               WRITE (LUPRI,'(/A)')
     &          ' Imaginary det-CI property matrix'
            END IF
            CALL OUTPUT(WORK(KPR),1,NZCONF,1,NZCONF,
     &           NZCONF,NZCONF,1,LUPRI)
            IF (NZ_in_CI .GT. 1) THEN
               WRITE (LUPRI,'(/A)')
     &          ' Imag part of det-CI property matrix'
               CALL OUTPUT(WORK(KPI),1,NZCONF,1,NZCONF,
     &              NZCONF,NZCONF,-1,LUPRI)
            END IF
         END IF
C
C        Transform to eigenvector basis
C

         CALL QTRANS('AOMO','S',0.0D0,NZCONF,NZCONF,MULTIPESR,MULTIPESR,
     &        WORK(KPR),NZCONF,NZCONF,NZ_in_CI,IQINDX(JZ),
     &        ESRVAL, MULTIPESR,MULTIPESR, 2,IQINDX,
     &        CIEVEC(1,STATE_ESR,1), NZCONF,N_CIESR, NZ_in_CI,IQINDX,
     &        CIEVEC(1,STATE_ESR,1), NZCONF,N_CIESR, NZ_in_CI,IQINDX,
     &        WORK(KFREE),LFREE,-1)
       END IF
C
C      Multiply with phase factors
C
       IF (IPRESR .GE. 5) THEN
          WRITE(LUPRI,'(/A)') 'ESRTRANS: ESRVAL before phase factors'
          CALL PRQMAT(ESRVAL,MULTIPESR,MULTIPESR,
     &                MULTIPESR,MULTIPESR,NZCMPLX,IQINDX,LUPRI)
       END IF

       CALL ESR_PHASE_TRANS(PHASE_factors, ESRVAL, MULTIPESR)
C
       IF (IPRESR .GE. 0) THEN
          CALL HEADER('ESRTRANS results in CI eigenvector basis',1)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(A,3X,A/A,3I5/A,1P,D10.2)')
     &       '* ESRTRANS: CI property matrix for operator :',
     &         PRPNAM(INDXPR),
     &       '  Boson and fermion symmetry, time reversal :',
     &         ISYM,IOPSY,ITIM,
     &       '  Requested relative accuracy in CI property matrix :',
     &         ACCURACY
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(1,INDXPR)
          CALL PRQMAT(ESRVAL,MULTIPESR,MULTIPESR,
     &                MULTIPESR,MULTIPESR,NZCMPLX,IQINDX,LUPRI)
      END IF
      IF (ITIM.EQ.1) THEN
         WRITE (LUPRI,*)' NOTE: Inactive contribution is not included!'
C TODO
c        DO I = 1,MULTIPESR
c           ESRVAL(I,I,1) = ESRVAL(I,I,1) + PCORE
c        END DO
      END IF
C
      CALL MEMREL('ESRTRANS',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('ESRTRANS')
      RETURN
      END

      SUBROUTINE ESR_PHASE_TRANS(PHASE_factors, ESRVAL, MULTIPESR)

!     ESRVAL_new(i,j) = P*(i) ESRVAL(i,j) P(j)

      implicit none
      REAL*8  PHASE_factors(MULTIPESR,2), ESRVAL(MULTIPESR,MULTIPESR,2)
      INTEGER MULTIPESR

      INTEGER I,J
      REAL*8  VAL_IJ_R, VAL_IJ_I

#if ESR_DEBUG > 2
      WRITE(6,*) 'ESR_PHASE_TRANS Phase factors:'
      do i = 1, multipesr
         write(6,*) phase_factors(i,1:2)
      end do
#endif

      DO j = 1, MULTIPESR
         DO i = 1, MULTIPESR

            VAL_IJ_R =
     &         PHASE_factors(i,1)*PHASE_factors(j,1)*ESRVAL(I,J,1)
     &       + PHASE_factors(i,2)*PHASE_factors(j,2)*ESRVAL(I,J,1)
     &       - PHASE_factors(i,1)*PHASE_factors(j,2)*ESRVAL(I,J,2)
     &       + PHASE_factors(i,2)*PHASE_factors(j,1)*ESRVAL(I,J,2)
            VAL_IJ_I =
     &         PHASE_factors(i,1)*PHASE_factors(j,1)*ESRVAL(I,J,2)
     &       + PHASE_factors(i,1)*PHASE_factors(j,2)*ESRVAL(I,J,1)
     &       - PHASE_factors(i,2)*PHASE_factors(j,1)*ESRVAL(I,J,1)
     &       + PHASE_factors(i,2)*PHASE_factors(j,2)*ESRVAL(I,J,2)

            ESRVAL(I,J,1) = VAL_IJ_R
            ESRVAL(I,J,2) = VAL_IJ_I

         END DO
      END DO
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesrodd */
      SUBROUTINE PRPESRODD(ESRVAL,IREP,ITIM,PMLF)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
C dgroup.h : NZ, NZ_in_CI
#include "dcbesr.h"
#include "dgroup.h"
#include "dcborb.h"
C
      REAL*8    ESRVAL(N_CIESR,N_CIESR,2),PMLF(NASHT,2,NASHT,2,*)
      LOGICAL   REAL_PMLF
C
      if (2*NASHT .ne. N_CIESR) then
         print *, 'PRPESRODD error'
         print *,'N_CIESR,NASHT',N_CIESR,NASHT
         stop 'PRPESRODD error'
      end if
      N = N_CIESR*N_CIESR*NZ_in_CI
      IF (NZ_in_CI .EQ. 2) THEN
         call dcopy(N,PMLF,1,ESRVAL,1)
      ELSE
C        ... real group, NZ_in_CI = NZ = 1
C        Find out if matrix is pure real or pure imaginary
C        (See QFC2MFC)
         JQ = IPQTOQ(1,IREP)
         REAL_PMLF = JQ.EQ.1 .OR. JQ.EQ.3
         IF (ITIM.EQ.-1) REAL_PMLF = .NOT. REAL_PMLF
         IF (REAL_PMLF) THEN
            call dcopy(N,PMLF,1,ESRVAL(1,1,1),1)
            call dzero(ESRVAL(1,1,2),N)
         ELSE
            call dzero(ESRVAL(1,1,1),N)
            call dcopy(N,PMLF,1,ESRVAL(1,1,2),1)
         END IF
      END IF
#ifdef OLDCODE
          if (nash(1) .eq. 1 .and. nash(2) .eq. 0) then
             JMOSY = 1
          else if (nash(1) .eq. 0 .and. nash(2) .eq. 1) then
             JMOSY = 2
          else
             call quit('Only one open shell implemented yet!')
          end if
          DO J = 1,NZ
             PVAL = PMAT(1+I2ASHX(JMOSY,JMOSY)+N2ASHX*(J-1))
             JQ = IPQTOQ(J,IREP)
             IF (JQ .EQ. 2) THEN
                ESRVAL(1,1,1) =  PVAL
                ESRVAL(2,2,1) = -PVAL
             ELSE IF (JQ .EQ. 3) THEN
                ESRVAL(1,2,2) = -PVAL
                ESRVAL(2,1,2) =  PVAL
             ELSE IF (JQ .EQ. 4) THEN
                ESRVAL(1,2,1) =  PVAL
                ESRVAL(2,1,1) =  PVAL
             ELSE
C               ... JQ .eq. 1
                ESRVAL(1,1,2) =  PVAL
                ESRVAL(2,2,2) =  PVAL
                WRITE(LUPRI, *) 'test PVAL zero?: JQ,J,IREP,JMOSY,PVAL',
     &             JQ,J,IREP,JMOSY,PVAL
             END IF
          END DO
#endif
          RETURN
          END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesreven */
      SUBROUTINE PRPESREVEN(ESRVAL,IREP,ITIM,PMLF)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
C dgroup.h : NZ, NZ_in_CI
#include "dcbesr.h"
#include "dgroup.h"
#include "dcborb.h"
C
      REAL*8    ESRVAL(N_CIESR,N_CIESR,2),PMLF(NASHT,2,NASHT,2,*)
      LOGICAL   REAL_PMLF
C
      CALL QUIT(
     &'ESR only implemented for odd number of active electrons')
      N = N_CIESR*N_CIESR*NZ_in_CI
      IF (NZ_in_CI .EQ. 2) THEN
         call quit('PRPESREVEN: NZ_in_CI=2 not implemented.')
C
      ELSE
C        ... real group, NZ_in_CI = NZ = 1
C        Find out if matrix is pure real or pure imaginary
C        (See QFC2MFC)
         JQ = IPQTOQ(1,IREP)
         REAL_PMLF = JQ.EQ.1 .OR. JQ.EQ.3
         IF (ITIM.EQ.-1) REAL_PMLF = .NOT. REAL_PMLF
         IF (REAL_PMLF) THEN
C
         ELSE
C
         END IF
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#ifdef WORK_IN_PROGRESS
C  /* Deck PRPESR2 */
      SUBROUTINE PRPESR2(LU1INT,INDXPR,PVAL,LESRST,
     &           DMAT,PMAT,PTRI,OP1INT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIESR/
C
C     A. Noerager - Mar 2000
C     Based on prpex2 by Trond Saue May 27 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbesr.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
C
      LOGICAL   FIRST(4),LESRST(2), RODHF
      REAL*8    DMAT(N2BBASX,NZ,*),PMAT(N2BBASX,NZ),
     &          PTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
C     Check if single open shell calculation (.OPEN 1 1/2,0 or .OPEN 1 1/0,2)
C
      RODHF = NOPEN.EQ.1 .AND. NAELEC.EQ.1 .AND. NASHT.EQ.1
C
C     Check symmetry
C
      PVAL = 0.0D0
C
C     1. Molecular point group symmetry
C
      IREP = IPRPSYM(INDXPR)-1
      IF(IREP.NE.0) LESRST(1) = .TRUE.
C
C     2. Time reversal symmetry
C
      ITIM = IPRPTIM(INDXPR)
      IF(ITIM.NE.1) LESRST(2) = .TRUE.
C
C     If not totally symmetric, return with zero value
C
      IF(LESRST(1) .OR. (LESRST(2) .AND. .NOT.RODHF)) RETURN
C
C     Totally symmetric operator ==> proceed
C
      CALL LSET(NZ,.TRUE.,FIRST)
      CALL PRPMAO(LU1INT,INDXPR,.TRUE.,WORK,.TRUE.,NNBBASX,
     &             PTRI,OP1INT,FIRST,IPRINT)
C
C     Get expectation value, using only the matrices containing
C     integrals
C
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        IQ = IPQTOQ(IZ,IREP)
        IH = IHQMAT(IQ,ITIM)
        IF(IH.EQ.1) THEN
          CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
        ELSEIF(IH.EQ.2) THEN
          CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
        ENDIF
        CALL BUTOBS(PMAT(1,IZ),1,WORK,LWORK)
          IF (LESRST(2) .AND. RODHF) THEN
C           ... < p | O | p > for ITIM = -1 :
            PVAL = PVAL +
     &         DF(1)*DDOT(N2BBASX,DMAT(1,IZ,2),1,PMAT(1,IZ),1)
          ELSE
            DO IOPEN = 0,NOPEN
              PVAL = PVAL + DF(IOPEN)
     &           *DDOT(N2BBASX,DMAT(1,IZ,IOPEN+1),1,PMAT(1,IZ),1)
            END DO
          ENDIF
      ENDIF
      ENDDO
      PVAL = PVAL + PVAL
C
      RETURN
      END
#endif /* WORK_IN_PROGRESS */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ESR_HFC */
      SUBROUTINE PRP_ESR_HFC(ESRVAL,SWITCH_XY,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for ESR properties
C     This routine is modelled after PRP_SHIELD
C
C     Written by A.Noerager APR 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbesr.h"
#include "nuclei.h"
C
      REAL*8    ESRVAL(MULTIPESR,MULTIPESR,2,*),WORK(*)
      LOGICAL   SWITCH_XY
      logical, allocatable :: DOATOM(:)

      KFRSAV = KFREE
      NCMAT =  NUCDEP*9
      CALL MEMGET2('REAL','TMAT',KTMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TISO',KTISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBUF',KIBUF,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IATM',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','RATM',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','CSTRA',KCSTRA,
     &   NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SCTRA',KSCTRA,
     &   NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TMPMAT',KTMPMAT,NCMAT,WORK,KFREE,LFREE)
C
C     SETATM code -29: accept .SELECT specification. This code is also
C     used in DEF_SHIELD, which is also called for esr
      allocate( DOATOM(NUCIND) )
      CALL SETATM(DOATOM,NATOM,-29)

      CALL PRP_ESR_HFC_1(ESRVAL,SWITCH_XY,
     &     WORK(KTMPMAT),WORK(KTMAT),WORK(KTISO),
     &     WORK(KIATM),WORK(KRATM),WORK(KIBUF),
     &     DOATOM,WORK(KCSTRA),WORK(KSCTRA) )
      CALL MEMREL('PRP_ESR',WORK,1,KFRSAV,KFREE,LFREE)
C
      deallocate( DOATOM )
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ESR1 */
      SUBROUTINE PRP_ESR_HFC_1(ESRVAL,SWITCH_XY,TMPMAT,TMAT,TISO,
     &                    IATINF,RATINF,IBUFA,DOATOM,CSTRA,SCTRA)
C***********************************************************************
C
C     Output routine for ESR properties
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
C Used from common blocks:
C  dcbprp: ABUND,?
C
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbesr.h"
#include "dcborb.h"
C
      LOGICAL   SWITCH_XY,TEST,DOATOM(NUCIND)
      REAL*8    ESRVAL(MULTIPESR,MULTIPESR,2,*),
     &          TMPMAT(3,*),TMAT(3,3,NUCDEP),
     &          TISO(NUCDEP),RATINF(3,5,NUCDEP),
     &          CSTRA(NUCDEP*NUCDEP*9),SCTRA(NUCDEP*NUCDEP*9)
      INTEGER   IATINF(3,NUCDEP), IBUFA(3*NUCDEP),IBUFB(3)
#include "nuclei.h"
#include "symmet.h"
C
      CALL HEADER('ESR nuclear hyperfine couplings',1)
C
C
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &'@ name  charge  isotope     mass       abundance    g_N factor'
      CALL PRSYMB(LUPRI,'-',61,1)
      NSPEC = 0
      IATOM = 0
      DO INUC = 1,NUCIND
      IDEG = NUCDEG(INUC)
      IF(DOATOM(INUC)) THEN
        TEST = .FALSE.
        IZ   = IZATOM(INUC)
        NISO = 0
        DO ISO = 1,5
          GVAL = DISOTP(IZ,ISO,'GVAL')
          IF(GVAL.NE.0) THEN
            ABND = DISOTP(IZ,ISO,'ABUNDANCE')
            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = GVAL
            ENDIF
          ENDIF
        ENDDO
        IF(NISO.GT.0) THEN
          IATINF(1,NSPEC+1) = IATOM + 1
          IATINF(2,NSPEC+1) = IZ
          IATINF(3,NSPEC+1) = NISO
          WRITE(LUPRI,'
     &      (A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &       '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
          DO J = 2,NISO
            WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &         '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
          ENDDO
          DO K = 2,IDEG
            IATINF(1,NSPEC+K) = IATOM + K
            IATINF(2,NSPEC+K) = IZ
            IATINF(3,NSPEC+K) = NISO
            DO L = 1,NISO
              RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
              RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
              RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
            ENDDO
          ENDDO
          NSPEC = NSPEC + IDEG
        ENDIF
      ENDIF
      IATOM = IATOM + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        GO TO 9999
      ENDIF
C
      NSCOOR = 3*NUCDEP
      NCMAT  = 9*NUCDEP
      CALL DZERO(TMPMAT,NCMAT)
C
C     *** Extract requested nuclear coordinates and B fields ***
C
      NA = 0
      DO I = 1,NSCOOR
         IF(IPESR(I,1).GT.0) THEN
            NA = NA + 1
            IBUFA(NA) = I
         ENDIF
      END DO
C
      NB = 0
      DO I = 1, 3
         IF(IPESR(I+MXCOOR,1).GT.0) THEN
            NB = NB + 1
            IBUFB(NB) = I
         ENDIF
      END DO
C
C     Fill ESR hfcc matrices Seff_p * TMPMAT * I_N,q
C
C     The factors in front of ESRVAL will cancel the corresponding
C     factors in the Seff_p matrices
C
C
      CALL DZERO(TMPMAT,NCMAT)
C
      IF (MULTIPESR .EQ. 2) THEN
C     ....doublet molecules (Seff = 1/2)
C
         jmhlf = 2
         jphlf = 1
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            IF (SWITCH_XY) THEN
               TMPMAT(IPTAX(2,2),JA) =  2.0D0*ESRVAL(jphlf,jmhlf,1,KA)
               TMPMAT(IPTAX(1,2),JA) = -2.0D0*ESRVAL(jphlf,jmhlf,2,KA)
            ELSE
               TMPMAT(IPTAX(1,2),JA) =  2.0D0*ESRVAL(jphlf,jmhlf,1,KA)
               TMPMAT(IPTAX(2,2),JA) = -2.0D0*ESRVAL(jphlf,jmhlf,2,KA)
            END IF
            TMPMAT(IPTAX(3,2),JA)    =  2.0D0*ESRVAL(jphlf,jphlf,1,KA)
         ENDDO
C
      ELSE IF (MULTIPESR .EQ. 3) THEN
C     .... triplet molecules (Seff = 1)
         jm1 = 3
         j0  = 2
         jp1 = 1
C
         SQRT2 = SQRT(2.0D0)
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            IF (SWITCH_XY) THEN
               TMPMAT(IPTAX(2,2),JA) =  SQRT2*ESRVAL(j0 ,jm1,1,KA)
               TMPMAT(IPTAX(1,2),JA) = -SQRT2*ESRVAL(j0 ,jm1,2,KA)
            ELSE
               TMPMAT(IPTAX(1,2),JA) =  SQRT2*ESRVAL(j0 ,jm1,1,KA)
               TMPMAT(IPTAX(2,2),JA) = -SQRT2*ESRVAL(j0 ,jm1,2,KA)
            END IF
            TMPMAT(IPTAX(3,2),JA)    =        ESRVAL(jp1,jp1,1,KA)
         ENDDO
C
      ELSE IF (MULTIPESR .EQ. 4) THEN
C     .... quartet molecules (Seff = 3/2)
C
C        Take mean value of M_Seff=3/2 and 1/2 values
         fac3hlf = 2.0d0/sqrt(3.0d0) ! fac1hlf is 1.0d0
C
         jm3hlf = 4
         jm1hlf = 3
         jp1hlf = 2
         jp3hlf = 1
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            IF (SWITCH_XY) THEN
               TMPMAT(IPTAX(2,2),JA) =
     &            0.5D0*( ESRVAL(jp3hlf,jp1hlf,1,KA)*fac3hlf
     &                  + ESRVAL(jp1hlf,jm1hlf,1,KA) )
               TMPMAT(IPTAX(1,2),JA) =
     &           -0.5D0*( ESRVAL(jp3hlf,jp1hlf,2,KA)*fac3hlf
     &                  + ESRVAL(jp1hlf,jm1hlf,2,KA) )
            ELSE
               TMPMAT(IPTAX(1,2),JA) =
     &            0.5D0*( ESRVAL(jp3hlf,jp1hlf,1,KA)*fac3hlf
     &                  + ESRVAL(jp1hlf,jm1hlf,1,KA) )
               TMPMAT(IPTAX(2,2),JA) =
     &           -0.5D0*( ESRVAL(jp3hlf,jp1hlf,2,KA)*fac3hlf
     &                  + ESRVAL(jp1hlf,jm1hlf,2,KA) )
            END IF
!           TMPMAT(IPTAX(3,2),JA) = 2.0D0*ESRVAL(jp1hlf,jp1hlf,1,KA)
            TMPMAT(IPTAX(3,2),JA) = ESRVAL(jp3hlf,jp3hlf,1,KA)/3.0D0
     &                            + ESRVAL(jp1hlf,jp1hlf,1,KA)
         ENDDO
      ELSE
         WRITE(LUPRI,'(//A/A)') '@ WARNING:'//
     &   ' Only doublet, triplet, and quartet hfcc implemented yet',
     &   '@ Returning without doing anything ...'
         GO TO 9999
      END IF
C
      IF (IPRPRP.GE.4) THEN
         CALL HEADER('ESR hfcc'//
     &         ' in symmetry coordinates (B~I,au)',1)
         CALL FCPRI (TMPMAT,'SIGMANO',CSTRA,SCTRA)
      ENDIF
C
C     Transform to non-symmetry basis
C
         CALL DZERO(TMAT(1,1,1),NCMAT)
         CALL TRADIP(TMPMAT,TMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C
C     Analyze results
C
         NSPEC = 0
         IATOM = 0
         DO 100 I = 1, NUCIND
            DO 200 ISYMOP = 0, MAXOPR
               IF (IAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(I)) THEN
                     NSPEC = NSPEC + 1
                     CALL ESR_HFC_ANA(TMAT(1,1,IATOM), TISO(IATOM),
     &                 IATINF(1,NSPEC),RATINF(1,1,NSPEC),
     &                 NAMDEP(IATOM),NAMDPX(3*(IATOM-1)+1), IPRPRP)
                  END IF
               END IF
 200        CONTINUE
 100     CONTINUE
C
C
 9999 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck ESR_HFC_ANA */
      SUBROUTINE ESR_HFC_ANA(TMAT,AVETOT,IATINF,RATINF,
     &                  NAME,NAMEX,IPRINT)
C***********************************************************************
C
C     Analyze routine for ESR properties
C     This routine is modelled after SHIANA in DALTON
C
C     April 2000 A. Noerager
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
C
      REAL*8    TMAT(3,3), AVETOT, RATINF(3,*)
      INTEGER   IATINF(3), IPRINT
      CHARACTER NAME*6, NAMEX(3)*6
#include "chrxyz.h"
      REAL*8    TMATNY(3,3), TMATX(3,3), HFC_PVAL(3), HFC_PVEC(3,3)
C
C     factor = 1/c^2 * \mu_N = 1/c^2 * 1/(2 m_p)
      FACTOR = ALPHAC*ALPHAC/(2.0D0*XFMP)
C
C     Multiply with FACTOR
C
      CALL DSCAL(9,FACTOR,TMAT,1)
C
C     Multiply with appropriate g_N factor for each isotope
C     to get final value (in au)
C
      NISO=IATINF(3)

      DO L=1,NISO
         DO J=1,3
           DO I=1,3
              TMATNY(I,J)=RATINF(3,L)*TMAT(I,J)
              TMATX(I,J)=(XTHZ*1.0D-6) * TMATNY(I,J)
           ENDDO
         ENDDO
         AVETOT = (TMATNY(1,1)+TMATNY(2,2)+TMATNY(3,3))/3.0D0
         AVETOTX= (XTHZ*1.0D-6) * AVETOT

C
         CALL TITLER(
     &      'ESR hyperfine coupling constants for '//NAME//':','=',1)
C
         WRITE (LUPRI,'(A,I5//A,F20.5,A)')
     &   '  Isotope number:',NINT(RATINF(1,L)),
     &   '  Isotropic hyperfine coupling constant',AVETOTX, ' MHz'
         IF (IPRINT .GT. 0) WRITE (LUPRI,'(A,1P,D20.10,A)')
     &   '  Isotropic hyperfine coupling constant',AVETOT,  ' au'
C
         CALL HEADER('Total hyperfine coupling matrix (MHz):',1)
         WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
         DO J = 1, 3
            WRITE (LUPRI,'(2X,A6,3F15.5)')
     &             NAMEX(J), (TMATX(I,J),I=1,3)
         END DO

         IF (IPRINT .GT. 0) THEN
            CALL HEADER('Total hyperfine coupling matrix (au):',1)
            WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
            DO J = 1, 3
               WRITE (LUPRI,'(2X,A6,1P,3D20.10)')
     &                NAMEX(J), (TMATNY(I,J),I=1,3)
            END DO
         END IF

C        Calculate principal values from symmetrized A matrix (in MHz)

         TMP = 0.5D0*(TMATX(1,2) + TMATX(2,1))
         TMATX(1,2) = TMP; TMATX(2,1) = TMP
         TMP = 0.5D0*(TMATX(1,3) + TMATX(3,1))
         TMATX(1,3) = TMP; TMATX(3,1) = TMP
         TMP = 0.5D0*(TMATX(2,3) + TMATX(3,2))
         TMATX(2,3) = TMP; TMATX(3,2) = TMP

         CALL RSJACO(3,3,3,TMATX,HFC_PVAL,1,1,0,HFC_PVEC)

         WRITE (LUPRI,'(/A//,(F15.5,A,3F15.5))')
     &      '  Principal values and principal axes'//
     &         ' from symmetrized hyperfine coupling matrix:',
     &      (HFC_PVAL(I),' MHz :',HFC_PVEC(1:3,I), I = 1,3)

      END DO ! L=1,NISO
C
C     END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE  ESR_READ_CMO(CMO,IPRINT,WORK,KFREE,LFREE)
C
C     Hans Joergen Aa. Jensen Feb. 2015
C
C     Read MO coefficents for ESR-CI calculation in this priority order:
C        1. Orbitals on KRMCSCF
C        2. MP2NO on KRMCOLD
C        3. DFCOEF with Hartree-Fock or DFT or MVO orbitals
C
#include "implicit.h"
#include "priunit.h"
      REAL*8  CMO(*), WORK(*)
      INTEGER IPRINT, KFREE, LFREE
      LOGICAL FILEX, FNDLAB
C
C Used from COMMON blocks:
C  dcbgen.h : LUCOEF, LUKRMC
C  dcbham.h : SPINFR; LEVYLE
C  dcbbas.h : N2BBASXQ, ...
C  dcborb.h : NORBT, ..., NCMOTQ
C  dgroup.h : IPQTOQ
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

C
      KFREE_start = KFREE
      CALL DZERO(CMO,N2BBASXQ)

C        1. Orbitals on KRMCSCF ?

      INQUIRE(FILE='KRMCSCF',EXIST=FILEX)
      IF (FILEX) THEN
         CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','PRPESR')

         REWIND (LUKRMC)
         IF (FNDLAB('NEWNATOB', LUKRMC)) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (ESR)  Reading "NEWNATOB" orbitals from file KRMCSCF'
            GO TO 101
         END IF

         REWIND (LUKRMC)
         IF (FNDLAB('MCCINATO', LUKRMC)) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (ESR)  Reading "MCCINATO" orbitals from file KRMCSCF'
            GO TO 101
         END IF

         REWIND (LUKRMC)
         IF (FNDLAB('NEWORB  ', LUKRMC)) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (ESR)  Reading "NEWORB  " orbitals from file KRMCSCF'
            GO TO 101
         END IF

         REWIND (LUKRMC)
         IF (FNDLAB('OLDORB  ', LUKRMC)) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (ESR)  Reading "OLDORB  " orbitals from file KRMCSCF'
            GO TO 101
         END IF

         ! no orbitals found on KRMCSCF, try next file
         GO TO 200

  101    CALL READT(LUKRMC,NCMOTQ,CMO)
         CLOSE(LUKRMC,STATUS='KEEP')
         GO TO 9000
      END IF

C        2. MP2NO on KRMCOLD ?

  200 INQUIRE(FILE='KRMCOLD',EXIST=FILEX)
      IF (FILEX) THEN
         WRITE(LUPRI,'(/A)')
     &      ' (ESR)  Reading MP2 natural orbitals from file KRMCOLD'
         call flshfo(lupri)
         CALL OPNFIL(LUKRMC,'KRMCOLD','OLD','PRPESR')
         CALL MEMGET2('REAL','OCCNO',KOCC ,NORBT   ,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','OCCNAT',KEIG_NAT,NESHT,WORK,KFREE,LFREE)
         CALL MEMGET2('INTE','IBEIG',KBEIG,NORBT   ,WORK,KFREE,LFREE)
C        ... initialize
         CALL DZERO(WORK(KOCC),NORBT)
         CALL DZERO(WORK(KEIG_NAT),NESHT)
         CALL IZERO(WORK(KBEIG),NORBT)
C        ... read from file
         CALL READMP2NAT(CMO,WORK(KOCC),WORK(KEIG_NAT),
     &                   WORK(KBEIG))
         CLOSE(LUKRMC,STATUS='KEEP')

!        WRITE(LUPRI,'(/A)')
!    &      ' (ESR)  MP2 natural orbitals read from file KRMCOLD'
C
         GO TO 9000
      END IF

C     MP2NO not found. Then:
C     3. DFCOEF with Hartree-Fock or DFT or MVO orbitals ?

      CALL MEMGET2('REAL','EIG',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      IF (SPINFR.OR.LEVYLE) THEN
         CALL MEMGET2('INTE','IBPRP',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
         CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,WORK(KIBRP),
     &               TOTERG,10)
      ELSE
         KIBRP = KFREE
         CALL REACMO(LUCOEF,'DFCOEF',CMO,WORK(KEIG),IDUM,
     &               TOTERG,2)
      END IF
      WRITE(LUPRI,'(/A)')
     &   ' (ESR)  molecular orbitals read from CHECKPOINT'
      GO TO 9000

C     DFCOEF not found. Then:
C     4. no more options (yet).

      WRITE(LUPRI,'(//A)') 'Fatal error for ESR: MO coefficients not'//
     &   ' found on KRMCOLD or CHECKPOINT'
      CALL QUIT('ESR fatal error: no MO coefficients found')

 9000 CALL MEMREL('ESR_READ_CMO',WORK,1,KFREE_start,KFREE,LFREE)

      IF ( IPRINT .gt. 20 ) THEN
         CALL HEADER('ESR: input molecular orbitals:',-1)
         DO I = 1,NFSYM
            CALL HEADER('Coefficients for irrep '//CHAR(I+48),-1)
            WRITE(LUPRI,*) ' printing offset ICMOQ(I)+1',ICMOQ(I)+1
            WRITE(LUPRI,*) ' NFBAS(I,0), NORB(I)',NFBAS(I,0),NORB(I)
            CALL PRQMAT(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF

      flush(lupri)
      RETURN
      END
!  -- end of pamesr.F --
