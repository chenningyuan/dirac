#define notESR_DEBUG
!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
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
C     Last revision:
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcborb.h"
#include "orgcom.h"
#include "dgroup.h"
      PARAMETER (NTABLE = 16)
      LOGICAL SET, NEWDEF, RESET, GORG_BY_USER
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.GAUGEO','.PRINT ','.OPERAT','.MULTIP',
     &            '.SPINPO','.SDCI  ','.MAX CI','.CI ROO',
     &            '.THRCVE','.THRCI ','.SDTCI ','.SDT*CI',
     &            '.COREHO','.MY CI ','.THRPCI','.KR-CON'/
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
      MAX_CIESR_IT = 0
      THR_CVEC_ESR = 1.0D-4
      THRPCI_ESR   = 1.0D-3
      THR_CIESR    = 1.0D-6 ! based on 1 ppm = 10**(-6)
C     variables for spin polarization in CI
      NCSH_ESR(:)  = -1
      NVSH_ESR(:)  = -1
      NSSH_ESR(:)  = -1
      LVL_CIESR(:) = 0
      USE_KRAMERS_CONJ = .FALSE.
C
C     Initialize orgcom.h
C     ===================
C
      GORG_BY_USER = .FALSE.
      CALL DZERO(GAGORG,3)
C
C     Process input
C     =========================
C
      NEWDEF = (WORD .EQ. '*ESR   ')
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
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ESRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ESRINP.')
C Read user specified gauge origin
    1          CONTINUE ! .GAUGEO
                  READ(LUCMD,*) (GAGORG(ICOOR), ICOOR = 1,3)
                  GORG_BY_USER = .TRUE.
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
                  DOESR = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
C Read multiplicity
    4          CONTINUE ! .MULTIPLICITY
                  READ(LUCMD,*) MULTIPESR
               GO TO 100
C&&&& SPINPOLARIZATION (= SCI)
    5          CONTINUE ! .SPINPOL
                  LVL_CIESR(1:2) = 1
                  READ(LUCMD,*) (NVSH_ESR(I),I=1,NFSYM)
                  READ(LUCMD,*) (NSSH_ESR(I),I=1,NFSYM)
               GO TO 100
C&&&& Singles and doubles CI
    6          CONTINUE ! .SDCI
                  LVL_CIESR(1:2) = 2
                  READ(LUCMD,*) (NVSH_ESR(I),I=1,NFSYM)
                  READ(LUCMD,*) (NSSH_ESR(I),I=1,NFSYM)
               GO TO 100
C&&&& MAX CI ITERATIONS
    7          CONTINUE ! .MAX CI ITERATIONS
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
                  LVL_CIESR(1:2) = 3
                  READ(LUCMD,*) (NVSH_ESR(I),I=1,NFSYM)
                  READ(LUCMD,*) (NSSH_ESR(I),I=1,NFSYM)
               GO TO 100
C&&&& MR-CI-SDT/SD
   12          CONTINUE ! .SDT*CI
                  LVL_CIESR(1) = 3
                  LVL_CIESR(2) = 2
                  READ(LUCMD,*) (NVSH_ESR(I),I=1,NFSYM)
                  READ(LUCMD,*) (NSSH_ESR(I),I=1,NFSYM)
               GO TO 100
C&&&& number of holes in non-frozen NISH(1:2) - NVSH_ESR(1:2) orbitals
   13          CONTINUE ! .COREHOLES
                  READ(LUCMD,*) LVL_CIESR(0)
               GO TO 100
C&&&& flexible input of the CI you want
   14          CONTINUE ! .MY CI
                  READ(LUCMD,*) LVL_CIESR(0:2)
                  READ(LUCMD,*) (NCSH_ESR(I),I=1,NFSYM)
                  READ(LUCMD,*) (NVSH_ESR(I),I=1,NFSYM)
                  READ(LUCMD,*) (NSSH_ESR(I),I=1,NFSYM)
               GO TO 100
C&&&& limit for printing CI coefficients in output
   15          CONTINUE ! .THRPCI
                  READ(LUCMD,*) THRPCI_ESR
               GO TO 100
C&&&& use Kramers conjugation on CI vector
   16          CONTINUE ! .KR_CON
                  USE_KRAMERS_CONJ = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized in ESRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ESRINP.')
            END IF
      END IF
  300 CONTINUE
      LVL_CIESR_T = LVL_CIESR(0) + LVL_CIESR(1) + LVL_CIESR(2)
      IF (LVL_CIESR_T == 0) THEN
         NCSH_ESR(:) = 0
         NVSH_ESR(:) = 0
         NSSH_ESR(:) = 0
      ELSE
         IF (NSSH_ESR(1) < 0) THEN
            NSSH_ESR(1:2) = NSSH(1:2)
         END IF
         IF (NCSH_ESR(1) < 0) THEN
            NCSH_ESR(1:2) = 0
            IF (NVSH_ESR(1) < 0) THEN
               NVSH_ESR(1:2) = NISH(1:2) - NFRO(1:2)
            ELSE IF (LVL_CIESR(0) > 0) THEN
               NCSH_ESR(1:2) = NISH(1:2) - NFRO(1:2) - NVSH_ESR(1:2)
            END IF
         ELSE
            IF (NVSH_ESR(1) < 0) THEN
               NVSH_ESR(1:2) = NISH(1:2) - NFRO(1:2) - NCSH_ESR(1:2)
            END IF
         END IF
      END IF
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
      IF (LVL_CIESR_T .gt. 0 .OR. MULTIPESR .ne. 2) THEN
      !  LVL_CIESR_T test because o.k. with no CI if just doublet one open-shell DHF
         IF (N_CIESR < 0) N_CIESR = MULTIPESR
         WRITE(LUPRI,'(/A,I5)')
     &   '* Number of roots to converge in CI:', N_CIESR
         IF (N_CIESR .lt. MULTIPESR) THEN
            WRITE(LUPRI,*) 'INPUT ERROR: .MULTIP > .CI ROOTS'
            INPERR = INPERR + 1
         END IF
      END IF
      IF (LVL_CIESR_T .GT. 0) THEN
         WRITE (LUPRI,'(/A)')
     &      ' * Spin polarization and correlation included with'

         IF (LVL_CIESR(0) .LT. 0) THEN
            WRITE(LUPRI,*) 'INPUT ERROR: number of core holes < 0'
            INPERR = INPER + 1
         ELSE IF (LVL_CIESR(0) .EQ. 0) THEN
         ELSE
            WRITE (LUPRI,'(A,I2,A,T50,2I5)') '   - 0 to',LVL_CIESR(0),
     &         ' holes in non-frozen core orbitals:', NCSH_ESR(1:NFSYM)
         END IF

         IF (LVL_CIESR(1) .LT. 0) THEN
            WRITE(LUPRI,*) 'INPUT ERROR: number of valence holes < 0'
            INPERR = INPER + 1
         ELSE IF (LVL_CIESR(1) .EQ. 0) THEN
         ELSE
            WRITE (LUPRI,'(A,I2,A,T50,2I5)') '   - 0 to',LVL_CIESR(1),
     &         ' holes in valence orbitals:', (NVSH_ESR(I),I=1,NFSYM)
         END IF

         WRITE (LUPRI,'(A,T50,2I5)') '   - active/open shell orbitals:',
     &      (NASH(I),I=1,NFSYM)

         IF (LVL_CIESR(2) .LT. 0) THEN
            WRITE(LUPRI,*)
     &         'INPUT ERROR: number of virtual space electrons < 0'
            INPERR = INPER + 1
         ELSE IF (LVL_CIESR(2) .EQ. 0) THEN
         ELSE
            WRITE (LUPRI,'(A,I2,A,T50,2I5)') '   - 0 to',LVL_CIESR(2),
     &        ' electrons in virtual orbitals:', (NSSH_ESR(I),I=1,NFSYM)
         END IF
      END IF

      IF (GORG_BY_USER) THEN
         CALL PRSYMB(LUPRI,'-',75,0)
         WRITE(LUPRI,'(/A,3F12.6)') '* User defined magnetic '//
     &           'gauge origin used:',(GAGORG(I), I = 1, 3)
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
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
C     Calculate expectation values defined in /CBIESR/
C
C     Written by Anette Noerager - Mar 2000
C
C***********************************************************************
       use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C     for PRQMAT print of complex matrices:
      PARAMETER (NZCMPLX = 2)
      DIMENSION IQINDX(NZCMPLX)
      SAVE IQINDX
      DATA IQINDX /1,2/
C
      DIMENSION ESRVAL(MULTIPESR,MULTIPESR,NZCMPLX,*),WORK(*)
      CHARACTER WF*4
      LOGICAL   FILEX
      DIMENSION GTENSOR(3,3),GCORR(3,3),GMATRIX(3,3),GEIG(3),GSHIFT(3),
     &          GEVEC(3,3), GCOR(3),JMS(-5:5)
C
      PARAMETER (D2=2.0D0, D3 = 3.0D0, D4 = 4.0D0)
C
C Used from common blocks:
C  dcbgen : LUCOEF
C
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbesr.h"
#include "gfac.h"
#include "dcbprp.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dcbopt.h"
C
      CALL QENTER('PRPESR')
      MZ = MIN(2,NZ)
!d    write (lupri,*) 'DEBUG 1 NZ, MZ',NZ,MZ
C
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
C
C **  Retrieve molecular orbitals
C
      CALL MEMGET2('REAL','CMO',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL ESR_READ_CMO(WORK(KCMO),WORK,KFREE,LFREE)

!     write (lupri,*) 'hj1 N_CIESR MAX_CIESR_IT NZCONF',
!    &                     N_CIESR,MAX_CIESR_IT,NZCONF
      IF (N_CIESR .GE. 0) THEN
         !
      IF (MULTIPESR .gt. NZCONF) THEN
         CALL QUIT('.MULTIPESR specification > # determinants')
      END IF
      IF (MAX_CIESR_IT .GT. 0) THEN
         ! Davidson iterations
!        N_CIESR = MIN(NZCONF, MULTIPESR + 1)
         N_CIESR = MULTIPESR
!      write(lupri,*) 'N_CIESR,MULTIPESR is', N_CIESR,MULTIPESR
      ELSE
         ! full CI (now primarily for debugging Davidson iterations code)
         N_CIESR = NZCONF
      END IF
      END IF
!d    write (lupri,*) 'hj2 N_CIESR MAX_CIESR_IT NZCONF',
!d   &                     N_CIESR,MAX_CIESR_IT,NZCONF
C
C **  Solve ESR CI problem
C     on exit: WORK(KEIGVL) - N_CIESR CI eigenvalues
C              WORK(KEIGVC) - N_CIESR CI eigenvectors (NZCONF,N_CIESR,MZ)
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
     &   MZ*MULTIPESR*MULTIPESR,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHCI),MZ*MULTIPESR*MULTIPESR)
      DO I = 1, MULTIPESR
         II = KHCI + (I-1)*MULTIPESR + (I-1)
         WORK(II) = WORK(KEIGVL-1 + I)
      END DO

      IF (IPRESR .GT. 3) THEN
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
C
#ifdef ESR_DEBUG
       WRITE(LUPRI,*)
     &   'PRPESR: Real part of the N_CIESR CI vectors'
       CALL OUTPUT(WORK(KEIGVC),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       IF (MZ .EQ. 2) THEN
          WRITE(LUPRI,*)
     &      'PRPESR: Imag part of the N_CIESR CI vectors'
          KEIGVC_I = KEIGVC + NZCONF*N_CIESR
          CALL OUTPUT(WORK(KEIGVC_I),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       END IF
#endif
      CALL ESRPHASE(ESRVAL,WORK(KCMO), WORK(KPMAT),
     &              WORK(KHCI),WORK(KEIGVC),MZ,WORK,KFREE,LFREE)
C     CALL ESRPHASE(ESRVAL,CMO,PMAT,
C    &              HCI,CIEVEC,MZ,WORK,KFREE,LFREE)
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
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      END IF
C
      CALL PRPESR1(ESRVAL,WORK(KESRST),WORK(KCMO),WORK(KPMAT),
     &             WORK(KEIGVC),MZ,WORK,KFREE,LFREE)
C     CALL PRPESR1(ESRVAL,LESRST,CMO,PMAT,
C    &             CIEVEC,MZ,WORK,KFREE,LFREE)
C
      IF (IPRESR .GT. 1) THEN
        CALL HEADER('PRPESR results in CI eigenvector basis'//
     &           ' (test output)',1)
        DO I = 1,NESRP
          INDXPR = LESRP(I)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
cq        WRITE(LUPRI,'(4X,A,4F15.6)')
cq   &       'Quaternion value:',(ESRVAL(J,I),J=1,4)
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
        jmhlf = 1
        jphlf = 2
        DO K = 1,3
C       ... for B field in direction k (x,y,z)
          I = IPESR(MXCOOR+K,1)
          GTENSOR(K,1) =  D4*ESRVAL(jphlf,jmhlf,1,I)
          GTENSOR(K,2) = -D4*ESRVAL(jphlf,jmhlf,2,I)
          GTENSOR(K,3) =  D4*ESRVAL(jphlf,jphlf,1,I)
          J = IPESR(MXCOOR+K,2)
          GCORR(K,1)   =  D4*ESRVAL(jphlf,jmhlf,1,J)
          GCORR(K,2)   = -D4*ESRVAL(jphlf,jmhlf,2,J)
          GCORR(K,3)   =  D4*ESRVAL(jphlf,jphlf,1,J)
        ENDDO
C
      ELSE IF (MULTIPESR .EQ. 3) THEN
C     .... triplet molecules (S_eff = 1)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
C  BRANCH CODE
        jm1 = 1
        j0  = 2
        jp1 = 3
        D2SQRT2 = D2*SQRT(D2)
        DO K = 1,3
C       ... for B field in direction k (x,y,z)
          I = IPESR(MXCOOR+K,1)
          GTENSOR(K,1) = D2SQRT2*ESRVAL(j0 ,jm1,1,I)
          GTENSOR(K,2) = D2SQRT2*ESRVAL(j0 ,jm1,2,I)
          GTENSOR(K,3) =      D2*ESRVAL(jp1,jp1,1,I)
c??       GTENSOR(K,3) = ESRVAL(jp1,jp1,1,I) - ESRVAL(jm1,jm1,1,I)
c         WRITE(lupri,*)'K =', K
c         WRITE(lupri,*)'GTENSOR(K,1)', GTENSOR(K,1),
c    &                  'GTENSOR(K,2)', GTENSOR(K,2),
c    &                  'GTENSOR(K,3)', GTENSOR(K,3)
          J = IPESR(MXCOOR+K,2)
          GCORR(K,1)   = D2SQRT2*ESRVAL(j0 ,jm1,1,J)
          GCORR(K,2)   = D2SQRT2*ESRVAL(j0 ,jm1,2,J)
          GCORR(K,3)   =      D2*ESRVAL(jp1,jp1,1,J)
c         WRITE(lupri,*)'K =', K
c         WRITE(lupri,*)'GCORR(K,1)', GCORR(K,1),
c    &      'GCORR(K,2)', GCORR(K,2),
c    &      'GCORR(K,3)', GCORR(K,3)
        ENDDO
C   END BRANCH CODE
C  TRUNK CODE
C        Find Sigma_Z expectation values for the three states:
C         I = IPSIGMA(3)
C         WRITE (LUPRI,*) 'Analysis to find triplet components'
C         JMS(-1) = -1
C         JMS( 0) = -1
C         JMS( 1) = -1
C         DO K = 1,3
C            SZK  = ESRVAL(K,K,1,I)
C            ISZK = NINT(SZK)
C            WRITE (LUPRI,*)'State, <Sigma_Z>, S_Z(eff)',K,SZK,ISZK
C            IF (ISZK .GT. 1 .OR. ISZK .LT. -1) THEN
C               CALL QUIT('error: S_Z(eff) too big to triplet')
C            ELSE IF (JMS(ISZK) .NE. -1) THEN
C               CALL QUIT('error: same S_Z(eff) twice!!!')
C            END IF
C            JMS(ISZK) = K
C         END DO
C         jm1 = jMS(-1)
C         j0  = jMS( 0)
C         jp1 = jMS( 1)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
C      D2SQRT2 = D2*SQRT(D2)
C      DO K = 1,3
C        ... for B field in direction k (x,y,z)
C         I = IPESR(MXCOOR+K,1)
C         GTENSOR(K,1) = -D2SQRT2*ESRVAL(jp1,j0 ,1,I)
C         GTENSOR(K,2) =  D2SQRT2*ESRVAL(jm1,j0 ,2,I)
C         GTENSOR(K,3) =  D2*ESRVAL(jp1,jp1,1,I)
c??      GTENSOR(K,3) =  ESRVAL(jp1,jp1,1,I) - ESRVAL(jm1,jm1,1,I)
C         J = IPESR(MXCOOR+K,2)
C         GCORR(K,1)   =  -D2SQRT2*ESRVAL(jp1,j0,1,J)
c         GCORR(K,2)   =   D2SQRT2*ESRVAL(jm1,j0,2,J)
c         GCORR(K,3)   =   D2*ESRVAL(jp1,jp1,1,J)
c      ENDDO

C  END TRUNK CODE
C
C     ***************************************
C     Zero-field-splitting parameters D and E
C     ***************************************
C
C  BRANCH CODE
        EM1 = WORK(KHCI)
        E0  = WORK(KHCI+  MULTIPESR+1)
        EP1 = WORK(KHCI+2*MULTIPESR+2)
        ZFS_D = -(E0-(EP1+EM1)/D2)
        ZFS_E = (EP1-EM1)/D2
        ZFS_DCM = ZFS_D*XTKAYS
        ZFS_ECM = ZFS_E*XTKAYS

        CALL HEADER(' Zero field splitting parameters : ',1)
        WRITE (LUPRI,'(2(A,F25.8,A/))')
     &   '  ZFS parameter D:      ',ZFS_DCM, ' cm-1 ',
     &   '  ZFS parameter E:      ',ZFS_ECM, ' cm-1 '
C   END BRANCH CODE
C
C   TRUNK CODE
C      ZFS_D = -(WORK(KEIGVL) - (WORK(KEIGVL+1)+WORK(KEIGVL+2))/D2)
C..   ZFS_D = EIGVAL(1) - (EIGVAL(2)+EIGVAL(3))/D2
C      ZFS_E = (WORK(KEIGVL+1)-WORK(KEIGVL+2))/D2
C      ZFS_DCM = ZFS_D*XTKAYS
C      ZFS_ECM = ZFS_E*XTKAYS
C      CALL HEADER('Zero field splitting parameters:',1)
C         WRITE (LUPRI,'(5(/A,F15.8,A))')
C     &   ' ZFS parameter D:                ',ZFS_DCM,  '(cm-1) '
C         WRITE (LUPRI,'(5(/A,F15.8,A))')
C     &   ' ZFS parameter E:                ',ZFS_ECM,  '(cm-1) '
C     print *, 'ZERO-FIELD D= ',ZFS_DCM, '(cm-1)'
C   END TRUNK CODE

C
      ELSE IF (MULTIPESR .EQ. 4) THEN
C     .... quartet molecules (S_eff = 3/2)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
        DO K = 1,3
C       ... for B field in direction k (x,y,z)
          I = IPESR(MXCOOR+K,1)
          GTENSOR(K,1) =  D2*ESRVAL(2,3,1,I)
          GTENSOR(K,2) = -D2*ESRVAL(2,3,2,I)
          GTENSOR(K,3) =  D4*ESRVAL(3,3,1,I)
c??       GTENSOR(K,3) = ESRVAL(jp1,jp1,1,I) - ESRVAL(jm1,jm1,1,I)
          J = IPESR(MXCOOR+K,2)
          GCORR(K,1)   =  D2*ESRVAL(2,3,1,J)
          GCORR(K,2)   = -D2*ESRVAL(2,3,2,J)
          GCORR(K,3)   =  D4*ESRVAL(3,3,1,J)
        ENDDO
C
      ELSE
         call quit('Only doublet, triplet, and quartet implemented yet')
      END IF
C
      CALL HEADER(' ESR g shifts : ',1)
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

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC,WORK(KFREE),LFREE)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - 2.0D0)*1.D6
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,'(/A//A//5X,4F18.8//A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The Dirac equation g values corresponding to g_e = 2',
     &   '  The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   '  The g shifts in ppm (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO,
     &   ' The corresponding eigenvectors:',
     &   ((GEVEC(I,K),K=1,3),I=1,3)
C
C     (2) Get the corrected g values by scaling
C
      DO I = 1,3
        GCOR(I) = (GEIG(I)/2)*GFAC
        GSHIFT(I) = (GCOR(I) - GFAC)*1.D6
      ENDDO
      GISO = (GCOR(1) + GCOR(2) + GCOR(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,'(/A,F13.10,A//A//5X,4F18.8//A//5X,4F18.1))')
     &   '* The corrected g values scaled with',GFAC,' / 2',
     &   '  The corrected g values (x, y, z, iso):',
     &   (GCOR(I),I =1,3), GISO,
     &   '  The corrected g shifts in ppm (x, y, z, iso):',
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

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC,WORK(KFREE),LFREE)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - GFAC)*1.D6
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,
     &   '(/A,F13.10/A//5X,4F18.8//A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The g_1 corrected Dirac eq. g values corresponding to g_e ='
     &   ,GFAC,
     &   '  The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   '  The g shifts in ppm (x, y, z, iso):',
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
     &   CALL PRP_ESR(ESRVAL,WORK,KFREE,LFREE)
C
C     Memory deallocation
C
      CALL MEMREL('PRPESR',WORK,1,KFRSAV,KFREE,LFREE)
 9999 CALL QEXIT('PRPESR')
      RETURN
C
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPESR: Coefficient file not found !'
      CALL QUIT('PRPESR: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ESRPHASE(ESRVAL,CMO,PMAT,HCI,CIEVEC,MZ,
     &                    WORK,KFREE,LFREE)
C*****************************************************************************
C
C     --- Find state vectors diagonalizing S_z (for a "standard" S_eff basis)
C         Modify phases acc. to convention so <M_S + 1 | S_x | M_S > .gt. 0
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C     for PRQMAT print of complex matrices:
      PARAMETER (NZCMPLX = 2)
      DIMENSION IQINDX(NZCMPLX)
      SAVE IQINDX
      DATA IQINDX /1,2/
C
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C
      LOGICAL   REAL_PMLF
      DIMENSION ESRVAL(MULTIPESR,MULTIPESR,2,4),CMO(*),PMAT(*),
     &          HCI(MULTIPESR,MULTIPESR,MZ), CIEVEC(NZCONF,N_CIESR,MZ),
     &          WORK(*)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   IPCMO(2)
C
C Used from common blocks:
C  dcbgen : LU1INT,?
C  dcbopt : KZCONF,NZCONF
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
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &      FILE='AOPROPER')
C
C     Initialize
C     ==========
C     (using first NPSIGMA elements of ESRVAL for IPSIGMA properties,
C     later reused for the NESRP requested properties)
C
      IF (N_CIESR .LT. MULTIPESR) THEN
         WRITE(LUPRI,*) 'ERROR: N_CIESR.lt.MULTIPESR',N_CIESR,MULTIPESR
         CALL QUIT('N_CIESR .lt. MULTIPESR')
      END IF
      IF (NPSIGMA.NE. 3)      CALL QUIT('NPSIGMA .ne. 3')
      N2MULT  = MULTIPESR ** 2
      CALL DZERO(ESRVAL,N2MULT*2*NPSIGMA)
C
      I = IPSIGMA(3)
      INDXPR = LESRP(I)

#ifdef ESR_DEBUG
       WRITE(LUPRI,*)
     &   'ESRPHASE: Real part of the N_CIESR CI vectors'
       CALL OUTPUT(CIEVEC(1,1,1),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       IF (MZ .EQ. 2) THEN
          WRITE(LUPRI,*)
     &      'ESRPHASE: Imag part of the N_CIESR CI vectors'
          CALL OUTPUT(CIEVEC(1,1,2),1,NZCONF,1,N_CIESR,
     &             NZCONF,N_CIESR,-1,LUPRI)
       END IF
#endif
      CALL ESRTRANS(INDXPR,ESRVAL,CMO,PMAT,CIEVEC,MZ,
     &              WORK,KFREE,LFREE)
C       temp. use ESRVAL(...,2:3) for eigenvalues and eigenvectors
C       note: using MZ in QDIAG and not 2 because we know Sigma_z
C             is a real operator
      CALL QDIAG(MZ,MULTIPESR,ESRVAL,MULTIPESR,MULTIPESR,
     &           ESRVAL(1,1,1,2),1,ESRVAL(1,1,1,3),
     &           MULTIPESR,MULTIPESR,WORK(KFREE),LFREE,IERR)

      IF (IPRESR .GT. 4) THEN
        DO I = 2,3
          write (lupri,'(A,I2/A)')
     &     ' test output after qdiag for S_z; ESRVAL matrix',I,
     &     ' (2 ~ eigenvalues, 3 ~ eigenvectors)'
          CALL PRQMAT(ESRVAL(1,1,1,I),MULTIPESR,MULTIPESR,
     &                MULTIPESR,MULTIPESR,MZ,IQINDX,LUPRI)
        END DO
      END IF

C
C     Transform the first HCI matrix to the new CI basis
C
      CALL DCOPY(N2MULT*MZ,HCI,1,ESRVAL,1)
      CALL QTRANS('AOMO','S',D0,MULTIPESR,MULTIPESR,MULTIPESR,MULTIPESR,
     &        ESRVAL,         MULTIPESR,MULTIPESR,MZ,IQINDX,
     &        HCI,            MULTIPESR,MULTIPESR,MZ,IQINDX,
     &        ESRVAL(1,1,1,3),MULTIPESR,MULTIPESR,MZ,IQINDX,
     &        ESRVAL(1,1,1,3),MULTIPESR,MULTIPESR,MZ,IQINDX,
     &        WORK(KFREE),LFREE,-1)
      IF (IPRESR .GT. 3) THEN
        write (lupri,*) ' H_CI matrix after transform. to S_z basis'
        CALL PRQMAT(HCI,MULTIPESR,MULTIPESR,
     &              MULTIPESR,MULTIPESR,MZ,IQINDX,LUPRI)
      END IF
C
C     Transform the first MULTIPESR CI eigenvectors to S_z basis
C
      L_CIVEC2 = NZCONF * N_CIESR * MZ
      CALL MEMGET2('REAL','CIVEC2',KCIVEC2 ,L_CIVEC2,WORK,KFREE,LFREE)
      CALL DCOPY(L_CIVEC2,CIEVEC,1,WORK(KCIVEC2),1)
      CALL QGEMM(NZCONF,MULTIPESR,MULTIPESR,D1,
     &           'N','N',IQINDX,
     &           WORK(KCIVEC2),NZCONF,N_CIESR,MZ,
     &           'N','N',IQINDX,
     &           ESRVAL(1,1,1,3),MULTIPESR,MULTIPESR,MZ,
     &           D0,IQINDX, CIEVEC,NZCONF,N_CIESR,MZ)
      CALL MEMREL('ESRPHASE-1',WORK,1,KCIVEC2,KFREE,LFREE)
      IF (IPRESR .GT. 10) THEN
        write (lupri,*) ' CIEVEC after transform. to S_z basis'
        CALL PRQMAT(CIEVEC,NZCONF,MULTIPESR,
     &              NZCONF,N_CIESR,MZ,IQINDX,LUPRI)
      END IF
C
C     Calculate S_x, S_y, S_z matrices in CIEVEC basis,
C     using ESRVAL(...,1:3)
C
        DO I = 1,NPSIGMA
          I1 = IPSIGMA(I)
          INDXPR = LESRP(I1)
C
C         =============================
C
          IF (IPRESR.GT.1) THEN
            CALL PRSYMB(LUPRI,'=',75,0)
            WRITE(LUPRI,'(A,I3,A,3X,A16)')
     &        '* ESRPHASE: calculating property no.',I,' :',
     &        PRPNAM(INDXPR)
            CALL PRSYMB(LUPRI,'.',75,0)
            CALL WRIXPR(I,INDXPR)
          END IF
C
          CALL ESRTRANS(INDXPR,ESRVAL(1,1,1,I),CMO,PMAT,
     &                  CIEVEC,MZ,WORK,KFREE,LFREE)
        ENDDO
C
      CLOSE(LU1INT,STATUS = 'KEEP')
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
         K   = NINT(SZJ+SEFF+D1)
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
C       Use Sigma_X to determine phase of J.gt.1
C       using the requirement that <J+1|Seff_+|J> .gt. 0
C       (which is same as <J+1|Sigma_X|J> .gt. 0) /tow+hjaaj
C
      PHASE_JM1_R = D1
      PHASE_JM1_I = D0
      DO J = 2,MULTIPESR
        JM1 = J-1
        WRITE (LUPRI,'(/A,2I5,F15.10,A,F15.10)')
     &     " J,JM1,<J|S_x|JM1>  :",J,JM1,ESRVAL(J,JM1,1,1),
     &                            ' + i',ESRVAL(J,JM1,2,1)
        WRITE (LUPRI,'(/A,2I5,F15.10,A,F15.10)')
     &     " J,JM1,<J|S_y|JM1>  :",J,JM1,ESRVAL(J,JM1,1,2),
     &                            ' + i',ESRVAL(J,JM1,2,2)
        SPLUS_R = ESRVAL(J,JM1,1,1) - ESRVAL(J,JM1,2,2)
        SPLUS_I = ESRVAL(J,JM1,2,1) + ESRVAL(J,JM1,1,2)
        WRITE (LUPRI,'(/A,2I5,F15.10,A,F15.10)')
     &     " J,JM1,<J|S_+|JM1>  :",J,JM1,SPLUS_R,
     &                            ' + i',SPLUS_I
        PHASE_J_R = PHASE_JM1_R*SPLUS_R
     &            - PHASE_JM1_I*SPLUS_I
        PHASE_J_I = PHASE_JM1_R*SPLUS_I
     &            + PHASE_JM1_I*SPLUS_R
C       ... here PHASE is new phase of vector J-1
        WRITE (LUPRI,'(/A,2I5,F15.10,A,F15.10)')
     &     " J,JM1,<J|S_x|JM1'> :",J,JM1,PHASE_J_R,' + i',PHASE_J_I
        PHASE_J_NORM = SQRT(PHASE_J_R**2 + PHASE_J_I**2)
        PHASE_J_R  = PHASE_J_R / PHASE_J_NORM
        PHASE_J_I  = PHASE_J_I / PHASE_J_NORM
           Seff_M = (JM1 - 1) - Seff
           CPLUS = SQRT( SEFF*(SEFF+1.0D0)-Seff_M*(Seff_M+1.0D0) ) / 2
           write(lupri,*) 'CPLUS / 2 ', CPLUS
           WRITE (LUPRI,'(/A,2I5,F15.10,A,F15.10,A,2F15.10)')
     &     " J,JM1,<J|S_x|JM1'> :",J,JM1,PHASE_J_R,' + i',PHASE_J_I,
     &     '; PHASE_J_NORM',PHASE_J_NORM, PHASE_J_NORM - CPLUS
C
        IF (MZ .EQ. 2) THEN
           DO M = 1, NZCONF
              C_M_R = CIEVEC(M,J,1)
              C_M_I = CIEVEC(M,J,2)
              CIEVEC(M,J,1) = PHASE_J_R*C_M_R - PHASE_J_I*C_M_I
              CIEVEC(M,J,2) = PHASE_J_R*C_M_I + PHASE_J_I*C_M_R
           END DO
           PHASE_JM1_R = PHASE_J_R
           PHASE_JM1_I = PHASE_J_I
        ELSE
           PHASCHK = PHASE_JM1_R * ESRVAL(J,JM1,1,1)
           IF (PHASCHK .LT. D0) THEN
             PHASE_JM1_R = -PHASE_JM1_R
C         ... here PHASE is updated to new phase of vector J
             WRITE (LUPRI,'(A)')
     &         ' Wrong phase, we scale this eigenvector with -1 ...'
             DO IZ = 1,MZ
               CALL DSCAL(NZCONF,DM1,CIEVEC(1,J,IZ),1)
             END DO
           END IF
        END IF
      END DO
C
C **  Form S**2 matrix in ESRVAL(...,4) in CIEVEC basis
C
      FAC  = D0
      DO I = 1,3
        CALL QGEMM(MULTIPESR,MULTIPESR,MULTIPESR,D1,
     &           'N','N',IQINDX, ESRVAL(1,1,1,I),MULTIPESR,MULTIPESR,2,
     &           'N','N',IQINDX, ESRVAL(1,1,1,I),MULTIPESR,MULTIPESR,2,
     &           FAC,    IQINDX, ESRVAL(1,1,1,4),MULTIPESR,MULTIPESR,2)
        FAC = D1
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
      SUBROUTINE PRPESR1(ESRVAL,LESRST,CMO,PMAT,
     &                   CIEVEC,MZ,WORK,KFREE,LFREE)
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
      PARAMETER(D0 = 0.0D0)
C
      LOGICAL LESRST(2,*),REAL_PMLF
      DIMENSION ESRVAL(MULTIPESR,MULTIPESR,2,*),CMO(*),PMAT(*),
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

      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &      FILE='AOPROPER')
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
C      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPESR')
C      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6,
C     &            WORK,KFREE,LFREE)
C      CLOSE(LUCOEF,STATUS='KEEP')

C     Expectation values normal procedure
C    ************************************
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
C        MZ = MIN(2,NZ)
C        CALL MEMGET('REAL',KPMLF  ,MZ*(2*NASHT)**2,WORK,KFREE,LFREE)
C        CALL MEMGET('REAL',KMUUUU ,0              ,WORK,KFREE,LFREE)
C        LPCI = N2CIESR
C        IF (MZ .EQ. 2) LPCI = 2*LPCI
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
             WRITE(LUPRI,'(A,I3,A,A16/A,3I5)')
     &       '* PRPESR1: calculating property no.',I,' : ',
     &         PRPNAM(INDXPR),
     &       '  Boson symmetry, g/u symmetry, time reversal :',
     &         ISYM,IOPSY,ITIM
          END IF
C
          CALL ESRTRANS(INDXPR,ESRVAL(1,1,1,I),CMO,PMAT,
     &                  CIEVEC,MZ,WORK,KFREE,LFREE)
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
C          CALL PRPMAT(INDXPR,IOPSY,PMAT,CMO,IBEIG,IPCMO,NASH,
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
C             IF (MZ .EQ. 2) THEN
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
C             IF (MZ .EQ. 2) THEN
C                WRITE (LUPRI,'(/A)')
C     &          ' Imag part of det-CI property matrix'
C                CALL OUTPUT(WORK(KPI),1,N_CIESR,1,N_CIESR,
C     &              N_CIESR,N_CIESR,1,LUPRI)
C             END IF
c         END IF
C
C         Transform to eigenvector basis
C
C          IF (MZ .EQ. 2) THEN
C             JZ = 1
C          ELSE IF (REAL_PMLF) THEN
C             JZ = 1
Cunoedvendig CALL DZERO(ESRVAL(1,1,2,I),N2CIESR)
C          ELSE
C             JZ = 2
Cunoedvendig CALL DZERO(ESRVAL(1,1,1,I),N2CIESR)
C          END IF
C
C          CALL QTRANS('AOMO','S',D0,N_CIESR,N_CIESR,N_CIESR,N_CIESR,
C     &        WORK(KPR),N_CIESR,N_CIESR,MZ,IQINDX(JZ),
C     &        ESRVAL(1,1,1,I),N_CIESR,N_CIESR,2,IQINDX,
C     &        CIEVEC,N_CIESR,N_CIESR,MZ,IQINDX,
C     &        CIEVEC,N_CIESR,N_CIESR,MZ,IQINDX,
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
c         WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
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
      CLOSE(LU1INT,STATUS = 'KEEP')
      call memchk('prpesr1 4',work,1)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ESRTRANS(INDXPR,ESRVAL,CMO,PMAT,CIEVEC,MZ,
     &                    WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate one-electron property INDXPR in CI basis (in CIEVEC)
C     and return result in ESRVAL.
C
C     Input:
C       INDXPR : index of property operator
C       CMO    : MO coefficients
C       CIEVEC : CI vectors
C       MZ     : = 1 if CI vectors real, =2 if complex
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
      PARAMETER (NZCMPLX = 2)
      DIMENSION IQINDX(NZCMPLX)
      SAVE IQINDX
      DATA IQINDX /1,2/
C
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C
      LOGICAL   REAL_PMLF
      DIMENSION ESRVAL(MULTIPESR,MULTIPESR,2),CMO(*),PMAT(*),
     &          CIEVEC(NZCONF,N_CIESR,*), WORK(*)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   IPCMO(2)
C
C Used from common blocks:
C  dcbgen : LU1INT,?
C  dcbopt : KZCONF,NZCONF
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
!d    print *,'ESRTRANS 1 KFREE,LFREE',KFREE,LFREE
!d    call memchk('esrTRANS 1',work,1)
C
c LU1INT must be opened before calling ESRTRANS
c     OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
c    &      FILE='AOPROPER')
C
C     Initialize
C     ==========
C
      IF (N_CIESR .lt. MULTIPESR) CALL QUIT('N_CIESR .Lt. MULTIPESR')
C
      CALL DZERO(ESRVAL,(MULTIPESR**2)*2)
C     Allocate memory for property matrix in molfdir format:
      CALL MEMGET2('REAL','PMLF ',
     &   KPMLF  ,MZ*(2*NASHT)**2,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','MUUUU',
     &   KMUUUU ,0,WORK,KFREE,LFREE)
C
C     Extract symmetry information
C     =============================
C
      ISYM = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      ITIM = IPRPTIM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
chj   IF (IPRESR.GT.1 .OR. ITIM.GE.0) THEN
      IF (IPRESR.GT.1) THEN
         CALL PRSYMB(LUPRI,'-',75,0)
         WRITE(LUPRI,'(A,A16/A,3I5)')
     &       '* ESRTRANS: calculating property  : ',
     &         PRPNAM(INDXPR),
     &       '  Boson symmetry, g/u symmetry, time reversal :',
     &         ISYM,IOPSY,ITIM
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
!d    call memchk('esrTRANS 2',work,1)
      CALL MEMGET2('REAL','IBEIG',KIBEIG ,NORBT,WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','ESRPRP')
      CALL REACMO_no_work(LUCOEF,'DFCOEF',DUMMY,DUMMY,WORK(KIBEIG),
     &                    TOTERG,8)
      CLOSE(LUCOEF,status='KEEP')
      CALL PRPMAT(INDXPR,IOPSY,PMAT,CMO,WORK(KIBEIG),IPCMO,NASH,
     &            WORK,KFREE,LFREE,IPRESR)
      CALL MEMREL('after PRPMAT',WORK,1,KIBEIG,KFREE,LFREE)
c    &            WORK,KFREE,LFREE,10)
!d    call memchk('esrTRANS 3',work,1)
      CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,IPRESR)
chjdbg    CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,40)
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

      IF (MZ .EQ. 2) THEN
         JZ = 1
      ELSE IF (REAL_PMLF) THEN
         JZ = 1
      ELSE
         JZ = 2
      END IF

      IF (MAX_CIESR_IT .GT. 0) THEN
C
C         Make property matrix in ESRVAL directly from the MULTIPESR first CI vectors
C           CIEVEC(NZCONF, N_CIESR, MZ)
C           ESRVAL(MULTIPESR, MULTIPESR, 2)
C

#ifdef ESR_DEBUG
         write(lupri,*)
     &   'Check hermiticity before GASCIP_MAKE_P1MAT ', IHRM
         write(lupri,*)
     &    'N_CIESR,NZCONF,KZCONF,MULTIPESR,IHRM,KPMLF',
     &    N_CIESR,NZCONF,KZCONF,MULTIPESR,IHRM,KPMLF
       WRITE(LUPRI,*)
     &   'GASCIP_MAKE_P1MAT: Real part of the N_CIESR CI vectors'//
     &   ' BEFORE MAKE_P1MAT'
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
         CALL GASCIP_MAKE_P1MAT(CIEVEC,N_CIESR,NZCONF,WORK(KZCONF),
     &                         ESRVAL(1,1,JZ),MULTIPESR,IHRM,
     &                         WORK(KPMLF),WORK(KFREE),LFREE,IPRESR)
C       CALL GASCIP_MAKE_P1MAT(CIVEC,N_CIVEC,NDET,IDET,
C    &                         P1MAT_CI,N_P1MAT_CI,IHRM_P1,
C    &                         P1MAT_MO,  WORK,LWORK, IPRINT)
      ELSE  ! here if MAX_CIESR_IT .le. 0
C
C         Our first non-iterative code with full, explicit CI matrices.
C         Construct complete property matrix over determinants
C         and transform to the desired property matrix in ESRVAL from
C         the MULTIPESR first CI vectors with QTRANS afterwards.
C
         LPCI = MZ*NZCONF*NZCONF
         CALL MEMGET2('REAL','PR',KPR,LPCI,WORK,KFREE,LFREE)
         KPI  = KPR + NZCONF*NZCONF
         CALL GASCIP_MAKEH(.FALSE.,IHRM,NZCONF,WORK(KZCONF),
     &         WORK(KPR),WORK(KPI),WORK(KPMLF),WORK(KMUUUU))
         IF (IPRESR .GT. 10) THEN
            IF (MZ .EQ. 2) THEN
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
            IF (MZ .GT. 1) THEN
               WRITE (LUPRI,'(/A)')
     &          ' Imag part of det-CI property matrix'
               CALL OUTPUT(WORK(KPI),1,NZCONF,1,NZCONF,
     &              NZCONF,NZCONF,1,LUPRI)
            END IF
         END IF
C
C        Transform to eigenvector basis
C

         CALL QTRANS('AOMO','S',D0,NZCONF,NZCONF,MULTIPESR,MULTIPESR,
     &        WORK(KPR),NZCONF,NZCONF,MZ,IQINDX(JZ),
     &        ESRVAL,   MULTIPESR,MULTIPESR, 2,IQINDX,
     &        CIEVEC,   NZCONF,N_CIESR,     MZ,IQINDX,
     &        CIEVEC,   NZCONF,N_CIESR,     MZ,IQINDX,
     &        WORK(KFREE),LFREE,-1)
       END IF
C
       IF (IPRESR .GT. 2) THEN
          CALL HEADER('ESRTRANS results in CI eigenvector basis'//
     &                ' (test output)',1)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,3X,A16)')
     &                'Operator no. 1 :',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(1,INDXPR)
          CALL PRQMAT(ESRVAL,MULTIPESR,MULTIPESR,
     &                MULTIPESR,MULTIPESR,NZCMPLX,IQINDX,LUPRI)
C
          CALL MEMREL('ESRTRANS',WORK,1,KFRSAV,KFREE,LFREE)
C
       END IF
       IF (ITIM.EQ.1) THEN
          WRITE (LUPRI,*)' WARNING: Inactive contribution not added!'
C TODO
c         DO I = 1,MULTIPESR
c            ESRVAL(I,I,1) = ESRVAL(I,I,1) + PCORE
c         END DO
       END IF
      CALL QEXIT('ESRTRANS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesrodd */
      SUBROUTINE PRPESRODD(ESRVAL,IREP,ITIM,PMLF)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dcbesr.h"
#include "dgroup.h"
#include "dcborb.h"
C
      DIMENSION ESRVAL(N_CIESR,N_CIESR,2),PMLF(NASHT,2,NASHT,2,*)
      LOGICAL   REAL_PMLF
C
      if (2*NASHT .ne. N_CIESR) then
         print *, 'PRPESRODD error'
         print *,'N_CIESR,NASHT',N_CIESR,NASHT
         stop 'PRPESRODD error'
      end if
      MZ = MIN(2,NZ)
      N = N_CIESR*N_CIESR*MZ
      IF (MZ .EQ. 2) THEN
         call dcopy(N,PMLF,1,ESRVAL,1)
      ELSE
C        ... real group, MZ = NZ = 1
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
      PARAMETER(D0 = 0.0D0)
C
#include "dcbesr.h"
#include "dgroup.h"
#include "dcborb.h"
C
      DIMENSION ESRVAL(N_CIESR,N_CIESR,2),PMLF(NASHT,2,NASHT,2,*)
      LOGICAL   REAL_PMLF
C
      CALL QUIT(
     &'ESR only implemented for odd number of active electrons')
      MZ = MIN(2,NZ)
      N = N_CIESR*N_CIESR*MZ
      IF (MZ .EQ. 2) THEN
C
      ELSE
C        ... real group, MZ = NZ = 1
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
      PARAMETER(D0 = 0.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbesr.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
C
      LOGICAL   FIRST(4),LESRST(2), RODHF
      DIMENSION DMAT(N2BBASX,NZ,*),PMAT(N2BBASX,NZ),
     &          PTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
C     Check if single open shell calculation (.OPEN 1 1/2,0 or .OPEN 1 1/0,2)
C
      RODHF = NOPEN.EQ.1 .AND. NAELEC.EQ.1 .AND. NASHT.EQ.1
C
C     Check symmetry
C
      PVAL = D0
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
CMI   CALL PRPMAO(LU1INT,INDXPR,PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
      CALL PRPMAO(LU1INT,INDXPR,NNBBASX,
     &             PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
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
C     /* Deck PRP_ESR */
      SUBROUTINE PRP_ESR(ESRVAL,WORK,KFREE,LFREE)
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
      DIMENSION ESRVAL(MULTIPESR,MULTIPESR,2,*),WORK(*)
      KFRSAV = KFREE
      NCMAT =  NUCDEP*9
      CALL MEMGET2('REAL','TMAT',KTMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TISO',KTISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBUF',KIBUF,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET2('LOGI','ATOM',KATOM,NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IATM',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','RATM',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','CSTRA',KCSTRA,
     &   NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SCTRA',KSCTRA,
     &   NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TMPMAT',KTMPMAT,NCMAT,WORK,KFREE,LFREE)
C
      CALL SETATM(WORK(KATOM),NATOM,29)
      CALL PRP_ESR1(ESRVAL,WORK(KTMPMAT),WORK(KTMAT),WORK(KTISO),
     &     WORK(KIATM),WORK(KRATM),WORK(KIBUF),
     &     WORK(KATOM),WORK(KCSTRA),WORK(KSCTRA) )
      CALL MEMREL('PRP_ESR',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ESR1 */
      SUBROUTINE PRP_ESR1(ESRVAL,TMPMAT,TMAT,TISO,IATINF,RATINF,
     &                    IBUFA,DOATOM,CSTRA,SCTRA)
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
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0)
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION ESRVAL(MULTIPESR,MULTIPESR,2,*),
     &          TMPMAT(3,*),TMAT(3,3,NUCDEP),
     &          TISO(NUCDEP),IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          IBUFA(3*NUCDEP),IBUFB(3),CSTRA(NUCDEP*NUCDEP*9),
     &          SCTRA(NUCDEP*NUCDEP*9)
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
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
        IZ   = NINT(CHARGE(INUC))
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
         jmhlf = 1
         jphlf = 2
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            TMPMAT(IPTAX(1,2),JA) =  D2*ESRVAL(jphlf,jmhlf,1,KA)
            TMPMAT(IPTAX(2,2),JA) = -D2*ESRVAL(jphlf,jmhlf,2,KA)
            TMPMAT(IPTAX(3,2),JA) =  D2*ESRVAL(jphlf,jphlf,1,KA)
         ENDDO
C
      ELSE IF (MULTIPESR .EQ. 3) THEN
C     .... triplet molecules (Seff = 1)
         jm1 = 1
         j0  = 2
         jp1 = 3
C
         SQRT2 = SQRT(D2)
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            TMPMAT(IPTAX(1,2),JA) = SQRT2*ESRVAL(j0 ,jm1,1,KA)
            TMPMAT(IPTAX(2,2),JA) = SQRT2*ESRVAL(j0 ,jm1,2,KA)
            TMPMAT(IPTAX(3,2),JA) =       ESRVAL(jp1,jp1,1,KA)
         ENDDO
C
      ELSE IF (MULTIPESR .EQ. 4) THEN
C     .... quartet molecules (Seff = 3/2)
C
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            TMPMAT(IPTAX(1,2),JA) =    ESRVAL(2,3,1,KA)
            TMPMAT(IPTAX(2,2),JA) =   -ESRVAL(2,3,2,KA)
            TMPMAT(IPTAX(3,2),JA) = D2*ESRVAL(3,3,1,KA)
         ENDDO
      ELSE
         WRITE(LUPRI,'(//A/A)') '@ WARNING:'//
     &   ' Only doublet and triplet hfcc implemented yet',
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
               IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(IATOM)) THEN
                     NSPEC = NSPEC + 1
                     CALL ESRANA(TMAT(1,1,IATOM), TISO(IATOM),
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
C     /* Deck ESRANA */
      SUBROUTINE ESRANA(TMAT,AVETOT,IATINF,RATINF,
     &                  NAME,NAMEX,IPRINT)
C***********************************************************************
C
C     Analyze routine for ESR properties
C     This routine is modelled after SHIANA in DALTON
C
C     April 2000 A. Noerager
C
C***********************************************************************
      use codata
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (FACTOR = ALPHAC*ALPHAC/(D2*XFMP))
C                factor = 1/c^2 * \mu_N = 1/c^2 * 1/(2 m_p)
      CHARACTER NAME*6, NAMEX(3)*6
      DIMENSION TMAT(3,3), TMATNY(3,3), IATINF(3), RATINF(3,*),
     &          TMATX(3,3)
#include "chrxyz.h"
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
         AVETOT = (TMATNY(1,1)+TMATNY(2,2)+TMATNY(3,3))/D3
         AVETOTX= (XTHZ*1.0D-6) * AVETOT

C
         CALL TITLER(
     &      'ESR hyperfine coupling constants for '//NAME//':','=',1)
C
         WRITE (LUPRI,'(A,I5//A,F20.2,A)')
     &   '  Isotope number:',NINT(RATINF(1,L)),
     &   '  Isotropic hyperfine coupling constant',AVETOTX, ' MHz'
         IF (IPRINT .GT. 0) WRITE (LUPRI,'(A,1P,D20.10,A)')
     &   '  Isotropic hyperfine coupling constant',AVETOT,  ' au'
C
         CALL HEADER('Total hyperfine coupling matrix (MHz):',1)
         WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
         DO J = 1, 3
            WRITE (LUPRI,'(2X,A6,3F15.2)')
     &             NAMEX(J), (TMATX(I,J),I=1,3)
         END DO
C
         IF (IPRINT .GT. 0) THEN
            CALL HEADER('Total hyperfine coupling matrix (au):',1)
            WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
            DO J = 1, 3
               WRITE (LUPRI,'(2X,A6,1P,3D20.10)')
     &                NAMEX(J), (TMATNY(I,J),I=1,3)
            END DO
         END IF
      ENDDO
C
C     END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE  ESR_READ_CMO(CMO,WORK,KFREE,LFREE)
C
C     Hans Joergen Aa. Jensen Feb. 2015
C
C     Read MO coefficents for ESR-CI calculation in this priority order:
C        1. MP2NO on KRMCOLD
C        2. DFCOEF with Hartree-Fock or DFT or MVO orbitals
C
#include "implicit.h"
#include "priunit.h"
      REAL*8  CMO(*), WORK(*)
      LOGICAL FILEX
C
C Used from COMMON blocks:
C  dcbgen.h : LUCOEF, LUKRMC
C  dcbham.h : SPINFR; LEVYLE
C  dcbbas.h : N2BBASXQ, ...
C  dcborb.h : NORBT, ...
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"

C
      KFREE_start = KFREE
      CALL DZERO(CMO,N2BBASXQ)

C        1. MP2NO on KRMCOLD ?

      INQUIRE(FILE='KRMCOLD',EXIST=FILEX)
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

         WRITE(LUPRI,'(/A)')
     &      ' (ESR)  MP2 natural orbitals read from file KRMCOLD'
C
         GO TO 9000
      END IF

C        MP2NO not found. Then:
C        2. DFCOEF with Hartree-Fock or DFT or MVO orbitals ?

      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(FILEX) THEN
         CALL MEMGET2('REAL','EIG',KEIG ,NORBT   ,WORK,KFREE,LFREE)
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPESR')
         IF (SPINFR.OR.LEVYLE) THEN
            CALL MEMGET2('INTE','IBPRP',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
            CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,WORK(KIBRP),
     &                  TOTERG,10,WORK,KFREE,LFREE)
         ELSE
            KIBRP = KFREE
            CALL REACMO(LUCOEF,'DFCOEF',CMO,WORK(KEIG),IDUM,
     &                  TOTERG,2,WORK,KFREE,LFREE)
         END IF
         CLOSE(LUCOEF,STATUS='KEEP')
         WRITE(LUPRI,'(/A)')
     &      ' (ESR)  molecular orbitals read from file DFCOEF'
         GO TO 9000
      END IF

C        DFCOEF not found. Then:
C        3. no more options (yet).

      WRITE(LUPRI,'(//A)') 'Fatal error for ESR: MO coefficients not'//
     &   ' found on KRMCOLD or DFCOEF'
      CALL QUIT('ESR fatal error: no MO coefficients found')

 9000 CALL MEMREL('ESR_READ_CMO',WORK,1,KFREE_start,KFREE,LFREE)
      RETURN
      END
!  -- end of pamesr.F --
