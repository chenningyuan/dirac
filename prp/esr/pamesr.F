!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "orgcom.h"
      PARAMETER (NTABLE = 8)
      LOGICAL SET, NEWDEF, RESET, GORG
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.GAUGEO','.PRINT ','.OPERAT','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
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
      IPRESR  = 0
      NESRP   = 0
      NCIESR  = -1
C
C     Initialize orgcom.h
C     ===================
C
      GORG = .FALSE.
      CALL DZERO(GAGORG,3)
C
C     Process input 
C     =========================
C
      ICHANG = 0
      INPERR = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
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
    1          CONTINUE
C Read user specified gauge origin
                  READ(LUCMD,*) (GAGORG(ICOOR), ICOOR = 1,3)
                  GORG = .TRUE.
               GO TO 100
    2          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRESR
               GO TO 100
    3          CONTINUE
C&&&& OPERATOR: Define operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPRESR)
                  CALL OP1IND('ESRINP',INDESR,LESRP,NESRP,INDXPR,
     &                        MAXESR)
                  DOESR = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ESRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ESRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF (GORG) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A)') 'ESRINP: '
         CALL PRSYMB(LUPRI,'=',75,0)
         IF (GORG) THEN
            WRITE(LUPRI,'(/,A,3F12.6,/)') 'User defined gauge '//
     &           'origin used:',(GAGORG(I), I = 1, 3)
         END IF
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
      END IF
  999 CONTINUE
      IF (INPERR.GT.0) CALL QUIT('Input error in *ESR')
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
      DIMENSION ESRVAL(NCIESR,NCIESR,NZCMPLX,*),WORK(*)
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
C
      CALL QENTER('PRPESR')
C
      KFRSAV = KFREE
      print *,'PRPESR  KFREE,LFREE',KFREE,LFREE
      call memchk('prpesr 1',work,1)
C
      IF (WF(1:3) .NE. 'DHF' .AND. WF(1:3) .NE. 'DFT') THEN
         WRITE (LUPRI,'(/2A/A)')
     &   ' WARNING: ESR properties are only'//
     &   ' implemented for open shell DHF; not for ',WF,
     &   ' WARNING: .ESR request is ignored !'
         GO TO 9999
      END IF
C
C     Retrieve molecular orbitals
C
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPESR')
      IF(SPINFR.OR.LEVYLE) THEN
        CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
        CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,WORK(KIBRP),
     &              TOTERG,10,WORK,KFREE,LFREE)
      ELSE
        KIBRP = KFREE
        CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,IDUM,TOTERG,2,
     &              WORK,KFREE,LFREE)
      ENDIF
CTROND      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),IDUM,TOTERG,6,
CTROND     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Solve ESR CI problem
C     on exit: WORK(KEIGVL) - NZCONF CI eigenvalues
C              WORK(KEIGVC) - NZCONF CI eigenvectors (NZCONF,NZCONF,MZ)
C
      CALL ESR_RESOLVE(WORK(KCMO),IPRESR,WORK,KEIGVL,KEIGVC,KFREE,LFREE)
C
C     Memory allocation for PRPESR1
C
      CALL MEMGET('LOGI',KESRST ,2*NESRP   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT  ,NZ*N2ASHX ,WORK,KFREE,LFREE)
C
      call memchk('prpesr 2',work,1)
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
      CALL PRPESR1(ESRVAL,WORK(KESRST),WORK(KCMO),WORK(KIBRP),
     &             WORK(KPMAT),WORK(KEIGVL),WORK(KEIGVC),
     &             WORK,KFREE,LFREE)
      call memchk('prpesr 3',work,1)
C
      IF (IPRESR .GT. 2) THEN
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
          CALL PRQMAT(ESRVAL(1,1,1,I),NCIESR,NCIESR,
     &                NCIESR,NCIESR,NZCMPLX,IQINDX,LUPRI)
        ENDDO
      END IF
C
      IF (ESRGTENS) THEN
C
      IF (NAELEC .EQ. 1) THEN
C     ....doublet molecules
C      
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
      DO K = 1,3
C        ... for B field in direction k (x,y,z)
         I = IPESR(MXCOOR+K,1)
         GTENSOR(K,1) =  D4*ESRVAL(1,2,1,I)
         GTENSOR(K,2) = -D4*ESRVAL(1,2,2,I)
         GTENSOR(K,3) =  D4*ESRVAL(1,1,1,I)
         J = IPESR(MXCOOR+K,2)
         GCORR(K,1)   =  D4*ESRVAL(1,2,1,J)
         GCORR(K,2)   = -D4*ESRVAL(1,2,2,J)
         GCORR(K,3)   =  D4*ESRVAL(1,1,1,J)
      ENDDO
C
      ELSE IF (NAELEC .EQ. 2) THEN
C     .... triplet molecules
C        Find Sigma_Z expectation values for the three states:
         I = IPSIGMA(3)
         WRITE (LUPRI,*) 'Analysis to find triplet components'
         JMS(-1) = -1
         JMS( 0) = -1
         JMS( 1) = -1
         DO K = 1,3
            SZK  = ESRVAL(K,K,1,I)
            ISZK = NINT(SZK)
            WRITE (LUPRI,*)'State, <Sigma_Z>, S_Z(eff)',K,SZK,ISZK
            IF (ISZK .GT. 1 .OR. ISZK .LT. -1) THEN
               CALL QUIT('error: S_Z(eff) too big to triplet')
            ELSE IF (JMS(ISZK) .NE. -1) THEN
               CALL QUIT('error: same S_Z(eff) twice!!!')
            END IF
            JMS(ISZK) = K
         END DO
         jm1 = jMS(-1)
         j0  = jMS( 0)
         jp1 = jMS( 1)
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C     (2) Get the QED correction in GCORR
C
      D2SQRT2 = D2*SQRT(D2)
      DO K = 1,3
C        ... for B field in direction k (x,y,z)
         I = IPESR(MXCOOR+K,1)
         GTENSOR(K,1) = -D2SQRT2*ESRVAL(jp1,j0 ,1,I)
         GTENSOR(K,2) =  D2SQRT2*ESRVAL(jm1,j0 ,2,I)
         GTENSOR(K,3) =  D2*ESRVAL(jp1,jp1,1,I)
c??      GTENSOR(K,3) =  ESRVAL(jp1,jp1,1,I) - ESRVAL(jm1,jm1,1,I)
         J = IPESR(MXCOOR+K,2)
         GCORR(K,1)   =  -D2SQRT2*ESRVAL(jp1,j0,1,J)
         GCORR(K,2)   =   D2SQRT2*ESRVAL(jm1,j0,2,J)
         GCORR(K,3)   =   D2*ESRVAL(jp1,jp1,1,J)
      ENDDO
C
C     ***************************************
C     Zero-field-splitting parameters D and E
C     ***************************************
C
      ZFS_D = -(WORK(KEIGVL) - (WORK(KEIGVL+1)+WORK(KEIGVL+2))/D2)
C..   ZFS_D = EIGVAL(1) - (EIGVAL(2)+EIGVAL(3))/D2
      ZFS_E = (WORK(KEIGVL+1)-WORK(KEIGVL+2))/D2
      ZFS_DCM = ZFS_D*XTKAYS
      ZFS_ECM = ZFS_E*XTKAYS
      CALL HEADER('Zero field splitting parameters:',1)
         WRITE (LUPRI,'(5(/,1X,A,F15.8,A))')
     &   'ZFS parameter D:                ',ZFS_DCM,  '(cm-1) '
         WRITE (LUPRI,'(5(/,1X,A,F15.8,A))')
     &   'ZFS parameter E:                ',ZFS_ECM,  '(cm-1) '
c     print *, 'ZERO-FIELD D= ',ZFS_DCM, '(cm-1)'
      ELSE
         call quit('Only doublet and triplet implemented yet')
      END IF
C
      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The g-tensor for H(Zeeman) = 1/2c B g Seff :',
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


      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The G-matrix = g * g^T :', ((GMATRIX(I,K),K=1,3),I=1,3)

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC,WORK(KFREE),LFREE)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - 2.0D0)*1.D6
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,'(/A//A//5X,4F18.8//A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The Dirac equation g values corresponding to g_e = 2',
     &   ' The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   ' The g shifts in ppm (x, y, z, iso):',
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
     &   ' The corrected g values (x, y, z, iso):',
     &   (GCOR(I),I =1,3), GISO,
     &   ' The corrected g shifts in ppm (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO
C
C     (3) Get the corrected g values by adding correction term
C
      DO I = 1,3
      DO K = 1,3
C        ... for B field in direction k (x,y,z)
         GTENSOR(K,I) =  GTENSOR(K,I) + GCORR(K,I)
      ENDDO
      ENDDO

      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The g_1 corrected g-tensor for H(Zeeman) = 1/2c B g Seff :',
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


      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The G-matrix = g * g^T :', ((GMATRIX(I,K),K=1,3),I=1,3)

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC,WORK(KFREE),LFREE)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - GFAC)*1.D6
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,
     &   '(/A,F13.10//A//5X,4F18.8//A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The g_1 corrected Dirac eq. g values corresponding to g_e ='
     &   ,GFAC,
     &   ' The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   ' The g shifts in ppm (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO,
     &   ' The corresponding eigenvectors:',
     &   ((GEVEC(I,K),K=1,3),I=1,3)

      call memchk('prpesr 4',work,1)
      ENDIF
C     .. end if (ESRGTENS) then
C
C     Print ESR anisotropies from hyperfine Hamiltonian
C
      IF (ESRANISO)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesr1 */
      SUBROUTINE PRPESR1(ESRVAL,LESRST,CMO,IBEIG,PMAT,
     &                   CIENER,CIEVEC,WORK,KFREE,LFREE)
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
      DIMENSION ESRVAL(NCIESR,NCIESR,2,*),CMO(*),IBEIG(*),PMAT(*),
     &          CIENER(NCIESR), CIEVEC(NCIESR,NCIESR,*), WORK(*)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   GETIOPEN, IPCMO(2), IPCMPL(2)
      DATA IPCMPL/1,2/
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
      IF (NCIESR .NE. NZCONF) CALL QUIT('NCIESR .ne. NZCONF')
      N2CIESR = NCIESR*NCIESR
      CALL DZERO(ESRVAL,N2CIESR*2*NESRP)
      N2ESR = 2*NESRP
      CALL LSET(N2ESR,.FALSE.,LESRST)
C
C     Generate density matrix
C     =======================
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPESR')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
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
        MZ = MIN(2,NZ)
        CALL MEMGET('REAL',KPMLF  ,MZ*(2*NASHT)**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KMUUUU ,0              ,WORK,KFREE,LFREE)
        LPCI = N2CIESR
        IF (MZ .EQ. 2) LPCI = 2*LPCI
        CALL MEMGET('REAL',KPR    ,LPCI           ,WORK,KFREE,LFREE)
        KPI = KPR + N2CIESR
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
chj       IF (ITIM .GE. 0) CALL QUIT(
chj  &    'PRPESR1 only implemented for time-antisymmetric operators')
C

          IPCMO(1) = icmoq(1) + nfbas(1,0)*(npsh(1)+nish(1))
          IPCMO(2) = icmoq(2) + nfbas(2,0)*(npsh(2)+nish(2))
Chj-start test output
c         iprtst = 10
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
          CALL PRPMAT(INDXPR,IOPSY,PMAT,CMO,IBEIG,IPCMO,NASH,
     &                WORK,KFREE,LFREE,IPRESR)
c    &                WORK,KFREE,LFREE,10)
      call memchk('prpesr1 3',work,1)
          CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,40)
chjdbg    CALL QFC2MFC(PMAT,WORK(KPMLF),ISYM,ITIM,IPRESR)
C
C         Calculate CI property matrix over determinants
C
          IF (NZ.GT.1) THEN
             IHRM = 1
             REAL_PMLF = .TRUE.
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
          CALL GASCIP_MAKEH(.FALSE.,IHRM,WORK(KZCONF),
     &         WORK(KPR),WORK(KPI),WORK(KPMLF),WORK(KMUUUU))
          IF (ITIM.EQ.1) THEN
             WRITE (LUPRI,*)' WARNING: Inactive contribution not added!'
C TODO
c            DO I = 1,NZCONF
c               JPR = KPR + (I-1)*NCIESR + (I-1)
c               WORK(JPR) = WORK(JPR) + PCORE
c            END DO
          END IF
c         IF (IPRESR .GT. 5) THEN
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
             CALL OUTPUT(WORK(KPR),1,NCIESR,1,NCIESR,
     &           NCIESR,NCIESR,1,LUPRI)
             IF (MZ .EQ. 2) THEN
                WRITE (LUPRI,'(/A)')
     &          ' Imag part of det-CI property matrix'
                CALL OUTPUT(WORK(KPI),1,NCIESR,1,NCIESR,
     &              NCIESR,NCIESR,1,LUPRI)
             END IF
c         END IF
C
C         Transform to eigenvector basis
C
          IF (MZ .EQ. 2) THEN
             JZ = 1
          ELSE IF (REAL_PMLF) THEN
             JZ = 1
Cunoedvendig CALL DZERO(ESRVAL(1,1,2,I),N2CIESR)
          ELSE
             JZ = 2
Cunoedvendig CALL DZERO(ESRVAL(1,1,1,I),N2CIESR)
          END IF
      
          CALL QTRANS('AOMO','S',D0,NCIESR,NCIESR,NCIESR,NCIESR,
     &        WORK(KPR),NCIESR,NCIESR,MZ,IPCMPL(JZ),
     &        ESRVAL(1,1,1,I),NCIESR,NCIESR,2,IPCMPL,
     &        CIEVEC,NCIESR,NCIESR,MZ,IPCMPL,
     &        CIEVEC,NCIESR,NCIESR,MZ,IPCMPL,
     &        WORK(KFREE),LFREE,-1)
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
      DIMENSION ESRVAL(NCIESR,NCIESR,2),PMLF(NASHT,2,NASHT,2,*)
      LOGICAL   REAL_PMLF
C
      if (2*NASHT .ne. NCIESR) then
         print *, 'PRPESRODD error'
         print *,'NCIESR,NASHT',NCIESR,NASHT
         stop 'PRPESRODD error'
      end if
      MZ = MIN(2,NZ)
      N = NCIESR*NCIESR*MZ
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
#include <dcbesr.h>
#include <dgroup.h>
#include <dcborb.h>
C
      DIMENSION ESRVAL(NCIESR,NCIESR,2),PMLF(NASHT,2,NASHT,2,*)
      LOGICAL   REAL_PMLF
C
      CALL QUIT(
     &'ESR only implemented for odd number of active electrons')
      MZ = MIN(2,NZ)
      N = NCIESR*NCIESR*MZ
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
            write (6,*) 'iz, df(1),esrval',iz,df(1),esrval
            write (6,*) 'dmat(1,iz,2)'
            nbast = ntbas(0)
            call output(dmat(1,iz,2),1,nbast,1,nbast,nbast,nbast,1,6)
            write (6,*) 'pmat(1,iz)'
            call output(pmat(1,iz),1,nbast,1,nbast,nbast,nbast,1,6)
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
      DIMENSION ESRVAL(NCIESR,NCIESR,2,*),WORK(*)
      KFRSAV = KFREE
      NCMAT =  NUCDEP*9
      CALL MEMGET('REAL',KTMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KATOM,NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMPMAT,NCMAT,WORK,KFREE,LFREE)
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
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION ESRVAL(NCIESR,NCIESR,2,*), TMPMAT(3,*),TMAT(3,3,NUCDEP),
     &          TISO(NUCDEP),IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          IBUFA(3*NUCDEP),IBUFB(3),CSTRA(NUCDEP*NUCDEP*9),
     &          SCTRA(NUCDEP*NUCDEP*9)
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C     
      CALL TITLER('DIRAC  - ESR PROPERTIES','*',116)
C
C
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &'@ name  charge  isotope   mass         abundance    g_N factor'
      CALL PRSYMB(LUPRI,'-',61,0)
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
        RETURN
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
C     Fill ESR hfs anisotropy matrices Seff_p * TMPMAT * I_N,q
C
C     
         CALL PRSYMB(LUPRI,'-',72,2)

         CALL DZERO(TMPMAT,NCMAT)
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            TMPMAT(IPTAX(1,2),JA) =  ESRVAL(1,2,1,KA)
            TMPMAT(IPTAX(2,2),JA) = -ESRVAL(1,2,2,KA)
            TMPMAT(IPTAX(3,2),JA) =  ESRVAL(1,1,1,KA)
         ENDDO
         IF (IPRPRP.GE.4) THEN
            CALL HEADER('ESR hfs anisotropies'//
     &           ' in symmetry coordinates (B~I,au)',1)
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
      PARAMETER (FACTOR = ALPHAC*ALPHAC/(XFMP))
C                factor = 1/c^2 * \mu_N = 1/c^2 * 1/(2 m_p)
      CHARACTER NAME*6, NAMEX(3)*6
      DIMENSION TMAT(3,3), TMATNY(3,3), IATINF(3), RATINF(3,*),
     &          TMATX(3,3)
#include "chrxyz.h"
C #include <dcbxlr.h>
C
      CALL DSCAL(9,FACTOR,TMAT,1)
C
      NISO=IATINF(3)
  
      DO L=1,NISO	
         DO I=1,3
           DO J=1,3
              TMATNY(I,J)=RATINF(3,L)*TMAT(I,J)
              TMATX(I,J)=XTHZ*TMATNY(I,J) * 1.0D-6
           ENDDO
         ENDDO
      AVETOT = (TMATNY(1,1)+TMATNY(2,2)+TMATNY(3,3))/D3
      AVETOTX= XTHZ*AVETOT * 1.0D-6
       
       
C
C     IF (IPRINT.GE.1) THEN
      CALL TITLER('ESR anisotropy for '//NAME//':','=',1)
C
         WRITE (LUPRI,'(5(/,1X,A,1P,D20.10,A))')
     &   ' Isotropic anisotropy:   ',AVETOT,  ' au'
         WRITE (LUPRI,'(5(/,1X,A,F15.2,A))')
     &   ' Isotropic anisotropy:        ',AVETOTX,  ' MHz'
         WRITE (LUPRI,'(5(/,1X,A,F15.8,A))')
     &   ' Isotope number:                ',RATINF(1,L),  ' ' 
      CALL HEADER('Total anisotropy tensor (au):',1)
      WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
      DO 400 J = 1, 3
         WRITE (LUPRI,'(2X,A6,1P,3D20.10)')
     &          NAMEX(J), (TMATNY(I,J),I=1,3)
  400 CONTINUE
      CALL HEADER('Total anisotropy tensor (MHz):',1)
      WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
      DO 500 J = 1, 3
         WRITE (LUPRI,'(2X,A6,3F15.2)')
     &          NAMEX(J), (TMATX(I,J),I=1,3)
  500 CONTINUE
      ENDDO
C
C     END IF
      RETURN
      END
