      SUBROUTINE QUIT(TEXT)
      CHARACTER TEXT*(*)
      stop
      end subroutine

      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
      LHEAD = LNBLNK(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(80 - LHEAD)/2 + 1)
      END IF
      WRITE (6, '(//,150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
      WRITE (6, '(   150A)') (' ',I=1,INDENT), ('-',I=1,LHEAD)
      WRITE (6, '()')
      RETURN
      END

      SUBROUTINE MXMA(A,ICA,IRA,B,ICB,IRB,C,ICC,IRC,NROW,NSUM,NCOL)
#include "implicit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION A(*),B(*),C(*)
#if !defined (VAR_ESSL) 
      DO 100 IROW=0,NROW-1
      DO 100 ICOL=0,NCOL-1
         SUM=D0
         DO 110 ISUM=0,NSUM-1
            SUM=SUM +
     &          A(1+IROW*ICA+ISUM*IRA) *
     &          B(1+ISUM*ICB+ICOL*IRB)
110      CONTINUE
         INDEX=1+ICC*IROW+IRC*ICOL
         C(INDEX)=SUM
100   CONTINUE
#else
      IF(NSUM.EQ.0) GOTO 600
      IF((NROW.EQ.0).OR.(NCOL.EQ.0)) RETURN
C  TRY TO FIT TO A STANDARD CASE AND SET SELECTOR:
      ISEL=1
      IF ((ICA.EQ.1).AND.(IRA.GE.NROW)) GOTO 101
      IF ((IRA.EQ.1).AND.(ICA.GE.NSUM)) GOTO 102
      GOTO 400
 101  ISEL=ISEL+4
 102  IF ((ICB.EQ.1).AND.(IRB.GE.NSUM)) GOTO 201
      IF ((IRB.EQ.1).AND.(ICB.GE.NCOL)) GOTO 202
      GOTO 400
 201  ISEL=ISEL+2
 202  IF ((ICC.EQ.1).AND.(IRC.GE.NROW)) GOTO 301
      IF ((IRC.EQ.1).AND.(ICC.GE.NCOL)) GOTO 302
      GOTO 400
 301  ISEL=ISEL+1
C  NOTE -- CASES 2 OR 7 IN COMPUTED GOTO STATEMENT BELOW ARE OF
C  THE T-T TYPE, WHICH SOMETIMES GOES WRONG IN DGEMUL.
C  FOR THE TIME BEING, THIS MUST BE TREATED AS A GENERAL CASE.
C  CHANGE BACK AS SOON AS ESSL HAS BEEN PATCHED!! P-AA M 87-03-23.
C302  GOTO (1,2,3,4,5,6,7,8) ISEL
 302  GOTO (1,400,3,4,5,6,400,8) ISEL
      WRITE(*,*)' SOMETHING IS VERY WRONG -- THIS SPOT SHOULD NOT BE'
      WRITE(*,*)' REACHABLE!! MXMA STOPS.'

   1  CALL DGEMUL(B,ICB,'N',A,ICA,'N',C,ICC,NCOL,NSUM,NROW)
      RETURN

C  2  CALL DGEMUL(A,ICA,'T',B,ICB,'T',C,IRC,NROW,NSUM,NCOL)
C     RETURN

   3  CALL DGEMUL(B,IRB,'T',A,ICA,'N',C,ICC,NCOL,NSUM,NROW)
      RETURN

   4  CALL DGEMUL(A,ICA,'T',B,IRB,'N',C,IRC,NROW,NSUM,NCOL)
      RETURN

   5  CALL DGEMUL(B,ICB,'N',A,IRA,'T',C,ICC,NCOL,NSUM,NROW)
      RETURN

   6  CALL DGEMUL(A,IRA,'N',B,ICB,'T',C,IRC,NROW,NSUM,NCOL)
      RETURN

C  7  CALL DGEMUL(B,IRB,'T',A,IRA,'T',C,ICC,NCOL,NSUM,NROW)
C     RETURN

   8  CALL DGEMUL(A,IRA,'N',B,IRB,'N',C,IRC,NROW,NSUM,NCOL)
      RETURN

 400  DO 520 IROW=0,NROW-1
        DO 520 ICOL=0,NCOL-1
          T=D0
          DO 510 ISUM=0,NSUM-1
            T=T+A(1+IROW*ICA+ISUM*IRA)*B(1+ISUM*ICB+ICOL*IRB)
 510      CONTINUE
          C(1+IROW*ICC+ICOL*IRC)=T
 520  CONTINUE
      RETURN

 600  DO 620 IROW=0,NROW-1
        DO 620 ICOL=0,NCOL-1
          C(1+IROW*ICC+ICOL*IRC)=0.0D00
 620  CONTINUE
#endif
      RETURN
      END

      SUBROUTINE TIMTXT(TEXT,TIMUSD,LUPRIN)
C
C TIMTXT based on TIMER by TUH //900709-hjaaj
C
#include "implicit.h"
      CHARACTER*(*) TEXT
      CHARACTER AHOUR*6, ASEC*8, AMIN*8
C
      ISECND = NINT(TIMUSD)
      IF (ISECND .GE. 60) THEN
         MINUTE = ISECND/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = ISECND - 3600*IHOURS - 60*MINUTE
         IF (IHOURS .EQ. 1) THEN
            AHOUR = ' hour '
         ELSE
            AHOUR = ' hours'
         END IF
         IF (MINUTE .EQ. 1) THEN
            AMIN = ' minute '
         ELSE
            AMIN = ' minutes'
         END IF
         IF (ISECND .EQ. 1) THEN
            ASEC = ' second '
         ELSE
            ASEC = ' seconds'
         END IF
         IF (IHOURS .GT. 0) THEN
            WRITE(LUPRIN,100)
     *            TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
         ELSE
            WRITE(LUPRIN,200) TEXT, MINUTE, AMIN, ISECND, ASEC
         END IF
      ELSE
         WRITE(LUPRIN,300) TEXT,TIMUSD
      END IF
  100 FORMAT(1X,A,I4,A,I3,A,I3,A)
  200 FORMAT(1X,A,     I3,A,I3,A)
  300 FORMAT(1X,A,F7.2,' seconds')
      RETURN
      END

      SUBROUTINE TITLER(HEAD,A,IN)
#include "priunit.h"
      CHARACTER HEAD*(*), A*(*)
C
      LHEAD  = LNBLNK(HEAD)
      LENGTH = LHEAD
      IF (IN .EQ. 200) THEN
         LENGTH = LENGTH + 2
      ELSE IF (IN .GE. 100) THEN
         MARG = (80 - LENGTH)/2 - 1
         MARG = MIN(MARG, IN - 100)
         IF (MARG .GT. 0) MARG = MARG + 1
         LENGTH = LENGTH + 2*MARG
      END IF
      IF (IN .GE. 0 .AND. IN .LT. 100) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(80 - LENGTH)/2 + 1)
      END IF
      IF (IN .EQ. 200) THEN
         WRITE (LUPRI, '(//150A)')
     *      (' ',I=1,INDENT),'.', ('-',I=1,LENGTH),'.'
         WRITE (LUPRI, '(150A)')
     &      (' ',I=1,INDENT),'| ', HEAD(1:LHEAD), ' |'
         WRITE (LUPRI, '(150A)')
     *      (' ',I=1,INDENT),'`', ('-',I=1,LENGTH),"'"
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD(1:LHEAD), ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//150A)') (' ',I=1,INDENT), HEAD(1:LHEAD)
         WRITE (LUPRI, '(150A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
!     CALL FLSHFO(LUPRI)
      RETURN
      END

C  /* Deck spdcar */
      FUNCTION SPDCAR(IORDER)
#include "implicit.h"
      CHARACTER*1 SPDCAR
      IF (IORDER .EQ. 0) THEN
         SPDCAR = 's'
      ELSE IF (IORDER .EQ. 1) THEN
         SPDCAR = 'p'
      ELSE IF (IORDER .EQ. 2) THEN
         SPDCAR = 'd'
      ELSE IF (IORDER .LE. 6) THEN
         SPDCAR = CHAR(ICHAR('f') + IORDER - 3)
      ELSE
C
C        Note: j is not used - therefore this special case
C
         SPDCAR = CHAR(ICHAR('k') + IORDER - 7)
      END IF
      RETURN
      END

      FUNCTION FACULT(N)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1=1.D0)
#include "sdpre.h"
      IF (N .LT. 0) THEN
         WRITE (LUPRI,'(/,A,I10,/A)')
     &         ' Argument less than zero in FACULT:',N,
     &         ' Program cannot continue.'
         CALL QUIT('Illegal argument in FACULT')
      ELSE
        FACULT = D1
        DO 100 I = 1, N
           FACULT = FACULT*SDPRE(I)
  100   CONTINUE
      END IF
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck q2bphase */
      SUBROUTINE Q2BPHASE(TYP,IM,JFAC,AMAT)
C***********************************************************************
C
C     Quaternion phase insertion in square matrix
C       TYP = F:  (e_K1)(e_IM)(e_K2*)   as Fock matrix
C       TYP = D:  (e_K1*)(e_IM)(e_K2)   as density matrix
C
C     Written by T.Saue, January 1995
C     Last revision: June 25 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "dcbbas.h"
#include "dcbqal.h"
#include "dgroup.h"
      CHARACTER TYP*1
      DIMENSION AMAT(*)
      EXTERNAL BDQALGEBRA
C
      IF(TYP.NE.'F'.AND.TYP.NE.'D') THEN
        CALL QUIT('Q2BPHASE: Unknown TYP !')
      ENDIF
      DO 10 IRP1 = 0,NBSYM-1
        ISY1 = IRP1 + 1
        DO 20 IC1 = 1,2
          NB1 = NBBAS(IRP1,IC1)
          IF(NB1.EQ.0) GOTO 20
          K1  = JQBAS(IRP1,IC1)
          DO 30 IRP2 = 0,NBSYM-1
            ISY2 = IRP2 + 1
            DO 40 IC2 = 1,2
              NB2 = NBBAS(IRP2,IC2)
              IF(NB2.EQ.0) GOTO 40
              K2   = JQBAS(IRP2,IC2)
              IF    (TYP.EQ.'F') THEN
                KFAC = JFAC*IQPHASE(K1,K2,IM)*IQPH(K1,1)*IQPH(K2,2)
              ELSEIF(TYP.EQ.'D') THEN
                IQ   = IQMULT(K1,K2,IM)
                KFAC = JFAC*IQPHASE(K1,K2,IQ)*IQPH(K1,2)*IQPH(K2,1)
              ENDIF
              IF(KFAC.EQ.-1) THEN
CTROND                IOFF = I2COSX(ISY1,ISY2,IC1,IC2)
                IOFF = I2BBASX(IRP1,IRP2,IC1,IC2)
                DO 50 ICOL = 1,NB2
                  CALL DSCAL(NB1,DM1,AMAT(IOFF+1),1)
                  IOFF = IOFF + NTBAS(0)
   50           CONTINUE
              ENDIF
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      END

      SUBROUTINE WTAOSO(AAO,ASO,NBAST,NVEC,IPRINT)
C***********************************************************************
C
C     Vector Transform:  AO-basis to SO-basis
C     VTAOSO/WTAOSO are covariant/contarvariant versions
C       VTAOSO is good for vectors and WTAOSO for functions
C     This routine is based on HERMIT routine DAOTSO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sept 10 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NVEC), AAO(NBAST,NVEC)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine WTAOSO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(ASO,NBAST*NVEC)
      DO 100 IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO 200 I = ISTRA,ISTRA + NORBA - 1
          IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
          NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
          NHKTA  = NHKT(IA)
          KHKTA  = KHKT(IA)
          MULA   = ISTBAO(IA)
          INDA   = KSTRT(IA) + NA - KHKTA
C
C         Loop over symmetry dependent centers
C
          DO 300 ISYMA = 0, MAXOPR
          IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
            INDA = INDA + KHKTA
            FAC  = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
            CALL DAXPY(NVEC,FAC,AAO(INDA,1),NBAST,ASO(I,1),NBAST)
          ENDIF
  300     CONTINUE
  200   CONTINUE
        ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Vectors in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
         CALL HEADER('OUTPUT: Vectors in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NVEC,NBAST,NVEC,-1,LUPRI)
      END IF
      RETURN
      END

      SUBROUTINE OUTPUT (AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,
     *                   NCTL,LUPRI)
C.......................................................................
C Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
#include <implicit.h>
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
C
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
C
      AMAX = ZERO
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
   10 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 RETURN
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
      SUBROUTINE VTBUBS(VBU,VBS,NVEC)
C***********************************************************************
C
C     Vector transform:  from unsorted basis to sorted basis
C
C     Written by T.Saue Mar 2 1995
C     Last revision : Jan 7 1998 - jth LINUX
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBU(NBAST,NVEC),VBS(NBAST,NVEC)
#else
      DIMENSION VBU(NTBAS(0),NVEC),VBS(NTBAS(0),NVEC)
#endif
C
      DO J = 1,NVEC
        DO IC = 1,2
          DO IREP = 0,MAXREP
            IUOFF = ICOS (IREP+1,IC)
            ISOFF = IBBAS(IREP  ,IC)
            IF(NBBAS(IREP,IC).GT.0)
     &        CALL DCOPY(NBBAS(IREP,IC),VBU(IUOFF+1,J),1,
     &                                  VBS(ISOFF+1,J),1)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
C
      END

         SUBROUTINE DZERO(DX,LENGTH)
#include "implicit.h"
C
C Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
C
C   Subroutine DZERO sets a real array of length *LENGTH*
C   to zero.
C...................................................................
C     DIMENSION DX(*)
      DIMENSION DX(LENGTH)
C
      IF (LENGTH.LE.0) RETURN
C
      DO 100 I = 1,LENGTH
  100    DX(I) = 0.0D00
C
      RETURN
      END

      SUBROUTINE BSTOBU_no_work(AMAT,MZ)
C***********************************************************************
C
C     Matrix transform:
C     from class/symmetry sorted basis ("sorted basis" = Dirac sorted)
C     to symmetry/class sorted basis ("unsorted basis" = Hermit sorted)
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(*)
      real(8), allocatable :: buf(:)
C
C     Memory allocation
C
      allocate(buf(mz*n2bbasx))
      CALL DCOPY(MZ*N2BBASX,AMAT,1,BUF,1)
      IZOFF = 0
      DO IZ = 1,MZ
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = IZOFF + (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = IZOFF + (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF(NBBAS(IREP,IC).GT.0)
     &              CALL DCOPY(NBBAS(IREP,IC),buf(1+ISOFF),1,
     &                                        AMAT(IUOFF+1)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IZOFF = IZOFF + N2BBASX
      ENDDO
C
C     Memory deallocation
C
      deallocate(buf)
      END
      SUBROUTINE BUTOBS_no_work(AMAT,MZ)
C***********************************************************************
C
C     Matrix transform:
C     from symmetry/class sorted basis ("unsorted basis" = Hermit sorted)
C     to class/symmetry sorted basis ("sorted basis" = Dirac sorted)
C
C     Written by T.Saue Mar 2 1995
C     Last revision : May 6 1996 - tsaue
C
C***********************************************************************

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
C
      DIMENSION AMAT(*)
      real(8), allocatable :: buf(:)
C

      allocate(buf(mz*n2bbasx))
      CALL DCOPY(MZ*N2BBASX,AMAT,1,buf,1)
      IZOFF = 0
      DO IZ = 1,MZ
        DO JC = 1,2
          DO JREP = 0,MAXREP
            DO J = 0,NBBAS(JREP,JC)-1
              JSOFF = IZOFF + (IBBAS(JREP  ,JC)+J)*NTBAS(0)
              JUOFF = IZOFF + (ICOS (JREP+1,JC)+J)*NTBAS(0)
              DO IC = 1,2
                DO IREP = 0,MAXREP
                  IUOFF = JUOFF + ICOS (IREP+1,IC)
                  ISOFF = JSOFF + IBBAS(IREP  ,IC)
                  IF(NBBAS(IREP,IC).GT.0)
     &              CALL DCOPY(NBBAS(IREP,IC),buf(1+IUOFF),1,
     &                                        AMAT(ISOFF+1)   ,1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IZOFF = IZOFF + N2BBASX
      ENDDO
C
      deallocate(buf)
      RETURN
C
      END
C  /* Deck lmnval */
      SUBROUTINE LMNVAL(NHKTA,KHKTA,LVALUE,MVALUE,NVALUE)
C
C     tuh March 87
C
#include "implicit.h"
#include "maxmom.h"
#include "xyzpow.h"
      DIMENSION LVALUE(KHKTA), MVALUE(KHKTA), NVALUE(KHKTA)
#include "ibtfun.h"
C
      ICOMP = 0
      DO 100 I = 1, KHKTA
         ICOMP = ICOMP + 1
         LVALUE(ICOMP) = NHKTA - ISTEP(I)
         MVALUE(ICOMP) = MVAL(I)
         NVALUE(ICOMP) = NVAL(I)
  100 CONTINUE
      RETURN
      END

      FUNCTION POLVAL(NORDER,B,XVAL)
#include "implicit.h"
      DIMENSION B(NORDER+1)
      POLVAL = B(1)
      XBUF = 1
      DO 10 I = 2,(NORDER+1)
        XBUF = XBUF*XVAL
        POLVAL = POLVAL + XBUF*B(I)
   10 CONTINUE
      END

      SUBROUTINE DTSOAO(ASO,AAO,NBAST,IREPDM,IPRINT)
C***********************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "pgroup.h"
#include "ibtfun.h"
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DTSOAO',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         CALL HEADER('OUTPUT:Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
      RETURN
      END

      SUBROUTINE GATMAT(IOPT,LDA,N,A,B,IND,NIND)
C***********************************************************************
C
C     This routine will collect columns elements of square matrix A
C     in groups based on pointer index IND.
C     IOPT discerns two options:
C       IOPT = 0        pick out the elements of largest absolute value
C       IOPT = 1        sum up elements
C
C       Written by T.Saue, June 1994
C       Last revision: October 9 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(LDA,N),B(NIND,NIND),IND(N)
C
C
CTROND No initialization !!!!      CALL DZERO(B,NIND*NIND)
C
C     Pick out largest elements for each group
C     ========================================
C
      IF    (IOPT.EQ.0) THEN
        DO 10 J = 1,N
          JJ = IND(J)
          DO 20 I = 1,N
            II = IND(I)
            IF(ABS(A(I,J)).GT.ABS(B(II,JJ))) B(II,JJ) = ABS(A(I,J))
   20     CONTINUE
   10   CONTINUE
C
C     Sum up elements within each group
C     =================================
C
      ELSEIF(IOPT.EQ.1) THEN
        DO 30 J = 1,N
          JJ = IND(J)
          DO 40 I = 1,N
            II = IND(I)
            B(II,JJ) = B(II,JJ) + A(I,J)
   40     CONTINUE
   30   CONTINUE
      ELSE
        CALL QUIT('GATMAT:Unknown option IOPT!')
      ENDIF
C
      RETURN
C
      END

      SUBROUTINE REACMO_no_work(IUNIT,CMOFIL,CMO,EIG,IBEIG,TOTERG,IOPT)
C***********************************************************************
C
C     Read DHF-coefficients and eigenvalues from unformatted file
C     Read option is provided from bit-packed IOPT:
C       0001 - give restart info
C       0010 - read coefficients
C       0100 - read eigenvalues
C       1000 - read boson irrep identification
C
C     IDIM(1,IFRP) - number of positronic solutions
C     IDIM(2,IFRP) - number of electronic solutions
C     IDIM(3,IFRP) - number of AO-basis functions
C
C     Written by T.Saue Sept 1 1995
C
C***********************************************************************

#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL INFO,LBIT,TEST
      CHARACTER TEXT*74,CMOFIL*6,FMT*6,MXFORM*6
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT)
      DIMENSION IDIM(3,2)
      real(8), allocatable :: kc(:)
      real(8), allocatable :: ko(:)
      integer, allocatable :: ki(:)
C
      IF(IOPT.EQ.0) RETURN
      REWIND IUNIT
C
C     Read title line
C     ===============
C
      READ (IUNIT,END=10,ERR=20) TEXT,NSYM,
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG

      IF(IOPT.EQ.-1) RETURN
C
C     Compatibility check for NFSYM and NFBAS
C     =======================================
C
      IERR = 0
      IF(NSYM .NE. NFSYM) THEN
         IERR = IERR + 1
         INFO = .TRUE.
         WRITE(LUPRI,'(/A/A,2I5)')
     &      'REACMO FATAL ERROR: Incompatible number of fermion ircops',
     &      ' from file and in this run:',NSYM,NFSYM
      ELSE
         DO IFRP = 1,NFSYM
           IF(IDIM(3,IFRP).NE.NFBAS(IFRP,0)) IERR = IERR + 1
         END DO
         IF (IERR.GT.0) THEN
           INFO = .TRUE.
           WRITE(LUPRI,'(/A/A,2(/A,I2,2I5))')
     &     'REACMO FATAL ERROR. Incompatible number of basis functions',
     &     ' from file and in this run:', (' Fermion ircop no.',IFRP,
     &     IDIM(3,IFRP),NFBAS(IFRP,0),IFRP=1,NFSYM)
         END IF
      END IF
C
C     File info
C     =========
C
      IF(LBIT(IOPT,1)) THEN
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/3A,'//FMT//')')
     &   '* REACMO: Coefficients read from file ',CMOFIL,
     &   ' - Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
C
      IF (IERR.GT.0) CALL QUIT(
     &'REACMO: old DFCOEF file is not consistent with this calculation')
C
      IF(IOPT.LT.2) RETURN
C
C     Compatibility tests for number of orbitals
C     ==========================================
C
      TEST = .FALSE.
      NC = 0
      NO = 0
      DO IFRP = 1,NFSYM
        IDIM(3,IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
        IF(IDIM(3,IFRP).NE.NORB(IFRP)) TEST=.TRUE.
        NC = NC + IDIM(3,IFRP)*NFBAS(IFRP,0)*NZ
        NO = NO + IDIM(3,IFRP)
      ENDDO
C
C     Read coefficients
C     =================
C
      IF(LBIT(IOPT,2)) THEN
        IF(TEST) THEN
          allocate(kc(nc))
          CALL READT(IUNIT,NC, kc)
          JC = 1
          DO IFRP = 1,NFSYM
            CALL ADACMO(NFBAS(IFRP,0),NZ,
     &         CMO(ICMOQ(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         kc(JC),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            JC = JC + IDIM(3,IFRP)*NFBAS(IFRP,0)*NZ
          ENDDO
          deallocate(kc)
        TEST = .FALSE.
        ELSE
          READ (IUNIT,END=30,ERR=40) CMO
        ENDIF
        IF(IOPT.LT.4) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
C
C     Read eigenvalues
C     ================
C
      IF(LBIT(IOPT,3)) THEN
        IF(TEST) THEN
          allocate(ko(no))
          CALL READT(IUNIT,NO,ko)
          JO = 1
          DO IFRP = 1,NFSYM
            CALL ADAEIG(
     &         EIG(IORB(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         ko(JO),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            JO = JO + IDIM(1,IFRP)+IDIM(2,IFRP)
          ENDDO
          deallocate(ko)
        ELSE
          READ (IUNIT,END=50,ERR=60) EIG
        ENDIF
        IF(IOPT.LT.8) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
C
C     Read single group irrep (spinfree and Levy-Leblond calcs)
C     =========================================================
C
      IF(LBIT(IOPT,4)) THEN
        IF(TEST) THEN
          allocate(ki(no))
          CALL READI(IUNIT,NO,ki)
          CALL ADAIBP(IBEIG,ki,IDIM)
          deallocate(ki)
        ELSE
          READ (IUNIT,END=70,ERR=80) IBEIG
        ENDIF
        IF(IOPT.LT.16) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('REACMO: ERROR reading TEXT')
 30   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading coefficients')
 40   CONTINUE
      CALL QUIT('REACMO: ERROR reading coefficients')
 50   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading eigenvalues')
 60   CONTINUE
c     CALL QUIT('REACMO: ERROR reading eigenvalues')
      WRITE(LUPRI,'(/A)') 'REACMO: WARNING !!! No eigenvalues read'
      RETURN
 70   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading irrep identification')
 80   CONTINUE
      CALL QUIT('REACMO: ERROR reading irrep identification')
      RETURN
      END
      SUBROUTINE INDEXX(N,ARR,INDX)
C***********************************************************************
C
C     This subroutine goes through array ARR and returns the ordering
C     index array INDX giving indices from smallest to largest
C
C     The array ARR is not touched.
C
C***********************************************************************
      INTEGER N,INDX(N),M,NSTACK
      REAL*8 ARR(N)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,INDXT,IR,ITEMP,J,JSTACK,K,L,ISTACK(NSTACK)
      REAL*8 A
      DO J=1,N
         INDX(J)=J
      END DO
      JSTACK=0
      L=1
      IR=N
 1    IF(IR-L.GE.M)THEN
         K=(L+IR)/2
         ITEMP=INDX(K)
         INDX(K)=INDX(L+1)
         INDX(L+1)=ITEMP
         IF(ARR(INDX(L+1)).GT.ARR(INDX(IR)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ARR(INDX(L)).GT.ARR(INDX(IR)))THEN
            ITEMP=INDX(L)
            INDX(L)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ARR(INDX(L+1)).GT.ARR(INDX(L)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(L)
            INDX(L)=ITEMP
         ENDIF
         I=L+1
         J=IR
         INDXT=INDX(L)
         A=ARR(INDXT)
 3       CONTINUE
         I=I+1
         IF(ARR(INDX(I)).LT.A)GOTO 3
 4       CONTINUE
         J=J-1
         IF(ARR(INDX(J)).GT.A)GOTO 4
         IF(J.LT.I)GOTO 5
         ITEMP=INDX(I)
         INDX(I)=INDX(J)
         INDX(J)=ITEMP
         GOTO 3
 5       INDX(L)=INDX(J)
         INDX(J)=INDXT
         JSTACK=JSTACK+2
         IF(JSTACK.GT.NSTACK)
     $        CALL QUIT('NSTACK TOO SMALL IN INDEXX')
         IF(IR-I+1.GE.J-L)THEN
            ISTACK(JSTACK)=IR
            ISTACK(JSTACK-1)=I
            IR=J-1
         ELSE
            ISTACK(JSTACK)=J-1
            ISTACK(JSTACK-1)=L
            L=I
         ENDIF
      ELSE
         DO J=L+1,IR
            INDXT=INDX(J)
            A=ARR(INDXT)
            DO I=J-1,1,-1
               IF(ARR(INDX(I)).LE.A)GOTO 2
               INDX(I+1)=INDX(I)
            END DO
            I=0
 2          INDX(I+1)=INDXT
         END DO
         IF(JSTACK.EQ.0)RETURN
         IR=ISTACK(JSTACK)
         L=ISTACK(JSTACK-1)
         JSTACK=JSTACK-2
      ENDIF
      GOTO 1
      END

      LOGICAL FUNCTION LBIT(I,N)
C       Based on an analogous DISCO routine(Jan Almloef)
#include "implicit.h"
#include "ibtfun.h"
      LBIT = .FALSE.
      IF (IBTAND(IBTSHR(I,N-1),1).EQ.1) LBIT = .TRUE.
      END

      SUBROUTINE READT (IT,N,X)
#include "implicit.h"
      CHARACTER*30 STRING
      DIMENSION X(N)
      READ (IT,END=10,ERR=20) X
      RETURN
 10   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: END reading unit  ',IT
      CALL QUIT(STRING)
 20   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: Error reading unit',IT
      CALL QUIT(STRING)
      END

      SUBROUTINE QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  BUF,LBUF,IPRINT)
C***********************************************************************
C
C     This routine performs unitary transformations indicated by TYP:
C
C       TYP = AOMO : AO-to-MO-transformation  FMO = (C+)FAOC + FADD*FMO
C       TYP = MOAO : MO-to-AO-transformation  FAO = CFMO(C+) + FADD*FAO
C
C     FADD gives the possibility of adding the results
C
C     TREV indicates symmetry of F-matrix under timereversal
C
C       TREV = 'S' - symmetric     (t =  1)
C       TREV = 'A' - anti-symmtric (t = -1)
C
C     The unitary C-matrix has a time-symmetric structure.
C
C     The AO-to-MO-transformation can be set up quaternionically as
C
C       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
C
C       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
C
C     and the MO-to-Ao-transformation can be set up as:
C
C       FAO = [(Ca)+t(Cb)j][(FMOa)+(FMOb)j][(Ca+)-(Cb)j]
C
C       F(NORB,NORB,NZ) --->  F(NBAS,NBAS,NZ)
C
C     Written by T.Saue,January 1995
C     Last revision: 29 January 1997 : Luuk Visscher
C**********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      CHARACTER TREV*1,TRA*1,TYP*4
      DIMENSION FAO(LRAO,LCAO,*),IQAO(*),
     &          FMO(LRMO,LCMO,*),IQMO(*),
     &          TM1(LR1,LC1,*),IQTM1(*),
     &          TM2(LR2,LC2,*),IQTM2(*),
     &          BUF(LBUF)
      DIMENSION IQW(4)
CMI   .... IPRINT threshold 
      DATA IPRLIM / 200 /
!     DATA IPRLIM /   4 /

C
C
C     Add the results on to the contents of the MO matrix ?
C
#ifdef DEBUG
      print*, 'QTRANS: The AO matrix'
      CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
      print*, 'QTRANS: The transformation matrix TM1'
      CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)
      print*, 'QTRANS: The transformation matrix TM2'
      CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
#endif
C
C     Time reversal symmetry of F-matrix
C
      IF    (TREV.EQ.'S') THEN
        TRA = 'N'
      ELSEIF(TREV.EQ.'A') THEN
        TRA = 'I'
      ELSE
        CALL QUIT('QTRANS: Unknown keyword TREV '//TREV)
      ENDIF
C
C     *******************************************************
C     *****  AO - to - MO -transformation : F' = (C+)FC *****
C     *******************************************************
C
      IF(TYP.EQ.'AOMO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF (NCAO.GE.NRAO) THEN
C
C         Perform the two steps:
C           1. W  = FC
C           2. F' = (C+)W
C
          CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
         
          IF (IPRINT.GE.IPRLIM) THEN
            CALL HEADER('QTRANS: AO-to-MO-transformation',-1)

            WRITE(LUPRI,'(2X,A,4I2)')
     &      'QTRANS: The AO matrix; IQAO=',(IQAO(IZ),IZ=1,NZAO)
            CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)

            WRITE(LUPRI,'(/2X,A,4I2)')
     &      'QTRANS: The transformation matrix TM1; IQTM1=',
     &      (IQTM1(IZ),IZ=1,NZTM1)
            CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)

            WRITE(LUPRI,'(/2X,A,4I2)')
     &      'QTRANS: The transformation matrix TM2; IQTM2=',
     &      (IQTM2(IZ),IZ=1,NZTM2)
            CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
          ENDIF
C
C         First part of transformation: W = FC
C         ------------------------------------
C
          CALL QGEMM(NRAO,NCMO,NCAO,D1,
     &              'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &              'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &               D0,IQW,BUF,NRAO,NCMO,NZW)

          IF (IPRINT.GE.IPRLIM) THEN
            WRITE(LUPRI,'(/2X,A,4I2)')
     &      'QTRANS: First part transformation matrix W, W=FC;'//
     &      ' IQW=',(IQW(IZ),IZ=1,NZW)
            CALL PRQMAT(BUF,NRAO,NCMO,NRAO,NCMO,NZW,IQW,LUPRI)
          ENDIF
C
C         Second part of transformation: H = (C+)W
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCMO,NRAO,D1,
     &               'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &               'N','N',IQW,BUF,NRAO,NCMO,NZW,
     &               FADD,IQMO,FMO,LRMO,LCMO,NZMO)

          IF (IPRINT.GE.IPRLIM) THEN
            WRITE(LUPRI,'(/2X,A,4I2)')
     &     'QTRANS: Second part of transformation H,H=(C+).W; '//
     &     'IQMO=',(IQMO(IZ),IZ=1,NZMO)
            CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IQMO,LUPRI)
          ENDIF

        ELSE
C
C         Perform the two steps:
C           1. W  = (C+)F
C           2. F' = WC
C
          CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF

C
C         First part of transformation: W = (C+)F
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCAO,NRAO,D1,
     &       'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &       'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &        D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: H = WC
C         -------------------------------------
C
          CALL QGEMM(NRMO,NCMO,NCAO,D1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
C
        ENDIF
        IF(IPRINT.GE.5) THEN
          CALL HEADER('QTRANS: AO-to-MO-transformed Fock matrix',-1)
          CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IQMO,LUPRI)
        ENDIF
C
C     *******************************************************
C     *****  MO - to - AO -transformation : F' = CF(C+) *****
C     *******************************************************
C
      ELSEIF(TYP.EQ.'MOAO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF(NCMO.GE.NRMO) THEN
          CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = F(C+)
C
          CALL QGEMM(NRMO,NCAO,NCMO,D1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL QGEMM(NRAO,NCAO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ELSE
          CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = CF
C
          CALL QGEMM(NRAO,NCMO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       D0,IQW,BUF,NRAO,NCMO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL QGEMM(NRAO,NCAO,NCMO,D1,
     &       'N','N',IQW ,BUF,NRAO,NCMO,NZW,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ENDIF
       IF(IPRINT.GE.5) THEN
          CALL HEADER('QTRANS: MO-to-AO-transformed Fock matrix',-1)
          CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
        ENDIF
      ELSE
        CALL QUIT('FTRANS: Unknown keyword TYP '//TYP)
      ENDIF
C
      RETURN
C
      END

      SUBROUTINE QTRANS90(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,IPRINT)
c     Like QTRANS, but uses fortran 90 memory allocation internally
#include "implicit.h"
      real(8), allocatable :: work(:)
      CHARACTER TREV*1,TYP*4
      integer IQW(4)
C     find out how much buffer is needed in qtrans
      IF(TYP.EQ.'AOMO') THEN
         IF (NCAO.GE.NRAO) THEN
            CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
            NBUF = NRAO*NCMO*NZW
         else
            CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
            NBUF = NRMO*NCAO*NZW
         endif
      else
         IF(NCMO.GE.NRMO) THEN
            CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
            NBUF = NRMO*NCAO*NZW
         else
            CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
            NBUF = NRAO*NCMO*NZW
         endif
      endif

      if(nbuf.le.0)then
        call quit(' *** error in QTRANS90: buffer space allocation not
     & possible for length <= 0.***')
      end if

      allocate(work(nbuf))
!     allocate( work(nbuf), stat=ierr )
!     if (ierr.ne.0) call quit('QTRANS90: error in allocation !')
      work = 0
      call QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  work,nbuf,IPRINT)
      deallocate(work)
CMI    ... necessary to eliminate gfortran memory bug...
!     deallocate(work,  stat=ierr)
!     if (ierr.ne.0) call quit('QTRANS90: error in dallocation !')
      
      end
      SUBROUTINE QGEMM(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA,
     &                             FB,TB,IQPB,B,LRB,LCB,NZB,
     &                              BETA,IQPC,C,LRC,LCC,NZC)
C*****************************************************************************
C
C     Performs matrix-matrix operations on general matrices:
C
C     C := alpha * op(A) * op(B) + beta * C
C
C     where alpha and beta are real constants. A,B and C
C     are matrices that may be real(NZ = 1), complex(NZ = 2)
C     or quaternionic (NZ = 4). Calls BLAS routine DGEMM.
C
C     LRX and LCX are leading rows and columns of matrix X.
C     op(A) is an M by K matrix, op(B) is a K by N matrix
C     and C an M by N matrix.
C
C     op(X) is determined by the corresponding character variable
C     FX:
C       FX = 'N'        Normal matrix
C       FX = 'T'        Transpose..
C       FX = 'C'        Complex conjugate
C       FX = 'H'        Hermitian conjugate
C     and TX
C       TX = 'N'        No transform    a+bj
C       TX = 'I'        i-transform   -i(a+bj)i = a - bj
C       TX = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TX = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     AFUL and BFUL indicates on input what components of
C     matrices A and B are non-zero, whereas CFUL on output
C     indicates what components of C are non-zero.
C
C     Written by T.Saue, University of Oslo, Nov. 1994
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D1=1.0D0)
C
C     Global variables
C
      DIMENSION A(LRA*LCA,NZA),B(LRB*LCB,NZB),C(LRC*LCC,NZC),
     &          IQPA(NZA),IQPB(NZB),IQPC(NZC)
      CHARACTER  FA*1,FB*1,TA*1,TB*1,OA*1,OB*1
C
C     Local variables
C
      DIMENSION IPB(4)
      LOGICAL   LBUF(4)      
C
      EXTERNAL BDQALGEBRA
C
C     Initialize LBUF
C     ===============
C
      CALL LSET(4,.FALSE.,LBUF)
C
C     Make pointer IPB
C     ================
C
      IPB = 0
      DO IZB = 1,NZB
        IPB(IQPB(IZB)) = IZB
      ENDDO
C
C     MATRIX A
C

C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FA.EQ.'N') THEN
        OA  = 'N'
        IFA = 1
C     Transpose...
      ELSEIF(FA.EQ.'T') THEN
        OA  = 'T'
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        OA  = 'N'
        IFA = 2
C     Hermitian conjugate...
      ELSEIF(FA.EQ.'H') THEN
        OA  = 'T'
        IFA = 2
      ELSE
        CALL QUIT('QGEMM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
C     MATRIX B
C
C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FB.EQ.'N') THEN
        OB  = 'N'
        IFB = 1
C     Transpose...
      ELSEIF(FB.EQ.'T') THEN
        OB  = 'T'
        IFB = 1
C     Complex conjugate....
      ELSEIF(FB.EQ.'C') THEN
        OB  = 'N'
        IFB = 2
C     Hermitian conjugate...
      ELSEIF(FB.EQ.'H') THEN
        OB  = 'T'
        IFB = 2
      ELSE
        CALL QUIT('QGEMM:Unknown format '//FA//' of matrix B!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TB.EQ.'N') THEN
        ITB = 1
C     i-transform
      ELSEIF(TB.EQ.'I') THEN
        ITB = 2
C     j-transform
      ELSEIF(TB.EQ.'J') THEN
        ITB = 3
C     k-transform
      ELSEIF(TB.EQ.'K') THEN
        ITB = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TB '//TB//' of matrix B!')
      ENDIF
C
      DO 10 IZC = 1,NZC
        IQC = IQPC(IZC)
        DO 20 IZA = 1,NZA
          IQA = IQPA(IZA)
          IQB = IQMULT(IQA,IQC,1)
          IZB = IPB(IQB)
          IF(IZB.EQ.0) GOTO 20
          FACA = ALPHA*(IQSIGN(IQA,IFA,ITA)*IQSIGN(IQB,IFB,ITB)
     &                *IQPHASE(IQA,IQB,1))
          IF(LBUF(IZC)) THEN
            FACB = D1
          ELSE
            FACB = BETA
          ENDIF
          CALL DGEMM(OA,OB,M,N,K,FACA,A(1,IZA),LRA,B(1,IZB),LRB,
     &               FACB,C(1,IZC),LRC)
          LBUF(IZC) = .TRUE.
   20   CONTINUE
   10 CONTINUE
C
C     Check that matrix C is filled
C
      DO IZC = 1,NZC
      IF(.NOT.LBUF(IZC)) 
     &    WRITE(LUPRI,'(A,I5)') '* WARNING(QGEMM): Not filled: ',IZC
      ENDDO
      RETURN
C
      END

      SUBROUTINE IQPACK(IQPA,NZA,IQPB,NZB,IQPC,NZC)
C***********************************************************************
C
C     IQPA and IQPB are pointers to quaternion units 1,i,j and k
C     IQPACK considers what quaternion units span the quaternion product
C     A * B and makes packs these in pointer IQPC as generated.
C
C     Written by T.Saue - June 27 1996
C     Last revision: June 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "dcbqal.h"
C
      LOGICAL LBUF(4)
      DIMENSION IQPA(NZA),IQPB(NZB),IQPC(*)
      EXTERNAL BDQALGEBRA
C
C     Initialize LBUF
C     ===============
C
      LBUF(1) = .FALSE.
      LBUF(2) = .FALSE.
      LBUF(3) = .FALSE.
      LBUF(4) = .FALSE.
C
C     Determine quaternion units that contribute to A*B
C     =================================================
C
      NZC = 0
      DO IZA = 1,NZA
        IQA = IQPA(IZA)
        DO IZB = 1,NZB
          IQB = IQPB(IZB)
          IQC = IQMULT(IQA,IQB,1)
          IF(.NOT.LBUF(IQC)) THEN
            NZC = NZC + 1
            IQPC(NZC) = IQC
            LBUF(IQC) = .TRUE.
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
C
      END
      SUBROUTINE PRQMAT(QMAT,NROW,NCOL,LRQ,LCQ,NZ,IQP,IUNIT)
#include "implicit.h"
      DIMENSION QMAT(LRQ,LCQ,NZ),IQP(*)
      CHARACTER M(4)*6
      M(1) = 'A real'
      M(2) = 'A imag'
      M(3) = 'B real'
      M(4) = 'B imag'
      DO 10 IZ = 1,NZ
        IQ = IQP(IZ)
        WRITE(IUNIT,'(/A,A6,A,4I4/)')
     &  '*** ',M(IQ),' part ***; NROW,NCOL,LRQ,LCQ:',NROW,NCOL,LRQ,LCQ
        CALL OUTPUT(QMAT(1,1,IZ),1,NROW,1,NCOL,LRQ,LCQ,-1,IUNIT)
   10 CONTINUE
      END
      SUBROUTINE OPNFIL(IUNIT,NAME,STATUS,SUB)
C***********************************************************************
C
C       Open sequential, unformatted file.
C       Based on OPENDX
C       T.Saue, University of Tromsoe, Norway, Feb.-94
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) NAME,STATUS,SUB
      IF (STATUS .EQ. 'NEW') GO TO 10
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 100
C
C       Open old file
C
      OPEN(IUNIT,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',
     +     ACCESS='SEQUENTIAL',ERR=10)
C
C     Just in case....
      REWIND IUNIT
      GOTO 20
   10 CONTINUE
      IF(STATUS.EQ.'OLD') GOTO 110
C
C       Open new file
C
      OPEN(IUNIT,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     +     ACCESS='SEQUENTIAL')
   20 CONTINUE
      RETURN
C
C       Error messages
C
  100 CONTINUE
      WRITE(LUPRI,'(6A,I2,A)') SUB,': Invalid STATUS keyword ',STATUS,
     +     ' for file ',NAME,', (unit ',IUNIT,').'
      CALL QUIT('OPNFIL:Invalid STATUS keyword')
  110 CONTINUE
      WRITE(LUPRI,'(4A,I2,A)') SUB,': Old file ',NAME,' (unit',IUNIT,
     +    ') not found.'
      CALL QUIT('OPNFIL: Old file not found')
      END
      SUBROUTINE LSET(N,LVAL,LVEC)
#include "implicit.h"
      LOGICAL LVEC(N),LVAL
      DO 10 I = 1,N
        LVEC(I) = LVAL
   10 CONTINUE
      END

      SUBROUTINE ADACMO(NB,NZ,CF1,NO1,NP1,NE1,CF2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt coefficients CF2 to current dimensions CF1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CF1(NB,NO1,NZ),CF2(NB,NO2,NZ)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        NN = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)
      ENDDO
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)
      ENDDO
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        NN  = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
      ENDIF
C
      RETURN
      END
      SUBROUTINE ADAEIG(EIG1,NO1,NP1,NE1,EIG2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt eigenvalues EIG2 to current dimensions EIG1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG1(NO1),EIG2(NO2)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        CALL DZERO(EIG1(IV1),NR)
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        CALL DZERO(EIG1(IV1),NR)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adaibp */
      SUBROUTINE ADAIBP(IBO1,IBO2,IDIM)
C***********************************************************************
C
C     Pre-routine to adaibo
C     Written by T.Saue Feb 1 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION IBO1(*),IBO2(*),IDIM(3,2)
      I2 = 1
      DO IFRP = 1,NFSYM
        CALL ADAIBO(
     &         IBO1(IORB(IFRP)+1),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         IBO2(I2),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
          I2 = I2 + IDIM(2,IFRP)+IDIM(1,IFRP)
      ENDDO
C
      RETURN
      END

      SUBROUTINE ADAIBO(IBO1,NO1,NP1,NE1,IBO2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt boson irrep array IBO2 to current dimensions IBO1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION IBO1(NO1),IBO2(NO2)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        ibo1(iv1:iv1+nr) = 0
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        ibo1(iv1:iv1+nr) = 0
      ENDIF
C
      RETURN
      END


      SUBROUTINE ICOPY(N,IX,INCX,IY,INCY)
C
C     COPY integer IX TO integer IY.
C     FOR I = 0 TO N-1, COPY IX(LX+I*INCX) TO IY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
C     (860516 - hjaaj - based on BLAS DCOPY)
C
      INTEGER IX(*),IY(*)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      JX = 1
      JY = 1
      IF(INCX.LT.0)JX = (-N+1)*INCX + 1
      IF(INCY.LT.0)JY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        IY(JY) = IX(JX)
        JX = JX + INCX
        JY = JY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        IY(I) = IX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        IY(I) = IX(I)
        IY(I + 1) = IX(I + 1)
        IY(I + 2) = IX(I + 2)
        IY(I + 3) = IX(I + 3)
        IY(I + 4) = IX(I + 4)
        IY(I + 5) = IX(I + 5)
        IY(I + 6) = IX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          IY(I) = IX(I)
   70     CONTINUE
      RETURN
      END

C/* Deck Nucbas */
      SUBROUTINE NUCBAS(NUCORB,AA,IPRINT)
C***********************************************************************
C
C     Extract basis information for all centers
C
C     Written by T.Saue March 12 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "primit.h"
      CHARACTER SPDCAR*1
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)
C
C     Initialize
C
      NUCORB = 0
C
      JCENT = 0
      JPRIM = -1
      JC = -1
      DO ISHELL = 1,KMAX
        ICENT = NCENT(ISHELL)
        IF(ICENT.NE.JCENT) THEN
          JCENT = ICENT
          JLVAL = 0
        ENDIF
        IC = LCLASS(ISHELL)
        IF(IC.NE.JC) THEN
          JC    = IC
          JLVAL = 0
        ENDIF
        ILVAL = NHKT(ISHELL)
        IF(ILVAL.NE.JLVAL) THEN
          JLVAL = ILVAL
          NUCORB(ILVAL,IC,ICENT) = 0
          AA(1,ILVAL,IC,ICENT)=D0
          AA(2,ILVAL,IC,ICENT)=DUMMY
        ENDIF
        NUCORB(ILVAL,IC,ICENT)=NUCORB(ILVAL,IC,ICENT)+1
        IPRIM = JSTRT(ISHELL)
        IF(IPRIM.NE.JPRIM) THEN
          JPRIM = IPRIM
          NPRIM = NUCO(ISHELL)
          DO IEXP = 1,NPRIM
            A=PRIEXP(IPRIM+IEXP)
            AA(1,ILVAL,IC,ICENT)=MAX(AA(1,ILVAL,IC,ICENT),A)
            AA(2,ILVAL,IC,ICENT)=MIN(AA(2,ILVAL,IC,ICENT),A)
          ENDDO
        ENDIF
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         CALL HEADER('NUCORB:Basis set information:',-1)
        DO I = 1,NUCIND
          WRITE(LUPRI,'(/A,A4,A/)') '*** Center: ',NAMN(I),' ***'
          WRITE(LUPRI,'(2X,A)') '* Large components:'
          IC = 1
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H :',AA(1,LL,IC,I),
     &      'Alpha_L :',AA(2,LL,IC,I)
          ENDIF
          ENDDO
          WRITE(LUPRI,'(2X,A)') '* Small components:'
          IC = 2
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H: ',AA(1,LL,IC,I),
     &      'Alpha_L: ',AA(2,LL,IC,I)
          ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END

      SUBROUTINE READI (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
      DIMENSION INTX(N)
      READ (IT) INTX
      RETURN
      END

      SUBROUTINE LSQMA1(ICLAB)
C***********************************************************************
C
C     ... Make list of type (L or S) that maps each basis function
C     to a nuclei.
C
C     Written by Jesper K. Pedersen sep. 2001
C     
C***********************************************************************

#include <implicit.h>
#include <priunit.h>
C
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <frame.h>
C
#include <dcblab.h>
#include <symmet.h>
#include <aosotr.h>
#include <dgroup.h>
#include <pgroup.h>
#include <nuclei.h>
#include <shells.h>
#include <dcbbas.h>
#include <ccom.h>
      DIMENSION ICLAB(*)
#include <dcbibt.h>
C
C     Get array ICLAB which maps AO basis functions (large and small)
C     onto atomic centers.
C
      DO 100 I = 1, NTBAS(0)
         ITEMP = IATTR(IPLAB(I,1),1)
         ICENT = JGET(ITEMP)
         IDEG  = KGET(ITEMP)
         IF (IDEG .EQ. 1) THEN
            ICLAB(I) = NUCNUM(ICENT,1)
         ELSE
            JDEG = 0
            DO J = 0, MAXREP
               IF (NUCNUM(ICENT,J+1) .NE. 0) THEN
                  JDEG = JDEG + 1
                  IF (JDEG .EQ. IDEG) THEN
                     ICLAB(I) = NUCNUM(ICENT,J+1)
                     GO TO 100
                  END IF
               END IF
            END DO
         END IF
  100 CONTINUE
C
      RETURN
C      
      END

      FUNCTION MXFORM(X,NDIG)
C***********************************************************************
C
C       Find the format that gives optimum precision for real
C       variable X
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER MXFORM*6
C
      MXFORM = '      '
      IF(X.EQ.D0) THEN
        NDEC = NDIG
      ELSE
        NDEC = INT(LOG10(ABS(X)))
      ENDIF
      IF((NDEC.LT.0).OR.(NDEC.GT.(NDIG-3))) THEN
        MXFORM(1:1) = 'D'
        NDEC = NDIG - 7
      ELSE
        MXFORM(1:1) = 'F'
        NDEC = NDIG - NDEC - 3
      ENDIF
      IF(NDEC.LE.0) THEN
        WRITE(LUPRI,*) 'MXFORM out of boundaries...'
        WRITE(LUPRI,*) 'X=',X
        WRITE(LUPRI,*) 'NDEC,NDIG',NDEC,NDIG
        CALL QUIT('MXFORM: NDEC.LE.0')
      ENDIF
      IF(NDIG.LT.10) THEN
        IDIG = ICHAR('0') + NDIG
        MXFORM(2:2) = CHAR(IDIG)
        IOFF= 0
      ELSE
        IDIG = ICHAR('0') + INT(NDIG/10)
        MXFORM(2:2) = CHAR(IDIG)
        IDIG = ICHAR('0') + MOD(NDIG,10)
        MXFORM(3:3) = CHAR(IDIG)
        IOFF = 1
      ENDIF
      IND = 3+IOFF
      MXFORM(IND:IND) = '.'
      IF(NDEC.LT.10) THEN
        IDIG = ICHAR('0') + NDEC
        IND  = 4+IOFF
        MXFORM(IND:IND) = CHAR(IDIG)
        IOFF = IOFF+ 0
      ELSE
        IDIG = ICHAR('0') + 1
        IND = 4 + IOFF
        MXFORM(IND:IND) = CHAR(IDIG)
        IND = IND + 1
        IDIG = ICHAR('0') + MOD(NDEC,10)
        MXFORM(IND:IND) = CHAR(IDIG)
        IOFF = IOFF + 1
      ENDIF
      RETURN
      END

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

      subroutine fortwrt(str, slength)

      implicit none

#include "priunit.h"

#ifdef INT_STAR8
      integer(8) :: slength
#else
      integer    :: slength
#endif /* ifdef INT_STAR8 */
      character  :: str(slength)

      write(lupri, *) str

      end subroutine

      SUBROUTINE DFTSETCAM(W,BETA)
#include "implicit.h"
#include "dcbham.h"
      HFXMU = W
      HFXATT = BETA
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GTDOAV(mat_dim,
     &                  nz,
     &                  nr_dmat,
     &                  dmat,
     &                  ifac)
!***********************************************************************
!
!     Get average density from open shells; only real part !
!     Written by T.Saue March 2004
!
!***********************************************************************
      implicit none
#include"dcbdhf.h"
      integer :: mat_dim, nz, nr_dmat
      integer :: ifac, ishell
      real(8) :: DMAT(mat_dim, mat_dim, nz, nr_dmat)
      real(8) :: fac
      IF(AOC) THEN
        DO ISHELL = 1, NOPEN
          FAC = IFAC*DF(ISHELL)
          CALL DAXPY(mat_dim*mat_dim,FAC,
     &        DMAT(1,1,1,ISHELL+1),1,DMAT,1)
        ENDDO
      ENDIF
!
      end subroutine
