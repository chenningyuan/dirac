CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck dftexc */                                                C
      SUBROUTINE GRTEST(WORK,LWORK,JPRINT)
C                                                                      C
C     T. Helgaker sep 99                                               C
C     DIRAC adaption by T.Saue                                         C
C     Modified for grid testing by Olav Fossgaard July 2002.           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
 
      PARAMETER (DP5 = 0.5D0)
      DIMENSION WORK(LWORK)

#include "dcbgen.h"
#include "dcbham.h"
#include "inforb.h"
#include "nuclei.h"
#include "dftcom.h"
#include "memint.h"

CTROND      LUQUAD = 29
      LUQUAD = 44
C
C     Allocate SMAT and send to  DFTDRV
C
      CALL MEMGET('REAL',KSMAT,NBAST*NBAST*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO,NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA1,3*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGA2,6*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGB1,3*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCNT,NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDST,NATOMS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,NBAST*4,WORK,KFREE,LFREE)
      CALL GRTES1(WORK(KSMAT),WORK(KGAO),WORK(KGA1),WORK(KGA2),
     &     WORK(KGB1),WORK(KCNT),WORK(KDST),WORK(KBUF),
     &     WORK(KFREE),LFREE,LUQUAD,IPRDFT)
      CALL MEMREL('GRTES1',WORK,KWORK,KWORK,KFREE,LFREE)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C /* Deck dftdrv */
      SUBROUTINE GRTES1(SMAT,GAO,GAO1,GAO2,GAB1,
     &                  NCNT,DST,BUF,WORK,LWORK,LUQUAD,IPRINT)
C                                                                      C
C     T. Helgaker sep 99                                               C
C     Modified for grid testing by Olav Fossgaard July 2002.           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)

      CHARACTER*12 MXLB1,MXLB2,RMXLB1,RMXLB2
      INTEGER MAXI, MAXJ, MAXIR, MAXJR, NLSOS, NSSOS        
      DOUBLE PRECISION ERRAV, ERRAB, ERRMX, RERRAV, RERRAB, RERRMX,
     &                 SMATDF, SMATAN, SMATFR

#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "mxorb.h" 
#include "dcblab.h"
#include "maxaqn.h" 
#include "symmet.h"
#include "nuclei.h"
#include "dfterg.h"
#include "dftcom.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "functionals.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST,2)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0),2)
#endif

      DIMENSION GAO(*), GAO1(*), 
     &          GAO2(*), GAB1(*), WORK(LWORK), 
     &          NCNT(*), DST(*), BUF(*)

#include "memint.h"

C
C     Get analytical SMAT(I,J,2)
C
      NELM = NTBAS(0)*NTBAS(0)
      CALL DZERO(SMAT(1,1,1),NELM)
      CALL GTOVLX(SMAT(1,1,2),SSMTRC,WORK(KFREE),LFREE)
C
C     integrate
C
      REWIND LUQUAD
      R2TOT = D0
      NPNTS = 0
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL MEMGET('REAL',KX,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KY,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KZ,NPOINT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KW,NPOINT,WORK,KFREE,LFREE)
         CALL GRTES3(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),LUQUAD,NPOINT)

      CALL GRTES2(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NPOINT,
     &               SMAT,DFTASC,
     &               GAO,GAO1,GAO2,GAB1,
     &               NCNT,DST,BUF,WORK(KFREE),LFREE,IPRINT)

         CALL MEMREL('GRTES2',WORK,KWORK,KWORK,KFREE,LFREE)
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      ELSE
         GO TO 200
      END IF
  200 CONTINUE

C
C     Subtract the analytical SMAT(I,J,2) from the numerical SMAT(I,J,1)
C     The numerical SMAT is overwritten.
C
      CALL DAXPY(NELM,-1.0D0,SMAT(1,1,2),1,SMAT(1,1,1),1)
C
C     Analyze the difference matrix:
C
C     ERRAV = average error
C     ERRAB = absolute average werror
C     ERRMX = largest error
C     (MAXI, MAXJ) = overlap integral with error ERRMAX

      ERRAV = 0.0D0
      ERRAB = 0.0D0
      ERRMX = 0.0D0
      MAXI = 0
      MAXJ = 0
C
C     Relative values (only calculated for diagonal elements)
C
      RERRAV = 0.0D0
      RERRAB = 0.0D0
      RERRMX = 0.0D0
      MAXIR = 0
      MAXJR = 0
C
      SMATDF = 0.0D0
      SMATAN = 0.0D0
      SMATFR = 0.0D0
C
C     Loop over SO's
C
      IF (LEVYLE) THEN
         MC = 1
      ELSE
         MC = 2
      END IF
C
      DO IC = 1,MC
         DO IRP = 0,MAXREP
            DO I = 1,NBBAS(IRP,IC)
               IS = IBBAS(IRP,IC) + I
               DO J = 1,NBBAS(IRP,IC)
                  JS = IBBAS(IRP,IC) + J
                  SMATDF = SMAT(IS,JS,1)
                  SMATAN = SMAT(IS,JS,2)
C
                  ERRAV = ERRAV + SMATDF
                  ERRAB = ERRAB + ABS(SMATDF)
                  IF(ABS(SMATDF) .GT. ABS(ERRMX)) THEN
                     ERRMX = SMATDF
                     MAXI = IS
                     MAXJ = JS
                  END IF
C   
                  IF(JS .EQ. IS) THEN
                     SMATFR = SMATDF/SMATAN
                     RERRAV = RERRAV + SMATFR
                     RERRAB = RERRAB + ABS(SMATFR)
C
COLAV                     MXLB1 = PLABEL(IPLAB(is,2),2)
COLAV                     MXLB2 = PLABEL(IPLAB(js,2),2) 
COLAV                     write(lupri,'(2(I4,A12),E30.20,2F10.4)') 
COLAV     &                    IS, MXLB1, JS, MXLB2,SMATAN, 
COLAV     &                    LOG10(ABS(SMATFR)),
COLAV     &                    LOG10(ABS(D1-SMATAN))
C                    
                     IF (ABS(SMATFR) .GT. ABS(RERRMX)) THEN
                        RERRMX = SMATFR
                        MAXIR = IS
                        MAXJR = JS
                     END IF   
                  END IF
C   
               END DO
            END DO
         END DO
      END DO
C
C     Get labels of maximum error overlap SO's
C
      MXLB1 = PLABEL(IPLAB(MAXI,2),2)
      MXLB2 = PLABEL(IPLAB(MAXJ,2),2) 
      RMXLB1 = PLABEL(IPLAB(MAXIR,2),2)
      RMXLB2 = PLABEL(IPLAB(MAXJR,2),2)
C
C     Calculate average errors:
C
      NLSOS = 0
      NSSOS = 0

      DO I = 1, MAXREP+1
         NLSOS = NLSOS + NCOS(I,1)*NCOS(I,1)
         NSSOS = NSSOS + NCOS(I,2)*NCOS(I,2) 
      END DO
      
      IF(LEVYLE) THEN
         ERRAV = ERRAV/DBLE(NLSOS)
         ERRAB = ERRAB/DBLE(NLSOS)
         RERRAV = RERRAV/DBLE(NTBAS(1))
         RERRAB = RERRAB/DBLE(NTBAS(1))
      ELSE
         ERRAV = ERRAV/(DBLE(NLSOS) + DBLE(NSSOS))
         ERRAB = ERRAB/(DBLE(NLSOS) + DBLE(NSSOS))
         RERRAV = RERRAV/DBLE(NTBAS(0))
         RERRAB = RERRAB/DBLE(NTBAS(0))
      END IF
C
C     Print section
C
      CALL HEADER('Output from GRTEST', -1)
C
      WRITE(LUPRI,'(/A)') '* Errors: '
C
      WRITE (LUPRI,'( 2X,A,E10.4)')
     &      ' Average SMAT error               : ', ERRAV
      WRITE (LUPRI,'( 2X,A,E10.4)')
     &      ' Absolute average SMAT error      : ', ERRAB      
      WRITE (LUPRI,'( 2X,A,E10.4)')
     &      ' Largest error                    : ', ERRMX
      WRITE (LUPRI,'( 2X,A,A,I4,A,I4,A)')
     &      ' Position of largest error        : ', 
     &      '(',MAXI,',',MAXJ,')'
      WRITE (LUPRI,'( 2X,A,A,A12,A,A12,A)')
     &      ' Labels for largest error integral: ', 
     &      '(',MXLB1,',',MXLB2,')'
C
      WRITE(LUPRI,'(/A)')'* Relative errors (diagonal elements only) : '
C
      WRITE (LUPRI,'( 2X,A,E10.4)')
     &      ' Average SMAT error               : ', RERRAV
      WRITE (LUPRI,'( 2X,A,E10.4)')
     &      ' Absolute average SMAT error      : ', RERRAB      
      WRITE (LUPRI,'( 2X,A,E10.4)')
     &      ' Largest error                    : ', RERRMX
      WRITE (LUPRI,'( 2X,A,A,I4,A,I4,A)')
     &      ' Position of largest error        : ', 
     &      '(',MAXIR,',',MAXJR,')'
      WRITE (LUPRI,'( 2X,A,A,A12,A,A12,A)')
     &      ' Labels for largest error integral: ', 
     &      '(',RMXLB1,',',RMXLB2,')'

      IF (INTCHK.GE.2) THEN
         CALL HEADER('Errors in numerical SMAT ', -1)
         CALL OUTPUT(SMAT(1,1,1),1,NTBAS(0),1,NTBAS(0),NTBAS(0),
     &        NTBAS(0),1,LUPRI)
      END IF

      RETURN 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck dftdr1 */                                                C
      SUBROUTINE GRTES2(CORX,CORY,CORZ,WEIGHT,NPOINT,
     &               SMAT,DFTASC,
     &               GAO,GAO1,GAO2,GAB1,
     &               NCNT,DST,BUF,
     &               WORK,LWORK,IPRINT)
C                                                                      C
C     Written by Trond Saue and T.Helgaker                             C
C     Modified for grid testing by Olav Fossgaard July 2002.           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)

#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "dftbuf.h"
#include "dfterg.h"
#include "functionals.h"
      
      LOGICAL DFTASC
      DIMENSION SMAT(*)
      DIMENSION CORX(NPOINT),CORY(NPOINT),CORZ(NPOINT),
     &            WEIGHT(NPOINT),
     &            GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &            GAB1(NBAST,3), 
     &            NCNT(NBAST), DST(NATOMS),BUF(NBAST,4),
     &            VX(DERVS1),
     &            WORK(LWORK)
      DIMENSION RHG(3)

      DO 300 IPNT = 1, NPOINT
         WGHT = WEIGHT(IPNT)
C
C       AOs
C       ===
C
        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,
     &              CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &              BUF,WORK,LWORK,NBAST,NDER,.FALSE.,IPRINT) 
C
C       Test accuracy of numerical integration
C       ======================================
C
        CALL GRTES4(WGHT,0,SMAT,GAO)

  300 CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck reaqua */
      SUBROUTINE GRTES3(X,Y,Z,WT,LUQUAD,NPOINT)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck doexc1 */
      SUBROUTINE GRTES4(VXC,IREPOP,EXCMAT,GAO)      
C
C     Construct Kohn-Sham matrix                                       C
C                                                                      C
C     Written by T.Saue and T.Helgaker                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"

      PARAMETER(D0=0.0D0,D2=2.0D0)

#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))

      DIMENSION GAO(NBAST),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0))
#endif
      
#include "ibtfun.h"

      MC = 2
      IF(LEVYLE) MC = 1      
      DO 10 IRP1 = 0,NBSYM-1
        DO 20 IRP2 = 0,NBSYM-1
          IF(IBTXOR(IRP1,IRP2).NE.IREPOP) GOTO 20
          DO 30 IC = 1,MC
            NB1 = NBBAS(IRP1,IC)
            IF(NB1.EQ.0) GOTO 30
             NB2 = NBBAS(IRP2,IC)
            IF(NB2.EQ.0) GOTO 30
            JB1 = IBBAS(IRP1,IC)
            JB2 = IBBAS(IRP2,IC)
            DO J = 1, NBBAS(IRP2,IC)
              GVXC = VXC*GAO(JB2+J)
              CALL DAXPY(NBBAS(IRP1,IC),GVXC,
     &                   GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
            ENDDO
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      END


