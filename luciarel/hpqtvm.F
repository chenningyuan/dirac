      SUBROUTINE HPQTVM(VECIN,VECUT,IDUM,JDUM)
*
* Outer routine for Preconditioner times vector
*
*. LUCIA version
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "glbbas.inc"
#include "crun.inc"
#include "wrkspc.inc"
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,
C    &            ICISTR,NOCSF,IDIAG




      COMMON/SHFT/SHIFT
*
      NPDM = MXP1+MXP2
      NQDM = MXQ
      KLPHP = KH0
      KLPHQ = KH0 + NPDM*(NPDM+1)/2
      KLQHQ = KLPHQ + MXP1*MXQ
*
C          HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
      CALL HPQTV(MXP1,MXP2,MXQ,
     &           WORK(KLPHP),WORK(KLPHQ),WORK(KLQHQ),
     &           VECIN,VECUT)
*
      IF(SHIFT.NE.0.0D0) THEN
       CALL VECSUM(VECUT,VECUT,VECIN,1.0D0,SHIFT,NPDM+NQDM)
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTVM'
        CALL WRTMAT(VECIN,1,NPDM+NQDM,1,NPDM+NQDM)
        CALL WRTMAT(VECUT,1,NPDM+NQDM,1,NPDM+NQDM)
      END IF
*
      RETURN
      END
