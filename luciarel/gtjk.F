      SUBROUTINE GTJK(RJ,RK,NTOOB,SCR,IREOTS)
*
* Interface routine for obtaining Coulomb (RJ) and
* Exchange integrals (RK)
*
* Ordering of intgrals is in the internal order
      IMPLICIT REAL*8(A-H,O-Z)
*
*.CRUN
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,NP1,NP2,NQ,INCORE,MXCIV,ICISTR,
C    &            NOCSF,IDIAG
#include "mxpdim.inc"
#include "crun.inc"
*.Input
      DIMENSION IREOTS(*)
*.Output
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*.Scratch
      DIMENSION SCR(2*NTOOB ** 2)
* Commented out
C     IF(INTIMP.EQ.1) THEN
      IF(INTIMP.EQ.11) THEN
*. Interface to MOTRA integrals
        CALL GTJKM(SCR(1),SCR(1+NTOOB**2))
*. Reorder to LUCIA order
        DO 100 J = 1, NTOOB
          JJ = IREOTS(J)
          DO 50 I = 1, NTOOB
            II = IREOTS(I)
C?          WRITE(6,*) ' I II J JJ ',I,II,J, JJ
            RJ(I,J) = SCR((JJ-1)*NTOOB+II)
            RK(I,J) = SCR((JJ-1)*NTOOB+II+ NTOOB **2)
   50     CONTINUE
  100  CONTINUE
      ELSE IF(INTIMP.EQ.1.OR.INTIMP.EQ.5) THEN
*. Interface to SIRIUS
        CALL GTJKS(RJ,RK,NTOOB)
      ELSE
*. Interface to LUCAS integrals
        CALL GTJKL(RJ,RK,NTOOB)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' RJ and RK from GTJK '
        CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
        CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
