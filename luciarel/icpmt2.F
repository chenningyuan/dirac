      SUBROUTINE ICPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
*
* Copy INTEGER matrix AIN to AOUT . Dimensions can differ
*
* If IZERO .ne. 0 , AOUT is zeroed  first
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input
      INTEGER AIN(NINR,NINC)
*. Output
      INTEGER AOUT(NOUTR,NOUTC)
*
      IF(IZERO.NE.0) CALL ISETVC(AOUT,0,NOUTR*NOUTC)
      DO 100 J = 1, NINC
       CALL ICOPVE(AIN(1,J),AOUT(1,J),NINR)
  100 CONTINUE
*
      RETURN
      END
