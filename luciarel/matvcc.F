      SUBROUTINE MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
*
* ITRNS = 0 : VOUT(I) = A(I,J)*VIN(J)
* ITRNS = 1 : VOUT(I) = A(J,I)*VIN(J)
*
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(NROW,NCOL)
       DIMENSION VIN(*),VOUT(*)
*
      IF(ITRNS.EQ.0) THEN
        DO 10 I = 1, NROW
   10   VOUT(I) = 0.0D0
        DO 100 J = 1, NCOL
         VINJ = VIN(J)
         DO 60 I = 1, NROW
           VOUT(I) = VOUT(I) + A(I,J)*VINJ
   60    CONTINUE
  100   CONTINUE
      ELSE IF( ITRNS.EQ.1) THEN
        DO 200 I = 1, NCOL
          X = 0.0D0
          DO 160 J = 1, NROW
            X = X + A(J,I)*VIN(J)
  160     CONTINUE
          VOUT(I) = X
  200   CONTINUE
      END IF
*
      RETURN
      END
