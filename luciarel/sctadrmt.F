      SUBROUTINE SCTADRMT(ISCT,SSCT,AIN,AOUT,NROWIN,NROWOUT,NCOL,FACTOR)
*
* AOUT(ISCT(I),J) = AOUT(ISCT(I),J)+SSCT(I)*AIN(I,J)*FACTOR if ISCT(I) .NE.0
* No update                                 if ISCT(I) .EQ.0
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION ISCT(*),SSCT(*),AIN(NROWIN,NCOL)
*.Output
      DIMENSION AOUT(NROWOUT,NCOL)
*
*. Currently programmed as stright loops,
* could be improved by processing batches of columns simultaneously
      DO J = 1, NCOL
        DO I =1, NROWIN
          IF(ISCT(I).NE.0) THEN
           AOUT(ISCT(I),J) = AOUT(ISCT(I),J) + SSCT(I)*FACTOR*AIN(I,J)
          END IF
        END DO
      END DO
*
      RETURN
      END
