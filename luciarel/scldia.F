      SUBROUTINE SCLDIA(A,FACTOR,NDIM,IPACK)
*
* scale diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*                assumed packed columnwise !!!!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      IF( IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          II = (I-1)*NDIM + I
          A(II) = A(II) * FACTOR
  100   CONTINUE
      ELSE
        II = 1
        DO 200 I = 1, NDIM
          A(II) = A(II) * FACTOR
          II = II + NDIM - I + 1
  200   CONTINUE
      END IF
*
      RETURN
      END
