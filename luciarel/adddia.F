      SUBROUTINE ADDDIA(A,FACTOR,NDIM,IPACK)
*
* add factor to diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      DO 100 I = 1,NDIM
        IF(IPACK .EQ. 0 ) THEN
          II = (I-1)*NDIM + I
        ELSE
          II = I*(I+1)/2
        END IF
        A(II) = A(II) + FACTOR
  100 CONTINUE
*
      RETURN
      END
