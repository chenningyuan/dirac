      SUBROUTINE COPDIA(A,VEC,NDIM,IPACK)
*
* Copy diagonal of matrix A into vector VEC
*
*   IPACK = 0 : Full matrix
*   IPACK = 1 : Lower triangular matrix
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION A(*),VEC(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          VEC(I) = A((I-1)*NDIM+I)
  100   CONTINUE
      ELSE
        DO 200 I = 1, NDIM
          VEC(I) = A(I*(I+1)/2)
  200   CONTINUE
      END IF
*
      RETURN
      END
