      SUBROUTINE ICOPVE2(IIN,IOFF,NDIM,IOUT)
*
* IOUT(I) = IIN(IOFF-1+I),I = 1, NDIM
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION IIN(*)
*. Output
      DIMENSION IOUT(*)
*
      DO I = 1, NDIM
        IOUT(I) = IIN(IOFF-1+I)
      END DO
*
      RETURN
      END
