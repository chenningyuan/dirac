      SUBROUTINE VECSUMC(C,A,B,FACAR,FACAI,FACBR,FACBI,NDIM)
*
* Multiply two complex vectors with scalars and add
*
*C(I) = FACA*A(I) + FACB*B(I)
*
* NDIM is length of real/imag part
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION C(NDIM,2)
*. Input
      DIMENSION A(NDIM,2),B(NDIM,2)
*
*
      DO I = 1, NDIM
*. Real part of C
        C(I,1) = FACAR * A(I,1)  + FACBR * B(I,1)
     &         - FACAI * A(I,2) -  FACBI * B(I,2)
*. Imaginary part of C
        C(I,2) = FACAR * A(I,2) + FACBR * B(I,2)
     &         + FACAI * A(I,1) + FACBI * B(I,1)
      END DO
*
      RETURN
      END
