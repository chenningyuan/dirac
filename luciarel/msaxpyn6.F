      SUBROUTINE MSAXPYN6(AX,A,X,NDIM,NVEC,INDEX)
*
* AX(I) = AX(I) + SUM(L=1,NVEC) A(L)*X(I,INDEX(L))
*
* Initial RISC/6000 version
*
* Jeppe Olsen
*

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AX(*),X(NDIM,*)
      DIMENSION A(*) ,INDEX(*)
*
* level of unrolling of column loop
      NROL = 5
      NRES = MOD(NVEC,NROL)
*
*. part that cannot be column unrolled
*
      IF(NRES.EQ.0.OR.NVEC.EQ.0) THEN
*. Nothing to do !
      ELSE IF (NRES.EQ.1) THEN
        I1 = INDEX(1)
        A1 = A(1)
        DO 20 I = 1, NDIM
          AX(I) = AX(I) + A1*X(I,I1)
   20   CONTINUE
      ELSE IF (NRES.EQ.2) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        DO 30 I = 1, NDIM
         AX(I) = AX(I) + A1*X(I,I1)+A2*X(I,I2)
   30   CONTINUE
      ELSE IF (NRES.EQ.3) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        I3 = INDEX(3)
        A3 = A(3)
        DO 40 I = 1, NDIM
         AX(I) =
     &   AX(I) + A1*X(I,I1)+A2*X(I,I2)+A3*X(I,I3)
   40   CONTINUE
      ELSE IF (NRES.EQ.4) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        I3 = INDEX(3)
        A3 = A(3)
        I4 = INDEX(4)
        A4 = A(4)
        DO 50 I = 1, NDIM
         AX(I) =
     &   AX(I) + A1*X(I,I1)+A2*X(I,I2)+A3*X(I,I3)+A4*X(I,I4)
   50   CONTINUE
      ELSE IF (NRES.GE.5) THEN
        WRITE(6,*) ' NRES to large in MSAXPY so STOP '
        Call Abend2( ' MSAXPY ' )
      END IF
*. unrolled column loops
      DO 285 L = 1+NRES,NVEC,NROL
        AL= A(L)
        IL= INDEX(L)
        AL1= A(L+1)
        IL1 = INDEX(L+1)
        AL2= A(L+2)
        IL2 = INDEX(L+2)
        AL3= A(L+3)
        IL3 = INDEX(L+3)
        AL4= A(L+4)
        IL4 = INDEX(L+4)
        DO 205 I = 1, NDIM
          AX(I) = AX(I) + AL*X(I,IL)+AL1*X(I,IL1)
     &          + AL2*X(I,IL2)+AL3*X(I,IL3)+AL4*X(I,IL4)
  205   CONTINUE
  285 CONTINUE
*
      RETURN
      END
