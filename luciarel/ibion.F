      INTEGER FUNCTION IBION(M,N)
C
C BIONOMIAL COEFFICIENT (M / N ) = IFAC(M)/(IFAC(M-N)*IFAC(N))
C
      IB = 1
      IF(M-N.GE.N) THEN
         DO 100 K = (M-N+1), M
           IB = IB * K
  100    CONTINUE
         IB = IB/IFAC(N)
      ELSE
         DO 200 K = N+1,M
           IB = IB * K
  200    CONTINUE
         IB = IB/IFAC(M-N)
      END IF
*
      IBION = IB
C
      RETURN
      END
