      SUBROUTINE IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
*
* Add two (scaled) integer vectors
*
*        IA(*) = IFACB*IB(*) + IFACC*IC(*)
*
      DIMENSION IA(*),IB(*),IC(*)
*
      DO 100 I = 1, NDIM
        IA(I) = IFACB * IB(I) + IFACC * IC(I)
  100 CONTINUE
*
      RETURN
      END
