      FUNCTION IFREQ(IVEC,IVAL,NDIM)
*
* Number of times IVAL occurs in IVEC
*
      DIMENSION IVEC(*)
*
      NTIME = 0
      DO 100 I = 1, NDIM
        IF(IVEC(I).EQ.IVAL) NTIME = NTIME + 1
  100 CONTINUE
*
      IFREQ = NTIME
*
      RETURN
      END
