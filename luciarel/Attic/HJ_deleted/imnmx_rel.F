      FUNCTION IMNMX_REL(IVEC,NDIM,MINMAX)
*
*     Find smallest (MINMAX=1) or largest (MINMAX=2)
*     absolute value of elements in integer vector IVEC
*
      DIMENSION IVEC(NDIM)
*
      IX = 0
      IF(MINMAX.EQ.1) THEN
        IF(NDIM.GT.0) IX=ABS(IVEC(1))
        DO 100 I=2,NDIM
          IX=MIN(IX,ABS(IVEC(I)))
  100   CONTINUE
      END IF
*
      IF(MINMAX.EQ.2) THEN
        IF(NDIM.GT.0) IX=ABS(IVEC(1))
        DO 200 I=2,NDIM
          IX=MAX(IX,ABS(IVEC(I)))
  200   CONTINUE
      END IF
*
      IMNMX_REL = IX
*
      RETURN
      END
