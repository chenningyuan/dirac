      SUBROUTINE SKPRCD2(NDIM,MBLOCK,IFILE)
C
C     Skip record in file IFILE
C
*. Version allowing zero and packed blocks
*
* Dos not work with FASTIO - I expect
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)

C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
*. Loop over packed records of dimension LPBLK
*. The next LPBLK elements
  999   CONTINUE
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE)
            READ(IFILE)
          END IF
          READ(IFILE) ISTOP
        IF(ISTOP.EQ.0) GOTO 999
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
C       IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
C        CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
C       END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
