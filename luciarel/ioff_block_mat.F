      FUNCTION IOFF_BLOCK_MAT(LROW,LCOL,IBLK)
*
* Offset of block IBLK in blocked matrix
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION LROW(*), LCOL(*)
*
      IOFF = 1
      DO JBLK = 1, IBLK-1
        IOFF = IOFF + LROW(JBLK)*LCOL(JBLK)
      END DO
*
      IOFF_BLOCK_MAT = IOFF
C?    WRITE(6,*) ' Offset for block ', IBLK, ' is ', IOFF
*
      RETURN
      END
