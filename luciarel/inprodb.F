      REAL*8 FUNCTION INPRODB(VEC1,VEC2,NBLK,LBLK,I0BLK)
*
* Inner products between blocked vectors with check of
* zero blocks. Zero blocks are flagged by a unit entry
* in I0BLK or a negative block length
*
*Jeppe Olsen, May 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC1(*),VEC2(*), LBLK(NBLK),I0BLK(NBLK)
*
      REAL*8 INPROD
*
      X = 0.0D0
      IOFF = 1
      DO IBLK = 1, NBLK
        LLBLK = LBLK(IBLK)
C?      WRITE(6,*) ' INPRODB IBLK LLBLK IOFF ',IBLK,LLBLK,IOFF
        IF(I0BLK(IBLK).EQ.0.AND.LLBLK.GT.0) THEN
          X = X + INPROD(VEC1(IOFF),VEC2(IOFF),LLBLK)
C?        WRITE(6,*) ' Vec1 and Vec2 blocks '
C?        CALL WRTMAT(VEC1(IOFF),1,LLBLK,1,LLBLK)
C?        CALL WRTMAT(VEC2(IOFF),1,LLBLK,1,LLBLK)
C?        write(6,*) ' Updated x', X
        END IF
        IF(LLBLK.GT.0) IOFF = IOFF + LLBLK
      END DO
*
      INPRODB = X
*
      RETURN
      END
