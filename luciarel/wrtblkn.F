
      SUBROUTINE WRTBLKN(VEC,NBLOCK,LBLOCK)
*
* Write the NBLOCK blocks of VEC
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*)
      DIMENSION LBLOCK(NBLOCK)
*
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
        IF( LBLOCK(IBLOCK).GT.0) THEN
          WRITE(6,*) ' Block : ', IBLOCK
          WRITE(6,*) ' ==================='
          WRITE(6,*)
          WRITE(6,*) ' Length : ', LBLOCK(IBLOCK)
          WRITE(6,*)
          CALL WRITVE(VEC(IOFF),LBLOCK(IBLOCK))
          IOFF = IOFF + LBLOCK(IBLOCK)
          WRITE(6,*)
        END IF
      END DO
*
      RETURN
      END
