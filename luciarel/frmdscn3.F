      SUBROUTINE FRMDSCN3(VEC,NREC,LBLK,LU,NO_ZEROING,I_AM_ZERO,
     &                    LBLOCK)
*
* Read  VEC as multiple record file, NREC records read
* If NO_ZEROING.EQ.0 then zero blocks are not
* set to zero;  a 1 is instead flagged in the relevant block
* of I_AM_ZERO
*
* IF LBLOCK(IBLK) is less than zero, then no space for this
* block is made
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION LBLOCK(*)
*. Output
      DIMENSION VEC(*),I_AM_ZERO(*)
*
      IOFF = 1
      DO IREC = 1, NREC
        CALL IFRMDS(LREC,1,LBLK,LU)
        IF(LBLOCK(IREC).GE.0) THEN
          CALL FRMDSC2(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          I_AM_ZERO(IREC) = IMZERO
          IOFF = IOFF + LREC
        ELSE
          I_AM_ZERO(IREC) = 1
          CALL SKPRCD2(LBLK,-1,LU)
        END IF

      END DO
*
      RETURN
      END
