      SUBROUTINE FIND_ACTIVE_BLOCKS
     &           (LUIN,LBLK,BLK_A,SEGMNT)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      CALL REWINE(LUIN  ,LBLK)
*
      IBLK = 0
      NBLK_A = 0
*. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
      NBLK =  IBLK-1
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,'(A,I8,I8)')
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active blocks '
        CALL WRTMAT(BLK_A,1,NBLK,1,NBLK)
      END IF
*
      RETURN
      END
