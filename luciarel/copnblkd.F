      SUBROUTINE COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*
* Copy Nblocks from FILE LUIN to file LUOUT
*
*
* Type of file LUOUT is inherited from LUIN
*
* If less than NBLOCK blocks are left in the given vector
* on LUIN a warning is printed,
* but normal return
*
*
* NO EOV mark written
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF

      DO IIBLK = 1, NBLOCK
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          IF(LBL.GE.0) CALL ITODS (LBL,1,-1,LUOUT)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IAMPACK.EQ.0) THEN
            CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
          ELSE
            CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
          END IF
        END IF
        IF(LBL.EQ.-1) THEN
           WRITE(6,*) ' END OF VECTOR mark for BLOCK'
           write(6,*) 'IBLOCK not defined!!'
           stop
           GOTO 3333
        END IF
      END DO
*     ^End of loop over blocks
 3333 CONTINUE
*
      RETURN
      END
