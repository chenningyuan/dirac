      SUBROUTINE ZAP_BLOCK_VEC(LUIN,LBLK,IBLKS_A,SEGMNT,LUSCR)
*
* Zap blocks in vector in file LUIN for which IBLKS_A is zero
*
* Vector is initially constructed on LUSCR, and is copied back to
* LUIN after use
*
* Note : Files are always rewinded
* Packed version
*
* LBLK defines file type
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
      INTEGER IBLKS_A(*)
*
      CALL REWINE(LUIN  ,LBLK)
      CALL REWINE(LUSCR ,LBLK)
*
      IBLK = 0
*. Loop over blocks
C?      write(6,*) ' ZAP_BLOCK_VEC :  LBLK = ', LBLK
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUSCR) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUSCR)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IBLKS_A(IBLK).EQ.0) THEN
            ZERO = 0.0D0
            CALL SETVEC(SEGMNT,ZERO,LBL)
          END IF
          CALL TODSCP(SEGMNT,LBL,KBLK,LUSCR)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
*. And then copy back to LUIN
C          COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      IREW = 1
      CALL COPVCDP(LUSCR,LUIN,SEGMNT,IREW,LBLK)
*
      RETURN
      END
