      SUBROUTINE FIND_ACTIVE_CLASSES(LUIN,
     &           LBLK,I_BLK_TO_CLS,ICLS_A,NCLS,SEGMNT)
*
*. Find the active (nonvanishing classes ) on LUIN)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION I_BLK_TO_CLS(*)
*. Output
      DIMENSION ICLS_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      CALL REWINE(LUIN  ,LBLK)
      IZERO = 0
      CALL ISETVC(ICLS_A,IZERO,NCLS)
*
      IBLK = 0
      NCLS_A = 0
*. Loop over blocks
C?      write(6,*) ' ZAP_BLOCK_VEC :  LBLK = ', LBLK
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
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN
           ICLS_A(I_BLK_TO_CLS(IBLK)) = 1
C          write(6,*) ' Active block and class ',IBLK,
C    &     I_BLK_TO_CLS(IBLK)
           NCLS_A = NCLS_A + 1
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
      NTEST = 1
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' FIND_A.... Number of active classes ',NCLS_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active classes '
        CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
      END IF
*
      RETURN
      END
