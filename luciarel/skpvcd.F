      SUBROUTINE SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
C
C SKIP OVER NVEC VECTORS ON FILE LUIN
C
C LBLK DEFINES STRUCTURE OF FILE
C (see note on structure of files )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE(LU ,LBLK)
      END IF
      DO 1001 IVEC = 1, NVEC
C
C LOOP OVER BLOCKS OF GIVEN VECTOR
C
 1000   CONTINUE
C
          IF( LBLK .GT. 0 ) THEN
            LBL = LBLK
          ELSE IF (LBLK .EQ. 0 ) THEN
            READ(LU) LBL
          ELSE IF (LBLK .LT. 0 ) THEN
            CALL IFRMDS(LBL,1,-1,LU)
          END IF
C
          IF( LBL .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = LBLK
            ELSE
              KBLK = -1
            END IF
            CALL FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,IAMPACK)
          END IF
        IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
 1001 CONTINUE
C
      RETURN
      END
