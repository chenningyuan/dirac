      SUBROUTINE DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
C
C  IF( INV .NE. 0 ) THEN
C    V3(I) = (V1(I)+FAC)-1 * V2(I)
C    LU3      LU1            LU2
C  IF( INV .EQ. 0 ) THEN
C    V3(I) = (V1(I)+FAC) * V2(I)
C    LU3         LU1        LU2
C WHERE V1 AND V2 ARE VECTORS ON FILES LU1 AND LU2,
C AND LU3 IS WRITTEN ON FILE LU3
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  VEC1(*),VEC2(*)
C
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE(LU1,LBLK)
          CALL REWINE(LU2,LBLK)
          CALL REWINE(LU3,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          CALL IFRMDS(LBL2,1,-1,LU2)
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVCD : '
     &                     , LBL1,LBL2
          Call Abend2( ' DIFFERENT BLOCKSIZES IN DMTVCD ' )
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(VEC1,LBL1,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,LBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
             CALL DIAVC2(VEC2,VEC2,VEC1,FAC,LBL1)
            ELSE
             CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
             CALL VECSUM(VEC2,VEC1,VEC2,1.0D0,FAC,LBL1)
           END IF
C          CALL TODSC(VEC2,LBL1,KBLK,LU3)
         END IF
         CALL TODSC(VEC2,LBL1,KBLK,LU3)
      END IF
C
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
