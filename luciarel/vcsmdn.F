      FUNCTION  VCSMDN(VEC1,VEC2,FAC1,FAC2,LU1,LU2,IREW,LBLK)
*
* Norm of sum of two vectors residing on disc
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPROD
*
      XNORM = 0.0D0
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
      END IF
*
* LOOP OVER BLOCKS OF VECTOR
*
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &  NBL1,NBL2
        Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
      END IF
C
      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        CALL FRMDSC(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
        CALL FRMDSC(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
        IF( NBL1 .GT. 0 ) THEN
          CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          XNORM = XNORM + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
*
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      VCSMDN = XNORM
      RETURN
      END
