      REAL*8 FUNCTION INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
C
C DISC VERSION OF INPROD
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
      LOGICAL DIFVEC
C
      X = 0.0D0
      IF( LU1 .NE. LU2 ) THEN
        DIFVEC = .TRUE.
      ELSE
        DIFVEC =  .FALSE.
      END IF
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          IF(DIFVEC) REWIND LU2
         ELSE
          CALL REWINE( LU1,LBLK)
          IF( DIFVEC ) CALL REWINE( LU2,LBLK)
         END IF
      END IF
C
C LOOP OVER BLOCKS OF VECTORS
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          IF( DIFVEC) READ(LU2) NBL2
        ELSE IF ( LBLK .LT. 0 ) THEN
          CALL IFRMDS(NBL1,1,-1,LU1)
          IF( DIFVEC)CALL IFRMDS(NBL2,1,-1,LU2)
        END IF
C
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
          IF( DIFVEC) THEN
            CALL FRMDSC(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
            IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
          ELSE
          IF(NBL1 .GT. 0 )
     &    X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      INPRDD = X
C
      RETURN
      END
