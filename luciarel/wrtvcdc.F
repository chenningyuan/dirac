      SUBROUTINE WRTVCDC(SEGMNT,LU,IREW,LBLK)
C
C PRINT COMPLEX VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU
        ELSE
          CALL REWINE(LU,LBLK)
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
      DO IRI = 1, 2
        IF(IRI.EQ.1) THEN
           WRITE(6,*) ' Real Part '
        ELSE
           WRITE(6,*) ' Imaginary part '
        END IF
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU) LBL
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,I_AM_PACKED)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      END DO
*. End of loop over real/imaginary parts
      RETURN
      END
