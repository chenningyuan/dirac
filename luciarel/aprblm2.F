      SUBROUTINE APRBLM2(A,LROW,LCOL,NBLK,ISYM)
C
C PRINT BLOCKED MATRIX
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
      DIMENSION LROW(NBLK),LCOL(NBLK)
C
      IBASE = 1
      WRITE(6,*) ' Blocked matrix '
      WRITE(6,*) '================'
      WRITE(6,*)
      DO 100 IBLK = 1, NBLK
        WRITE(6,'(A,I3)') ' Block ... ',IBLK
        IF(ISYM.EQ.0) THEN
        IF(IBLK .NE. 1 ) IBASE = IBASE + LROW(IBLK-1)*LCOL(IBLK-1)
        CALL WRTMAT(A(IBASE),LROW(IBLK),LCOL(IBLK),
     &              LROW(IBLK),LCOL(IBLK) )
       ELSE
        IF(IBLK .NE. 1 )
     &  IBASE = IBASE + LROW(IBLK-1)*(LCOL(IBLK-1)+1)/2
        CALL PRSYM(A(IBASE),LROW(IBLK))
       END IF

  100 CONTINUE
      RETURN
      END
