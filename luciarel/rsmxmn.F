      SUBROUTINE RSMXMN(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     &                  MIN1_T,MAX1,MIN3,MAX3_T,NTEST)
*
* Construct accumulated MAX and MIN arrays for a RAS set of strings
*
      IMPLICIT REAL*8           ( A-H,O-Z)
      DIMENSION  MINEL(*),MAXEL(*)
*
      NORB = NORB1 + NORB2 + NORB3
*. accumulated max and min in each of the three spaces
*. ( required max and min at final orbital in each space )
COLD  MIN1A = MIN1_T
      MIN1A = MAX(MIN1_T,NEL-MAX3_T-NORB2)
      MAX1A = MAX1
*
      MIN2A = NEL - MAX3_T
      MAX2A = NEL - MIN3
*
      MIN3A = NEL
      MAX3A = NEL
*
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF(NORB1 .GT. 0 )
     &    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF(NORB1+NORB2 .GT. 0 )
     &    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
  100 CONTINUE
*
      IF( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from RSMXMN '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
