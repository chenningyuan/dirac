      SUBROUTINE NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* A set of numbers INUM(I),I=1,NELMNT is given.
* Find next compound number.
* Digit I must be in the range MINVAL(I),MAXVAL(I).
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen, Oct 1994
*
*. Input
      DIMENSION MINVAL(*),MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input number in NXTNUM '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      IF(NELMNT.EQ.0) THEN
        NONEW = 1
        GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            INUM(JPLACE) = MINVAL(JPLACE)
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
