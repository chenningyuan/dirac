* Output
      SUBROUTINE NXTIJ(I,J,NI,NJ,IJSM,NONEW)
*
* An ordered pair (I,J) is given ,I.LE.NI,J.LE.NJ
*
* Find next pair, if IJSM .ne. 0 ,I .ge. J
*
      NONEW = 0
  100 CONTINUE
      IF(I.LT.NI) THEN
        I = I + 1
      ELSE
        IF(J.LT.NJ) THEN
          I = 1
          J = J+1
        ELSE
          NONEW = 1
          GOTO 101
        END IF
      END IF
      IF(IJSM.NE.0.AND.I.LT.J) GOTO 100
  101 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' next (i,j) pair ', I,J
      END IF
*
      RETURN
      END
