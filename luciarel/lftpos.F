
      SUBROUTINE LFTPOS(CARD,LENGTH)
*
* Eliminate blanks and left position chararcter string CARD
*
      CHARACTER*72 CARD
*
      IEFF = 0
C     WRITE(6,*) ' INPUT string to LFTPOS '
C     WRITE(6,'(1H ,A)') CARD

      DO 100 IPOS = 1, LENGTH
       IF(CARD(IPOS:IPOS).NE.' ') THEN
         IEFF = IEFF + 1
         IF(IEFF.NE.IPOS) CARD(IEFF:IEFF) = CARD(IPOS:IPOS)
       END IF
  100 CONTINUE
*.Fill end with trailing blanks
      DO 200 IPOS = IEFF+1,LENGTH
        CARD(IPOS:IPOS) = ' '
  200 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Left adjusted character string '
        WRITE(6,'(1H ,A)') CARD
      END IF
*
      RETURN
      END
