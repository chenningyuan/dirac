      SUBROUTINE REWINE(LU,LBLK)
*
* LBLK .LT. 0 :  REWIND SEQ FILE LU WITH FASTIO ROUTINES
* LBLK .GE. 0 :  rewinf seq file LU with normal REWIND
      ICRAY = 1
      IF ( ICRAY.EQ.0.AND.LBLK .LT. 0 ) THEN
        IDUM = 1
        CALL SQFILE(LU,5,IDUM,IDUM)
      ELSE
        REWIND LU
      END IF
*
      RETURN
      END
