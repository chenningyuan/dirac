      SUBROUTINE ZERORC_REL(MBLOCK,IFIL)
*
* A record was known to be identical  zero
*
*
* Write corresponding info to file IFIL
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISCR(2)
* Zero record
      ISCR(1) = 1
*. Packed form
      ISCR(2) = 1
*
      CALL ITODS(ISCR,2,2,IFIL)
*
      RETURN
      END
