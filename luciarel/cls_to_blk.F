      SUBROUTINE CLS_TO_BLK(NBLOCK,IBLK_TO_CLS,ICLS_A,IBLK_A)
*
* an array ICLS_A is given for each class as well as
* a block to class array IBLK_TO_CLS.
* Obtain ICLS_A in block form (Well, I am bit tired and my
* pedagogical explanations can be pretty lousy even when I
* am awake !A)
*
* Jeppe Olsen, Jan. 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*.input
      DIMENSION ICLS_A(*),IBLK_TO_CLS(*)
*.output
      INTEGER IBLK_A(*)
*
      DO IBLK = 1, NBLOCK
        IBLK_A(IBLK) = ICLS_A(IBLK_TO_CLS(IBLK))
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CLS_TO BLK : IBLKS_a array '
        CALL IWRTMA(IBLK_A,1,NBLOCK,1,NBLOCK)
      END IF
*
      RETURN
      END
