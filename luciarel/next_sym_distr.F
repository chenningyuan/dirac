      SUBROUTINE NEXT_SYM_DISTR(NGAS,MINVAL,MAXVAL,
     &           ISYM,ISYM_TOT,IFIRST,NONEW)
*
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION MINVAL(NGAS),MAXVAL(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*. Symmetry of first NGAS -1 spaces
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 1, NGAS-1
          ISYM(IGAS) = MINVAL(IGAS)
        END DO
        IFIRST = 0
        NONEW = 0
      ELSE
C            NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISYM,NGAS-1,MINVAL,MAXVAL,NONEW)
      END IF
*. Symmetry of last space
      IF(NONEW.EQ.0) THEN
        JSYM = 1
        DO IGAS = 1, NGAS-1
          CALL SYMCOM(3,0,JSYM,ISYM(IGAS),KSYM)
          JSYM = KSYM
        END DO
        CALL SYMCOM(2,0,JSYM,ISYM(NGAS),ISYM_TOT)
C      write(6,*) ' after SYMCOM : JSYM, ISYM(NGAS), ISYM_TOT'
C      write(6,*) JSYM,ISYM(NGAS),ISYM_TOT
*
C       IF(MINVAL(NGAS).GT.ISYM(NGAS).OR.
C    &     MAXVAL(NGAS).LT.ISYM(NGAS)    )NONEW = 1
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
