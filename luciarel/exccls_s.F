      SUBROUTINE EXCCLS_S(NGAS,ICLS,NCLS,IACTIN,IACTOUT,IEXC)
*
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* Slave routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS(NGAS,NCLS)
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*
C     write(6,*) ' ICLS in EXCCLS '
C     call iwrtma(icls,ngas,ncls,ngas,ncls)

      IZERO = 0
      CALL ISETVC(IACTOUT,IZERO,NCLS)
      DO ICLSIN = 1, NCLS
        IF(IACTIN(ICLSIN).NE.0) THEN
        DO ICLSOUT = 1, NCLS
*. Number of anihilations and creations required to connect classes
          NANNI = 0
          NCREA = 0
          DO IGAS = 1, NGAS
*
            IDIF = ICLS(IGAS,ICLSOUT)-ICLS(IGAS,ICLSIN)
            IF(IDIF .GT. 0 ) THEN
              NCREA = NCREA + IDIF
            ELSE IF (IDIF .LT. 0 ) THEN
              NANNI = NANNI - IDIF
            END IF
          END DO
*
          IF(NCREA.LE.IEXC) THEN
            IACTOUT(ICLSOUT) = 1
          END IF
*
        END DO
        END IF
      END DO
*
      NACTOUT = 0
      DO ICLSOUT = 1, NCLS
        IF(IACTOUT(ICLSOUT).GT.0) NACTOUT = NACTOUT + 1
      END DO
*
      NTEST = 1
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Output from EXCCLS '
         WRITE(6,*) ' ==================='
         WRITE(6,*)
         WRITE(6,*) ' Number of active output classes ',NACTOUT
         WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Active output classes '
         CALL IWRTMA(IACTOUT,NCLS,1,NCLS,1)
      END IF
*
      RETURN
      END
