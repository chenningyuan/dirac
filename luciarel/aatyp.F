      SUBROUTINE AATYP(NAATP,ITYP,JTYP,LOCC,ROCC,NGAS)
*
* Obtain types of I,J so
* <L!aI aJ !R> is nonvanishing
* only combinations with type(I) .ge. type(J) * are included
*
*
      INTEGER LOCC(NGAS),ROCC(NGAS)
*. First ensure that the types differs by annihilation of two elecs
      ICREA = 0
      IANNI = 0
      DO IGAS = 1, NGAS
        IF(LOCC(IGAS).LT.ROCC(IGAS)) THEN
         IANNI = IANNI + ROCC(IGAS)-LOCC(IGAS)
        ELSE IF (LOCC(IGAS).GT.ROCC(IGAS) ) THEN
         ICREA = ICREA + LOCC(IGAS) - ROCC(IGAS)
        END IF
      END DO
*
      IF(ICREA.EQ.0 .AND. IANNI.EQ.2) THEN
       NAATP = 1
      ELSE
       NAATP = 0
      END IF
*
      IANNI = 0
      IF(NAATP.EQ.1) THEN
        DO IGAS = 1, NGAS
          IF(LOCC(IGAS)-ROCC(IGAS) .EQ. -2 ) THEN
            ITYP = IGAS
            JTYP = IGAS
          ELSE IF (LOCC(IGAS)-ROCC(IGAS) .EQ. -1 ) THEN
            IF(IANNI.EQ.1) THEN
             ITYP = IGAS
             IANNI = 2
            ELSE IF( IANNI.EQ. 0 ) THEN
             JTYP = IGAS
             IANNI = 1
            END IF
          END IF
        END DO
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Double Annihilation connecting  '
        WRITE(6,*) ' Left string '
        CALL IWRTMA(LOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Right string '
        CALL IWRTMA(ROCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        IF(NAATP.EQ.0) THEN
          WRITE(6,*) ' No connecting double annihilation'
        ELSE IF(NAATP.EQ.1) THEN
          WRITE(6,*) ' Types of connecting annihilation ops', ITYP,JTYP
        END IF
      END IF
*
      RETURN
      END
