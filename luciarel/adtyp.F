      SUBROUTINE ADTYP(NATP,JTYP,NGAS,LOCC,ROCC)
*
* Obtain type of J so
* <L! a+J !R> is nonvanishing
*
      INTEGER LOCC(NGAS),ROCC(NGAS)
*. First ensure that the types differ by creation of one electron
      ICREA = 0
      IANNI = 0
      DO IGAS = 1,NGAS,1
        IF(LOCC(IGAS).LT.ROCC(IGAS)) THEN
         IANNI = IANNI + ROCC(IGAS)-LOCC(IGAS)
        ELSE IF (LOCC(IGAS).GT.ROCC(IGAS) ) THEN
         ICREA = ICREA + LOCC(IGAS) - ROCC(IGAS)
        END IF
      END DO
*
      IF(ICREA.EQ.1 .AND. IANNI.EQ.0) THEN
       NATP = 1
      ELSE
       NATP = 0
      END IF
*
      IF(NATP.EQ.1) THEN
        DO IGAS = 1,NGAS,1
          IF (LOCC(IGAS)-ROCC(IGAS) .EQ. 1) THEN
            IF(ICREA.EQ.1) THEN
             JTYP = IGAS
             ICREA = 1
            end if
          END IF
        END DO
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Single Creation connecting  '
        WRITE(6,*) ' Left string '
        CALL IWRTMA(LOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Right string '
        CALL IWRTMA(ROCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        IF(NATP.EQ.0) THEN
          WRITE(6,*) ' No connecting single creation'
        ELSE IF(NATP.EQ.1) THEN
          WRITE(6,*) ' Type of connecting creation op', JTYP
        END IF
      END IF
*
      RETURN
      END
