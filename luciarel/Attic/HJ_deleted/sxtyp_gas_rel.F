      SUBROUTINE SXTYP_GAS_REL(NSXTYP,ITP,JTP,NGAS,MXPNGAS,ILTP,IRTP)
*
* Two supergroups are given. Find single excitations that connect
* these supergroups
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ILTP(NGAS),IRTP(NGAS)
*. Output
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
*. Differences between occupations :
      NCREA = 0
      NANNI = 0
      DO IAS = 1, NGAS
        IF(ILTP(IAS).GT.IRTP(IAS)) THEN
         NCREA = NCREA + ILTP(IAS) - IRTP(IAS)
         ICREA = IAS
        ELSE IF(IRTP(IAS).GT.ILTP(IAS)) THEN
         NANNI = NANNI + IRTP(IAS) - ILTP(IAS)
         IANNI = IAS
        END IF
      END DO
*
      IF(NCREA.GT.1) THEN
*. Sorry : No single excitation connects
        NSXTYP = 0
      ELSE IF(NCREA .EQ. 1 ) THEN
*. Supergroups differ by one single excitation.
        NSXTYP = 1
        ITP(1) = ICREA
        JTP(1) = IANNI
      ELSE IF (NCREA.EQ.0 ) THEN
*. Supergroups are identical, connects with all
*  diagonal excitations.
        NSXTYP = 0
        DO IAS = 1, NGAS
          IF(IRTP(IAS).NE.0) THEN
            NSXTYP = NSXTYP + 1
            ITP(NSXTYP) = IAS
            JTP(NSXTYP) = IAS
          END IF
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SXTYP_GAS : '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input left supergroup '
        CALL IWRTMA(ILTP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Input right supergroup '
        CALL IWRTMA(IRTP,1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of connecting single excitations ', NSXTYP
        IF(NSXTYP.NE.0) THEN
          WRITE(6,*) ' Connecting single excitations '
          DO ISX = 1, NSXTYP
            WRITE(6,*) ITP(ISX),JTP(ISX)
          END DO
        END IF
      END IF
*
      RETURN
      END
