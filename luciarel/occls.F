      SUBROUTINE OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and
* IOCCLS = allowed distributions of electrons
*
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS)
*. Output
      DIMENSION  IOCCLS(NGAS,*)
*. Local scratch
#include "mxpdim.inc"
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
*
      ISKIP = 1
      NOCCLS = 0
*. start with smallest allowed number
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: '
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space
        NEGA=0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: '
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : '
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
            CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes '
           CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      RETURN
      END
