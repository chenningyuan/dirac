      SUBROUTINE OP_T_OCC_REL(IOPOCC,ITOCC,IOPTOCC,IMZERO)
*
* The occupation of  general operator OP and an CC cluster operator 
* are given. Find occupation of resulting operator after all 
* contractions have been performed. All deexcitation operators
* (anni of particles and creation of holes) are contracted. If 
*  this is not possible, Imzero = 1 is returned
*
* Jeppe Olsen, still May 2000
*
*. The only allowed open shell reference allowed is 
*  the high spin state with all singly occupied electrons having alpha spin
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*. Input
      INTEGER IOPOCC(NGAS,4),ITOCC(NGAS,4)
*. Output
      INTEGER IOPTOCC(NGAS,4)
*
      CALL ICOPVE(ITOCC,IOPTOCC,4*NGAS)
*. Loop over parts of operator
      IMZERO = 0
      DO IGAS = 1, NGAS
       DO IAB = 1, 2
        IF(IHPVGAS_AB(IGAS,IAB).EQ.1)THEN 
*. Contraction of hole creations in operator with annihilations in T
          IOPTOCC(IGAS,2+IAB) = IOPTOCC(IGAS,2+IAB)-IOPOCC(IGAS,IAB)
          IF( IOPTOCC(IGAS,2+IAB) .LT. 0 ) IMZERO = 1
*. Addition of hole annihilation operators in Op and T
          IOPTOCC(IGAS,2+IAB) = IOPTOCC(IGAS,2+IAB)+IOPOCC(IGAS,2+IAB)
        ELSE IF (IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
*. Contraction of particle annihilations in operator with creations in T
          IOPTOCC(IGAS,IAB) = IOPTOCC(IGAS,IAB) - IOPOCC(IGAS,2+IAB)
          IF( IOPTOCC(IGAS,IAB) .LT. 0 ) IMZERO = 1
*. Addition of particle creation operators in Op and T
          IOPTOCC(IGAS,IAB) = IOPTOCC(IGAS,IAB) + IOPOCC(IGAS,IAB)
        END IF
       END DO 
      END DO
*. Ensure that number of creations/annihilations do not exceed 
*. number of orbitals
      DO ICAAB = 1, 4
        DO IGAS = 1, NGAS
          IF(IOPTOCC(IGAS,ICAAB).GT.NOBPT(IGAS)) IMZERO = 1 
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Combination of general operator and CC operator'
        WRITE(6,*) ' Input general and CC operators '
        CALL WRT_SPOX_TP_CC_REL(IOPOCC,1)
        CALL WRT_SPOX_TP_CC_REL(ITOCC,1)
        WRITE(6,*) ' Output operator '
        CALL WRT_SPOX_TP_CC_REL(IOPTOCC,1)
        IF(IMZERO.EQ.1) THEN
          WRITE(6,*) ' Vanishing operator '
        END IF
      END IF
*
      RETURN
      END
