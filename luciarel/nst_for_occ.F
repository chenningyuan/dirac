      FUNCTION NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
*
* A supergroup is defined by NELEC_PER_GAS. Find the 
* Total number of strings of this supergroup
*
* Jeppe Olsen, July 2000 (HNIE)
*
#include "implicit.inc"
*
      INTEGER NELEC_PER_GAS(*),NORB_PER_GAS(*)
*
      NSTR = 1
      DO IGAS = 1, NGAS
       IF(NELEC_PER_GAS(IGAS).NE.0) THEN
        NSTR = NSTR*IBION(NORB_PER_GAS(IGAS),NELEC_PER_GAS(IGAS))
       END IF
      END DO
*
      NST_FOR_OCC = NSTR
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of electrons and orbitals per gasspace '
        CALL IWRTMA(NELEC_PER_GAS,1,NGAS,1,NGAS)
        CALL IWRTMA(NORB_PER_GAS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of strings = ', NSTR
      END IF
*
      RETURN
      END
