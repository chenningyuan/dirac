      SUBROUTINE WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
*
* Write occupation in each GASspace for general operator
*
* Jeppe Olsen, August 1999
*
#include "implicit.inc"
*
      INTEGER ICA(NGAS),ICB(NGAS),IAA(NGAS),IAB(NGAS)
*
      WRITE(6,*) ' Occupation of alpha-creation string '
      CALL IWRTMA(ICA,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of beta -creation string '
      CALL IWRTMA(ICB,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of alpha-annihilation string '
      CALL IWRTMA(IAA,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of beta -annihilation string '
      CALL IWRTMA(IAB,1,NGAS,1,NGAS)
*
      RETURN
      END
