      FUNCTION LEN_TT_BLOCK(ISM,NIA,NIB,NSMST)
*
* Length of TT block with total sym ISM
*
* Jeppe Olsen, May 99 in Aarhus
* Modified for relativistic case, Timo Fleig, June 2001
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Specific input
      INTEGER NIA(*), NIB(*)
*. General input
#include "symm.inc"
*
      LEN = 0
      DO IASM = 1, NSMST
        IBSM = IDBGMULT(INVELM(IASM),ISM)
        LEN = LEN + NIA(IASM)*NIB(IBSM)
      END DO
*
      LEN_TT_BLOCK = LEN
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' LEN_TT_BLOCK = ', LEN_TT_BLOCK 
      END IF
*
      RETURN
      END
