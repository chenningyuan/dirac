      SUBROUTINE IOFF_SYMBLK_MAT(NSMST,NA,NB,ITOTSM,IOFF,IRESTRICT)
*
* Offset for symmetry blocked matrix. 
* Row-symmetry used as primary index
*
* Jeppe Olsen, Summer of 99
*               July 2000 : IRESTRICT added
* Modified for relativistic case, Timo Fleig, July 2001
*
* IRESTRICT = 0 : No restrictions of symmetry blocks
*           = 1 : Only blocks with ASM .GE. BSM are included
*           =-1 : Only blocks with ASM .LT. BSM are included
#include "implicit.inc"
*. General input  
#include "mxpdim.inc"
#include "symm.inc"
*. Specific input 
      INTEGER NA(NSMST),NB(NSMST)
*. Output
      INTEGER IOFF(NSMST)
*
      IB = 1
      DO IASM = 1, NSMST
        IBSM = IDBGMULT(INVELM(IASM),ITOTSM)
        INCLUDED = 0
        IF(IRESTRICT.EQ.0) INCLUDED = 1
        IF(IRESTRICT.EQ.1.AND.IASM.GE.IBSM) INCLUDED = 1
        IF(IRESTRICT.EQ.-1.AND.IASM.LT.IBSM) INCLUDED = -1
        IF(INCLUDED.EQ.1) THEN
          IOFF(IASM) = IB
          IB = IB + NA(IASM)*NB(IBSM)
        ELSE
          IOFF(IASM) = 0
        END IF
      END DO
*
      RETURN
      END
