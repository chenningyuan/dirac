      FUNCTION I_GET_ADJSYM(ISYM)
*
* Obtain conjugate symmetry of ISYM
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "symm.inc"
#include "lucinp.inc"
*
      IF(PNTGRP.EQ.1) THEN
*. Usual pointgroup, symmetry and complex symmetry identical
        I_GET_ADJSYM = ISYM
      ELSE
*. Currently only working alternative is double group symmetry
        I_GET_ADJSYM = IADJSYM(ISYM)
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' I_GET_ADJSYM , ISYM, IADJSYM',
     &  ISYM, IADJSYM(ISYM)
      END IF
*
      RETURN
      END
