      SUBROUTINE KRA_SYM(ISYM,IBU,IKRSM,IBUK)
* A spinor symmetry ISYM is given for a spinor, Barred/Unbarred
* given by IBU
*
* Find corresponding symmetry of spinor with opposite Barred/Unbarred symmetry
* as well as opposite Kramers index IBUK
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "symm.inc"
*
      IPNTGP_SYM = ISPSF_MO(ISYM,IBU)
      IF(IBU.EQ.1) THEN
        IKRSM = IMOSF_SP(IPNTGP_SYM,2)
        IBUK = 2
      ELSE
        IKRSM = IMOSF_SP(IPNTGP_SYM,1)
        IBUK = 1
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Kramers connected symmetries',ISYM,IKRSM
      END IF
*
      RETURN
      END
