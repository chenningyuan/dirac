      FUNCTION IOFF_DBG_LIST(ILIST,ISYM,JSYM,KSYM,LSYM)
*
* Offset to given integral block of double group integrals
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "lucinp.inc"
*
      IOFF_DBG_LIST =
     &IOFF_DBG_LIST_S(WORK(KPINT2D(ILIST)),ISYM,JSYM,KSYM,LSYM,
     &                NSMOB)
*
      RETURN
      END
