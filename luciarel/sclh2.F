      SUBROUTINE SCLH2(XLAMBDA)
*
*. Scale two electron integrals
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      CALL SCALVE(WORK(KINT2),XLAMBDA,NINT2)
*
      RETURN
      END
