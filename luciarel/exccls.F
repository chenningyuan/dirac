      SUBROUTINE EXCCLS(NCLS,IACTIN,IACTOUT,IEXC)
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* Master routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*. From the common blocks
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "cgas.inc"
*
      CALL EXCCLS_S(NGAS,WORK(KLOCCLS),NCLS,IACTIN,IACTOUT,IEXC)
*
      RETURN
      END
