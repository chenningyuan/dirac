      FUNCTION IFRMR_REL(WORK,IROFF,IELMNT)
*
* An integer array is stored in real array WORK,
* starting from WORK(IROFF). Obtain element
* IELMNT of this array
*
      INTEGER WORK(*)
*
#include "irat.inc"
C     PARAMETER(IRAT = 2 )
*. offset when work is integer array
      IIOFF = 1 + IRAT * (IROFF-1)
      IFRMR_REL = WORK(IIOFF-1+IELMNT)
*
      RETURN
      END
