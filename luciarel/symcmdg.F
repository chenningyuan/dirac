      SUBROUTINE SYMCMDG(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1*I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* Double group version, written for compatibility with general symmetry
*
#include "mxpdim.inc"
#include "symm.inc"
*
      IF(ITASK.EQ.1) THEN
        I1  = IDBGMULT(INVELM(I2),I12)
      ELSE IF(ITASK.EQ.2) THEN
        I2  = IDBGMULT(INVELM(I1),I12)
      ELSE IF (ITASK.EQ.3) THEN
        I12 = IDBGMULT(I1,I2)
      END IF
*
      RETURN
      END
