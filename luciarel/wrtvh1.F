      SUBROUTINE WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
*
* Write one-electron integrals with symmetry IVSM
* ISYM = 1 => Only lower triangular matrix included
*
* Jeppe Olsen, Jan. 1999
* Timo Fleig, June 2001: relativistic case
*
#include "implicit.inc"
#include "mxpdim.inc"
*. General input
#include "symm.inc"
      INTEGER NRPSM(NSMOB),NCPSM(NSMOB)
*. Specific input
      DIMENSION H(*)
*
      IOFF = 1
      DO ISM = 1, NSMOB
        JSM = IDBGMULT(INVELM(ISM),IHSM)
        NI = NRPSM(ISM)
        NJ = NCPSM(JSM)
        print*,'JSM,NJ',JSM,NJ
        IF(ISYM.EQ.0.OR.ISM.GT.JSM) THEN
*. Complete block
          WRITE(6,*) ' Block with symmetry ISM, JSM ',ISM,JSM
          CALL WRTMAT(H(IOFF),NI,NJ,NI,NJ)
          IOFF = IOFF + NI*NJ
        ELSE IF (ISYM.EQ.1.AND.ISM.EQ.JSM) THEN
          CALL PRSYM(H(IOFF),NI)
          IOFF = IOFF * NI*(NI+1)/2
        END IF
      END DO
*
      RETURN
      END
