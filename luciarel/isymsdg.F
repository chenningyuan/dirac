      FUNCTION ISYMSDG(STRING,NEL)
*
* Symmmetry of string, Double group  version
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "symm.inc"
*

*. Specific input
      INTEGER STRING(*)
*. Note : How to switch between up and down spinors
      ISYM = 1
      DO IEL = 1, NEL
        ISYM = IDBGMULT(ISYM,ISMFTO(STRING(IEL)))
      END DO
*
      ISYMSDG = ISYM
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
*
      RETURN
      END
