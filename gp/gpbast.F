      SUBROUTINE IFILL(IX,INSIDE,LENGTH)
C
#include "implicit.h"
C
      INTEGER IX(LENGTH)
C
C
      IF(LENGTH.LE.0) RETURN
C
      DO I=1,LENGTH
        IX(I) = INSIDE
      ENDDO
C
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE DFILL(DX,DINSIDE,LENGTH)
C
#include "implicit.h"
C
      DIMENSION DX(LENGTH)
C
C
      IF(LENGTH.LE.0) RETURN
C
      DO I=1,LENGTH
        DX(I) = DINSIDE
      ENDDO
C
C
      RETURN
      END
