      SUBROUTINE IFILL(IX,IPUT,LENGTH)
#include "implicit.h"
      DIMENSION IX(LENGTH)

      IF(LENGTH .LE. 0) RETURN
      DO I = 1,LENGTH
        IX(I) = IPUT
      ENDDO

      RETURN
      END


      SUBROUTINE DFILL(DX,DPUT,LENGTH)
#include "implicit.h"
      DIMENSION DX(LENGTH)

      IF(LENGTH .LE. 0) RETURN
      DO I = 1,LENGTH
        DX(I) = DPUT
      ENDDO

      RETURN
      END


      SUBROUTINE SET_LBIT(LBIT_INTEGER,IDIGIT)
#include "implicit.h"
      LOGICAL LBIT

      IF(.NOT. LBIT(LBIT_INTEGER,IDIGIT)) THEN
        LBIT_INTEGER = LBIT_INTEGER + 2**(IDIGIT - 1)
      ENDIF

      RETURN
      END
