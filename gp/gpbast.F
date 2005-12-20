      SUBROUTINE MATRIXIO(QMATRIX,NROW,NCOL,NZ,DOWRITE,FILENAME)
C
#include "implicit.h"
C
      LOGICAL TOBE,DOWRITE
      DIMENSION QMATRIX(NROW,NCOL,NZ)
      CHARACTER FILENAME*(*)
C
C
      IUNIT = 66
      INQUIRE(FILE = FILENAME,EXIST = TOBE)
C
      IF(TOBE.AND.DOWRITE) THEN
C       GOTO 1
        OPEN(IUNIT,FILE = FILENAME,STATUS = 'OLD',FORM = 'FORMATTED')
        CLOSE(IUNIT,STATUS = 'DELETE')
      ENDIF
C
      IF(.NOT.(TOBE).AND..NOT.(DOWRITE)) THEN
C       panic message
        GOTO 1
      ENDIF
C
      IF(.NOT.(DOWRITE)) THEN
        CALL DZERO(QMATRIX,NROW*NCOL*NZ)
      ENDIF
C
      OPEN(IUNIT,FILE = FILENAME,STATUS = 'UNKNOWN',FORM = 'FORMATTED')
C
      DO IZ = 1,NZ
        DO IROW = 1,NROW
          DO ICOL = 1,NCOL
            IF(DOWRITE) THEN
              WRITE(IUNIT,'(D15.6)') QMATRIX(IROW,ICOL,IZ)
            ELSE
              READ(IUNIT,'(D15.6)') QMATRIX(IROW,ICOL,IZ)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      CLOSE(IUNIT,STATUS = 'KEEP')
 1    CONTINUE
C
C
      RETURN
      END
C
C
C
C
C.......................................................................
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
C.......................................................................
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
