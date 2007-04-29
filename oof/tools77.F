C
C
C   Routines and functions which could be only in Fortran77 format
C
C

      SUBROUTINE PRINT_REAL_VECTOR(N,VEC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Called from "vect_mod.F90"/print_vector FOTRAN90 routine
C
C Written by Miro Ilias, April 2007, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION VEC(N)
      CALL QENTER('PRINT_REAL_VECTOR')

      DO I = 1, N
        WRITE(LUPRI,*) I,'...',VEC(I)
      ENDDO

      CALL QEXIT('PRINT_REAL_VECTOR')
      RETURN
      END

      SUBROUTINE PRINT_INT_VECTOR(N,IVEC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Print the integer vector. 
C
C Called from "vect_mod.F90"/print_vector FOTRAN90 routine. The F77 routine
C is due to data typa change from REAL*8 to INTEGER, which is NOT premitted
C in standard FOTRAN90 language.
C
C Written by Miro Ilias, April 2007, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION IVEC(N)
      CALL QENTER('PRINT_INT_VECTOR')
      DO I = 1, N
        WRITE(LUPRI,'(1X,I3,A,I8,A,$)') I,'->',IVEC(I),'|'
        IF (MOD(I,5).EQ.0) WRITE(LUPRI,*)
      ENDDO
      CALL QEXIT('PRINT_INT_VECTOR')
      RETURN
      END

