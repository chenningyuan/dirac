C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck FRAC2I */
      SUBROUTINE FRAC2I(A,IN,ID,IERR)
C*****************************************************************************
C
C     Transform a fraction (in the string A) to the numerator 
C     and denumerator
C     E.g. A = ' 3/4'  --> IN = 3, ID = 4
C     
C*****************************************************************************
#include <implicit.h>
      CHARACTER A*(*),B*80
C
      IERR = 0
      IDX = INDEX(A,'/')
      IF (IDX.EQ.0) THEN
         IERR = -1
         RETURN
      END IF
C
      WRITE(B,'(A)') A(1:IDX-1)
      READ (B,'(I)') IN
      I = LEN(A)
 20   CONTINUE
         IF (INDEX('0123456789',A(I:I)).EQ.0.AND.I.GT.IDX) THEN
            I = I - 1
            GOTO 20
         END IF
C
      WRITE(B,'(A)') A(IDX+1:I)
      READ(B,'(I)') ID
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SUMMMH */
      SUBROUTINE SUMMMH(AMAT,NDIM,NZ,LRA,LCA)
C*****************************************************************************
C
C     Take sum of a general quaternionic matrix AMAT and
C     its Hermitian conjugate; return in AMAT. Note that this gives
C     an Hermitian result !
C
C     Written by J. Thyssen, Dec 1, 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0, D2 = 2.00D00)
      DIMENSION AMAT(LRA,LCA,NZ)
C
C     Take difference
C
      DO J = 1,NDIM
        DO I = 1,(J-1)
          AMAT(I,J,1) = AMAT(I,J,1)+AMAT(J,I,1)
          AMAT(J,I,1) = AMAT(I,J,1)
        ENDDO
        AMAT(J,J,1) = D2*AMAT(J,J,1)
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(I,J,IZ) = AMAT(I,J,IZ)-AMAT(J,I,IZ)
            AMAT(J,I,IZ) = -AMAT(I,J,IZ)
          ENDDO
          AMAT(J,J,IZ) = D0
        ENDDO
      ENDDO
C
      RETURN
      END
