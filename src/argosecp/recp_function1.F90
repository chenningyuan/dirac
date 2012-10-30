MODULE RECP_FUNCTION1
CONTAINS

SUBROUTINE RECP_FUNCTION_FACAB (binom,An,Bn,crda,crdb,xab)
  IMPLICIT NONE
! variables : global
  INTEGER :: An, Bn
! dimension binom(*), crda(*), crdb(*), xab(*)
  REAL(8) :: binom,crda,crdb,xab
  DIMENSION binom(*), crda(*), crdb(*), xab(*)
! variables : local
  INTEGER :: I,Ai,Bi,SUM_A,SUM_B

  do I = 1, (An+Bn-1) 
     xab(I) = 0.d0
  enddo
 
  SUM_A=(An*(An-1))/2  ! sum from 1 to An
  SUM_B=(Bn*(Bn-1))/2
  do Ai=1,An
    do Bi=1,Bn
      xab((Ai-1)+Bi) = xab((Ai-1)+Bi)        &
           + binom(SUM_A+Ai)*crda((An+1)-Ai) &
           * binom(SUM_B+Bi)*crdb((Bn+1)-Bi)
    enddo
  enddo
END SUBROUTINE RECP_FUNCTION_FACAB

SUBROUTINE RECP_SETZERO_I1(A,D1)
  IMPLICIT NONE
  INTEGER :: A(*)
  INTEGER :: I,D1
  DO I = 1,D1
     A(I) = 0
  ENDDO
END SUBROUTINE RECP_SETZERO_I1

SUBROUTINE RECP_SETZERO_I2(A,D1,D2)
  IMPLICIT NONE
  INTEGER :: I,J,D1,D2,A(D1,D2)
  DO I = 1,D1
     DO J = 1,D2
        A(I,J) = 0
     ENDDO
  ENDDO
END SUBROUTINE RECP_SETZERO_I2


SUBROUTINE RECP_SETZERO_I3(A,D1,D2,D3)
  IMPLICIT NONE
  INTEGER :: I,J,K,D1,D2,D3,A(D1,D2,D3)
  DO I = 1,D1
     DO J = 1,D2
        DO K = 1,D3
           A(I,J,K) = 0
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE RECP_SETZERO_I3


SUBROUTINE RECP_SETZERO_R1(A,D1)
  IMPLICIT NONE
  REAL(8) :: A(:)
  INTEGER :: I,D1
  DO I = 1,D1
     A(I) = 0.0d0
  ENDDO
END SUBROUTINE RECP_SETZERO_R1

END MODULE RECP_FUNCTION1
