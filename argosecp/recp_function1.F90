MODULE RECP_FUNCTION
CONTAINS

SUBROUTINE RECP_FUNCTION_FACAB (binom,An,Bn,crda,crdb,xab)
  IMPLICIT NONE
! variables : global
  INTEGER :: An, Bn
! dimension binom(*), crda(*), crdb(*), xab(*)
  REAL(8) :: binom
  REAL(8) :: crda
  REAL(8) :: crdb
  REAL(8) :: xab
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

END MODULE RECP_FUNCTION
