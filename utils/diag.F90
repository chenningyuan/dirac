!************************************************************************ 
! Standalone program for comprehensive testing of all diagonalization
! routines used by DIRAC. 

! These are Householder, RSJACOBI, QJACOBI and
! DSYEVR.
!
!
! The program either reads the matrix file, or generates random matrix.
!
!
! Written by Miro Ilias, July 2014 Banska Bystrica, Slovakia
!   with the help of P.S. Bagus, T.Saue
!
!
!************************************************************************ 

module  prgdata

  integer :: N,NZ
  real*8, allocatable :: A(:,:)

end module prgdata


Program Test_DIRAC_Diagonalization_Routines
  use prgdata
 
  print *, 'hello!'
  

End Program Test_DIRAC_Diagonalization_Routines
