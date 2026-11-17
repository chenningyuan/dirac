!-----------------------------------------
! Module with variables previously defined
! at ihm.inc include file (COMMON /CCIHI/, /CCIHR/, /CCIHL/)
!-----------------------------------------
!     Small common block to pass on information used in the
!     intermediate Hamiltonian formalism.
!     For the Intermediate Hamiltonian method we need to identify spinors as either
!     element of Pm or Pi space. We use an array that gives this info that has
!     the same structure as the EPS array (that is passed as argument).
   module ihm_mod
   implicit none
   integer, parameter :: MAXSPIM = 5000
   real*8 :: CIH, ER_IH(2,2), AIH, SHIFT_IH(2,6)
   integer :: IPIORB(MAXSPIM), IHSCHEME, NIH
   logical :: DOIH
   end module ihm_mod
