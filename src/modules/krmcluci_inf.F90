!-----------------------------------------
! Module with variables previously defined
! at krmcluci_inf.h include file (COMMON /CIVCINF/)
!-----------------------------------------
!  Information used for KRMC and LUCIAREL interface
!  ------------------------------------------------
!  1) CI energy
!  2) IRIQ for real/imaginary mismatch when running real groups
!     in LUCIAREL (preliminary implementation)
   module krmcluci_inf_mod
   implicit none
   real*8 :: ECI
   integer :: IRIQ
   end module krmcluci_inf_mod
