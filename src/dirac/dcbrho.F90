!-----------------------------------------
! Module with variables previously defined
! at dcbrho.h include file (COMMON /CBRRHO/)
!-----------------------------------------
!*** 3-dimensional density plot ***
!     NCUBORB - number of g/u orbitals to make cubes of
!     ICUBORB - which orbitals to make cubes from
   module dcbrho_mod
   implicit none
   real*8 :: CUBADJ(2)
   integer :: ILSCUB, IPRRHO, NCUBE(0:3), ICUBORB(1000,2), NCUBORB(2)
   end module dcbrho_mod
