!-----------------------------------------
! Module with variables previously defined
! at dcbibn.h include file
!-----------------------------------------
!  Labels for vector types:
!    JBCNDX for CI-type,
!    JBENDX for electron orb-type,
!    JBPNDX for positron orb-type
   module dcbibn_mod
   implicit none
   integer, parameter :: JBCNDX = -1
   integer, parameter :: JBENDX = 1
   integer, parameter :: JBPNDX = 2
   end module dcbibn_mod
