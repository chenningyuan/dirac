!-----------------------------------------
! Module with variables previously defined
! at dcbtrp.h include file (COMMON /DCBTRP/)
!-----------------------------------------
!*****************************************************************************
!  DCBTRP:  Property opearators to be transformed to MO-basis
   module dcbtrp_mod
   implicit none
   integer, parameter :: MAXTRP = 350
   integer :: IPRTRP, NTRPP, LTRPP(MAXTRP)
   end module dcbtrp_mod
