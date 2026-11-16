!-----------------------------------------
! Module with variables previously defined
! at dcb_ashift.h include file (COMMON /DCC_ASHIFT/, /DCR_AHIFT/, /DCI_ASHIFT/)
!-----------------------------------------
!** Atomic projectors ***
   module dcb_ashift_mod
   implicit none
   integer, parameter :: MAXAREF = 10
   integer, parameter :: MAXAVAL = 25
   real*8 :: ASHIFT_PAR(MAXAVAL,MAXAREF)
   integer :: ASHIFT_NVAL(MAXAREF)
   character(len=6) :: AREFFIL(MAXAREF)
   character(len=72) :: AVECREF(MAXAVAL,MAXAREF)
   end module dcb_ashift_mod
