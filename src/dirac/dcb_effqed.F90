!-----------------------------------------
! Module with variables previously defined
! at dcb_effqed.h include file (COMMON /DCC_EFFQED/, /DCR_EFFQED/, /DCI_EFFQED/)
!-----------------------------------------
!** Effective QED potentials ***
   module dcb_effqed_mod
   implicit none
   integer, parameter :: MAXREF = 10
   integer, parameter :: MAXKAP = 7
   real*8 :: EFFQED_PAR(2,MAXKAP,MAXREF)
   integer :: EFFQED_NKAP(MAXREF), EFFQED_IPAR(MAXKAP,MAXREF)
   character(len=6) :: REFFIL(MAXREF)
   character(len=72) :: VECREF(MAXKAP,MAXREF)
   end module dcb_effqed_mod
