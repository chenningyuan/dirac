!-----------------------------------------
! Module with variables previously defined
! at freeze.inc include file (COMMON /RELCCSD_FREEZE/)
!-----------------------------------------
   module freeze_mod
   implicit none
   integer, parameter :: MAXFROZ = 100
   integer, parameter :: MAXIREP = 32
   integer :: NFROZEN, IFROZ(MAXFROZ,MAXIREP)
   logical :: IFROZ_SET_IN_INPUT
   end module freeze_mod
