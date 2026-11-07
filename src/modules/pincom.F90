!-----------------------------------------
! Module with variables previously defined
! at pincom.h include file (COMMON /PINCOM/)
!-----------------------------------------
   module pincom_mod
   use maxorb_mod, only: MXCORB
   implicit none
   integer :: IPIND(MXCORB), INDGEN(MXCORB)
   end module pincom_mod
