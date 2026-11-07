!-----------------------------------------
! Module with variables previously defined
! at aosotr.h include file (COMMON /AOSOTR/)
!-----------------------------------------
   module aosotr_mod
   use maxorb_mod, only: MXCORB
   implicit none
   real*8 :: CTRAN(MXCORB,8)
   integer :: ITRAN(MXCORB,8), IAOAO(MXCORB), IPTYP(MXCORB), IPCEN(MXCORB), JTRAN(MXCORB)
   end module aosotr_mod
