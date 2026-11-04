!-----------------------------------------
! Module with variables previously defined
! at cotabs.h include file (COMMON /COTABS/)
!-----------------------------------------
   module cotabs_mod
   use maxaqn_mod, only: MXAQN
   use maxorb_mod, only: MXCORB, MXSHEL
   implicit none
   integer :: ISHLCO(MXCORB), ICMPCO(MXCORB), ICSHCM(0:MXSHEL,0:MXAQN), ICOMID, ICOMAX
   end module cotabs_mod
