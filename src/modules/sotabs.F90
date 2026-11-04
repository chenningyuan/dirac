!-----------------------------------------
! Module with variables previously defined
! at sotabs.h include file (COMMON /SPOTAB/)
!-----------------------------------------
   module sotabs_mod
   use maxaqn_mod, only: MXAQN
   use maxorb_mod, only: MXSHEL
   implicit none
   integer :: ISHLDF(3,2,MXAQN,MXSHEL), ICMPDF(3,2,MXAQN,MXSHEL), KMAX1
   end module sotabs_mod
