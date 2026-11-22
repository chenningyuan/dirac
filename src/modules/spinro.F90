!-----------------------------------------
! Module with variables previously defined
! at spinro.h include file (COMMON /SPINRO/)
!-----------------------------------------
   module spinro_mod
   use mxcent_mod, only: MXCENT
   implicit none
   real*8 :: SPRNUC(3,3,MXCENT), SPRNUC2(3,3,MXCENT), SPRNUC3(3,3,MXCENT), &
             SPRDNL(3,3,MXCENT), GTRANT(3,3,MXCENT)
   logical :: LINDET
   end module spinro_mod
