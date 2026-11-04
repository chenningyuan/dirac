!-----------------------------------------
! Module with variables previously defined
! at oneadr.h include file (COMMON /ONEADR/)
!-----------------------------------------
   module oneadr_mod
   use mxcent_mod, only: MXCENT
   implicit none
   integer :: LASTAD(2,3*MXCENT,0:7), MXABUF, MXABF2, LABUFD, LABUFI, LABUF1
   end module oneadr_mod
