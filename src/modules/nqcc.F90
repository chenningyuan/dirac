!-----------------------------------------
! Module with variables previously defined
! at nqcc.h include file (COMMON /NQCC/)
!-----------------------------------------
   module nqcc_mod
   use mxcent_mod, only: MXCENT
   implicit none
   real*8 :: EFGEL(3,3,mxcent), EFGNUC(3,3,mxcent), EFGCORR(3,3,mxcent), &
             EFGTOT(3,3,mxcent), EFG_METHOD(3,3,mxcent)
   end module nqcc_mod
