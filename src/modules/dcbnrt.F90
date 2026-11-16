!-----------------------------------------
! Module with variables previously defined
! at dcbnrt.h include file (COMMON /DCINRT/, /DCLNRT/)
!-----------------------------------------
! For freezing ("NO ROTATION") orbitals in MCSCF
!     Dependencies: maxorb.h
   module dcbnrt_mod
   use maxorb_mod, only: MXCORB
   implicit none
   integer :: NOROT(MXCORB), NOROTC(MXCORB)
   logical :: LNOROT
   end module dcbnrt_mod
