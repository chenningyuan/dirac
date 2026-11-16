!-----------------------------------------
! Module with variables previously defined
! at dcbrot.h include file (COMMON /CBLMIX/, /CBRMIX/, /CBCMIX/)
!-----------------------------------------
!** Modified virtual orbitals ***
   module dcbrot_mod
   implicit none
   real*8 :: ROTANG
   logical :: ROTORB
   character(len=72) :: VECROT(2)
   end module dcbrot_mod
