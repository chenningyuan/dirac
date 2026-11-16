!-----------------------------------------
! Module with variables previously defined
! at dcbvec.h include file (COMMON /CBIVEC/, /CBCVEC/)
!-----------------------------------------
!*** Vector print ***
   module dcbvec_mod
   implicit none
   integer :: IPRCMP, ILABDF, IPRVEC
   character(len=72) :: VECPRI(2)
   end module dcbvec_mod
