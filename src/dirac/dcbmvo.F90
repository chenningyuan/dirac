!-----------------------------------------
! Module with variables previously defined
! at dcbmvo.h include file (COMMON /CBLMVO/, /CBIMVO/, /CBCMVO/, /CBDMVO/)
!-----------------------------------------
!** Modified virtual orbitals ***
   module dcbmvo_mod
   implicit none
   real*8 :: MVOWEIGHT
   integer :: IPRMVO, INTMVO, IF2MVO
   logical :: IONMVO, ADDREP
   character(len=72) :: VECMVO(2)
   end module dcbmvo_mod
