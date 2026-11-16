!-----------------------------------------
! Module with variables previously defined
! at dcblap.h include file (COMMON /DCBLAP/)
!-----------------------------------------
!  DCBLAP:  Parameters for Laplace module
   module dcblap_mod
   use lapdim_mod, only: MXLAP
   implicit none
   real*8 :: STEPMX, TOLRNG, TOLPAR, TOLLAP, XPNTS(MXLAP), WGHTS(MXLAP)
   integer :: NLAP, MXITER, IPRLAP
   logical :: INILAP, FNDLAP
   end module dcblap_mod
