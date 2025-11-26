!-----------------------------------------
! Module with variables previously defined
! at fitprimit.h include file
!-----------------------------------------
   module fitprimit_mod
!
   use primitives_mod, only: MXPRIM_FIT
   use aovec_mod, only: MXCONT
   implicit none
!
   real*8 :: PRIEXP_FIT(MXPRIM_FIT)
   real*8 :: PRICCF_FIT(MXPRIM_FIT, MXCONT)
   real*8 :: PRICRX_FIT(MXPRIM_FIT)
   real*8 :: PRICRY_FIT(MXPRIM_FIT)
   real*8 :: PRICRZ_FIT(MXPRIM_FIT)
!
   integer :: NPSHEL_FIT
   integer :: NPLSH_FIT
   integer :: NPSSH_FIT
!
   end module fitprimit_mod

