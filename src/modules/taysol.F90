!-----------------------------------------
! Module with variables previously defined
! at taysol.h include file (COMMON /TAYSOL/)
!-----------------------------------------
   module taysol_mod
   use mxcent_mod, only: MXCOOR
   implicit none
   real*8 :: ESOLTT, GSOLTT(MXCOOR), GSOLNN(MXCOOR), ESOLNN
   end module taysol_mod
