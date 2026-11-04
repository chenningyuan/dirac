!-----------------------------------------
! Module with variables previously defined
! at dirprt.h include file (COMMON /DIRPRT/)
!-----------------------------------------
   module dirprt_mod
   implicit none
   real*8 :: SIGNDR(4,3)
   integer :: NOATMS, IATOMS(4), ISOPDR(4), NCNTDR(4), ICNTDR(4)
   logical :: TRANS
   end module dirprt_mod
