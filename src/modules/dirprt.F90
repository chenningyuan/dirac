!-----------------------------------------
!
! Module with variables previously defined
!
! at dirprt.h include file
!
!-----------------------------------------

module dirprt_mod

   implicit none

   real*8 :: SIGNDR(4,3)

   logical :: TRANS

   integer :: NOATMS
   integer :: IATOMS(4)
   integer :: ISOPDR(4)
   integer :: NCNTDR(4)
   integer :: ICNTDR(4)

end module dirprt_mod
