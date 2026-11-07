!-----------------------------------------
! Module with variables previously defined
! at cbihr2.h include file (COMMON /CBIHR2/, /CBRHR2/)
!-----------------------------------------
   module cbihr2_mod
   implicit none
   real*8 :: THRTWO, THRFAC(2), SCRFCK
   integer :: IPRTWO, IPRNTA, IPRNTB, IPRNTC, IPRNTD, ICEDIF, IFTHRS
   logical :: RUNTWO, RTNTWO, TKTIME, SOFOCK, USRSCR, NCLHER
   end module cbihr2_mod
