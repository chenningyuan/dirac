!-----------------------------------------
! Module with variables previously defined
! at dummy.h include file
!-----------------------------------------
   module dummy_mod
   implicit none
!  make DUMMY write protected (on most computers)
   real*8,  parameter :: DUMMY = 1.0D20
   integer, parameter :: IDUMMY = -9999999
   integer, parameter :: IDUMMYARR(1:1) = -9999999
   end module dummy_mod
