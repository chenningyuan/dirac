!-----------------------------------------
! Module with variables previously defined
! at maxaqn.h include file
!-----------------------------------------
   module maxaqn_mod
   implicit none
!
   integer, parameter :: MXQN  = 10
   integer, parameter :: MXAQN = MXQN * (MXQN + 1) / 2
!
   end module maxaqn_mod
