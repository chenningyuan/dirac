!-----------------------------------------
! Module with variables previously defined
! at maxmom.h include file
!-----------------------------------------
   module maxmom_mod
   implicit none
   integer, parameter :: MXQNM  = 21
   integer, parameter :: MXAQNM = MXQNM*(MXQNM + 1)/2
   end module maxmom_mod
