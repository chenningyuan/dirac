!-----------------------------------------
! Module with variables previously defined
! at para.h include file
!-----------------------------------------
   module para_mod
   implicit none
   integer, parameter :: LMAX = 4
   integer, parameter :: LMAX_OCC = 3
   integer, parameter :: MXPRIML = 38
   integer, parameter :: MXCONTL = 38
   integer, parameter :: MXCART = 520
   integer, parameter :: NDFMX = 4*Lmax+4
   end module para_mod
