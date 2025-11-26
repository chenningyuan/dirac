!-----------------------------------------
! Module with variables previously defined
! at mxcent.h include file
!-----------------------------------------
   module mxcent_mod
   implicit none

   integer, parameter :: MXATOM   = 20
   integer, parameter :: MXCENT   = 300
   integer, parameter :: MXCOOR   = 3 * MXCENT
   integer, parameter :: MXBSETS  = 5
   integer, parameter :: MAX_NONT = 150

   end module mxcent_mod
