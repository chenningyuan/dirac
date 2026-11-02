!-----------------------------------------
! Module with variables previously defined
! at aovec.h include file
!-----------------------------------------
!
! MXCONT = maximum number of contracted functions in an AO input block
!          (in "MOLECULE.INP" input file)
! MXAOVC = used for MXCONT in a few places
!
! IF you change MXCONT you should make sure
! MAXPRD in cbirea.h is at least the same value and
! then rebuild the program using the command "make".
!
   module aovec_mod

   implicit none

   integer, parameter :: MXCONT = 40
   integer, parameter :: MXAOVC = MXCONT

   end module aovec_mod
