!-----------------------------------------
! Module with variables previously defined
! at maxorb.h include file
!-----------------------------------------
!
! hjaaj July 2000:
! MXSHEL = max number of symmetry independent contracted shells
! MXPRIM = max number of symmetry independent primitive shells
!  (e.g. a p shell contains 3 components)
! MXCORB = max total number of contracted functions = max NBASIS
!
   module maxorb_mod

   implicit none

   integer, parameter :: MXSHEL = 15000
   integer, parameter :: MXPRIM = 15000
   integer, parameter :: MXCORB = 40000

!  hjaaj Jan 2019: fitting is not implemented in Dirac, so
!  to save static memory I introduced the following three
!  variables for the fitting common blocks
   integer, parameter :: MXSHEL_FIT = 1
   integer, parameter :: MXPRIM_FIT = 1
   integer, parameter :: MXCORB_FIT = 1

   end module maxorb_mod
