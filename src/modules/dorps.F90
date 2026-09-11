!-----------------------------------------
! Module with variables previously defined
! at dorps.h include file
!-----------------------------------------
   module dorps_mod
!
   use mxcent_mod, only: MXCENT
   implicit none
   save
!
   logical :: DOREPS(0:7)
   logical :: DOCOOR(3, MXCENT)
   logical :: DCORD(MXCENT, 3, 2)
   logical :: DCORGD(0:MXCENT, 3, 2)
   logical :: DOPERT(0:3*MXCENT, 2)
!
   integer :: NDCORD(2)
!
 end module dorps_mod
