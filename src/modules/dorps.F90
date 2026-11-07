!-----------------------------------------
! Module with variables previously defined
! at dorps.h include file (COMMON /DORPS/)
!-----------------------------------------
   module dorps_mod
   use mxcent_mod, only: MXCENT
   implicit none
   integer :: NDCORD(2)
   logical :: DOREPS(0:7), DOCOOR(3,MXCENT), DCORD(MXCENT,3,2), DCORGD(0:MXCENT,3,2), &
              DOPERT(0:3*MXCENT,2)
   end module dorps_mod
