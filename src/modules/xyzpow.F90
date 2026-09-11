!-----------------------------------------
! Module with variables previously defined
! at xyzpow.h include file
!-----------------------------------------
   module xyzpow_mod
!
   use maxmom_mod, only: MXAQNM
   implicit none
   save
!
   integer :: ISTEP(MXAQNM)
   integer :: MVAL(MXAQNM)
   integer :: NVAL(MXAQNM)
!
   end module xyzpow_mod
