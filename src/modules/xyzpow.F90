!-----------------------------------------
! Module with variables previously defined
! at xyzpow.h include file (COMMON /XYZPOW/)
!-----------------------------------------
   module xyzpow_mod
   use maxmom_mod, only: MXAQNM
   implicit none
   integer :: ISTEP(MXAQNM), MVAL(MXAQNM), NVAL(MXAQNM)
   end module xyzpow_mod
