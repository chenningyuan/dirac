!-----------------------------------------
! Module with variables previously defined
! at dipole.h include file
!-----------------------------------------
   module dipole_mod
   use mxcent_mod, only: MXCOOR
   use kinds_mod, only: dp
   implicit none
!
   REAL(dp) DIPMN(3), DIPME(3), DIPMCORR(3), DDIPN(3, MXCOOR)
   REAL(dp) DDIPE(3, MXCOOR), DDIPS(3, MXCOOR), DDIPR(3, MXCOOR)
!
   end module dipole_mod