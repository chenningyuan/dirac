!-----------------------------------------
! Module with variables previously defined
! at dipole.h include file (COMMON /DIPOLE/)
!-----------------------------------------
   module dipole_mod
   use mxcent_mod, only: MXCOOR
   implicit none
   real*8 :: DIPMN(3), DIPME(3), DIPMCORR(3), DDIPN(3,MXCOOR), DDIPE(3,MXCOOR), &
             DDIPS(3,MXCOOR), DDIPR(3,MXCOOR)
   end module dipole_mod
