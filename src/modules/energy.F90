!-----------------------------------------
! Module with variables previously defined
! at energy.h include file (COMMON /ENERGY/)
!-----------------------------------------
! Prefix -
!  ENER : energy contribution
!  GRAD : gradient contribution
! Suffix -
!  KE : kinetic energy
!  NA : nuclear attraction
!  EE : electron-electron repulsion
!  NN : nucleus-nucleus repulsion
!  FS : reorthonormalization ( trace(Fmat Smat) )
   module energy_mod
   use mxcent_mod, only: MXCOOR
   implicit none
   real*8 :: ENERKE, ENERNA, ENEREE, ENERNN, GRADKE(MXCOOR), GRADNA(MXCOOR), &
             GRADEE(MXCOOR), GRADNN(MXCOOR), GRADFS(MXCOOR)
   end module energy_mod
