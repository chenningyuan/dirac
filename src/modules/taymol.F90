!-----------------------------------------
! Module with variables previously defined
! at taymol.h include file (COMMON /TAYMOL/)
!-----------------------------------------
   module taymol_mod
   use mxcent_mod, only: MXCOOR
   implicit none
   real*8 :: ERGMOL, GRDMOL(MXCOOR), HESMOL(MXCOOR,MXCOOR)
   end module taymol_mod
