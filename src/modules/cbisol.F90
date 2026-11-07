!-----------------------------------------
! Module with variables previously defined
! at cbisol.h include file (COMMON /CBISOL/)
!-----------------------------------------
   module cbisol_mod
   implicit none
   real*8 :: ORICAV(3), RCAV(3), EPDIEL
   integer :: LCAVMX, LMTOT, LMNTOT, NCNTCV
   logical :: SOLVNT
   end module cbisol_mod
