!-----------------------------------------
! Module with variables previously defined
! at dcbres.h include file (COMMON /DCIRES/, /DCLRES/, /DCRRES/)
!-----------------------------------------
!  DCBRES:  Parameters for resolution of open-shell manifold
! Updated by Miro Ilias, August 2016 for DFT-COSCI
   module dcbres_mod
   implicit none
   real*8 :: SCRRES, RESOUT
   integer :: IPRRES, ISTRES, INTRES
   logical :: NOVXC, COREDENS
   end module dcbres_mod
