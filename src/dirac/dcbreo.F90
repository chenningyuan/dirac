!-----------------------------------------
! Module with variables previously defined
! at dcbreo.h include file (COMMON /CBIREO/, /CBLREO/, /CBI2RE/, /CBL2RE/)
!-----------------------------------------
!     LMOORD: Reorder orbitals?
!     MXREORD: maximum number of reorderings
!     IREORD, IMOORD: orbital reordering
   module dcbreo_mod
   implicit none
   integer, parameter :: MXREORD = 250
   integer :: IMOORD(MXREORD,2), IREORD(2), IMOORF(MXREORD,2), IREORF(2)
   logical :: LMOORD, LMOORF
   end module dcbreo_mod
