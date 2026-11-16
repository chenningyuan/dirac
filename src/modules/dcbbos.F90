!-----------------------------------------
! Module with variables previously defined
! at dcbbos.h include file (COMMON /DCIBOS/)
!-----------------------------------------
!     Real or approximate boson symmetry of orbitals.
!     Depends on maxorb.h
   module dcbbos_mod
   use maxorb_mod, only: MXCORB
   implicit none
   integer :: IBOSYM(MXCORB)
   end module dcbbos_mod
