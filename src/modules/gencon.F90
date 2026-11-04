!-----------------------------------------
! Module with variables previously defined
! at gencon.h include file (COMMON /GENCON/)
!-----------------------------------------
!   - Parameter for maximum number of functions in a block:
!     MXCON_PARAM = NPRIM * l(l+1)/2 * NCONTR
!   - DIRCON true if contracted basis set is used
   module gencon_mod
   implicit none
   integer, parameter :: MXCON_PARAM = 300
   logical :: DIRCON(4)
   end module gencon_mod
