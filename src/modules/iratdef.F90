!-----------------------------------------
! Module with variables previously defined
! at iratdef.h include file
!-----------------------------------------
!     IRAT  = (real word length) / (integer word length)
!     IRAT2 = (real word length) / (half-integer word length)
!             if available and used, otherwise IRAT2 = IRAT
!     LRAT  = (real word length) / (logical word length)
   module iratdef_mod
   implicit none
#if defined (INT_STAR8)
   integer, parameter :: IRAT = 1, LRAT = 1
#else
   integer, parameter :: IRAT = 2, LRAT = 2
#endif
#if defined (VAR_STAR2)
   integer, parameter :: IRAT2 = 4
#else
   integer, parameter :: IRAT2 = IRAT
#endif
   end module iratdef_mod
