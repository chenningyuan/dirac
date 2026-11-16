!-----------------------------------------
! Module with variables previously defined
! at SysDef.inc include file
! (the header's Arch_Type, a DATA-initialized local, was never used)
!-----------------------------------------
   module SysDef_mod
   implicit none
#if defined(INT_STAR8)
   integer, parameter :: ITOB = 8, RTOB = 8, RTOI = 1
#else
   integer, parameter :: ITOB = 4, RTOB = 8, RTOI = 2
#endif
   end module SysDef_mod
