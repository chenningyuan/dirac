!-----------------------------------------
! Module with variables previously defined
! at pincom.h include file
!-----------------------------------------
   module pincom_mod
!
   use primitives_mod, only: MXCORB
   implicit none
   private
!
   integer :: IPIND(MXCORB)
   integer :: INDGEN(MXCORB)
!
   public :: IPIND, INDGEN
!
   end module pincom_mod
