!-----------------------------------------
! Module with variables previously defined
! at chrsgn.h include file
!-----------------------------------------
   module chrsgn_mod
   implicit none
   character(len=1), parameter :: CHRSGN(-1:1) = ['-', ' ', '+']
   end module chrsgn_mod
