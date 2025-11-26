!-----------------------------------------
! Module with variables previously defined
! at nqcc.h include file
!-----------------------------------------
   module nqcc_mod
!
   use MXCENT_mod, only: MXCENT
   implicit none
!
   real(8) :: EFGEL     (3, 3, MXCENT)
   real(8) :: EFGNUC    (3, 3, MXCENT)
   real(8) :: EFGTOT    (3, 3, MXCENT)
   real(8) :: EFG_METHOD(3, 3, MXCENT)
!
   end module nqcc_mod
