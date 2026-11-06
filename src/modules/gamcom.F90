!-----------------------------------------
! Module with variables previously defined
! at gamcom.h include file (COMMON /GAMCOM/)
!-----------------------------------------
   module gamcom_mod
   use maxaqn_mod, only: MXQN
   implicit none
   integer, parameter :: MAXJ = 4*(MXQN - 1) + 2
   integer :: JMAX0
   real(kind=8) :: WVAL
   real(kind=8) :: FJW(0:MAXJ), TABFJW(121*(MAXJ+7))
   end module gamcom_mod
