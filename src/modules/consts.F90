!-----------------------------------------
! Module with variables previously defined
! at consts.h include file
!-----------------------------------------
   module consts_mod
   use kinds_mod, only: dp
   implicit none
   real(dp), parameter :: DM1  = -1.0_dp
   real(dp), parameter :: DMP5 = -0.5_dp
   real(dp), parameter :: DP5  =  0.5_dp
   real(dp), parameter :: D0   =  0.0_dp
   real(dp), parameter :: D1   =  1.0_dp
   real(dp), parameter :: D2   =  2.0_dp
   real(dp), parameter :: D3   =  3.0_dp
   real(dp), parameter :: D4   =  4.0_dp
   real(dp), parameter :: D6   =  6.0_dp
   real(dp), parameter :: D7   =  7.0_dp
   real(dp), parameter :: D8   =  8.0_dp
   end module consts_mod
