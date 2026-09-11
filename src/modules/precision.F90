!-----------------------------------------
! Module to define kinds of variables
!-----------------------------------------
   module kinds_mod
   use iso_fortran_env, only: int8, int16, int32, int64, real32, real64
   implicit none
    ! Integers
   integer, parameter :: i8  = int8
   integer, parameter :: i16 = int16
   integer, parameter :: i32 = int32
   integer, parameter :: i64 = int64
!
    ! Reals
   integer, parameter :: sp = real32   ! "single precision"
   integer, parameter :: dp = real64   ! "double precision"
   end module kinds_mod