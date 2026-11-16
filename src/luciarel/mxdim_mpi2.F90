!-----------------------------------------
! Module with variables previously defined
! at mxdim_mpi2.h include file
!-----------------------------------------
!     MPI-2 memory window dimensioning
!     ======================================
   module mxdim_mpi2_mod
   implicit none
   integer, parameter :: MAX_WINDOW_SIZE = 1900000000
   integer, parameter :: IDUMMY_BUFF_SZ = 1000000
   end module mxdim_mpi2_mod
