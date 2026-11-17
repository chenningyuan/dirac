!-----------------------------------------
! Module with variables previously defined
! at typesz_mpi2.h include file (COMMON /TYPE_SZMPI2/)
! Byte sizes of MPI-2 data types (set by MPIXCALC_TYPESZ); the header
! declared them only with VAR_MPI2, the module always does.
!-----------------------------------------
   module typesz_mpi2_mod
   implicit none
   integer(kind=8) :: ISIZE_DP, ISIZE_MY_INT, ISIZE_INTDEF, ISIZE_INT8, ISIZE_LOG
   end module typesz_mpi2_mod
