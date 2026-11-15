!-----------------------------------------
! Module with variables previously defined
! at infpar.h include file (COMMON /DIRACINFPAR/)
! MPARID and NODES are also read from C (src/cfun/general.c), hence bind(C).
!-----------------------------------------
!  specific information for parallel calculations
!     .. mpi variables:
!        NUMNOD = number of slaves = number of nodes - 1 (from MPI_Comm_rank)
!        MYTID  = mpi ID of this process (from MPI_Comm_rank)
!
!        PARHER true when MPI parallel Hermit calculation
!                    (only set true in subroutines HERNOD and NODTRA)
!
!     .. openMP variables:
!         MAX_OMP_THREADS
   module infpar_mod
   use, intrinsic :: iso_c_binding, only: c_int, c_long
   implicit none
   private :: c_int, c_long, CINT
!  C kind of the default integer (C side: int, or long with INT_STAR8)
#if defined (INT_STAR8)
   integer, parameter :: CINT = c_long
#else
   integer, parameter :: CINT = c_int
#endif
   real*8 :: TOTWAL
   integer :: IPRPAR, NTASK, MYTID, NUMNOD, ICOUNT, NODE, IPAREX, MAX_OMP_THREADS
   integer(CINT), bind(C, name='diracinfpar_mparid') :: MPARID
   integer(CINT), bind(C, name='diracinfpar_nodes') :: NODES
   logical :: PARHER, SLAVE, TIMING
   end module infpar_mod
