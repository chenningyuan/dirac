!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

module xc_mpi

#ifdef VAR_MPI
   use mpi
#endif

   implicit none

   public xc_mpi_is_master
   public xc_mpi_get_nr_proc
#ifdef VAR_MPI
   public xc_mpi_bcast
   public xc_mpi_reduce
   public xc_mpi_distribute_points
#endif

   private

#ifdef VAR_MPI
   interface xc_mpi_bcast
      module procedure xc_mpi_bcast_i0
      module procedure xc_mpi_bcast_i1
      module procedure xc_mpi_bcast_i2
      module procedure xc_mpi_bcast_i3
      module procedure xc_mpi_bcast_i4
      module procedure xc_mpi_bcast_r0
      module procedure xc_mpi_bcast_r1
      module procedure xc_mpi_bcast_r2
      module procedure xc_mpi_bcast_r3
      module procedure xc_mpi_bcast_r4
      module procedure xc_mpi_bcast_l0
      module procedure xc_mpi_bcast_l1
      module procedure xc_mpi_bcast_l2
   end interface

   interface xc_mpi_reduce
      module procedure xc_mpi_reduce_i0
      module procedure xc_mpi_reduce_i1
      module procedure xc_mpi_reduce_r0
      module procedure xc_mpi_reduce_r1
      module procedure xc_mpi_reduce_r2
      module procedure xc_mpi_reduce_r4
   end interface
#endif

   integer            :: istat, ierr
   integer, parameter :: irank_master = 0

contains

#ifdef VAR_MPI
   function get_my_rank()
      integer :: get_my_rank
      integer :: irank
      call mpi_comm_rank(mpi_comm_world, irank, ierr)
      get_my_rank = irank
   end function
#endif

   function xc_mpi_is_master()
      logical :: xc_mpi_is_master
#ifdef VAR_MPI
      xc_mpi_is_master = (get_my_rank() == irank_master)
#else
      xc_mpi_is_master = .true.
#endif
   end function

   function xc_mpi_get_nr_proc()
      integer :: xc_mpi_get_nr_proc
      integer :: nr_proc
#ifdef VAR_MPI
      call mpi_comm_size(mpi_comm_world, nr_proc, ierr)
      xc_mpi_get_nr_proc = nr_proc
#else
      xc_mpi_get_nr_proc = 1
#endif
   end function

#ifdef VAR_MPI
   subroutine xc_mpi_bcast_i0(x)
      integer :: x
      call mpi_bcast(x, 1, mpi_integer, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_i1(x)
      integer :: x(:)
      call mpi_bcast(x, size(x), mpi_integer, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_i2(x)
      integer :: x(:, :)
      call mpi_bcast(x, size(x), mpi_integer, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_i3(x)
      integer :: x(:, :, :)
      call mpi_bcast(x, size(x), mpi_integer, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_i4(x)
      integer :: x(:, :, :, :)
      call mpi_bcast(x, size(x), mpi_integer, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_r0(x)
      real(8) :: x
      call mpi_bcast(x, 1, mpi_double_precision, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_r1(x)
      real(8) :: x(:)
      call mpi_bcast(x, size(x), mpi_double_precision, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_r2(x)
      real(8) :: x(:, :)
      call mpi_bcast(x, size(x), mpi_double_precision, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_r3(x)
      real(8) :: x(:, :, :)
      call mpi_bcast(x, size(x), mpi_double_precision, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_r4(x)
      real(8) :: x(:, :, :, :)
      call mpi_bcast(x, size(x), mpi_double_precision, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_l0(x)
      logical :: x
      call mpi_bcast(x, 1, mpi_logical, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_l1(x)
      logical :: x(:)
      call mpi_bcast(x, size(x), mpi_logical, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_bcast_l2(x)
      logical :: x(:, :)
      call mpi_bcast(x, size(x), mpi_logical, irank_master, mpi_comm_world, ierr)
   end subroutine

   subroutine xc_mpi_reduce_i0(x)
      integer :: x
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_i1(x)
      integer :: x(:)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r0(x)
      real(8) :: x
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r1(x)
      real(8) :: x(:)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r2(x)
      real(8) :: x(:, :)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      end if
   end subroutine

   subroutine xc_mpi_reduce_r4(x)
      real(8) :: x(:, :, :, :)
      if (xc_mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
      end if
   end subroutine

   subroutine xc_mpi_scatter_r1(x, n)
      real(8)              :: x(:)
      real(8), allocatable :: buffer(:)
      integer              :: n
   end subroutine

   subroutine xc_mpi_distribute_points(rx, ry, rz, rw, nr_points_batch, nr_points_on_this_proc)
      real(8)              :: rx(*)
      real(8)              :: ry(*)
      real(8)              :: rz(*)
      real(8)              :: rw(*)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable :: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable :: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/xc_mpi_get_nr_proc()
      if (xc_mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call mpi_scatter(rx, nr_points_on_this_proc, mpi_real8, buffer, &
              nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         call mpi_scatter(ry, nr_points_on_this_proc, mpi_real8, buffer, &
              nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         call mpi_scatter(rz, nr_points_on_this_proc, mpi_real8, buffer, &
              nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         call mpi_scatter(rw, nr_points_on_this_proc, mpi_real8, buffer, &
              nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         deallocate(buffer)
      else
         call mpi_scatter((/0.0d0/), 0, mpi_real8, rx, nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         call mpi_scatter((/0.0d0/), 0, mpi_real8, ry, nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         call mpi_scatter((/0.0d0/), 0, mpi_real8, rz, nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
         call mpi_scatter((/0.0d0/), 0, mpi_real8, rw, nr_points_on_this_proc, mpi_real8, irank_master, mpi_comm_world, ierr)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, xc_mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (xc_mpi_is_master()) then
            i = nr_points_batch - nr_points_rest + 1
            call mpi_send(rx(i), nr_points_rest, mpi_real8, xc_mpi_get_nr_proc()-1, 90, mpi_comm_world, ierr)
            call mpi_send(ry(i), nr_points_rest, mpi_real8, xc_mpi_get_nr_proc()-1, 91, mpi_comm_world, ierr)
            call mpi_send(rz(i), nr_points_rest, mpi_real8, xc_mpi_get_nr_proc()-1, 92, mpi_comm_world, ierr)
            call mpi_send(rw(i), nr_points_rest, mpi_real8, xc_mpi_get_nr_proc()-1, 93, mpi_comm_world, ierr)
         end if
         if (get_my_rank() == (xc_mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc + 1
            call mpi_recv(rx(i), nr_points_rest, mpi_real8, irank_master, 90, mpi_comm_world, istat, ierr)
            call mpi_recv(ry(i), nr_points_rest, mpi_real8, irank_master, 91, mpi_comm_world, istat, ierr)
            call mpi_recv(rz(i), nr_points_rest, mpi_real8, irank_master, 92, mpi_comm_world, istat, ierr)
            call mpi_recv(rw(i), nr_points_rest, mpi_real8, irank_master, 93, mpi_comm_world, istat, ierr)
            nr_points_on_this_proc = nr_points_on_this_proc + nr_points_rest
         end if
      end if

      allocate(dist_overview(0:xc_mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call xc_mpi_reduce(dist_overview)
      if (xc_mpi_is_master()) then
         do i = 0, xc_mpi_get_nr_proc()-1
!           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call xc_mpi_reduce(isum)
      if (xc_mpi_is_master()) then
         if (isum /= nr_points_batch) then
            call quit('programming error in distribution of points in parallel DFT')
         end if
      end if
   end subroutine
#endif /* ifdef VAR_MPI */

end module
