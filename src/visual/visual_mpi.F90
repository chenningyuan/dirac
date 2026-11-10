module visual_mpi

#ifdef VAR_MPI
  use interface_to_mpi
  use integer_kind_mpilib
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif
#else
  implicit none
#endif

  public mpi_is_master
  public mpi_get_nr_proc
  public is_parallel
#ifdef VAR_MPI
  public get_my_rank
  public visual_interface_mpi_gatherv
  public visual_interface_mpi_distribute_points
  public visual_interface_mpi_reduce

  interface visual_interface_mpi_distribute_points
    module procedure mpi_distribute_r1
    module procedure mpi_distribute_r3
  end interface

  interface visual_interface_mpi_reduce
    module procedure mpi_reduce_r0
  end interface
#endif

  save
  private

  integer, parameter, public:: irank_master = 0

contains

  function is_parallel()
    logical:: is_parallel
    is_parallel = (mpi_get_nr_proc() > 1)
  end function

  function mpi_is_master()
    logical:: mpi_is_master
    mpi_is_master = .true.
    if (is_parallel()) then
#ifdef VAR_MPI
      mpi_is_master = (get_my_rank() == irank_master)
#endif
    end if
  end function

  function mpi_get_nr_proc()
#if defined (VAR_MPI)
      use dcbgen_mod, only: PARCAL
#endif
    integer:: mpi_get_nr_proc
    integer :: nr_proc
    mpi_get_nr_proc = 1
#ifdef VAR_MPI
    if (parcal) then
      call interface_mpi_comm_size(global_communicator, nr_proc)
      mpi_get_nr_proc = nr_proc
    end if
#endif
  end function

#ifdef VAR_MPI
  function get_my_rank()
    integer:: get_my_rank
    call interface_mpi_comm_rank(global_communicator,get_my_rank)
  end function
#endif

#ifdef VAR_MPI
  subroutine mpi_reduce_r0(x, y)
    real(8):: x, y

    call interface_mpi_reduce(x,y,1,op_mpi_sum,&
                              irank_master,global_communicator)
  end subroutine

  subroutine visual_interface_mpi_gatherv(x, y, n)
    real(8):: x(:), y(:)
    integer, intent(in):: n
    integer:: i, ierr
    integer, allocatable:: rcounts(:)
    integer, allocatable:: displs(:)
    integer, allocatable:: rcounts_node(:)
    integer:: nproc

!   set ncouts and displs arrays
!
    nproc = mpi_get_nr_proc()

    allocate(rcounts_node(nproc))
    allocate(rcounts(nproc))
    allocate(displs(nproc))

    rcounts_node = 0
    rcounts = 0
    rcounts_node(1) = n

!
!   gosia 08/06/2023:
!   in the calls to generic mpi subroutine, we use mpi_double_precision
!   which is different than mpi_real8 used in interface_to_mpi.F90
!   (see https://stackoverflow.com/questions/65579359/mpi-double-precision-vs-mpi-real8
!   and https://stackoverflow.com/questions/19726509/how-do-i-retain-precision-for-a-fortran-mpi-program-in-a-portable-way)
!   i didn't test it here, i kept how things were done in embedding/emb_mpi.F90 subroutines
!
    !call interface_mpi_gather(rcounts_node,1,rcounts,1,irank_master,global_communicator)
    call mpi_gather(rcounts_node, 1, df_my_mpi_integer, rcounts, 1, df_my_mpi_integer, &
                    irank_master, global_communicator, ierr)

    displs(1) = 0
    do i = 2, nproc
          displs(i) = displs(i-1) + rcounts(i-1)
    end do

    !call interface_mpi_gatherv(x,n,y,rcounts,displs,irank_master,global_communicator)
    call mpi_gatherv(x, n, mpi_double_precision, y, rcounts, displs, mpi_double_precision, &
                     irank_master, global_communicator, ierr)

    deallocate(rcounts_node)
    deallocate(rcounts)
    deallocate(displs)
  end subroutine


  subroutine mpi_distribute_r1(a, nr_points_batch, nr_points_on_this_proc)
    real(8)              :: a(*)
    integer              :: nr_points_batch
    integer              :: nr_points_on_this_proc
    real(8), allocatable:: buffer(:)
    integer              :: nr_points_rest
    integer              :: i
    integer              :: isum=0

    nr_points_on_this_proc = nr_points_batch/mpi_get_nr_proc()
    if (mpi_is_master()) then
      allocate(buffer(nr_points_on_this_proc))
      call interface_mpi_scatter_r1_work_f77(a, nr_points_on_this_proc, &
                                              buffer, nr_points_on_this_proc, &
                                              irank_master, global_communicator)
      deallocate(buffer)
    else
      call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, &
                                             a, nr_points_on_this_proc, &
                                             irank_master, global_communicator)
    end if

!   if there is a rest of points, send them to the last proc
    nr_points_rest = mod(nr_points_batch, mpi_get_nr_proc())
    if (nr_points_rest > 0) then
      if (mpi_is_master()) then
        i = nr_points_batch-nr_points_rest+1
        call interface_mpi_send(a(i), nr_points_rest, mpi_get_nr_proc()-1, 111, global_communicator)
      end if
      if (get_my_rank() == (mpi_get_nr_proc()-1)) then
        i = nr_points_on_this_proc+1
        call interface_mpi_recv(a(i), nr_points_rest, irank_master, 111, global_communicator)
        nr_points_on_this_proc = nr_points_on_this_proc+nr_points_rest
      end if
    end if

!   verify that sum matches nr_points_batch
    isum = nr_points_on_this_proc
    call interface_mpi_reduce(nr_points_on_this_proc,isum,1,op_mpi_sum,&
                              irank_master,global_communicator)
    if (mpi_is_master()) then
      if (isum /= nr_points_batch) then
        write(*,*) 'programming error in distribution of points'
        stop
      end if
    end if
  end subroutine

  subroutine mpi_distribute_r3(ax, ay, az, nr_points_batch, nr_points_on_this_proc)
    real(8)              :: ax(*)
    real(8)              :: ay(*)
    real(8)              :: az(*)
    integer              :: nr_points_batch
    integer              :: nr_points_on_this_proc
    real(8), allocatable:: buffer(:)
    integer              :: nr_points_rest
    integer              :: i
    integer              :: isum=0

    nr_points_on_this_proc = nr_points_batch/mpi_get_nr_proc()
    if (mpi_is_master()) then
      allocate(buffer(nr_points_on_this_proc))
      call interface_mpi_scatter_r1_work_f77(ax, nr_points_on_this_proc, &
                                              buffer, nr_points_on_this_proc, &
                                              irank_master, global_communicator)
      call interface_mpi_scatter_r1_work_f77(ay, nr_points_on_this_proc, &
                                              buffer, nr_points_on_this_proc, &
                                              irank_master, global_communicator)
      call interface_mpi_scatter_r1_work_f77(az, nr_points_on_this_proc, &
                                              buffer, nr_points_on_this_proc, &
                                              irank_master, global_communicator)
      deallocate(buffer)
    else
      call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, &
                                             ax, nr_points_on_this_proc, &
                                             irank_master, global_communicator)
      call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, &
                                             ay, nr_points_on_this_proc, &
                                             irank_master, global_communicator)
      call interface_mpi_scatter_r1_work_f77((/0.0d0/), 0, &
                                             az, nr_points_on_this_proc, &
                                             irank_master, global_communicator)
    end if

!   if there is a rest of points, send them to the last proc
    nr_points_rest = mod(nr_points_batch, mpi_get_nr_proc())
    if (nr_points_rest > 0) then
      if (mpi_is_master()) then
        i = nr_points_batch-nr_points_rest+1
        call interface_mpi_send(ax(i), nr_points_rest, mpi_get_nr_proc()-1, 117, global_communicator)
        call interface_mpi_send(ay(i), nr_points_rest, mpi_get_nr_proc()-1, 118, global_communicator)
        call interface_mpi_send(az(i), nr_points_rest, mpi_get_nr_proc()-1, 119, global_communicator)
      end if
      if (get_my_rank() == (mpi_get_nr_proc()-1)) then
        i = nr_points_on_this_proc+1
        call interface_mpi_recv(ax(i), nr_points_rest, irank_master, 117, global_communicator)
        call interface_mpi_recv(ay(i), nr_points_rest, irank_master, 118, global_communicator)
        call interface_mpi_recv(az(i), nr_points_rest, irank_master, 119, global_communicator)
        nr_points_on_this_proc = nr_points_on_this_proc+nr_points_rest
      end if
    end if

!   verify that sum matches nr_points_batch
    isum = nr_points_on_this_proc
    call interface_mpi_reduce(nr_points_on_this_proc,isum,1,op_mpi_sum,&
                              irank_master,global_communicator)
    if (mpi_is_master()) then
      if (isum /= nr_points_batch) then
        write(*,*) 'programming error in distribution of points'
        stop
      end if
    end if
  end subroutine


#endif /* ifdef VAR_MPI */
end module
