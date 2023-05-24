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
  public visual_interface_mpi_init
  public visual_interface_mpi_finalize
  public visual_interface_mpi_send
  public visual_interface_mpi_recv
  public visual_interface_mpi_bcast
  public visual_interface_mpi_reduce
  public visual_interface_mpi_gather
  public visual_interface_mpi_gatherv
  public visual_interface_mpi_distribute_points
  public visual_interface_mpi_wait

  interface visual_interface_mpi_send
    module procedure visual_interface_mpi_send_i0
    module procedure visual_interface_mpi_send_r0
    module procedure visual_interface_mpi_send_r1
  end interface

  interface visual_interface_mpi_recv
    !module procedure visual_interface_mpi_recv_i0
    module procedure visual_interface_mpi_recv_r0
    module procedure visual_interface_mpi_recv_r1
  end interface

  interface visual_interface_mpi_reduce
    module procedure visual_interface_mpi_reduce_i0
    module procedure visual_interface_mpi_reduce_i1
  end interface

  interface visual_interface_mpi_gather
    module procedure visual_interface_mpi_gather_r
  end interface

  interface visual_interface_mpi_bcast
    module procedure visual_interface_mpi_bcast_i0
    module procedure visual_interface_mpi_bcast_i1
    module procedure visual_interface_mpi_bcast_r2
    module procedure visual_interface_mpi_bcast_l0
    module procedure visual_interface_mpi_bcast_cx
    module procedure visual_interface_mpi_bcast_c1
  end interface

  interface visual_interface_mpi_distribute_points
    module procedure mpi_distribute_r1
  end interface

  interface visual_interface_mpi_scatter
    module procedure visual_interface_mpi_scatter_r1
  end interface
#endif

  save
  private

  integer, parameter, public:: irank_master = 0
#ifdef VAR_MPI
  integer, parameter, public:: my_MPI_INTEGER = MPI_INTEGER
  integer, parameter, public:: my_mpi_real8   = MPI_REAL8
  integer, parameter, public:: my_mpi_status_size  = mpi_status_size
  logical, public           :: MPI_INIT_called     = .false.
#endif

contains

  function is_parallel()
    logical:: is_parallel
    is_parallel = (mpi_get_nr_proc() > 1)
  end function

  function mpi_is_master()
    logical:: mpi_is_master
    mpi_is_master = .true.
    if (is_parallel()) then
      mpi_is_master = (get_my_rank() == irank_master)
    end if
  end function

  function mpi_get_nr_proc()
    integer:: mpi_get_nr_proc
    mpi_get_nr_proc = 1
#ifdef VAR_MPI
    call interface_mpi_comm_size(global_communicator,mpi_get_nr_proc)
#endif
   end function

#ifdef VAR_MPI
  function get_my_rank()
    integer:: get_my_rank
    call interface_mpi_comm_rank(global_communicator,get_my_rank)
  end function
#endif

#ifdef VAR_MPI
  subroutine visual_interface_mpi_init()
    integer:: ierr
    call MPI_Init(ierr)
    MPI_INIT_called = .true.
    call logprint
  end subroutine

  subroutine visual_interface_mpi_finalize()
    integer:: ierr
    if (.not. MPI_INIT_called) return
    call mpi_finalize(ierr)
  end subroutine


!-------------------------------------------------------------------------------
! mpi_send
!-------------------------------------------------------------------------------
  subroutine visual_interface_mpi_send_i0(x, ndim, receiver, message_id, communicator)
    integer:: x
    integer:: ndim
    integer:: receiver
    integer:: message_id
    integer:: communicator
    integer:: ierr

    if (.not. MPI_INIT_called) return

    call mpi_send(x, ndim   , my_MPI_integer, receiver   , message_id   , communicator   , ierr   )

  end subroutine


  subroutine visual_interface_mpi_send_r0(x, ndim, receiver, message_id, communicator)
    real(8)                    :: x
    integer:: ndim
    integer:: receiver
    integer:: message_id
    integer:: communicator
    integer:: ierr

    if (.not. MPI_INIT_called) return
    call mpi_send(x, ndim   , MPI_REAL8, receiver   , message_id   , communicator   , ierr   )

  end subroutine


  subroutine visual_interface_mpi_send_r1(x, ndim, receiver, message_id, communicator)
    real(8)                    :: x(:)
    integer:: ndim
    integer:: receiver
    integer:: message_id
    integer:: communicator
    integer:: ierr

    if (.not. MPI_INIT_called) return

    call mpi_send(x, ndim   , MPI_REAL8, receiver   , message_id   , communicator   , ierr   )

  end subroutine


!-------------------------------------------------------------------------------
! mpi_recv
!-------------------------------------------------------------------------------

  subroutine visual_interface_mpi_recv_r0(x, ndim, sender, message_id, communicator)
    real(8)                    :: x
    integer:: ndim
    integer:: sender
    integer:: message_id
    integer:: communicator
    integer:: ierr
    integer:: status_array(MPI_STATUS_SIZE)

    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_REAL8, sender   , message_id   , communicator   , status_array   , ierr   )

  end subroutine


  subroutine visual_interface_mpi_recv_r1(x, ndim, sender, message_id, communicator)
    real(8)                    :: x(:)
    integer:: ndim
    integer:: sender
    integer:: message_id
    integer:: communicator
    integer:: ierr
    integer:: status_array(MPI_STATUS_SIZE)

    if (.not. MPI_INIT_called) return
    call mpi_recv(x, ndim   , MPI_REAL8, sender   , message_id   , communicator   , status_array   , ierr   )

  end subroutine


!-------------------------------------------------------------------------------
! mpi_reduce
!-------------------------------------------------------------------------------

  subroutine visual_interface_mpi_reduce_i0(x)
    integer:: x
    integer:: ierr

    if (.not. MPI_INIT_called) return
      if (mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER, MPI_SUM, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, 1, MPI_INTEGER, MPI_SUM, 0, global_communicator, ierr)
      end if

  end subroutine


  subroutine visual_interface_mpi_reduce_i1(x)
    integer:: x(:)
    integer:: ierr
    if (.not. MPI_INIT_called) return
      if (mpi_is_master()) then
         call mpi_reduce(MPI_IN_PLACE, x, size(x), mpi_integer, mpi_sum, 0, global_communicator, ierr)
      else
         call mpi_reduce(x, MPI_IN_PLACE, size(x), mpi_integer, mpi_sum, 0, global_communicator, ierr)
      end if
  end subroutine



!-------------------------------------------------------------------------------
! mpi_gather
!-------------------------------------------------------------------------------

 subroutine visual_interface_mpi_gather_r(x, y)
    real(8):: x(:)
    real(8):: y(:)
    integer:: n
    integer:: disp
    integer:: ierr

    if (.not. MPI_INIT_called) return
    n = size(x)/mpi_get_nr_proc()
    call mpi_gather(x, n, MPI_REAL8, y, n, MPI_REAL8, irank_master, global_communicator   , ierr   )

  end subroutine


    subroutine visual_interface_mpi_gatherv(x, y, n)
      real(8):: x(:), y(:)
       integer, intent(in):: n
       integer:: i, ierr
       integer, allocatable:: rcounts(:)
       integer, allocatable:: displs(:)
       integer, allocatable:: rcounts_node(:)
       integer:: nproc

    if (.not. MPI_INIT_called) return

! set ncouts and displs arrays
!
       nproc = mpi_get_nr_proc()

       allocate(rcounts_node(nproc))
       allocate(rcounts(nproc))
       allocate(displs(nproc))

       rcounts_node = 0
       rcounts = 0
       rcounts_node(1) = n

       call mpi_gather(rcounts_node, 1, my_MPI_INTEGER, rcounts, 1, my_MPI_INTEGER, irank_master, global_communicator, ierr)

       displs(1) = 0
       do i = 2, nproc
             displs(i) = displs(i-1) + rcounts(i-1)
       end do

       call mpi_gatherv(x, n, mpi_double_precision, y, rcounts, displs, mpi_double_precision, &
    &                  irank_master, global_communicator, ierr)

       deallocate(rcounts_node)
       deallocate(rcounts)
       deallocate(displs)
    end subroutine




!-------------------------------------------------------------------------------
! mpi_bcast
!-------------------------------------------------------------------------------

  subroutine visual_interface_mpi_bcast_i0(x)
    integer:: x
    integer:: ierr

    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, 1, MPI_integer, irank_master, global_communicator, ierr   )

  end subroutine


  subroutine visual_interface_mpi_bcast_i1(x)
    integer:: x(:)
    integer:: ierr

    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, size(x), MPI_integer, irank_master, global_communicator, ierr   )

  end subroutine


  subroutine visual_interface_mpi_bcast_r2(x)
    real(8):: x(:,:)
    integer:: ierr

    if (.not. MPI_INIT_called) return
    call mpi_bcast(x, size(x), MPI_REAL8, irank_master, global_communicator, ierr)

  end subroutine


  subroutine visual_interface_mpi_bcast_l0(x)
    logical:: x
    integer:: ierr

     if (.not. MPI_INIT_called) return
     call mpi_bcast(x, 1, MPI_INTEGER, irank_master, global_communicator, ierr)

  end subroutine


  subroutine visual_interface_mpi_bcast_cx(x)
    character*(*):: x
    integer:: ierr

     if (.not. MPI_INIT_called) return
     call mpi_bcast(x, len(x), MPI_CHARACTER, irank_master, global_communicator, ierr)

  end subroutine


  subroutine visual_interface_mpi_bcast_c1(x)
    character*(*):: x(:)
    integer:: ierr

     if (.not. MPI_INIT_called) return
     !todo: test with size(x)*len(x)
     call mpi_bcast(x, size(x), MPI_CHARACTER, irank_master, global_communicator, ierr)

  end subroutine


!-------------------------------------------------------------------------------
! other subroutines
!-------------------------------------------------------------------------------


   subroutine mpi_distribute_r1(a, nr_points_batch, nr_points_on_this_proc)
      real(8)              :: a(:)
      integer              :: nr_points_batch
      integer              :: nr_points_on_this_proc
      real(8), allocatable:: buffer(:)
      integer              :: nr_points_rest
      integer              :: i
      integer, allocatable:: dist_overview(:)
      integer              :: isum

      nr_points_on_this_proc = nr_points_batch/mpi_get_nr_proc()
      if (mpi_is_master()) then
         allocate(buffer(nr_points_on_this_proc))
         call visual_interface_mpi_scatter(a, nr_points_on_this_proc, &
                                    buffer, nr_points_on_this_proc, &
                                    irank_master, global_communicator)
         deallocate(buffer)
      else
         call visual_interface_mpi_scatter((/0.0d0/), 0, &
                                    a, nr_points_on_this_proc, &
                                    irank_master, global_communicator)
      end if

!     if there is a rest of points, send them to the last proc
      nr_points_rest = mod(nr_points_batch, mpi_get_nr_proc())
      if (nr_points_rest > 0) then
         if (mpi_is_master()) then
            i = nr_points_batch-nr_points_rest+1
            call visual_interface_mpi_send(a(i), nr_points_rest, mpi_get_nr_proc()-1, 90, global_communicator)
         end if
         if (get_my_rank() == (mpi_get_nr_proc()-1)) then
            i = nr_points_on_this_proc+1
            call visual_interface_mpi_recv(a(i), nr_points_rest, irank_master, 90, global_communicator)
            nr_points_on_this_proc = nr_points_on_this_proc+nr_points_rest
         end if
      end if

      allocate(dist_overview(0:mpi_get_nr_proc()-1))
      dist_overview = 0
      dist_overview(get_my_rank()) = nr_points_on_this_proc
      call visual_interface_mpi_reduce(dist_overview)

      if (mpi_is_master()) then
         do i = 0, mpi_get_nr_proc()-1
           write(*, *) 'proc', i+1, 'will calculate', dist_overview(i), 'points'
         end do
      end if
      deallocate(dist_overview)

!     verify that sum matches nr_points_batch
      isum = nr_points_on_this_proc
      call visual_interface_mpi_reduce(isum)
      if (mpi_is_master()) then
         if (isum /= nr_points_batch) then
            write(*,*) 'programming error in distribution of points'
            stop
         end if
      end if
   end subroutine




  subroutine visual_interface_mpi_wait(request, status_container)
    integer:: request
    integer:: status_container(MPI_STATUS_SIZE)
    integer:: ierr

    if (.not. MPI_INIT_called) return
    call mpi_wait(request   , status_container   , ierr   )

  end subroutine


!-------------------------------------------------------------------------------
! private
!-------------------------------------------------------------------------------

  subroutine logprint

  if (MPI_INIT_called) then

    if (mpi_is_master()) then
      write(*,*) 'mpi is initialized'
      if (is_parallel()) then
        write(*,*) 'this is a parallel run'
      end if
    end if

    if (get_my_rank() == irank_master) then
      write(*,*) 'size of communicator', mpi_get_nr_proc()
    else
      write(*,*) 'rank of the process on node ', get_my_rank()
    end if

  end if

  end subroutine logprint


  subroutine visual_interface_mpi_scatter_r1(x, ndim1, y, ndim2, root_proc, communicator)
    real(8)                    :: x(:)
    real(8)                    :: y(:)
    integer:: ndim1
    integer:: ndim2
    integer:: root_proc
    integer:: communicator
    integer:: ierr

    if (.not. MPI_INIT_called) return
    call mpi_scatter(x, ndim1   , MPI_REAL8, y, ndim2   , MPI_REAL8, root_proc   , communicator   , ierr   )

  end subroutine

#endif /* ifdef VAR_MPI */
end module
