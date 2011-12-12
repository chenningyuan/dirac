!dirac_copyright_start
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
!dirac_copyright_end
#define VAR_MPI_32bit_int
!
! purpose: determine integer size in MPI library to 
!
module interface_to_mpi
#ifdef VAR_MPI

  use integer_model_in_mpi
  use mpi

  implicit none

#include "priunit.h"

#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#define my_MPI_LOGICAL MPI_LOGICAL8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#define my_MPI_LOGICAL MPI_LOGICAL4
#endif


  public interface_mpi_bcast
  public interface_mpi_gather
  public interface_mpi_comm_rank
  public interface_mpi_comm_size
  public interface_mpi_iprobe

  interface interface_mpi_bcast
    module procedure interface_mpi_bcast_i0
    module procedure interface_mpi_bcast_i1
    module procedure interface_mpi_bcast_i2
    module procedure interface_mpi_bcast_i3
    module procedure interface_mpi_bcast_i4
    module procedure interface_mpi_bcast_r0
    module procedure interface_mpi_bcast_r1
    module procedure interface_mpi_bcast_r2
    module procedure interface_mpi_bcast_r3
    module procedure interface_mpi_bcast_r4
    module procedure interface_mpi_bcast_cx
    module procedure interface_mpi_bcast_l0
    module procedure interface_mpi_bcast_l1
  end interface

  interface interface_mpi_gather
    module procedure interface_mpi_gather_i
    module procedure interface_mpi_gather_r
  end interface
  
  interface interface_mpi_iprobe
    module procedure interface_mpi_iprobe_l0
  end interface

  private

#ifdef INT_STAR8
  integer, parameter :: integer_kind = 8
#else
  integer, parameter :: integer_kind = 4
#endif

  save

contains

!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i0(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter i0: ndim == ',ndim 
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter i1: ndim == ',ndim 
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i2(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter i2: ndim == ',ndim 
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i3(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:,:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter i3: ndim == ',ndim 
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_i4(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x(:,:,:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter i4: ndim == ',ndim 
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r0(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter r0: ndim == ',ndim 
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter r1: ndim == ',ndim 
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r2(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter r2: ndim == ',ndim 
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r3(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter r3: ndim == ',ndim 
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_r4(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x(:,:,:,:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter r4: ndim == ',ndim 
    call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_cx(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    character*(*)              :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_bcast(x, ndim_i4, MPI_CHARACTER, root_proc_i4, communicator_i4, ierr_i4)
#else
    call mpi_bcast(x, ndim   , MPI_CHARACTER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_l0(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    logical(kind=integer_kind) :: x
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    logical(kind=4)            :: x_l4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    x_l4            = x
    call mpi_bcast(x_l4, ndim_i4, my_MPI_LOGICAL, root_proc_i4, communicator_i4, ierr_i4)
    x = x_l4
    write(lupri,*) ' check l0 - this might not necessarily work!!!: ndim == ',ndim
#else
    write(lupri,*) ' enter l0: ndim == ',ndim 
    call mpi_bcast(x, ndim   , my_MPI_LOGICAL, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_l1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    logical                      :: x(:)
    integer(kind=integer_kind)   :: ndim
    integer(kind=integer_kind)   :: root_proc
    integer(kind=integer_kind)   :: communicator
    integer(kind=integer_kind)   :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)              :: ndim_i4
    integer(kind=4)              :: root_proc_i4
    integer(kind=4)              :: communicator_i4
    integer(kind=4)              :: ierr_i4
    logical(kind=4), allocatable :: x_l4(:)
    integer                      :: i
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    allocate(x_l4(ndim_i4))
    do i = 1, ndim_i4
      x_l4(i) = x(i)
    end do
    call mpi_bcast(x_l4, ndim_i4, my_MPI_LOGICAL, root_proc_i4, communicator_i4, ierr_i4)
    write(lupri,*) ' check l1 - this might not necessarily work!!!: ndim == ',ndim
    do i = 1, ndim_i4
      x(i) = x_l4(i)
    end do
    deallocate(x_l4)
#else
    call mpi_bcast(x, ndim   , my_MPI_LOGICAL, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_iprobe_l0(rank,message_id,x,communicator,status_container)

!-------------------------------------------------------------------------------
    logical(kind=integer_kind) :: x
    integer(kind=integer_kind) :: message_id
    integer(kind=integer_kind) :: status_container(MPI_STATUS_SIZE)
    integer(kind=integer_kind) :: rank
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    logical(kind=4)            :: x_l4
    integer(kind=4)            :: message_id_i4
    integer(kind=4)            :: status_container_i4(MPI_STATUS_SIZE)
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: rank_i4
    integer(kind=4)            :: i
!-------------------------------------------------------------------------------
    communicator_i4 = communicator
    message_id_i4   = message_id
    rank_i4         = rank

    call mpi_iprobe(rank_i4, message_id_i4, communicator_i4, x_l4, status_container_i4, ierr_i4)

    do i = 1, MPI_STATUS_SIZE
      status_container(i) = status_container_i4(i)
    end do
    x               = x_l4
#else
    call mpi_iprobe(rank   , message_id   , communicator   , x   , status_container   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! mpi_send interfaces
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! mpi_recv interfaces
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! mpi_reduce interfaces
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! mpi_gather interfaces
!-------------------------------------------------------------------------------

  subroutine interface_mpi_gather_i(x,ndim,y,disp,root_proc,communicator)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: x
    integer(kind=integer_kind) :: y(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: disp
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: disp_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    disp_i4         = disp
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_gather(x, ndim_i4, my_MPI_INTEGER, y, disp_i4, my_MPI_INTEGER, root_proc_i4, communicator_i4, ierr_i4)
#else
    call mpi_gather(x, ndim   , my_MPI_INTEGER, y, disp   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_gather_r(x,ndim,y,disp,root_proc,communicator)

!-------------------------------------------------------------------------------
    real(8)                    :: x
    real(8)                    :: y(:)
    integer(kind=integer_kind) :: ndim
    integer(kind=integer_kind) :: disp
    integer(kind=integer_kind) :: root_proc
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: ndim_i4
    integer(kind=4)            :: disp_i4
    integer(kind=4)            :: root_proc_i4
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
!-------------------------------------------------------------------------------
    ndim_i4         = ndim
    disp_i4         = disp
    root_proc_i4    = root_proc
    communicator_i4 = communicator
    call mpi_gather(x, ndim_i4, MPI_REAL8, y, disp_i4, MPI_REAL8, root_proc_i4, communicator_i4, ierr_i4)
#else
    call mpi_gather(x, ndim   , MPI_REAL8, y, disp   , MPI_REAL8, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! interfaces to non-buffer MPI routines
!-------------------------------------------------------------------------------

  subroutine interface_mpi_comm_rank(communicator,rank)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: rank
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: rank_i4
!-------------------------------------------------------------------------------
    communicator_i4 = communicator
    call mpi_comm_rank(communicator_i4, rank_i4, ierr_i4)
    rank            = rank_i4
#else
    call mpi_comm_rank(communicator   , rank   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_comm_size(communicator,is_size)

!-------------------------------------------------------------------------------
    integer(kind=integer_kind) :: is_size
    integer(kind=integer_kind) :: communicator
    integer(kind=integer_kind) :: ierr
!-------------------------------------------------------------------------------
#ifdef VAR_MPI_32bit_int
    integer(kind=4)            :: communicator_i4
    integer(kind=4)            :: ierr_i4
    integer(kind=4)            :: is_size_i4
!-------------------------------------------------------------------------------
    communicator_i4 = communicator
    call mpi_comm_rank(communicator_i4, is_size_i4, ierr_i4)
    is_size         = is_size_i4
#else
    call mpi_comm_rank(communicator   , is_size      , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------


#else
  implicit none
  integer :: my_dummy_interface_to_mpi
#endif
end module interface_to_mpi
