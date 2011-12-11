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
!
! purpose: determine integer size in MPI library to 
!
module interface_to_mpi
#ifdef VAR_MPI

  use integer_model_in_mpi
  use mpi

  implicit none

#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif


  public interface_mpi_bcast

  interface interface_mpi_bcast
    module procedure interface_mpi_bcast_i0
    module procedure interface_mpi_bcast_i1
    module procedure interface_mpi_bcast_i2
    module procedure interface_mpi_bcast_i3
    module procedure interface_mpi_bcast_i4
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
    call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

#else
  implicit none
  integer :: my_dummy_interface_to_mpi
#endif
end module interface_to_mpi
