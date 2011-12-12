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

#include "priunit.h"

#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif


  public interface_mpi_bcast
  public interface_mpi_gather

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
    logical                    :: x
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
    call mpi_bcast(x, ndim_i4, MPI_LOGICAL, root_proc_i4, communicator_i4, ierr_i4)
#else
    write(lupri,*) ' enter i1: ndim == ',ndim 
    call mpi_bcast(x, ndim   , MPI_LOGICAL, root_proc   , communicator   , ierr   )
#endif

  end subroutine
!-------------------------------------------------------------------------------

  subroutine interface_mpi_bcast_l1(x,ndim,root_proc,communicator)

!-------------------------------------------------------------------------------
    logical                    :: x(:)
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
    call mpi_bcast(x, ndim_i4, MPI_LOGICAL, root_proc_i4, communicator_i4, ierr_i4)
#else
    call mpi_bcast(x, ndim   , MPI_LOGICAL, root_proc   , communicator   , ierr   )
#endif

  end subroutine
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
#else
  implicit none
  integer :: my_dummy_interface_to_mpi
#endif
end module interface_to_mpi
