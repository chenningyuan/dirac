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
! purpose: interface to MPI library via f77 interface
!
#ifdef VAR_MPI

      subroutine interface_mpi_bcast_l0_f77(x,ndim,root_proc,
     &                                      communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
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
!-------------------------------------------------------------------------------
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, my_MPI_LOGICAL, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      print *, ' enter l0 special: ndim == ',ndim
      call mpi_bcast(x, ndim   , my_MPI_LOGICAL, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_bcast_l1_f77(x,ndim,root_proc,
     &                                      communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      logical(kind=integer_kind) :: x(ndim)
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
      call mpi_bcast(x, ndim_i4, my_MPI_LOGICAL, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      print *, ' enter l1 special: ndim == ',ndim
      call mpi_bcast(x, ndim   , my_MPI_LOGICAL, root_proc   , 
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! special bcast routines for communicating data on the work array
!-------------------------------------------------------------------------------
    
      subroutine interface_mpi_bcast_i1_work_f77(x,ndim,root_proc,
     &                                           communicator)
    
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      integer(kind=integer_kind) :: x(ndim)
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
!------------------------------------------------------------------------------
      ndim_i4         = ndim
      root_proc_i4    = root_proc
      communicator_i4 = communicator
      call mpi_bcast(x, ndim_i4, my_MPI_INTEGER, root_proc_i4, 
     &               communicator_i4, ierr_i4)
#else
      print *, ' enter i1 special: ndim == ',ndim
      call mpi_bcast(x, ndim   , my_MPI_INTEGER, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------

      subroutine interface_mpi_bcast_r1_work_f77(x,ndim,root_proc,
     &                                           communicator)
!-------------------------------------------------------------------------------
      use integer_kind_mpilib
      implicit none
#include "mpif.h"
#include "interface_to_mpi_f77_typedef.h"
      real(8)                    :: x(ndim)
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
      call mpi_bcast(x, ndim_i4, MPI_REAL8, root_proc_i4,
     &               communicator_i4, ierr_i4)
#else
      print *, ' enter r1 special: ndim == ',ndim
      call mpi_bcast(x, ndim   , MPI_REAL8, root_proc   ,
     &               communicator   , ierr   )
#endif

      end subroutine
!-------------------------------------------------------------------------------
#else
      subroutine dummy_mpi_f77_interface
      end subroutine
#endif
