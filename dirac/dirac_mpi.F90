!      Copyright (c) 2010 by the authors of DIRAC.
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

module dirac_mpi

  implicit none

#ifdef VAR_MPI
  public dirac_mpi_bcast

  private

  interface dirac_mpi_bcast
    module procedure dirac_mpi_bcast_i0
    module procedure dirac_mpi_bcast_i1
    module procedure dirac_mpi_bcast_i2
    module procedure dirac_mpi_bcast_i3
    module procedure dirac_mpi_bcast_i4
    module procedure dirac_mpi_bcast_r0
    module procedure dirac_mpi_bcast_r1
    module procedure dirac_mpi_bcast_r2
    module procedure dirac_mpi_bcast_r3
    module procedure dirac_mpi_bcast_r4
    module procedure dirac_mpi_bcast_l0
    module procedure dirac_mpi_bcast_l1
    module procedure dirac_mpi_bcast_l2
  end interface

#include "infpar.h"
#include "mpif.h"

  integer :: ierr

contains

  subroutine dirac_mpi_bcast_i0(x)

!   ----------------------------------------------------------------------------
    integer :: x
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, 1, my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_i1(x)

!   ----------------------------------------------------------------------------
    integer :: x(:)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_i2(x)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_i3(x)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_i4(x)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_r0(x)

!   ----------------------------------------------------------------------------
    real(8) :: x
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, 1, mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_r1(x)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_r2(x)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_r3(x)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_r4(x)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_l0(x)

!   ----------------------------------------------------------------------------
    logical :: x
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, 1, mpi_logical, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_l1(x)

!   ----------------------------------------------------------------------------
    logical :: x(:)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_logical, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dirac_mpi_bcast_l2(x)

!   ----------------------------------------------------------------------------
    logical :: x(:, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_logical, mparid, mpi_comm_world, ierr)

  end subroutine
#endif /* ifdef VAR_MPI */

end module
