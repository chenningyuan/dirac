module dft_mpi

  implicit none

#ifdef VAR_MPI
  public dft_mpi_bcast
#else
  public dft_mpi_empty
#endif

  private

#ifdef VAR_MPI
  interface dft_mpi_bcast
    module procedure dft_mpi_bcast_n
    module procedure dft_mpi_bcast_i1
    module procedure dft_mpi_bcast_i2
    module procedure dft_mpi_bcast_r1
    module procedure dft_mpi_bcast_r2
    module procedure dft_mpi_bcast_l
  end interface

#include "infpar.h"
#include "mpif.h"

  integer :: ierr
#endif

contains

#ifdef VAR_MPI
  subroutine dft_mpi_bcast_n(n)

!   ----------------------------------------------------------------------------
    integer :: n
!   ----------------------------------------------------------------------------

    call mpi_bcast(n, 1, my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dft_mpi_bcast_i1(x)

!   ----------------------------------------------------------------------------
    integer :: x(:)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dft_mpi_bcast_i2(x)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dft_mpi_bcast_r1(x)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dft_mpi_bcast_r2(x)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, mparid, mpi_comm_world, ierr)

  end subroutine

  subroutine dft_mpi_bcast_l(x)

!   ----------------------------------------------------------------------------
    logical :: x(:)
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_logical, mparid, mpi_comm_world, ierr)

  end subroutine
#endif /* ifdef VAR_MPI */

#ifndef VAR_MPI */
  subroutine dft_mpi_empty()
  end subroutine
#endif

end module
