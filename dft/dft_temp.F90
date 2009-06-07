module dft_temp

! radovan: (1) if you modify this module (add new stuff)
!              don't forget to update
!              routines init_dft_temp and sync_dft_temp (below)
!          (2) the difference between dft_temp and dft_cfg is ...

  implicit none

  type dft_options_t
!   see also routine reset_dft_options
!   and sync routine(s) for parallel dft
    logical :: energy
    logical :: gradient
    logical :: london
    logical :: lr
    logical :: qr
    integer :: aor_xc_order
  end type

  save

  integer :: dft_temp_file_unit
  integer :: dft_temp_nr_grid_points
  integer :: dft_temp_nr_grid_points_node

  logical :: dft_temp_grid_is_done
  logical :: dft_temp_gga_hs
  logical :: dft_temp_gga_ha
  logical :: dft_temp_need_gao1
  logical :: dft_temp_reload_batches

! these are not broadcast, only mpi_reduced
  real(8) :: dft_temp_xc_energy
  real(8) :: dft_temp_xc_mat_energy
  real(8) :: dft_temp_nr_electrons

contains

  subroutine reset_dft_options(dft_options)

!   ----------------------------------------------------------------------------
    type(dft_options_t) :: dft_options
!   ----------------------------------------------------------------------------

    dft_options%energy       = .false.
    dft_options%gradient     = .false.
    dft_options%london       = .false.
    dft_options%lr           = .false.
    dft_options%qr           = .false.
    dft_options%aor_xc_order = 0

  end subroutine

  subroutine init_dft_temp()

    dft_temp_grid_is_done   = .false.
    dft_temp_gga_hs         = .false.
    dft_temp_gga_ha         = .false.
    dft_temp_need_gao1      = .false.
    dft_temp_reload_batches = .false.

  end subroutine

#ifdef VAR_MPI
  subroutine sync_dft_temp()

#include "mpif.h"
#include "infpar.h"

!   ----------------------------------------------------------------------------
    logical             :: i_am_the_master
    integer             :: idx, ierr
    logical             :: pack_logical(5)
!                                   1000 is a phantasy number
    integer, parameter  :: nbytes = 1000
    character           :: scratch(nbytes)
!   ----------------------------------------------------------------------------

    if (mytid == mparid) then
      i_am_the_master = .true.
    else
      i_am_the_master = .false.
    end if

    if (i_am_the_master) then

      pack_logical(1) = dft_temp_grid_is_done
      pack_logical(2) = dft_temp_gga_hs
      pack_logical(3) = dft_temp_gga_ha
      pack_logical(4) = dft_temp_need_gao1
      pack_logical(5) = dft_temp_reload_batches

      idx = 1
      call mpi_pack(pack_logical,       &
                    size(pack_logical), &
                    mpi_logical,        &
                    scratch,            &
                    nbytes,             &
                    idx,                &
                    mpi_comm_world,     &
                    ierr)

    end if

    call mpi_bcast(scratch,        &
                   nbytes,         &
                   mpi_packed,     &
                   0,              &
                   mpi_comm_world, &
                   ierr)

    if (.not. i_am_the_master) then

      idx = 1
      call mpi_unpack(scratch,            &
                      nbytes,             &
                      idx,                &
                      pack_logical,       &
                      size(pack_logical), &
                      mpi_logical,        &
                      mpi_comm_world,     &
                      ierr)

      dft_temp_grid_is_done   = pack_logical(1)
      dft_temp_gga_hs         = pack_logical(2)
      dft_temp_gga_ha         = pack_logical(3)
      dft_temp_need_gao1      = pack_logical(4)
      dft_temp_reload_batches = pack_logical(5)

    end if

  end subroutine

  subroutine sync_dft_options(dft_options)

#include "mpif.h"
#include "infpar.h"

!   ----------------------------------------------------------------------------
    type(dft_options_t), intent(inout) :: dft_options
!   ----------------------------------------------------------------------------
    logical             :: i_am_the_master
    integer             :: idx, ierr
    logical             :: pack_logical(5)
    integer             :: pack_integer(1)
!                                   1000 is a phantasy number
    integer, parameter  :: nbytes = 1000
    character           :: scratch(nbytes)
!   ----------------------------------------------------------------------------

    if (mytid == mparid) then
      i_am_the_master = .true.
    else
      i_am_the_master = .false.
    end if

    if (i_am_the_master) then

      pack_logical(1) = dft_options%energy
      pack_logical(2) = dft_options%gradient
      pack_logical(3) = dft_options%london
      pack_logical(4) = dft_options%lr
      pack_logical(5) = dft_options%qr
      pack_integer(1) = dft_options%aor_xc_order

      idx = 1
      call mpi_pack(pack_logical,       &
                    size(pack_logical), &
                    mpi_logical,        &
                    scratch,            &
                    nbytes,             &
                    idx,                &
                    mpi_comm_world,     &
                    ierr)
      call mpi_pack(pack_integer,       &
                    size(pack_integer), &
                    my_MPI_INTEGER,     &
                    scratch,            &
                    nbytes,             &
                    idx,                &
                    mpi_comm_world,     &
                    ierr)

    end if

    call mpi_bcast(scratch,        &
                   nbytes,         &
                   mpi_packed,     &
                   0,              &
                   mpi_comm_world, &
                   ierr)

    if (.not. i_am_the_master) then

      idx = 1
      call mpi_unpack(scratch,            &
                      nbytes,             &
                      idx,                &
                      pack_logical,       &
                      size(pack_logical), &
                      mpi_logical,        &
                      mpi_comm_world,     &
                      ierr)

      dft_options%energy       = pack_logical(1)
      dft_options%gradient     = pack_logical(2)
      dft_options%london       = pack_logical(3)
      dft_options%lr           = pack_logical(4)
      dft_options%qr           = pack_logical(5)
      dft_options%aor_xc_order = pack_integer(1)

    end if

  end subroutine
#endif /* ifdef VAR_MPI */

end module
