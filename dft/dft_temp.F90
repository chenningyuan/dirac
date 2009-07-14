module dft_temp

! radovan: (1) if you modify this module (add new stuff)
!              don't forget to update
!              routines init_dft_temp and sync_dft_temp (below)
!          (2) the difference between dft_temp and dft_cfg is ...

#ifdef VAR_MPI
  use dft_mpi
#endif

  implicit none

  public reset_dft_options
#ifdef VAR_MPI
  public sync_dft_options
  public sync_dft_temp
#endif

  public dft_options_t

  type dft_options_t
!   see also routine reset_dft_options
!   and sync routine(s) for parallel dft
    logical :: energy       = .false.
    logical :: gradient     = .false.
    logical :: london       = .false.
    logical :: lr           = .false.
    logical :: qr           = .false.
    integer :: aor_xc_order = 0
  end type

  save

! only master needs this
  integer, public :: dft_temp_file_unit = 0

  logical, public :: dft_temp_grid_is_done   = .false.
  logical, public :: dft_temp_gga_hs         = .false.
  logical, public :: dft_temp_gga_ha         = .false.
  logical, public :: dft_temp_reload_batches = .false.

! these are not broadcast, only mpi_reduced
  real(8), public :: dft_temp_xc_energy     = 0.0d0
  real(8), public :: dft_temp_xc_mat_energy = 0.0d0
  real(8), public :: dft_temp_nr_electrons  = 0.0d0

  private

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

#ifdef VAR_MPI
  subroutine sync_dft_temp()

    call dft_mpi_bcast(dft_temp_grid_is_done  )
    call dft_mpi_bcast(dft_temp_gga_hs        )
    call dft_mpi_bcast(dft_temp_gga_ha        )
    call dft_mpi_bcast(dft_temp_reload_batches)

  end subroutine

  subroutine sync_dft_options(dft_options)

!   ----------------------------------------------------------------------------
    type(dft_options_t), intent(inout) :: dft_options
!   ----------------------------------------------------------------------------

    call dft_mpi_bcast(dft_options%energy      )
    call dft_mpi_bcast(dft_options%gradient    )
    call dft_mpi_bcast(dft_options%london      )
    call dft_mpi_bcast(dft_options%lr          )
    call dft_mpi_bcast(dft_options%qr          )
    call dft_mpi_bcast(dft_options%aor_xc_order)

  end subroutine
#endif /* ifdef VAR_MPI */

end module
