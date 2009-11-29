module dft_temp

! radovan: (1) if you modify this module (add new stuff)
!              don't forget to update
!              routines init_dft_temp and sync_dft_temp (below)
!          (2) the difference between dft_temp and dft_cfg is ...

#ifdef VAR_MPI
  use dirac_mpi
#endif

  implicit none

#ifdef VAR_MPI
  public sync_dft_temp
#endif

  save

! only master needs this
  integer, public :: dft_temp_file_unit = 0

  logical, public :: dft_temp_grid_is_done   = .false.
  logical, public :: dft_temp_reload_batches = .false.

! these are not broadcast, only mpi_reduced
  real(8), public :: dft_temp_xc_energy     = 0.0d0
  real(8), public :: dft_temp_xc_mat_energy = 0.0d0
  real(8), public :: dft_temp_nr_electrons  = 0.0d0

  private

contains

#ifdef VAR_MPI
  subroutine sync_dft_temp()

    call dirac_mpi_bcast(dft_temp_grid_is_done  )
    call dirac_mpi_bcast(dft_temp_reload_batches)

  end subroutine
#endif /* ifdef VAR_MPI */

end module
