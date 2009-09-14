module dft_cfg

! radovan: (1) if you modify this module (add new stuff)
!              don't forget to update
!              sync_dft_cfg (below)
!          (2) this module should contain things set by the user
!              or in principle accessible by the user
!              temporary parameters used by the program
!              (like integrated electron count etc.)
!              should not be here, put them in dft_temp

#ifdef VAR_MPI
  use dft_mpi
#endif

  implicit none

#ifdef VAR_MPI
  public sync_dft_cfg
#endif

  save

  character(30), public :: dft_cfg_gridfile

  logical, public :: dft_cfg_keep_ao                 = .false.
  logical, public :: dft_cfg_explicit_potential      = .false.
  logical, public :: dft_cfg_no_pruning              = .false.
  logical, public :: dft_cfg_estimate_radii          = .false.
  logical, public :: dft_cfg_import_grid             = .false.
  logical, public :: dft_cfg_export_grid             = .false.
  logical, public :: dft_cfg_nozip                   = .false.
  logical, public :: dft_cfg_alda_hs                 = .false.
  logical, public :: dft_cfg_alda_ha                 = .false.
  logical, public :: dft_cfg_xalda                   = .false.
  logical, public :: dft_cfg_force_4c_grid           = .false.
  logical, public :: dft_cfg_old_dft                 = .false.
  logical, public :: dft_cfg_sdft                    = .false.
  logical, public :: dft_cfg_sdft_collinear          = .false.
  logical, public :: dft_cfg_sdft_beta               = .false.
  logical, public :: dft_cfg_cdft_vrg                = .false.
  logical, public :: dft_cfg_autodiff                = .false.
  logical, public :: dft_cfg_grac                    = .false.
  logical, public :: dft_cfg_grac_is_active          = .false.
  logical, public :: dft_cfg_saop                    = .false.
  logical, public :: dft_cfg_saop_is_active          = .false.
  logical, public :: dft_cfg_saop_with_response_part = .false.
  logical, public :: dft_cfg_asymptote_is_lb94       = .false.
  logical, public :: dft_cfg_asymptote_is_lbalpha    = .false.
  logical, public :: dft_cfg_is_gga_autodiff         = .false.

  real(8), public :: dft_cfg_radint                          = 1.0d-13
  real(8), public :: dft_cfg_accepted_error                  = 1.0d-3
  real(8), public :: dft_cfg_dfthr0                          = 1.0d-9
  real(8), public :: dft_cfg_dfthrl                          = 1.0d-10
  real(8), public :: dft_cfg_dfthri                          = 1.0d-15
  real(8), public :: dft_cfg_ac_ip                           = 0.0d0
  real(8), public :: dft_cfg_ac_threshold                    = 0.0d0
  real(8), public :: dft_cfg_homo_eigenvalue                 = 0.0d0
  real(8), public :: dft_cfg_grac_alpha                      = 0.0d0
  real(8), public :: dft_cfg_grac_beta                       = 0.0d0
  real(8), public :: dft_cfg_lda_functional_weight_array(20) = 0.0d0
  real(8), public :: dft_cfg_gga_functional_weight_array(20) = 0.0d0

  integer, public :: dft_cfg_angmin                  = 15
  integer, public :: dft_cfg_angint                  = 41
  integer, public :: dft_cfg_ao_mem_limit            = 0
  integer, public :: dft_cfg_print_level             = 2
  integer, public :: dft_cfg_debug_level             = 0
  integer, public :: dft_cfg_integration_check_level = 0

  private

#ifdef VAR_MPI
contains

  subroutine sync_dft_cfg()

    call dft_mpi_bcast(dft_cfg_keep_ao                    )
    call dft_mpi_bcast(dft_cfg_explicit_potential         )
    call dft_mpi_bcast(dft_cfg_no_pruning                 )
    call dft_mpi_bcast(dft_cfg_estimate_radii             )
    call dft_mpi_bcast(dft_cfg_import_grid                )
    call dft_mpi_bcast(dft_cfg_export_grid                )
    call dft_mpi_bcast(dft_cfg_nozip                      )
    call dft_mpi_bcast(dft_cfg_alda_hs                    )
    call dft_mpi_bcast(dft_cfg_alda_ha                    )
    call dft_mpi_bcast(dft_cfg_xalda                      )
    call dft_mpi_bcast(dft_cfg_force_4c_grid              )
    call dft_mpi_bcast(dft_cfg_old_dft                    )
    call dft_mpi_bcast(dft_cfg_sdft                       )
    call dft_mpi_bcast(dft_cfg_sdft_collinear             )
    call dft_mpi_bcast(dft_cfg_sdft_beta                  )
    call dft_mpi_bcast(dft_cfg_cdft_vrg                   )
    call dft_mpi_bcast(dft_cfg_autodiff                   )
    call dft_mpi_bcast(dft_cfg_grac                       )
    call dft_mpi_bcast(dft_cfg_grac_is_active             )
    call dft_mpi_bcast(dft_cfg_saop                       )
    call dft_mpi_bcast(dft_cfg_saop_is_active             )
    call dft_mpi_bcast(dft_cfg_saop_with_response_part    )
    call dft_mpi_bcast(dft_cfg_asymptote_is_lb94          )
    call dft_mpi_bcast(dft_cfg_asymptote_is_lbalpha       )
    call dft_mpi_bcast(dft_cfg_is_gga_autodiff            )
    call dft_mpi_bcast(dft_cfg_radint                     )
    call dft_mpi_bcast(dft_cfg_accepted_error             )
    call dft_mpi_bcast(dft_cfg_dfthr0                     )
    call dft_mpi_bcast(dft_cfg_dfthrl                     )
    call dft_mpi_bcast(dft_cfg_dfthri                     )
    call dft_mpi_bcast(dft_cfg_ac_ip                      )
    call dft_mpi_bcast(dft_cfg_ac_threshold               )
    call dft_mpi_bcast(dft_cfg_homo_eigenvalue            )
    call dft_mpi_bcast(dft_cfg_grac_alpha                 )
    call dft_mpi_bcast(dft_cfg_grac_beta                  )
    call dft_mpi_bcast(dft_cfg_lda_functional_weight_array)
    call dft_mpi_bcast(dft_cfg_gga_functional_weight_array)
    call dft_mpi_bcast(dft_cfg_angmin                     )
    call dft_mpi_bcast(dft_cfg_angint                     )
    call dft_mpi_bcast(dft_cfg_ao_mem_limit               )
    call dft_mpi_bcast(dft_cfg_print_level                )
    call dft_mpi_bcast(dft_cfg_debug_level                )
    call dft_mpi_bcast(dft_cfg_integration_check_level    )

  end subroutine
#endif /* ifdef VAR_MPI */

end module
