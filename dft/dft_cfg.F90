!dirac_copyright_start
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
!dirac_copyright_end

module dft_cfg

! radovan: (1) if you modify this module (add new stuff)
!              don't forget to update
!              sync_dft_cfg (below)

#ifdef VAR_MPI
   use xc_mpi
#endif
   use block_config
#ifdef MOD_XCFUN
   use xcfun
#endif

  implicit none

  public mxbllen
#ifdef VAR_MPI
  public sync_dft_cfg
  public sync_functional
#endif

  public functional
  type functional
#ifdef MOD_XCFUN
     real(8) :: w(XC_NR_PARAMS)
#else
     real(8) :: w(1)
#endif
     integer :: id
     logical :: xcfun_id_is_set
     real(8), allocatable :: input(:, :)
     real(8), allocatable :: output(:, :)
  end type

  character(30), public :: dft_cfg_gridfile

  character(30), public, allocatable :: dft_cfg_multigrid_file(:)
  integer,       public, allocatable :: dft_cfg_multigrid_file_unit(:)
  integer,       public              :: dft_cfg_multigrid_max

  logical, public :: dft_cfg_screening               = .false.
  logical, public :: dft_cfg_keep_ao                 = .false.
  logical, public :: dft_cfg_explicit_potential      = .false.
  logical, public :: dft_cfg_no_pruning              = .false.
  logical, public :: dft_cfg_estimate_radii          = .false.
  logical, public :: dft_cfg_import_grid             = .false.
  logical, public :: dft_cfg_export_grid             = .false.
  logical, public :: dft_cfg_multigrid               = .false.
  logical, public :: dft_cfg_zipgrid                 = .true.
  logical, public :: dft_cfg_alda_hs                 = .false.
  logical, public :: dft_cfg_alda_ha                 = .false.
  logical, public :: dft_cfg_xalda                   = .false.
  logical, public :: dft_cfg_force_4c_grid           = .false.
  logical, public :: dft_cfg_no_sdft                 = .false.
  logical, public :: dft_cfg_sdft_collinear          = .false.
  logical, public :: dft_cfg_cdft_vrg                = .false.
  logical, public :: dft_cfg_xcfun                   = .false.
  logical, public :: dft_cfg_grac                    = .false.
  logical, public :: dft_cfg_grac_is_active          = .false.
  logical, public :: dft_cfg_saop                    = .false.
  logical, public :: dft_cfg_saop_is_active          = .false.
  logical, public :: dft_cfg_saop_with_response_part = .false.
  logical, public :: dft_cfg_asymptote_is_lb94       = .false.
  logical, public :: dft_cfg_asymptote_is_lbalpha    = .false.
  logical, public :: dft_cfg_betasigma               = .false.
  logical, public :: dft_cfg_blocked                 = .false.

  real(8), public :: dft_cfg_radint          = 1.0d-13
  real(8), public :: dft_cfg_accepted_error  = 1.0d-3

!                    for high accuracy high order response this is not tight enough
  real(8), public :: dft_cfg_dfthr0          = 1.0d-9

!                    this threshold is not used in default dft calculations
!                    only when screening is on
!                    currently not read by input
  real(8), public :: dft_cfg_dfthri          = 2.0d-12

  real(8), public :: dft_cfg_ac_ip           = 0.0d0
  real(8), public :: dft_cfg_ac_threshold    = 0.0d0
  real(8), public :: dft_cfg_homo_eigenvalue = 0.0d0
  real(8), public :: dft_cfg_grac_alpha      = 0.0d0
  real(8), public :: dft_cfg_grac_beta       = 0.0d0

  integer, public :: dft_cfg_angmin                  = 15
  integer, public :: dft_cfg_angint                  = 41
  integer, public :: dft_cfg_ao_mem_limit            = 0
  integer, public :: dft_cfg_print_level             = 2
  integer, public :: dft_cfg_debug_level             = 0
  integer, public :: dft_cfg_integration_check_level = 0
  integer, public :: dft_cfg_alda_order              = 99

  type(functional), public :: xc_fun
  type(functional), public :: xc_fun_alda
  type(functional), public :: xc_fun_xalda
  type(functional), public :: kin_fun
  type(functional), public :: kin_fun_alda

  private

#ifdef VAR_MPI
contains

   subroutine sync_functional(f)

!     --------------------------------------------------------------------------
      type(functional) :: f
!     --------------------------------------------------------------------------

      call xc_mpi_bcast(f%w)
      call xc_mpi_bcast(f%id)
#ifdef MOD_XCFUN
      call xc_mpi_bcast(f%xcfun_id_is_set)
#endif

   end subroutine

  subroutine sync_dft_cfg()

    call xc_mpi_bcast(dft_cfg_screening              )
    call xc_mpi_bcast(dft_cfg_keep_ao                )
    call xc_mpi_bcast(dft_cfg_explicit_potential     )
    call xc_mpi_bcast(dft_cfg_no_pruning             )
    call xc_mpi_bcast(dft_cfg_estimate_radii         )
!   only master imports and exports
!   call xc_mpi_bcast(dft_cfg_import_grid            )
!   call xc_mpi_bcast(dft_cfg_export_grid            )
    call xc_mpi_bcast(dft_cfg_multigrid              )
    call xc_mpi_bcast(dft_cfg_zipgrid                )
    call xc_mpi_bcast(dft_cfg_alda_hs                )
    call xc_mpi_bcast(dft_cfg_alda_ha                )
    call xc_mpi_bcast(dft_cfg_xalda                  )
    call xc_mpi_bcast(dft_cfg_force_4c_grid          )
    call xc_mpi_bcast(dft_cfg_no_sdft                )
    call xc_mpi_bcast(dft_cfg_sdft_collinear         )
    call xc_mpi_bcast(dft_cfg_cdft_vrg               )
    call xc_mpi_bcast(dft_cfg_xcfun                  )
    call xc_mpi_bcast(dft_cfg_grac                   )
    call xc_mpi_bcast(dft_cfg_grac_is_active         )
    call xc_mpi_bcast(dft_cfg_saop                   )
    call xc_mpi_bcast(dft_cfg_saop_is_active         )
    call xc_mpi_bcast(dft_cfg_saop_with_response_part)
    call xc_mpi_bcast(dft_cfg_asymptote_is_lb94      )
    call xc_mpi_bcast(dft_cfg_asymptote_is_lbalpha   )
    call xc_mpi_bcast(dft_cfg_betasigma              )
    call xc_mpi_bcast(dft_cfg_blocked                )
    call xc_mpi_bcast(dft_cfg_radint                 )
    call xc_mpi_bcast(dft_cfg_accepted_error         )
    call xc_mpi_bcast(dft_cfg_dfthr0                 )
    call xc_mpi_bcast(dft_cfg_dfthri                 )
    call xc_mpi_bcast(dft_cfg_ac_ip                  )
    call xc_mpi_bcast(dft_cfg_ac_threshold           )
    call xc_mpi_bcast(dft_cfg_homo_eigenvalue        )
    call xc_mpi_bcast(dft_cfg_grac_alpha             )
    call xc_mpi_bcast(dft_cfg_grac_beta              )
    call xc_mpi_bcast(dft_cfg_angmin                 )
    call xc_mpi_bcast(dft_cfg_angint                 )
    call xc_mpi_bcast(dft_cfg_ao_mem_limit           )
    call xc_mpi_bcast(dft_cfg_print_level            )
    call xc_mpi_bcast(dft_cfg_debug_level            )
    call xc_mpi_bcast(dft_cfg_integration_check_level)
    call xc_mpi_bcast(dft_cfg_alda_order             )

  end subroutine
#endif /* ifdef VAR_MPI */

end module
