module dft_cfg

  use dft_cb_interface

! radovan: (1) if you modify this module (add new stuff)
!              don't forget to update
!              routines init_dft_cfg and sync_dft_cfg (below)
!          (2) this module should contain things set by the user
!              or in principle accessible by the user
!              temporary parameters used by the program
!              (like integrated electron count etc.)
!              should not be here, put them in dft_temp

  implicit none

  save

  character(30) :: dft_cfg_gridfile

  logical       :: dft_cfg_set
  logical       :: dft_cfg_is_dft_calculation

  logical       :: dft_cfg_keep_ao
  logical       :: dft_cfg_explicit_potential
  logical       :: dft_cfg_no_pruning
  logical       :: dft_cfg_estimate_radii
  logical       :: dft_cfg_import_grid
  logical       :: dft_cfg_export_grid
  logical       :: dft_cfg_zipgrid
  logical       :: dft_cfg_alda_hs
  logical       :: dft_cfg_alda_ha
  logical       :: dft_cfg_xalda
  logical       :: dft_cfg_force_4c_grid
  logical       :: dft_cfg_old_dft
  logical       :: dft_cfg_sdft
  logical       :: dft_cfg_sdft_collinear
  logical       :: dft_cfg_sdft_beta
  logical       :: dft_cfg_cdft_vrg
  logical       :: dft_cfg_autodiff
  logical       :: dft_cfg_grac
  logical       :: dft_cfg_grac_is_active
  logical       :: dft_cfg_saop
  logical       :: dft_cfg_saop_is_active
  logical       :: dft_cfg_saop_with_response_part
  logical       :: dft_cfg_asymptote_is_lb94
  logical       :: dft_cfg_asymptote_is_lbalpha

  real(8)       :: dft_cfg_radint
  real(8)       :: dft_cfg_accepted_error
  real(8)       :: dft_cfg_dfthr0
  real(8)       :: dft_cfg_dfthrl
  real(8)       :: dft_cfg_dfthri
  real(8)       :: dft_cfg_ac_ip
  real(8)       :: dft_cfg_ac_threshold
  real(8)       :: dft_cfg_homo_eigenvalue
  real(8)       :: dft_cfg_grac_alpha
  real(8)       :: dft_cfg_grac_beta

  integer       :: dft_cfg_angmin
  integer       :: dft_cfg_angint
  integer       :: dft_cfg_ao_mem_limit
  integer       :: dft_cfg_print_level
  integer       :: dft_cfg_debug_level
  integer       :: dft_cfg_integration_check_level

contains

  subroutine init_dft_cfg()

    dft_cfg_keep_ao                 = .false.
    dft_cfg_explicit_potential      = .false.
    dft_cfg_no_pruning              = .false.
    dft_cfg_estimate_radii          = .false.
    dft_cfg_import_grid             = .false.
    dft_cfg_export_grid             = .false.
    dft_cfg_zipgrid                 = .false.
    dft_cfg_alda_hs                 = .false.
    dft_cfg_alda_ha                 = .false.
    dft_cfg_xalda                   = .false.
    dft_cfg_force_4c_grid           = .false.
    dft_cfg_old_dft                 = .false.
    dft_cfg_sdft                    = .false.
    dft_cfg_sdft_collinear          = .false.
    dft_cfg_sdft_beta               = .false.
    dft_cfg_cdft_vrg                = .false.
    dft_cfg_autodiff                = .false.
    dft_cfg_grac                    = .false.
    dft_cfg_grac_is_active          = .false.
    dft_cfg_saop                    = .false.
    dft_cfg_saop_is_active          = .false.
    dft_cfg_saop_with_response_part = .false.
    dft_cfg_asymptote_is_lb94       = .false.
    dft_cfg_asymptote_is_lbalpha    = .false.

    dft_cfg_radint                  = 1.0d-13
    dft_cfg_accepted_error          = 1.0d-3
    dft_cfg_dfthr0                  = 1.0d-9
    dft_cfg_dfthrl                  = 1.0d-10
    dft_cfg_dfthri                  = 1.0d-15
    dft_cfg_ac_ip                   = 0.0d0
    dft_cfg_ac_threshold            = 0.0d0
    dft_cfg_homo_eigenvalue         = 0.0d0
    dft_cfg_grac_alpha              = 0.0d0
    dft_cfg_grac_beta               = 0.0d0

    dft_cfg_angmin                  = 15
    dft_cfg_angint                  = 41
!radovan: don't reset here, set with pam
!   dft_cfg_ao_mem_limit            = 0
    dft_cfg_print_level             = 2
    dft_cfg_debug_level             = 0
    dft_cfg_integration_check_level = 0

#ifdef MOD_AUTODIFF
    call init_factorial_common_block()
#endif
    call init_dft_cb()

  end subroutine

#ifdef VAR_MPI
  subroutine sync_dft_cfg()

#include "mpif.h"
#include "infpar.h"

!   ----------------------------------------------------------------------------
    logical             :: i_am_the_master
    integer             :: idx, ierr
    logical             :: pack_logical(24)
    real(8)             :: pack_double_precision(10)
    integer             :: pack_integer(6)
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

      pack_logical( 1) = dft_cfg_keep_ao
      pack_logical( 2) = dft_cfg_explicit_potential
      pack_logical( 3) = dft_cfg_no_pruning
      pack_logical( 4) = dft_cfg_estimate_radii
      pack_logical( 5) = dft_cfg_import_grid
      pack_logical( 6) = dft_cfg_export_grid
      pack_logical( 7) = dft_cfg_zipgrid
      pack_logical( 8) = dft_cfg_alda_hs
      pack_logical( 9) = dft_cfg_alda_ha
      pack_logical(10) = dft_cfg_xalda
      pack_logical(11) = dft_cfg_force_4c_grid
      pack_logical(12) = dft_cfg_old_dft
      pack_logical(13) = dft_cfg_sdft
      pack_logical(14) = dft_cfg_sdft_collinear
      pack_logical(15) = dft_cfg_sdft_beta
      pack_logical(16) = dft_cfg_cdft_vrg
      pack_logical(17) = dft_cfg_autodiff
      pack_logical(18) = dft_cfg_grac
      pack_logical(19) = dft_cfg_grac_is_active
      pack_logical(20) = dft_cfg_saop
      pack_logical(21) = dft_cfg_saop_is_active
      pack_logical(22) = dft_cfg_saop_with_response_part
      pack_logical(23) = dft_cfg_asymptote_is_lb94
      pack_logical(24) = dft_cfg_asymptote_is_lbalpha

      pack_double_precision( 1) = dft_cfg_radint
      pack_double_precision( 2) = dft_cfg_accepted_error
      pack_double_precision( 3) = dft_cfg_dfthr0
      pack_double_precision( 4) = dft_cfg_dfthrl
      pack_double_precision( 5) = dft_cfg_dfthri
      pack_double_precision( 6) = dft_cfg_ac_ip
      pack_double_precision( 7) = dft_cfg_ac_threshold
      pack_double_precision( 8) = dft_cfg_homo_eigenvalue
      pack_double_precision( 9) = dft_cfg_grac_alpha
      pack_double_precision(10) = dft_cfg_grac_beta

      pack_integer(1) = dft_cfg_angmin
      pack_integer(2) = dft_cfg_angint
      pack_integer(3) = dft_cfg_ao_mem_limit
      pack_integer(4) = dft_cfg_print_level
      pack_integer(5) = dft_cfg_debug_level
      pack_integer(6) = dft_cfg_integration_check_level

      idx = 1
      call mpi_pack(pack_logical,       &
                    size(pack_logical), &
                    mpi_logical,        &
                    scratch,            &
                    nbytes,             &
                    idx,                &
                    mpi_comm_world,     &
                    ierr)
      call mpi_pack(pack_double_precision,       &
                    size(pack_double_precision), &
                    mpi_double_precision,        &
                    scratch,                     &
                    nbytes,                      &
                    idx,                         &
                    mpi_comm_world,              &
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
      call mpi_unpack(scratch,                     &
                      nbytes,                      &
                      idx,                         &
                      pack_double_precision,       &
                      size(pack_double_precision), &
                      mpi_double_precision,        &
                      mpi_comm_world,              &
                      ierr)
      call mpi_unpack(scratch,            &
                      nbytes,             &
                      idx,                &
                      pack_integer,       &
                      size(pack_integer), &
                      my_MPI_INTEGER,     &
                      mpi_comm_world,     &
                      ierr)

      dft_cfg_keep_ao                 = pack_logical( 1)
      dft_cfg_explicit_potential      = pack_logical( 2)
      dft_cfg_no_pruning              = pack_logical( 3)
      dft_cfg_estimate_radii          = pack_logical( 4)
      dft_cfg_import_grid             = pack_logical( 5)
      dft_cfg_export_grid             = pack_logical( 6)
      dft_cfg_zipgrid                 = pack_logical( 7)
      dft_cfg_alda_hs                 = pack_logical( 8)
      dft_cfg_alda_ha                 = pack_logical( 9)
      dft_cfg_xalda                   = pack_logical(10)
      dft_cfg_force_4c_grid           = pack_logical(11)
      dft_cfg_old_dft                 = pack_logical(12)
      dft_cfg_sdft                    = pack_logical(13)
      dft_cfg_sdft_collinear          = pack_logical(14)
      dft_cfg_sdft_beta               = pack_logical(15)
      dft_cfg_cdft_vrg                = pack_logical(16)
      dft_cfg_autodiff                = pack_logical(17)
      dft_cfg_grac                    = pack_logical(18)
      dft_cfg_grac_is_active          = pack_logical(19)
      dft_cfg_saop                    = pack_logical(20)
      dft_cfg_saop_is_active          = pack_logical(21)
      dft_cfg_saop_with_response_part = pack_logical(22)
      dft_cfg_asymptote_is_lb94       = pack_logical(23)
      dft_cfg_asymptote_is_lbalpha    = pack_logical(24)

      dft_cfg_radint          = pack_double_precision( 1)
      dft_cfg_accepted_error  = pack_double_precision( 2)
      dft_cfg_dfthr0          = pack_double_precision( 3)
      dft_cfg_dfthrl          = pack_double_precision( 4)
      dft_cfg_dfthri          = pack_double_precision( 5)
      dft_cfg_ac_ip           = pack_double_precision( 6)
      dft_cfg_ac_threshold    = pack_double_precision( 7)
      dft_cfg_homo_eigenvalue = pack_double_precision( 8)
      dft_cfg_grac_alpha      = pack_double_precision( 9)
      dft_cfg_grac_beta       = pack_double_precision(10)

      dft_cfg_angmin                  = pack_integer(1)
      dft_cfg_angint                  = pack_integer(2)
      dft_cfg_ao_mem_limit            = pack_integer(3)
      dft_cfg_print_level             = pack_integer(4)
      dft_cfg_debug_level             = pack_integer(5)
      dft_cfg_integration_check_level = pack_integer(6)

    end if

  end subroutine
#endif /* ifdef VAR_MPI */

end module
