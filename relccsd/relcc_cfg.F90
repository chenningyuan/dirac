module relcc_cfg

! First step in removing all common blocks in the RELCCSD programme.
! This one contains all input set by the user

  implicit none

  save

! General control flags
  logical, public :: relcc_debug                   = .false.
  logical, public :: relcc_timing                  = .false.
  logical, public :: relcc_do_sort                 = .true.
  logical, public :: relcc_do_energy               = .false.
  logical, public :: relcc_do_gradient             = .false.
  logical, public :: relcc_do_fspc                 = .false.
  integer, public :: relcc_print                   = 0
 
! Information relevant for the integral sorting and setup stage
  logical, public :: relcc_breit                   = .false.
  logical, public :: relcc_add_finite_field        = .false.
  logical, public :: relcc_use_orbital_energies    = .true.
  logical, public :: relcc_no_recompute            = .false.
  logical, public :: relcc_low_memory              = .false.
  integer, public :: relcc_nelec(16)               = 0
  integer, public :: relcc_nelec_f1(16)            = 0
  integer, public :: relcc_nelec_f2(16)            = 0
  integer, public :: relcc_memory_in_mw            = 0
  character(10), public :: relcc_integral_interface
  character( 8), public :: relcc_ff_property_names

! Information relevant for the energy calculation
  logical, public :: relcc_do_mp2                  = .true.
  logical, public :: relcc_do_ccsd                 = .true.
  logical, public :: relcc_do_ccsd_t               = .true.
  logical, public :: relcc_no_zg                   = .false.
  integer, public :: relcc_max_dimension_diis      = 0
  integer, public :: relcc_max_iterations          = 0
  real(8), public :: relcc_desired_convergence     = 1.0d-13

! Information relevant for the gradient calculation
  logical, public :: relcc_do_mp2gradient          = .false.
  logical, public :: relcc_do_oldmp2gradient       = .false.
  logical, public :: relcc_do_ccsdgradient         = .false.
  logical, public :: relcc_do_naturalorbitals      = .false.

! Information relevant for the (Intermediate Hamiltonian) Fock space calculations
  logical, public :: relcc_fs_do_IntHamiltonian    = .false.
  integer, public :: relcc_fs_sectors(6)           = 0
  integer, public :: relcc_fs_n_inactive(16)       = 0
  integer, public :: relcc_fs_n_active(16)         = 0
  integer, public :: relcc_max_iterations_fs(6)    = 0
  real(8), public :: relcc_desired_convergence_in_fs       = 1.0d-13
  real(8), public :: relcc_ff_property_strength(2,10)      = 0.0d+00

  private

end module
