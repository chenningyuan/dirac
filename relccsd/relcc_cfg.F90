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

module relcc_cfg

! First step in removing all common blocks in the RELCCSD programme.
! This one contains all input set by the user

  implicit none
!Pawel new relcc input 1.03.2011
  integer, dimension(1:16), public :: relcc_nelec   = 0  
  integer, dimension(1:16), public :: relcc_nelec_f1= 0 
  integer, dimension(1:16), public :: relcc_nelec_f2= 0 
  integer, dimension(1:16), public :: relcc_nfroz   = 0
  save

! General control flags
  logical, public :: nelec_input                   = .false.
  logical, public :: relcc_debug                   = .false.
  logical, public :: relcc_timing                  = .false.
  logical, public :: relcc_do_sort                 = .true.
  logical, public :: relcc_do_energy               = .false.
  logical, public :: relcc_do_gradient             = .false.
  logical, public :: relcc_do_fspc                 = .false.
  logical, public :: relcc_carith                  = .false.
  integer, public :: relcc_print                   = 0
  integer, public :: relcc_max_dimension_diis      = 8
  integer, public :: relcc_max_iterations          = 30
  integer, public :: relcc_desired_convergence     = 12
 
! Information relevant for the integral sorting and setup stage
  logical, public :: relcc_breit                    = .false.
  logical, public :: relcc_add_finite_field         = .false.
  logical, public :: relcc_use_orbital_energies     = .true.
  logical, public :: relcc_no_recompute             = .false.
  logical, public :: relcc_low_memory               = .false.
  integer, public :: relcc_memory_in_mw             = 0
  character(10), public :: relcc_fs_interface
  character( 8), public :: relcc_ff_property_names

! Information relevant for the energy calculation
  logical, public :: relcc_do_mp2                  = .true.
  logical, public :: relcc_do_singles              = .false.
  logical, public :: relcc_do_doubles              = .false.
  logical, public :: relcc_do_ccsd                 = .true.
  logical, public :: relcc_do_ccsd_t               = .true.
  integer, public :: relcc_ccener_max_dimension_diis      = 8
  integer, public :: relcc_ccener_max_iterations          = 30
  integer, public :: relcc_ccener_desired_convergence     = 12

! Information relevant for the gradient calculation
  logical, public :: relcc_do_mp2gradient          = .false.
  logical, public :: relcc_do_oldmp2gradient       = .false.
  logical, public :: relcc_do_ccsdgradient         = .false.
  logical, public :: relcc_do_ccsdtgradient        = .false.
  logical, public :: relcc_do_naturalorbitals      = .false.
  logical, public :: relcc_do_zipgrid              = .false.
  integer, public :: relcc_ne_oper                 = 1
  integer, public :: relcc_fopr_max_dimension_diis = 8
  integer, public :: relcc_fopr_max_iterations     = 30
  integer, public :: relcc_fopr_desired_convergence= 12

! Information relevant for the (Intermediate Hamiltonian) Fock space calculations
  logical, public :: relcc_fs_do_ea                   = .false.
  logical, public :: relcc_fs_do_ea2                  = .false.
  logical, public :: relcc_fs_do_exc                  = .false.
  logical, public :: relcc_fs_do_ie                   = .false.
  logical, public :: relcc_fs_do_ie2                  = .false.
  logical, public :: relcc_fs_do_ih                   = .false.
  integer, public, dimension(1:16) :: relcc_fs_nacth  = 0
  integer, public, dimension(1:16) :: relcc_fs_nactp  = 0
  integer, public, dimension(1:6)  :: relcc_fs_fssect = (/1,0,0,0,0,0/)
  integer, public :: relcc_fs_max_dimension_diis      = 8
  integer, public :: relcc_fs_max_iterations          = 30
  integer, public :: relcc_fs_desired_convergence     = 12
  integer, public :: relcc_fs_max02_iterations        = -1
  integer, public :: relcc_fs_max01_iterations        = -1
  integer, public :: relcc_fs_max20_iterations        = -1
  integer, public :: relcc_fs_max10_iterations        = -1
  integer, public :: relcc_fs_max11_iterations        = -1
  integer, public :: relcc_fs_max00_iterations        = -1
  integer, public :: relcc_fs_tshold                  = 0
! Intermediate Hamiltonian  
  integer, public, dimension(1:16) :: relcc_ih_nacthi = 0
  integer, public, dimension(1:16) :: relcc_ih_nactpi = 0
  logical, public :: relcc_ih_do_zipgrid              = .false.
  integer, public :: relcc_ih_scheme                  = 2
  integer, public :: relcc_ih_nih                     = 1
  real(8), public :: relcc_ih_shift_h11               = 0.0d+00
  real(8), public :: relcc_ih_shift_h12               = 0.0d+00
  real(8), public :: relcc_ih_shift_p11               = 0.0d+00
  real(8), public :: relcc_ih_shift_p12               = 0.0d+00
  real(8), public :: relcc_ih_shift_h2                = 0.0d+00
  real(8), public :: relcc_ih_shift_p2                = 0.0d+00
  real(8), public :: relcc_ih_shift_hp                = 0.0d+00
  real(8), public :: relcc_ih_aih                     = 0.0d+00
  real(8), public :: relcc_ih_eh_min                  = 1000.0d+00
  real(8), public :: relcc_ih_eh_max                  = 100.0d+00
  real(8), public :: relcc_ih_ep_min                  = -1d+00
  real(8), public :: relcc_ih_ep_max                  = -1000.0d+00

end module
