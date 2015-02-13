!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

module relcc_cfg

! First step in removing all common blocks in the RELCCSD programme.
! This one contains all input set by the user

  implicit none

! stefan: sync co-workers with new input-style variables
#ifdef VAR_MPI
  public relcc_sync_cw
#endif

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
  logical, public :: relcc_do_energy               = .true.
  logical, public :: relcc_do_gradient             = .false.
  logical, public :: relcc_do_hessian              = .false. ! this currently only relates to the old Hartree-Fock MO response code 
  logical, public :: relcc_do_fspc                 = .false.
  logical, public :: relcc_carith                  = .false.
  integer, public :: relcc_print                   = 0
  integer, public :: relcc_max_dimension_diis      = 8
  integer, public :: relcc_max_iterations          = 30
  integer, public :: relcc_desired_convergence     = 12
 
! Information relevant for the integral sorting and setup stage
  logical, public :: relcc_breit                    = .false.
  logical, public :: relcc_add_finite_field         = .false.
  logical, public :: relcc_use_orbital_energies     = .false.
  logical, public :: relcc_no_recompute             = .false.
  logical, public :: relcc_low_memory               = .false.
  integer, public :: relcc_memory_in_mw             = 0
  character(10), public :: relcc_integral_interface = 'DIRAC6    '
  character( 8), public :: relcc_ff_property_names

! Information relevant for the energy calculation
  logical, public :: relcc_do_mp2                  = .true. 
  logical, public :: relcc_no_singles              = .false.
  logical, public :: relcc_no_doubles              = .false.
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
  logical, public :: relcc_do_unrelaxed            = .false.
  logical, public :: relcc_ifroz_input             = .false.
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
  integer, public :: relcc_fs_select_state_for_numgrad_energy = -1
! Intermediate Hamiltonian  
  integer, public, dimension(1:16) :: relcc_ih_nacthi = 0
  integer, public, dimension(1:16) :: relcc_ih_nactpi = 0
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

#ifdef VAR_MPI

contains
  subroutine relcc_sync_cw()
!
! stefan: please add bcast here if a new variable has been added to the list
!         above
!
  use interface_to_mpi

#include "freeze.inc"
  
  call interface_mpi_BCAST(relcc_nelec,   size(relcc_nelec),   0,global_communicator)
  call interface_mpi_BCAST(relcc_nelec_f1,size(relcc_nelec_f1),0,global_communicator)
  call interface_mpi_BCAST(relcc_nelec_f2,size(relcc_nelec_f2),0,global_communicator)
  call interface_mpi_BCAST(relcc_nfroz,   size(relcc_nfroz),   0,global_communicator)

  call interface_mpi_bcast_l0(nelec_input              ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_debug              ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_timing             ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_sort            ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_energy          ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_gradient        ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_hessian         ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_fspc            ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_carith             ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_print              ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_max_dimension_diis ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_max_iterations     ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_desired_convergence,1,0,global_communicator)

  call interface_mpi_bcast_l0(relcc_breit               ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_add_finite_field    ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_use_orbital_energies,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_no_recompute        ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_low_memory          ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_memory_in_mw        ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_integral_interface, 10,0,global_communicator)
  call interface_mpi_bcast(relcc_ff_property_names,   8,0,global_communicator)

! Information relevant for the energy calculation
  call interface_mpi_bcast_l0(relcc_do_mp2                    ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_no_singles                ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_no_doubles                ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_ccsd                   ,1,0,global_communicator) 
  call interface_mpi_bcast_l0(relcc_do_ccsd_t                 ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ccener_max_dimension_diis ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ccener_max_iterations     ,1,0,global_communicator) 
  call interface_mpi_bcast(relcc_ccener_desired_convergence,1,0,global_communicator)

! Information relevant for the gradient calculation
  call interface_mpi_bcast_l0(relcc_do_mp2gradient       ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_oldmp2gradient    ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_ccsdgradient      ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_ccsdtgradient     ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_naturalorbitals   ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_do_unrelaxed         ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_ifroz_input          ,1,0,global_communicator)
  if (relcc_ifroz_input) then  ! transfer IFROZ array in freeze.inc
     IFROZ_set_in_input = .true.
     call interface_mpi_bcast(IFROZ        ,MAXFROZ*MAXIREP,0,global_communicator)
  end if
  call interface_mpi_bcast(relcc_ne_oper                 ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fopr_max_dimension_diis ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fopr_max_iterations     ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fopr_desired_convergence,1,0,global_communicator)

! Information relevant for the (Intermediate Hamiltonian) Fock space calculations
  call interface_mpi_bcast_l0(relcc_fs_do_ea             ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_fs_do_ea2            ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_fs_do_exc            ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_fs_do_ie             ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_fs_do_ie2            ,1,0,global_communicator)
  call interface_mpi_bcast_l0(relcc_fs_do_ih             ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_nacth, size(relcc_fs_nacth),   0,global_communicator)
  call interface_mpi_bcast(relcc_fs_nactp, size(relcc_fs_nactp),   0,global_communicator)
  call interface_mpi_bcast(relcc_fs_fssect,size(relcc_fs_fssect),  0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max_dimension_diis         ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max_iterations             ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_desired_convergence        ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max02_iterations           ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max01_iterations           ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max20_iterations           ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max10_iterations           ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max11_iterations           ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_max00_iterations           ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_tshold                     ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_fs_select_state_for_numgrad_energy,1,0,global_communicator)
! Intermediate Hamiltonian  
  call interface_mpi_bcast(relcc_ih_nacthi,size(relcc_ih_nacthi),0,global_communicator)
  call interface_mpi_bcast(relcc_ih_nactpi,size(relcc_ih_nactpi),0,global_communicator)
  call interface_mpi_bcast(relcc_ih_scheme                  ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_nih                     ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_h11               ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_h12               ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_p11               ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_p12               ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_h2                ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_p2                ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_shift_hp                ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_aih                     ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_eh_min                  ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_eh_max                  ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_ep_min                  ,1,0,global_communicator)
  call interface_mpi_bcast(relcc_ih_ep_max                  ,1,0,global_communicator)

  end subroutine relcc_sync_cw
#endif

end module
