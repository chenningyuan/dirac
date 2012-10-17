!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!
!
! this module contains all the functionality required to add one- and two-electron
! spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis. 
! the spin-orbit corrections will be obtained as the spin-dependent parts of the 
! (converged) 4c-Fock operator in an atomic Dirac run.
!
! written by sknecht june 2011
!
module x2c_dirac_in_house_so_correction

  use memory_allocator 
  use x2c_utils

  implicit none

  public prepare_so_correction_driver

  private

contains

!**********************************************************************
  subroutine prepare_so_correction_driver(nr_fsym,    &
                                          nr_quat,    &
                                          norb_dim,   &
                                          print_lvl,  &
                                          xmat)
!**********************************************************************
     real(8), optional, intent(in)      :: xmat(*)
     integer, intent(in)                :: nr_fsym
     integer, intent(in)                :: nr_quat
     integer, intent(in)                :: norb_dim(nr_fsym)
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     integer                            :: i
     integer                            :: j
     integer                            :: l
     integer                            :: offset_boson_info
     integer                            :: offset_mo_orbitals
     integer                            :: file_unit
     integer, allocatable               :: boson_info(:)
     integer, allocatable               :: scratch1(:)
     real(8), allocatable               :: molecular_orbitals_onmo(:)
     real(8), allocatable               :: scratch2(:)
     real(8), allocatable               :: scratch3(:)
     real(8), allocatable               :: scratch4(:)
     integer                            :: nzsave
     logical                            :: dojaco_save
     logical                            :: sbsave
!**********************************************************************

!     out: 4c-fock operator in SA-AO basis -> h1_saao_bas
!     call get_4cdf_saao_basis(h1_saao_bas,        &
!                              scrmat1,            &
!                              scrmat2,            &
!                              naosh_ls,           &
!                              naosh_all,          &
!                              norb_dim,           &
!                              nfsym,              &
!                              nz,                 &
!                              ipq_off_in,         &
!                              linear_sym,         &
!                              lwrk,               &
!                              x2c_file_scr,       &
!                              print_lvl)

!     call x2c_export_mat2C1


  end subroutine prepare_so_correction_driver

!**********************************************************************

end module x2c_dirac_in_house_so_correction
