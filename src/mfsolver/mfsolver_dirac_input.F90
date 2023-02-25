!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module mfsolver_dirac_input

  implicit none

#ifdef VAR_MPI
  public mfsolver_dirac_input_sync_cw
#endif

! Information relevant for matrix-free diagonalization (Davidson etc)
  real*8, public  :: mfsolver_options_convergence_threshold = 1.0d-8
  real*8, public  :: mfsolver_options_eigenvalues_energy_shift(2) = (/0.0d0, 0.0d0/)
  integer, public :: mfsolver_options_max_subspace_size = 128
  integer, public :: mfsolver_options_max_iterations = 80
  integer, public :: mfsolver_options_refresh_rate = 0
  logical, public :: mfsolver_options_overlap_sorting = .false.
  logical, public :: mfsolver_options_verbose = .false.
  logical, public :: mfsolver_options_trial_diagonal = .true.
  logical, public :: mfsolver_options_trial_ccs = .false. 
  logical, public :: mfsolver_options_trial_full_matrix = .false.
  logical, public :: mfsolver_options_trial_restart = .false.
  logical, public :: mfsolver_options_trial_lhs_use_rhs = .true.
  logical, public :: mfsolver_options_trial_double_mfold = .false.

  logical, public :: mfsolver_projectors_do_restricted_excitation_window = .false. 
  logical, public :: mfsolver_projectors_do_core_valence_separation = .false. 
  logical, public :: mfsolver_projectors_rew_strict = .false.
  logical, public :: mfsolver_projectors_rew_remove_double_occupied = .false.
  logical, public :: mfsolver_projectors_frozen_core = .false.
  real*8, public  :: mfsolver_projectors_rew_occ_min_energy = -huge(mfsolver_projectors_rew_occ_min_energy)
  real*8, public  :: mfsolver_projectors_rew_occ_max_energy =  huge(mfsolver_projectors_rew_occ_max_energy)
  real*8, public  :: mfsolver_projectors_rew_virt_min_energy = -huge(mfsolver_projectors_rew_virt_min_energy)
  real*8, public  :: mfsolver_projectors_rew_virt_max_energy =  huge(mfsolver_projectors_rew_virt_max_energy)
  real*8, public  :: mfsolver_projectors_frozen_core_max_energy =  huge(mfsolver_projectors_frozen_core_max_energy)

#ifdef VAR_MPI

contains
  subroutine mfsolver_dirac_input_sync_cw()
!
! please add bcast here if a new variable has been added to the list above
!
  use interface_to_mpi

! Information relevant to the davidson or other type of matrix-free diagonalization methods

  call interface_mpi_bcast(mfsolver_options_convergence_threshold ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_options_eigenvalues_energy_shift,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_options_max_subspace_size     ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_options_max_iterations        ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_options_refresh_rate          ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_overlap_sorting    ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_verbose            ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_trial_diagonal     ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_trial_ccs          ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_trial_full_matrix  ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_trial_restart      ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_trial_lhs_use_rhs  ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_options_trial_double_mfold ,1,0,global_communicator)

! information relevant to projectors that modify the matrix-free trial vectores

  call interface_mpi_bcast_l0(mfsolver_projectors_do_restricted_excitation_window,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_projectors_do_core_valence_separation,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_projectors_rew_strict                     ,1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_projectors_rew_remove_double_occupied     , 1,0,global_communicator)
  call interface_mpi_bcast_l0(mfsolver_projectors_frozen_core     , 1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_projectors_rew_occ_min_energy    ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_projectors_rew_occ_max_energy    ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_projectors_rew_virt_min_energy   ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_projectors_rew_virt_max_energy   ,1,0,global_communicator)
  call interface_mpi_bcast(mfsolver_projectors_frozen_core_max_energy    ,1,0,global_communicator)

  end subroutine mfsolver_dirac_input_sync_cw
#endif

end module
