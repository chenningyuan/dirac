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

module dirac_MO

!
! Contains all the necessary information for a range of molecular orbitals
! In particular:  - Orbital energy
!                 - Boson characterization (for spinfree runs)
!                 - Coefficients of AO-functions
!                 - Position in the full list of MO vectors

   implicit none
  
!  Private

  Public :: create_qmo, delete_qmo, sync_qmo

  type MO
     
    integer ::  &
      fermion_irrep, &     ! fermion irrep (1=gerade, 2=ungerade) to which the spinors belongs
      nao,      &          ! number of aos
      nmo,      &          ! number of mos
      nz                   ! number of active quaternion units

    real(8) :: &
      total_energy         ! total (SCF) energy as stored on the coefficient file

    real(8), pointer ::  &
      energy    (:),     & ! orbital energies
      coeff (:,:,:)        ! coefficients (AO,MO,IZ)
     
    integer, pointer ::  &
      index(:),          & ! index of the spinor in the full list of spinors for this fermion symmetry
      boson_irrep(:)       ! for spinfree calculations: boson irrep to which the spinor belongs

  end type MO

  interface create_qmo
    module procedure AllocMO
  end interface  

  interface delete_qmo
    module procedure DeAllocMO
  end interface  
  
  interface sync_qmo
    module procedure SyncMO
  end interface  

contains

  subroutine AllocMO (M, nao, nmo, nz)
   
!   Prepare to work with a block of Mos

    type (MO), intent (out) :: M
    integer, intent(in)           :: nao
    integer, intent(in)           :: nmo
    integer, intent(in), optional :: nz

    M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
    M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)
    if (present(nz)) then
       M%nz = nz
    else
       M%nz = 1  ! We default to one, easy for spinfree applications
    endif

!   Reserve the memory needed to work with a block of Mos
    nullify  (M%energy)
    allocate (M%energy(M%nmo))
    M%energy = 0.0

    nullify  (M%coeff)
    allocate (M%coeff(M%nao,M%nmo,M%nz))
    M%coeff = 0.0

    nullify  (M%boson_irrep)
    allocate (M%boson_irrep(M%nmo))
    M%boson_irrep = 0

    nullify  (M%index)
    allocate (M%index(M%nmo))
    M%index = 0
    
  end subroutine AllocMO

  subroutine DeAllocMO (M)

!   Clean up workspace

    type (MO), intent (inout) :: M

!   Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
    M%nao = 0
    M%nmo = 0
    M%nz = 0

!   Free the memory associated with a block of Mos
    deallocate (M%energy)
    deallocate (M%coeff)
    deallocate (M%boson_irrep)
    deallocate (M%index)
    
  end subroutine DeAllocMO

#if defined (VAR_MPI)

  subroutine SyncMO (M,my_MPI_master)
   
    use interface_to_mpi

!   Synchronize and allocate (if not on master) a block of MOcoefficients
 
    type (MO), intent (inout)   :: M
    integer, intent(in)         :: my_MPI_master
    integer                     :: my_MPI_rank
    integer                     :: n_coef

    call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
    call interface_mpi_bcast (M%nao,1,my_MPI_master,global_communicator)
    call interface_mpi_bcast (M%nmo,1,my_MPI_master,global_communicator)
    call interface_mpi_bcast (M%nz ,1,my_MPI_master,global_communicator)
    if (my_MPI_rank /= my_MPI_master) then
       call AllocMO (M,M%nao,M%nmo,M%nz)
    end if
    call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
    call interface_mpi_bcast (M%boson_irrep,  M%nmo,my_MPI_master,global_communicator)
    call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
    call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
    n_coef = M%nao * M%nmo * M%nz
    call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)

  end subroutine SyncMO

#else

  subroutine SyncMO (M,my_MPI_master)
 
!   do nothing, makes it possible to call sync also in serial runs

    type (MO), intent (inout)   :: M
    integer, intent(in)         :: my_MPI_master

  end subroutine SyncMO

#endif

end module dirac_MO
