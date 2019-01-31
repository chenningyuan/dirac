module exacorr_mo 


! Contains all the necessary information for a range of molecular orbitals
! In particular:  - Orbital energy
!                 - Boson characterization (for spinfree runs)
!                 - Coefficients of AO-functions
!                 - Position in the full list of MO vectors


!!!!!!!!!!!!!!!!!LOSE DEPENDENCY ON DIRAC_MO -> PUT EVERYTHING IN HERE
!!!!!!!!!!!!!!!!!HAVE BASIC MO TYPE, EXTEND WITH EITHER QMO OR CMO
!!!!!!!!!!!!!!!!!FOR NOW ONLY CMO


    implicit none

    private

    public :: cmo, create_cmo, delete_cmo, sync_cmo

    type cmo

        integer             :: &
         fermion_irrep, &     ! fermion irrep (1=gerade, 2=ungerade) to which the spinors belongs
         nao, &               ! number of aos
         nmo                  ! number of mos

        real(8)             :: &
         total_energy         ! total (SCF) energy as stored on the coefficient file

        real(8), pointer    :: &
         energy (:)           ! orbital energies


        complex(8), pointer :: &
         coeff (:,:)          ! coefficients (AO,MO)
         
        integer, pointer    :: &
         index(:), &          ! index of the spinor in the full list of spinors for this fermion symmetry
         boson_irrep(:)       ! for spinfree calculations: boson irrep to which the spinor belongs

    end type cmo

    interface create_cmo
      module procedure alloc_cmo
    end interface

    interface delete_cmo
      module procedure dealloc_cmo
    end interface

    interface sync_cmo
      module procedure sync_cmo
    end interface

    contains

     subroutine alloc_cmo (M, nao, nmo)

!     Prepare to work with a block of Mos

      type(cmo), intent(out) :: M
      integer,   intent(in ) :: nao
      integer,   intent(in ) :: nmo

      M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
      M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)

!     Reserve the memory needed to work with a block of Mos
      nullify  (M%energy)
      allocate (M%energy(M%nmo))
      M%energy = 0.0

      nullify  (M%coeff)
      allocate (M%coeff(M%nao,M%nmo))
      M%coeff = 0.0

      nullify  (M%boson_irrep)
      allocate (M%boson_irrep(M%nmo))
      M%boson_irrep = 0

      nullify  (M%index)
      allocate (M%index(M%nmo))
      M%index = 0

     end subroutine alloc_cmo

     subroutine dealloc_cmo (M)

!     Clean up workspace

      type(cmo), intent(inout) :: M

!     Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
      M%nao = 0
      M%nmo = 0

!     Free the memory associated with a block of Mos
      deallocate (M%energy)
      deallocate (M%coeff)
      deallocate (M%boson_irrep)
      deallocate (M%index)

     end subroutine dealloc_cmo

#if defined (VAR_MPI)

     subroutine sync_cmo (M,my_MPI_master)

      use interface_to_mpi

!     Synchronize and allocate (if not on master) a block of MOcoefficients

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master
      integer                  :: my_MPI_rank
      integer                  :: n_coef

      call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
      call interface_mpi_bcast (M%nao,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%nmo,1,my_MPI_master,global_communicator)
      if (my_MPI_rank /= my_MPI_master) then
         call alloc_cmo (M,M%nao,M%nmo)
      end if
      call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%boson_irrep,  M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
      n_coef = 2*M%nao * 2*M%nmo
      call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)

     end subroutine sync_cmo

#else

     subroutine sync_cmo (M,my_MPI_master)

!     do nothing, makes it possible to call sync also in serial runs

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master

     end subroutine sync_cmo

#endif

end module exacorr_mo
