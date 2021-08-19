module exacorr_mo 

!   Contains all the necessary information for a range of molecular orbitals
!   In particular:  - Orbital energy
!                   - Symmetry characterization
!                   - Coefficients of AO-functions
!                   - Position in the full list of MO vectors

!   Extended from dirac_mo by Stan Papadopoulos, January 2019

    implicit none

    private

    public :: alloc_mo, dealloc_mo, sync_mo, compress_mo, copy_mo
    public :: convert_mo_to_quaternion,convert_mo_to_complex,apply_Kramers_operator

    type, public  :: cmo
        integer             :: &
         nao, &               ! number of aos
         nmo, &               ! number of mos
         nz=4                 ! number of active quaternion units

        logical             :: &
         allocated=.false.,    &        ! coefficient array is allocated
         restricted=.true.,    &        ! spinors (or spin-orbitals) are paired
         complex_values=.true.          ! coefficients are complex

        character           :: &
         algebra=' '                    ! c=complex, r=real, q=quaternion, value can be set explicitly or by alloc_mo

        real(8)             :: &
         total_energy         ! total (SCF) energy as stored on the coefficient file

        real(8), pointer    :: &
         energy (:)           ! orbital energies

        complex(8), pointer :: &
         coeff (:,:,:)        ! coefficients (AO,MO,spin) (default complex representation which is always valid)

        real(8), pointer :: &
         coeff_r (:,:,:)      ! coefficients (AO,MO,spin) (alternate real representation, can be used for spinfree calculations)

        real(8), pointer :: &
         coeff_q (:,:,:)      ! coefficients (AO,MO,IZ)  (alternate quaternion representation, can be used for Kramers-restricted calculations)

        integer, pointer    :: &
         index(:), &          ! index of the spinor in the full list of spinors for this fermion symmetry
         symmetry(:)          ! possibility for an additional symmetry label (e.g. in case of spinfree calculations)

        integer             :: &
         fermion_irrep        ! fermion irrep (1=gerade, 2=ungerade) to which the spinors belongs

    end type cmo

    interface alloc_mo
      module procedure alloc_mo
    end interface

    interface dealloc_mo
      module procedure dealloc_mo
    end interface

    interface copy_mo
      module procedure copy_mo
    end interface

    interface sync_mo
      module procedure sync_mo
    end interface

    interface compress_mo
      module procedure compress_mo_complex
    end interface

    contains

     subroutine alloc_mo(M,nao,nmo)

!     Prepare to work with a block of Mos

      type(cmo), intent(inout) :: M
      integer,   intent(in )   :: nao
      integer,   intent(in )   :: nmo

      M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
      M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)

      ! Only set algebra in case it is not yet set, in this way we can allocate in a non-default format (coeff_q instead of coeff)
      if (M%algebra == ' ') M%algebra = 'c'

!     Reserve the memory needed to work with a block of Mos
      nullify  (M%energy)
      allocate (M%energy(M%nmo))
      M%energy = 0.0

      if (M%algebra == 'c') then
         nullify  (M%coeff)
         allocate (M%coeff(M%nao,M%nmo,2))
         M%coeff = 0.0
      else if (M%algebra == 'q') then
         nullify  (M%coeff_q)
         allocate (M%coeff_q(M%nao,M%nmo,M%nz))
         M%coeff_q = 0.0
      end if

      nullify  (M%symmetry)
      allocate (M%symmetry(M%nmo))
      M%symmetry = 0

      nullify  (M%index)
      allocate (M%index(M%nmo))
      M%index = 0

      M%allocated = .true.

     end subroutine alloc_mo

     subroutine dealloc_mo(M)

!     Clean up workspace

      type(cmo), intent(inout) :: M

!     Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
      M%nao = 0
      M%nmo = 0

!     Free the memory associated with a block of Mos
      deallocate (M%energy)
      deallocate (M%coeff)
      deallocate (M%symmetry)
      deallocate (M%index)
      M%allocated = .false.

     end subroutine dealloc_mo

     subroutine copy_mo(M,M_copy)

      type(cmo), intent(inout) :: M, M_copy
      integer                  :: n_coef

      M_copy%algebra        = M%algebra
      call alloc_mo (M_copy,M%nao,M%nmo)

      M_copy%restricted     = M%restricted
      M_copy%complex_values = M%complex_values
      M_copy%total_energy   = M%total_energy
      M_copy%fermion_irrep  = M%fermion_irrep
      M_copy%symmetry       = M%symmetry
      M_copy%index          = M%index
      M_copy%energy         = M%energy
      if (M%algebra == 'c') M_copy%coeff = M%coeff
      if (M%algebra == 'q') M_copy%coeff_q = M%coeff_q

     end subroutine copy_mo

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))

     subroutine sync_mo(M,my_MPI_master)

      use interface_to_mpi

!     Synchronize and allocate (if not on master) a block of MOcoefficients

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master
      integer                  :: my_MPI_rank
      integer                  :: n_coef

      call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
      call interface_mpi_bcast (M%nao,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%nmo,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%total_energy,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%algebra,1,my_MPI_master,global_communicator)
      if (my_MPI_rank /= my_MPI_master) then
         call alloc_mo (M,M%nao,M%nmo)
         M%allocated = .true.
      end if
      call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%symmetry,     M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
      n_coef = M%nao * M%nmo * 4 ! (AO,MO,spin) * 2 for complex or (AO,MO,4) for quaternion
      if (M%algebra == 'c') call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)
      if (M%algebra == 'q') call interface_mpi_bcast (M%coeff_q,n_coef,my_MPI_master,global_communicator)

     end subroutine sync_mo

#else

     subroutine sync_mo(M,my_MPI_master)

!     do nothing, makes it possible to call sync also in serial runs

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master

     end subroutine sync_mo

#endif

     subroutine compress_mo_complex(M,ao_map,mo_map)

!     Reduce or reorder MO coefficients based on AO and MO compression maps

      type(cmo), intent(inout) :: M
      integer,   intent(in), optional   :: ao_map(:,:)
      integer,   intent(in), optional   :: mo_map(:,:)

      integer, allocatable :: map_ao(:,:), map_mo(:,:)
      integer :: i, j, nao_new, nmo_new
      real(8), allocatable    :: energy(:)
      complex(8), allocatable :: coeff(:,:,:)
      integer, allocatable    :: index(:), symmetry(:)

      ! start by allocating local copies of the ao and mo maps (with a sanity check)
      if (present(ao_map)) then
         nao_new = size(ao_map,2)
         if (size(ao_map,1) /= 2 .or. nao_new > M%nao) stop 'compress_mo: inconsistent ao_map'
         allocate (map_ao(2,nao_new))
         map_ao = ao_map
      else
         nao_new = M%nao
         allocate (map_ao(2,nao_new))
         forall (i=1:nao_new) map_ao(1:2,i) = i
      end if

      if (present(mo_map)) then
         nmo_new = size(mo_map,2)
         if (size(mo_map,1) /= 2 .or. nmo_new > M%nmo) stop 'compress_mo: inconsistent mo_map'
         allocate (map_mo(2,nmo_new))
         map_mo = mo_map
      else
         nmo_new = M%nmo
         allocate (map_mo(2,nmo_new))
         forall (i=1:nmo_new) map_mo(1:2,i) = i
      end if

      ! additional sanity checks
      if (any(map_ao<1))            stop 'compress_mo: corrupt ao_map'
      if (any(map_ao(1,:)>M%nao))   stop 'compress_mo: corrupt ao_map'
      if (any(map_ao(2,:)>nao_new)) stop 'compress_mo: corrupt ao_map'
      if (any(map_mo<1))            stop 'compress_mo: corrupt mo_map'
      if (any(map_mo(1,:)>M%nmo))   stop 'compress_mo: corrupt mo_map'
      if (any(map_mo(2,:)>nmo_new)) stop 'compress_mo: corrupt mo_map'

      ! resort and store the MO information in the object
      ! we only need to change this if we reorder MOs

      if (present(mo_map)) then

         allocate (energy(nmo_new))
         allocate (symmetry(nmo_new))
         allocate (index(nmo_new))

         do i = 1, nmo_new
            energy(map_mo(2,i))      = M%energy(map_mo(1,i))
            symmetry(map_mo(2,i)) = M%symmetry(map_mo(1,i))
            index(map_mo(2,i))       = M%index(map_mo(1,i))
         enddo

         deallocate (M%energy)
         allocate (M%energy(nmo_new))
         M%energy = energy

         deallocate  (M%symmetry)
         allocate (M%symmetry(nmo_new))
         M%symmetry = symmetry

         deallocate (M%index)
         allocate (M%index(nmo_new))
         M%index = index

      end if

      ! resort the coefficients
      allocate (coeff(nao_new,nmo_new,2))
      do j = 1, nmo_new
         do i = 1, nao_new
            coeff(map_ao(2,i),map_mo(2,j),1) = M%coeff(map_ao(1,i),map_mo(1,j),1)
            coeff(map_ao(2,i),map_mo(2,j),2) = M%coeff(map_ao(1,i),map_mo(1,j),2)
         end do
      end do

      ! store the resorted coefficients in the spinor object
      deallocate  (M%coeff)
      allocate (M%coeff(nao_new,nmo_new,2))
      M%coeff = coeff

      ! we can now update the dimension information
      M%nao = nao_new
      M%nmo = nmo_new

     end subroutine compress_mo_complex

     subroutine convert_mo_to_complex(cm)

!     switch to Kramer's unrestricted picture

      implicit none
      type(cmo),  intent(inout) :: cm
      integer                   :: i, j
      integer, allocatable      :: symmetry(:),index(:)
      real(8), allocatable      :: energy(:)

      ! Start by sanity check
      if (cm%algebra == 'c') return ! nothing to be done, is already complex
      if (cm%algebra == 'r') stop 'real to complex conversion is not yet implemented'

      ! Adjust dimensions to reflect the new format
      cm%nmo = cm%nmo*2  ! the number of spinors is twice the number of Kramers pairs

      ! Adjust arrays that depend on the number of spinors to spinors format (twice the size)
      allocate(energy(cm%nmo))
      allocate(symmetry(cm%nmo))
      allocate(index(cm%nmo))
      do j = 1, cm%nmo/2
        energy(2*j-1)      = cm%energy(j)
        energy(2*j)        = cm%energy(j)
        symmetry(2*j-1) = cm%symmetry(j)
        symmetry(2*j)   = cm%symmetry(j)
        index(2*j-1)       = cm%index(j)
        index(2*j)         = cm%index(j)
      end do
      deallocate (cm%energy,cm%symmetry,cm%index)
      allocate (cm%energy(cm%nmo))
      allocate (cm%symmetry(cm%nmo))
      allocate (cm%index(cm%nmo))
      cm%energy = energy
      cm%symmetry = symmetry
      cm%index = index
      deallocate (energy,symmetry,index)

      ! Allocate the array to store coefficients in complex format and copy these
      nullify  (cm%coeff)
      allocate (cm%coeff(cm%nao,cm%nmo,2))
      cm%coeff = 0.0
      do j = 1, cm%nmo/2
         do i = 1, cm%nao
          cm%coeff(i,2*j-1,1) = dcmplx( cm%coeff_q(i,j,1), cm%coeff_q(i,j,2)) ! unbar-alpha
          cm%coeff(i,2*j-1,2) = dcmplx(-cm%coeff_q(i,j,3), cm%coeff_q(i,j,4)) ! unbar-beta
          cm%coeff(i,2*j  ,1) = dcmplx( cm%coeff_q(i,j,3), cm%coeff_q(i,j,4)) ! bar-alpha
          cm%coeff(i,2*j  ,2) = dcmplx( cm%coeff_q(i,j,1),-cm%coeff_q(i,j,2)) ! bar-beta
         end do
      end do

      ! we can remove the quaternion representation now and end by adjusting the state of this object
      deallocate  (cm%coeff_q)
      cm%algebra = 'c'

     end subroutine convert_mo_to_complex

     subroutine convert_mo_to_quaternion(cm)

!     switch to Kramer's restricted format, taking unbarred spinors (odd indices) as reference
!     NB: will not work for unrestricted calculations in which unbar/bar are not intertwined exactly
!         and even if they are one should be careful with degenerate Kramers pairs as applying a
!         Kramers operator is then not guaranteed to produce an orthogonal set of spinors

      implicit none
      type(cmo),  intent(inout) :: cm
      integer                   :: i, j
      integer, allocatable      :: symmetry(:),index(:)
      real(8), allocatable      :: energy(:)

      ! Start by sanity checks
      if (cm%algebra == 'q') return ! nothing to be done, is already quaternion
      if (cm%algebra == 'r') stop 'real to quaternion conversion is not yet implemented'
      if (.not.cm%restricted) stop "can not make quaternion representation for unrestricted spinors"
      if (mod(cm%nmo,2) /= 0) stop "can not make quaternion representation for odd number of spinors"

      ! Adjust dimensions to reflect the new format
      cm%nmo = cm%nmo/2  ! the number of Kramers pairs is half the number of spinors
      cm%nz = 4          ! nz=4 as we have nosym, this should be generalized when symmetry is implemented.

      ! Adjust arrays that depend on the number of spinors to Kramers format (half the size)
      allocate(energy(cm%nmo))
      allocate(symmetry(cm%nmo))
      allocate(index(cm%nmo))
      do j = 1, cm%nmo
        energy(j)      = cm%energy(2*j-1)
        symmetry(j) = cm%symmetry(2*j-1)
        index(j)       = cm%index(2*j-1)
      end do
      deallocate (cm%energy,cm%symmetry,cm%index)
      allocate (cm%energy(cm%nmo))
      allocate (cm%symmetry(cm%nmo))
      allocate (cm%index(cm%nmo))
      cm%energy = energy
      cm%symmetry = symmetry
      cm%index = index
      deallocate (energy,symmetry,index)

      ! Allocate the array to store coefficients in quaternion format and copy these
      nullify  (cm%coeff_q)
      allocate (cm%coeff_q(cm%nao,cm%nmo,cm%nz))
      cm%coeff_q = 0.0
      do j = 1, cm%nmo
         do i = 1, cm%nao
          cm%coeff_q(i,j,1) =  real(cm%coeff(i,2*j-1,1))
          cm%coeff_q(i,j,2) = aimag(cm%coeff(i,2*j-1,1))
          cm%coeff_q(i,j,3) = -real(cm%coeff(i,2*j-1,2))
          cm%coeff_q(i,j,4) = aimag(cm%coeff(i,2*j-1,2))
         end do
      end do

      ! we can remove the complex representation now and end by adjusting the state of this object
      deallocate  (cm%coeff)
      cm%algebra = 'q'

     end subroutine convert_mo_to_quaternion

     subroutine apply_Kramers_operator(cm,cmbar)

!     apply Kramers operator: returns time-reversed set of spinors

      implicit none
      type(cmo),  intent(inout) :: cm    ! original set of spinors
      type(cmo), intent(inout)  :: cmbar ! time-reversed set of spinors
      integer                   :: i, j

      if (cmbar%allocated) call dealloc_mo(cmbar)
      call alloc_mo (cmbar,cm%nao,cm%nmo)

      cmbar%total_energy  = cm%total_energy
      cmbar%fermion_irrep = cm%fermion_irrep

      do j = 1, cm%nmo
        cmbar%energy(j)      = cm%energy(j)
        cmbar%index(j)       = cm%index(j)
        cmbar%symmetry(j) = cm%symmetry(j)
      end do

      do j = 1, cm%nmo
         do i = 1, cm%nao
          cmbar%coeff(i,j,1) = - dconjg(cm%coeff(i,j,2)) ! c_alpha <- - c_beta*
          cmbar%coeff(i,j,2) =   dconjg(cm%coeff(i,j,1)) ! c_beta  <- + c_alpha*
         end do
      end do

     end subroutine apply_Kramers_operator


end module exacorr_mo
