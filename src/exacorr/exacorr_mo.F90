module exacorr_mo 

!   Contains all the necessary information for a range of molecular orbitals
!   In particular:  - Orbital energy
!                   - Boson characterization (for spinfree runs)
!                   - Coefficients of AO-functions
!                   - Position in the full list of MO vectors

!   Extended from dirac_mo by Stan Papadopoulos, January 2019


    implicit none

    private

    public :: alloc_mo, dealloc_mo, sync_mo
    public :: convert_quaternion_to_complex,convert_complex_to_quaternion

    type  :: mo_root
        integer             :: &
         fermion_irrep, &     ! fermion irrep (1=gerade, 2=ungerade) to which the spinors belongs
         nao, &               ! number of aos
         nmo                  ! number of mos

        real(8)             :: &
         total_energy         ! total (SCF) energy as stored on the coefficient file

        real(8), pointer    :: &
         energy (:)           ! orbital energies

        integer, pointer    :: &
         index(:), &          ! index of the spinor in the full list of spinors for this fermion symmetry
         boson_irrep(:)       ! for spinfree calculations: boson irrep to which the spinor belongs

    end type mo_root

    type, public, extends(mo_root) :: cmo

        complex(8), pointer :: &
         coeff (:,:,:)          ! coefficients (AO,MO,spin)

    end type cmo

    type, public, extends(mo_root) :: qumo

        integer :: &
         nz                   ! number of active quaternion units

        real(8), pointer :: &
         coeff (:,:,:)        ! coefficients (AO,MO,IZ)

    end type qumo

    interface alloc_mo
      module procedure alloc_mo_complex
      module procedure alloc_mo_quaternion
    end interface

    interface dealloc_mo
      module procedure dealloc_mo_complex
      module procedure dealloc_mo_quaternion
    end interface

    interface sync_mo
      module procedure sync_mo_complex
      module procedure sync_mo_quaternion
    end interface

    contains

     subroutine alloc_mo_complex(M,nao,nmo)

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
      allocate (M%coeff(M%nao,M%nmo,2))
      M%coeff = 0.0

      nullify  (M%boson_irrep)
      allocate (M%boson_irrep(M%nmo))
      M%boson_irrep = 0

      nullify  (M%index)
      allocate (M%index(M%nmo))
      M%index = 0

     end subroutine alloc_mo_complex

     subroutine alloc_mo_quaternion(M,nao,nmo,nz)

!      Prepare to work with a block of Mos

       type(qumo), intent(out)           :: M
       integer,    intent(in )           :: nao
       integer,    intent(in )           :: nmo
       integer,    intent(in ), optional :: nz

       M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
       M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)
       if (present(nz)) then
          M%nz = nz
       else
          M%nz = 1  ! We default to one, easy for spinfree applications
       endif

!      Reserve the memory needed to work with a block of Mos
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

     end subroutine alloc_mo_quaternion

     subroutine dealloc_mo_complex(M)

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

     end subroutine dealloc_mo_complex

     subroutine dealloc_mo_quaternion(M)

!      Clean up workspace

       type (qumo), intent(inout) :: M

!      Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
       M%nao = 0
       M%nmo = 0
       M%nz = 0

!      Free the memory associated with a block of Mos
       deallocate (M%energy)
       deallocate (M%coeff)
       deallocate (M%boson_irrep)
       deallocate (M%index)

     end subroutine dealloc_mo_quaternion

#if defined (VAR_MPI)

     subroutine sync_mo_complex(M,my_MPI_master)

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
         call alloc_mo (M,M%nao,M%nmo)
      end if
      call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%boson_irrep,  M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
      n_coef = M%nao * M%nmo * 4 ! (AO,MO,spin) * 2 for complex
      call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)

     end subroutine sync_mo_complex

     subroutine sync_mo_quaternion(M,my_MPI_master)

      use interface_to_mpi

!     Synchronize and allocate (if not on master) a block of MOcoefficients
 
      type(qumo), intent(inout) :: M
      integer,    intent(in   ) :: my_MPI_master
      integer                   :: my_MPI_rank
      integer                   :: n_coef

      call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
      call interface_mpi_bcast (M%nao,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%nmo,1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%nz ,1,my_MPI_master,global_communicator)
      if (my_MPI_rank /= my_MPI_master) then
         call alloc_mo(M,M%nao,M%nmo,M%nz)
      end if
      call interface_mpi_bcast (M%fermion_irrep,    1,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%boson_irrep,  M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%index,        M%nmo,my_MPI_master,global_communicator)
      call interface_mpi_bcast (M%energy,       M%nmo,my_MPI_master,global_communicator)
      n_coef = M%nao * M%nmo * M%nz
      call interface_mpi_bcast (M%coeff,n_coef,my_MPI_master,global_communicator)

     end subroutine sync_mo_quaternion

#else

     subroutine sync_mo_complex(M,my_MPI_master)

!     do nothing, makes it possible to call sync also in serial runs

      type(cmo), intent(inout) :: M
      integer,   intent(in   ) :: my_MPI_master

     end subroutine sync_mo_complex

     subroutine sync_mo_quaternion(M,my_MPI_master)

!     do nothing, makes it possible to call sync also in serial runs

      type(qumo), intent(inout) :: M
      integer,    intent(in   ) :: my_MPI_master

     end subroutine sync_mo_quaternion

#endif

     subroutine convert_quaternion_to_complex(coeff,scoeff,eig,seig,nao,nmo,npo)

!     switch to Kramer's unrestricted picture

!     input variables
      real(8),    intent(inout) :: coeff(:,:,:), eig(:)
      real(8),    intent(  out) :: seig(:)
      complex(8), intent(  out) :: scoeff(:,:,:)
      integer,    intent(in   ) :: nao, nmo, npo
!     summation and dummy variables
      integer                 :: i, j, spin

!     create new arrays grouped by Kramer's partners + drop "positronic" solutions

      do i = 1, npo+nmo
        seig(2*i-1) = eig(npo+i)
        seig(2*i  ) = eig(npo+i)
      end do

      do spin = 1, 2
        if (spin == 1) then
            do j = 1, nmo
                do i = 1, nao
                scoeff(i,2*j-1,spin) = dcmplx(coeff(i,j+npo,1),coeff(i,j+npo,2))
                scoeff(i,2*j  ,spin) = dcmplx(coeff(i,j+npo,3),coeff(i,j+npo,4))
                end do
            end do
        else if (spin == 2) then
            do j = 1, nmo
                do i = 1, nao
                scoeff(i,2*j-1,spin) = dcmplx(-coeff(i,j+npo,3), coeff(i,j+npo,4))
                scoeff(i,2*j  ,spin) = dcmplx( coeff(i,j+npo,1),-coeff(i,j+npo,2))
                end do
            end do
        end if
      end do

     end subroutine convert_quaternion_to_complex

     subroutine convert_complex_to_quaternion(cm,qm)

!     switch to Kramer's restricted format, taking unbarred spinors (odd indices) as reference
!     NB: will not work for unrestricted calculations in which unbar/bar are not intertwined exactly

      implicit none
      type(cmo),  intent(inout) :: cm
      type(qumo), intent(inout) :: qm
      integer                   :: i, j

      call alloc_mo (qm,cm%nao,cm%nmo/2,4)
    
      qm%total_energy  = cm%total_energy
      qm%fermion_irrep = cm%fermion_irrep

      do j = 1, qm%nmo
        qm%energy(j) = cm%energy(2*j-1)
      end do

      do j = 1, qm%nmo
         do i = 1, qm%nao
          qm%coeff(i,j,1) =  real(cm%coeff(i,2*j-1,1))
          qm%coeff(i,j,2) = aimag(cm%coeff(i,2*j-1,1))
          qm%coeff(i,j,3) = -real(cm%coeff(i,2*j-1,2))
          qm%coeff(i,j,4) = aimag(cm%coeff(i,2*j-1,2))
         end do
      end do

     end subroutine convert_complex_to_quaternion


end module exacorr_mo
