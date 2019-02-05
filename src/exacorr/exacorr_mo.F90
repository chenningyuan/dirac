module exacorr_mo 

!   Contains all the necessary information for a range of molecular orbitals
!   In particular:  - Orbital energy
!                   - Boson characterization (for spinfree runs)
!                   - Coefficients of AO-functions
!                   - Position in the full list of MO vectors

!   Extended from dirac_mo by Stan Papadopoulos, January 2019

!!!!!!!!!!!!!!!!!LOSE DEPENDENCY ON DIRAC_MO -> PUT EVERYTHING IN HERE

    implicit none

    private

    public :: alloc_mo, dealloc_mo, sync_mo
    public :: convert_quaternion_to_complex

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
         coeff (:,:)          ! coefficients (AO,MO)

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
      allocate (M%coeff(2*M%nao,M%nmo))
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
      !n_coef = M%nao * M%nmo
      n_coef = 2* M%nao * M%nmo
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

!-----input variables
      real(8),    intent(inout) :: coeff(:,:,:), eig(:)
      real(8),    intent(  out) :: seig(:)
      complex(8), intent(  out) :: scoeff(:,:,:)
      integer,    intent(in   ) :: nao, nmo, npo
!-----summation and dummy variables
      integer                 :: i, j, spin
      complex(8), allocatable :: t_cmplx(:,:)

!-----create new arrays grouped by spin + drop "positronic" solutions

      seig(    1:  nmo) = eig(npo+1:npo+nmo)
      seig(nmo+1:2*nmo) = eig(npo+1:npo+nmo)

!     initial structure of scoeff
!            nmo  nmo
!     nao (  un-   -  )

      do spin = 1, 2
        if (spin == 1) then
            do j = 1, npo+nmo
                do i = 1, nao
                scoeff(i,j    ,spin) = dcmplx(coeff(i,j+npo,1),coeff(i,j+npo,2))
                scoeff(i,j+nmo,spin) = dcmplx(coeff(i,j+npo,3),coeff(i,j+npo,4))
                end do
            end do
        else if (spin == 2) then
            do j = 1, npo+nmo
                do i = 1, nao
                scoeff(i,j    ,spin) = dcmplx(-coeff(i,j+npo,3), coeff(i,j+npo,4))
                scoeff(i,j+nmo,spin) = dcmplx( coeff(i,j+npo,1),-coeff(i,j+npo,2))
                end do
            end do
        end if
      end do

!-----sort output arrays in terms of ascending energy

      allocate(t_cmplx(nao,2))

      call sort_by_energy(scoeff,seig,t_cmplx,1,2*nmo)

      deallocate(t_cmplx)

     end subroutine convert_quaternion_to_complex

     recursive subroutine sort_by_energy(coeff,eig,t_cmplx,first,last)

!     Quicksort algorithm, taken from https://gist.github.com/t-nissie/479f0f16966925fa29ea

!-----input variables
      complex(8), intent(inout) :: coeff(:,:,:)
      complex(8), intent(inout) :: t_cmplx(:,:)
      real(8),    intent(inout) :: eig(:)
      integer,    intent(in   ) :: first, last
!-----tmp storage and summation
      integer    :: i, j
      real(8)    :: x, t

      x = eig( (first+last) / 2 )
      i = first
      j = last

      do
        do while (eig(i) < x)
            i = i + 1
        end do
        do while (x < eig(j))
            j = j - 1
        end do
        if (i >= j) exit

!       sort eigenvalues
        t      = eig(i)
        eig(i) = eig(j)
        eig(j) = t
!       sort coefficients
        t_cmplx      = coeff(:,i,:)
        coeff(:,i,:) = coeff(:,j,:)
        coeff(:,j,:) = t_cmplx

        i = i + 1
        j = j - 1
      end do
      if (first < i - 1) call sort_by_energy(coeff,eig,t_cmplx,first,i-1 )
      if (j + 1 < last ) call sort_by_energy(coeff,eig,t_cmplx,j+1  ,last)

     end subroutine sort_by_energy

end module exacorr_mo
