module rose_mo 

!   Defines all the necessary information to manipulate molecular orbitals
!   Minimally required - Dimensions (size of AO basis, number of MOs)
!                      - Orbital energies
!                      - Coefficients of MOs in terms of AO-functions
!                      - Algebra of the coefficients (real or complex)
!                      - Pairing of orbitals (restricted or unrestricted)
!                      - Mixing of alpha and beta spin functions (spin-orbitals or M)
!   Optional           - Reference total molecular energy obtained with this MO set
!                      - Information about the symmetry species of these orbitals

!   Taken from earlier DIRAC code but modified for general usage. Lucas Visscher, 2019.

    implicit none

    private

    public :: mo_basis, alloc_mo_basis, dealloc_mo_basis, fill_mo_basis, copy_mo_basis
    public :: reduce_mo_basis, insert_mo_basis, convert_mo_basis, imo_to_mobasis, apply_Kramers_operator
    public :: mo_type, scale_spinor

    type  :: quaternion
        real(8),dimension(4) :: val
    end type quaternion

    type  :: mo_basis
        integer             :: &
         nao, &               ! number of aos
         nmo, &               ! number of mos
         type=0               ! type of MOs, 0 means that data is not yet intialized

        real(8)             :: &
         total_energy=0.D0    ! total (SCF) energy as taken from the coefficient file

        real(8), pointer    :: &
         energy (:)           ! orbital energies as taken from the coefficient file

        real(8), pointer :: &
         coeff_r (:,:,:)      ! coefficients (AO,MO,spin)

        complex(8), pointer :: &
         coeff_c (:,:,:)      ! coefficients (AO,MO,spin)

        real(8), pointer :: &
         coeff_q (:,:,:)      ! coefficients (AO,MO,IZ)

        integer, pointer    :: &
         irrep(:)             ! symmetry species (irreducible representation number) of the orbital

    end type mo_basis

    interface alloc_mo_basis
      module procedure alloc_mo_basis
    end interface

    interface dealloc_mo_basis
      module procedure dealloc_mo_basis
    end interface

    interface reduce_mo_basis
      module procedure reduce_mo_basis
    end interface

    interface insert_mo_basis
      module procedure insert_mo_basis
    end interface

    interface fill_mo_basis
      module procedure fill_mo_basis
    end interface

    interface copy_mo_basis
      module procedure copy_mo_basis
    end interface

    interface scale_spinor
      module procedure scale_spinor
    end interface

    contains

     subroutine alloc_mo_basis(M,nao,nmo,mo_type)

!     Prepare to work with a block of Mos

      type(mo_basis), intent(inout) :: M
      integer,   intent(in ) :: nao
      integer,   intent(in ) :: nmo
      integer,   intent(in ), optional :: mo_type
      integer                          :: l_type

      M%nao = nao  ! Size of AO basis.
      M%nmo = nmo  ! Size of MO basis.

!     Reserve the memory needed to work with a block of Mos
      nullify  (M%energy)
      allocate (M%energy(M%nmo))
      M%energy = 0.0

      ! make local copy of the optional type argument (can not be tested directly if not specified)
      if (present(mo_type)) then
         l_type = mo_type
      else
         l_type = 1
      end if 

      select case (l_type)
         case (1)
            ! type 1: standard RHF case
            nullify  (M%coeff_r)
            allocate (M%coeff_r(M%nao,M%nmo,1)) 
            M%coeff_r = 0.0
            M%type = 1
         case (2)
            ! type 2: standard UHF case
            nullify  (M%coeff_r)
            allocate (M%coeff_r(M%nao,M%nmo,2))
            M%coeff_r = 0.0
            M%type = 2
         case (3)
            ! type 3: spin-orbit UHF case, we work with complex orbitals
            nullify  (M%coeff_c)
            allocate (M%coeff_c(M%nao,M%nmo,2))
            M%coeff_c = 0.0
            M%type = 3
         case (4)
            ! type 4: spin-orbit KR-HF case, we work with quaternion orbitals
            nullify  (M%coeff_q)
            allocate (M%coeff_q(M%nao,M%nmo,4))
            M%coeff_q = 0.0
            M%type = 4
         case default
            stop "unsupported type of spinors requested"
      end select

      nullify  (M%irrep)
      allocate (M%irrep(M%nmo))
      M%irrep = 0

     end subroutine alloc_mo_basis

     subroutine dealloc_mo_basis(M)

!     Clean up workspace

      type(mo_basis), intent(inout) :: M

!     Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
      M%nao = 0
      M%nmo = 0

!     Free the memory associated with a block of Mos
      if (M%type > 0) then
         deallocate (M%energy)
         deallocate (M%irrep)
         if (M%type <= 2)  deallocate (M%coeff_r)
         if (M%type == 3)  deallocate (M%coeff_c)
         if (M%type == 4)  deallocate (M%coeff_q)
      end if
      M%type = 0

     end subroutine dealloc_mo_basis

     subroutine reduce_mo_basis (M,nmo_new,start,mo_list)

      ! Reduce the size of the mo_basis
      type(mo_basis), intent(inout) :: M
      integer, intent(in)           :: nmo_new
      integer, intent(in), optional :: start
      integer, intent(in), optional :: mo_list(:)
      integer                       :: offset,i

      real(8), allocatable          :: coeff(:,:,:)
      complex(8), allocatable       :: coeff_c(:,:,:)
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)
      integer, allocatable          :: mo_index(:)

      if (present(start)) then
         offset = start
      else
         offset = 1
      end if

      ! do a sanity check first, size should match
      if (M%nmo < offset+nmo_new-1) then
          print*,"we have ",M%nmo,'orbitals, but the range ',offset,':',offset+nmo_new-1,' is adressed in reduce_mo_basis'
          error stop "inconsistency in call to reduce_mo_basis"
      end if

      allocate(mo_index(nmo_new))
      if (present(mo_list)) then
         ! sanity checks before copying into the local array
         if (present(start)) error stop 'both start and mo_list specified in call to reduce_mo_basis'
         if (nmo_new .ne. size(mo_list) .or. maxval(mo_list) > M%nmo) then
             print*, M%nmo,"orbitals,are to be reduced to ",nmo_new," orbitals, but the list of indices is:"
             print*, mo_list
             error stop "inconsistency in call to reduce_mo_basis"
         end if
         mo_index = mo_list
      else
         ! contiguous copy of spinors
         forall (i=1:nmo_new) mo_index(i) = offset+i-1
      end if

      allocate(energy(nmo_new))
      forall (i=1:nmo_new) energy(i) = M%energy(mo_index(i))
      !energy = M%energy(offset:offset+nmo_new-1)
      deallocate(M%energy)
      allocate(M%energy(nmo_new))
      M%energy = energy
      deallocate(energy)

      allocate(irrep(nmo_new))
      forall (i=1:nmo_new) irrep(i) = M%irrep(mo_index(i))
      !irrep = M%irrep(offset:offset+nmo_new-1)
      deallocate(M%irrep)
      allocate(M%irrep(nmo_new))
      M%irrep = irrep
      deallocate(irrep)

      select case (M%type)
         case (1,2) ! real (un)restricted
            allocate(coeff(M%nao,nmo_new,size(M%coeff_r,3)))
            forall (i=1:nmo_new) coeff(:,i,:) = M%coeff_r(:,mo_index(i),:)
            !coeff = M%coeff_r(:,offset:offset+nmo_new-1,:)
            deallocate(M%coeff_r)
            allocate(M%coeff_r(M%nao,nmo_new,size(coeff,3)))
            M%coeff_r = coeff
            deallocate(coeff)
         case (3) ! complex unrestricted
            allocate(coeff_c(M%nao,nmo_new,2))
            !coeff_c = M%coeff_c(:,offset:offset+nmo_new-1,:)
            forall (i=1:nmo_new) coeff_c(:,i,:) = M%coeff_c(:,mo_index(i),:)
            deallocate(M%coeff_c)
            allocate(M%coeff_c(M%nao,nmo_new,2))
            M%coeff_c = coeff_c
            deallocate(coeff_c)
         case (4) ! quaternion restricted
            allocate(coeff(M%nao,nmo_new,4))
            !coeff = M%coeff_q(:,offset:offset+nmo_new-1,:)
            forall (i=1:nmo_new) coeff(:,i,:) = M%coeff_q(:,mo_index(i),:)
            deallocate(M%coeff_q)
            allocate(M%coeff_q(M%nao,nmo_new,4))
            M%coeff_q = coeff
            deallocate(coeff)
         case default
      end select

      M%nmo = nmo_new

     end subroutine reduce_mo_basis

     subroutine insert_mo_basis (M,M1,after)

      ! Insert the MOs from M1 into M at position after
      type(mo_basis), intent(inout) :: M, M1
      integer, intent(in), optional :: after
      integer                       :: offset, nmo_new

      real(8), allocatable          :: coeff(:,:,:)
      complex(8), allocatable       :: coeff_c(:,:,:)
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      ! do a sanity check first, size and type should match
      if (M%nao.ne.M1%nao .or. M%type.ne.M1%type) stop "inconsistency in  insert_mo_basis"

      if (present(after)) then
         offset = after
         if (offset < 0 .or. offset > M%nmo) stop "wrong offset in insert_mo_basis"
      else
         offset = 0
      end if

      nmo_new = M%nmo + M1%nmo

      allocate(energy(nmo_new))
      if (offset > 0) energy(1:offset) = M%energy(1:offset)
      energy(offset+1:offset+M1%nmo)   = M1%energy
      if (offset < M%nmo) energy(offset+M1%nmo+1:nmo_new)  = M%energy(offset+1:M%nmo)
      deallocate(M%energy)
      allocate(M%energy(nmo_new))
      M%energy = energy
      deallocate(energy)

      allocate(irrep(nmo_new))
      if (offset > 0) irrep(1:offset) = M%irrep(1:offset)
      irrep(offset+1:offset+M1%nmo)   = M1%irrep
      if (offset < M%nmo) irrep(offset+M1%nmo+1:nmo_new)  = M%irrep(offset+1:M%nmo)
      deallocate(M%irrep)
      allocate(M%irrep(nmo_new))
      M%irrep = irrep
      deallocate(irrep)

      select case (M%type)
         case (1,2) ! real (un)restricted
            allocate(coeff(M%nao,nmo_new,size(M%coeff_r,3)))
            if (offset > 0) coeff(:,1:offset,:) = M%coeff_r(:,1:offset,:)
            coeff(:,offset+1:offset+M1%nmo,:) = M1%coeff_r
            if (offset < M%nmo) coeff(:,offset+M1%nmo+1:nmo_new,:) = M%coeff_r(:,offset+1:M%nmo,:)
            deallocate(M%coeff_r)
            allocate(M%coeff_r(M%nao,nmo_new,size(coeff,3)))
            M%coeff_r = coeff
            deallocate(coeff)
         case (3) ! complex unrestricted
            allocate(coeff_c(M%nao,nmo_new,2))
            if (offset > 0) coeff_c(:,1:offset,:) = M%coeff_c(:,1:offset,:)
            coeff_c(:,offset+1:offset+M1%nmo,:) = M1%coeff_c
            if (offset < M%nmo) coeff(:,offset+M1%nmo+1:nmo_new,:) = M%coeff_c(:,offset+1:M%nmo,:)
            deallocate(M%coeff_c)
            allocate(M%coeff_c(M%nao,nmo_new,2))
            M%coeff_c = coeff_c
            deallocate(coeff_c)
         case (4) ! quaternion restricted
            allocate(coeff(M%nao,nmo_new,4))
            if (offset > 0) coeff(:,1:offset,:) = M%coeff_q(:,1:offset,:)
            coeff(:,offset+1:offset+M1%nmo,:) = M1%coeff_q
            if (offset < M%nmo) coeff(:,offset+M1%nmo+1:nmo_new,:) = M%coeff_q(:,offset+1:M%nmo,:)
            deallocate(M%coeff_q)
            allocate(M%coeff_q(M%nao,nmo_new,4))
            M%coeff_q = coeff
            deallocate(coeff)
         case default
      end select

      M%nmo = nmo_new

     end subroutine insert_mo_basis

     subroutine copy_mo_basis (M1,M2)

      ! Copy the contents of M1 to M2
      type(mo_basis), intent(in) :: M1
      type(mo_basis), intent(out):: M2

      call dealloc_mo_basis(M2) ! does nothing if it was not allocated
      call alloc_mo_basis(M2,M1%nao,M1%nmo,M1%type)
     
      M2%total_energy = M1%total_energy
      M2%energy       = M1%energy
      M2%irrep        = M1%irrep
      select case (M1%type)
         case (1,2) ! real (un)restricted
            M2%coeff_r = M1%coeff_r
         case (3) ! complex unrestricted
            M2%coeff_c = M1%coeff_c
         case (4) ! quaternion restricted
            M2%coeff_q = M1%coeff_q
          case default
      end select

     end subroutine copy_mo_basis

     subroutine fill_mo_basis(M,start,coeff,energy,total_energy)

!     Copy data into this type

      type(mo_basis), intent(inout) :: M

      integer, intent(in), optional    :: start
      real(8), intent(in), optional    :: coeff(:,:,:)
      real(8), intent(in), optional    :: energy(:)
      real(8), intent(in), optional    :: total_energy

      integer :: offset, nao, nmo_act, n

!     Simplest copy first
      if (present(total_energy)) then
         M%total_energy = total_energy
      end if

!     We may work with an offset and leave the first MOs untouched
      if (present(start)) then
         offset = start
      else
         offset = 1
      end if

      if (present(energy)) then
         nmo_act = size(energy)
         ! do a sanity check first, size should match this MO type
         if (M%nmo < offset+nmo_act-1) stop "inconsistency in call to fill_mo_basis"
         M%energy(offset:offset+nmo_act-1) = energy
      end if

      if (present(coeff)) then
         nao     = size(coeff,1)
         nmo_act = size(coeff,2)
         n       = size(coeff,3)
         ! do sanity checks, size should match this MO type
         if (M%nao /= nao)             stop "inconsistency in call to fill_mo_basis"
         if (M%nmo < offset+nmo_act-1) stop "inconsistency in call to fill_mo_basis"
         select case (M%type)
            case (1) ! real restricted
               if (n /= 1) stop "inconsistency in call to fill_mo_basis"
               M%coeff_r(:,offset:offset+nmo_act-1,:) = coeff
            case (2) ! real unrestricted
               if (n /= 2) stop "inconsistency in call to fill_mo_basis"
               M%coeff_r(:,offset:offset+nmo_act-1,:) = coeff
            case (3) ! complex unrestricted
               if (n /= 4) stop "inconsistency in call to fill_mo_basis"
               M%coeff_c(:,offset:offset+nmo_act-1,:) = cmplx(coeff(:,:,1:2),coeff(:,:,3:4))
            case (4) ! quaternion restricted
               if (n /= 4) stop "inconsistency in call to fill_mo_basis"
               M%coeff_q(:,offset:offset+nmo_act-1,:) = coeff
         case default
         end select
      end if

     end subroutine fill_mo_basis

     subroutine convert_mo_basis(M,newtype)

!     Convert MO basis to a different type

      type(mo_basis), intent(inout) :: M
      integer, intent(in)           :: newtype
      integer                       :: oldtype

      oldtype = M%type

      if (oldtype == newtype) return ! nothing to be done

      ! call the right converter
      ! second argument is the variable that indicates whether input or output real coefficients are restricted
      ! note that the converter between the restricted and unrestricted real types is still to be written
      if     (oldtype == 1 .and. newtype == 3) then
         call convert_real_to_complex(M,.true.)
      elseif (oldtype == 1 .and. newtype == 4) then
         call convert_real_to_quaternion(M)
      elseif (oldtype == 2 .and. newtype == 3) then
         call convert_real_to_complex(M,.false.)
      elseif (oldtype == 2 .and. newtype == 4) then
         call convert_real_to_quaternion(M)
      elseif (oldtype == 3 .and. newtype == 1) then
         call convert_complex_to_real(M,.true.)
      elseif (oldtype == 3 .and. newtype == 2) then
         call convert_complex_to_real(M,.false.)
      elseif (oldtype == 3 .and. newtype == 4) then
         call convert_complex_to_quaternion(M)
      elseif (oldtype == 4 .and. newtype == 1) then
         call convert_quaternion_to_real(M,.true.)
      elseif (oldtype == 4 .and. newtype == 2) then
         call convert_quaternion_to_real(M,.false.)
      elseif (oldtype == 4 .and. newtype == 3) then
         call convert_quaternion_to_complex(M)
      else
         write (*,'(A,i2,A,i2,A)') "Conversion from type",oldtype," to type",newtype," is not possible."
         stop "error in convert_mo_basis"
      end if

     end subroutine convert_mo_basis

     subroutine IMO_to_mobasis(IMO,spinor_original,spinor_new,energies,last_mo)

        ! Takes intrinsic molecular orbital coefficients to construct a new mo_basis object

        use linear_algebra

        implicit none
        real(8), intent(in)               :: IMO(:,:,:)           ! (part of) unitary transformation in MO basis
        type(mo_basis),intent(in)         :: spinor_original      ! original MOs in AO basis B1
        type(mo_basis),intent(out)        :: spinor_new           ! transformed MOs in AO basis B1
        real(8), optional, intent(in)     :: energies(:)
        integer, optional, intent(in)     :: last_mo              ! possibility to use fewer MOs (truncation of virtual space)

        integer                           :: nmo1,nao1,nmo2,nmo2r
        real(8), allocatable              :: IMO_aoB1(:,:,:)      ! IMOs in terms of AO basis B1
        complex(8), allocatable           :: IMO_aoB1_c(:,:,:)    ! temporary complex format (needed for type 3)

        nao1 = spinor_original%nao
        nmo1 = size(IMO,1)
        nmo2 = size(IMO,2)

        if (present(last_mo)) then
           nmo2r = last_mo
           if (nmo2r > nmo2) stop "IMO_to_mobasis: last mo cannot be larger than nmo2"
        else
           nmo2r = nmo2
        end if

        !-----Express spinors in terms of the original AO basis B1----- (nao1 x nmo2)
        select case (spinor_original%type)
           case (1)
              allocate(IMO_aoB1(nao1,nmo2r,1))
              call matrix_mul(spinor_original%coeff_r,IMO(:,1:nmo2r,:),IMO_aoB1)
           case (2)
              allocate(IMO_aoB1(nao1,nmo2r,2)) ! third dimension is the spin
              IMO_aoB1(:,:,1) = matmul(spinor_original%coeff_r(:,:,1),IMO(:,1:nmo2r,1))
              IMO_aoB1(:,:,2) = matmul(spinor_original%coeff_r(:,:,2),IMO(:,1:nmo2r,1))
           case (3)
              !multiplication needs to be done with complex arithmetic
              allocate(IMO_aoB1_c(nao1,nmo2r,2)) ! third dimension is just the spin in the complex format
              IMO_aoB1_c(:,:,1) = matmul(spinor_original%coeff_c(:,:,1),dcmplx(IMO(:,1:nmo2r,1),IMO(:,1:nmo2r,2)))
              IMO_aoB1_c(:,:,2) = matmul(spinor_original%coeff_c(:,:,2),dcmplx(IMO(:,1:nmo2r,1),IMO(:,1:nmo2r,2)))
              !copy into the real array
              allocate(IMO_aoB1(nao1,nmo2r,4))   ! third dimension is spin and real / imaginary
              IMO_aoB1(:,:,1:2) = real(IMO_aoB1_c)
              IMO_aoB1(:,:,3:4) = imag(IMO_aoB1_c)
              deallocate(IMO_aoB1_c)
           case (4)
              allocate(IMO_aoB1(nao1,nmo2r,4)) ! third dimension is the quaternion unit
              call matrix_mul(spinor_original%coeff_q,IMO(:,1:nmo2r,:),IMO_aoB1)
        end select

        ! Write the coefficient matrix into the spinor type.
        call alloc_mo_basis(spinor_new,nao1,nmo2r,spinor_original%type)
        spinor_new%total_energy = spinor_original%total_energy
        call fill_mo_basis(spinor_new,coeff=IMO_aoB1)
        if (present(energies)) call fill_mo_basis(spinor_new,energy=energies(1:nmo2r))
        deallocate(IMO_aoB1)

     end subroutine IMO_to_mobasis

     subroutine convert_quaternion_to_complex(M)

!     switch to Kramer's unrestricted picture
!     we will work with twice the number of spinors instead of Kramers pairs

      type(mo_basis), intent(inout) :: M
      integer                       :: i, j, spin
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      if (M%type .ne. 4) stop 'error in convert_quaternion_to_complex'

      ! we have degeneracy in the energy, irrep number may need modification when we really start using symmetry
      ! start by copying the energy and symmetry info to the longer list that is to be used
      allocate (energy(2*M%nmo))
      allocate (irrep (2*M%nmo))
      do j = 1, M%nmo
         energy(2*j-1) = M%energy(j)
         energy(2*j  ) = M%energy(j)
         irrep(2*j-1)  = M%irrep(j)
         irrep(2*j  )  = M%irrep(j)
      end do
      deallocate (M%energy)
      deallocate (M%irrep)
      allocate(M%energy(2*M%nmo))
      allocate(M%irrep(2*M%nmo))
      M%energy = energy
      M%irrep  = irrep
      deallocate (energy)
      deallocate (irrep)

      allocate (M%coeff_c(M%nao,2*M%nmo,2)) ! same procedure now for the coefficients, where we use Kramers symmetry
      do spin = 1, 2
        if (spin == 1) then
            do j = 1, M%nmo
                do i = 1, M%nao
                M%coeff_c(i,2*j-1,spin) = dcmplx(M%coeff_q(i,j,1),M%coeff_q(i,j,2)) ! j_alpha 
                M%coeff_c(i,2*j  ,spin) = dcmplx(M%coeff_q(i,j,3),M%coeff_q(i,j,4)) ! jbar_alpha = -j_beta*
                end do
            end do
        else if (spin == 2) then
            do j = 1, M%nmo
                do i = 1, M%nao
                M%coeff_c(i,2*j-1,spin) = dcmplx(-M%coeff_q(i,j,3), M%coeff_q(i,j,4)) ! j_beta
                M%coeff_c(i,2*j  ,spin) = dcmplx( M%coeff_q(i,j,1),-M%coeff_q(i,j,2)) ! jbar_beta = j_alpha*
                end do
            end do
        end if
      end do

      ! set to new type and clean up
      M%type = 3
      M%nmo  = 2 * M%nmo
      deallocate (M%coeff_q)
      nullify    (M%coeff_q)

     end subroutine convert_quaternion_to_complex

     subroutine convert_complex_to_quaternion(M)

      implicit none

!     reverse of convert_quaternion_to_complex
!     note that we assume that the first half of the spinors will form a suitable set of Kramers pairs
!     we only test that we have an even number of spinors, this is only a minimal requirement !

      type(mo_basis), intent(inout) :: M
      integer                       :: i, j
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      if (M%type .ne. 3 .or. mod(M%nmo,2) .ne. 0) stop 'error in convert_complex_to_quaternion'

      ! we have degeneracy in the energy, irrep number may need modification when we really start using symmetry
      ! start by copying the energy and symmetry info to the shorter list that is to be used
      allocate (energy(M%nmo/2))
      allocate (irrep (M%nmo/2))
      do j = 1, M%nmo/2
         energy(j) = M%energy(2*j-1)
         irrep(j)  = M%irrep(2*j-1)
      end do
      deallocate (M%energy)
      deallocate (M%irrep)
      allocate(M%energy(M%nmo/2))
      allocate(M%irrep(M%nmo/2))
      M%energy = energy
      M%irrep  = irrep
      deallocate (energy)
      deallocate (irrep)

      allocate (M%coeff_q(M%nao,M%nmo/2,4)) ! same procedure now for the coefficients, where we use Kramers symmetry
      do j = 1, M%nmo/2
         do i = 1, M%nao
            M%coeff_q(i,j,1) =  real(M%coeff_c(i,2*j-1,1)) ! j_alpha
            M%coeff_q(i,j,2) = aimag(M%coeff_c(i,2*j-1,1)) ! j_alpha
            M%coeff_q(i,j,3) = -real(M%coeff_c(i,2*j-1,2)) ! j_beta
            M%coeff_q(i,j,4) = aimag(M%coeff_c(i,2*j-1,2)) ! j_beta
         end do
      end do

      ! set to new type and clean up
      M%type = 4
      M%nmo  = M%nmo/2
      deallocate (M%coeff_c)
      nullify    (M%coeff_c)

     end subroutine convert_complex_to_quaternion

     subroutine convert_complex_to_real(M,restricted)

     ! Dirac read quaterions, that can then be converted into complex by the convert_quaternion_to_complex subroutine.
     ! For non relativistic Hamiltonian, coefficients are real and there is no need for complex algebra. Therefore,
     ! we can use this subroutine to get rid of the (zero) imaginary part.

      implicit none

      logical, intent(in)           :: restricted
      type(mo_basis), intent(inout) :: M
      integer                       :: i, j
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      if (M%type .ne. 3 .or. mod(M%nmo,2) .ne. 0) stop 'error in convert_complex_to_real'

      if (restricted) then
         allocate (M%coeff_r(M%nao,M%nmo/2,1))
         M%coeff_r = 0.D0
         M%type = 1
         allocate (energy(M%nmo/2))
         allocate (irrep (M%nmo/2))
         do j = 1, M%nmo/2
            energy(j) = M%energy(2*j-1)
            irrep(j)  = M%irrep(2*j-1)
         end do
         deallocate (M%energy)
         deallocate (M%irrep)
         allocate(M%energy(M%nmo/2))
         allocate(M%irrep(M%nmo/2))
         M%energy = energy
         M%irrep  = irrep
         deallocate (energy)
         deallocate (irrep)
         do i=1,M%nmo/2
            do j=1,M%nao
               M%coeff_r(j,i,1) = real(M%coeff_c(j,2*i-1,1)) ! j_alpha
            end do
         end do
         M%nmo = M%nmo/2
      else
         allocate (M%coeff_r(M%nao,M%nmo,2))
         M%coeff_r = 0.D0
         M%type = 2
         do i=1,M%nmo
            do j=1,M%nao
               M%coeff_r(j,i,1) = real(M%coeff_c(j,i,1)) ! j_alpha
               M%coeff_r(j,i,2) = real(M%coeff_c(j,i,2)) ! j_beta
            end do
         end do
      end if


      deallocate (M%coeff_c)
      nullify    (M%coeff_c)

     end subroutine convert_complex_to_real

     subroutine convert_quaternion_to_real(M,restricted)

     ! Dirac reads quaterions, that can then be converted into complex by the convert_quaternion_to_complex subroutine.
     ! For non relativistic Hamiltonian, coefficients are real and there is no need for complex algebra.
      
      implicit none
     
      logical                       :: restricted ! output type can be restricted or unrestricted
      type(mo_basis), intent(inout) :: M
      integer                       :: i, j
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      if (M%type .ne. 4) stop 'error in convert_quaternion_to_real'

      if (restricted) then
         allocate (M%coeff_r(M%nao,M%nmo,1))
         M%coeff_r = 0.D0
         M%type = 1
         M%coeff_r(:,:,1) = M%coeff_q(:,:,1) ! j_alpha
      else
         allocate (M%coeff_r(M%nao,2*M%nmo,2))
         M%coeff_r = 0.D0
         allocate (energy(2*M%nmo))
         allocate (irrep (2*M%nmo))
         do j = 1, M%nmo
            energy(2*j-1) = M%energy(j)
            energy(2*j  ) = M%energy(j)
            irrep(2*j-1)  = M%irrep(j)
            irrep(2*j  )  = M%irrep(j)
         end do
         deallocate (M%energy)
         deallocate (M%irrep)
         allocate(M%energy(2*M%nmo))
         allocate(M%irrep(2*M%nmo))
         M%energy = energy
         M%irrep  = irrep
         deallocate (energy)
         deallocate (irrep)
         M%type = 2
         do i=1,M%nao
            do j=1,M%nmo
               M%coeff_r(i,2*j-1,1) =  M%coeff_q(i,j,1) ! j_alpha
!               M%coeff_r(i,2*j,2)   = -M%coeff_q(i,j,3) ! j_beta
               ! coeff_q(:,:,3) is 0 for nonrel Hamiltonian.
               M%coeff_r(i,2*j,2)   = -M%coeff_q(i,j,1) ! j_beta
            end do
         end do
         M%nmo = 2*M%nmo
      end if

           
      deallocate (M%coeff_q)
      nullify    (M%coeff_q)
        
     end subroutine convert_quaternion_to_real

     subroutine convert_real_to_complex(M,restricted)

     ! Dirac read quaterions, that can then be converted into complex by the convert_quaternion_to_complex subroutine.
     ! For non relativistic Hamiltonian, coefficients are real and there is no need for complex algebra. Therefore,
     ! we can use this subroutine to get rid of the (zero) imaginary part.

      implicit none

      logical                       :: restricted
      type(mo_basis), intent(inout) :: M
      integer                       :: i, j
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      if (M%type .ne. 1 .and. M%type .ne. 2) stop 'error in convert_real_to_complex'

      if (restricted) then
         allocate (M%coeff_c(M%nao,2*M%nmo,2))
         M%coeff_c = (0.D0,0.D0)
         allocate (energy(M%nmo*2))
         allocate (irrep (M%nmo*2))
         do j = 1, M%nmo
            energy(2*j-1) = M%energy(j)
            energy(2*j)   = M%energy(j)
            irrep(2*j-1)  = M%irrep(j)
            irrep(2*j)    = M%irrep(j)
         end do
         deallocate (M%energy)
         deallocate (M%irrep)
         allocate(M%energy(M%nmo*2))
         allocate(M%irrep(M%nmo*2))
         M%energy = energy
         M%irrep  = irrep
         deallocate (energy)
         deallocate (irrep)
         do i=1,M%nmo
            do j=1,M%nao
               M%coeff_c(j,2*i-1,1) = dcmplx(M%coeff_r(j,i,1),0.D0) ! j_alpha
               M%coeff_c(j,2*i,2) = dcmplx(M%coeff_r(j,i,1),0.D0) ! j_beta same as j_alpha
            end do
         end do
         M%nmo = 2*M%nmo
      else
         allocate (M%coeff_c(M%nao,M%nmo,2))
         M%coeff_c = (0.D0,0.D0)
         M%coeff_c(:,:,1) = dcmplx(M%coeff_r(:,:,1),0.D0) ! j_alpha
         M%coeff_c(:,:,2) = dcmplx(M%coeff_r(:,:,2),0.D0) ! j_alpha
      end if

      ! set to new type and clean up
      M%type = 3
      deallocate (M%coeff_r)
      nullify    (M%coeff_r)

     end subroutine convert_real_to_complex

     subroutine convert_real_to_quaternion(M)

     ! Dirac read quaterions, that can then be converted into complex by the convert_quaternion_to_complex subroutine.
     ! For non relativistic Hamiltonian, coefficients are real and there is no need for complex algebra. Therefore,
     ! we can use this subroutine to get rid of the (zero) imaginary part.

      implicit none

      logical                       :: restricted
      type(mo_basis), intent(inout) :: M
      integer                       :: i, j
      real(8), allocatable          :: energy(:)
      integer, allocatable          :: irrep(:)

      if (M%type .ne. 1 .and. M%type .ne. 2) stop 'error in convert_real_to_complex'

      if (M%type==1) then
         ! restricted case
         allocate (M%coeff_q(M%nao,M%nmo,4)) ! same procedure now for the coefficients, where we use Kramers symmetry
         M%coeff_q = 0.D0
         M%coeff_q(:,:,1) = M%coeff_r(:,:,1) ! j_alpha
         M%coeff_q(:,:,2) = 0.D0 ! j_alpha
         M%coeff_q(:,:,3) = 0.D0
         M%coeff_q(:,:,4) = 0.D0
      else
         ! unrestricted case
         allocate (M%coeff_q(M%nao,M%nmo/2,4)) ! same procedure now for the coefficients, where we use Kramers symmetry
         M%coeff_q = 0.D0
         allocate (energy(M%nmo/2))
         allocate (irrep (M%nmo/2))
         do j = 1, M%nmo/2
            energy(j) = M%energy(2*j-1)
            irrep(j)  = M%irrep(2*j-1)
         end do
         deallocate (M%energy)
         deallocate (M%irrep)
         allocate(M%energy(M%nmo/2))
         allocate(M%irrep(M%nmo/2))
         M%energy = energy
         M%irrep  = irrep
         deallocate (energy)
         deallocate (irrep)
         do j = 1, M%nmo/2
            do i = 1, M%nao
               M%coeff_q(i,j,1) = M%coeff_r(i,2*j-1,1) ! j_alpha
               M%coeff_q(i,j,2) = 0.D0 ! j_alpha
!               M%coeff_q(i,j,3) = -M%coeff_r(i,2*j,2) ! j_beta. Actually, once we have real algebra, quaternion should be 0.D0 here... ?
               M%coeff_q(i,j,3) = 0.D0
               M%coeff_q(i,j,4) = 0.D0 ! j_beta
            end do
         end do
         M%nmo = M%nmo/2
      end if

      ! set to new type and clean up
      M%type = 4
      deallocate (M%coeff_r)
      nullify    (M%coeff_r)

     end subroutine convert_real_to_quaternion

     subroutine apply_Kramers_operator(cm,cmbar)

!     apply Kramers operator: returns time-reversed set of spinors

      implicit none
      type(mo_basis),  intent(inout) :: cm    ! original set of spinors
      type(mo_basis), intent(inout)  :: cmbar ! time-reversed set of spinors
      integer                   :: i, j

      if (cm%type .ne. 3) stop 'intended for complex spinors only'

      if (cmbar%type >0) call dealloc_mo_basis(cmbar)
      call alloc_mo_basis (cmbar,cm%nao,cm%nmo,3)

      cmbar%total_energy  = cm%total_energy

      do j = 1, cm%nmo
        cmbar%energy(j)      = cm%energy(j)
      end do

      do j = 1, cm%nmo
         do i = 1, cm%nao
          cmbar%coeff_c(i,j,1) = - dconjg(cm%coeff_c(i,j,2)) ! c_alpha <- - c_beta*
          cmbar%coeff_c(i,j,2) =   dconjg(cm%coeff_c(i,j,1)) ! c_beta  <- + c_alpha*
         end do
      end do

     end subroutine apply_Kramers_operator

     subroutine scale_spinor(cm,factor,index_spinor)

!     right multiply spinors by a constant factor (to renormalize or scale)

      implicit none
      type(mo_basis),  intent(inout) :: cm           ! set of spinors
      real(8), intent(in)            :: factor(:)    ! factor
      integer, optional, intent(in)  :: index_spinor ! index of spinor to be changed (default is to change all spinors)
      real(8), allocatable           :: coeff(:)
      complex(8)                     :: coeff_c(2)
      integer                        :: i, j, ind, nz
      type(quaternion)               :: a, b, c

      if (present(index_spinor)) then
         ind = index_spinor
      else
         ind = 0
      endif

      ! Perform sanity check and allocate coeff with proper algebra
      select case (cm%type)
      case (1)
         if (size(factor) .ne. 1) stop 'this factor is not real'
         allocate(coeff(1))
      case (2)
         if (size(factor) .ne. 1) stop 'this factor is not real'
         allocate(coeff(2))
      case (3)
         if (size(factor) .ne. 2) stop 'this factor is not complex'
      case (4)
         if (size(factor) .ne. 4) stop 'this factor is not quaternion'
      case default
         stop 'input spinors are corrupted'
      end select

      do j = 1, cm%nmo
         if (ind.ne.0 .and. ind.ne.j) cycle ! check whether we only want to scale one spinor
         select case (cm%type)
         case(1)
            do i = 1, cm%nao
             coeff(1) = cm%coeff_r(i,j,1) ! real coefficient for spatial orbital
             cm%coeff_r(i,j,1) = coeff(1) * factor(1)
            end do
         case(2)
            do i = 1, cm%nao
             coeff(1:2) = cm%coeff_r(i,j,1:2) ! real coefficient for two spinorbitals
             cm%coeff_r(i,j,1) = coeff(1) * factor(1)
             cm%coeff_r(i,j,2) = coeff(2) * factor(1)
            end do
         case(3)
            do i = 1, cm%nao
             coeff_c(1:2) = cm%coeff_c(i,j,1:2) ! complex coefficient for two spinors
             cm%coeff_c(i,j,1) = coeff(1) * dcmplx(factor(1),factor(2))
             cm%coeff_c(i,j,2) = coeff(2) * dcmplx(factor(1),factor(2))
            end do
         case(4)
            do i = 1, cm%nao
             a%val = cm%coeff_q(i,j,1:4) ! quaternion coefficient for spinor
             b%val = factor
             c = quaternion_multiplication(a,b)
             cm%coeff_q(i,j,1:4) = c%val(1:4)
            end do
         end select
      end do

     end subroutine scale_spinor

     function mo_type (restricted, spatial_orbitals) result(spinor_type)

      logical, intent(in)    :: restricted, spatial_orbitals
      integer                :: spinor_type

      if (spatial_orbitals) then
         if (restricted) then
            spinor_type = 1 ! one set of real orbitals
         else
            spinor_type = 2 ! two sets of complex orbitals
         end if
      else
         if (restricted) then
            spinor_type = 4 ! one set of quaternion orbitals
         else
            spinor_type = 3 ! two sets of complex orbitals
         end if
      end if

     end function mo_type

     function quaternion_multiplication(a,b) result(c)

         type(quaternion) :: a,b,c

         c%val(1) = a%val(1)*b%val(1) - a%val(2)*b%val(2) - a%val(3)*b%val(3) - a%val(4)*b%val(4)
         c%val(2) = a%val(1)*b%val(2) + a%val(2)*b%val(1) + a%val(3)*b%val(4) - a%val(4)*b%val(3)
         c%val(3) = a%val(1)*b%val(3) - a%val(2)*b%val(4) + a%val(3)*b%val(1) + a%val(4)*b%val(2)
         c%val(4) = a%val(1)*b%val(4) + a%val(2)*b%val(3) - a%val(3)*b%val(2) + a%val(4)*b%val(1)

     end function quaternion_multiplication

end module rose_mo
