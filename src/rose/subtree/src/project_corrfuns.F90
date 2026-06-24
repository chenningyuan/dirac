module rose_correlating_virtuals

        private

        public :: generate_corrfuns

contains

subroutine generate_corrfuns(molecule,ao_basis_1,spinor_complement,Cadd,nadd_frag)

        use rose_ao
        use rose_mo
        use read_xyzfile
        use rose_input_data
        use rose_global

        implicit none

        type(molecule_info),intent(in)        :: molecule
        type(basis_set_info_t),intent(in)     :: ao_basis_1
        type(mo_basis),intent(in)             :: spinor_complement
        real(8),allocatable,intent(out)       :: Cadd(:,:,:)
        integer,allocatable,intent(out)       :: nadd_frag(:)

        type(basis_set_info_t),allocatable    :: ao_basis_2(:)
        type(mo_basis),allocatable            :: spinor_fragment(:)
        integer                               :: motype, i

        motype = spinor_complement%type
        if (n_fragments == 0) then
            n_fragments = molecule%natoms
            if (.not. file_for_MOs .eq. "adf" .and. .not. file_for_MOs .eq. "ams") atom_fragments = .true.
        end if
        call build_fragments(molecule,motype,ao_basis_2,spinor_fragment,nadd_frag)

        call get_correlating_space(ao_basis_1,ao_basis_2,spinor_complement,spinor_fragment,Cadd)

        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#      Correlating virtuals        #"
        write(luout,*) "####################################"
        write(luout,*)
        write(luout,'(A,I6)') " Available molecular virtuals for correlation: ",spinor_complement%nmo
        write(luout,'(A,I6)') " Basis functions intended for correlation:     ",size(Cadd,2)
        write(luout,*)
        do i = 1, n_fragments
           write(luout,'(//A,I3,A)') " Fragment",i,": correlating functions"
           call print_basis(luout,ao_basis_2(i))
        end do
        write(luout,*)

end subroutine generate_corrfuns

subroutine build_fragments(molecule,motype,ao_basis_2,spinor_fragment,nadd_frag)

        use rose_ao
        use rose_mo
        use rose_utils
        use read_xyzfile
        use rose_global
        use rose_input_data
        use cartesian_to_spherical
        use correlating_basisfunctions, only: read_correlation_basis

        implicit none

        type(molecule_info),intent(in)        :: molecule
        integer,intent(in)                    :: motype
        type(basis_set_info_t),allocatable,intent(out) :: ao_basis_2(:)
        type(mo_basis),allocatable,intent(out) :: spinor_fragment(:)
        integer,allocatable,intent(inout)     :: nadd_frag(:)

        integer                               :: i,j
        integer                               :: nxyz,nlm,ixyz,ilm
        real(8),allocatable                   :: coeff(:,:,:)

        allocate (ao_basis_2(n_fragments))
        allocate (spinor_fragment(n_fragments))

        ! Read the correlating functions of all fragments
        call read_correlation_basis(molecule,ao_basis_2)

        if (.not. allocated(nadd_frag)) allocate(nadd_frag(n_fragments))

        ! Define coefficients for the spherical combinations
        do i = 1, n_fragments
           ! AO basis is in Cartesians, but we want sphericals, first count how many there are
           ilm = 0
           do j=1,ao_basis_2(i)%nshells
              ilm = ilm + nfunctions(ao_basis_2(i)%gtos(j)%orb_momentum,2)
           enddo
           nadd_frag(i) = ilm
           call alloc_mo_basis (spinor_fragment(i),ao_basis_2(i)%nao,ilm,motype)
           allocate (coeff(ao_basis_2(i)%nao,ilm,rcq))
           ! Now get the spherical combinations and fill the coefficient matrix with these
           coeff = 0.D0
           ixyz = 0
           ilm = 0
           do j=1,ao_basis_2(i)%nshells
            nxyz = nfunctions(ao_basis_2(i)%gtos(j)%orb_momentum,1)
            nlm  = nfunctions(ao_basis_2(i)%gtos(j)%orb_momentum,2)
            call get_sphcom(ao_basis_2(i)%gtos(j)%orb_momentum,coeff(ixyz+1:ixyz+nxyz,ilm+1:ilm+nlm,1))
            ixyz = ixyz + nxyz
            ilm  = ilm  + nlm
           end do
           call fill_mo_basis (spinor_fragment(i),coeff=coeff)
           deallocate(coeff)
        end do

end subroutine build_fragments

subroutine get_correlating_space (ao_basis_1,ao_basis_2,spinor_complement,spinor_fragment,ICO)

        use rose_ao
        use rose_mo
        use get_overlap
        use rose_global
        use linear_algebra

        implicit none

        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis),intent(in)          :: spinor_complement  ! Canonical(ized) MO's of the virtual space to be projected upon
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO's representing the correlating virtuals
        real(8), allocatable,intent(out)   :: ICO(:,:,:) ! coefficient matrices to transform to projected virtual space

        integer                            :: nmo1, nmo2
        integer                            :: motype, rcq
        real(8), allocatable               :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        integer                            :: i

        nmo1 = spinor_complement%nmo
        nmo2 = sum(spinor_fragment(:)%nmo)
        motype  = spinor_complement%type

        if (motype < 3) then
           rcq = 1 ! real
        elseif (motype == 3) then
           rcq = 2 ! complex
        elseif (motype == 4) then
           rcq = 4 ! quaternion
        else
           stop 'unknown motype encountered'
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate (S11(nmo1,nmo1,rcq))
        allocate (S22(nmo2,nmo2,rcq))
        allocate (S12(nmo1,nmo2,rcq))
        allocate (S21(nmo2,nmo1,rcq))
        call compute_all_overlaps (ao_basis_1,ao_basis_2,spinor_complement,spinor_fragment,S11,S12,S21,S22)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Projection Matrix     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ICO = S12

        deallocate(S22,S21,S12)

        !-----Symmetrically orthonormalize the projected orbitals with S11-----
        call SymOrth(ICO,S11)

end subroutine get_correlating_space

end module rose_correlating_virtuals
