module rose_correlating_virtuals

        private

        public :: project_corrfuns
        public :: generate_corrfuns

contains

subroutine project_corrfuns(molecule,ao_basis_1,spinor_original,spinor_frozen)

        use rose_ao
        use rose_mo
        use rose_utils
        use read_xyzfile 
        use rose_input_data
        use rose_global
        use make_ibo
        use hdf5_interface
        use fchk_interface
        use ADF_interface


        implicit none

        type(molecule_info),intent(in)     :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(in)          :: spinor_original ! Canonical virtual MO-coefficients
        type(mo_basis),intent(in)          :: spinor_frozen ! Frozen occupied MO-coefficients
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis)                     :: spinor_ico ! MO-coefficients that span the space of atomic correlating functions (orthogonal to occupied)
        type(mo_basis), allocatable        :: spinor_fragment(:) ! MO-coefficients of the atomic correlating functions (non-orthogonal)
        real(8), allocatable               :: ICO(:,:,:) ! coefficient matrices to transform from B1 to projected virtual space
        real(8), allocatable               :: energies_com(:) ! complementary orbitals
        integer                            :: nmo1, nao1, nmo2, nocc, i, motype
        integer                            :: nhard(1) ! Number of hard virtuals
        integer, allocatable               :: nadd_frag(:)! Number of additional correlating virtual orbitals for each fragment
 
        motype = spinor_original%type
        nmo1 = spinor_original%nmo
        nao1 = spinor_original%nao
        nocc = spinor_frozen%nmo

        ! Define the list of atoms for which we want to add correlation functions
        if (n_fragments == 0) then
            n_fragments = molecule%natoms
            if (.not. file_for_MOs .eq. "adf" .and. .not. file_for_MOs .eq. "ams") atom_fragments = .true.
            allocate(nadd_frag(n_fragments))
            nadd_frag = -1 ! signal that this is not set by the user
        end if

        call build_fragments(molecule,motype,ao_basis_2,spinor_fragment,nadd_frag)
        nmo2 = sum(spinor_fragment(:)%nmo)

        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#         Basis set size           #"
        write(luout,*) "#               and                #"
        write(luout,*) "#      Canonical MO Energies       #"
        write(luout,*) "####################################"
        write(luout,*)
        write(luout,'(A,I6)') " No. basis functions:       ",nao1
        write(luout,'(A,I6)') " No. occupied  orbitals:    ",nocc
        write(luout,'(A,I6)') " No. virtual orbitals:      ",nmo1
        write(luout,'(A,I6)') " No. correlating functions: ",nmo2
        write(luout,'(//A)') " Details of the primary basis:"
        call print_basis(luout,ao_basis_1)
        do i = 1, n_fragments
           write(luout,'(//A,I3,A)') " Fragment",i,": correlating functions"
           call print_basis(luout,ao_basis_2(i))
        end do
        write(luout,*)

        call get_correlating_space(ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,ICO)
        nhard(1) = size(ICO,2)
        allocate(energies_com(size(ICO,2)))
        call recanonize(ICO,spinor_original%energy,energies_com,nhard)

        write(luout,*)
        if (trim(virtual_mode) .eq. "canonical") then
           write(luout,*)"----- Globally recanonicalized virtual energies -----"
           write(luout,'(2A12,A20)')"Fragment","MO","Energy"
           do i=1,nhard(1)
             write(luout,'(2I12,F20.9)') 1,i, energies_com(i)
           end do
        else
           error stop 'local recanonicalization not yet implemented for this pathway'
        end if
 
        call IMO_to_mobasis(ICO,spinor_original,spinor_ico,energies_com,size(ICO,2))
        deallocate(ICO,energies_com)

        ! Insert the occupied orbitals before these virtuals
        call insert_mo_basis (spinor_ico,spinor_frozen)

        if (trim(file_for_MOs) .eq. "dirac") then
           call write_mos_to_hdf5(spinor_ico)
        elseif (trim(file_for_MOs) .eq. "h5") then
           call write_mos_to_hdf5(spinor_ico,"CHECKPOINT","ico")
        else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
           call write_mos_to_fchk(spinor_ico,"CHECKPOINT","ico")
        else if (trim(file_for_MOs) .eq. "adf") then
           call write_mos_to_ADF ("ico",spinor_ico)
        else
           print*,"**WARNING** error in the extension of file? I cannot write the resulting IAO/IBO to file."
        endif
        
        write(luout,'(//A)') " Projection complete: written orbitals to CHECKPOINT."

end subroutine project_corrfuns

subroutine generate_corrfuns(molecule,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)

        use rose_ao
        use rose_mo
        use linear_algebra
        use read_xyzfile
        use rose_input_data
        use rose_global

        implicit none

        type(molecule_info),intent(in)        :: molecule
        type(basis_set_info_t),intent(in)     :: ao_basis_1
        type(mo_basis),intent(in)             :: spinor_original
        real(8),intent(inout)                 :: ICO(:,:,:)
        real(8),allocatable,intent(out)       :: Ladd(:,:,:)
        integer,allocatable,intent(out)       :: nadd_frag(:)

        type(basis_set_info_t),allocatable    :: ao_basis_2(:)
        type(mo_basis),allocatable            :: spinor_fragment(:)
        type(mo_basis)                        :: spinor_complement
        integer                               :: motype, i
        real(8),allocatable                   :: Cadd(:,:,:)

        motype = spinor_original%type
        call build_fragments(molecule,motype,ao_basis_2,spinor_fragment,nadd_frag)

        call IMO_to_mobasis(ICO,spinor_original,spinor_complement)

        call get_correlating_space(ao_basis_1,ao_basis_2,spinor_complement,spinor_fragment,Cadd)
        allocate(Ladd(size(ICO,1),size(Cadd,2),size(ICO,3)))
        call matrix_mul(ICO,Cadd,Ladd)

        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#      Correlating virtuals        #"
        write(luout,*) "####################################"
        write(luout,*)
        write(luout,'(A,I6)') " Available molecular virtuals for correlation: ",size(ICO,2)
        write(luout,'(A,I6)') " Basis functions intended for correlation:     ",size(Cadd,2)
        write(luout,*)
        do i = 1, n_fragments
           write(luout,'(//A,I3,A)') " Fragment",i,": correlating functions"
           call print_basis(luout,ao_basis_2(i))
        end do
        write(luout,*)

        deallocate(Cadd)

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

subroutine get_correlating_space (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,ICO)

        use rose_ao
        use rose_mo
        use get_overlap
        use rose_global
        use linear_algebra

        implicit none

        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis),intent(in)          :: spinor_original    ! Canonical MO-coefficients of the virtual space
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO-coefficients representing the correlating virtuals
        real(8), allocatable,intent(out)   :: ICO(:,:,:) ! coefficient matrices to transform from B1 to projected virtual space

        integer                            :: nmo1, nmo2
        integer                            :: motype, rcq
        real(8), allocatable               :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        integer                            :: i

        nmo1 = spinor_original%nmo
        nmo2 = sum(spinor_fragment(:)%nmo)
        motype  = spinor_original%type

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
        call compute_all_overlaps (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,S11,S12,S21,S22)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Projection Matrix     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ICO = S12

        deallocate(S22,S21,S12)

        !-----Symmetrically orthonormalize the projected orbitals with S11-----
        call SymOrth(ICO,S11)

end subroutine get_correlating_space

end module rose_correlating_virtuals
