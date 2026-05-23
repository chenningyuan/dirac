subroutine project_corrfuns

        use rose_ao
        use rose_mo
        use rose_utils
        use linear_algebra
        use get_overlap
        use make_ibo
        use read_xyzfile 
        use rose_import_fragments
        use rose_input_data
        use rose_global
        use rose_intrinsic_mos
        use correlating_functions
        use cartesian_to_spherical
        use hdf5_interface
        use fchk_interface
        use ADF_interface
        use correlating_basisfunctions, only: read_correlation_basis


        implicit none

        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis)                     :: spinor_original ! Canonical MO-coefficients
        type(mo_basis)                     :: spinor_occupied ! Canonical occupied MO-coefficients
        type(mo_basis)                     :: spinor_ico ! MO-coefficients that span the space of atomic correlating functions (orthogonal to occupied)
        type(mo_basis), allocatable        :: spinor_fragment(:) ! MO-coefficients of the atomic correlating functions (non-orthogonal)
        integer                            :: ierr, nmo1, nao1, nmo2, nocc, rcq, i, j, motype, nxyz,nlm, ixyz, ilm
        integer, allocatable               :: nadd_frag(:)! Number of additional correlating virtual orbitals for each fragment
        real(8), allocatable               :: coeff(:,:,:)
 
             
        ! We first define the type of algebra (real/complex/quaternion) that will be used
        call define_algebra_type(spatial_orbitals,restricted,rcq)
        motype = mo_type(restricted,spatial_orbitals)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ INFORMATION FROM FILES (MO COEFF, BASIS SET)         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! For the full molecule:
        if (file_for_MOs .eq. "dirac") then
           call read_from_hdf5 ('CHECKPOINT',ierr,ao_basis_1,spinor_original,molecule)
           call convert_mo_basis(spinor_original,mo_type(restricted,spatial_orbitals))
        elseif (file_for_MOs .eq. "h5") then
           call read_from_hdf5 ('CHECKPOINT',ierr,ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "fchk") then
           call read_from_fchk ('MOLECULE',ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "ams") then
           molecule%name = 'adf'
           call read_in_ADF (ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "adf") then
           call read_from_ADF ("Molecule_B1_MO1",ao_basis_1,spinor_original,molecule)
        end if

        ! Keep only the virtual orbitals
        nmo1 = spinor_original%nmo
        nao1 = spinor_original%nao
        nocc = molecule%znumber - charge
        if (restricted) nocc = nocc/2
        nmo1 = nmo1 - nocc

        ! Make a copy of the occupied spinors, so that we can re-insert these later
        call copy_mo_basis (spinor_original,spinor_occupied)
        ! Remove the virtuals from the occupied set
        call reduce_mo_basis (spinor_occupied,nocc)
        ! Remove the occupied from the original set
        call reduce_mo_basis (spinor_original,nmo1,nocc+1)

        ! Define the list of atoms for which we want to add correlation functions
        if (n_fragments == 0) then
            n_fragments = molecule%natoms
            if (.not. file_for_MOs .eq. "adf" .and. .not. file_for_MOs .eq. "ams") atom_fragments = .true.
            allocate(nadd_frag(n_fragments))
            nadd_frag = -1 ! signal that this is not set by the user
        end if
        allocate (ao_basis_2(n_fragments))
        allocate (spinor_fragment(n_fragments))

        ! Read the correlating functions of all fragments
        call read_correlation_basis(molecule,ao_basis_2)

        ! Define coefficients for the spherical combinations
        nmo2 = 0
        do i = 1, n_fragments
           ! AO basis is in Cartesians, but we want sphericals, first count how many there are
           ilm = 0
           do j=1,ao_basis_2(i)%nshells
              ilm = ilm + nfunctions(ao_basis_2(i)%gtos(j)%orb_momentum,2)
           enddo
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
           nmo2 = nmo2 + ilm
        end do

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

        call get_projected_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,spinor_ico)

        ! Insert the occupied orbitals before these virtuals
        call insert_mo_basis (spinor_ico,spinor_occupied)

        if (trim(file_for_MOs) .eq. "dirac") then
           call write_mos_to_hdf5(spinor_ico)
        elseif (trim(file_for_MOs) .eq. "h5") then
           call write_mos_to_hdf5(spinor_ico,"CHECKPOINT","ico")
        else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
           call write_mos_to_fchk(spinor_ico,ao_basis_1,molecule,"ico")
        else if (trim(file_for_MOs) .eq. "adf") then
           call write_mos_to_ADF ("ico",spinor_ico)
        else
           print*,"**WARNING** error in the extension of file? I cannot write the resulting IAO/IBO to file."
        endif
        

         write(luout,'(//A)') " Projection complete: written orbitals to CHECKPOINT."

        contains
   
   subroutine get_projected_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,spinor_ico)

        use rose_ao
        use rose_mo
        use linear_algebra
        use rose_global

 
        implicit none

        type(mo_basis),intent(in)          :: spinor_original    ! Canonical MO-coefficients of the virtual space
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO-coefficients representing the correlating virtuals
        type(mo_basis),intent(out)         :: spinor_ico         ! Projected and recanonicalized MOs
        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments

        real(8), allocatable               :: ICO(:,:,:) ! coefficient matrices to transform from B1 to projected virtual space
        integer                            :: nmo1,nmo2,nmo3
        integer                            :: motype, rcq
        real(8), allocatable               :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        real(8), allocatable               :: energies_com(:) ! complementary orbitals
        integer                            :: nhard(1) ! Number of hard virtuals
        integer                            :: i, j
        real(8), allocatable               :: pnorm(:)

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
        allocate (ICO(nmo1,nmo2,rcq))

        !-----Symmetrically orthonormalize the projected orbitals with S11-----
        !call SymOrth(ICO,S11)
        call get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,1,ICO)

        nhard(1) = nmo2
        allocate(energies_com(nmo2))
        call recanonize (ICO,spinor_original%energy,energies_com,nhard)

        allocate(pnorm(nmo1))
        pnorm = 0.0
        write (luout,'(/A)') " Recanonization of projected virtuals: energies and quaternion coefficients for expansion in the original MOs."
        do i = 1, nmo2
          write(luout,'(/A,i7,A,f16.8)') "  index ",i," new energy: ",energies_com(i)
          do j = 1, nmo1
            if (dot_product(ICO(j,i,:),ICO(j,i,:))>1e-2) write (luout,'(i4,5x,f15.6,4f10.3)') j,spinor_original%energy(j),ICO(j,i,1:rcq)
            pnorm(j) = pnorm(j) + dot_product(ICO(j,i,:),ICO(j,i,:))
          end do
        enddo

        write (luout,'(/A)') " Contribution of original MOs: energies and norms."
        do j = 1, nmo1
           write (*,'(i4,5x,2f15.6)') j,spinor_original%energy(j),pnorm(j)
        end do

        ! Express everything in terms of the original AOs.
        call IMO_to_mobasis(ICO,spinor_original,spinor_ico,energies_com,nmo2)
        deallocate(ICO)

   end subroutine get_projected_mos

end subroutine project_corrfuns
