   program project_corrfuns

        use rose_ao
        use rose_mo
        use rose_utils
        use correlating_functions
        use linear_algebra
        use read_xyz
        use get_overlap
        use import_fragments
        use hdf5_interface

        implicit none

        complex(8), parameter              :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the correlation functions
        type(molecule_info)                :: molecule ! store informations : charge, nelec, nalpha, nbeta, multiplicit. Used for fchk file..
        type(atom), pointer                :: first_atom,current_atom ! linked list to store coordindate information
        type(mo_basis)                     :: cspinor_original ! MO-coefficients before localization 
        type(mo_basis),allocatable         :: cspinor_fragment(:) ! MO-coefficients of the correlation functions
        type(mo_basis)                     :: cspinor_projected   ! MO-coefficients after projecting the virtuals on the correlation space
        integer                            :: luxyz, ierr, i, j, k, l
        integer                            :: nmo1, nmo2, imo, jmo, nmoi, nmoj, nao1, ncore, nvale, rcq
        integer                            :: n_fragments ! the number of fragments used to project on to
        integer                            :: nocc ! Number of occupied spinors
        integer                            :: nvir ! Number of virtual spinors
        integer                            :: charge ! Charge of the molecule
        integer                            :: restrict ! restrict = 1 means that the logical variable "restricted" will be set to .true.
        real(8), allocatable               :: S_block(:,:,:) ! sub block of the overlap matrix
        real(8), allocatable               :: S_block_bar(:,:,:) ! sub block of the overlap matrix, transpose conjugate.
        real(8), allocatable               :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) !  the overlap matrices
        real(8), allocatable               :: P12(:,:,:),P21(:,:,:)         !  the projection operators
        real(8), allocatable               :: E(:,:,:)                    !  the projected virtuals
        real(8), allocatable               :: E_vec(:,:,:)                !  eigenvectors after projection
        real(8), allocatable               :: E_val(:)                  !  eigenvalues after projection
        character(len=18)                  :: filename_projected        !  name of the file to write the projected virtuals on to
        logical                            :: restricted = .true.       !  maintain (Kramers)-restriction for the localized orbital set
        logical                            :: spatial_orbitals = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ Input Files         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_free_fileunit(luxyz)
        open (luxyz,file="INPUT",form='formatted')
        read(luxyz,'(I2)') restrict
        close(luxyz)

        if (restrict .ne. 1) restricted = .false.

        call define_algebra_type(spatial_orbitals,restricted,rcq)

        call get_free_fileunit(luxyz)
        open (luxyz,file='MOLECULE.XYZ',form='formatted')
        call get_atoms(luxyz,first_atom)
        close (luxyz,status='keep')

        call read_from_hdf5 ('MOLECULE',ierr,ao_basis_1,cspinor_original,molecule)
        call convert_mo_basis (cspinor_original,3)
        nmo1 = cspinor_original%nmo
        nao1 = cspinor_original%nao

        current_atom => first_atom
        print*,"Molecular structure taken from the MOLECULE.XYZ file"
        ! For now we will just consider atomic fragments, which should be available from a library of DFCOEFs
        n_fragments = 0
        nocc = 0
        do while(associated(current_atom))
           n_fragments = n_fragments + 1
           print*,current_atom%name,current_atom%bascharge,current_atom%x
           current_atom=>current_atom%next
        end do
        print*,"Locations of basis functions"
        print*, ao_basis_1%gtos(1)%coord
        do i = 2, ao_basis_1%nshells
           if (ao_basis_1%gtos(i)%atom_number /= ao_basis_1%gtos(i-1)%atom_number) print*, ao_basis_1%gtos(i)%coord
        end do
        allocate (ao_basis_2(n_fragments))
        allocate (cspinor_fragment(n_fragments))
        nmo2 = 0
        current_atom => first_atom
        do i = 1, n_fragments
           ! The call below should give back the ao basis consisting of just the correlating functions
           !call read_corrfuns(nint(current_atom%bascharge),current_atom%x,2,-1,0,ao_basis_2(i))
           ! This is for the moment replaced by a call that retrives the full basis by reading from a dfcoef file
           ! the also read cspinor_fragment structure is overwritten by init_fragment to be able to check that routine
           call get_atomic_fragment (nint(current_atom%bascharge),current_atom%x,ao_basis_2(i),cspinor_fragment(i),&
                                     mo_type(restricted,spatial_orbitals),file_for_MOs="dfcoef")
           ! Determine the size of the (mo) basis that should be kept
           call atom_minbas_size (nint(current_atom%bascharge),ncore,nvale)
           ! Keep only the spinors that we need
           call reduce_mo_basis (cspinor_fragment(i),(ncore+nvale)/2)
           call init_fragment (ao_basis_2(i),cspinor_fragment(i))
           call convert_mo_basis (cspinor_fragment(i),3)
           nmo2 = nmo2 + cspinor_fragment(i)%nmo
           current_atom=>current_atom%next
        end do
        print*,"Size of original and projection basis ",nmo1,nmo2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !-----Compute S11----- (nmo1 x nmo1) (Should be equal to Identity)
        allocate (S11(nmo1,nmo1,rcq))
        call get_mo_overlap (ao_basis_1,cspinor_original,ao_basis_1,cspinor_original,S11)

        !-----Compute S22----- (nmo2 x nmo2)
        allocate (S22(nmo2,nmo2,rcq))
        imo = 0
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           jmo = imo
           do j = i, n_fragments
              nmoj = cspinor_fragment(j)%nmo
              allocate (S_block(nmoi,nmoj,rcq))
              call get_mo_overlap (ao_basis_2(i),cspinor_fragment(i),ao_basis_2(j),cspinor_fragment(j),S_block)
              S22(imo+1:imo+nmoi,jmo+1:jmo+nmoj,:) = S_block
              if (i /= j) then
                 allocate (S_block_bar(nmoj,nmoi,rcq))
                 call transpose_conjg(S_block,S_block_bar)
                 S22(jmo+1:jmo+nmoj,imo+1:imo+nmoi,:) = S_block_bar
                 deallocate(S_block_bar)
              end if
              jmo = jmo + nmoj
              deallocate (S_block)
           enddo
           imo = imo + nmoi
        enddo

        !-----Compute S21----- (nmo2 x nmo1)
        allocate (S21(nmo2,nmo1,rcq))
        imo = 0
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           allocate (S_block(nmoi,nmo1,rcq))
           call get_mo_overlap (ao_basis_2(i),cspinor_fragment(i),ao_basis_1,cspinor_original,S_block)
           S21(imo+1:imo+nmoi,:,:) = S_block
           imo = imo + nmoi
           deallocate (S_block)
        enddo

        !-----Compute S12----- (nmo1 x nmo2)
        allocate (S12(nmo1,nmo2,rcq))
        call transpose_conjg(S21,S12)

        allocate (P21(nmo2,nmo1,rcq))
        call solve_lineq(S22,S21,P21)
        allocate (E(nmo1,nmo1,rcq))
        call matrix_mul(S12,P21,E)
        allocate (E_vec(nmo1,nmo1,rcq))
        allocate (E_val(nmo1))
        E_vec = E
        call diag_matrix (E_vec,E_val)
        do i = 1, nmo1
           print*,i,E_val(i)
        end do
        print*,"Trace of E ",sum(E_val)

   end program project_corrfuns
