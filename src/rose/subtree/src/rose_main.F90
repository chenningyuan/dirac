   subroutine rose

        use rose_utils
        use rose_input_data
        use rose_global
        use linear_algebra
        use rose_ao
        use rose_mo
        use rose_intrinsic_mos
        use rose_correlating_virtuals, only: generate_corrfuns
        use read_xyzfile
        use rose_import_fragments
        use check

        implicit none

        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info), allocatable   :: fragments(:) ! same information for the fragment
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis)                     :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis)                     :: spinor_frozen ! Frozen part of spinor_original
        type(mo_basis)                     :: spinor_iao ! MO-coefficients for the iao basis 
        type(mo_basis)                     :: spinor_complement ! MO-coefficients spanning the virtual complement
        integer                            :: nocc ! number of occupied spinors / orbitals
        integer                            :: ncore ! number of core spinors / orbitals
        real(8), allocatable               :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), allocatable               :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable               :: Cocc_IAO(:,:,:) ! occupied MOs in terms of IAOs
        real(8), allocatable               :: Cvir_IAO(:,:,:) ! valence virtual MOs in terms of IAOs
        real(8), allocatable               :: Cadd(:,:,:) ! correlating virtuals in terms of the current complement
        real(8), allocatable               :: Locc(:,:,:),Lvir(:,:,:),Ladd(:,:,:)! Localized MOs in terms of original MOs of B1
        real(8), allocatable               :: energies_occ(:) ! Recanonicalized occupied energies
        real(8), allocatable               :: energies_vir(:) ! Recanonicalized virtual energies
        integer, allocatable               :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        integer, allocatable               :: nadd_frag(:) ! Number of molecular orbitals (minimal complementary basis) for each fragment.
        integer, allocatable               :: nocc_frag(:)! Number of localized occupied orbitals assigned to each fragment
        integer, allocatable               :: nvir_frag(:)! Number of localized valence virtual orbitals assigned to each fragment
        integer                            :: i
        logical                            :: need_iao

        write(luout,*)
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*) "@                                                @"
        write(luout,*) "@       Reduction of Orbital Space Extent        @"
        write(luout,*) "@                     (ROSE)                     @"
        write(luout,*) "@                                                @"
        write(luout,*) "@ By Bruno Senjean, Souloke Sen, Michal Repisky, @"
        write(luout,*) "@        Gerald Knizia and Lucas Visscher.       @"
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*)


        ! Check the input (should already be read)
        call rose_check_input

        molecule%name = "CHECKPOINT"
        
        call rose_read_structure  (molecule)
        call rose_print_structure (molecule)
        call rose_read_mos ('CHECKPOINT',file_for_MOs,mo_type(restricted,spatial_orbitals),ao_basis_1,spinor_original)

        need_iao = trim(occupied_mode) .eq. "localized" .or. &
                   (trim(virtual_mode) .ne. "original" .and. trim(virtual_source) .ne. "correlating")
        if (need_iao) then
           call rose_define_fragments (molecule,fragments)
           call rose_count_mospace (molecule,fragments,spinor_original,nocc,ncore,nmo_frag)
           call rose_split_mospace (ncore,spinor_original,spinor_frozen)
           call generate_iao (nocc,nmo_frag,fragments,ao_basis_1,spinor_original,IAO,ICO)
           call express_in_intrinsic_mos(IAO,nocc,Cocc_IAO,Cvir_IAO)
           call IMO_to_mobasis(IAO,spinor_original,spinor_iao)
           if (avas_orbitals%array_type /= ' ') then
              call rose_generate_avas(ao_basis_1,spinor_original,spinor_frozen,spinor_iao,nmo_frag,nocc)
           end if
        else
           call rose_count_occupied(molecule,nocc)
           if (trim(virtual_source) .eq. "correlating") then
              if (n_fragments == 0) then
                 n_fragments = molecule%natoms
                 if (.not. file_for_MOs .eq. "adf" .and. .not. file_for_MOs .eq. "ams") atom_fragments = .true.
              end if
              call rose_split_mospace(nocc,spinor_original,spinor_frozen)
              allocate(ICO(spinor_original%nmo,spinor_original%nmo,rcq))
              ICO = 0.D0
              do i = 1, spinor_original%nmo
                 ICO(i,i,1) = 1.D0
              end do
              nocc = 0
           else
              call rose_split_mospace(0,spinor_original,spinor_frozen)
           end if
        end if

        select case (trim(occupied_mode))
        case ("original")
           call copy_original_occ(nocc,spinor_original%energy,Locc,energies_occ,nocc_frag)
        case ("localized")
           call generate_ibos( 1,Cocc_IAO,IAO,spinor_original%energy,nmo_frag, &
                               Locc,energies_occ,nocc_frag)
           if (test) call test_IAO_and_IBO(IAO,Cocc_IAO,Locc)
        end select

        select case (trim(virtual_mode))
        case ("original")
           call copy_original_vir(nocc,spinor_original%energy,Lvir,energies_vir,nvir_frag)
        case ("localized","canonical")
           if (trim(virtual_source) .eq. "valence" .or. &
               trim(virtual_source) .eq. "valcorr" .or. &
               trim(virtual_source) .eq. "valhard") then
              call generate_ibos( nocc+1,Cvir_IAO,IAO,spinor_original%energy,nmo_frag, &
                                  Lvir,energies_vir,nvir_frag)
           end if
           if (trim(virtual_source) .eq. "correlating" .or. &
               trim(virtual_source) .eq. "valcorr") then
              if (need_iao) then
                 call IMO_to_mobasis(ICO,spinor_original,spinor_complement)
              else
                 call copy_mo_basis(spinor_original,spinor_complement)
              end if
              call generate_corrfuns(molecule,ao_basis_1,spinor_complement,Cadd,nadd_frag)
              allocate(Ladd(size(ICO,1),size(Cadd,2),size(ICO,3)))
              call matrix_mul(ICO,Cadd,Ladd)
              deallocate(Cadd)
           else if (trim(virtual_source) .eq. "valhard") then
              call generate_fragvirs(fragments,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)
           else
              allocate(nadd_frag(n_fragments))
              nadd_frag = 0
           end if
           call rose_merge_virtuals(nvir_frag,nadd_frag,Lvir,Ladd,spinor_original%energy,energies_vir)
        end select

        if (need_iao) then
           call rose_write_mos(nocc_frag,nvir_frag,Locc,Lvir,energies_occ,energies_vir, &
                               spinor_frozen,spinor_original,spinor_iao)
        else
           call rose_write_mos(nocc_frag,nvir_frag,Locc,Lvir,energies_occ,energies_vir, &
                               spinor_frozen,spinor_original)
        end if
        call delete_adf_objects(ao_basis_1)

!!!!!!!!!!!!!!!!!!!!
!       END        !
!!!!!!!!!!!!!!!!!!!!
        write(luout,*)
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*) "@                                                @"
        write(luout,*) "@       Reduction of Orbital Space Extent        @"
        write(luout,*) "@                                                @"
        write(luout,*) "@                      END                       @"
        write(luout,*) "@                                                @"
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

   end
