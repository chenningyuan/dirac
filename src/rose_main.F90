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
        integer                            :: nocc ! number of occupied spinors / orbitals
        integer                            :: nocc_active ! number of occupied spinors / orbitals in the active MO space
        integer                            :: ncore ! number of core spinors / orbitals
        real(8), allocatable               :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), allocatable               :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable               :: Cocc_IAO(:,:,:) ! occupied MOs in terms of IAOs
        real(8), allocatable               :: Cvir_IAO(:,:,:) ! valence virtual MOs in terms of IAOs
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
        write(luout,*) "@                                                @"
        write(luout,*) "@ Description of algorithms:                     @"
        write(luout,*) "@ - Senjean et al., DOI:10.1021/acs.jctc.0c00964 @"
        write(luout,*) "@ - Sen et al., DOI:10.1063/5.0137729            @"
        write(luout,*) "@ - Khandavilli et al., DOI:10.1063/5.0311698    @"
        write(luout,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(luout,*)


        ! Check the input (should already be read)
        call rose_check_input

        molecule%name = 'CHECKPOINT'
        
        call rose_read_structure  (molecule)
        call rose_print_structure (molecule)
        ! The MOs are read before the fragments are defined because the atomic fragments
        ! are enumerated over the centers of the supermolecular AO basis, not over the
        ! nucleus list. rose_read_mos does not depend on any fragment information.
        call rose_read_mos (molecule%name,file_for_MOs,ao_basis_1,spinor_original)
        call rose_define_fragments (molecule,fragments,ao_basis_1)

        need_iao = trim(occupied_mode) .eq. "localized" .or. &
                   (trim(virtual_mode) .ne. "original" .and. trim(virtual_source) .ne. "correlating")
        if (need_iao) then
           ! We need intrinsic fragment orbitals for either the occupied or the virtual space generation
           call rose_count_mospace (molecule,fragments,spinor_original,nocc_active,ncore,nmo_frag)
           call rose_split_mospace (ncore,spinor_original,spinor_frozen)
           call generate_iao (nocc_active,nmo_frag,fragments,ao_basis_1,spinor_original, &
                              IAO,ICO,Cocc_IAO,Cvir_IAO)
        else
           call rose_count_occupied(molecule,nocc)
           if (trim(virtual_source) .eq. "correlating") then
             ! No intrinsic orbitals required, but for uniform treatment of virtual orbital generation
             ! we do need to initialize the ICO array.
              call rose_split_mospace(nocc,spinor_original,spinor_frozen)
              nocc_active = 0
              allocate(ICO(spinor_original%nmo,spinor_original%nmo,rcq))
              ICO = 0.D0
              forall (i = 1:spinor_original%nmo) ICO(i,i,1) = 1.D0
           else
              ! Original occupied and original virtual orbitals: keep the full MO space,
              ! but initialize spinor_frozen for uniform downstream writing.
              call rose_split_mospace(0,spinor_original,spinor_frozen)
              nocc_active = nocc
           end if
        end if

        select case (trim(occupied_mode))
        case ("original")
           call copy_original_occ(nocc_active,spinor_original%energy,Locc,energies_occ,nocc_frag)
        case ("localized")
           call generate_ibos( 1,Cocc_IAO,IAO,spinor_original%energy,nmo_frag, &
                               Locc,energies_occ,nocc_frag)
           if (test) call test_IAO_and_IBO(IAO,Cocc_IAO,Locc)
        end select

        select case (trim(virtual_mode))
        case ("original")
           call copy_original_vir(nocc_active,spinor_original%energy,Lvir,energies_vir,nvir_frag)
        case ("localized","canonical")
           if (trim(virtual_source) .eq. "valence" .or. &
               trim(virtual_source) .eq. "valcorr" .or. &
               trim(virtual_source) .eq. "valhard") then
              ! The first set of virtuals to use is the set of valence virtuals
              call generate_ibos( nocc_active+1,Cvir_IAO,IAO,spinor_original%energy,nmo_frag, &
                                  Lvir,energies_vir,nvir_frag)
           end if
           if (trim(virtual_source) .eq. "correlating" .or. &
               trim(virtual_source) .eq. "valcorr") then
              ! The second set of virtuals to use is a set of designated correlating functions read from file
              call generate_corrfuns(molecule,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)
           else if (trim(virtual_source) .eq. "valhard") then
              ! The second set of virtuals to use is a set of newly generated intrinsic "hard" virtuals
              call generate_fragvirs(fragments,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)
           else
              ! No additional virtuals are added, just initialize for correct processing further on
              allocate(nadd_frag(n_fragments))
              nadd_frag = 0
           end if
           ! Merge available virtual sets and recanonize
           call rose_merge_virtuals(nvir_frag,nadd_frag,Lvir,Ladd,spinor_original%energy,energies_vir)
        end select

        if (need_iao) then
           ! We have intrinsic orbitals and potentially also avas orbitals to write as well
           if (avas_orbitals%array_type /= ' ') then
              call rose_generate_avas(ao_basis_1,spinor_original,spinor_frozen,IAO,nmo_frag,nocc_active)
           end if
           call rose_write_mos(nocc_frag,nvir_frag,Locc,Lvir,energies_occ,energies_vir, &
                               spinor_frozen,spinor_original,IAO)
        else
           ! We only write the produced occupied and virtual orbitals
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
