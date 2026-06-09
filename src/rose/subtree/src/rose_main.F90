   subroutine rose

        use rose_utils
        use rose_input_data
        use rose_global
        use linear_algebra
        use rose_ao
        use rose_mo
        use rose_intrinsic_mos
        use read_xyzfile
        use rose_import_fragments

        implicit none

        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info), allocatable   :: fragments(:) ! same information for the fragment
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis)                     :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis)                     :: spinor_original_core ! Core part of spinor_original
        type(mo_basis)                     :: spinor_iao ! MO-coefficients for the iao basis 
        type(mo_basis)                     :: spinor_ibo ! occupied MO-coefficients after localization 
        type(mo_basis)                     :: spinor_ico ! MO-coefficients for the space complementary to the iaos (and the core)
        type(mo_basis)                     :: spinor_ribo ! occupied MO-coefficients after localization and recanonization
        type(mo_basis)                     :: spinor_rido ! valence(+hard) virtual MO-coefficients after localization and recanonization
        type(mo_basis)                     :: spinor_avas ! MO coefficients for the AVAS approach to define an active space
        integer                            :: nocc ! number of occupied spinors / orbitals
        integer                            :: ncore ! number of core spinors / orbitals
        integer                            :: nvir_val ! number of localized valence virtual orbitals
        integer                            :: nmo4,nvir_add
        real(8), allocatable               :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), allocatable               :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable               :: Cocc_IAO(:,:,:),Cvir_IAO(:,:,:) ! orbital coefficient matrices in terms of IAOs
        real(8), allocatable               :: Cvir_IBO(:,:,:) ! Localized valence virtual MOs in terms of IAOs
        real(8), allocatable               :: Locc(:,:,:),Lvir(:,:,:),Ladd1(:,:,:)! Localized MOs in terms of original MOs of B1
        integer, allocatable               :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        integer, allocatable               :: nmo_frag2(:) ! Number of molecular orbitals (minimal complementary basis) for each fragment.
        integer, allocatable               :: nocc_frag(:)! Number of localized occupied orbitals assigned to each fragment

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
        call rose_define_fragments (molecule,fragments)
        call rose_read_mos ('CHECKPOINT',ao_basis_1,spinor_original,mo_type(restricted,spatial_orbitals),file_for_MOs)
        call rose_split_mospace (molecule,fragments,spinor_original,spinor_original_core,nocc,ncore)
        call generate_iao (nocc,fragments,ao_basis_1,spinor_original,spinor_original_core, &
                           spinor_ico,IAO,ICO,Cocc_IAO,Cvir_IAO,nmo_frag)
        call generate_ibo (nocc,nmo_frag,IAO,Cocc_IAO,Cvir_IAO,ao_basis_1,spinor_original, &
                            spinor_original_core,spinor_ribo,spinor_ico,ICO,Cvir_IBO,Locc,Lvir,nocc_frag,nvir_val)
        call generate_ico (fragments,ao_basis_1,spinor_original,spinor_ico,ICO,Ladd1,nmo_frag2,nmo4,nvir_add)
        call rose_write_mos (nocc,ncore,nmo_frag,nmo_frag2,nocc_frag,nvir_val,nmo4,nvir_add, &
                             IAO,Cvir_IBO,Locc,Lvir,Ladd1,ao_basis_1,spinor_original, &
                             spinor_original_core,spinor_ribo)
        if (add_correlating)   call project_corrfuns

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
