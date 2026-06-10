   subroutine rose

        use rose_utils
        use rose_input_data
        use rose_global
        use linear_algebra
        use rose_ao
        use rose_mo
        use rose_intrinsic_mos
        use rose_correlating_virtuals, only: generate_corrfuns, project_corrfuns
        use read_xyzfile
        use rose_import_fragments

        implicit none

        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info), allocatable   :: fragments(:) ! same information for the fragment
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis)                     :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis)                     :: spinor_original_core ! Core part of spinor_original
        type(mo_basis)                     :: spinor_iao ! MO-coefficients for the iao basis 
        type(mo_basis)                     :: spinor_rido ! valence(+hard) virtual MO-coefficients after localization and recanonization
        integer                            :: nocc ! number of occupied spinors / orbitals
        integer                            :: ncore ! number of core spinors / orbitals
        real(8), allocatable               :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), allocatable               :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable               :: Locc(:,:,:),Lvir(:,:,:),Ladd(:,:,:)! Localized MOs in terms of original MOs of B1
        real(8), allocatable               :: energies_occ(:) ! Recanonicalized occupied energies
        real(8), allocatable               :: energies_vir(:) ! Recanonicalized virtual energies
        integer, allocatable               :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        integer, allocatable               :: nadd_frag(:) ! Number of molecular orbitals (minimal complementary basis) for each fragment.
        integer, allocatable               :: nocc_frag(:)! Number of localized occupied orbitals assigned to each fragment
        integer, allocatable               :: nvir_frag(:)! Number of localized valence virtual orbitals assigned to each fragment

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

        if (localize_orbitals) then
           call rose_define_fragments (molecule,fragments)
           call rose_count_mospace (molecule,fragments,spinor_original,nocc,ncore,nmo_frag)
           call rose_split_mospace (ncore,spinor_original,spinor_original_core)
           call generate_iao (nocc,nmo_frag,fragments,ao_basis_1,spinor_original,IAO,ICO)
           call IMO_to_mobasis(IAO,spinor_original,spinor_iao)
           if (avas_orbitals%array_type /= ' ') then
              call rose_generate_avas(ao_basis_1,spinor_original,spinor_original_core,spinor_iao,nmo_frag,nocc)
           end if
           call generate_ibo   (nocc,nmo_frag,IAO,spinor_original%energy, &
                                Locc,Lvir,energies_occ,nocc_frag,nvir_frag)
           if (add_correlating) then
              call generate_corrfuns(molecule,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)
           else
              call generate_fragvirs(fragments,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)
           end if
           call rose_merge_virtuals(nvir_frag,nadd_frag,Lvir,Ladd,spinor_original%energy,energies_vir)
           call rose_write_mos(nocc_frag,nvir_frag,Locc,Lvir,energies_occ,energies_vir, &
                               spinor_original_core,spinor_original,spinor_iao)
        elseif (add_correlating) then
           call project_corrfuns
        else
           error stop "Non-localized ROSE runs currently require correlating functions."
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
