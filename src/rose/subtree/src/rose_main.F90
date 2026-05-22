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

        if (localize_occupied .or. localize_virtuals) call generate_ibo(molecule,fragments)
        if (add_correlating) call project_corrfuns

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
