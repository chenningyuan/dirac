   subroutine rose_read_structure (molecule)

        use rose_ao
        use rose_mo
        use read_xyzfile
        use read_xyzfile, only: atom
        use hdf5_interface
        use fchk_interface
        use ADF_interface
        use rose_input_data
        use rose_global

        implicit none

        type(molecule_info), intent(inout)  :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..

        ! next two variables are not used, just present as the readers from adf and fchck cannot work with optional arguments yet
        type(basis_set_info_t) :: ao_basis_1 
        type(mo_basis)         :: spinor_original

        integer                            :: ierr

        ! Read geometry of the full molecule:
        if (file_for_MOs .eq. "dirac" .or. file_for_MOs .eq. "h5") then
           call read_from_hdf5 (molecule%name,ierr,molecule=molecule)
        else if (file_for_MOs .eq. "fchk") then
           call read_from_fchk (molecule%name,ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "ams") then
           molecule%name = 'adf'
           call read_in_ADF (ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "adf") then
           call read_from_ADF (molecule%name,ao_basis_1,spinor_original,molecule)
        end if
 
   end subroutine

   subroutine rose_print_structure (molecule)

        use rose_global
        use read_xyzfile
        use read_xyzfile, only: atom

        implicit none

        type(molecule_info), intent(in)  :: molecule
        type(atom), pointer              :: current_atom
        character(len=10)                :: atom_label
 
        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#      Molecular  Coordinates      #"
        write(luout,*) "####################################"
        write(luout,*)
        write(luout,'(A10,A6,3A12)') "Atom","Z","X","Y","Z"

        current_atom => molecule%firstatom
        do while(associated(current_atom))
           if (abs(current_atom%charge) < ghost_threshold .and. current_atom%bascharge > ghost_threshold) then
              atom_label = 'Gh(' // trim(current_atom%name) // ')'
           else
              atom_label = current_atom%name
           end if
           write(luout,'(A10,F6.2,3F12.8)') trim(atom_label),current_atom%bascharge,current_atom%x
           current_atom => current_atom%next
        end do

   end subroutine
