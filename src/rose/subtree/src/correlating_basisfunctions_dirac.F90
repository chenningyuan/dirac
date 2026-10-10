submodule (correlating_basisfunctions) correlating_basisfunctions_dirac

contains

    module procedure read_correlation_basis

       use correlating_functions, only: read_corrfuns
       use hermodules, only: atomtype, basistype, firstatomtype
       use read_xyzfile, only: atom
       use rose_basis_locations, only: get_basis_file_location
       use rose_input_data, only: get_corr_levels

       implicit none

       type(atom), pointer      :: current_atom
       type(basistype), pointer :: baspointer
       type(atomtype), pointer  :: typepointer
       character(len=32)        :: basis_name
       character(len=300)       :: basis_file
       logical                  :: found_basis_file
       integer                  :: i, j

       j = 0
       typepointer=>firstatomtype
       do while (associated(typepointer))
          baspointer=>typepointer%basis
          basis_name = trim(baspointer%file)
          if (len_trim(basis_name) == 0) basis_name = trim(baspointer%name)
          current_atom => typepointer%first
          do i = 1, typepointer%natom
             j = j + 1
             call get_basis_file_location(trim(basis_name),basis_file,found_basis_file, &
                                          basis_role='correlating',fallback_role='orbital')
             if (.not. found_basis_file) then
                print*, "ROSE could not find the registered basis file location for basis ", &
                        trim(basis_name)
                stop
             end if
             call read_corrfuns(nint(current_atom%bascharge),current_atom%x,get_corr_levels(current_atom%name),ao_basis(j), &
                                trim(basis_file))
             current_atom => current_atom%next
          end do
          typepointer=>typepointer%next
       end do

    end procedure read_correlation_basis

end submodule correlating_basisfunctions_dirac
