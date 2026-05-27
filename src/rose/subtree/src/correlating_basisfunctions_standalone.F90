submodule (correlating_basisfunctions) correlating_basisfunctions_standalone

contains

    module procedure read_correlation_basis

       use correlating_functions, only: read_corrfuns
       use read_xyzfile, only: atom

       implicit none

       type(atom), pointer :: current_atom
       integer             :: i

       current_atom => molecule%firstatom
       do i = 1, size(ao_basis)
          call read_corrfuns(nint(current_atom%bascharge),current_atom%x,(/0/),ao_basis(i), &
                             '/Users/visscher/sources/Dirac1/basis_dalton/cc-pVTZ')
          current_atom => current_atom%next
       end do

    end procedure read_correlation_basis

end submodule correlating_basisfunctions_standalone
