submodule (correlating_basisfunctions) correlating_basisfunctions_dirac

contains

    module procedure read_correlation_basis

       use correlating_functions, only: read_corrfuns
       use hermodules, only: atomtype, basistype, firstatomtype
       use read_xyzfile, only: atom

       implicit none

       type(atom), pointer      :: current_atom
       type(basistype), pointer :: baspointer
       type(atomtype), pointer  :: typepointer
       integer                  :: i, j

       j = 0
       typepointer=>firstatomtype
       do while (associated(typepointer))
          baspointer=>typepointer%basis
          current_atom => typepointer%first
          do i = 1, typepointer%natom
             j = j + 1
             call read_corrfuns(nint(current_atom%bascharge),current_atom%x,(/0/),ao_basis(j), &
                                trim(baspointer%name))
             current_atom => current_atom%next
          end do
          typepointer=>typepointer%next
       end do

    end procedure read_correlation_basis

end submodule correlating_basisfunctions_dirac
