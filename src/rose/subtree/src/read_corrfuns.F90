module correlating_functions

   implicit none
   private

!  Linked list storage type for reading blocks of basis functions
   type basis_func_block
      integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
      integer:: n_primitives=0        ! number of primitives
      real, pointer    :: exponent(:) ! exponents
      integer, pointer :: designation(:) 
      type(basis_func_block), pointer :: next
   end type basis_func_block

   ! designation: 
   ! -3: inner core correlation
   ! -2: outer core correlation
   ! -1: valence correlation
   !  0: valence 
   !  1: diffuse 
   !  2: very diffuse

   public read_corrfuns

contains

   subroutine read_corrfuns (element,coord,nzeta,correlation,diffuse,ao_basis)

        use rose_ao
        use rose_utils

        implicit none

        integer, parameter  :: basis_angular = 1  ! Basis set is defined in terms of Cartesians for DIRAC / Respect
        integer, intent(in) :: element, nzeta, correlation, diffuse
        real(8), intent(in) :: coord(3)
        type(basis_set_info_t),intent(out) :: ao_basis ! the ao basis containing the correlating and diffuse functions
        type(basis_func_block), pointer :: first_block, current_block
        integer :: i,ierr, ishell, nshells, nao, file_unit,l_block,l_value
        character(len=10) basisfile

        ! Should include path information, for now assume this is copied into the working dir
        write (basisfile,'(a6,i1,a1)') 'dyall.',nzeta,'z'
        call get_free_fileunit(file_unit)
        open (file_unit,file=basisfile,form='formatted')

        ierr = find_element (file_unit,element) ! this routine (private to this module) is still to be written
        if (ierr /= 0) then
           write (*,*) 'Element ',element,' not found on file ',basisfile
           stop 'element not found'
        end if
        
        allocate(first_block)
        first_block%next => NULL()

        current_block => first_block
        do while (read_block(current_block) >= 0)  ! also the private read_block function is still to be written
             allocate(current_block%next)
             current_block => current_block%next
        end do

        ! count the number of shells that we will retain and the number of correlating+diffuse functions
        nshells = 0
        nao     = 0
        current_block => first_block
        do while (associated(current_block))
           do i = 1, current_block%n_primitives
              if ( (current_block%designation(i) < 0 .AND. current_block%designation(i) >= correlation) .OR. &
                   (current_block%designation(i) > 0 .AND. current_block%designation(i) <= diffuse) ) then
                 nshells = nshells + 1
                 nao = nao + nfunctions(current_block%orb_momentum,1)
              end if
           end do
           current_block => current_block%next
        end do
        
        ! allocate basis set and fill in the exponents
        ao_basis%nshells       = nshells
        ao_basis%nao           = nao
        ao_basis%basis_angular = basis_angular
        allocate(ao_basis%gtos(nshells))
        ishell = 0
        nao = 0
        current_block => first_block
        do while (associated(current_block))
           do i = 1, current_block%n_primitives
              if ( (current_block%designation(i) < 0 .AND. current_block%designation(i) >= correlation) .OR. &
                   (current_block%designation(i) > 0 .AND. current_block%designation(i) <= diffuse) ) then
                 ishell = ishell + 1
                 ao_basis%gtos(ishell)%orb_momentum = current_block%orb_momentum
                 ao_basis%gtos(ishell)%atom_element = element
                 ao_basis%gtos(ishell)%coord = coord
                 ao_basis%gtos(ishell)%atom_number  = 0      ! initialize to zero, can be set later if we want to order atoms
                 ao_basis%gtos(ishell)%n_primitives = 1      ! basis functions are uncontracted
                 allocate (ao_basis%gtos(ishell)%exponent(1))
                 allocate (ao_basis%gtos(ishell)%coefficient(1))
                 ao_basis%gtos(ishell)%exponent(1) = current_block%exponent(i)
                 ao_basis%gtos(ishell)%coefficient(1) = 1.D0
              end if
           end do
           current_block => current_block%next
        end do

        ! Everything is read, end by deallocating the memory used for the blocks

   end subroutine read_corrfuns

   integer function find_element (file_unit,element)

       integer,intent(in) :: file_unit,element

       find_element = -1 ! not found

   end function find_element

   integer function read_block(current_block)

      type(basis_func_block), pointer :: current_block

      read_block = -1 ! not read

   end function read_block

end module
