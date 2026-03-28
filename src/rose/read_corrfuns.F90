module correlating_functions

   implicit none
   private

   type, public :: corr_func
      integer :: atomic_number ! atomic number of the atom
      integer :: nlevels ! number of levels to consider
      character(30) :: basis_name ! basis name for the correlation functions
      integer, pointer :: levels(:) ! levels to consider
   end type corr_func

!  Linked list storage type for reading blocks of basis functions
   type basis_func_block
      integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
      integer:: n_primitives=0        ! number of primitives
      integer:: designation
      real(8), pointer :: exponent(:) ! exponents
      type(basis_func_block), pointer :: next
   end type basis_func_block

   ! designation, Levels:
   ! -2: second shell below valence
   ! -1: shell below valence
   !  0: valence correlation
   !  1: first diffuse augmentation shell
   

   public read_corrfuns

contains

   subroutine read_corrfuns (atomic_number,coord,levels,ao_basis,basisfile)

        use rose_ao
        use rose_utils

        implicit none

        integer, parameter  :: basis_angular = 1  ! Basis set is defined in terms of Cartesians for DIRAC / Respect
        integer, intent(in) :: atomic_number
        integer, intent(in) :: levels(:) ! assumed shape array
        real(8), intent(in) :: coord(3)
        character(len=*), intent(in) :: basisfile
        type(basis_set_info_t),intent(out) :: ao_basis ! the ao basis containing the correlating and diffuse functions
        type(basis_func_block), pointer :: first_block, current_block
        integer :: i,ierr, ishell, nshells, nao, file_unit

        ! Should include path information, for now assume this is copied into the working dir
        call get_free_fileunit(file_unit)
        open (file_unit,file=basisfile,form='formatted')

        ierr = find_element (file_unit,atomic_number)
        if (ierr /= 0) then
           write(*,*) 'Element with atomic number ', atomic_number, ' not found in file ',basisfile
           stop
        end if

        allocate(first_block)
        first_block%next => NULL()

        current_block => first_block
        do while ( read_block(file_unit,current_block) >= 0)  ! read_block read severals blocks of a given level
        end do

        ! count the number of shells that we will retain and the number of basis functions
        nshells = 0
        nao     = 0
        current_block => first_block
        do while (associated(current_block))
           do i = 1, current_block%n_primitives
              if (ANY (levels == current_block%designation )) then
                 nshells = nshells + 1
                 nao = nao + nfunctions(current_block%orb_momentum,1)
              end if
           end do
           current_block => current_block%next
        end do

        allocate (ao_basis%shell_indices(nao))

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
              if (ANY (levels == current_block%designation )) then
                 ishell = ishell + 1
                 ao_basis%gtos(ishell)%orb_momentum = current_block%orb_momentum
                 ao_basis%gtos(ishell)%atom_element = atomic_number
                 ao_basis%gtos(ishell)%coord = coord
                 ao_basis%gtos(ishell)%atom_number  = 0 ! this is set to 0 for now but can be generalized.
                 ao_basis%gtos(ishell)%n_primitives = 1      ! basis functions are uncontracted
                 allocate (ao_basis%gtos(ishell)%exponent(1))
                 allocate (ao_basis%gtos(ishell)%coefficient(1))
                 ao_basis%gtos(ishell)%exponent(1) = current_block%exponent(i)
                 ao_basis%gtos(ishell)%coefficient(1) = &
                 & normalization_constant_psi4(current_block%orb_momentum-1,current_block%exponent(i))
                 ao_basis%shell_indices(nao+1:nao+nfunctions(current_block%orb_momentum,1)) = ishell
                 nao = nao + nfunctions(current_block%orb_momentum,1)
              end if
           end do
           current_block => current_block%next
        end do

        close(file_unit)

        ! Everything is read, end by deallocating the memory used for the blocks

   end subroutine read_corrfuns

      integer function find_element ( file_unit, atom_number )
      !=================================================================
      !  Function to check if the atomic element is present in the file
      !=================================================================
         integer         :: file_unit
         integer         :: atom_number
         integer         :: buff_1
         character(40)   :: buffer
         character(7)    :: buffer_string1,buffer_string2
         integer         :: l, element, stat

         find_element = -1 ! Case when we do not find the element
         l = 0
         do ! Infinite loop to read the file
           l = l + 1 ! Number of the line that is treated
           read(file_unit, '(A)', iostat=stat) buffer
           if (stat == -1) exit ! If the end of the file is reached, stat = -1.
           if (buffer(1:7) == 'Element') then
              read(buffer,*) buffer_string1,buffer_string2,element
              if (element == atom_number) then
                find_element = 0
                exit
              endif
           endif
         enddo
         !write(*,*) 'Element detected at line :', l ! TEST : to see where the element is located in the file

      end function find_element

      integer function read_block ( file_unit, current_block )
      !==========================================================================
      !  Function to read a block of data associated to a specific atom and level
      !==========================================================================
         type(basis_func_block), pointer :: current_block
         integer                         :: file_unit
         integer                         :: line, buffer_number
         character(40)                   :: buffer
         character(7)                    :: buffer_string1,buffer_string2
         integer                         :: stat, designation, orb_mom, max_l, l, ll, n_primitives

         ! First line contains the designation number
         read(file_unit, '(A)', iostat=stat) buffer
         if (buffer(2:6) /= 'Level' .or. stat == -1) then
            read_block = -1
            return
         else
            read(buffer,*) buffer_string1, designation, buffer_string2, max_l
            read_block = 0
         end if
        
         ! Loop for the blocks of a given level
         
         do l = 0, max_l

           read(file_unit, *, iostat=stat) buffer_string1, ll, n_primitives
           if (ll /= l) stop 'correlating functions file is ill-formatted'

             current_block%orb_momentum = l + 1 ! basis file starts at 0 while ROSE starts with 1.
             current_block%designation  = designation
             current_block%n_primitives = n_primitives
             allocate( current_block%exponent( n_primitives ) )

             ! We read the values of the associated exponents
             do line = 1, n_primitives
               read(file_unit,*) current_block%exponent(line)
             enddo

             !print*,"found block with designation ",current_block%designation
             allocate(current_block%next)
             current_block => current_block%next
             current_block%next => NULL()
          
         enddo

      end function read_block

end module
