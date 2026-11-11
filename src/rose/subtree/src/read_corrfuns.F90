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
      real(8), pointer :: exponent(:) => null() ! exponents
      type(basis_func_block), pointer :: next => null()
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
        integer :: i,k,ierr, ishell, nshells, nao, file_unit, nsel
        integer, allocatable :: sel_l(:)
        real(8), allocatable :: sel_exp(:)
        logical :: duplicate
        character(len=64) :: level_string

        ! Should include path information, for now assume this is copied into the working dir
        call get_free_fileunit(file_unit)
        open (file_unit,file=basisfile,form='formatted',status='old',iostat=ierr)
        if (ierr /= 0) stop 'correlating functions: cannot open basis file'

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

        ! Select the shells of the requested levels. The same (l, exponent) pair is listed in neighbouring
        ! levels on purpose (shared functions); keep it only once to avoid linearly dependent shells.
        nsel = 0
        do k = 1, size(levels)
           if (.not. level_present(first_block,levels(k))) then
              write(*,*) 'ERROR: level ',levels(k),' requested for element ',atomic_number, &
                       & ' is not present in file ',trim(basisfile)
              stop 'correlating functions: requested level not found'
           end if
        end do
        current_block => first_block
        do while (associated(current_block))
           if (current_block%n_primitives > 0 .and. ANY (levels == current_block%designation)) &
              & nsel = nsel + current_block%n_primitives
           current_block => current_block%next
        end do
        allocate (sel_l(nsel), sel_exp(nsel))
        nshells = 0
        nao     = 0
        current_block => first_block
        do while (associated(current_block))
           if (current_block%n_primitives > 0 .and. ANY (levels == current_block%designation)) then
              do i = 1, current_block%n_primitives
                 duplicate = .false.
                 do k = 1, nshells
                    if (sel_l(k) == current_block%orb_momentum .and. &
                      & abs(sel_exp(k)-current_block%exponent(i)) <= 1.d-8*abs(sel_exp(k))) then
                       duplicate = .true.
                       exit
                    end if
                 end do
                 if (duplicate) cycle
                 nshells = nshells + 1
                 sel_l(nshells)   = current_block%orb_momentum
                 sel_exp(nshells) = current_block%exponent(i)
                 nao = nao + nfunctions(current_block%orb_momentum,1)
              end do
           end if
           current_block => current_block%next
        end do

        allocate (ao_basis%shell_indices(nao))

        ! allocate basis set and fill in the exponents
        ao_basis%nshells       = nshells
        ao_basis%nao           = nao
        ao_basis%basis_angular = basis_angular
        allocate(ao_basis%gtos(nshells))
        nao = 0
        do ishell = 1, nshells
           ao_basis%gtos(ishell)%orb_momentum = sel_l(ishell)
           ao_basis%gtos(ishell)%atom_element = atomic_number
           ao_basis%gtos(ishell)%coord = coord
           ao_basis%gtos(ishell)%atom_number  = 0 ! this is set to 0 for now but can be generalized.
           ao_basis%gtos(ishell)%n_primitives = 1      ! basis functions are uncontracted
           allocate (ao_basis%gtos(ishell)%exponent(1))
           allocate (ao_basis%gtos(ishell)%coefficient(1))
           ao_basis%gtos(ishell)%exponent(1) = sel_exp(ishell)
           ao_basis%gtos(ishell)%coefficient(1) = &
           & normalization_constant_dirac(sel_l(ishell)-1,sel_exp(ishell))
           ao_basis%shell_indices(nao+1:nao+nfunctions(sel_l(ishell),1)) = ishell
           nao = nao + nfunctions(sel_l(ishell),1)
        end do
        deallocate (sel_l, sel_exp)

        close(file_unit)

        write(level_string,'(*(I0,1X))') levels
        write(*,'(A,I0,A,A,A,I0,A,I0,A)') ' ROSE correlating functions for element ',atomic_number, &
           ', levels ',trim(level_string),': ',nshells,' shells (',nao,' Cartesian functions)'

        ! Everything is read, end by deallocating the memory used for the blocks
        call free_blocks(first_block)

   end subroutine read_corrfuns

      integer function find_element ( file_unit, atom_number )
      !=================================================================
      !  Function to check if the atomic element is present in the file
      !=================================================================
         integer         :: file_unit
         integer         :: atom_number
         character(256)  :: buffer
         character(16)   :: buffer_string1,buffer_string2
         integer         :: l, element, stat

         find_element = -1 ! Case when we do not find the element
         l = 0
         do ! Infinite loop to read the file
           l = l + 1 ! Number of the line that is treated
           read(file_unit, '(A)', iostat=stat) buffer
           if (stat /= 0) exit ! end of file or read error
           if (buffer(1:7) == 'Element') then
              read(buffer,*,iostat=stat) buffer_string1,buffer_string2,element
              if (stat == 0 .and. element == atom_number) then
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
      !  Header: ' Level k Lmax n ...' (n = highest angular momentum) or
      !          ' Level k Symmetries n ...' (n = number of angular momenta, i.e. Lmax+1)
      !==========================================================================
         type(basis_func_block), pointer :: current_block
         integer                         :: file_unit
         integer                         :: line
         character(256)                  :: buffer
         character(16)                   :: buffer_string1,buffer_string2,shape_string
         integer                         :: stat, designation, max_l, l, ll, n_primitives, n

         ! First line contains the designation number
         read(file_unit, '(A)', iostat=stat) buffer
         if (stat /= 0) then
            read_block = -1
            return
         end if
         if (buffer(2:6) /= 'Level') then
            read_block = -1
            return
         end if
         read(buffer,*,iostat=stat) buffer_string1, designation, buffer_string2, n
         if (stat /= 0) stop 'correlating functions file: cannot parse Level header'
         shape_string = buffer_string2
         call to_upper(shape_string)
         select case (trim(shape_string))
         case ('LMAX')
            max_l = n
         case ('SYMMETRIES')
            max_l = n - 1
         case default
            write(*,*) 'Unknown keyword in level header: ',trim(buffer)
            stop 'correlating functions file: expected Lmax or Symmetries'
         end select
         read_block = 0

         ! Loop for the blocks of a given level
         do l = 0, max_l

           read(file_unit, *, iostat=stat) buffer_string1, ll, n_primitives
           if (stat /= 0) stop 'correlating functions file is ill-formatted (angular momentum line)'
           if (ll /= l) stop 'correlating functions file is ill-formatted'

             current_block%orb_momentum = l + 1 ! basis file starts at 0 while ROSE starts with 1.
             current_block%designation  = designation
             current_block%n_primitives = n_primitives
             allocate( current_block%exponent( n_primitives ) )

             ! We read the values of the associated exponents
             do line = 1, n_primitives
               read(file_unit,*,iostat=stat) current_block%exponent(line)
               if (stat /= 0) stop 'correlating functions file is ill-formatted (exponent)'
             enddo

             allocate(current_block%next)
             current_block => current_block%next
             current_block%next => NULL()

         enddo

      end function read_block

      subroutine to_upper (string)
         character(*), intent(inout) :: string
         integer :: i
         do i = 1, len(string)
            if (string(i:i) >= 'a' .and. string(i:i) <= 'z') string(i:i) = achar(iachar(string(i:i))-32)
         end do
      end subroutine to_upper

      logical function level_present (first_block, level)
         type(basis_func_block), pointer :: first_block, b
         integer, intent(in) :: level
         level_present = .false.
         b => first_block
         do while (associated(b))
            if (b%n_primitives > 0 .and. b%designation == level) level_present = .true.
            b => b%next
         end do
      end function level_present

      subroutine free_blocks (first_block)
         type(basis_func_block), pointer :: first_block, b, nxt
         b => first_block
         do while (associated(b))
            nxt => b%next
            if (associated(b%exponent)) deallocate(b%exponent)
            deallocate(b)
            b => nxt
         end do
         first_block => null()
      end subroutine free_blocks

end module
