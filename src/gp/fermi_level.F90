module fermi_level

!  Get energy and location of Fermi level, HOMO and LUMO
!  Written by Lucas Visscher, 2026

  use checkpoint

  implicit none
  private

  character*(:),allocatable :: group

  public fermi_energy, mo_position
  public to_upper, set_mosection

  contains

   subroutine set_mosection(section)

      character(len=*), optional, intent(in)  :: section

      ! By default we consider SCF orbitals, but this can be modified via the optional input argument
      if (present(section)) then
          group = section
      else
          group = '/result/wavefunctions/scf/'
      endif

   end subroutine set_mosection

   integer function mo_position (mo)

      implicit none

      character(len=*), intent(in)  :: mo

      real(8)                       :: fermi_level
      integer                       :: homo, lumo, adjust

      call orbital_information (fermi_level, homo, lumo)

      if (to_upper(trim(mo(1:1))) == 'H') then
         mo_position = homo
      elseif (to_upper(trim(mo(1:1))) == 'L') then
         mo_position = lumo
      else
         print*, " ERROR: position of "//mo//" can not be established."
         stop ' Fatal error in module fermi_level'
      end if

      if (len(trim(mo))>1) then ! allow for lumo+2, homo-3, etc.
         read (mo(2:len(mo)),*) adjust
         mo_position = mo_position + adjust
      end if

   end function mo_position

   real(8) function fermi_energy()

      implicit none

      integer  :: homo, lumo

      call orbital_information (fermi_energy, homo, lumo)

   end function fermi_energy

! generic subroutine, not to be exposed to module users to avoid abuse
! to programmers: please keep the interfaces simple, we may want to revise
! this in the future and extract data from other sources than checkpoint

   subroutine orbital_information (fermi_level, homo, lumo)

      implicit none

      real(8), intent(out)                    :: fermi_level
      integer, intent(out)                    :: homo, lumo

      real(8), allocatable                    :: eig(:)
      integer, allocatable                    :: shell_id(:)
      logical                                 :: data_present
      integer                                 :: i, j, nmo, npo, n_mo(2), n_po(2)

      if (.not.allocated(group)) call set_mosection()

      ! Check whether the necessary data section is indeed present and get number of MOs
      call checkpoint_query(group//'mobasis/n_mo',data_present)
      if (data_present) then
         call checkpoint_read (group//'mobasis/n_mo',idata=n_mo)
      else
         print*, " ERROR: no orbitals present on "//group
         stop ' Fatal error in module fermi_level'
      end if

      ! 4C calculations will have negative energy solutions as well, check that
      call checkpoint_query(group//'mobasis/n_po',data_present)
      if (data_present) then
         call checkpoint_read (group//'mobasis/n_po',idata=n_po)
      else
         n_po = 0
      end if

      nmo = n_mo(1)
      npo = n_po(1)
      if (n_mo(2) > 0 .or. n_po(2) > 0) then
         print*, " Not yet possible to use inversion symmetry in module fermi_level"
         stop ' Fatal error in module fermi_level'
      end if

      call checkpoint_query(group//'mobasis/shell_id',data_present)
      if (data_present) then
         allocate(shell_id(nmo))
         allocate(eig(nmo))
         call checkpoint_read (group//'mobasis/shell_id',idata=shell_id)
         call checkpoint_read (group//'mobasis/eigenvalues',rdata=eig)
      else
         print*, " ERROR: no shell ids present on "//group
         stop ' Fatal error in module fermi_level'
      end if

      homo = 1
      lumo = nmo
      do i = 1, nmo
         if (i <= npo) cycle ! positronic orbitals
         if (shell_id(i) == -1) then ! occupied orbitals
            homo = i
            cycle
         elseif (shell_id(i) > 0) then ! open shell orbitals
            do j = i+1, nmo ! find the lumo (should come after the open shell orbitals)
               if (shell_id(j) == 0) then
                  lumo = j
                  exit
               end if
            end do
         elseif (i /= 1) then ! virtual orbitals
            lumo = i          
         else ! special case: only unoccupied orbitals
            print*, " WARNING: No occupied orbitals, setting Fermi level to first virtual orbital"
            homo = 0
            lumo = 1
            fermi_level = eig(i)
            exit
         end if
         fermi_level = 0.5d0*(eig(homo)+eig(lumo))
         exit ! we know all we need to know and should exit
       end do

   end subroutine orbital_information

   pure function to_upper(str) result(out)
         implicit none
         character(len=*), intent(in) :: str
         character(len=len(str)) :: out
         integer :: i, c

          out = str
          do i = 1, len(str)
             c = iachar(str(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) then
               out(i:i) = achar(c - 32)
            end if
          end do
   end function to_upper

end module fermi_level

