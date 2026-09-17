module rose_basis_locations

  implicit none
  private

  public :: set_basis_file_location
  public :: get_basis_file_location
  public :: has_basis_file_location
  public :: clear_basis_file_locations

  integer, parameter :: basis_name_len = 80
  integer, parameter :: basis_role_len = 32
  integer, parameter :: basis_file_len = 300

  type :: basis_location_t
     character(len=basis_name_len) :: name = ''
     character(len=basis_role_len) :: role = ''
     character(len=basis_file_len) :: file = ''
  end type basis_location_t

  type(basis_location_t), allocatable :: basis_locations(:)

contains

  subroutine set_basis_file_location(basis_name,basis_file,basis_role)

    character(len=*), intent(in) :: basis_name
    character(len=*), intent(in) :: basis_file
    character(len=*), intent(in), optional :: basis_role

    character(len=basis_role_len) :: role
    integer :: i

    if (present(basis_role)) then
       role = trim(basis_role)
    else
       role = 'orbital'
    end if

    if (.not. allocated(basis_locations)) allocate(basis_locations(0))

    do i = 1, size(basis_locations)
       if (trim(basis_locations(i)%name) == trim(basis_name) .and. &
           trim(basis_locations(i)%role) == trim(role)) then
          basis_locations(i)%file = trim(basis_file)
          return
       end if
    end do

    basis_locations = [basis_locations, basis_location_t(trim(basis_name),role,trim(basis_file))]

  end subroutine set_basis_file_location

  recursive subroutine get_basis_file_location(basis_name,basis_file,found,basis_role,fallback_role)

    character(len=*), intent(in) :: basis_name
    character(len=*), intent(out) :: basis_file
    logical, intent(out), optional :: found
    character(len=*), intent(in), optional :: basis_role
    character(len=*), intent(in), optional :: fallback_role

    character(len=basis_role_len) :: role
    integer :: i

    if (present(basis_role)) then
       role = trim(basis_role)
    else
       role = 'orbital'
    end if
    basis_file = ''

    if (allocated(basis_locations)) then
       do i = 1, size(basis_locations)
          if (trim(basis_locations(i)%name) == trim(basis_name) .and. &
              trim(basis_locations(i)%role) == trim(role)) then
             basis_file = trim(basis_locations(i)%file)
             if (present(found)) found = .true.
             return
          end if
       end do
    end if

    if (present(fallback_role)) then
       if (trim(fallback_role) /= trim(role)) then
          call get_basis_file_location(basis_name,basis_file,found,basis_role=fallback_role)
          return
       end if
    end if

    if (present(found)) found = .false.

  end subroutine get_basis_file_location

  logical function has_basis_file_location(basis_name,basis_role)

    character(len=*), intent(in) :: basis_name
    character(len=*), intent(in), optional :: basis_role

    character(len=basis_file_len) :: basis_file

    call get_basis_file_location(basis_name,basis_file,has_basis_file_location,basis_role)

  end function has_basis_file_location

  subroutine clear_basis_file_locations()

    if (allocated(basis_locations)) deallocate(basis_locations)

  end subroutine clear_basis_file_locations

end module rose_basis_locations
