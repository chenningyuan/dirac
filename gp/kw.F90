module kw

  use character_processing

  implicit none

  public kw_matches
  public reset_available_kw_list
  public check_whether_kw_found

  private

  integer, parameter   :: kw_length = 7
  integer, parameter   :: max_nr_kw = 200

  character(kw_length) :: available_kw_list(max_nr_kw)
  integer              :: nr_available_kw
  logical              :: kw_found

contains

  subroutine reset_available_kw_list()

    nr_available_kw = 0
    kw_found        = .false.

  end subroutine

  subroutine check_whether_kw_found(kw_input, kw_section)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: kw_input
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i
!   ----------------------------------------------------------------------------

    if (.not. kw_found) then

      write(*, *) 'illegal keyword '//kw_input//' in section '//kw_section

      write(*, *) 'list of available keywords in section '//kw_section//':'
      do i = 1, nr_available_kw
        write(*, *) available_kw_list(i)
      end do

      call quit('illegal keyword '//kw_input//' in section '//kw_section)

    end if

  end subroutine

  function kw_matches(kw_input, kw_option)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: kw_input
    character(kw_length), intent(in) :: kw_option
!   ----------------------------------------------------------------------------
    character(kw_length)             :: kw_input_lower
    character(kw_length)             :: kw_option_lower
    logical                          :: kw_matches
!   ----------------------------------------------------------------------------

    kw_input_lower  = kw_input
    kw_option_lower = kw_option

    call lowercase(kw_input_lower)
    call lowercase(kw_option_lower)

    if (kw_input_lower == kw_option_lower) then
      kw_matches = .true.
      kw_found   = .true.
    else
      kw_matches = .false.
      nr_available_kw = nr_available_kw + 1
      available_kw_list(nr_available_kw) = kw_option
    end if

  end function

end module
