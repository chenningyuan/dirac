module character_processing

  implicit none

  public lowercase
  public word_count
  public word_contains

  private

contains

  subroutine lowercase(s)

!   ----------------------------------------------------------------------------
    character(*), intent(inout) :: s
!   ----------------------------------------------------------------------------
    integer                     :: off, i, ia
!   ----------------------------------------------------------------------------

    off = iachar('a') - iachar('A')

    do i = 1, len(s)
      ia = iachar(s(i:i))
      if (ia >= iachar('A') .and. ia <= iachar('Z')) then
        s(i:i) = achar(ia + off)
      end if
    enddo

  end subroutine

  function word_count(s)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: s
    integer                  :: word_count
!   ----------------------------------------------------------------------------
    integer                  :: i
    logical                  :: is_blank
!   ----------------------------------------------------------------------------

    word_count = 0

    if (len(s) <= 0) return

    is_blank = .true.

    do i = 1, len(s)
      if (s(i:i) == ' ') then
        is_blank = .true.
      else if (is_blank) then
        word_count = word_count + 1
        is_blank = .false.
      end if
    end do

  end function

  function word_contains(word, substring)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: word
    character(*), intent(in) :: substring
!   ----------------------------------------------------------------------------
    logical                  :: word_contains
!   ----------------------------------------------------------------------------

    word_contains = .false.
    if (index(word, substring) > 0) then
      word_contains = .true.
    end if

  end function

end module
