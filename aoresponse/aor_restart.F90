module aor_restart

  use dirac_interface

  implicit none

  public open_restart_file
  public close_restart_file
  public find_restart_tensor
  public read_restart_tensor
  public write_restart_tensor

  private

  character(*), parameter :: restart_file = 'aoresponse.restart'
  logical                 :: restart_file_exists
  logical                 :: restart_file_is_open
  integer,      parameter :: file_unit    = 137
  integer,      parameter :: label_length = 16
  character(5), parameter :: label_format = '(a16)'

contains

  subroutine open_restart_file()

    inquire(file  = restart_file, &
            exist = restart_file_exists)
    if (restart_file_exists) then
      inquire(file   = restart_file, &
              opened = restart_file_is_open)
      if (.not. restart_file_is_open) then
        open(file_unit,             &
             file   = restart_file, &
             status = 'unknown',    &
             form   = 'formatted',  &
             access = 'sequential')
      end if
    else
      open(file_unit,             &
           file   = restart_file, &
           status = 'new',        &
           form   = 'formatted',  &
           access = 'sequential')
      rewind(file_unit)
    end if
    
  end subroutine

  subroutine close_restart_file()

    close(file_unit, status = 'keep')

  end subroutine

  function find_restart_tensor(l, tensor_size)

!   ----------------------------------------------------------------------------
    character(label_length), intent(in) :: l
    integer,                 intent(in) :: tensor_size
    logical                             :: find_restart_tensor
!   ----------------------------------------------------------------------------
    character(label_length)             :: l_file
    integer                             :: tensor_size_file
    integer                             :: i
!   ----------------------------------------------------------------------------

    rewind(file_unit)
   
    do while (.true.) 

      read(file_unit, label_format, end = 99) l_file
      read(file_unit, *) tensor_size_file

      if ((l == l_file) .and. (tensor_size == tensor_size_file)) then
        find_restart_tensor = .true.
        return
      else
        do i = 1, tensor_size_file
          read(file_unit, *)
        end do
      end if

    end do

 99 find_restart_tensor = .false.

  end function

  subroutine read_restart_tensor()

! rewind
! find label
! read label

  end subroutine

  subroutine write_restart_tensor()

! goto end
! write

  end subroutine

end module
