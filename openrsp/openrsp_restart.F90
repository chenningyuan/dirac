module openrsp_restart

!radovan: this module is common to Dalton and DIRAC
!         and is used to restart property calculations
!         (perhaps more interesting in DIRAC where calculations are long)

  implicit none

  public open_restart_file
  public close_restart_file
  public find_restart_tensor
  public read_restart_tensor
  public write_restart_tensor

  private

  interface read_restart_tensor
    module procedure read_restart_tensor_2
    module procedure read_restart_tensor_3
  end interface

  interface write_restart_tensor
    module procedure write_restart_tensor_2
    module procedure write_restart_tensor_3
  end interface

  character(*), parameter :: restart_file = 'openrsp.restart'
  logical                 :: restart_file_exists
  logical                 :: restart_file_is_open
  integer,      parameter :: file_unit    = 137
  integer,      parameter :: l_length = 8
  character(4), parameter :: l_format = '(a8)'
  character(9), parameter :: c_format = '(2e24.16)'

  real(8),      parameter :: tiny = 1.0d-6

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

  function find_restart_tensor(l, freq, tensor_size)

!   return .true. if record found where
!   label .and. frequencies .and. tensor size match

!   ----------------------------------------------------------------------------
    logical                         :: find_restart_tensor
!   ----------------------------------------------------------------------------
    character(l_length), intent(in) :: l
    complex(8),          intent(in) :: freq(:)
    integer,             intent(in) :: tensor_size
!   ----------------------------------------------------------------------------
    character(l_length)             :: l_file
    integer                         :: nr_freq_file
    integer                         :: tensor_size_file
    logical                         :: l_match
    logical                         :: freq_match
    logical                         :: tensor_size_match
    integer                         :: i
    real(8)                         :: cr, ci
!   ----------------------------------------------------------------------------

    rewind(file_unit)

    do while (.true.)

      l_match           = .false.
      freq_match        = .false.
      tensor_size_match = .false.

      read(file_unit, l_format, end = 99) l_file
      if (l == l_file) l_match = .true.

      read(file_unit, *, end = 99) nr_freq_file
      if (size(freq) == nr_freq_file) then
        freq_match = .true.
        do i = 1, nr_freq_file
          read(file_unit, c_format) cr, ci
          if (dabs(cr - real(freq(i)))  > tiny) freq_match = .false.
          if (dabs(ci - aimag(freq(i))) > tiny) freq_match = .false.
        end do
      else
        do i = 1, nr_freq_file
          read(file_unit, c_format)
        end do
      end if

      read(file_unit, *, end = 99) tensor_size_file
      if (tensor_size == tensor_size_file) tensor_size_match = .true.

      if (l_match .and. freq_match .and. tensor_size_match) then
        find_restart_tensor = .true.
        return
      else
        do i = 1, tensor_size_file
          read(file_unit, c_format)
        end do
      end if

    end do

 99 find_restart_tensor = .false.

  end function

  subroutine read_restart_tensor_2(tensor)

!   ----------------------------------------------------------------------------
    complex(8) :: tensor(:, :)
!   ----------------------------------------------------------------------------
    real(8)    :: cr, ci
    integer    :: i, j
!   ----------------------------------------------------------------------------

    do j = 1, size(tensor, 2)
      do i = 1, size(tensor, 1)
        read(file_unit, c_format) cr, ci
        tensor(i, j) = cmplx(cr, ci)
      end do
    end do

  end subroutine

  subroutine read_restart_tensor_3(tensor)

!   ----------------------------------------------------------------------------
    complex(8) :: tensor(:, :, :)
!   ----------------------------------------------------------------------------
    real(8)    :: cr, ci
    integer    :: i, j, k
!   ----------------------------------------------------------------------------

    do k = 1, size(tensor, 3)
      do j = 1, size(tensor, 2)
        do i = 1, size(tensor, 1)
          read(file_unit, c_format) cr, ci
          tensor(i, j, k) = cmplx(cr, ci)
        end do
      end do
    end do

  end subroutine

  subroutine write_restart_tensor_2(l, freq, tensor)

!   move to the end of the file and write record

!   ----------------------------------------------------------------------------
    character(l_length), intent(in) :: l
    complex(8),          intent(in) :: freq(:)
    complex(8),          intent(in) :: tensor(:, :)
!   ----------------------------------------------------------------------------
    integer                         :: i, j
!   ----------------------------------------------------------------------------

    rewind(file_unit)
    do while (.true.)
      read(file_unit, *, end = 1)
    end do

 1  write(file_unit, l_format) l

    write(file_unit, *) size(freq)
    do i = 1, size(freq)
      write(file_unit, c_format) real(freq(i)), aimag(freq(i))
    end do

    write(file_unit, *) size(tensor)
    do j = 1, size(tensor, 2)
      do i = 1, size(tensor, 1)
        write(file_unit, c_format) real(tensor(i, j)), aimag(tensor(i, j))
      end do
    end do

  end subroutine

  subroutine write_restart_tensor_3(l, freq, tensor)

!   move to the end of the file and write record

!   ----------------------------------------------------------------------------
    character(l_length), intent(in) :: l
    complex(8),          intent(in) :: freq(:)
    complex(8),          intent(in) :: tensor(:, :, :)
!   ----------------------------------------------------------------------------
    integer                         :: i, j, k
!   ----------------------------------------------------------------------------

    rewind(file_unit)
    do while (.true.)
      read(file_unit, *, end = 1)
    end do

 1  write(file_unit, l_format) l

    write(file_unit, *) size(freq)
    do i = 1, size(freq)
      write(file_unit, c_format) real(freq(i)), aimag(freq(i))
    end do

    write(file_unit, *) size(tensor)
    do k = 1, size(tensor, 3)
      do j = 1, size(tensor, 2)
        do i = 1, size(tensor, 1)
          write(file_unit, c_format) real(tensor(i, j, k)), aimag(tensor(i, j, k))
        end do
      end do
    end do

  end subroutine

end module
