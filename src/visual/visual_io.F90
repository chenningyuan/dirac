module visual_io

  use visual_txt_storage
  use visual_labeled_storage

  use labeled_storage
  use visual_cfg
  use visual_utils
  use visual_datatypes

  implicit none

  public store_grid
  public store_data
  public get_grid_info_from_file
  public get_grid_points_from_file
  public get_grid_function_info_from_file

  private

contains

! -----------------------------------------------------------------------------
  subroutine store_grid(grid)
    type(grid_t), intent(in):: grid
    type(file_info_t):: data_file
    ! todo:
    ! - add support for other formats, now - hdf5 format only
    ! - check what file%status == -1 means ('never used' = 'file not exists'?)
    !   now - hardcoded (after labeled_storage)

    data_file%name = grid%file_out
    data_file%status = -1
    if (data_file%status == -1) then
      write(*,*) " Opening new HDF5 file for grid export: ", data_file%name
      data_file%type = 2
      data_file%status = 0
    end if

    select case (data_file%type)
      case (2)
        call write_grid_to_checkpoint(data_file, grid)
      case default
        call visual_quit('programming error in store_grid; unsupported file type')
    end select

  end subroutine


! -----------------------------------------------------------------------------
  subroutine store_data(grid, grid_function)
    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(in):: grid_function
    type(file_info_t):: data_file
    integer::i
    character(len=:), allocatable:: fmt
    ! todo:
    ! - add support for other formats, now - hdf5 format only
    ! - check what file%status == -1 means ('never used' = 'file not exists'?)
    !   now - hardcoded (after labeled_storage)

    data_file%name = grid_function%file_out
    data_file%status = -1
    if (data_file%status == -1) then
      write(*,*) " Opening new HDF5 file for grid function export: ", data_file%name
      data_file%type = 2
      data_file%status = 0
    end if

    select case (data_file%type)
      case (2)
        call write_grid_to_checkpoint(data_file, grid)
        call write_grid_function_to_checkpoint(data_file, grid_function)
      case default
        call visual_quit('programming error in store_data; unsupported file type')
    end select

!   test
    !fmt = '(*(es18.8e3, :))'
    !do i = 1, grid%nr_points_total
    !  write(*, fmt) grid%points(1,i), grid%points(2,i), grid%points(3,i)
    !end do

  end subroutine


! -----------------------------------------------------------------------------
  subroutine get_grid_info_from_file(grid)
    type(grid_t), intent(inout):: grid
    type(file_info_t):: data_file
    ! todo:
    ! - add support for other formats, now - hdf5 format only
    ! - check what file%status == -1 means ('never used' = 'file not exists'?)
    !   now - hardcoded (after labeled_storage)

    if (grid%input  /=  1) return

    !if (data_file%status == -1) return ! file not exists
    data_file%status = 0 ! assuming file present and closed
    data_file%name = grid%file_inp
    data_file%type = 2

    select case (data_file%type)
      case (2)
        call get_grid_info_from_labeled_storage(data_file, grid)
      case default
        call visual_quit('programming error in get_grid_info_from_file; unsupported file type')
    end select

  end subroutine get_grid_info_from_file


! -----------------------------------------------------------------------------
  subroutine get_grid_points_from_file(grid)
    type(grid_t), intent(inout):: grid
    type(file_info_t):: data_file

    if (grid%input  /=  1) return

    !if (data_file%status == -1) return ! file not exists
    data_file%status = 0
    data_file%name = grid%file_inp
    data_file%type = 2

    select case (data_file%type)
      case (2)
        call get_grid_points_from_labeled_storage(data_file, grid)
      case default
        call visual_quit('programming error in get_grid_points_from_file; unsupported file type')
    end select
  end subroutine get_grid_points_from_file


! -----------------------------------------------------------------------------
  subroutine get_grid_function_info_from_file(grid_function)
    type(grid_function_t), intent(inout):: grid_function
    character*(:), allocatable  :: group
    type(file_info_t):: data_file
    type(rs_label_t)::label
    integer:: i

    if (grid_function%input  /=  1) return
    if (grid_function%status /= -1) return

    data_file%name = grid_function%file_inp
    data_file%type = 2
    data_file%status = 0

    select case (data_file%type)
      case (2)
        call get_grid_function_info_from_labeled_storage(data_file, grid_function)
      case default
        call visual_quit('programming error in get_grid_function_info_from_file; unsupported file type')
    end select
  end subroutine get_grid_function_info_from_file

end module
