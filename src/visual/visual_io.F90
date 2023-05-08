module visual_io

  use visual_txt_storage
  use visual_labeled_storage

  use labeled_storage, only: file_info_t
  use visual_cfg
  use visual_utils
  use visual_datatypes

  implicit none

  public store_grid
  public store_data
  public get_grid_info_from_file
  public get_grid_points_from_file
  public get_grid_function_info_from_file

!  public visual_plot_file

  private

!! to do: extends file_info_t
!  type:: visual_plot_file
!    integer:: unit = -1
!    character(len=:), allocatable:: name
!    character(len=:), allocatable:: format
!    character(len=:), allocatable:: domain
!    character(len=:), allocatable:: data_type
!    character(len=:), allocatable:: float_format
!    logical:: exists = .false.                                                                           
!    logical:: is_opened = .false.                                                                        
!  !contains
!    !    procedure:: assign, get, length
!    !    generic:: assignment(=) => assign
!    !    generic:: operator(-) => get
!    !    generic:: l => length
!  end type visual_plot_file

  !interface store_data
  !    !module procedure store_array_data
  !    module procedure store_labeled_data
  !end interface

contains

  subroutine store_grid(grid)
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check what grid%status == -1 means ('never used' a.k.a 'file not exists'?)
    type(grid_t), intent(in):: grid
    type(file_info_t):: data_file

    data_file%name = grid%file_out
    ! gosia fixme - hardcoded
    data_file%status = -1
    if (data_file%status == -1) then
      write(*,*) " Opening new HDF5 file for grid export", data_file%name
      data_file%type = 2
      data_file%status = 0
    end if

    select case (data_file%type)
      case (2)
        call write_grid_to_checkpoint(data_file, grid)
    end select

  end subroutine


  subroutine store_data(grid, grid_function)
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check what grid%status == -1 means ('never used' a.k.a 'file not exists'?)
    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(in):: grid_function
    type(file_info_t):: data_file

    data_file%name = grid_function%file_out
    ! gosia fixme - hardcoded
    data_file%status = -1
    if (data_file%status == -1) then
      write(*,*) " Opening new HDF5 file for full data export", data_file%name
      data_file%type = 2
      data_file%status = 0
    end if

    select case (data_file%type)
      case (2)
        call write_grid_to_checkpoint(data_file, grid)
        call write_grid_function_to_checkpoint(data_file, grid_function)
    end select
  end subroutine


  subroutine get_grid_info_from_file(grid)
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check what grid%status == -1 means ('never used' a.k.a 'file not exists'?)
    type(grid_t), intent(inout):: grid
    type(file_info_t):: data_file

    if (grid%input  /=  1) return
    if (grid%status /= -1) return

    data_file%name = grid%file_inp
    ! gosia fixme - hardcoded
    data_file%type = 2
    data_file%status = 0

    select case (data_file%type)
      case (2)
        call get_grid_info_from_labeled_storage(data_file, grid)
    end select
  end subroutine get_grid_info_from_file


  subroutine get_grid_points_from_file(grid)
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check what grid%status == -1 means ('never used' a.k.a 'file not exists'?)
    type(grid_t), intent(inout):: grid
    type(file_info_t):: data_file

    if (grid%input  /=  1) return
    if (grid%status /= -1) return

    data_file%name = grid%file_inp
    ! gosia fixme - hardcoded
    data_file%type = 2
    data_file%status = 0

    select case (data_file%type)
      case (2)
        call get_grid_points_from_labeled_storage(data_file, grid)
    end select
  end subroutine get_grid_points_from_file


  subroutine get_grid_function_info_from_file(grid_function)
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check what grid%status == -1 means ('never used' a.k.a 'file not exists'?)
    type(grid_function_t), intent(inout):: grid_function
    character*(:), allocatable  :: group
    type(file_info_t):: data_file
    type(rs_label_t)::label
    integer:: i

    if (grid_function%input  /=  1) return
    if (grid_function%status /= -1) return

    data_file%name = grid_function%file_inp
    ! gosia fixme - hardcoded
    data_file%type = 2
    data_file%status = 0

    select case (data_file%type)
      case (2)
        call get_grid_function_info_from_file_from_labeled_storage(data_file, grid_function)
    end select
  end subroutine get_grid_function_info_from_file

end module
