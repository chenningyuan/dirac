module visual_labeled_storage

  use labeled_storage
  use visual_datatypes, only: rs_label_t, grid_t, grid_function_t
  use visual_utils, only: visual_quit

  implicit none
  private

  public write_grid_to_checkpoint
  public write_grid_function_to_checkpoint
  public get_grid_info_from_labeled_storage
  public get_grid_points_from_labeled_storage
  public get_grid_function_info_from_labeled_storage

contains

  subroutine assign_label_grid_function(grid_function,label,top_section,level)
!   ----------------------------------------------------------------------------
    type(grid_function_t), intent(in):: grid_function
    type(rs_label_t), intent(out):: label
    character(len=*), intent(in) :: top_section
    integer, intent(in) :: level
    character*(:),allocatable :: group_label
    character*(:),allocatable :: l0, l1, l2, l3, l4
    character(len = 200):: errmsg
    character(len=1)::temp
!   ----------------------------------------------------------------------------
    if (top_section /= 'input' .and. top_section /= 'result') return

    l0 = trim(top_section)
    l1 = l0 // '/grid_function/' // trim(grid_function%name)

    select case (grid_function%purpose)
      case (0)
        l2 = l1 // "/visualize"
      case (1)
        l2 = l1 // "/integrate"
    end select

    l3 = l2 // '/values'

    write(temp, '(i1)') grid_function%tensor_order

    l4 = l3 // '/' // '/t' // trim(temp)

    select case (level)
      case (1)
        group_label = trim(l1)
      case (2)
        group_label = trim(l2)
      case (3)
        group_label = trim(l3)
      case (4)
        group_label = trim(l4)
    end select

    if (len(group_label) > 80) then
      ! todo - check whether this limitation is due to DIRAC or HDF5
      write (errmsg, '(a)') 'STOP. Grid function label is too long. ', group_label
      call visual_quit(errmsg)
    else
      label%name = trim(group_label)
    end if

  end subroutine assign_label_grid_function


  subroutine assign_label_grid(grid, label, top_section, level)
!   ----------------------------------------------------------------------------
    type(grid_t), intent(in):: grid
    type(rs_label_t), intent(out):: label
    character(len=*), intent(in) :: top_section
    integer, intent(in) :: level
    character*(:),allocatable :: group_label
    character*(:),allocatable :: l0, l1, l2
    character(len = 200):: errmsg
!   ----------------------------------------------------------------------------
    if (top_section /= 'input' .and. top_section /= 'result') return

    l0 = trim(top_section)
    l1 = l0//'/grid'
    l2 = l0//'/grid/points'

    select case (level)
      case (1)
        group_label = trim(l1)
      case (2)
        group_label = trim(l2)
    end select

    if (len(group_label) > 80) then
      !  todo - check whether this limitation is due to DIRAC or HDF5
      write (errmsg, '(a)') 'STOP. Grid function label is too long. ', group_label
      call visual_quit(errmsg)
    else
      label%name = trim(group_label)
    end if
  end subroutine assign_label_grid


  subroutine write_grid_to_checkpoint(data_file, grid)
!   ----------------------------------------------------------------------------
    type(file_info_t), intent(inout):: data_file
    type(grid_t), intent(in):: grid
    type(rs_label_t):: label
    character*(:), allocatable  :: group
    character(len=1):: ich
    integer:: i
!   ----------------------------------------------------------------------------
    call assign_label_grid(grid, label, top_section='result', level=1)
    group = label%name
    call lab_write(data_file, trim(group)//'/nr_points_total',idata = grid%nr_points_total)
    call lab_write(data_file, trim(group)//'/ndim',idata = grid%ndim)

    call assign_label_grid(grid, label, top_section='result', level=2)
    group = label%name
    do i = 1, grid%ndim
      write(ich, '(i1)') i
      call lab_write(data_file, trim(group)//'/coor_'//trim(ich), rdata = grid%points(i, :))
    end do
  end subroutine


  subroutine write_grid_function_to_checkpoint(data_file, grid_function)
!   ----------------------------------------------------------------------------
    type(file_info_t), intent(inout):: data_file
    type(grid_function_t), intent(in):: grid_function
    type(rs_label_t):: label
    character*(:), allocatable  :: group
    character(len=1):: ich
    integer:: i
!   ----------------------------------------------------------------------------
    call assign_label_grid_function(grid_function, label, top_section='result', level=2)
    group = label%name
    call lab_write(data_file, group)

    call assign_label_grid_function(grid_function, label, top_section='result', level=4)
    group = label%name
    do i = 1, grid_function%get_ndim()
      write(ich, '(i1)') i
      call lab_write(data_file, group//'_icomp_'//ich, rdata = grid_function%values(i, :))
    end do
  end subroutine


  subroutine get_grid_info_from_labeled_storage(data_file, grid)
!   ----------------------------------------------------------------------------
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check whether "grid" group is under "input" or "result"
    type(file_info_t), intent(inout):: data_file
    type(grid_t), intent(inout):: grid
    type(rs_label_t)::label
    integer:: i
    character*(:), allocatable  :: group
!   ----------------------------------------------------------------------------
    call assign_label_grid(grid, label, top_section='result', level=1)
    group = label%name
    call lab_read(data_file, trim(group)//'/ndim',idata = grid%ndim)
    call lab_read(data_file, trim(group)//'/nr_points_total',idata = grid%nr_points_total)
    grid%ndim_known = .true.
    grid%nr_points_total_known = .true.
  end subroutine


  subroutine get_grid_points_from_labeled_storage(data_file, grid)
!   ----------------------------------------------------------------------------
    ! gosia fixme - unused thus far
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check whether "grid" group is under "input" or "result"
    type(file_info_t), intent(inout):: data_file
    type(grid_t), intent(inout):: grid
    type(rs_label_t)::label
    integer:: i
    character*(:), allocatable  :: group
    character(len=1):: ich
!   ----------------------------------------------------------------------------
    call assign_label_grid(grid, label, top_section='result', level=2)
    group = label%name
    do i = 1, grid%ndim
      write(ich, '(i1)') i
      write(*,*) 'will read ', trim(group)//'/coor_'//trim(ich), ' from ', data_file%name
      call lab_read(data_file, trim(group)//'/coor_'//trim(ich), rdata = grid%points(i,:))
    end do
  end subroutine


  subroutine get_grid_function_info_from_labeled_storage(data_file,grid_function)
!   ----------------------------------------------------------------------------
    ! gosia fixme - unused thus far
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check whether "grid function" group is under "input" or "result"
    type(grid_function_t), intent(inout):: grid_function
    character*(:), allocatable  :: group
    type(file_info_t):: data_file
    type(rs_label_t)::label
    integer:: i
!   ----------------------------------------------------------------------------
    call assign_label_grid_function(grid_function, label, top_section='result', level=1)
  end subroutine

end module
