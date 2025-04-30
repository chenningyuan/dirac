module visual_labeled_storage

  use labeled_storage
  use visual_utils
  use visual_datatypes, only: rs_label_t, grid_t, grid_function_t

  implicit none

  public write_grid_to_checkpoint
  public write_grid_function_to_checkpoint
  public get_grid_info_from_labeled_storage
  public get_grid_points_from_labeled_storage
  public get_grid_weights_from_labeled_storage
  public get_grid_function_info_from_labeled_storage

  private

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
        l2 = l1 // "/visualization"
        l3 = l2 // '/values'
      case (1)
        l2 = l1 // "/integration"
        l3 = l2 // '/integrated_values'
    end select

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
    if (grid%typ == 2) then
      l2 = l0//'/grid/radial_points'
    else
      l2 = l0//'/grid/points'
    end if

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


! ------------------------------------------------------------------------------
  subroutine write_grid_to_checkpoint(data_file, grid)
    type(file_info_t), intent(inout):: data_file
    type(grid_t), intent(in):: grid
    type(rs_label_t):: label
    character*(:), allocatable:: group
    character(len=1):: ich
    integer:: i

    call assign_label_grid(grid, label, top_section='result', level=1)
    group = label%name
    if (grid%typ == 2) then
      call lab_write(data_file, trim(group)//'/radial_nr_steps',&
                     idata = grid%radial_nr_steps)
      call lab_write(data_file, trim(group)//'/radial_nr_points_angular',&
                     idata = grid%radial_nr_points_angular)
    else
      call lab_write(data_file, trim(group)//'/nr_points_total',&
                     idata = grid%nr_points_total)
      call lab_write(data_file, trim(group)//'/ndim',idata = grid%ndim)
      do i = 1, grid%ndim
        write(ich, '(i1)') i
        if (grid%nr_points_known(i)) &
         call lab_write(data_file, trim(group)//'/nr_points_dim_'//trim(ich), idata = grid%nr_points(i))
        if (grid%spacing_known(i)) &
         call lab_write(data_file, trim(group)//'/spacing_dim_'//trim(ich), rdata = grid%spacing(i))
        if (grid%origin_known) &
         call lab_write(data_file, trim(group)//'/origin_dim_'//trim(ich), rdata = grid%origin(i))
      end do
    end if

    call assign_label_grid(grid, label, top_section='result', level=2)
    group = label%name
    if (grid%typ == 2) then
      if (grid%allocated_radial_coord()) then
        !call lab_write(data_file, trim(group)//'/radial_coord', rdata = grid%radial_coord(:))
        call lab_write(data_file, trim(group), rdata = grid%radial_coord(:))
      end if
    else
      do i = 1, grid%ndim
        write(ich, '(i1)') i
        call lab_write(data_file, trim(group)//'/coor_'//trim(ich), rdata = grid%points(i, :))
      end do
      if (grid%allocated_weights()) then
        call lab_write(data_file, trim(group)//'/weights', rdata = grid%weights(:))
      end if
    end if
  end subroutine


! ------------------------------------------------------------------------------
  subroutine write_grid_function_to_checkpoint(data_file, grid_function)
    type(file_info_t), intent(inout):: data_file
    type(grid_function_t), intent(in):: grid_function
    type(rs_label_t):: label
    character*(:), allocatable  :: group
    character(len=1):: ich
    integer:: i

    call assign_label_grid_function(grid_function, label, top_section='result', level=2)
    group = label%name
    call lab_write(data_file, group)

    call assign_label_grid_function(grid_function, label, top_section='result', level=4)
    group = label%name
    do i = 1, grid_function%get_ndim()
      write(ich, '(i1)') i
      if (grid_function%purpose == 0 .and. grid_function%allocated_values()) then
        call lab_write(data_file, group//'_icomp_'//ich, &
             rdata = grid_function%values(i, :))
      end if
      if (grid_function%purpose == 1 .and. grid_function%allocated_integrated_values()) then
        call lab_write(data_file, group//'_icomp_'//ich, &
             rdata = grid_function%integrated_values(i))
      end if
    end do
  end subroutine

! ------------------------------------------------------------------------------
  subroutine get_grid_info_from_labeled_storage(data_file, grid)
    ! gosia todo:
    ! - add support for other formats (now, hdf5 format only: data_file%type = 2)
    ! - check whether "grid" group is under "input" or "result"
    type(file_info_t), intent(inout):: data_file
    type(grid_t), intent(inout):: grid
    type(rs_label_t)::label
    integer:: i
    character*(:), allocatable  :: group
    character(len=1):: ich
    logical::data_available

    call assign_label_grid(grid, label, top_section='result', level=1)
    group = label%name
    call lab_read(data_file, trim(group)//'/ndim',idata = grid%ndim)
    call lab_read(data_file, trim(group)//'/nr_points_total',idata = grid%nr_points_total)
    do i = 1, grid%ndim
      write(ich, '(i1)') i
      call lab_query(data_file, trim(group)//'/nr_points_dim_'//trim(ich), &
                     exist=data_available)
      if (data_available) then
        call lab_read(data_file, trim(group)//'/nr_points_dim_'//trim(ich), idata = grid%nr_points(i))
        grid%nr_points_known(i) = .true.
      end if
    end do
    if (grid%nr_points_total > 0) grid%nr_points_total_known = .true.
  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_grid_points_from_labeled_storage(data_file, grid)
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

! ------------------------------------------------------------------------------
  subroutine get_grid_weights_from_labeled_storage(data_file, grid)
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
    write(*,*) 'will read ', trim(group)//'/coor_'//trim(ich), ' from ', data_file%name
    call lab_read(data_file, trim(group)//'/weights', rdata = grid%weights(:))
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
