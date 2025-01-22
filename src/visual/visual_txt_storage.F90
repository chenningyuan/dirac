module visual_txt_storage

  use visual_utils
  use visual_datatypes, only: rs_label_t, grid_t, grid_function_t

  implicit none

  public write_plot_data_to_DIRAC_output
  public get_grid_info_from_txt_storage
  public get_grid_points_from_txt_storage

  type, public:: txt_file_info_t
    integer:: unit = -1
    character(len=:), allocatable:: name
    character(len=:), allocatable:: format
    character(len=:), allocatable:: domain
    character(len=:), allocatable:: data_type
    character(len=:), allocatable:: float_format
    logical:: exists = .false.
    logical:: is_opened = .false.
  end type

  character(len = 16):: float_fmt = '(*(es18.8e3, :))'
  private

contains

! -----------------------------------------------------------------------------
  subroutine write_plot_data_to_DIRAC_output(grid, grid_function)
    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(in):: grid_function
    integer:: g_ndim, gf_ndim, npoints, npoints_radial
    integer:: i, j

    npoints = grid%nr_points_total
    npoints_radial = grid%radial_nr_steps
    g_ndim  = grid%ndim
    gf_ndim = grid_function%get_ndim()

    write(*, '(a)') '**************************************'
    if (grid_function%purpose == 1) then 
      write(*, '(a)') 'Final grid_function values integrated:'
      write(*, float_fmt) (grid_function%integrated_values(i),i=1,gf_ndim)
    else if (grid_function%purpose == 0) then
      write(*, '(a)') 'Final grid_function values in grid points:'
!add radial!
      if (grid%allocated_radial_coord()) then
        write(*, '(a)') 'radial_step,grid_function values'
      else if (grid%allocated_weights()) then
        write(*, '(a)') 'x,y,z,w,grid_function values'
      else
        write(*, '(a)') 'x,y,z,grid_function values'
      end if
      if (grid%allocated_radial_coord()) then
        do j = 1, npoints_radial
          write(*, float_fmt) (grid%radial_coord(j)), &
                              (grid_function%values(i,j),i=1,gf_ndim)
        end do
      else
        if (grid%allocated_weights()) then
          do j = 1, npoints
            write(*, float_fmt) (grid%points(i,j),i=1,g_ndim), &
                                 grid%weights(j), &
                                (grid_function%values(i,j),i=1,gf_ndim)
          end do
        else
          do j = 1, npoints
            write(*, float_fmt) (grid%points(i,j),i=1,g_ndim), &
                                (grid_function%values(i,j),i=1,gf_ndim)
          end do
        end if
      end if
    end if
    write(*, '(a)') '**************************************'

   end subroutine

! ------------------------------------------------------------------------------
  subroutine open_txt_file(plot_file)
    type(txt_file_info_t), intent(inout):: plot_file
    character(len = 100):: errmsg
    logical:: exists, is_opened

    inquire(plot_file%unit, exist = exists, opened = is_opened)
    if (exists) then
      plot_file%exists    = .true.
      if (.not.is_opened) then
        open(unit   = plot_file%unit, &
             file   = plot_file%name, &
             status = 'unknown',      &
             access = 'sequential',   &
             form   = 'formatted')
        plot_file%is_opened = .true.
      end if
    else
      write(errmsg,*) 'file ',plot_file%name,' does not exist. Stop.'
      call visual_quit(errmsg)
    end if
 
  end subroutine

! ------------------------------------------------------------------------------
  subroutine close_txt_file(plot_file)
    type(txt_file_info_t),intent(inout):: plot_file
    if (plot_file%is_opened) close(plot_file%unit, status='keep')
  end subroutine

! ------------------------------------------------------------------------------
  subroutine get_grid_info_from_txt_storage(grid)
use interface_grid
    type(txt_file_info_t):: grid_file
    type(grid_t), intent(inout):: grid
    integer:: ival, ierror, nr_points

    if (grid%file_inp /= 'numerical_grid') return

    grid_file%name = 'numerical_grid'
    grid_file%unit = newunit()

    call open_txt_file(grid_file)
    rewind(grid_file%unit)

    nr_points = 0
    do
      read(grid_file%unit, '(i20)',iostat=ierror) ival
      if (ierror == 0) then
        if (ival > 0) then
          nr_points = nr_points + ival
        else
          exit
        end if
      else
        cycle
      end if
    end do
    grid%nr_points_total = nr_points
    if (grid%nr_points_total > 0) grid%nr_points_total_known = .true.

    call close_txt_file(grid_file)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_grid_points_from_txt_storage(grid)
use interface_grid
    type(txt_file_info_t):: grid_file
    type(grid_t), intent(inout):: grid
    integer:: ioff, ierror, nr_points

    if (grid%file_inp /= 'numerical_grid') return

    grid_file%name = 'numerical_grid'
    grid_file%unit = newunit()

    call open_txt_file(grid_file)
    rewind(grid_file%unit)

    ioff = 1
 1  CONTINUE
    READ(grid_file%unit, *) NR_POINTS

    IF(NR_POINTS .LE. 0) GOTO 2
    call num_grid_read(grid%points(1,ioff:ioff+nr_points), &
                       grid%points(2,ioff:ioff+nr_points), &
                       grid%points(3,ioff:ioff+nr_points), &
                       grid%weights(ioff:ioff+nr_points), &
                       grid_file%unit, nr_points)
    ioff = nr_points + ioff
    GOTO 1
 2  CONTINUE

    call close_txt_file(grid_file)

  end subroutine

end module
