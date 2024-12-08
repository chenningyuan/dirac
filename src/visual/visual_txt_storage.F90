module visual_txt_storage

  use visual_utils
  use visual_datatypes, only: rs_label_t, grid_t, grid_function_t

  implicit none

  public write_plot_data
  public get_grid_info_from_txt_storage
  !WIPpublic get_grid_points_from_txt_storage
  !WIPpublic get_grid_function_info_from_txt_storage

  type, public:: visual_plot_file_t
    integer:: unit = -1
    character(len=:), allocatable:: name
    character(len=:), allocatable:: format
    character(len=:), allocatable:: domain
    character(len=:), allocatable:: data_type
    character(len=:), allocatable:: float_format
    logical:: exists = .false.                                                                           
    logical:: is_opened = .false.                                                                        
  contains
    procedure, public:: display => display_visual_plot_file
  end type

  interface write_plot_data 
    module procedure write_plot_data_txt_by_line
      !TODOmodule procedure write_plot_data_cube
  end interface

  character(len = 16):: float_format_txt = '(*(es18.8e3, :))'
  private

contains

! -----------------------------------------------------------------------------
  subroutine write_plot_data_txt_by_line(plot_file, grid_point, grid_function,&
                                         outpri)
    type(visual_plot_file_t), intent(inout):: plot_file
    real(8), intent(in):: grid_point(:)
    real(8), intent(in):: grid_function(:)
    logical, intent(in), optional:: outpri
    character(len=:), allocatable:: fmt
    integer:: i, j

    if (present(outpri)) then
      if (outpri) then
        ! also write to DIRAC output
        fmt = '(*(es18.8e3, :))'
        write(*, fmt) (grid_point(j), j = 1, size(grid_point)), &
                      (grid_function(j), j = 1, size(grid_function))
      end if
    end if

         ! write to files
         !do i = 1, visual_cfg_nr_plot_files
         !
         !   fmt = '(*('//trim(plot_files(i)%float_format)//',:))'
         !
         !   if (plot_files(i)%is_opened) then
         !      write(plot_files(i)%unit, fmt) (grid_point(j), j = 1, size(grid_point)), &
         !                                     (grid_function(j), j = 1, size(grid_function))
         !   end if
         !
         !end do
 
   end subroutine

! ------------------------------------------------------------------------------
  subroutine open_txt_file(plot_file)
    type(visual_plot_file_t), intent(inout):: plot_file
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
    type(visual_plot_file_t),intent(inout):: plot_file
    if (plot_file%is_opened) close(plot_file%unit, status='keep')
  end subroutine

! ------------------------------------------------------------------------------
  subroutine init_txt_file(plot_file, plot_domain, file_name)
    type(visual_plot_file_t), intent(out):: plot_file
    character(*), intent(in):: plot_domain
    character(*), intent(in), optional:: file_name

    plot_file%unit = newunit()
    if (present(file_name)) then
      plot_file%name = trim(file_name)
    else
      plot_file%name = 'plot_file'
    end if
    if (plot_file%domain == 'cube') then
      plot_file%format = 'cube'
    else
      plot_file%format = 'txt'
    end if
    plot_file%name = plot_file%name//'.'//plot_file%format
    plot_file%data_type = 'scalar'
    plot_file%float_format= 'es20.10e3'

  end subroutine

! ------------------------------------------------------------------------------
  subroutine display_visual_plot_file(self)
    class(visual_plot_file_t), intent(in):: self

    write(*,*)     
    write(*,*) self%unit
    write(*,*) self%name
    write(*,*) self%format
    write(*,*) self%domain
    write(*,*) self%data_type
    write(*,*) self%float_format
    write(*,*) self%exists
    write(*,*) self%is_opened
    write(*,*)     
  end subroutine

! ------------------------------------------------------------------------------
  subroutine get_grid_info_from_txt_storage(plot_file, grid, &
                                            plot_domain, file_name)
    type(visual_plot_file_t), intent(inout):: plot_file
    type(grid_t), intent(inout):: grid
    character(*), intent(in):: plot_domain
    character(*), intent(in), optional:: file_name
    integer:: i

    call init_txt_file(plot_file, plot_domain, file_name)
    call open_txt_file(plot_file)
    call close_txt_file(plot_file)

    grid%nr_points_total_known = .true.
  end subroutine

end module
