module visual_txt_storage

  use visual_utils
  use visual_datatypes, only: rs_label_t, grid_t, grid_function_t, mol_t

  implicit none

  public write_plot_data_to_DIRAC_output
  !public write_grid_to_txt_storage
  public write_grid_function_to_txt_storage
  public write_grid_function_to_cube
  public get_grid_info_from_txt_storage
  public get_grid_points_from_txt_storage

  type, public:: txt_file_info_t
    integer:: unit = -1
    character(len=:), allocatable:: name
    character(len=:), allocatable:: format
    character(len=:), allocatable:: float_format
    logical:: exists = .false.
    logical:: is_opened = .false.
  end type

  character(len = 16):: float_fmt = '(*(es18.8e3, :))'
  private

contains


!! -----------------------------------------------------------------------------
!  subroutine write_grid_to_txt_storage(grid)
!    type(grid_t), intent(in):: grid
!    type(txt_file_info_t):: data_file
!    integer:: g_ndim, gf_ndim, npoints, npoints_radial
!    integer:: i, j, funit
!    character(len=:),allocatable:: float_fmt_local
!
!    npoints = grid%nr_points_total
!    npoints_radial = grid%radial_nr_steps
!    g_ndim  = grid%get_ndim()
!    data_file%name = grid%file_out
!    data_file%float_format = float_fmt
!
!    call open_txt_file(data_file, 'w')
!    funit = data_file%unit
!
!    if (grid%export == 'txt') then
!
!      ! this is for backward compatibility:
!      ! data is written in columns, separator=space, no header
!
!      if (grid%allocated_radial_coord()) then
!        do j = 1, npoints_radial
!          write(funit, float_fmt) (grid%radial_coord(j))
!        end do
!      else
!        if (grid%allocated_weights()) then
!          do j = 1, npoints
!            write(funit, float_fmt) (grid%points(i,j),i=1,g_ndim), &
!                                     grid%weights(j)
!          end do
!        else
!          do j = 1, npoints
!            write(funit, float_fmt) (grid%points(i,j),i=1,g_ndim)
!          end do
!        end if
!      end if
!
!    else if (grid%export == 'csv') then
!
!      ! data is written in columns, separator=',', with header 
!
!      if (grid%allocated_radial_coord()) then
!        write(funit, *) 'radial_step'
!      else
!        if (grid%allocated_weights()) then
!          write(funit, *) 'x,y,z,w'
!        else
!          write(funit, *) 'x,y,z'
!        end if
!      end if
!      if (grid%allocated_radial_coord()) then
!        do j = 1, npoints_radial
!          write(funit, float_fmt) (grid%radial_coord(j))
!        end do
!      else
!        if (grid%allocated_weights()) then
!          do j = 1, npoints
!            if (g_ndim == 3) then
!              float_fmt_local='(*(es18.8e3,:,","))'
!              write(funit, float_fmt_local) grid%points(1,j), &
!                                            grid%points(2,j), &
!                                            grid%points(3,j), &
!                                            grid%weights(j)
!            end if
!          end do
!        else
!          do j = 1, npoints
!            if (g_ndim == 3) then
!              float_fmt_local='(*(es18.8e3,:,","))'
!              write(funit, float_fmt_local) grid%points(1,j), &
!                                            grid%points(2,j), &
!                                            grid%points(3,j)
!            end if
!          end do
!        end if
!      end if
!    end if
!
!    call close_txt_file(data_file, 'keep')
!   end subroutine


! -----------------------------------------------------------------------------
  subroutine write_grid_function_to_txt_storage(grid,grid_function)
    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(in):: grid_function
    type(txt_file_info_t):: data_file
    integer:: g_ndim, gf_ndim, npoints, npoints_radial
    integer:: i, j, funit
    character(len=:),allocatable:: float_fmt_local

    npoints = grid%nr_points_total
    npoints_radial = grid%radial_nr_steps
    g_ndim  = grid%get_ndim()
    gf_ndim = grid_function%get_ndim()
    data_file%name = grid_function%file_out
    data_file%float_format = float_fmt

    call open_txt_file(data_file, 'w')
    funit = data_file%unit

    if (grid_function%export == 'txt') then

      ! this is for backward compatibility:
      ! data is written in columns, separator=space, no header

      if (grid%allocated_radial_coord()) then
        do j = 1, npoints_radial
          write(funit, float_fmt) (grid%radial_coord(j)), &
                                  (grid_function%values(i,j),i=1,gf_ndim)
        end do
      else
        if (grid%allocated_weights()) then
          do j = 1, npoints
            write(funit, float_fmt) (grid%points(i,j),i=1,g_ndim), &
                                     grid%weights(j), &
                                    (grid_function%values(i,j),i=1,gf_ndim)
          end do
        else
          do j = 1, npoints
            write(funit, float_fmt) (grid%points(i,j),i=1,g_ndim), &
                                    (grid_function%values(i,j),i=1,gf_ndim)
          end do
        end if
      end if

    else if (grid_function%export == 'csv') then

      ! data is written in columns, separator=',', with header 

      if (grid%allocated_radial_coord()) then
        write(funit, *) 'radial_step'
      else
        if (grid%allocated_weights()) then
          write(funit, *) 'x,y,z,w,', grid_function%name
        else
          write(funit, *) 'x,y,z,', grid_function%name
        end if
      end if
      if (grid%allocated_radial_coord()) then
        do j = 1, npoints_radial
          float_fmt_local='(*(es18.8e3,:,","))'
          write(funit, float_fmt_local) (grid%radial_coord(j)), &
                                        (grid_function%values(i,j),i=1,gf_ndim)
        end do
      else
        if (grid%allocated_weights()) then
          do j = 1, npoints
            if (g_ndim == 3) then
              float_fmt_local='(*(es18.8e3,:,","))'
              write(funit, float_fmt_local) grid%points(1,j), &
                                            grid%points(2,j), &
                                            grid%points(3,j), &
                                            grid%weights(j), &
                                            (grid_function%values(i,j),i=1,gf_ndim)
            end if
          end do
        else
          do j = 1, npoints
            if (g_ndim == 3) then
              float_fmt_local='(*(es18.8e3,:,","))'
              write(funit, float_fmt_local) grid%points(1,j), &
                                            grid%points(2,j), &
                                            grid%points(3,j), &
                                            (grid_function%values(i,j),i=1,gf_ndim)
            end if
          end do
        end if
      end if
    end if

    call close_txt_file(data_file, 'keep')
   end subroutine



! -----------------------------------------------------------------------------
  subroutine write_grid_function_to_cube(grid, grid_function, mol)
    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(in):: grid_function
    type(mol_t), intent(in):: mol
    type(txt_file_info_t):: data_file
    integer:: g_ndim, gf_ndim, npoints, npoints_radial
    integer:: funit, i, j, ix, iy, iz, ioff

!   copied and adapted from old VISUAL

    data_file%name = grid_function%file_out

    call open_txt_file(data_file, 'w')
    funit = data_file%unit

!   write header

    write(funit, '(a)') 'cube file'
    write(funit, '(a)') 'generated by DIRAC'

!   Origin and grid size
    write(funit, '(i5, 4f12.6)') mol%n_atoms, grid%origin(1:3)
    write(funit, '(i5, 4f12.6)') grid%nr_points(1),grid%spacing(1),0.d0,0.0d0
    write(funit, '(i5, 4f12.6)') grid%nr_points(2),0.0d0,grid%spacing(2),0.d0
    write(funit, '(i5, 4f12.6)') grid%nr_points(3),0.d0,0.0d0,grid%spacing(3)

!   Atom info
    do i = 1, mol%n_atoms
      write(funit, '(i5, 4f12.6)') &
        int(mol%nuc_charge(i)), mol%nuc_charge(i), (mol%geometry(j, i), j = 1, 3)
    end do

!   Data block: Write values in z-fastest order
    ioff = 0
    do ix = 1, grid%nr_points(1)
      do iy =1, grid%nr_points(2)
        write(funit, '(6ES13.5E3)') (grid_function%values(1,i), i=ioff+1,ioff+grid%nr_points(3))
        ioff = ioff+grid%nr_points(3)
      end do
    end do


    call close_txt_file(data_file, 'keep')
  end subroutine

! -----------------------------------------------------------------------------
  subroutine write_plot_data_to_DIRAC_output(grid, grid_function)
    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(in):: grid_function
    integer:: g_ndim, gf_ndim, npoints, npoints_radial
    integer:: i, j

    npoints = grid%nr_points_total
    npoints_radial = grid%radial_nr_steps
    g_ndim  = grid%get_ndim()
    gf_ndim = grid_function%get_ndim()

    if (grid_function%purpose == 1) then 
      write(*, '(2x,a,i3,a,i3,a)') &
        'Final grid function values (id=', grid_function%id,') integrated on grid (id=', grid%id,'): '
      write(*, float_fmt) (grid_function%integrated_values(i),i=1,gf_ndim)

    else if (grid_function%purpose == 0) then
      write(*, '(2x,a,i3,a,i3,a)') &
        'Final grid function values (id=', grid_function%id,') in grid points (id=', grid%id,'): '
      if (grid%allocated_radial_coord()) then
        write(*, '(4x,a)') 'radial step, grid function values'
      else
        if (grid%allocated_weights()) then
          write(*, '(4x,a)') 'x  y  z  w  grid_function values'
        else
          write(*, '(4x,a)') 'x  y  z  grid_function values'
        end if
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
    write(*, *) ''

   end subroutine

! ------------------------------------------------------------------------------
  subroutine open_txt_file(plot_file, mode)
    type(txt_file_info_t), intent(inout):: plot_file
    character(len = 1), intent(in):: mode
    character(len = 100):: errmsg
    logical:: exists, is_opened

    inquire(file=plot_file%name, exist=exists, opened=is_opened)

    if (mode == 'r') then
      ! open file to read
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
    else if (mode == 'w') then
      ! open new file to write
      plot_file%unit = newunit()
      open(unit   = plot_file%unit, &
           file   = plot_file%name, &
           status = 'new',      &
           access = 'sequential',   &
           form   = 'formatted')
      plot_file%is_opened = .true.
      plot_file%exists    = .true.
    else
      call visual_quit('open_txt_file: unsupported mode')
    end if
 
  end subroutine

! ------------------------------------------------------------------------------
  subroutine close_txt_file(plot_file, st)
    type(txt_file_info_t),intent(inout):: plot_file
    character(len=*),intent(in):: st

    if ((st == 'keep') .or. (st == 'delete')) then
      if (plot_file%is_opened) then
        close(plot_file%unit, status=st)
      end if
    else
      call visual_quit('close_txt_file: unsupported status')
    end if
  end subroutine

! ------------------------------------------------------------------------------
  subroutine get_grid_info_from_txt_storage(grid)
    type(txt_file_info_t):: grid_file
    type(grid_t), intent(inout):: grid
    integer:: ival, nr_points, i_point
    character*80 :: line


    if (grid%file_inp /= 'numerical_grid') return

    grid_file%name = 'numerical_grid'
    grid_file%unit = newunit()

    call open_txt_file(grid_file, 'r')
    rewind(grid_file%unit)

    nr_points = 0
    !we may have multiple grids on this file, loop over them to read all points
    do
       read(grid_file%unit, *) ival
       if (ival > 0) then
          !skip over lines with grid points
          do i_point = 1, ival
             read (grid_file%unit, *) line
          end do
          nr_points = nr_points + ival
       else ! -1 indicates the end of the list of grids
          exit
       end if
    end do

    grid%nr_points_total = nr_points
    if (grid%nr_points_total > 0) grid%nr_points_total_known = .true.
    
    call close_txt_file(grid_file, 'keep')

  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_grid_points_from_txt_storage(grid)
use interface_grid
    type(txt_file_info_t):: grid_file
    type(grid_t), intent(inout):: grid
    integer:: ioff, nr_points
    character(len = 100):: errmsg

    if (grid%file_inp /= 'numerical_grid') return

    grid_file%name = 'numerical_grid'
    grid_file%unit = newunit()

    call open_txt_file(grid_file, 'r')
    rewind(grid_file%unit)

    ioff = 0
    do
      read(grid_file%unit, *) nr_points
      if (nr_points < 0) exit
      call num_grid_read(grid%points(1,ioff+1:ioff+nr_points),&
                         grid%points(2,ioff+1:ioff+nr_points),&
                         grid%points(3,ioff+1:ioff+nr_points),&
                         grid%weights(ioff+1:ioff+nr_points),&
                         grid_file%unit, nr_points)
      ioff = ioff + nr_points
    end do

    call close_txt_file(grid_file, 'keep')

  end subroutine

end module
