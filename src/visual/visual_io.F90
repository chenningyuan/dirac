module visual_io


   use visual_cfg
   use visual_utils
   use visual_datatypes

   implicit none

!  public subroutines
   public open_plot_files_on_import
   public open_plot_files_for_export
   public close_plot_txt_files
   public print_plot_files_info

!  public types
   type, public:: visual_plot_file
    integer:: unit = -1
    character(len=:), allocatable:: name
    character(len=:), allocatable:: format                                                               
    character(len=:), allocatable:: domain                                                             
    character(len=:), allocatable:: data_type                                                            
    logical:: exists = .false.                                                                           
    logical:: is_opened = .false.                                                                        
   end type visual_plot_file

   private

   contains


!  -----------------------------------------------------------------------------
   subroutine open_plot_files_on_import(use_legacy, plot_files)
!  -----------------------------------------------------------------------------
      logical, intent(in):: use_legacy
      type(visual_plot_file), intent(out), allocatable:: plot_files(:)
      integer:: i
!  -----------------------------------------------------------------------------

      !if (use_legacy) then
      !   call visual_legacy_files_for_data_import()
      !else
      !   call prepare_plot_txt_file(plot_files)
      !end if

      !if (allocated(plot_files)) then
      !   do i = 1, size(plot_files)
      !      call open_plot_txt_file(plot_files(i))
      !   end do
      !end if

   end subroutine open_plot_files_on_import


!  -----------------------------------------------------------------------------
   subroutine open_plot_files_for_export(use_legacy, plot_files)
!  -----------------------------------------------------------------------------
      logical, intent(in):: use_legacy
      type(visual_plot_file), intent(out), allocatable:: plot_files(:)
      integer:: i
!  -----------------------------------------------------------------------------

      if (use_legacy) then
         call visual_legacy_files_for_data_export(plot_files)
      else
         call prepare_plot_txt_file(plot_files)
      end if

      if (allocated(plot_files)) then
         do i = 1, size(plot_files)
            call open_plot_txt_file(plot_files(i))
         end do
      end if

   end subroutine open_plot_files_for_export


!  -----------------------------------------------------------------------------
   subroutine open_plot_txt_file(plot_file)
!  -----------------------------------------------------------------------------
      type(visual_plot_file), intent(inout):: plot_file
!  -----------------------------------------------------------------------------

      open(unit   = plot_file%unit, &
           file   = plot_file%name, &
           status = 'new',          &
           access = 'sequential',   &
           form   = 'formatted')

      plot_file%exists    = .true.
      plot_file%is_opened = .true.
 
   end subroutine open_plot_txt_file


!  -----------------------------------------------------------------------------
   subroutine close_plot_txt_files(plot_files)
!  -----------------------------------------------------------------------------
      type(visual_plot_file), intent(inout), dimension(:), allocatable:: plot_files
      integer:: i
!  -----------------------------------------------------------------------------

      do i = 1, size(plot_files)
         if (plot_files(i)%is_opened) close(plot_files(i)%unit, status='keep')
      end do

      deallocate(plot_files)

   end subroutine close_plot_txt_files


!  -----------------------------------------------------------------------------
   subroutine prepare_plot_txt_file(plot_files)
!  -----------------------------------------------------------------------------
      type(visual_plot_file), intent(out), allocatable:: plot_files(:)
!  -----------------------------------------------------------------------------

      allocate(plot_files(1))

      plot_files(1)%unit = newunit()
      plot_files(1)%name = 'plot_file.txt'
      plot_files(1)%data_type = 'scalar'
      plot_files(1)%format    = 'txt'
      plot_files(1)%domain  = '2d'

   end subroutine prepare_plot_txt_file


!  -----------------------------------------------------------------------------
   subroutine visual_legacy_files_for_data_export(plot_files)
!  -----------------------------------------------------------------------------
      type(visual_plot_file), intent(out), allocatable:: plot_files(:)

      integer:: lunit_2d_scalar
      integer:: lunit_2d_vector
      integer:: lunit_2d_flow
      integer:: lunit_2d_tensor

      integer:: lunit_3d_cube
      integer:: lunit_3d_scalar
      integer:: lunit_3d_vector
      integer:: lunit_3d_tensor
      integer:: lunit_3d_vector_vtk

      integer:: lunit_line_scalar
      integer:: lunit_line_vector
      integer:: lunit_line_tensor

      integer:: lunit_list_scalar
      integer:: lunit_list_vector
      integer:: lunit_list_tensor

      integer:: lunit_radial_scalar
      integer:: lunit_radial_vector
      integer:: lunit_radial_tensor

      integer:: i, nfiles

!     prepare legacy plot files for export
!  -----------------------------------------------------------------------------

      lunit_2d_scalar     = 14
      lunit_2d_vector     = 15
      lunit_2d_flow       = 16
      lunit_2d_tensor     = 17

      lunit_3d_cube       = 18
      lunit_3d_scalar     = 19
      lunit_3d_vector     = 20
      lunit_3d_vector_vtk = 21
      lunit_3d_tensor     = 22

      lunit_line_scalar   = 23
      lunit_line_vector   = 24
      lunit_line_tensor   = 25

      lunit_list_scalar   = 26
      lunit_list_vector   = 27
      lunit_list_tensor   = 28

      lunit_radial_scalar = 29
      lunit_radial_vector = 30
      lunit_radial_tensor = 31

      if (visual_cfg_2d) then

         allocate(plot_files(size(visual_cfg_2d_data_types)))

         ! data_type = 'scalar'
         i = 1
         plot_files(i)%unit = lunit_2d_scalar 
         plot_files(i)%name = 'plot.2d.scalar'
         plot_files(i)%data_type = 'scalar'

         ! data_type = 'vector'
         i = i+1
         plot_files(i)%unit = lunit_2d_vector 
         plot_files(i)%name = 'plot.2d.vector'
         plot_files(i)%data_type = 'vector'

         ! data_type = 'tensor'
         i = i+1
         plot_files(i)%unit = lunit_2d_tensor
         plot_files(i)%name = 'plot.2d.tensor'
         plot_files(i)%data_type = 'tensor'

         ! data_type = 'flow'
         i = i+1
         plot_files(i)%unit = lunit_2d_flow
         plot_files(i)%name = 'plot.2d.flow'
         plot_files(i)%data_type = 'flow'

         nfiles = i

         ! file attributes common to all data types
         do i = 1, nfiles
            plot_files(i)%format = 'txt'
            plot_files(i)%domain = '2d'
         end do   

      end if

      if (visual_cfg_3d .or. visual_cfg_3d_fast .or. &
          visual_cfg_3d_gridimp .or. visual_cfg_3d_readjb) then

         allocate(plot_files(size(visual_cfg_3d_data_types)))

         ! data_type = 'scalar'
         i = 1
         plot_files(i)%unit = lunit_3d_scalar 
         plot_files(i)%name = 'plot.3d.scalar'
         plot_files(i)%data_type = 'scalar'

         if (.not. visual_cfg_london_skip_cube) then
            i = i+1
            plot_files(i)%unit = lunit_3d_cube
            plot_files(i)%name = 'plot.3d.cube'
            plot_files(i)%data_type = 'scalar'
         end if

         ! data_type = 'vector'
         i = i+1
         plot_files(i)%unit = lunit_3d_vector 
         plot_files(i)%name = 'plot.3d.vector'
         plot_files(i)%data_type = 'vector'

         if (.not. visual_cfg_london_skip_vtk) then
            i = i+1
            plot_files(i)%unit = lunit_3d_vector_vtk
            plot_files(i)%name = 'plot.3d.vector_vtk'
            plot_files(i)%data_type = 'vector'
         end if

         ! data_type = 'tensor'
         i = i+1
         plot_files(i)%unit = lunit_3d_tensor 
         plot_files(i)%name = 'plot.3d.tensor'
         plot_files(i)%data_type = 'tensor'

         nfiles = i

         ! file attributes common to all data types
         do i = 1, nfiles
            plot_files(i)%format = 'txt'
            plot_files(i)%domain = '3d'
         end do   

      end if

      if (visual_cfg_line) then
         allocate(plot_files(size(visual_cfg_line_data_types)))

         ! data_type = 'scalar'
         i = 1
         plot_files(i)%unit = lunit_line_scalar 
         plot_files(i)%name = 'plot.line.scalar'
         plot_files(i)%data_type = 'scalar'

         ! data_type = 'vector'
         i = i+1
         plot_files(i)%unit = lunit_line_vector 
         plot_files(i)%name = 'plot.line.vector'
         plot_files(i)%data_type = 'vector'

         ! data_type = 'tensor'
         i = i+1
         plot_files(i)%unit = lunit_line_tensor 
         plot_files(i)%name = 'plot.line.tensor'
         plot_files(i)%data_type = 'tensor'

         nfiles = i

         ! file attributes common to all data types
         do i = 1, nfiles
            plot_files(i)%format = 'txt'
            plot_files(i)%domain = 'line'
         end do   

      end if

      if (visual_cfg_list) then
         allocate(plot_files(size(visual_cfg_list_data_types)))

         ! data_type = 'scalar'
         i = 1
         plot_files(i)%unit = lunit_list_scalar 
         plot_files(i)%name = 'plot.list.scalar'
         plot_files(i)%data_type = 'scalar'

         ! data_type = 'vector'
         i = i+1
         plot_files(i)%unit = lunit_list_vector 
         plot_files(i)%name = 'plot.list.vector'
         plot_files(i)%data_type = 'vector'

         ! data_type = 'tensor'
         i = i+1
         plot_files(i)%unit = lunit_list_tensor 
         plot_files(i)%name = 'plot.list.tensor'
         plot_files(i)%data_type = 'tensor'

         nfiles = i

         ! file attributes common to all data types
         do i = 1, nfiles
            plot_files(i)%format = 'txt'
            plot_files(i)%domain = 'list'
         end do   

      end if

      if (visual_cfg_radial) then
         allocate(plot_files(size(visual_cfg_radial_data_types)))

         ! data_type = 'scalar'
         i = 1
         plot_files(i)%unit = lunit_radial_scalar 
         plot_files(i)%name = 'plot.radial.scalar'
         plot_files(i)%data_type = 'scalar'

         ! data_type = 'vector'
         i = i+1
         plot_files(i)%unit = lunit_radial_vector 
         plot_files(i)%name = 'plot.radial.vector'
         plot_files(i)%data_type = 'vector'

         ! data_type = 'tensor'
         i = i+1
         plot_files(i)%unit = lunit_radial_tensor 
         plot_files(i)%name = 'plot.radial.tensor'
         plot_files(i)%data_type = 'tensor'

         nfiles = i

         ! file attributes common to all data types
         do i = 1, nfiles
            plot_files(i)%format = 'txt'
            plot_files(i)%domain = 'radial'
         end do   

      end if


   end subroutine visual_legacy_files_for_data_export


   subroutine print_plot_files_info(plot_files)
      type(visual_plot_file), intent(in), dimension(:):: plot_files
      integer:: i

      write(*,*) 'size(plot_files): ', size(plot_files)
      do i = 1, size(plot_files)
         write(*,*)     
         write(*,*) plot_files(i)%unit
         write(*,*) plot_files(i)%name
         write(*,*) plot_files(i)%format
         write(*,*) plot_files(i)%domain
         write(*,*) plot_files(i)%data_type
         write(*,*) plot_files(i)%exists
         write(*,*) plot_files(i)%is_opened
      end do
   end subroutine print_plot_files_info
end module
