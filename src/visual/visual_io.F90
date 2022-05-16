module visual_io


   use visual_cfg
   use visual_utils
   use visual_datatypes

   implicit none

!  public subroutines
   public open_plot_files
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
   subroutine open_plot_files(use_legacy, plot_files)
!  -----------------------------------------------------------------------------
      logical, intent(in):: use_legacy
      type(visual_plot_file), intent(out), allocatable:: plot_files(:)
      integer:: i
!  -----------------------------------------------------------------------------

      if (use_legacy) then
         call visual_legacy_files(plot_files)
      else
         call prepare_plot_txt_file(plot_files)
      end if

      do i = 1, size(plot_files)
         call open_plot_txt_file(plot_files(i))
      end do

   end subroutine open_plot_files


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
   subroutine visual_legacy_files(plot_files)
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
      integer:: lunit_3d_jbfile

      integer:: lunit_line_scalar
      integer:: lunit_line_vector
      integer:: lunit_line_tensor
      integer:: lunit_3d_gridimp
      integer:: lunit_radial_scalar
      integer:: lunit_radial_vector
      integer:: lunit_radial_tensor
      integer:: lunit_outfil_custom_scalar
      integer:: lunit_outfil_custom_vector
      integer:: lunit_outfil_custom_2ndtensor

      integer:: i
!  -----------------------------------------------------------------------------

      lunit_2d_scalar     = 14
      lunit_2d_vector     = 15
      lunit_2d_flow       = 16
      lunit_2d_tensor     = 27

      lunit_3d_cube       = 17
      lunit_3d_scalar     = 18
      lunit_3d_vector     = 19
      lunit_3d_vector_vtk = 20
      lunit_3d_tensor     = 28
      lunit_3d_gridimp    = 23
      lunit_3d_jbfile     = 26

      lunit_line_scalar   = 21
      lunit_line_vector   = 22
      lunit_line_tensor   = 29

      lunit_radial_scalar = 24
      lunit_radial_vector = 25
      lunit_radial_tensor = 30

      lunit_outfil_custom_scalar    = 31
      lunit_outfil_custom_vector    = 32
      lunit_outfil_custom_2ndtensor = 33

      !filename_3d_cube       = 'plot.3d.cube'
      !filename_3d_scalar     = 'plot.3d.scalar'
      !filename_3d_vector     = 'plot.3d.vector'
      !filename_3d_tensor     = 'plot.3d.tensor'
      !filename_3d_vector_vtk = 'plot.3d.vector.vtk'

      !filename_line_scalar   = 'plot.line.scalar'
      !filename_line_vector   = 'plot.line.vector'
      !filename_line_tensor   = 'plot.line.tensor'

      !filename_radial_scalar = 'plot.radial.scalar'
      !filename_radial_vector = 'plot.radial.vector'
      !filename_radial_tensor = 'plot.radial.tensor'

      if (visual_cfg_2d) then

         allocate(plot_files(size(visual_cfg_2d_data_types)))

         ! data_type = 'scalar'
         plot_files(1)%unit = lunit_2d_scalar 
         plot_files(1)%name = 'plot.2d.scalar'
         plot_files(1)%data_type = 'scalar'

         ! data_type = 'vector'
         plot_files(2)%unit = lunit_2d_vector 
         plot_files(2)%name = 'plot.2d.vector'
         plot_files(2)%data_type = 'vector'

         ! data_type = 'tensor'
         plot_files(3)%unit = lunit_2d_tensor
         plot_files(3)%name = 'plot.2d.tensor'
         plot_files(3)%data_type = 'tensor'

         ! data_type = 'flow'
         plot_files(4)%unit = lunit_2d_flow
         plot_files(4)%name = 'plot.2d.flow'
         plot_files(4)%data_type = 'flow'

         ! file attributes common to all data types
         do i = 1, size(visual_cfg_2d_data_types)
            plot_files(i)%format    = 'txt'
            plot_files(i)%domain  = '2d'
         end do   

      end if

   end subroutine visual_legacy_files


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
