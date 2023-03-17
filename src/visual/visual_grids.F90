module visual_grids


   use visual_cfg
   use visual_utils
   use visual_datatypes

   implicit none

!  public subroutines
   public open_grid_files
   public close_grid_files
   public generate_grid_file
   public print_grid_files_info
   public read_grid_files
!  temp
   public open_old_grid_txt_file

!  public types
   type, public:: visual_grid_file
    integer:: unit = -1
    character(len=:), allocatable:: name
    character(len=:), allocatable:: format
    character(len=:), allocatable:: domain
    character(len=:), allocatable:: float_format
    character(len=:), allocatable:: grid_units
    logical:: exists = .false.                                                                           
    logical:: is_opened = .false.                                                                        
    !contains
    !    procedure:: assign, get, length
    !    generic:: assignment(=) => assign
    !    generic:: operator(-) => get
    !    generic:: l => length
   end type visual_grid_file

   interface generate_grid_file
      module procedure generate_grid_file_cube_regular
   end interface

   character(len = 16):: float_format_txt = '(*(es18.8e3, :))'

   private

#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbxlr.h"
#include "codata.h"
#include "orgcom.h"

   contains



!  -----------------------------------------------------------------------------
   subroutine generate_grid_file_cube_regular(grid_files, nr_points)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout), dimension(:):: grid_files
      integer, intent(in):: nr_points(3)
      real(8), allocatable:: grid_points(:,:)
      integer:: i, icoor, jcoor, kcoor, idx, n
      character(len = 16):: fmt
!  -----------------------------------------------------------------------------

      n = nr_points(1)*nr_points(2)*nr_points(3)
      allocate(grid_points(n, 3))
      grid_points = 0.0d0

      call get_grid_points(grid_points, nr_points)

      do i = 1, visual_cfg_nr_grid_files
      
         !fmt = '(*('//trim(grid_files(i)%float_format)//',:))'
         fmt = float_format_txt 

         if (grid_files(i)%is_opened) then
            call grid_file_txt_write(grid_files(i)%unit, & 
                                     fmt, &
                                     n, &
                                     grid_points(:,1), &
                                     grid_points(:,2), &
                                     grid_points(:,3))
         else
            write(*, *) 'error-grid file closed'
         end if
      
      end do

      deallocate(grid_points)

   end subroutine generate_grid_file_cube_regular


!  -----------------------------------------------------------------------------
   subroutine get_grid_points(grid_points, nr_points)
!  -----------------------------------------------------------------------------
      real(8), intent(inout), dimension(:,:):: grid_points
      integer, intent(in):: nr_points(3)
      real(8), allocatable:: atom_xyz(:, :)
      real(8)              :: step_3d(3, 0:3), r(3)
      integer              :: i, j, k, l, m, istep, n
      integer              :: icube, jcube, kcube
!  -----------------------------------------------------------------------------

         allocate(atom_xyz(4, nucdep))

!        step_3d(*, 0)  coordinates of initial point
!        step_3d(*, 1)  step vector for the "x"-direction
!        step_3d(*, 2)  step vector for the "y"-direction
!        step_3d(*, 3)  step vector for the "z"-direction

!        atom_xyz(1, i) x-coordinate of atom i
!        atom_xyz(2, i) y-coordinate of atom i
!        atom_xyz(3, i) z-coordinate of atom i
!        atom_xyz(4, i) charge of atom i

         step_3d = 0.0d0

         do i = 1, 3
            step_3d(i, 0) = cord(i, 1)
            step_3d(i, i) = cord(i, 1)
         enddo

         i = 0
         do j = 1, nucind
            k = istbnu(j)
            do l = 0, maxopr
               if (iand(l, k) == 0) then
                  i = i+1
                  do m = 1, 3
                     atom_xyz(m, i) = pt(iand(isymax(m, 1), l))*cord(m, j)

                     step_3d(m, 0) = min(step_3d(m, 0), atom_xyz(m, i))
                     step_3d(m, m) = max(step_3d(m, m), atom_xyz(m, i))
                  end do
                  atom_xyz(4, i) = charge(j)
               end if
            end do
         end do

!        add some space, calculate step sizes and convert step vectors to angstrom

         do i = 1, 3
            step_3d(i, 0) = (step_3d(i, 0) - (visual_cfg_add_3d*2.0d0))*xtang
            step_3d(i, i) = (step_3d(i, i) + (visual_cfg_add_3d*2.0d0))*xtang
            step_3d(i, i) = (step_3d(i, i) -  step_3d(i, 0))/(visual_cfg_ncube(i) - 1)
         end do

         i = 1
         do icube = 1, nr_points(1)
            do jcube = 1, nr_points(2)
               do kcube = 1, nr_points(3)
                  r(1) = step_3d(1, 0) + (icube-1)*step_3d(1, 1) &
                                       + (jcube-1)*step_3d(1, 2) &
                                       + (kcube-1)*step_3d(1, 3)
                  r(2) = step_3d(2, 0) + (icube-1)*step_3d(2, 1) &
                                       + (jcube-1)*step_3d(2, 2) &
                                       + (kcube-1)*step_3d(2, 3)
                  r(3) = step_3d(3, 0) + (icube-1)*step_3d(3, 1) &
                                       + (jcube-1)*step_3d(3, 2) &
                                       + (kcube-1)*step_3d(3, 3)
                  grid_points(i, 1) = r(1)
                  grid_points(i, 2) = r(2)
                  grid_points(i, 3) = r(3)
                  i = i+1
               end do
            end do
         end do

         deallocate(atom_xyz)
 
   end subroutine






!  -----------------------------------------------------------------------------
   subroutine open_grid_files(grid_files)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout):: grid_files(:)

      integer:: lunit_3d_gridimp
      integer:: lunit_3d_jbfile

      integer:: i
!  -----------------------------------------------------------------------------

      lunit_3d_gridimp = 32
      lunit_3d_jbfile  = 33

      if (visual_cfg_3d_gridimp) then
!        grid imported from file
         grid_files(1)%unit = lunit_3d_gridimp
         grid_files(1)%name = visual_cfg_3d_gridfil
         grid_files(1)%format = 'txt'
         grid_files(1)%domain = '3d'
         grid_files(1)%grid_units = visual_cfg_3d_grid_units
         ! open

      else if (visual_cfg_3d_readjb) then
!        grid and data imported from file
!gosia fixme-move this elsewhere
         grid_files(1)%unit = lunit_3d_jbfile
         grid_files(1)%name = visual_cfg_3d_jbfile
         grid_files(1)%format = 'txt'
         grid_files(1)%domain = '3d'
         grid_files(1)%grid_units  = visual_cfg_3d_grid_units  ! fixme
         ! open

      else
         do i = 1, visual_cfg_nr_grid_files
!           open a new file for grid generated in this run
            grid_files(i)%unit = get_new_file_unit()
            grid_files(i)%name = 'testme.grid'
            grid_files(i)%format = 'txt'
            grid_files(i)%domain = '3d'
            grid_files(i)%float_format = float_format_txt
            grid_files(i)%grid_units  = 'au' ! visual_cfg_3d_grid_units  ! fixme
            call open_new_grid_txt_file(grid_files(i))
         end do

      end if


   end subroutine


!  -----------------------------------------------------------------------------
   subroutine open_new_grid_txt_file(f)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout):: f
!  -----------------------------------------------------------------------------

      open(unit   = f%unit, &
           file   = f%name, &
           status = 'new',          &
           access = 'sequential',   &
           form   = 'formatted')
      f%is_opened = .true.
      f%exists    = .true.
 
   end subroutine


!  -----------------------------------------------------------------------------
   subroutine open_old_grid_txt_file(f)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout):: f
      logical:: ex, op
!  -----------------------------------------------------------------------------

      !inquire(file = f%name, exist = ex, opened = op)
      !if (ex .and. .not. op) then
      if (f%exists .and..not.f%is_opened) then
         f%unit = get_new_file_unit()
         open(unit   = f%unit, &
              file   = f%name, &
              status = 'OLD',          &
              access = 'sequential',   &
              form   = 'formatted')
         f%is_opened = .true.
      end if
      write(*,*) 'opened file: ', f%name, f%unit
   end subroutine


!  -----------------------------------------------------------------------------
   subroutine close_grid_files(grid_files)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout), dimension(:):: grid_files
      integer:: i
!  -----------------------------------------------------------------------------

      do i = 1, size(grid_files)
         if (grid_files(i)%is_opened) close(grid_files(i)%unit, status='keep')
      end do

   end subroutine


!  -----------------------------------------------------------------------------
   subroutine read_grid_files(x, y, z, grid_files, nr_points)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(in), dimension(:), allocatable:: grid_files
      real(8), intent(inout), dimension(:):: x, y, z
      integer, intent(in):: nr_points(3)
      integer:: i, j, n
!  -----------------------------------------------------------------------------

      n = nr_points(1)*nr_points(2)*nr_points(3)

      !do i = 1, visual_cfg_nr_grid_files
      !
      !   if (grid_files(i)%is_opened) then
      !      rewind(grid_files(i)%unit)
      !   else
      !      write(*, *) 'error-grid file closed'
      !   end if

      !   call grid_file_txt_read(grid_files(i)%unit, grid_files(i)%float_format, n, x, y, z)
      !end do
        i = 1
            rewind(grid_files(i)%unit)
         call grid_file_txt_read(grid_files(i)%unit, grid_files(i)%float_format, n, x, y, z)

   end subroutine


!  -----------------------------------------------------------------------------
   subroutine print_grid_files_info(grid_files)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(in), dimension(:):: grid_files
      integer:: i
!  -----------------------------------------------------------------------------

      write(*,*) 'Grid files used in VISUAL:'
      write(*,*) '--------------------------'

      do i = 1, size(grid_files)
         write(*,*) ' grid file #',i    
         write(*,*) grid_files(i)%unit
         write(*,*) grid_files(i)%name
         write(*,*) grid_files(i)%format
         write(*,*) grid_files(i)%domain
         write(*,*) grid_files(i)%float_format
         write(*,*) grid_files(i)%grid_units
         write(*,*) grid_files(i)%exists
         write(*,*) grid_files(i)%is_opened
      end do

   end subroutine print_grid_files_info



!XXXXXXXXXXXXXXXXXXXXXXXXX


   subroutine grid_file_txt_write(file_unit, float_format, nr_points, x, y, z)
!     --------------------------------------------------------------------------
      integer, intent(in):: file_unit
      character(*), intent(in):: float_format
      integer, intent(in):: nr_points
      real(8), intent(in):: x(:), y(:), z(:)
      character(len = 8):: xc, yc, zc
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      do i = 1, nr_points
         write(file_unit, float_format) x(i), y(i), z(i)
      enddo

   end subroutine


   subroutine grid_file_txt_read(file_unit, float_format, nr_points, x, y, z)
!     --------------------------------------------------------------------------
      real(8), intent(inout), dimension(:):: x, y, z
      integer, intent(in)  :: file_unit
      character(*), intent(in):: float_format
      integer, intent(in)  :: nr_points
!     --------------------------------------------------------------------------
      integer              :: i
!     --------------------------------------------------------------------------

      rewind(file_unit)

      i = 1
      do
         read(file_unit, *, end = 98) x(i), y(i), z(i)
         write(*, *) 'in grid_file_txt_read: ', i, x(i), y(i), z(i)
         i = i+1
      end do
 98   continue


   end subroutine








end module
