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
      type(visual_grid_file), intent(inout), dimension(:), allocatable:: grid_files
      integer, intent(in):: nr_points(3)
      real(8), allocatable:: grid_points(:,:)
      integer:: i, icoor, jcoor, kcoor, idx, n
!  -----------------------------------------------------------------------------

      n = maxval(nr_points)
      allocate(grid_points(n, 3))
      grid_points = 0.0d0

      call get_grid_points(grid_points, nr_points)

      do i = 1, visual_cfg_nr_grid_files
      
         !fmt = '(*('//trim(grid_files(i)%float_format)//',:))'

         if (grid_files(i)%is_opened) then
            do icoor = 1, nr_points(1)
               do jcoor = 1, nr_points(2)
                  do kcoor = 1, nr_points(3)
                     write(grid_files(i)%unit, float_format_txt) &
                        grid_points(icoor, 1), &
                        grid_points(jcoor, 2), &
                        grid_points(kcoor, 3)
                  end do
                end do
            end do
         end if
      
      end do

      deallocate(grid_points)

   end subroutine generate_grid_file_cube_regular


!  -----------------------------------------------------------------------------
   subroutine get_grid_points(grid_points, nr_points)
!  -----------------------------------------------------------------------------
!#include "nuclei.h" !nucdep, cord, nucind, istbnu, charge 
!#include "symmet.h" !maxopr, pt, isymax
!#include "codata.h" !xtang
      real(8), intent(inout), dimension(:,:):: grid_points
      integer, intent(in):: nr_points(3)
      real(8), allocatable:: atom_xyz(:, :)
      real(8)              :: step_3d(3, 0:3), r(3)
      integer              :: i, j, k, l, m, istep, n
      integer              :: icube, jcube, kcube
      integer:: nr_points_total
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
                  grid_points(icube, 1) = r(1)
                  grid_points(jcube, 2) = r(2)
                  grid_points(kcube, 3) = r(3)
               end do
            end do
         end do

         deallocate(atom_xyz)
 
   end subroutine






!  -----------------------------------------------------------------------------
   subroutine open_grid_files(grid_files)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(out), allocatable:: grid_files(:)

      integer:: lunit_3d_gridimp
      integer:: lunit_3d_jbfile

      integer:: i
!  -----------------------------------------------------------------------------

      lunit_3d_gridimp = 32
      lunit_3d_jbfile  = 33

      if (visual_cfg_3d_gridimp) then
!        grid imported from file
         allocate(grid_files(1))
         grid_files(1)%unit = lunit_3d_gridimp
         grid_files(1)%name = visual_cfg_3d_gridfil
         grid_files(1)%format = 'txt'
         grid_files(1)%domain = '3d'
         grid_files(1)%grid_units = visual_cfg_3d_grid_units

      else if (visual_cfg_3d_readjb) then
!        grid and data imported from file
!gosia fixme-move this elsewhere
         allocate(grid_files(1))
         grid_files(1)%unit = lunit_3d_jbfile
         grid_files(1)%name = visual_cfg_3d_jbfile
         grid_files(1)%format = 'txt'
         grid_files(1)%domain = '3d'
         grid_files(1)%grid_units  = visual_cfg_3d_grid_units  ! fixme

      else
!        open a new file for grid generated in this run
         allocate(grid_files(1))
         grid_files(1)%unit = get_new_file_unit()
         grid_files(1)%name = 'testme.grid'
         grid_files(1)%format = 'txt'
         grid_files(1)%domain = '3d'
         grid_files(1)%float_format = float_format_txt
         grid_files(1)%grid_units  = 'au' ! visual_cfg_3d_grid_units  ! fixme

      end if

      if (allocated(grid_files)) then
         do i = 1, size(grid_files)
            call open_inp_txt_file(grid_files(i))
         end do
      end if


   end subroutine

!  -----------------------------------------------------------------------------
   subroutine open_inp_txt_file(f)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout):: f
      character(len = 100):: errmsg
      logical:: exists, is_opened
!  -----------------------------------------------------------------------------

      inquire(f%unit, exist = exists, opened = is_opened)
      if (exists) then
         if (.not.is_opened) then
            open(unit   = f%unit, &
                 file   = f%name, &
                 status = 'new',          &
                 access = 'sequential',   &
                 form   = 'formatted')
            f%is_opened = .true.
         end if
         f%exists    = .true.
      else
          write (errmsg, *) 'file ', f%name, ' does not exist',&
                            ' Stopping now!'
          call visual_quit(errmsg)
      end if
 
   end subroutine open_inp_txt_file



!  -----------------------------------------------------------------------------
   subroutine close_grid_files(grid_files)
!  -----------------------------------------------------------------------------
      type(visual_grid_file), intent(inout), dimension(:), allocatable:: grid_files
      integer:: i
!  -----------------------------------------------------------------------------

      do i = 1, size(grid_files)
         if (grid_files(i)%is_opened) close(grid_files(i)%unit, status='keep')
      end do

      deallocate(grid_files)

   end subroutine


   subroutine print_grid_files_info(grid_files)
      type(visual_grid_file), intent(in), dimension(:):: grid_files
      integer:: i

      write(*,*) 'size(grid_files): ', size(grid_files)
      do i = 1, size(grid_files)
         write(*,*)     
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
end module
