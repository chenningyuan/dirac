module visual_grids


   use visual_io
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
   public get_grid_info
   public get_grid_points

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

   interface generate_grid_rectilinear
      module procedure generate_grid_rectilinear
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

  subroutine generate_grid_rectilinear(grid)
    type(grid_t), intent(inout):: grid
    integer:: i,j,k, ipoint
    real(8):: v(3)

    if (grid%ndim /= 3) return

    ipoint = 1
    do i = 1, grid%nr_points(1)
      v(1) = grid%origin(1) + (i - 1)*grid%spacing(1)
      do j = 1, grid%nr_points(2)
        v(2) = grid%origin(2) + (j - 1)*grid%spacing(2)
        do k = 1, grid%nr_points(3)
          v(3) = grid%origin(3) + (k - 1)*grid%spacing(3)

          grid%points(1:3, ipoint) = v(1:3)

          ipoint = ipoint + 1
        end do
      end do
    end do

  end subroutine generate_grid_rectilinear


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

      call get_grid_points_legacy(grid_points, nr_points)

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
   subroutine get_grid_points_legacy(grid_points, nr_points)
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
      !write(*,*) 'opened file: ', f%name, f%unit
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

! NEW

  subroutine grid_extents_from_mol(mol, grid, extents)
#include "codata.h"
    type(mol_t), intent(in):: mol
    type(grid_t), intent(in):: grid
    real(8), intent(out):: extents(0:1, 3)
    real(8):: atoms_min(1:3), atoms_max(1:3)
    integer::i, j

    ! mol%geometry is in Angstroms, so we rescale it
    atoms_min(1:3) = mol%geometry(1:3, 1)/xtang
    atoms_max(1:3) = mol%geometry(1:3, 1)/xtang

    do i = 1, mol%n_atoms
      do j = 1, 3
        atoms_min(j) = min(mol%geometry(j, i)/xtang, atoms_min(j))
        atoms_max(j) = max(mol%geometry(j, i)/xtang, atoms_max(j))
      end do
    end do

    ! add space around
    do i = 1, 3
      atoms_min(i) = atoms_min(i) - grid%margin
      atoms_max(i) = atoms_max(i) + grid%margin
      extents(0, i) = atoms_min(i)
      extents(1, i) = atoms_max(i)
      !write(*, *) 'grid_extents_from_mol : i, min, max ', i, extents(0, i), extents(1, i) 
    end do
  end subroutine grid_extents_from_mol


  subroutine grid_npoints_from_spacing(grid, extents)
    type(grid_t), intent(inout):: grid
    real(8), intent(in):: extents(0:1, 3)
    integer:: i

    do i = 1, grid%ndim
      if (grid%nr_points_known(i)) return
      if (grid%spacing(i) /= 0) then
        grid%nr_points(i) = int((extents(1, i) - extents(0, i))/grid%spacing(i)) + 1
        if (grid%nr_points(i) <= 0) then
          ! assign 1 point to this dimension;
          ! we do this because we work in grid%ndim dimensions
          ! and so we cannot have 0 points in one of the dimensions
          grid%nr_points(i) = 1
        end if
      end if
      grid%nr_points_known(i) = .true.
    end do
  end subroutine grid_npoints_from_spacing


  subroutine grid_spacing_from_npoints(grid, extents)
    type(grid_t), intent(inout):: grid
    real(8), intent(in):: extents(0:1, 3)
    integer:: i

    do i = 1, grid%ndim
      if (grid%spacing_known(i)) return
      grid%spacing(i) = 0.0d0
      if (grid%nr_points(i) > 1) then
        grid%spacing(i) = (extents(1, i) - extents(0, i))/(grid%nr_points(i) - 1)
      else
        grid%spacing(i) = extents(1, i) - extents(0, i)
      end if
      grid%spacing_known(i) = .true.
    end do
  end subroutine grid_spacing_from_npoints

  subroutine grid_origin_from_mol(grid, mol, extents)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    real(8), intent(inout), optional:: extents(0:1, 3)

    if (.not. present(extents)) call grid_extents_from_mol(mol, grid, extents)
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
  end subroutine grid_origin_from_mol


  subroutine center_mol_in_grid(grid)
#include "orgcom.h"
    type(grid_t), intent(inout):: grid
    real(8):: extents(0:1, 3)
    real(8):: cm_0, cm_1
    integer:: i

    do i = 1, 3
      ! get actual grid extents
      extents(0, i) = grid%origin(i)
      extents(1, i) = grid%origin(i) + (grid%nr_points(i) - 1)*grid%spacing(i)
      cm_0 = 0.5d0*(extents(1, i) + extents(0, i))
      cm_1 = cmxyz(i)  ! this is center-of-mass from "orgcom.h" in a.u.
      extents(0, i) = extents(0, i) + (cm_1 - cm_0)
      extents(1, i) = extents(1, i) + (cm_1 - cm_0)
    end do
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
  end subroutine center_mol_in_grid


  subroutine grid_get_nr_points_total(grid)
    type(grid_t), intent(inout):: grid
    integer:: i

    if (grid%nr_points_total_known) return
    grid%nr_points_total = 1
    do i = 1, grid%ndim
      grid%nr_points_total = grid%nr_points_total * grid%nr_points(i)
    end do
    grid%nr_points_total_known = .true.
  end subroutine grid_get_nr_points_total


  subroutine get_grid_essentials(grid, mol)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    real(8):: origin(3)
    real(8):: extents(0:1, 3)
    integer:: i

    call grid_extents_from_mol(mol, grid, extents)
    if (.not. grid%origin_known) call grid_origin_from_mol(grid, mol, extents)

    do i = 1, grid%ndim
      if (.not. grid%spacing_known(i)) call grid_spacing_from_npoints(grid, extents)
      if (.not. grid%nr_points_known(i)) call grid_npoints_from_spacing(grid, extents)
    end do

    if (.not. grid%nr_points_total_known) call grid_get_nr_points_total(grid)
  end subroutine get_grid_essentials


  subroutine get_grid_info(grid_id, grid, mol)
    integer, intent(in):: grid_id
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid

    if (grid%input  ==  1) then
      ! grid data to read from file
      call get_grid_info_from_file(grid)
      if (grid%id == 0) grid%id = grid_id
    else if (grid%input  ==  0) then
      ! grid to generate
      if (grid_id /= grid%id) then
        print *, 'ERROR: problem with grid id in get_grid_info'
        return
      end if
      call get_grid_essentials(grid, mol)
      ! center grid around molecule,
      ! modify grid extents and grid%origin (but do not move molecule)
      call center_mol_in_grid(grid)
    end if

    ! we should now have all information to work on the grid
    grid%status = 0

  end subroutine get_grid_info


  subroutine get_grid_points(grid)
    type(grid_t), intent(inout):: grid

    if (.not. grid%allocated_points()) return

    if (grid%input == 1) then
      call get_grid_points_from_file(grid)
      !write(*, *) 'Grid points for grid id = ', grid%id, ' are read from file ', grid%file_inp
    else if (grid%input == 0) then
      call generate_grid_points(grid)
      !write(*, *) 'Grid points for grid id = ', grid%id, ' are generated'
    end if
  end subroutine


  subroutine generate_grid_points(grid)
    type(grid_t), intent(inout):: grid
    integer :: i

    select case (grid%typ)
      case (1)
        call generate_grid_rectilinear(grid)
      case default
        call random_number(grid%points)
    end select
  end subroutine

end module
