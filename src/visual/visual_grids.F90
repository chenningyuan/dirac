module visual_grids

  use visual_datatypes, only: grid_t, mol_t
  use visual_io
  use visual_utils

  implicit none

  public get_grid_info
  public get_grid_points

  interface generate_grid_rectilinear
    module procedure generate_grid_rectilinear
  end interface

  private

contains

! ------------------------------------------------------------------------------
  subroutine get_grid_info(grid_id, grid, mol)
    integer, intent(in):: grid_id
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    character(len = 150):: errmsg

    if (grid%input  ==  1) then
!
!   read grid information from file
!
      call get_grid_info_from_file(grid)
      if (grid%id == 0) grid%id = grid_id

!
!   in other cases - set up information to generate grid points
!
    else if (grid%input==0 .or. grid%input==5) then
!     3d grid
      if (grid_id /= grid%id) then
        write (errmsg, *) 'ERROR: problem with grid id in get_grid_info: ',&
                          ' grid%id=', grid%id,&
                          ' grid_id=', grid_id
        call visual_quit(errmsg)
      end if
      if (grid%typ == 1) then
        ! rectilinear uniform 3d grid
        call get_3d_grid_essentials(grid, mol)
        ! center grid around molecule,
        ! modify grid extents and grid%origin (but do not move molecule)
        ! todo - maybe add a keyword + an option to control grid extents
        call center_mol_in_3d_grid(grid)
      else if (grid%typ == 2) then
        ! radial 3d grid
        call get_3d_radial_grid_essentials(grid, mol)
      end if

    else if (grid%input  ==  3) then
      ! line
      call get_line_grid_essentials(grid)

    end if

    ! we should now have all information to work on the grid
    grid%status = 0

  end subroutine get_grid_info


! ------------------------------------------------------------------------------
  subroutine get_grid_points(grid)
    type(grid_t), intent(inout):: grid

    if (.not. grid%allocated_points()) return

    if (grid%input == 1) then
      call get_grid_points_from_file(grid)
    else if (grid%input /= 2) then
      ! grid%input = 2 is for .LIST;
      ! in this case, the points should already be read from input
      call generate_grid_points(grid)
    end if
  end subroutine



! ------------------------------------------------------------------------------
  subroutine generate_grid_radial(grid)
use num_grid_cfg
    type(grid_t), intent(inout):: grid
    real(8),allocatable:: temp(:,:)
    real(8):: rad, radial_spacing
    integer:: nleb, i, m, n, istep, ipoint
    character(len = 150):: errmsg

    if (grid%typ /= 2) return

    if (grid%input == 0 .or. grid%input == 5) then

      ! Get angular points and weights
      m = NLEB(num_grid_cfg_angint)
      write(*, '(/,4x, a)')  '** Radial distribution **'
      write(*, '(6x, A, I3)')  '   Lebedev angular grid L = ',num_grid_cfg_angint
      write(*, '(6x, A, I10)') '   Angular points:   ',m
      allocate(temp(m, 4))
      temp=0.0d0
      CALL LEBEDEV(num_grid_cfg_angint,temp(1,1),temp(1,2),temp(1,3),temp(1, 4),n)
      CALL DSCAL(m, 4*ACOS(-1.0D0), temp(1, 4), 1)

      ! Loop over radial points
      if (grid%radial_nr_steps > 1) then
        radial_spacing = grid%radial_length/(grid%radial_nr_steps - 1)
      else
        radial_spacing = grid%radial_length
      end if

      ipoint = 1
      do istep = 1, grid%radial_nr_steps
        rad = (istep-1)*radial_spacing
        do i = 1, m
          grid%points(1:3, ipoint) = grid%origin(1:3) + rad*temp(i, 1:3)
          ipoint = ipoint + 1
        end do
      end do
      deallocate(temp)

    else
      write (errmsg, *) 'ERROR: generation of radial grids other than ',&
                        '3d is not implemented'
      call visual_quit(errmsg)
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine generate_grid_rectilinear(grid)
! rectilinear equidistant
    type(grid_t), intent(inout):: grid
    integer:: i,j,k, ipoint, N(grid%ndim)
    real(8):: v(3), vx(3), vy(3), vz(3)
    real(8):: vx_norm, vy_norm, vz_norm
    character(len = 150):: errmsg

    if (grid%typ /= 1) return

    if (grid%input == 0 .or. grid%input == 5) then
!     3d
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
    else if (grid%input == 4) then
!     2d
      ipoint = 1
      vx(1:3) = (grid%extents(1,1:3) - grid%extents(0,1:3))/(grid%nr_points(1)-1)
      vy(1:3) = (grid%extents(2,1:3) - grid%extents(0,1:3))/(grid%nr_points(2)-1)
      vz(1) = vx(2)*vy(3) - vx(3)*vy(2)
      vz(2) = vx(3)*vy(1) - vx(1)*vy(3)
      vz(3) = vx(1)*vy(2) - vx(2)*vy(1)
      vx_norm = dsqrt(vx(1)*vx(1) + vx(2)*vx(2) + vx(3)*vx(3))
      vy_norm = dsqrt(vy(1)*vy(1) + vy(2)*vy(2) + vy(3)*vy(3))
      vz_norm = dsqrt(vz(1)*vz(1) + vz(2)*vz(2) + vz(3)*vz(3))
      ipoint = 1
      do i = 1, grid%nr_points(1)
        do j = 1, grid%nr_points(2)
          grid%points(1:3, ipoint) = grid%origin(1:3) + (i-1)*vx(1:3) + (j-1)*vy
          ipoint = ipoint + 1
        end do
      end do
    else if (grid%input == 3) then
!     line
      ipoint = 1
      do i = 1, grid%nr_points(1)
        v(1:3) = grid%origin(1:3) + (i - 1)*grid%spacing(1:3)
        grid%points(1:3, ipoint) = v(1:3)
        ipoint = ipoint + 1
      end do
    else
      write (errmsg, *) 'ERROR: unsupported grid input '
      call visual_quit(errmsg)
    end if

  end subroutine

! ------------------------------------------------------------------------------
  subroutine grid_file_txt_write(file_unit, float_format, nr_points, x, y, z)
     integer, intent(in):: file_unit
     character(*), intent(in):: float_format
     integer, intent(in):: nr_points
     real(8), intent(in):: x(:), y(:), z(:)
     character(len = 8):: xc, yc, zc
     integer             :: i

     do i = 1, nr_points
        write(file_unit, float_format) x(i), y(i), z(i)
     enddo
  end subroutine


! ------------------------------------------------------------------------------
  subroutine grid_file_txt_read(file_unit, float_format, nr_points, x, y, z)
     real(8), intent(inout), dimension(:):: x, y, z
     integer, intent(in)  :: file_unit
     character(*), intent(in):: float_format
     integer, intent(in)  :: nr_points
     integer              :: i

     rewind(file_unit)

     i = 1
     do
        read(file_unit, *, end = 98) x(i), y(i), z(i)
        i = i+1
     end do
 98  continue

  end subroutine


! ------------------------------------------------------------------------------
  subroutine grid_origin_extents_from_mol(mol, grid)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    real(8):: extents(0:1, 3)
    real(8):: atoms_min(1:3), atoms_max(1:3)
    integer::i, j

    atoms_min(1:3) = mol%geometry(1:3, 1)
    atoms_max(1:3) = mol%geometry(1:3, 1)

    do i = 1, mol%n_atoms
      do j = 1, 3
        atoms_min(j) = min(mol%geometry(j, i), atoms_min(j))
        atoms_max(j) = max(mol%geometry(j, i), atoms_max(j))
      end do
    end do

    ! add space around
    do i = 1, 3
      atoms_min(i) = atoms_min(i) - grid%margin
      atoms_max(i) = atoms_max(i) + grid%margin
      extents(0, i) = atoms_min(i)
      extents(1, i) = atoms_max(i)
    end do

    ! update grid info
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
    grid%extents(0, 1:3) = extents(0, 1:3)
    grid%extents(1, 1:3) = extents(1, 1:3)
    grid%extents_known = .true.

  end subroutine grid_origin_extents_from_mol


! ------------------------------------------------------------------------------
  subroutine grid_3d_npoints_from_spacing(grid)
    type(grid_t), intent(inout):: grid
    integer:: i

    if (.not.grid%extents_known) return
    do i = 1, grid%ndim
      if (grid%nr_points_known(i) .or. &
          grid%spacing(i) == 0    .or. &
          .not.grid%spacing_known(i) .or. &
          grid%nr_points(i) <= 0) return

      grid%nr_points(i) = int((grid%extents(1,i) - grid%extents(0,i))/grid%spacing(i)) + 1
      grid%nr_points_known(i) = .true.
    end do
  end subroutine grid_3d_npoints_from_spacing


! ------------------------------------------------------------------------------
  subroutine grid_3d_spacing_from_npoints(grid)
    type(grid_t), intent(inout):: grid
    integer:: i
    character(len = 150):: errmsg

    do i = 1, grid%ndim
      if (grid%spacing_known(i)) return
      grid%spacing(i) = 0.0d0
      if (grid%nr_points_known(i)) then
        if (grid%nr_points(i) > 1) then
          grid%spacing(i) = (grid%extents(1,i) - grid%extents(0, i))/(grid%nr_points(i) - 1)
        else
          grid%spacing(i) = grid%extents(1,i) - grid%extents(0,i)
        end if
        grid%spacing_known(i) = .true.
      else
        write (errmsg, *) 'ERROR in grid_3d_spacing_from_npoints'
        call visual_quit(errmsg)
      end if
    end do
  end subroutine grid_3d_spacing_from_npoints

! ------------------------------------------------------------------------------
  subroutine grid_line_spacing_from_npoints(grid)
    type(grid_t), intent(inout):: grid
    integer:: i
    character(len = 150):: errmsg

    do i = 1, grid%ndim
      if (grid%spacing_known(i)) return
      grid%spacing(i) = 0.0d0
      if (grid%nr_points_total_known) then
        if (grid%nr_points_total > 1) then
          grid%spacing(i) = (grid%extents(1,i) - grid%extents(0, i))/(grid%nr_points_total - 1)
        else
          grid%spacing(i) = grid%extents(1,i) - grid%extents(0,i)
        end if
        grid%spacing_known(i) = .true.
      else
        write (errmsg, *) 'ERROR in grid_3d_spacing_from_npoints'
        call visual_quit(errmsg)
      end if
    end do
  end subroutine grid_line_spacing_from_npoints


! ------------------------------------------------------------------------------
  subroutine center_mol_in_3d_grid(grid)
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
      !cm_1 = cmxyz(i)  ! this is center-of-mass from "orgcom.h" in a.u.
      cm_1 = gagorg(i) ! this is a "magnetic gauge origin" (from "orgcom.h" in a.u.), was used in old DIRAC-VISUAL code 
      extents(0, i) = extents(0, i) + (cm_1 - cm_0)
      extents(1, i) = extents(1, i) + (cm_1 - cm_0)
    end do
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
  end subroutine center_mol_in_3d_grid


! ------------------------------------------------------------------------------
  subroutine grid_get_nr_points_total(grid)
    type(grid_t), intent(inout):: grid
    integer:: i

    if (grid%nr_points_total_known) return
    grid%nr_points_total = 1
    do i = 1, 3
      if (grid%nr_points(i) > 0) then
        grid%nr_points_total = grid%nr_points_total * grid%nr_points(i)
      end if
    end do
    grid%nr_points_total_known = .true.
  end subroutine grid_get_nr_points_total


! ------------------------------------------------------------------------------
  subroutine get_3d_grid_essentials(grid, mol)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    integer:: i

    if (grid%typ /= 1) return

    call grid_origin_extents_from_mol(mol, grid)

    do i = 1, 3
      if (.not. grid%spacing_known(i)) call grid_3d_spacing_from_npoints(grid)
      if (.not. grid%nr_points_known(i)) call grid_3d_npoints_from_spacing(grid)
    end do

    if (.not. grid%nr_points_total_known) call grid_get_nr_points_total(grid)
  end subroutine get_3d_grid_essentials


! ------------------------------------------------------------------------------
  subroutine get_3d_radial_grid_essentials(grid, mol)
use num_grid_cfg
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    integer:: m, nleb

    if (grid%typ /= 2) return

    m = NLEB(num_grid_cfg_angint)
    grid%nr_points_total = m*grid%radial_nr_steps
    grid%nr_points_total_known = .true.

  end subroutine get_3d_radial_grid_essentials


! ------------------------------------------------------------------------------
  subroutine get_line_grid_essentials(grid)
    type(grid_t), intent(inout):: grid
    real(8):: origin(3)
    real(8):: extents(0:1, 3)
    integer:: i

    if (grid%typ /= 1) return
    call grid_line_spacing_from_npoints(grid)

  end subroutine get_line_grid_essentials


! ------------------------------------------------------------------------------
  subroutine generate_grid_points(grid)
    type(grid_t), intent(inout):: grid
    integer :: i
    character(len = 150):: errmsg

    select case (grid%typ)
      case (1)
        call generate_grid_rectilinear(grid)
      case (2)
        call generate_grid_radial(grid)
      case default
        write (errmsg, *) 'ERROR: generation of grids other than ',&
                          '3d-rectilinear or radial is not implemented'
        call visual_quit(errmsg)
    end select
  end subroutine

end module
