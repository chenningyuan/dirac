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

  character(len = 16):: float_format_txt = '(*(es18.8e3, :))'

  private

contains

! ------------------------------------------------------------------------------
  subroutine generate_grid_rectilinear(grid)
    type(grid_t), intent(inout):: grid
    integer:: i,j,k, ipoint
    real(8):: v(3)

    if (grid%ndim /= 3) return
    write(*, *) 'HERENOw!'

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
  subroutine grid_extents_from_mol(mol, grid, extents)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(in):: grid
    real(8), intent(out):: extents(0:1, 3)
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
  end subroutine grid_extents_from_mol


! ------------------------------------------------------------------------------
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


! ------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------
  subroutine grid_origin_from_mol(grid, mol, extents)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    real(8), intent(inout), optional:: extents(0:1, 3)

    if (.not. present(extents)) call grid_extents_from_mol(mol, grid, extents)
    grid%origin(1:3) = extents(0, 1:3)
    grid%origin_known = .true.
  end subroutine grid_origin_from_mol


! ------------------------------------------------------------------------------
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


! ------------------------------------------------------------------------------
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


! ------------------------------------------------------------------------------
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


! ------------------------------------------------------------------------------
  subroutine get_grid_info(grid_id, grid, mol)
    integer, intent(in):: grid_id
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid

    !if (grid%input  ==  0) then ! .or. grid%input  ==  5) then
    if (grid%input==0 .or. grid%input==5) then
      ! generate 3d grid
      if (grid_id /= grid%id) then
        call visual_quit('ERROR: problem with grid id in get_grid_info')
        print *, 'ERROR: problem with grid id in get_grid_info'
        print *, ' * grid%id=', grid%id
        print *, ' * grid_id=', grid_id
        return
      end if
      call get_grid_essentials(grid, mol)
      ! center grid around molecule,
      ! modify grid extents and grid%origin (but do not move molecule)
      call center_mol_in_grid(grid)
    else if (grid%input  ==  1) then
      ! read grid data from file
      call get_grid_info_from_file(grid)
      if (grid%id == 0) grid%id = grid_id
    else if (grid%input  ==  2) then
      ! read grid data from user input: look for '.LINE'
      !call get_grid_info_from_txt_storage(grid, 'line')
    else if (grid%input  ==  3) then
    else if (grid%input  ==  4) then
    !else if (grid%input  ==  5) then
    else if (grid%input  ==  6) then
    end if

    ! we should now have all information to work on the grid
    grid%status = 0

  end subroutine get_grid_info


! ------------------------------------------------------------------------------
  subroutine get_grid_points(grid)
    type(grid_t), intent(inout):: grid

    if (.not. grid%allocated_points()) return

    !if (grid%input == 0) then ! .or. grid%input == 5) then
    if (grid%input== 0 .or. grid%input==5) then
      call generate_grid_points(grid)
    else if (grid%input == 1) then
      call get_grid_points_from_file(grid)
    end if
  end subroutine


! ------------------------------------------------------------------------------
  subroutine generate_grid_points(grid)
    type(grid_t), intent(inout):: grid
    integer :: i

    select case (grid%typ)
      case (1)
        call generate_grid_rectilinear(grid)
      case default
        stop
        !FIXME
        !call random_number(grid%points)
    end select
  end subroutine

end module
