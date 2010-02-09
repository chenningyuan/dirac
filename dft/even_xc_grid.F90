module even_xc_grid

  use memory_allocator

  implicit none

  public init_even_xc_grid
  public set_even_xc_grid

  integer, public :: even_xc_grid_nr_points

  logical, public :: plot_even_xc_grid    = .false.
  logical, public :: plot_regular_xc_grid = .false.

  private

  real(8) :: d(3), l(3)
  integer :: n(3)

contains

  subroutine set_even_xc_grid(nr_points, x, y, z, w)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: nr_points
    real(8), intent(inout) :: x(nr_points)
    real(8), intent(inout) :: y(nr_points)
    real(8), intent(inout) :: z(nr_points)
    real(8), intent(inout) :: w(nr_points)
!   ----------------------------------------------------------------------------
    integer                :: ipoint, ix, iy, iz
!   ----------------------------------------------------------------------------

    w = 1.0d0

    ipoint = 0
    do ix = 1, 2*n(1) + 1
      do iy = 1, 2*n(2) + 1
        do iz = 1, 2*n(3) + 1
          ipoint = ipoint + 1

          x(ipoint) = -l(1) + (ix - 1)*d(1)
          y(ipoint) = -l(2) + (iy - 1)*d(2)
          z(ipoint) = -l(3) + (iz - 1)*d(3)
        end do
      end do
    end do

  end subroutine

  subroutine init_even_xc_grid(nx, ny, nz, lx, ly, lz)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    real(8), intent(in) :: lx
    real(8), intent(in) :: ly
    real(8), intent(in) :: lz
!   ----------------------------------------------------------------------------
    integer             :: ixyz
!   ----------------------------------------------------------------------------

    n(1) = nx
    n(2) = ny
    n(3) = nz

    l(1) = lx
    l(2) = ly
    l(3) = lz

    do ixyz = 1, 3
      if (n(ixyz) > 0) then
        d(ixyz) = l(ixyz)/n(ixyz)
      else
        d(ixyz) = 0.0d0
        l(ixyz) = 0.0d0
      end if
    end do

    even_xc_grid_nr_points = (2*n(1) + 1)*(2*n(2) + 1)*(2*n(3) + 1)

  end subroutine

end module
