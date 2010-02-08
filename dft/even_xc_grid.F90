module even_xc_grid

  use memory_allocator

  implicit none

  real(8), allocatable, target :: even_xc_grid_x(:)
  real(8), allocatable, target :: even_xc_grid_y(:)
  real(8), allocatable, target :: even_xc_grid_z(:)
  real(8), allocatable, target :: even_xc_grid_w(:)

  public init_even_xc_grid
  public kill_even_xc_grid

  private

contains

  subroutine init_even_xc_grid(nx, ny, nz, lx, ly, lz)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    real(8), intent(in) :: lx
    real(8), intent(in) :: ly
    real(8), intent(in) :: lz
!   ----------------------------------------------------------------------------
    integer             :: ipoint, nr_points
    integer             :: ix, iy, iz
    real(8)             :: dx, dy, dz
!   ----------------------------------------------------------------------------

    dx = lx/nx
    dy = ly/ny
    dz = lz/nz

    nr_points = (2*nx + 1)*(2*ny + 1)*(2*nz + 1)

    call alloc(even_xc_grid_x, nr_points)
    call alloc(even_xc_grid_y, nr_points)
    call alloc(even_xc_grid_z, nr_points)
    call alloc(even_xc_grid_w, nr_points)

    even_xc_grid_w = 1.0d0

    ipoint = 0

    do ix = 1, 2*nx + 1
      do iy = 1, 2*ny + 1
        do iz = 1, 2*nz + 1

          ipoint = ipoint + 1

          even_xc_grid_x(ipoint) = -lx + (ix - 1)*dx
          even_xc_grid_y(ipoint) = -ly + (iy - 1)*dy
          even_xc_grid_z(ipoint) = -lz + (iz - 1)*dz

        end do
      end do
    end do

  end subroutine

  subroutine kill_even_xc_grid()

!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------

    call dealloc(even_xc_grid_x)
    call dealloc(even_xc_grid_y)
    call dealloc(even_xc_grid_z)
    call dealloc(even_xc_grid_w)

  end subroutine

end module
