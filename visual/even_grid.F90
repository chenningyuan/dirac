
!-------------------------------------------------------------------------------
!      this is just for experimentation ...
!
!      radovan bast (2008-04-06)
!      inspiration: discussions with jonas juselius
!-------------------------------------------------------------------------------

module even_grid_class

  use kinds_m
  use vector_algebra_m

  implicit none

  type grid_data_t
    real(kind_dp), dimension(:), pointer :: x, y, z, w
  end type grid_data_t

  type even_grid_t
    type(dp_vector_t) :: origin, u, v, w
    integer(kind_i4)  :: nr_steps(3)
    type(grid_data_t) :: grid_data
  end type even_grid_t

  public grid_data_t, even_grid_t
  public init_even_grid, del_even_grid, plot_even_grid

  private

! blub

contains

  subroutine init_even_grid(self)
    type(even_grid_t), intent(inout) :: self

    type(dp_vector_t)                :: a, b, c, ab, ac, avg
    real(kind_dp)                    :: ab_along_u, ac_along_u, ac_along_v, l, f

    self%nr_steps(1) = 50
    self%nr_steps(2) = 50
    self%nr_steps(3) = 0

    call set_dp_vector(a, 0.0d0, 4.0d0, 1.0d0)
    call set_dp_vector(b, 1.0d0, 2.0d0, 0.0d0)
    call set_dp_vector(c, 0.0d0, 0.0d0, 5.0d0)

    ab = b - a
    ac = c - a

    self%u = ab
    self%w = ab .x. ac
    self%v = self%w .x. self%u

    call normalize_dp_vector_to(self%u, 1.0d0)
    call normalize_dp_vector_to(self%v, 1.0d0)

    ab_along_u = ab*self%u
    ac_along_u = ac*self%u
    ac_along_v = ac*self%v

    l = max(ab_along_u, ac_along_v)

    f = l/self%nr_steps(1)

    call normalize_dp_vector_to(self%u, f)
    call normalize_dp_vector_to(self%v, f)

    avg = (a + b + c)/3.0d0

    self%origin = avg - self%nr_steps(1)*self%u - self%nr_steps(2)*self%v
  end subroutine init_even_grid

  subroutine plot_even_grid(self)
    type(even_grid_t), intent(inout) :: self

    type(dp_vector_t)                :: p
    real(kind_dp)                    :: quasi_x, quasi_y
    real(kind_dp)                    :: vector_density_along_u, vector_density_along_v
    real(kind_dp)                    :: scalar_density
    type(dp_vector_t)                :: vector_density
    integer(kind_i4)                 :: i, j

    do i = 0, 2*self%nr_steps(1)
      do j = 0, 2*self%nr_steps(2)

         p = self%origin + i*self%u + j*self%v

         scalar_density = fake_scalar_density(p)
         vector_density = fake_vector_density(p)

         quasi_x = i*dp_vector_norm(self%u)
         quasi_y = j*dp_vector_norm(self%v)

         vector_density_along_u = vector_density*self%u
         vector_density_along_v = vector_density*self%v

         write(*,'(5e18.10)') quasi_x, quasi_y, scalar_density, vector_density_along_u, vector_density_along_v

      enddo
    enddo
  end subroutine plot_even_grid

  subroutine del_even_grid(self)
    type(even_grid_t) :: self

  end subroutine del_even_grid

  function fake_scalar_density(p) result(density)
    type(dp_vector_t), intent(in) :: p
    real(kind_dp)                 :: density

    type(dp_vector_t)             :: a, b, c, ap, bp, cp
    real(kind_dp)                 :: norm_ap, norm_bp, norm_cp

    call set_dp_vector(a, 0.0d0, 4.0d0, 1.0d0)
    call set_dp_vector(b, 1.0d0, 2.0d0, 0.0d0)
    call set_dp_vector(c, 0.0d0, 0.0d0, 5.0d0)

    ap = p - a
    bp = p - b
    cp = p - c

    density = 1.0d0/dp_vector_norm(ap) + 1.0d0/dp_vector_norm(bp) + 1.0d0/dp_vector_norm(cp)
  end function fake_scalar_density

  function fake_vector_density(p) result(density)
    type(dp_vector_t), intent(in) :: p
    type(dp_vector_t)             :: density

    type(dp_vector_t)             :: a, ap, n

    call set_dp_vector(a, 0.0d0, 4.0d0, 1.0d0)

    ap = p - a

    call set_dp_vector(n, 0.0d0, 0.0d0, 1.0d0)

    density = ap .x. n

    call normalize_dp_vector_to(density, 1.0d0)
  end function fake_vector_density

end module even_grid_class
