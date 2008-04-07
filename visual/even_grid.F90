
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

  type even_grid_t
    type(dp_vector_t)   :: origin, u, v, w
    integer(kind_i4) :: nr_steps(3)
  end type even_grid_t

  public even_grid_t
  public init_even_grid, del_even_grid

  private

! blub

contains

  subroutine init_even_grid(self)
    type(even_grid_t) :: self

    type(dp_vector_t) :: a, b, c, ab, ac, avg
    real(kind_dp)     :: ab_along_u, ac_along_u, ac_along_v, l, f
    integer(kind_i4)  :: nr_steps = 50

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

    l = max(ac_along_u, ac_along_v) !change to l = max(ab_along_u, ac_along_v)

    f = l/nr_steps

    call normalize_dp_vector_to(self%u, f)
    call normalize_dp_vector_to(self%v, f)

    avg = (a + b + c)/3.0d0
 
    self%origin = avg - nr_steps*self%u - nr_steps*self%v

    write(*,*) 'origin%x ', self%origin%x
    write(*,*) 'origin%y ', self%origin%y
    write(*,*) 'origin%z ', self%origin%z



!for i in range(2*nr_steps + 1):
!    for j in range(2*nr_steps + 1):
!
!        p[0] = bottom_left[0] + i*u[0] + j*v[0]
!        p[1] = bottom_left[1] + i*u[1] + j*v[1]
!        p[2] = bottom_left[2] + i*u[2] + j*v[2]
!
!        scalar = fake_density(a, b, c, p)
!
!        quasi_x = i*scale - l
!        quasi_y = j*scale - l
!
!        j = fake_current(a, p)
!
!        uector_along_u = j[0]*u[0] + j[1]*u[1] + j[2]*u[2]
!        uector_along_v = j[0]*v[0] + j[1]*v[1] + j[2]*v[2]
!
!        print "%14.6e%14.6e%14.6e%14.6e%14.6e" % (quasi_x, quasi_y, scalar, uector_along_u, uector_along_v)








  end subroutine init_even_grid

  subroutine del_even_grid(self)
    type(even_grid_t) :: self

  end subroutine del_even_grid

end module even_grid_class
