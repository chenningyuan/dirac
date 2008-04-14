
!-------------------------------------------------------------------------------
!      collection of vector algebra operations
!
!      radovan bast (2008-04-06)
!      inspiration: the nice code example of A. Bay and T. Monteiro
!                   http://www.nag.co.uk/nagware/Examples/Vector.f90
!-------------------------------------------------------------------------------

module vector_algebra_m

  use kinds_m

  implicit none

  type dp_vector_t
    real(kind_dp) :: component(1:3)
  end type dp_vector_t

  interface operator(+)
    module procedure dp_vector_plus_dp_vector
  end interface

  interface operator(-)
    module procedure dp_vector_minus_dp_vector
  end interface

  interface operator(*)
    module procedure dp_scalar_times_dp_vector
    module procedure i4_scalar_times_dp_vector
    module procedure dp_vector_times_dp_scalar
    module procedure dp_vector_times_i4_scalar
    module procedure dp_vector_dot_dp_vector
  end interface

  interface operator(/)
    module procedure dp_vector_divided_by_dp_scalar
    module procedure dp_vector_divided_by_i4_scalar
  end interface

  interface operator(.x.)
    module procedure dp_vector_cross_dp_vector
  end interface

  interface assignment(=)
    module procedure copy_dp_vector_to_dp_vector
    module procedure copy_dp_array_to_dp_vector
    module procedure copy_dp_scalar_to_dp_vector_components
  end interface

contains

  subroutine set_dp_vector(v, a, b, c)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: a, b, c

    v%component(1) = a
    v%component(2) = b
    v%component(3) = c
  end subroutine set_dp_vector

  function dp_vector_plus_dp_vector(u, v) result(w)
    type(dp_vector_t), intent(in) :: u, v
    type(dp_vector_t)             :: w

    w%component(1) = u%component(1) + v%component(1)
    w%component(2) = u%component(2) + v%component(2)
    w%component(3) = u%component(3) + v%component(3)
  end function dp_vector_plus_dp_vector

  function dp_vector_minus_dp_vector(u, v) result(w)
    type(dp_vector_t), intent(in) :: u, v
    type(dp_vector_t)             :: w

    w%component(1) = u%component(1) - v%component(1)
    w%component(2) = u%component(2) - v%component(2)
    w%component(3) = u%component(3) - v%component(3)
  end function dp_vector_minus_dp_vector

  function dp_scalar_times_dp_vector(s, v) result(u)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp),     intent(in) :: s
    type(dp_vector_t)             :: u

    u%component(1) = s*v%component(1)
    u%component(2) = s*v%component(2)
    u%component(3) = s*v%component(3)
  end function dp_scalar_times_dp_vector

  function i4_scalar_times_dp_vector(s, v) result(u)
    type(dp_vector_t), intent(in) :: v
    integer(kind_i4),  intent(in) :: s
    type(dp_vector_t)             :: u

    u%component(1) = s*v%component(1)
    u%component(2) = s*v%component(2)
    u%component(3) = s*v%component(3)
  end function i4_scalar_times_dp_vector

  function dp_vector_times_dp_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp),     intent(in) :: s
    type(dp_vector_t)             :: u

    u%component(1) = s*v%component(1)
    u%component(2) = s*v%component(2)
    u%component(3) = s*v%component(3)
  end function dp_vector_times_dp_scalar

  function dp_vector_times_i4_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    integer(kind_i4),  intent(in) :: s
    type(dp_vector_t)             :: u

    u%component(1) = s*v%component(1)
    u%component(2) = s*v%component(2)
    u%component(3) = s*v%component(3)
  end function dp_vector_times_i4_scalar

  function dp_vector_dot_dp_vector(u, v) result(s)
    type(dp_vector_t), intent(in) :: u, v
    real(kind_dp)                 :: s

    s = u%component(1)*v%component(1) &
      + u%component(2)*v%component(2) &
      + u%component(3)*v%component(3)
  end function dp_vector_dot_dp_vector

  function dp_vector_divided_by_dp_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp),     intent(in) :: s
    type(dp_vector_t)             :: u

    u%component(1) = v%component(1)/s
    u%component(2) = v%component(2)/s
    u%component(3) = v%component(3)/s
  end function dp_vector_divided_by_dp_scalar

  function dp_vector_divided_by_i4_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    integer(kind_i4),  intent(in) :: s
    type(dp_vector_t)             :: u

    u%component(1) = v%component(1)/s
    u%component(2) = v%component(2)/s
    u%component(3) = v%component(3)/s
  end function dp_vector_divided_by_i4_scalar

  function dp_vector_cross_dp_vector(u, v) result(w)
    type(dp_vector_t), intent(in) :: u, v
    type(dp_vector_t)             :: w

    w%component(1) = u%component(2)*v%component(3) - u%component(3)*v%component(2)
    w%component(2) = u%component(3)*v%component(1) - u%component(1)*v%component(3)
    w%component(3) = u%component(1)*v%component(2) - u%component(2)*v%component(1)
  end function dp_vector_cross_dp_vector

  subroutine copy_dp_vector_to_dp_vector(u, v)
    type(dp_vector_t), intent(out) :: u
    type(dp_vector_t), intent(in)  :: v

    u%component(1) = v%component(1)
    u%component(2) = v%component(2)
    u%component(3) = v%component(3)
  end subroutine copy_dp_vector_to_dp_vector

  subroutine copy_dp_array_to_dp_vector(v, m)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: m(3)

    v%component(1) = m(1)
    v%component(2) = m(1)
    v%component(3) = m(1)
  end subroutine copy_dp_array_to_dp_vector

  subroutine copy_dp_scalar_to_dp_vector_components(v, s)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: s

    v%component(1) = s
    v%component(2) = s
    v%component(3) = s
  end subroutine copy_dp_scalar_to_dp_vector_components

  subroutine normalize_dp_vector_to(v, norm)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: norm

    real(kind_dp)                  :: actual_norm

    actual_norm = (v%component(1)**2.0d0 + v%component(2)**2.0d0 + v%component(3)**2.0d0)**0.5d0

    v%component(1) = norm*v%component(1)/actual_norm
    v%component(2) = norm*v%component(2)/actual_norm
    v%component(3) = norm*v%component(3)/actual_norm
  end subroutine normalize_dp_vector_to

  function dp_vector_norm(v) result(norm)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp)                 :: norm

    norm = (v%component(1)**2.0d0 + v%component(2)**2.0d0 + v%component(3)**2.0d0)**0.5d0
  end function dp_vector_norm

end module vector_algebra_m
