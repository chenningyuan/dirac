
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
    real(kind_dp) :: x
    real(kind_dp) :: y
    real(kind_dp) :: z
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

    v%x = a
    v%y = b
    v%z = c
  end subroutine set_dp_vector

  function dp_vector_plus_dp_vector(u, v) result(w)
    type(dp_vector_t), intent(in) :: u, v
    type(dp_vector_t)             :: w

    w%x = u%x + v%x
    w%y = u%y + v%y
    w%z = u%z + v%z
  end function dp_vector_plus_dp_vector

  function dp_vector_minus_dp_vector(u, v) result(w)
    type(dp_vector_t), intent(in) :: u, v
    type(dp_vector_t)             :: w

    w%x = u%x - v%x
    w%y = u%y - v%y
    w%z = u%z - v%z
  end function dp_vector_minus_dp_vector

  function dp_scalar_times_dp_vector(s, v) result(u)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp),     intent(in) :: s
    type(dp_vector_t)             :: u

    u%x = s*v%x
    u%y = s*v%y
    u%z = s*v%z
  end function dp_scalar_times_dp_vector

  function i4_scalar_times_dp_vector(s, v) result(u)
    type(dp_vector_t), intent(in) :: v
    integer(kind_i4),  intent(in) :: s
    type(dp_vector_t)             :: u

    u%x = s*v%x
    u%y = s*v%y
    u%z = s*v%z
  end function i4_scalar_times_dp_vector

  function dp_vector_times_dp_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp),     intent(in) :: s
    type(dp_vector_t)             :: u

    u%x = s*v%x
    u%y = s*v%y
    u%z = s*v%z
  end function dp_vector_times_dp_scalar

  function dp_vector_times_i4_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    integer(kind_i4),  intent(in) :: s
    type(dp_vector_t)             :: u

    u%x = s*v%x
    u%y = s*v%y
    u%z = s*v%z
  end function dp_vector_times_i4_scalar

  function dp_vector_dot_dp_vector(u, v) result(s)
    type(dp_vector_t), intent(in) :: u, v
    real(kind_dp)                 :: s

    s = u%x*v%x + u%y*v%y + u%z*v%z
  end function dp_vector_dot_dp_vector

  function dp_vector_divided_by_dp_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp),     intent(in) :: s
    type(dp_vector_t)             :: u

    u%x = v%x/s
    u%y = v%y/s
    u%z = v%z/s
  end function dp_vector_divided_by_dp_scalar

  function dp_vector_divided_by_i4_scalar(v, s) result(u)
    type(dp_vector_t), intent(in) :: v
    integer(kind_i4),  intent(in) :: s
    type(dp_vector_t)             :: u

    u%x = v%x/s
    u%y = v%y/s
    u%z = v%z/s
  end function dp_vector_divided_by_i4_scalar

  function dp_vector_cross_dp_vector(u, v) result(w)
    type(dp_vector_t), intent(in) :: u, v
    type(dp_vector_t)             :: w

    w%x = u%y*v%z - u%z*v%y
    w%y = u%z*v%x - u%x*v%z
    w%z = u%x*v%y - u%y*v%x
  end function dp_vector_cross_dp_vector

  subroutine copy_dp_vector_to_dp_vector(u, v)
    type(dp_vector_t), intent(out) :: u
    type(dp_vector_t), intent(in)  :: v

    u%x = v%x
    u%y = v%y
    u%z = v%z
  end subroutine copy_dp_vector_to_dp_vector

  subroutine copy_dp_array_to_dp_vector(v, m)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: m(3)

    v%x = m(1)
    v%y = m(2)
    v%z = m(3)
  end subroutine copy_dp_array_to_dp_vector

  subroutine copy_dp_scalar_to_dp_vector_components(v, s)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: s

    v%x = s
    v%y = s
    v%z = s
  end subroutine copy_dp_scalar_to_dp_vector_components

  subroutine normalize_dp_vector_to(v, norm)
    type(dp_vector_t), intent(out) :: v
    real(kind_dp),     intent(in)  :: norm

    real(kind_dp)                  :: actual_norm

    actual_norm = (v%x**2.0d0 + v%y**2.0d0 + v%z**2.0d0)**0.5d0

    v%x = norm*v%x/actual_norm
    v%y = norm*v%y/actual_norm
    v%z = norm*v%z/actual_norm
  end subroutine normalize_dp_vector_to

  function dp_vector_norm(v) result(norm)
    type(dp_vector_t), intent(in) :: v
    real(kind_dp)                 :: norm

    norm = (v%x**2.0d0 + v%y**2.0d0 + v%z**2.0d0)**0.5d0
  end function dp_vector_norm

end module vector_algebra_m
