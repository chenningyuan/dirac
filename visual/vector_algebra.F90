
!-------------------------------------------------------------------------------
!      collection of vector algebra operations
!      assumes dp scalars and vector components
!
!      radovan bast (2008-04-06)
!      inspiration: the nice code example of A. Bay and T. Monteiro
!                   http://www.nag.co.uk/nagware/Examples/Vector.f90
!-------------------------------------------------------------------------------

module vector_algebra_m

  use kinds_m

  implicit none

  type vector_t
    real(kind_dp) :: x
    real(kind_dp) :: y
    real(kind_dp) :: z
  end type vector_t

  interface operator(+)
    module procedure vector_plus_vector
  end interface

  interface operator(-)
    module procedure vector_minus_vector
  end interface

  interface operator(*)
    module procedure scalar_times_vector, vector_times_scalar, vector_dot_vector
  end interface

  interface operator(.cross.)
    module procedure vector_cross_vector
  end interface

  interface assignment(=)
    module procedure copy_vector_to_vector, copy_array_to_vector, copy_scalar_to_vector_components
  end interface

contains

  subroutine set_vector(v, a, b, c)
    type(vector_t), intent(out) :: v
    real(kind_dp),  intent(in)  :: a, b, c

    v%x = a
    v%y = b
    v%z = c
  end subroutine set_vector

  function vector_plus_vector(u, v)
    type(vector_t), intent(in) :: u, v
    type(vector_t) vector_plus_vector

    vector_plus_vector%x = u%x + v%x
    vector_plus_vector%y = u%y + v%y
    vector_plus_vector%z = u%z + v%z
  end function vector_plus_vector

  function vector_minus_vector(u, v)
    type(vector_t), intent(in) :: u, v
    type(vector_t) vector_minus_vector

    vector_minus_vector%x = u%x - v%x
    vector_minus_vector%y = u%y - v%y
    vector_minus_vector%z = u%z - v%z
  end function vector_minus_vector

  function scalar_times_vector(s, v)
    real(kind_dp),  intent(in) :: s
    type(vector_t), intent(in) :: v
    type(vector_t) scalar_times_vector

    scalar_times_vector%x = s*v%x
    scalar_times_vector%y = s*v%y
    scalar_times_vector%z = s*v%z
  end function scalar_times_vector

  function vector_times_scalar(v, s)
    real(kind_dp),  intent(in) :: s
    type(vector_t), intent(in) :: v
    type(vector_t) vector_times_scalar

    vector_times_scalar%x = s*v%x
    vector_times_scalar%y = s*v%y
    vector_times_scalar%z = s*v%z
  end function vector_times_scalar

  function vector_dot_vector(u, v)
    type(vector_t), intent(in) :: u, v
    real(kind_dp) vector_dot_vector

    vector_dot_vector = u%x*v%x + u%y*v%y + u%z*v%z
  end function vector_dot_vector

  function vector_cross_vector(u, v)
    type(vector_t), intent(in) :: u, v
    type(vector_t) vector_cross_vector

    vector_cross_vector%x = u%y*v%z - u%z*v%y
    vector_cross_vector%y = u%z*v%x - u%x*v%z
    vector_cross_vector%z = u%x*v%y - u%y*v%x
  end function vector_cross_vector

  subroutine copy_vector_to_vector(u, v)
    type(vector_t), intent(out) :: u
    type(vector_t), intent(in)  :: v

    u%x = v%x
    u%y = v%y
    u%z = v%z
  end subroutine copy_vector_to_vector

  subroutine copy_array_to_vector(v, m)
    type(vector_t), intent(out) :: v
    real(kind_dp),  intent(in)  :: m(3)

    v%x = m(1)
    v%y = m(2)
    v%z = m(3)
  end subroutine copy_array_to_vector

  subroutine copy_scalar_to_vector_components(v, s)
    type(vector_t), intent(out) :: v
    real(kind_dp),  intent(in)  :: s

    v%x = s
    v%y = s
    v%z = s
  end subroutine copy_scalar_to_vector_components

end module vector_algebra_m
