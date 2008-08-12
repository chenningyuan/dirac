module svwn_m

  use kinds_m

  implicit none

  public slater_recursive

  private

  real(kind_dp), parameter :: small = 1.0d-20
  real(kind_dp), parameter :: pi    = 3.141592653589793d0
  real(kind_dp), parameter :: f_13  = 1.0d0/3.0d0
  real(kind_dp), parameter :: g     = -(6.0d0/pi)**f_13

contains

  function slater_recursive(order, n, w) result(v)

!          closed-shell arbitrary-order Slater functional derivative
!          ============
!
!          w \frac{d^i}{dn^i} F(n) = slater_recursive(i, n, w)
!
!          note that for this functional:
!
!          \frac{d^2}{dn^2} F(n, s) = F_nn    = F_ss
!          \frac{d^3}{dn^3} F(n, s) = F_nnn   = F_nss
!          \frac{d^4}{dn^4} F(n, s) = F_nnnn  = F_nnss = F_ssss
!          \frac{d^5}{dn^5} F(n, s) = F_nnnnn = ...
!          ...
!
!          all derivatives with unpair number of "s" are zero

    integer,       intent(in) :: order
    real(kind_dp), intent(in) :: n, w
    real(kind_dp)             :: v
!   ----------------------------------------------------------------------------
    real(kind_dp)             :: na, p, e, f
    integer                   :: i

    na = n/2.0d0

    if (abs(na) > small) then
      p = 1.0d0
      e = f_13

      if (order > 1) then
        do i = 2, order

           p = p*e
           e = e - 1.0d0

        end do
      end if

      f = p*g*na**e
      v = w*f/2.0d0**(order-1)
    else
      v = 0.0d0
    end if

  end function

end module
