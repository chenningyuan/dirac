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

!          closed-shell arbitrary-order slater functional derivative
!
!          w*exc_n(n)    = slater_recursive(1, n, w)
!          w*exc_nn(n)   = slater_recursive(2, n, w)
!          w*exc_nnn(n)  = slater_recursive(3, n, w)
!          w*exc_nnnn(n) = ...
!
!          note that for this functional: exc_nn    = exc_ss
!                                         exc_nnn   = exc_nss
!                                         exc_nnnn  = exc_nnss = exc_ssss
!                                         exc_nnnnn = ...

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
