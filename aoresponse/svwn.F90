module svwn_m

  use kinds_m

  implicit none

  public slater_recursive
  public vwn_energy

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
!          (which is true for all LDA-type functionals at closed-shell reference)

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
      v = w*f/(2.0d0**(order-1))
    else
      v = 0.0d0
    end if

  end function

  function vwn_energy_sub_expanded(n, p) result(v)

    real(kind_dp), intent(in) :: n
    real(kind_dp), intent(in) :: p(4)
    real(kind_dp)             :: v
!   ----------------------------------------------------------------------------
    real(kind_dp)             :: a, b, c, s, x, y, z

    a = p(3)*p(1)/(p(1)*p(1) + p(1)*p(3) + p(4)) - 1.0d0
    b = 2.0d0*(p(3)*p(1)/(p(1)*p(1) + p(1)*p(3) + p(4)) - 1.0d0) + 2.0d0
    c = 2.0d0*p(3)*(1.0d0/sqrt(4.0d0*p(4) - p(3)*p(3)) - p(1)/((p(1)*p(1)   & 
      + p(1)*p(3) + p(4))*sqrt(4.0d0*p(4) - p(3)*p(3))/(p(3) + 2.0d0*p(1))))

    s = (3.0d0/(4.0d0*pi))**(1.0d0/6.0d0)/sqrt(n**(1.0d0/3.0d0))

    x = s*s + p(3)*s+p(4)
    y = s - p(1)
    z = sqrt(4.0d0*p(4) - p(3)*p(3))/(2.0d0*s + p(3))

    v = 0.5d0*p(2)*(2.0d0*log(s) + a*log(x) - b*log(y) + c*atan(z))

  end function

  function vwn_energy_sub(n, p) result(v)

    real(kind_dp), intent(in) :: n
    real(kind_dp), intent(in) :: p(4)
    real(kind_dp)             :: v
!   ----------------------------------------------------------------------------
    real(kind_dp)             :: a, b, c, d, f, g, h, i, r, s, t, u, z, e

    a = p(1)*p(1) + p(1)*p(3) + p(4)
    b = p(1)/a
    c = b*p(3) - 1.0d0
    d = 2.0d0*c + 2.0d0

    f = sqrt(4.0d0*p(4) - p(3)*p(3))
    g = f/(p(3) + 2.0d0*p(1))
    h = a*g
    i = 2.0d0*p(3)*(1.0d0/f - p(1)/h)

    r = (3.0d0/(4.0d0*pi))**(1.0d0/6.0d0)
    s = r/sqrt(n**(1.0d0/3.0d0))
    t = 2.0d0*s + p(3)
    u = f/t
    z = s*s + p(3)*s + p(4)

    e = 2.0d0*log(s) + c*log(z) - d*log(s - p(1)) + i*atan(u)
    v = 0.5d0*p(2)*e

  end function

  function vwn_energy(na, nb, w) result(v)

    real(kind_dp), intent(in) :: na, nb, w
    real(kind_dp)             :: v
!   ----------------------------------------------------------------------------
    real(kind_dp)             :: n
    real(kind_dp)             :: v_para, v_ferro, v_inter
    real(kind_dp)             :: para(4), ferro(4), inter(4)
    real(kind_dp)             :: z, z4, fz, d

    n = na + nb

    para(1)   = -0.1049800d0
    para(2)   =  0.0621814d0
    para(3)   =  3.72744d0
    para(4)   = 12.9352d0
              
    ferro(1)  = -0.3250000d0
    ferro(2)  =  0.0310907d0
    ferro(3)  =  7.06042d0
    ferro(4)  = 18.0578d0
              
    inter(1)  = -0.0047584d0
    inter(2)  = -0.0337737d0
    inter(3)  =  1.13107d0
    inter(4)  = 13.0045d0
              
    v_para    = vwn_energy_sub(n, para)
!   v_para    = vwn_energy_sub_expanded(n, para)
    v         = v_para

    if (abs(na - nb) > small) then
      v_ferro = vwn_energy_sub(n, ferro)
      v_inter = vwn_energy_sub(n, inter)
    
      z       = (na - nb)/n
      z4      = z**4.0d0
      fz      = 1.92366105093154d0*((1.0d0 + z)**(4.0d0/3.0d0) + (1.0d0 - z)**(4.0d0/3.0d0) - 2.0d0)
      d       = fz*((v_ferro - v_para)*z4 + v_inter*(1.0d0 - z4)*0.584822305543806d0)

      v       = v + d
    end if

    v         = w*n*v

  end function

end module
