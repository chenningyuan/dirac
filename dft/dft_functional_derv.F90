module dft_functional_derv

  use dft_cfg
#ifdef MOD_AUTODIFF
  use xc_fun
#endif

  implicit none

  public get_functional_derv
#ifdef MOD_AUTODIFF
  public xc_fun_set_functional
#endif

  integer, parameter, public :: d_0   =  0
  integer, parameter, public :: d_r   =  1
  integer, parameter, public :: d_z   =  2
  integer, parameter, public :: d_x   =  3
  integer, parameter, public :: d_rr  =  4
  integer, parameter, public :: d_ss  =  5
  integer, parameter, public :: d_zz  =  6
  integer, parameter, public :: d_yy  =  7
  integer, parameter, public :: d_rz  =  8
  integer, parameter, public :: d_sy  =  9
  integer, parameter, public :: d_rx  = 10
  integer, parameter, public :: d_zx  = 11
  integer, parameter, public :: d_rrr = 12
  integer, parameter, public :: d_zzz = 13
  integer, parameter, public :: d_rss = 14
  integer, parameter, public :: d_rzz = 15
  integer, parameter, public :: d_ryy = 16
  integer, parameter, public :: d_rrz = 17
  integer, parameter, public :: d_ssz = 18
  integer, parameter, public :: d_rsy = 19
  integer, parameter, public :: d_szy = 20
  integer, parameter, public :: d_zyy = 21

  integer, parameter, public :: nr_nonzero_derv(3) = (/2, 11, 21/)

  private

#include "functionals.h"

  real(8), parameter :: tiny    = 1.0d-20
  integer, parameter :: max_dim = 99

contains

#ifdef MOD_AUTODIFF
  subroutine xc_fun_set_functional(order)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: order
!   ----------------------------------------------------------------------------
    integer :: i
    integer :: gga_functionals(4) = (/XC_BECKE_EXCHANGE,  &
                                      XC_LYP_CORRELATION, &
                                      XC_PBE_EXCHANGE,    &
                                      XC_PBE_CORRELATION/)
!   ----------------------------------------------------------------------------

    do i = 1, size(gga_functionals)
      if (dabs(dft_cfg_functional_w(gga_functionals(i))) > tiny) then
        dft_cfg_is_gga_autodiff = .true.
      end if
    end do

    call xc_set_functional(XC_RSZXY, dft_cfg_functional_w)

  end subroutine
#endif /* ifdef MOD_AUTODIFF */

  subroutine get_functional_derv(order,        &
                                 block_length, &
                                 w,            &
                                 r_0,          &
                                 z_0,          &
                                 n,            &
                                 derv)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: order
    integer, intent(in)  :: block_length
    real(8), intent(in)  :: w(*)
    real(8), intent(in)  :: r_0(*)
    real(8), intent(in)  :: z_0(*)
    integer, intent(in)  :: n
    real(8), intent(out) :: derv(mxbllen, 0:n)
!   ----------------------------------------------------------------------------
    real(8)              :: t(max_dim, mxbllen)
    integer              :: k, l
    real(8), external    :: dftenergy
!   ----------------------------------------------------------------------------

    derv = 0.0d0

#ifdef MOD_AUTODIFF
    if (dft_cfg_autodiff) then

      do k = 1, block_length
        call xc_eval(t(1, k), order, (/r_0(k), 0.0d0, z_0(k), 0.0d0, 0.0d0/))
      end do

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = t(1, k)
          end do

        case (1)
          do k = 1, block_length
            derv(k, d_r) = t(xc_index(XC_RSZXY, (/1, 0, 0, 0, 0/)), k)
            derv(k, d_z) = t(xc_index(XC_RSZXY, (/0, 0, 1, 0, 0/)), k)
          end do

        case (2)
          do k = 1, block_length
            derv(k, d_rr) = t(xc_index(XC_RSZXY, (/2, 0, 0, 0, 0/)), k)
            derv(k, d_ss) = t(xc_index(XC_RSZXY, (/0, 2, 0, 0, 0/)), k)
            derv(k, d_rz) = t(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)), k)
            derv(k, d_zz) = t(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)), k)
            derv(k, d_sy) = t(xc_index(XC_RSZXY, (/0, 1, 0, 0, 1/)), k)
            derv(k, d_yy) = t(xc_index(XC_RSZXY, (/0, 0, 0, 0, 2/)), k)
            derv(k, d_z ) = t(xc_index(XC_RSZXY, (/0, 0, 1, 0, 0/)), k)
            derv(k, d_x ) = t(xc_index(XC_RSZXY, (/0, 0, 0, 1, 0/)), k)
          end do

        case (3)
          do k = 1, block_length
            derv(k, d_rrr) = t(xc_index(XC_RSZXY, (/3, 0, 0, 0, 0/)), k)
            derv(k, d_rss) = t(xc_index(XC_RSZXY, (/1, 2, 0, 0, 0/)), k)
            derv(k, d_rrz) = t(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)), k)
            derv(k, d_rsy) = t(xc_index(XC_RSZXY, (/1, 1, 0, 0, 1/)), k)
            derv(k, d_rx ) = t(xc_index(XC_RSZXY, (/1, 0, 0, 1, 0/)), k)
            derv(k, d_ryy) = t(xc_index(XC_RSZXY, (/1, 0, 0, 0, 2/)), k)
            derv(k, d_rz ) = t(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)), k)
            derv(k, d_rzz) = t(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)), k)
            derv(k, d_ssz) = t(xc_index(XC_RSZXY, (/0, 2, 1, 0, 0/)), k)
            derv(k, d_sy ) = t(xc_index(XC_RSZXY, (/0, 1, 0, 0, 1/)), k)
            derv(k, d_szy) = t(xc_index(XC_RSZXY, (/0, 1, 1, 0, 1/)), k)
            derv(k, d_yy ) = t(xc_index(XC_RSZXY, (/0, 0, 0, 0, 2/)), k)
            derv(k, d_zx ) = t(xc_index(XC_RSZXY, (/0, 0, 1, 1, 0/)), k)
            derv(k, d_zyy) = t(xc_index(XC_RSZXY, (/0, 0, 1, 0, 2/)), k)
            derv(k, d_zz ) = t(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)), k)
            derv(k, d_zzz) = t(xc_index(XC_RSZXY, (/0, 0, 3, 0, 0/)), k)
          end do

        case default
          call quit('get_functional_derv: order too high')

      end select

    else
#endif /* ifdef MOD_AUTODIFF */

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = dftenergy(r_0(k), dsqrt(z_0(k)))
          end do

        case (1)
          do k = 1, block_length
            call dftpot0(t(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case (2)
          do k = 1, block_length
            call sdftpot1(t(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case (3)
          do k = 1, block_length
            call sdftpot2(t(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case default
          call quit('get_functional_derv: order too high')

      end select

      do l = 1, n
        do k = 1, block_length
          derv(k, l) = derv(k, l) + t(l, k)
        end do
      end do

#ifdef MOD_AUTODIFF
    end if
#endif

!   multiply by weight
    do l = 0, n
      do k = 1, block_length
        derv(k, l) = w(k)*derv(k, l)
      end do
    end do

  end subroutine

end module
