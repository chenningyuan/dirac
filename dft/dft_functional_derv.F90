module dft_functional_derv

  use dft_cfg
  use dft_basis_info

  implicit none

  public get_functional_derv

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
#ifdef MOD_AUTODIFF
#include "../autodiff/xc_fortran_interface.h"
#endif

  real(8), parameter :: tiny    = 1.0d-20
  integer, parameter :: max_dim = 99

contains

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
    real(8)              :: t_lda(max_dim, mxbllen)
    real(8)              :: t_gga(max_dim, mxbllen)
    real(8)              :: f, taylor_deriv, dftenergy
    integer              :: i, k, l
!   ----------------------------------------------------------------------------

    derv = 0.0d0

#ifdef MOD_AUTODIFF
    if (dft_cfg_autodiff) then

!     lda
      t_lda = 0.0d0
      do i = 1, size(dft_cfg_lda_functional_weight_array)
        f = dft_cfg_lda_functional_weight_array(i)
        if (dabs(f) > tiny) then
          do k = 1, block_length
            call lda_eval(i, order, r_0(k), t(1, k))
            call daxpy(max_dim, f*w(k), t(1, k), 1, t_lda(1, k), 1)
          end do
        end if
      end do

!     gga
      t_gga = 0.0d0
      do i = 1, size(dft_cfg_gga_functional_weight_array)
        f = dft_cfg_gga_functional_weight_array(i)
        if (dabs(f) > tiny) then
          do k = 1, block_length
            call gga_eval(i, order, r_0(k), z_0(k), t(1, k))
            call daxpy(max_dim, f*w(k), t(1, k), 1, t_gga(1, k), 1)
          end do
        end if
      end do

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = t_lda(1, k) + t_gga(1, k)
          end do

        case (1)
          do k = 1, block_length
            derv(k, d_r) = taylor_deriv(t_lda(1, k), 2, (/1, 0/)) &
                         + taylor_deriv(t_gga(1, k), 5, (/1, 0, 0, 0, 0/))
            derv(k, d_z) = taylor_deriv(t_gga(1, k), 5, (/0, 0, 1, 0, 0/))
          end do

        case (2)
          do k = 1, block_length
            derv(k, d_rr) = taylor_deriv(t_lda(1, k), 2, (/2, 0/)) &
                          + taylor_deriv(t_gga(1, k), 5, (/2, 0, 0, 0, 0/))
            derv(k, d_ss) = taylor_deriv(t_lda(1, k), 2, (/0, 2/)) &
                          + taylor_deriv(t_gga(1, k), 5, (/0, 2, 0, 0, 0/))
            derv(k, d_rz) = taylor_deriv(t_gga(1, k), 5, (/1, 0, 1, 0, 0/))
            derv(k, d_zz) = taylor_deriv(t_gga(1, k), 5, (/0, 0, 2, 0, 0/))
            derv(k, d_sy) = taylor_deriv(t_gga(1, k), 5, (/0, 1, 0, 0, 1/))
            derv(k, d_yy) = taylor_deriv(t_gga(1, k), 5, (/0, 0, 0, 0, 2/))
            derv(k, d_z)  = taylor_deriv(t_gga(1, k), 5, (/0, 0, 1, 0, 0/))
            derv(k, d_x)  = taylor_deriv(t_gga(1, k), 5, (/0, 0, 0, 1, 0/))
          end do

        case (3)
          do k = 1, block_length
            derv(k, d_rrr) = taylor_deriv(t_lda(1, k), 2, (/3, 0/)) &
                           + taylor_deriv(t_gga(1, k), 5, (/3, 0, 0, 0, 0/))
            derv(k, d_rss) = taylor_deriv(t_lda(1, k), 2, (/1, 2/)) &
                           + taylor_deriv(t_gga(1, k), 5, (/1, 2, 0, 0, 0/))
            derv(k, d_rrz) = taylor_deriv(t_gga(1, k), 5, (/2, 0, 1, 0, 0/))
            derv(k, d_rsy) = taylor_deriv(t_gga(1, k), 5, (/1, 1, 0, 0, 1/))
            derv(k, d_rx)  = taylor_deriv(t_gga(1, k), 5, (/1, 0, 0, 1, 0/))
            derv(k, d_ryy) = taylor_deriv(t_gga(1, k), 5, (/1, 0, 0, 0, 2/))
            derv(k, d_rz)  = taylor_deriv(t_gga(1, k), 5, (/1, 0, 1, 0, 0/))
            derv(k, d_rzz) = taylor_deriv(t_gga(1, k), 5, (/1, 0, 2, 0, 0/))
            derv(k, d_ssz) = taylor_deriv(t_gga(1, k), 5, (/0, 2, 1, 0, 0/))
            derv(k, d_sy)  = taylor_deriv(t_gga(1, k), 5, (/0, 1, 0, 0, 1/))
            derv(k, d_szy) = taylor_deriv(t_gga(1, k), 5, (/0, 1, 1, 0, 1/))
            derv(k, d_yy)  = taylor_deriv(t_gga(1, k), 5, (/0, 0, 0, 0, 2/))
            derv(k, d_zx)  = taylor_deriv(t_gga(1, k), 5, (/0, 0, 1, 1, 0/))
            derv(k, d_zyy) = taylor_deriv(t_gga(1, k), 5, (/0, 0, 1, 0, 2/))
            derv(k, d_zz)  = taylor_deriv(t_gga(1, k), 5, (/0, 0, 2, 0, 0/))
            derv(k, d_zzz) = taylor_deriv(t_gga(1, k), 5, (/0, 0, 3, 0, 0/))
          end do

        case default
          call quit('get_functional_derv: order too high')

      end select

    else
#endif /* ifdef MOD_AUTODIFF */

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = w(k)*dftenergy(r_0(k), dsqrt(z_0(k)))
          end do

        case (1)
          do k = 1, block_length
            call dftpot0(t(1, k), w(k), r_0(k), dsqrt(z_0(k)))
          end do

        case (2)
          do k = 1, block_length
            call sdftpot1(t(1, k), w(k), r_0(k), dsqrt(z_0(k)))
          end do

        case (3)
          do k = 1, block_length
            call sdftpot2(t(1, k), w(k), r_0(k), dsqrt(z_0(k)))
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

  end subroutine

end module
