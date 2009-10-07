module dft_functional_derv

  use dft_cfg
#ifdef MOD_AUTODIFF
  use xc_fun
#endif

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
    integer              :: k, l
    integer, external    :: taylor_index
    real(8), external    :: dftenergy
!   ----------------------------------------------------------------------------

    derv = 0.0d0

#ifdef MOD_AUTODIFF
    if (dft_cfg_autodiff) then

      do k = 1, block_length
        call eval_gga(t(1, k), order, (/r_0(k), 0.0d0, z_0(k), 0.0d0, 0.0d0/))
      end do

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = t(1, k)
          end do

        case (1)
          do k = 1, block_length
            derv(k, d_r) = t(taylor_index(5, (/1, 0, 0, 0, 0/)), k)
            derv(k, d_z) = t(taylor_index(5, (/0, 0, 1, 0, 0/)), k)
          end do

        case (2)
          do k = 1, block_length
            derv(k, d_rr) = t(taylor_index(5, (/2, 0, 0, 0, 0/)), k)
            derv(k, d_ss) = t(taylor_index(5, (/0, 2, 0, 0, 0/)), k)
            derv(k, d_rz) = t(taylor_index(5, (/1, 0, 1, 0, 0/)), k)
            derv(k, d_zz) = t(taylor_index(5, (/0, 0, 2, 0, 0/)), k)
            derv(k, d_sy) = t(taylor_index(5, (/0, 1, 0, 0, 1/)), k)
            derv(k, d_yy) = t(taylor_index(5, (/0, 0, 0, 0, 2/)), k)
            derv(k, d_z ) = t(taylor_index(5, (/0, 0, 1, 0, 0/)), k)
            derv(k, d_x ) = t(taylor_index(5, (/0, 0, 0, 1, 0/)), k)
          end do

        case (3)
          do k = 1, block_length
            derv(k, d_rrr) = t(taylor_index(5, (/3, 0, 0, 0, 0/)), k)
            derv(k, d_rss) = t(taylor_index(5, (/1, 2, 0, 0, 0/)), k)
            derv(k, d_rrz) = t(taylor_index(5, (/2, 0, 1, 0, 0/)), k)
            derv(k, d_rsy) = t(taylor_index(5, (/1, 1, 0, 0, 1/)), k)
            derv(k, d_rx ) = t(taylor_index(5, (/1, 0, 0, 1, 0/)), k)
            derv(k, d_ryy) = t(taylor_index(5, (/1, 0, 0, 0, 2/)), k)
            derv(k, d_rz ) = t(taylor_index(5, (/1, 0, 1, 0, 0/)), k)
            derv(k, d_rzz) = t(taylor_index(5, (/1, 0, 2, 0, 0/)), k)
            derv(k, d_ssz) = t(taylor_index(5, (/0, 2, 1, 0, 0/)), k)
            derv(k, d_sy ) = t(taylor_index(5, (/0, 1, 0, 0, 1/)), k)
            derv(k, d_szy) = t(taylor_index(5, (/0, 1, 1, 0, 1/)), k)
            derv(k, d_yy ) = t(taylor_index(5, (/0, 0, 0, 0, 2/)), k)
            derv(k, d_zx ) = t(taylor_index(5, (/0, 0, 1, 1, 0/)), k)
            derv(k, d_zyy) = t(taylor_index(5, (/0, 0, 1, 0, 2/)), k)
            derv(k, d_zz ) = t(taylor_index(5, (/0, 0, 2, 0, 0/)), k)
            derv(k, d_zzz) = t(taylor_index(5, (/0, 0, 3, 0, 0/)), k)
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
