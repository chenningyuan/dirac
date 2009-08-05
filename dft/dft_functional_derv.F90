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

  real(8), parameter :: tiny = 1.0d-20

contains

  subroutine get_functional_derv(order, &
                                 w,     &
                                 r_0,   &
                                 z_0,   &
                                 n,     &
                                 derv)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: order
    real(8), intent(in)  :: w
    real(8), intent(in)  :: r_0
    real(8), intent(in)  :: z_0
    integer, intent(in)  :: n
    real(8), intent(out) :: derv(0:n)
!   ----------------------------------------------------------------------------
    real(8)              :: t(99), t_lda(99), t_gga(99)
    real(8)              :: f, taylor_deriv, dftenergy
    integer              :: i
!   ----------------------------------------------------------------------------

    derv = 0.0d0

#ifdef MOD_AUTODIFF
    if (dft_cfg_autodiff) then

!     lda
      t_lda = 0.0d0
      do i = 1, size(dft_cfg_lda_functional_weight_array)
        f = dft_cfg_lda_functional_weight_array(i)
        if (dabs(f) > tiny) then
          call lda_eval(i, order, r_0, t)
          t_lda = t_lda + f*t
        end if
      end do

!     gga
      t_gga = 0.0d0
      do i = 1, size(dft_cfg_gga_functional_weight_array)
        f = dft_cfg_gga_functional_weight_array(i)
        if (dabs(f) > tiny) then
          call gga_eval(i, order, r_0, z_0, t)
!         daxpy me
          t_gga = t_gga + f*t
        end if
      end do

      select case (order)

        case (0)
          derv(d_0) = t_lda(1) + t_gga(1)

        case (1)
          derv(d_r) = taylor_deriv(t_lda, 2, (/1, 0/)) &
                    + taylor_deriv(t_gga, 5, (/1, 0, 0, 0, 0/))
          derv(d_z) = taylor_deriv(t_gga, 5, (/0, 0, 1, 0, 0/))

        case default
          call quit('get_functional_derv: order too high')

      end select

    else
#endif /* ifdef MOD_AUTODIFF */

      select case (order)

        case (0)
          derv(d_0) = dftenergy(r_0, dsqrt(z_0))

        case (1)
          call dftpot0(derv(d_r), 1.0d0, r_0, dsqrt(z_0))

        case default
          call quit('get_functional_derv: order too high')

      end select

#ifdef MOD_AUTODIFF
    end if
#endif

    derv = w*derv

  end subroutine

end module
