module dft_functional_derv

  use dft_cfg
#ifdef MOD_XCFUN
  use xcfun
#endif

  implicit none

  public get_functional_derv
#ifdef MOD_XCFUN
  public xcfun_set_functional
  public xcfun_is_gga
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

  real(8), parameter   :: tiny    = 1.0d-20
  integer, parameter   :: max_dim = 99

  real(8), allocatable :: xcfun_in(:)
  real(8), allocatable :: xcfun_out(:, :)
  integer              :: xcfun_in_length  = 0
  integer              :: xcfun_out_length = 0

contains

#ifdef MOD_XCFUN

  function xcfun_is_gga()

!   ----------------------------------------------------------------------------
    logical :: xcfun_is_gga
!   ----------------------------------------------------------------------------

    xcfun_is_gga = .false.
    if (xc_get_type(dft_cfg_fun_id) > 0) then
      xcfun_is_gga = .true.
    end if

  end function

  subroutine xcfun_set_functional(response_order, response_order_ao)
!fixme: rename, functional is not set here, only room is made for in/out buffers

!   ----------------------------------------------------------------------------
    integer, intent(in) :: response_order
    integer, intent(in) :: response_order_ao
!   ----------------------------------------------------------------------------
    integer             :: i, alda_order
    integer             :: need_in_length
    integer             :: need_out_length
!   ----------------------------------------------------------------------------

    need_in_length  = xc_input_length(dft_cfg_fun_id)
    need_out_length = xc_output_length(dft_cfg_fun_id, response_order)

    if (need_in_length > xcfun_in_length) then
      xcfun_in_length = need_in_length
      if (allocated(xcfun_in)) then
        deallocate(xcfun_in)        
      end if
      allocate(xcfun_in(xcfun_in_length))
    end if

    if (need_out_length > xcfun_out_length) then
      xcfun_out_length = need_out_length
      if (allocated(xcfun_out)) then
        deallocate(xcfun_out)        
      end if
      allocate(xcfun_out(xcfun_out_length, mxbllen))
    end if

!   alda_order = dft_cfg_alda_order
!   if (dft_cfg_alda_hs .and. dft_cfg_alda_ha) then
!     alda_order = 1
!   end if

!   if ((alda_order <= response_order) .or. (alda_order <= response_order_ao)) then
!     if (dft_cfg_xalda) then
!       call xc_set_functional(XC_RSZXY, dft_cfg_functional_w_xalda)
!     else
!       call xc_set_functional(XC_RSZXY, dft_cfg_functional_w_alda)
!     end if
!     dft_cfg_is_gga_xcfun = .false.
!   else
!     call xc_set_functional(XC_RSZXY, dft_cfg_functional_w)
!radovan: fix later
!     do i = 1, size(gga_functionals)
!       if (dabs(dft_cfg_functional_w(gga_functionals(i))) > tiny) then
!         dft_cfg_is_gga_xcfun = .true.
!       end if
!     end do
!   end if

  end subroutine
#endif /* ifdef MOD_XCFUN */

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

#ifdef MOD_XCFUN
    if (dft_cfg_xcfun) then

      do k = 1, block_length
        xcfun_in = 0.0d0
        xcfun_in(1) = r_0(k)
        if (xcfun_in_length > 2) then
          xcfun_in(3) = z_0(k)
        end if
        call xc_eval(dft_cfg_fun_id, xcfun_out(1, k), order, xcfun_in)
      end do

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = xcfun_out(1, k)
          end do

        case (1)
          do k = 1, block_length
            derv(k, d_0) = xcfun_out(1, k)
            derv(k, d_r) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 0, 0, 0, 0/)), k)
            derv(k, d_z) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 1, 0, 0/)), k)
          end do

        case (2)
!fixme
!#define PARTIAL_ALDA_WITH_XCFUN
#ifdef PARTIAL_ALDA_WITH_XCFUN
          if      (dft_cfg_alda_hs .and. .not. dft_cfg_alda_ha) then
            call xc_set_functional(XC_RSZXY, dft_cfg_functional_w_alda)
            do k = 1, block_length
              call xc_eval(t(1, k), order, (/r_0(k), 0.0d0, z_0(k), 0.0d0, 0.0d0/))
              derv(k, d_rr) = t(xc_index(XC_RSZXY, (/2, 0, 0, 0, 0/)), k)
            end do
            call xc_set_functional(XC_RSZXY, dft_cfg_functional_w)
            do k = 1, block_length
              call xc_eval(t(1, k), order, (/r_0(k), 0.0d0, z_0(k), 0.0d0, 0.0d0/))
              derv(k, d_ss) = t(xc_index(XC_RSZXY, (/0, 2, 0, 0, 0/)), k)
              derv(k, d_sy) = t(xc_index(XC_RSZXY, (/0, 1, 0, 0, 1/)), k)
              derv(k, d_yy) = t(xc_index(XC_RSZXY, (/0, 0, 0, 0, 2/)), k)
              derv(k, d_x ) = t(xc_index(XC_RSZXY, (/0, 0, 0, 1, 0/)), k)
            end do
          else if (dft_cfg_alda_ha .and. .not. dft_cfg_alda_hs) then
            call xc_set_functional(XC_RSZXY, dft_cfg_functional_w_alda)
            do k = 1, block_length
              call xc_eval(t(1, k), order, (/r_0(k), 0.0d0, z_0(k), 0.0d0, 0.0d0/))
              derv(k, d_ss) = t(xc_index(XC_RSZXY, (/0, 2, 0, 0, 0/)), k)
            end do
            call xc_set_functional(XC_RSZXY, dft_cfg_functional_w)
            do k = 1, block_length
              call xc_eval(t(1, k), order, (/r_0(k), 0.0d0, z_0(k), 0.0d0, 0.0d0/))
              derv(k, d_rr) = t(xc_index(XC_RSZXY, (/2, 0, 0, 0, 0/)), k)
              derv(k, d_rz) = t(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)), k)
              derv(k, d_zz) = t(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)), k)
              derv(k, d_z ) = t(xc_index(XC_RSZXY, (/0, 0, 1, 0, 0/)), k)
            end do
          else
#endif /* ifdef PARTIAL_ALDA_WITH_XCFUN */
          do k = 1, block_length
            derv(k, d_rr) = xcfun_out(xc_index(dft_cfg_fun_id, (/2, 0, 0, 0, 0/)), k)
            derv(k, d_ss) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 2, 0, 0, 0/)), k)
            derv(k, d_rz) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 0, 1, 0, 0/)), k)
            derv(k, d_zz) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 2, 0, 0/)), k)
            derv(k, d_sy) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 1, 0, 0, 1/)), k)
            derv(k, d_yy) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 0, 0, 2/)), k)
            derv(k, d_z ) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 1, 0, 0/)), k)
            derv(k, d_x ) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 0, 1, 0/)), k)
          end do
#ifdef PARTIAL_ALDA_WITH_XCFUN
          end if
#endif /* ifdef PARTIAL_ALDA_WITH_XCFUN */

        case (3)
          do k = 1, block_length
            derv(k, d_rrr) = xcfun_out(xc_index(dft_cfg_fun_id, (/3, 0, 0, 0, 0/)), k)
            derv(k, d_rss) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 2, 0, 0, 0/)), k)
            derv(k, d_rrz) = xcfun_out(xc_index(dft_cfg_fun_id, (/2, 0, 1, 0, 0/)), k)
            derv(k, d_rsy) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 1, 0, 0, 1/)), k)
            derv(k, d_rx ) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 0, 0, 1, 0/)), k)
            derv(k, d_ryy) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 0, 0, 0, 2/)), k)
            derv(k, d_rz ) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 0, 1, 0, 0/)), k)
            derv(k, d_rzz) = xcfun_out(xc_index(dft_cfg_fun_id, (/1, 0, 2, 0, 0/)), k)
            derv(k, d_ssz) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 2, 1, 0, 0/)), k)
            derv(k, d_sy ) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 1, 0, 0, 1/)), k)
            derv(k, d_szy) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 1, 1, 0, 1/)), k)
            derv(k, d_yy ) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 0, 0, 2/)), k)
            derv(k, d_zx ) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 1, 1, 0/)), k)
            derv(k, d_zyy) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 1, 0, 2/)), k)
            derv(k, d_zz ) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 2, 0, 0/)), k)
            derv(k, d_zzz) = xcfun_out(xc_index(dft_cfg_fun_id, (/0, 0, 3, 0, 0/)), k)
          end do

        case default
          call quit('get_functional_derv: order too high')

      end select

    else
#endif /* ifdef MOD_XCFUN */

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = dftenergy(r_0(k), dsqrt(z_0(k)))
          end do

        case (1)
          do k = 1, block_length
            derv(k, d_0) = dftenergy(r_0(k), dsqrt(z_0(k)))
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

#ifdef MOD_XCFUN
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
