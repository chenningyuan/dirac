module xc_response

#ifdef MOD_AUTODIFF
  use dft_cfg
  use dft_density
  use dft_cb_interface
  use aor_cfg
  use xc_fun

  implicit none

  public xc_response_contribution

  private
    real(8), parameter     :: tiny = 1.0d-30
    integer, parameter     :: max_nr_dmat = 20

!    "the universe works on a math equation
!     that never even ever really ends in the end"
!                      modest mouse - never ending math equation - building nothing out of something


!     we have five density variables:
!       n
!       s
!       z = \nabla n \cdot nabla n
!       y = \nabla n \cdot nabla s
!       x = \nabla s \cdot nabla s
!
!     response contribution blocks:
!         n    s    z    y    x
!        _________________________
!     n | n    ns   nz   rest rest
!     s | ns   ns   sz   rest rest
!     z | nz   sz   nz   rest rest
!     y | rest rest rest rest rest
!     x | rest rest rest rest rest
!
!     \nabla n contribution | s contribution | \nabla s contribution | blocks
!     -----------------------------------------------------------------------
!     no                    | no             | no                    | n
!     no                    | yes            | no                    | n + ns
!     yes                   | no             | no                    | n + nz
!     yes                   | yes            | no                    | n + nz + ns + sz
!     yes                   | yes            | yes                   | n + nz + ns + sz + rest

  logical :: calculate_blocks_n    = .true.
  logical :: calculate_blocks_ns   = .false.
  logical :: calculate_blocks_nz   = .true.
  logical :: calculate_blocks_sz   = .false.
  logical :: calculate_blocks_rest = .false.

contains

  subroutine xc_response_contribution(order,   &
                         nr_dmat, &
                         use_gga, &
                         dmat,    &
                         fmat,    &
                         isymop,  &
                         ao_g0_m0,     &
                         ao_g1_m0,    &
                         n_0,     &
                         gn_0,    &
                         w,       &
                         buffer)

#include "dgroup.h"
#include "dcbbas.h"

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: order
    integer, intent(in)    :: nr_dmat
    logical, intent(in)    :: use_gga
    real(8), intent(in)    :: dmat(ntbas(0), ntbas(0), nz, nr_dmat)
    real(8), intent(inout) :: fmat(ntbas(0), ntbas(0), nz)
    integer, intent(in)    :: isymop(nr_dmat)
    real(8), intent(in)    :: ao_g0_m0(*)
    real(8), intent(in)    :: ao_g1_m0(*)
    real(8), intent(in)    :: n_0
    real(8), intent(in)    :: gn_0(3)
    real(8), intent(in)    :: w
    real(8), intent(in)    :: buffer(*)
!   ----------------------------------------------------------------------------
    real(8)                :: u
    real(8)                :: v(3)

    real(8)                :: derv(999)

    real(8)                :: n(0:max_nr_dmat)
    real(8)                :: gn(3, 0:max_nr_dmat)
    real(8)                :: z(0:max_nr_dmat, 0:max_nr_dmat)

    integer                :: i, j
    logical                :: use_alda
    logical                :: u_above_threshold
    logical                :: v_above_threshold

    real(8), external      :: ddot
!   ----------------------------------------------------------------------------

!   get perturbed densities

    n(0) = n_0
    do i = 1, nr_dmat
      call get_n(n(i),             &
                 isymop(i) - 1,    &
                 dmat(1, 1, 1, i), &
                 buffer,           &
                 ao_g0_m0)
    end do

!   get perturbed densitiy gradients and gradient norms

    if (use_gga) then

      gn(:, 0) = gn_0
      do i = 1, nr_dmat
        call get_gn(gn(1, i),         &
                    isymop(i) - 1,    &
                    dmat(1, 1, 1, i), &
                    buffer,           &
                    ao_g0_m0,         &
                    ao_g1_m0)
      end do

      do i = 0, nr_dmat
        do j = i, nr_dmat
!         the perturbed z carries factor 2
!         because derivative of z gives 2 identical terms
!         (one "left", one "right")
          z(i, j) = 2.0d0*ddot(3, gn(1, i), 1, gn(1, j), 1)
        end do
      end do
!     correct wrong prefactor 2 of the unperturbed z
      z(0, 0) = 0.5d0*z(0, 0)

    end if

    use_alda = .false.
    if (order >= aor_cfg_alda_order) then
      use_alda = .true.
    end if

    if (.not. use_gga) then
      use_alda = .true.
    end if

    derv = 0.0d0
    call xc_eval(derv, order, (/n(0), 0.0d0, z(0, 0), 0.0d0, 0.0d0/))

    u = 0.0d0
    v = 0.0d0

    select case (order)

      case (2)
        call xc_response_2(u, v, derv, n, gn, z, use_alda)

      case (3)
        call xc_response_3(u, v, derv, n, gn, z, use_alda)

      case (4)
        call xc_response_4(u, v, derv, n, gn, z, use_alda)

      case (5)
        call xc_response_5(u, v, derv, n, gn, z, use_alda)

      case default
        call quit('xc_response_contribution: order not implemented')

    end select

!#define WRITE_PERTURBED_DENSITIES
#ifdef WRITE_PERTURBED_DENSITIES
    if (n(0) > 0.1d0) then
      do i = 0, nr_dmat
        write(*, *) 'debug n', i, n(i)
      end do
      write(*, *) 'debug z(0)', z(0, 0)
      select case (order)

! - | 1 2 3  4  5  6  7   8  9  10 11  12  13  14  15
! 2 | i
! 3 | i j ij
! 4 | i j k  ij ik jk ijk
! 5 | i j k  l  ij ik il  jk jl kl ijk ijl ikl jkl ijkl

        case (2)
          write(*, *) 'debug z(1)', z(0, 1)

        case (3)
          write(*, *) 'debug z(1)', z(0, 1)

        case (4)
          write(*, *) 'debug z(1)', z(0, 1)
          write(*, *) 'debug z(2)', z(0, 4) + z(1, 2)

        case (5)
          write(*, *) 'debug z(1)', z(0, 1)
          write(*, *) 'debug z(2)', z(0, 4)  + z(1, 2)
          write(*, *) 'debug z(3)', z(0, 11) + z(1, 8)

      end select
      write(*, *) 'debug u', order, u
    end if
#endif

!   multiply by the weight
    u = u*w
    v = v*w

    u_above_threshold = .false.
    if (dabs(u) > tiny) then
      u_above_threshold = .true.
    end if

    v_above_threshold = .false.
    if (maxval((/dabs(v(1)), dabs(v(2)), dabs(v(3))/)) > tiny) then
      v_above_threshold = .true.
    end if

    if (v_above_threshold) then
      call gga_distribute(u, 2.0d0*v, 0, fmat, ao_g0_m0, ao_g1_m0)
    else
      if (u_above_threshold) then
        call lda_distribute(u, 0, fmat, ao_g0_m0)
      end if
    end if

  end subroutine

  subroutine xc_response_2(u, v, derv, n, gn, z, use_alda)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_alda
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tu = z(0, 1)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tu
    tv = gn(:, 1)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 1, 0, 0/)))*tv(:)
    end if

  end subroutine

  subroutine xc_response_3(u, v, derv, n, gn, z, use_alda)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_alda
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)*n(2)
    u  = u + derv(xc_index(XC_RSZXY, (/3, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tv = z(0, 1)*gn(:, 2) &
       + gn(:, 1)*z(0, 2)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2) &
       + z(0, 1)*n(2)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tu
    tu = z(0, 1)*z(0, 2)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tu
    tv = n(1)*gn(:, 2) &
       + gn(:, 1)*n(2)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 2)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tu
    end if

  end subroutine

  subroutine xc_response_4(u, v, derv, n, gn, z, use_alda)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_alda
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)*n(2)*n(3)
    u  = u + derv(xc_index(XC_RSZXY, (/4, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/3, 0, 1, 0, 0/)))*tu*gn(:, 0)
    end if
    tu = n(4)*n(3) &
       + n(5)*n(2) &
       + n(6)*n(1)
    u  = u + derv(xc_index(XC_RSZXY, (/3, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tv = z(0, 4)*gn(:, 3) &
       + z(1, 2)*gn(:, 3) &
       + z(0, 5)*gn(:, 2) &
       + z(1, 3)*gn(:, 2) &
       + z(0, 6)*gn(:, 1) &
       + z(2, 3)*gn(:, 1) &
       + gn(:, 4)*z(0, 3) &
       + gn(:, 5)*z(0, 2) &
       + gn(:, 6)*z(0, 1)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*n(2)*z(0, 3) &
       + n(1)*z(0, 2)*n(3) &
       + z(0, 1)*n(2)*n(3)
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/3, 0, 1, 0, 0/)))*tu
    tu = n(4)*z(0, 3) &
       + n(5)*z(0, 2) &
       + n(6)*z(0, 1) &
       + z(0, 4)*n(3) &
       + z(1, 2)*n(3) &
       + z(0, 5)*n(2) &
       + z(1, 3)*n(2) &
       + z(0, 6)*n(1) &
       + z(2, 3)*n(1)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tu
    tu = z(0, 1)*z(0, 2)*z(0, 3)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 3, 0, 0/)))*tu
    tu = z(0, 4)*z(0, 3) &
       + z(1, 2)*z(0, 3) &
       + z(0, 5)*z(0, 2) &
       + z(1, 3)*z(0, 2) &
       + z(0, 6)*z(0, 1) &
       + z(2, 3)*z(0, 1)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tu
    tv = n(4)*gn(:, 3) &
       + n(5)*gn(:, 2) &
       + n(6)*gn(:, 1) &
       + gn(:, 4)*n(3) &
       + gn(:, 5)*n(2) &
       + gn(:, 6)*n(1)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tv(:)
    tv = n(1)*z(0, 2)*gn(:, 3) &
       + n(1)*gn(:, 2)*z(0, 3) &
       + z(0, 1)*n(2)*gn(:, 3) &
       + z(0, 1)*gn(:, 2)*n(3) &
       + gn(:, 1)*n(2)*z(0, 3) &
       + gn(:, 1)*z(0, 2)*n(3)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2)*z(0, 3) &
       + z(0, 1)*n(2)*z(0, 3) &
       + z(0, 1)*z(0, 2)*n(3)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 2, 0, 0/)))*tu
    tv = n(1)*n(2)*gn(:, 3) &
       + n(1)*gn(:, 2)*n(3) &
       + gn(:, 1)*n(2)*n(3)
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 6) &
       + z(2, 5) &
       + z(3, 4)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tu
    tv = z(0, 1)*z(0, 2)*gn(:, 3) &
       + z(0, 1)*gn(:, 2)*z(0, 3) &
       + gn(:, 1)*z(0, 2)*z(0, 3)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 3, 0, 0/)))*tv(:)
    end if

  end subroutine

  subroutine xc_response_5(u, v, derv, n, gn, z, use_alda)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_alda
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)*n(2)*n(3)*n(4)
    u  = u + derv(xc_index(XC_RSZXY, (/5, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/4, 0, 1, 0, 0/)))*tu*gn(:, 0)
    end if
    tu = n(11)*n(4) &
       + n(12)*n(3) &
       + n(13)*n(2) &
       + n(14)*n(1) &
       + n(5)*n(10) &
       + n(6)*n(9) &
       + n(7)*n(8)
    u  = u + derv(xc_index(XC_RSZXY, (/3, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    end if
    tu = n(5)*n(3)*n(4) &
       + n(6)*n(2)*n(4) &
       + n(7)*n(2)*n(3) &
       + n(8)*n(1)*n(4) &
       + n(9)*n(1)*n(3) &
       + n(10)*n(1)*n(2)
    u  = u + derv(xc_index(XC_RSZXY, (/4, 0, 0, 0, 0/)))*tu
    if (.not. use_alda) then
    v  = v + derv(xc_index(XC_RSZXY, (/3, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tv = z(0, 11)*gn(:, 4) &
       + z(1, 8)*gn(:, 4) &
       + z(2, 6)*gn(:, 4) &
       + z(3, 5)*gn(:, 4) &
       + z(0, 12)*gn(:, 3) &
       + z(1, 9)*gn(:, 3) &
       + z(2, 7)*gn(:, 3) &
       + z(4, 5)*gn(:, 3) &
       + z(0, 13)*gn(:, 2) &
       + z(1, 10)*gn(:, 2) &
       + z(3, 7)*gn(:, 2) &
       + z(4, 6)*gn(:, 2) &
       + z(0, 14)*gn(:, 1) &
       + z(2, 10)*gn(:, 1) &
       + z(3, 9)*gn(:, 1) &
       + z(4, 8)*gn(:, 1) &
       + gn(:, 11)*z(0, 4) &
       + gn(:, 12)*z(0, 3) &
       + gn(:, 13)*z(0, 2) &
       + gn(:, 14)*z(0, 1) &
       + z(0, 5)*gn(:, 10) &
       + z(1, 2)*gn(:, 10) &
       + z(0, 6)*gn(:, 9) &
       + z(1, 3)*gn(:, 9) &
       + z(0, 7)*gn(:, 8) &
       + z(1, 4)*gn(:, 8) &
       + z(0, 8)*gn(:, 7) &
       + z(2, 3)*gn(:, 7) &
       + z(0, 9)*gn(:, 6) &
       + z(2, 4)*gn(:, 6) &
       + z(0, 10)*gn(:, 5) &
       + z(3, 4)*gn(:, 5)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tv(:)
    tv = n(11)*gn(:, 4) &
       + n(12)*gn(:, 3) &
       + n(13)*gn(:, 2) &
       + n(14)*gn(:, 1) &
       + gn(:, 11)*n(4) &
       + gn(:, 12)*n(3) &
       + gn(:, 13)*n(2) &
       + gn(:, 14)*n(1) &
       + n(5)*gn(:, 10) &
       + n(6)*gn(:, 9) &
       + n(7)*gn(:, 8) &
       + n(8)*gn(:, 7) &
       + n(9)*gn(:, 6) &
       + n(10)*gn(:, 5)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2)*z(0, 3)*z(0, 4) &
       + z(0, 1)*n(2)*z(0, 3)*z(0, 4) &
       + z(0, 1)*z(0, 2)*n(3)*z(0, 4) &
       + z(0, 1)*z(0, 2)*z(0, 3)*n(4)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 3, 0, 0/)))*tu
    tv = z(0, 5)*z(0, 3)*gn(:, 4) &
       + z(1, 2)*z(0, 3)*gn(:, 4) &
       + z(0, 6)*z(0, 2)*gn(:, 4) &
       + z(1, 3)*z(0, 2)*gn(:, 4) &
       + z(0, 7)*z(0, 2)*gn(:, 3) &
       + z(1, 4)*z(0, 2)*gn(:, 3) &
       + z(0, 8)*z(0, 1)*gn(:, 4) &
       + z(2, 3)*z(0, 1)*gn(:, 4) &
       + z(0, 9)*z(0, 1)*gn(:, 3) &
       + z(2, 4)*z(0, 1)*gn(:, 3) &
       + z(0, 10)*z(0, 1)*gn(:, 2) &
       + z(3, 4)*z(0, 1)*gn(:, 2) &
       + z(0, 5)*gn(:, 3)*z(0, 4) &
       + z(1, 2)*gn(:, 3)*z(0, 4) &
       + z(0, 6)*gn(:, 2)*z(0, 4) &
       + z(1, 3)*gn(:, 2)*z(0, 4) &
       + z(0, 7)*gn(:, 2)*z(0, 3) &
       + z(1, 4)*gn(:, 2)*z(0, 3) &
       + z(0, 8)*gn(:, 1)*z(0, 4) &
       + z(2, 3)*gn(:, 1)*z(0, 4) &
       + z(0, 9)*gn(:, 1)*z(0, 3) &
       + z(2, 4)*gn(:, 1)*z(0, 3) &
       + z(0, 10)*gn(:, 1)*z(0, 2) &
       + z(3, 4)*gn(:, 1)*z(0, 2) &
       + gn(:, 5)*z(0, 3)*z(0, 4) &
       + gn(:, 6)*z(0, 2)*z(0, 4) &
       + gn(:, 7)*z(0, 2)*z(0, 3) &
       + gn(:, 8)*z(0, 1)*z(0, 4) &
       + gn(:, 9)*z(0, 1)*z(0, 3) &
       + gn(:, 10)*z(0, 1)*z(0, 2)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 3, 0, 0/)))*tv(:)
    tu = z(0, 1)*z(0, 2)*z(0, 3)*z(0, 4)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 5, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 4, 0, 0/)))*tu
    tu = z(0, 11)*z(0, 4) &
       + z(1, 8)*z(0, 4) &
       + z(2, 6)*z(0, 4) &
       + z(3, 5)*z(0, 4) &
       + z(0, 12)*z(0, 3) &
       + z(1, 9)*z(0, 3) &
       + z(2, 7)*z(0, 3) &
       + z(4, 5)*z(0, 3) &
       + z(0, 13)*z(0, 2) &
       + z(1, 10)*z(0, 2) &
       + z(3, 7)*z(0, 2) &
       + z(4, 6)*z(0, 2) &
       + z(0, 14)*z(0, 1) &
       + z(2, 10)*z(0, 1) &
       + z(3, 9)*z(0, 1) &
       + z(4, 8)*z(0, 1) &
       + z(0, 5)*z(0, 10) &
       + z(0, 5)*z(3, 4) &
       + z(1, 2)*z(0, 10) &
       + z(1, 2)*z(3, 4) &
       + z(0, 6)*z(0, 9) &
       + z(0, 6)*z(2, 4) &
       + z(1, 3)*z(0, 9) &
       + z(1, 3)*z(2, 4) &
       + z(0, 7)*z(0, 8) &
       + z(0, 7)*z(2, 3) &
       + z(1, 4)*z(0, 8) &
       + z(1, 4)*z(2, 3)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tu
    tv = n(5)*n(3)*gn(:, 4) &
       + n(6)*n(2)*gn(:, 4) &
       + n(7)*n(2)*gn(:, 3) &
       + n(8)*n(1)*gn(:, 4) &
       + n(9)*n(1)*gn(:, 3) &
       + n(10)*n(1)*gn(:, 2) &
       + n(5)*gn(:, 3)*n(4) &
       + n(6)*gn(:, 2)*n(4) &
       + n(7)*gn(:, 2)*n(3) &
       + n(8)*gn(:, 1)*n(4) &
       + n(9)*gn(:, 1)*n(3) &
       + n(10)*gn(:, 1)*n(2) &
       + gn(:, 5)*n(3)*n(4) &
       + gn(:, 6)*n(2)*n(4) &
       + gn(:, 7)*n(2)*n(3) &
       + gn(:, 8)*n(1)*n(4) &
       + gn(:, 9)*n(1)*n(3) &
       + gn(:, 10)*n(1)*n(2)
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tv(:)
    tu = n(5)*z(0, 3)*z(0, 4) &
       + n(6)*z(0, 2)*z(0, 4) &
       + n(7)*z(0, 2)*z(0, 3) &
       + n(8)*z(0, 1)*z(0, 4) &
       + n(9)*z(0, 1)*z(0, 3) &
       + n(10)*z(0, 1)*z(0, 2) &
       + z(0, 5)*n(3)*z(0, 4) &
       + z(1, 2)*n(3)*z(0, 4) &
       + z(0, 6)*n(2)*z(0, 4) &
       + z(1, 3)*n(2)*z(0, 4) &
       + z(0, 7)*n(2)*z(0, 3) &
       + z(1, 4)*n(2)*z(0, 3) &
       + z(0, 8)*n(1)*z(0, 4) &
       + z(2, 3)*n(1)*z(0, 4) &
       + z(0, 9)*n(1)*z(0, 3) &
       + z(2, 4)*n(1)*z(0, 3) &
       + z(0, 10)*n(1)*z(0, 2) &
       + z(3, 4)*n(1)*z(0, 2) &
       + z(0, 5)*z(0, 3)*n(4) &
       + z(1, 2)*z(0, 3)*n(4) &
       + z(0, 6)*z(0, 2)*n(4) &
       + z(1, 3)*z(0, 2)*n(4) &
       + z(0, 7)*z(0, 2)*n(3) &
       + z(1, 4)*z(0, 2)*n(3) &
       + z(0, 8)*z(0, 1)*n(4) &
       + z(2, 3)*z(0, 1)*n(4) &
       + z(0, 9)*z(0, 1)*n(3) &
       + z(2, 4)*z(0, 1)*n(3) &
       + z(0, 10)*z(0, 1)*n(2) &
       + z(3, 4)*z(0, 1)*n(2)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 2, 0, 0/)))*tu
    tv = n(1)*z(0, 2)*z(0, 3)*gn(:, 4) &
       + n(1)*z(0, 2)*gn(:, 3)*z(0, 4) &
       + n(1)*gn(:, 2)*z(0, 3)*z(0, 4) &
       + z(0, 1)*n(2)*z(0, 3)*gn(:, 4) &
       + z(0, 1)*n(2)*gn(:, 3)*z(0, 4) &
       + z(0, 1)*z(0, 2)*n(3)*gn(:, 4) &
       + z(0, 1)*z(0, 2)*gn(:, 3)*n(4) &
       + z(0, 1)*gn(:, 2)*n(3)*z(0, 4) &
       + z(0, 1)*gn(:, 2)*z(0, 3)*n(4) &
       + gn(:, 1)*n(2)*z(0, 3)*z(0, 4) &
       + gn(:, 1)*z(0, 2)*n(3)*z(0, 4) &
       + gn(:, 1)*z(0, 2)*z(0, 3)*n(4)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 3, 0, 0/)))*tv(:)
    tu = z(0, 5)*z(0, 3)*z(0, 4) &
       + z(1, 2)*z(0, 3)*z(0, 4) &
       + z(0, 6)*z(0, 2)*z(0, 4) &
       + z(1, 3)*z(0, 2)*z(0, 4) &
       + z(0, 7)*z(0, 2)*z(0, 3) &
       + z(1, 4)*z(0, 2)*z(0, 3) &
       + z(0, 8)*z(0, 1)*z(0, 4) &
       + z(2, 3)*z(0, 1)*z(0, 4) &
       + z(0, 9)*z(0, 1)*z(0, 3) &
       + z(2, 4)*z(0, 1)*z(0, 3) &
       + z(0, 10)*z(0, 1)*z(0, 2) &
       + z(3, 4)*z(0, 1)*z(0, 2)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 3, 0, 0/)))*tu
    tu = n(1)*n(2)*n(3)*z(0, 4) &
       + n(1)*n(2)*z(0, 3)*n(4) &
       + n(1)*z(0, 2)*n(3)*n(4) &
       + z(0, 1)*n(2)*n(3)*n(4)
    v  = v + derv(xc_index(XC_RSZXY, (/3, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/4, 0, 1, 0, 0/)))*tu
    tv = z(0, 1)*z(0, 2)*z(0, 3)*gn(:, 4) &
       + z(0, 1)*z(0, 2)*gn(:, 3)*z(0, 4) &
       + z(0, 1)*gn(:, 2)*z(0, 3)*z(0, 4) &
       + gn(:, 1)*z(0, 2)*z(0, 3)*z(0, 4)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 4, 0, 0/)))*tv(:)
    tv = n(5)*z(0, 3)*gn(:, 4) &
       + n(6)*z(0, 2)*gn(:, 4) &
       + n(7)*z(0, 2)*gn(:, 3) &
       + n(8)*z(0, 1)*gn(:, 4) &
       + n(9)*z(0, 1)*gn(:, 3) &
       + n(10)*z(0, 1)*gn(:, 2) &
       + n(5)*gn(:, 3)*z(0, 4) &
       + n(6)*gn(:, 2)*z(0, 4) &
       + n(7)*gn(:, 2)*z(0, 3) &
       + n(8)*gn(:, 1)*z(0, 4) &
       + n(9)*gn(:, 1)*z(0, 3) &
       + n(10)*gn(:, 1)*z(0, 2) &
       + z(0, 5)*n(3)*gn(:, 4) &
       + z(1, 2)*n(3)*gn(:, 4) &
       + z(0, 6)*n(2)*gn(:, 4) &
       + z(1, 3)*n(2)*gn(:, 4) &
       + z(0, 7)*n(2)*gn(:, 3) &
       + z(1, 4)*n(2)*gn(:, 3) &
       + z(0, 8)*n(1)*gn(:, 4) &
       + z(2, 3)*n(1)*gn(:, 4) &
       + z(0, 9)*n(1)*gn(:, 3) &
       + z(2, 4)*n(1)*gn(:, 3) &
       + z(0, 10)*n(1)*gn(:, 2) &
       + z(3, 4)*n(1)*gn(:, 2) &
       + z(0, 5)*gn(:, 3)*n(4) &
       + z(1, 2)*gn(:, 3)*n(4) &
       + z(0, 6)*gn(:, 2)*n(4) &
       + z(1, 3)*gn(:, 2)*n(4) &
       + z(0, 7)*gn(:, 2)*n(3) &
       + z(1, 4)*gn(:, 2)*n(3) &
       + z(0, 8)*gn(:, 1)*n(4) &
       + z(2, 3)*gn(:, 1)*n(4) &
       + z(0, 9)*gn(:, 1)*n(3) &
       + z(2, 4)*gn(:, 1)*n(3) &
       + z(0, 10)*gn(:, 1)*n(2) &
       + z(3, 4)*gn(:, 1)*n(2) &
       + gn(:, 5)*n(3)*z(0, 4) &
       + gn(:, 6)*n(2)*z(0, 4) &
       + gn(:, 7)*n(2)*z(0, 3) &
       + gn(:, 8)*n(1)*z(0, 4) &
       + gn(:, 9)*n(1)*z(0, 3) &
       + gn(:, 10)*n(1)*z(0, 2) &
       + gn(:, 5)*z(0, 3)*n(4) &
       + gn(:, 6)*z(0, 2)*n(4) &
       + gn(:, 7)*z(0, 2)*n(3) &
       + gn(:, 8)*z(0, 1)*n(4) &
       + gn(:, 9)*z(0, 1)*n(3) &
       + gn(:, 10)*z(0, 1)*n(2)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tv(:)
    tu = n(11)*z(0, 4) &
       + n(12)*z(0, 3) &
       + n(13)*z(0, 2) &
       + n(14)*z(0, 1) &
       + z(0, 11)*n(4) &
       + z(1, 8)*n(4) &
       + z(2, 6)*n(4) &
       + z(3, 5)*n(4) &
       + z(0, 12)*n(3) &
       + z(1, 9)*n(3) &
       + z(2, 7)*n(3) &
       + z(4, 5)*n(3) &
       + z(0, 13)*n(2) &
       + z(1, 10)*n(2) &
       + z(3, 7)*n(2) &
       + z(4, 6)*n(2) &
       + z(0, 14)*n(1) &
       + z(2, 10)*n(1) &
       + z(3, 9)*n(1) &
       + z(4, 8)*n(1) &
       + n(5)*z(0, 10) &
       + n(5)*z(3, 4) &
       + n(6)*z(0, 9) &
       + n(6)*z(2, 4) &
       + n(7)*z(0, 8) &
       + n(7)*z(2, 3) &
       + n(8)*z(0, 7) &
       + n(8)*z(1, 4) &
       + n(9)*z(0, 6) &
       + n(9)*z(1, 3) &
       + n(10)*z(0, 5) &
       + n(10)*z(1, 2)
    v  = v + derv(xc_index(XC_RSZXY, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/2, 0, 1, 0, 0/)))*tu
    tu = n(5)*n(3)*z(0, 4) &
       + n(6)*n(2)*z(0, 4) &
       + n(7)*n(2)*z(0, 3) &
       + n(8)*n(1)*z(0, 4) &
       + n(9)*n(1)*z(0, 3) &
       + n(10)*n(1)*z(0, 2) &
       + n(5)*z(0, 3)*n(4) &
       + n(6)*z(0, 2)*n(4) &
       + n(7)*z(0, 2)*n(3) &
       + n(8)*z(0, 1)*n(4) &
       + n(9)*z(0, 1)*n(3) &
       + n(10)*z(0, 1)*n(2) &
       + z(0, 5)*n(3)*n(4) &
       + z(1, 2)*n(3)*n(4) &
       + z(0, 6)*n(2)*n(4) &
       + z(1, 3)*n(2)*n(4) &
       + z(0, 7)*n(2)*n(3) &
       + z(1, 4)*n(2)*n(3) &
       + z(0, 8)*n(1)*n(4) &
       + z(2, 3)*n(1)*n(4) &
       + z(0, 9)*n(1)*n(3) &
       + z(2, 4)*n(1)*n(3) &
       + z(0, 10)*n(1)*n(2) &
       + z(3, 4)*n(1)*n(2)
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/3, 0, 1, 0, 0/)))*tu
    tv = n(1)*n(2)*z(0, 3)*gn(:, 4) &
       + n(1)*n(2)*gn(:, 3)*z(0, 4) &
       + n(1)*z(0, 2)*n(3)*gn(:, 4) &
       + n(1)*z(0, 2)*gn(:, 3)*n(4) &
       + n(1)*gn(:, 2)*n(3)*z(0, 4) &
       + n(1)*gn(:, 2)*z(0, 3)*n(4) &
       + z(0, 1)*n(2)*n(3)*gn(:, 4) &
       + z(0, 1)*n(2)*gn(:, 3)*n(4) &
       + z(0, 1)*gn(:, 2)*n(3)*n(4) &
       + gn(:, 1)*n(2)*n(3)*z(0, 4) &
       + gn(:, 1)*n(2)*z(0, 3)*n(4) &
       + gn(:, 1)*z(0, 2)*n(3)*n(4)
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*n(2)*z(0, 3)*z(0, 4) &
       + n(1)*z(0, 2)*n(3)*z(0, 4) &
       + n(1)*z(0, 2)*z(0, 3)*n(4) &
       + z(0, 1)*n(2)*n(3)*z(0, 4) &
       + z(0, 1)*n(2)*z(0, 3)*n(4) &
       + z(0, 1)*z(0, 2)*n(3)*n(4)
    v  = v + derv(xc_index(XC_RSZXY, (/2, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/3, 0, 2, 0, 0/)))*tu
    tv = n(1)*n(2)*n(3)*gn(:, 4) &
       + n(1)*n(2)*gn(:, 3)*n(4) &
       + n(1)*gn(:, 2)*n(3)*n(4) &
       + gn(:, 1)*n(2)*n(3)*n(4)
    v  = v + derv(xc_index(XC_RSZXY, (/3, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 14) &
       + z(2, 13) &
       + z(3, 12) &
       + z(4, 11) &
       + z(5, 10) &
       + z(6, 9) &
       + z(7, 8)
    v  = v + derv(xc_index(XC_RSZXY, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(XC_RSZXY, (/1, 0, 1, 0, 0/)))*tu
    end if

  end subroutine
#endif /* ifdef MOD_AUTODIFF */

end module
