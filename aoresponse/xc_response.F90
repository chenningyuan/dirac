module xc_response

#ifdef MOD_AUTODIFF
  use dft_cfg
  use dft_density
  use aor_binomial
  use dft_cb_interface
  use xc_fun

  implicit none

  public xc_response_contribution

  private
    integer, parameter     :: max_nr_dmat = 20



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

#define SCALE_DENSITIES
#ifdef SCALE_DENSITIES
!     perturbed matrices are scaled by 2 in dftexc
!     correct for it here
      n(i) = 0.5d0*n(i)
#endif
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

#ifdef SCALE_DENSITIES
!       perturbed matrices are scaled by 2 in dftexc
!       correct for it here
        gn(:, i) = 0.5d0*gn(:, i)
#endif
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

    derv = 0.0d0
    call eval_gga(derv, order, (/n(0), 0.0d0, z(0, 0), 0.0d0, 0.0d0/))
 
!   multiply by the weight
    derv = derv*w

    u = 0.0d0
    v = 0.0d0

    select case (order)
    
      case (2)
        call xc_response_1(u, v, derv, n, gn, z)
 
      case (3)
        call xc_response_2(u, v, derv, n, gn, z)
 
      case (4)
        call xc_response_3(u, v, derv, n, gn, z)
 
!     case (5)
!       call xc_response_4(u, v, derv, n, gn, z)
    
      case default
        call quit('xc_response_contribution: order not implemented')
    
    end select

    u = u*2.0d0**(order - 2)
    v = v*2.0d0**(order - 2)

    if (maxval((/dabs(v(1)), dabs(v(2)), dabs(v(3))/)) > 0.0d0) then
      call gga_distribute(u, 2.0d0*v, 0, fmat, ao_g0_m0, ao_g1_m0)
    else
      call lda_distribute(u,          0, fmat, ao_g0_m0)
    end if

  end subroutine

  subroutine xc_response_1(u, v, derv, n, gn, z) 
    
!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
!   ----------------------------------------------------------------------------
    integer, external      :: taylor_index
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)
    v  = v + derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/2, 0, 0, 0, 0/)))*tu
    tu = z(0, 1)
    v  = v + derv(taylor_index(5, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*tu
    tv = gn(:, 1)
    v  = v + derv(taylor_index(5, (/0, 0, 1, 0, 0/)))*tv(:)
    
  end subroutine

  subroutine xc_response_2(u, v, derv, n, gn, z) 
    
!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
!   ----------------------------------------------------------------------------
    integer, external      :: taylor_index
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tv = z(0, 1)*gn(:, 2) &
       + gn(:, 1)*z(0, 2)
    v  = v + derv(taylor_index(5, (/0, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2) &
       + z(0, 1)*n(2)
    v  = v + derv(taylor_index(5, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/2, 0, 1, 0, 0/)))*tu
    tu = z(0, 1)*z(0, 2)
    v  = v + derv(taylor_index(5, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/1, 0, 2, 0, 0/)))*tu
    tv = n(1)*gn(:, 2) &
       + gn(:, 1)*n(2)
    v  = v + derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 2)
    v  = v + derv(taylor_index(5, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*tu
    tu = n(1)*n(2)
    v  = v + derv(taylor_index(5, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/3, 0, 0, 0, 0/)))*tu
    
  end subroutine

  subroutine xc_response_3(u, v, derv, n, gn, z) 
    
!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
!   ----------------------------------------------------------------------------
    integer, external      :: taylor_index
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tv = (z(0, 4) + 2.0d0*z(1, 2))*gn(:, 3) &
       + (z(0, 5) + 2.0d0*z(1, 3))*gn(:, 2) &
       + (z(0, 6) + 2.0d0*z(2, 3))*gn(:, 1) &
       + gn(:, 4)*z(0, 3) &
       + gn(:, 5)*z(0, 2) &
       + gn(:, 6)*z(0, 1)
    v  = v + 0.5d0*derv(taylor_index(5, (/0, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*n(2)*z(0, 3) &
       + n(1)*z(0, 2)*n(3) &
       + z(0, 1)*n(2)*n(3)
    v  = v + derv(taylor_index(5, (/2, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/3, 0, 1, 0, 0/)))*tu
    tu = n(4)*z(0, 3) &
       + n(5)*z(0, 2) &
       + n(6)*z(0, 1) &
       + (z(0, 4) + 2.0d0*z(1, 2))*n(3) &
       + (z(0, 5) + 2.0d0*z(1, 3))*n(2) &
       + (z(0, 6) + 2.0d0*z(2, 3))*n(1)
    v  = v + 0.5d0*derv(taylor_index(5, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + 0.5d0*derv(taylor_index(5, (/2, 0, 1, 0, 0/)))*tu
    tu = z(0, 1)*z(0, 2)*z(0, 3)
    v  = v + derv(taylor_index(5, (/0, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/1, 0, 3, 0, 0/)))*tu
    tu = (z(0, 4) + 2.0d0*z(1, 2))*z(0, 3) &
       + (z(0, 5) + 2.0d0*z(1, 3))*z(0, 2) &
       + (z(0, 6) + 2.0d0*z(2, 3))*z(0, 1)
    v  = v + 0.5d0*derv(taylor_index(5, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + 0.5d0*derv(taylor_index(5, (/1, 0, 2, 0, 0/)))*tu
    tv = n(4)*gn(:, 3) &
       + n(5)*gn(:, 2) &
       + n(6)*gn(:, 1) &
       + gn(:, 4)*n(3) &
       + gn(:, 5)*n(2) &
       + gn(:, 6)*n(1)
    v  = v + 0.5d0*derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*tv(:)
    tv = n(1)*z(0, 2)*gn(:, 3) &
       + n(1)*gn(:, 2)*z(0, 3) &
       + z(0, 1)*n(2)*gn(:, 3) &
       + z(0, 1)*gn(:, 2)*n(3) &
       + gn(:, 1)*n(2)*z(0, 3) &
       + gn(:, 1)*z(0, 2)*n(3)
    v  = v + derv(taylor_index(5, (/1, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*n(2)*n(3)
    v  = v + derv(taylor_index(5, (/3, 0, 1, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/4, 0, 0, 0, 0/)))*tu
    tu = n(1)*z(0, 2)*z(0, 3) &
       + z(0, 1)*n(2)*z(0, 3) &
       + z(0, 1)*z(0, 2)*n(3)
    v  = v + derv(taylor_index(5, (/1, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(taylor_index(5, (/2, 0, 2, 0, 0/)))*tu
    tv = n(1)*n(2)*gn(:, 3) &
       + n(1)*gn(:, 2)*n(3) &
       + gn(:, 1)*n(2)*n(3)
    v  = v + derv(taylor_index(5, (/2, 0, 1, 0, 0/)))*tv(:)
    tu = (z(1, 6) + z(2, 5) + z(3, 4))
    v  = v + 0.5d0*derv(taylor_index(5, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + 0.5d0*derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*tu
    tv = z(0, 1)*z(0, 2)*gn(:, 3) &
       + z(0, 1)*gn(:, 2)*z(0, 3) &
       + gn(:, 1)*z(0, 2)*z(0, 3)
    v  = v + derv(taylor_index(5, (/0, 0, 3, 0, 0/)))*tv(:)
    tu = n(4)*n(3) &
       + n(5)*n(2) &
       + n(6)*n(1)
    v  = v + 0.5d0*derv(taylor_index(5, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    u  = u + 0.5d0*derv(taylor_index(5, (/3, 0, 0, 0, 0/)))*tu
    
  end subroutine
#endif /* ifdef MOD_AUTODIFF */

end module
