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
    integer, parameter     :: max_order = 10



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
!#include "../dft/functionals.h"

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

    real(8)                :: n(0:max_order)
    real(8)                :: gn(3, 0:max_order)
    real(8)                :: z(0:max_order, 0:max_order)

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

    derv = 0.0d0
    call eval_gga(derv, order, (/n(0), 0.0d0, z(0, 0), 0.0d0, 0.0d0/))
 
!   multiply by the weight
    derv = derv*w

    u = 0.0d0
    v = 0.0d0

    if (calculate_blocks_n) then
      select case (order)
    
        case (2)
          call xc_response_1r_n(u, derv, n)
    
        case default
          call quit('xc_response_contribution: order not implemented')
    
      end select
    end if

    if (calculate_blocks_nz) then
      select case (order)
    
        case (2)
          call xc_response_1r_nz(u, v, derv, n, gn, z)
    
        case default
          call quit('xc_response_contribution: order not implemented')
    
      end select
    end if

    if (dabs(maxval(v)) > 0.0d0) then
      call gga_distribute(u, 2.0d0*v, 0, fmat, ao_g0_m0, ao_g1_m0)
    else
      call lda_distribute(u,          0, fmat, ao_g0_m0)
    end if

  end subroutine

  subroutine xc_response_1r_n(u, derv, n)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(in)  :: derv(*)
    real(8), intent(in)  :: n(0:max_order)
!   ----------------------------------------------------------------------------
    integer, external    :: taylor_index
!   ----------------------------------------------------------------------------

    u = u + derv(taylor_index(5, (/2, 0, 0, 0, 0/)))*n(1)

  end subroutine

  subroutine xc_response_1r_nz(u, v, derv, n, gn, z)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)  :: derv(*)
    real(8), intent(in)  :: n(0:max_order)
    real(8), intent(in)  :: gn(3, 0:max_order)
    real(8), intent(in)  :: z(0:max_order, 0:max_order)
!   ----------------------------------------------------------------------------
    integer, external    :: taylor_index
!   ----------------------------------------------------------------------------

    u = u + derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*z(0, 1)

    v = v + derv(taylor_index(5, (/0, 0, 2, 0, 0/)))*gn(:, 0)*z(0, 1)
    v = v + derv(taylor_index(5, (/1, 0, 1, 0, 0/)))*gn(:, 0)*n(1)
    v = v + derv(taylor_index(5, (/0, 0, 1, 0, 0/)))*gn(:, 1)

  end subroutine

#endif /* ifdef MOD_AUTODIFF */

end module
