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












contains

  subroutine xc_response_contribution(order,   &
                         nr_dmat, &
                         use_gga, &
                         dmat,    &
                         fmat,    &
                         isymop,  &
                         gao,     &
                         gao1,    &
                         r_0,     &
                         gr_0,    &
                         w,       &
                         buf)

#include "dgroup.h"
#include "dcbbas.h"
#include "../dft/functionals.h"

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: order
    integer, intent(in)    :: nr_dmat
    logical, intent(in)    :: use_gga
    real(8), intent(in)    :: dmat(ntbas(0), ntbas(0), nz, nr_dmat)
    real(8), intent(inout) :: fmat(ntbas(0), ntbas(0), nz)
    integer, intent(in)    :: isymop(nr_dmat)
    real(8), intent(in)    :: gao(nr_ao_total)
    real(8), intent(in)    :: gao1(nr_ao_total*3)
    real(8), intent(in)    :: r_0
    real(8), intent(in)    :: gr_0(3)
    real(8), intent(in)    :: w
    real(8), intent(in)    :: buf(*)
!   ----------------------------------------------------------------------------

    real(8)                :: u
!   real(8)                :: s_aux(max_order)
    real(8)                :: v(3)

    real(8)                :: derv(0:max_order, 0:max_order)
    real(8)                :: derv_i(999)
    real(8)                :: p
    integer                :: taylor_index

    real(8)                :: r(0:max_order)
    real(8)                :: gr(3, 0:max_order)
    real(8)                :: z(0:max_order, 0:max_order)

    integer                :: i, j, k, n
    integer                :: on, oz, sum_on_oz
    real(8)                :: ddot
!   ----------------------------------------------------------------------------

!   get perturbed densities

    r(0) = r_0

    do i = 1, nr_dmat
      call get_n(r(i),             &
                 isymop(i) - 1,    &
                 dmat(1, 1, 1, i), &
                 buf,              &
                 gao)
    end do

!   get perturbed densitiy gradients and gradient norms

    if (use_gga) then

      gr(:, 0) = gr_0

      do i = 1, nr_dmat

        call get_gn(gr(1, i),         &
                    isymop(i) - 1,    &
                    dmat(1, 1, 1, i), &
                    buf,              &
                    gao,              &
                    gao1)
      end do

      do i = 0, nr_dmat
!       do j = 0, nr_dmat
        do j = i, nr_dmat
!         the perturbed z carries factor 2
!         because derivative of z gives 2 identical terms
!         (one "left", one "right")
          z(i, j) = 2.0d0*ddot(3, gr(1, i), 1, gr(1, j), 1)
        end do
      end do
!     correct wrong prefactor 2 of the unperturbed z
      z(0, 0) = 0.5d0*z(0, 0)

    end if

!   get derivatives and put them in derv_i (intermediate)

    derv_i = 0.0d0

    call eval_gga(derv_i, order, (/r(0), 0.0d0, z(0, 0), 0.0d0, 0.0d0/))
 
!   multiply by the weight

    derv_i = derv_i*w

!   distribute derivatives into d

      do on = 0, order
        do oz = 0, order
          sum_on_oz = on + oz
          if ((sum_on_oz > 0) .and. (sum_on_oz <= order)) then

!           this order zeroes them out
            derv(on, oz) = derv_i(taylor_index(5, (/on, 0, oz, 0, 0/)))


          end if
        end do
      end do

    u = 0.0d0
    v = 0.0d0

    select case (order)

      case (2)
        call xc_response_1r_n(u, derv_i, r)

      case default
        call quit('xc_response_contribution: order not implemented')

    end select

    if (use_gga) then
      select case (order)
      
        case (2)
          call xc_response_1r_nz(u, v, derv_i, r, gr, z)
      
        case default
          call quit('xc_response_contribution: order not implemented')
      
      end select
    end if

    if (use_gga) then
      call gga_distribute(u, 2.0d0*v, 0, fmat, gao, gao1)
    else
      call lda_distribute(u,          0, fmat, gao)
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
