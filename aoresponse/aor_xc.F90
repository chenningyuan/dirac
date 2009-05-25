module aor_xc

  use dft_cfg
  use dft_density

  implicit none

contains

  subroutine aor_xc_empty()
  end subroutine

#ifdef MOD_AUTODIFF
  subroutine aor_xc_new(order,   &
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
#include "../autodiff/xc_fortran_interface.h"

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: order
    integer, intent(in)    :: nr_dmat
    logical, intent(in)    :: use_gga
    real(8), intent(in)    :: dmat(ntbas(0), ntbas(0), nz, nr_dmat)
    real(8), intent(inout) :: fmat(ntbas(0), ntbas(0), nz)
    integer, intent(in)    :: isymop(nr_dmat)
    real(8), intent(in)    :: gao(*)
    real(8), intent(in)    :: gao1(*)
    real(8), intent(in)    :: r_0
    real(8), intent(in)    :: gr_0(3)
    real(8), intent(in)    :: w
    real(8), intent(in)    :: buf(*)
!   ----------------------------------------------------------------------------
    integer, parameter     :: max_order = 10

    real(8)                :: s
    real(8)                :: v(3)

    real(8)                :: derv(max_order, max_order)
    real(8)                :: derv_i(999, 3)
!#define CHECK_DERV
#ifdef CHECK_DERV
    real(8)                :: f, ff
    real(8)                :: derv_ip(999, 3)
    real(8)                :: derv_im(999, 3)
    real(8)                :: derv_p(max_order, max_order)
    real(8)                :: derv_m(max_order, max_order)
#endif
    real(8)                :: p, taylor_deriv

    real(8)                :: r(0:max_order)
    real(8)                :: gr(3, 0:max_order)
    real(8)                :: z(0:max_order, 0:max_order)

    integer                :: i, j
    integer                :: on, oz, sum_on_oz
    real(8)                :: ddot
!   ----------------------------------------------------------------------------

!   get perturbed densities

    r(0) = r_0
#ifdef CHECK_DERV
    f = r(0)/1.0d5
#endif

    do i = 1, nr_dmat
      call get_n(r(i),             &
                 0,                &
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
#ifdef CHECK_DERV
    derv_ip = 0.0d0
    derv_im = 0.0d0
#endif

    if (use_gga) then

      call lda_eval(LDA_SLATER, order, r(0),          derv_i(1, 1))
      call gga_eval(GGA_BECKE,  order, r(0), z(0, 0), derv_i(1, 2))
      call gga_eval(GGA_LYP,    order, r(0), z(0, 0), derv_i(1, 3))

    else

      call lda_eval(LDA_SLATER, order, r(0), derv_i(1, 1))
      call lda_eval(LDA_VWN,    order, r(0), derv_i(1, 2))
#ifdef CHECK_DERV
      call lda_eval(LDA_SLATER, order, r(0)+f, derv_ip(1, 1))
      call lda_eval(LDA_VWN,    order, r(0)+f, derv_ip(1, 2))
      call lda_eval(LDA_SLATER, order, r(0)-f, derv_im(1, 1))
      call lda_eval(LDA_VWN,    order, r(0)-f, derv_im(1, 2))
#endif

    end if

!   multiply by the weight

    derv_i = derv_i*w
#ifdef CHECK_DERV
    derv_ip = derv_ip*w
    derv_im = derv_im*w
#endif

!   distribute derivatives into d

    if (use_gga) then

      do on = 0, order
        do oz = 0, order
          sum_on_oz = on + oz
          if ((sum_on_oz > 0) .and. (sum_on_oz <= order)) then

!           this order zeroes them out
            derv(on, oz) = taylor_deriv(derv_i(1, 2), 5, (/on, 0, oz, 0, 0/)) &
                         + taylor_deriv(derv_i(1, 3), 5, (/on, 0, oz, 0, 0/))

            if (oz == 0) then
              derv(on, oz) = derv(on, oz) &
                           + taylor_deriv(derv_i(1, 1), 2, (/on, 0/))
            end if

          end if
        end do
      end do

    else

      do on = 1, order
        derv(on, 0) = taylor_deriv(derv_i(1, 1), 2, (/on, 0/)) &
                    + taylor_deriv(derv_i(1, 2), 2, (/on, 0/))
#ifdef CHECK_DERV
        derv_p(on, 0) = taylor_deriv(derv_ip(1, 1), 2, (/on, 0/)) &
                      + taylor_deriv(derv_ip(1, 2), 2, (/on, 0/))
        derv_m(on, 0) = taylor_deriv(derv_im(1, 1), 2, (/on, 0/)) &
                      + taylor_deriv(derv_im(1, 2), 2, (/on, 0/))

        if (on > 1) then
          ff = (derv_p(on-1, 0) - derv_m(on-1, 0))/(2.0d0*f)
          write(*, '(a7, i2, 3e14.6)') 'compare',   &
                                       on,          &
                                       r(0),        &
                                       derv(on, 0), &
                                       ff
        end if
#endif
      end do

    end if

    s = 0.0d0
    v = 0.0d0


!   gradient contribution
!   =====================

    if (use_gga) then

      if      (order == 2) then

        s = s + derv(1, 1)*z(0, 1)

        v = v + gr(:, 0)*(                   &
                          derv(0, 2)*z(0, 1) &
                        + derv(1, 1)*r(1)    &
                         )
        v = v + gr(:, 1)*(           &
                          derv(0, 1) &
                         )

      else if (order == 3) then

        s = s + 2.0d0*derv(2, 1)*r(1)*z(0, 2)    &
              + 2.0d0*derv(2, 1)*r(2)*z(0, 1)    &
              + 2.0d0*derv(1, 2)*z(0, 1)*z(0, 2) &
              + 2.0d0*derv(1, 1)*z(1, 2)

        v = v + gr(:, 0)*(                           &
                          derv(0, 3)*z(0, 1)*z(0, 2) &
                        + derv(1, 2)*r(1)*z(0, 2)    &
                        + derv(1, 2)*r(2)*z(0, 1)    &
                        + derv(2, 1)*r(1)*r(2)       &
                        + derv(0, 2)*z(1, 2)         &
                         )
        v = v + gr(:, 1)*(                   &
                          derv(1, 1)*r(2)    &
                        + derv(0, 2)*z(0, 2) &
                         )
        v = v + gr(:, 2)*(                   &
                          derv(1, 1)*r(1)    &
                        + derv(0, 2)*z(0, 1) &
                         )

!       p (1, 2)
        v = 2.0d0*v

      else if (order == 4) then

!       fix prefactors
!       s = s +       derv(3, 1)*z(0, 1)*r(2)*r(3)
!                     derv(3, 1)*z(0, 2)*r(1)*r(3)
!                     derv(3, 1)*z(0, 3)*r(1)*r(2)
!
!                     derv(2, 1)*z(1, 2)*r(3)
!                     derv(2, 1)*z(1, 3)*r(2)
!                     derv(2, 1)*z(2, 3)*r(1)
!
!                     derv(0, 1)*z(1, 6)
!                     derv(0, 1)*z(2, 5)
!                     derv(0, 1)*z(3, 4)
!
!                     derv(2, 1)*z(0, 1)*r(6)
!                     derv(2, 1)*z(0, 2)*r(5)
!                     derv(2, 1)*z(0, 3)*r(4)
!                     derv(2, 1)*z(0, 4)*r(3)
!                     derv(2, 1)*z(0, 5)*r(2)
!                     derv(2, 1)*z(0, 6)*r(1)

      end if

    end if


!   pure density contribution
!   =========================

!   premultiply derivatives
    do i = 3, order
      derv(i, 0) = derv(i, 0)*2.0d0**(i - 2)
    end do

    if      (order == 2) then

      s = s + derv(2, 0)*r(1)

    else if (order == 3) then

!     the r(1)*r(2) here may be confusing
!     it's because from order 4 on i restrict
!     perturbations to be identical
!     here in order 3 it's still explicit
!     here both r(1) and r(2) are first order densities
      s = s + derv(3, 0)*r(1)*r(2)

    else if (order == 4) then

!     from here on r(1) is a first-order density
!     and r(2) is a second-order density etc.
      s = s +       derv(4, 0)*r(1)*r(1)*r(1) &
            + 3.0d0*derv(3, 0)*r(1)*r(2) !i jk, j ik, k ij

    else if (order == 5) then

      s = s +       derv(5, 0)*r(1)*r(1)*r(1)*r(1) &
            + 6.0d0*derv(4, 0)*r(1)*r(1)*r(2)      & !i j kl, i k jl, i l jk, j k il, j l ik, k l ij
            + 4.0d0*derv(3, 0)*r(1)*r(3)           & !i jkl, j ikl, k ijl, l ijk
            + 3.0d0*derv(3, 0)*r(2)*r(2)             !ij kl, ik jl, il jk

    else if (order == 6) then

!     4 1
!     3 1 1
!     3 2
!     2 1 1 1
!     2 2 1
!     2 3
!     1 1 1 1 1
!     1 2 1 1
!     1 3 1
!     1 4

      s = s +        derv(6, 0)*r(1)*r(1)*r(1)*r(1)*r(1) &
            + 10.0d0*derv(5, 0)*r(1)*r(1)*r(1)*r(2)      &
            + 10.0d0*derv(4, 0)*r(1)*r(1)*r(3)           &
            +  5.0d0*derv(3, 0)*r(1)*r(4)                &
            + 15.0d0*derv(4, 0)*r(1)*r(2)*r(2)           &
            + 10.0d0*derv(3, 0)*r(2)*r(3)

    else if (order == 7) then

!     5 1
!     4 1 1
!     4 2
!     3 1 1 1
!     3 2 1
!     3 3
!     2 1 1 1 1
!     2 2 1 1
!     2 3 1
!     2 4
!     1 1 1 1 1 1
!     1 2 1 1 1
!     1 3 1 1
!     1 4 1
!     1 5

      s = s +        derv(7, 0)*r(1)*r(1)*r(1)*r(1)*r(1)*r(1) & !i j k l m n
            + 15.0d0*derv(6, 0)*r(1)*r(1)*r(1)*r(1)*r(2)      & !i j k l mn, ...
            + 20.0d0*derv(5, 0)*r(1)*r(1)*r(1)*r(3)           & !i j k lmn, ...
            + 15.0d0*derv(4, 0)*r(1)*r(1)*r(4)                & !i j klmn, ...
            +  6.0d0*derv(3, 0)*r(1)*r(5)                     & !i jklmn, ...
            + 30.0d0*derv(5, 0)*r(1)*r(1)*r(2)*r(2)           & !i j kl mn, ...
            + 60.0d0*derv(4, 0)*r(1)*r(2)*r(3)                & !i jk lmn, ...
            + 15.0d0*derv(3, 0)*r(2)*r(4)                     & !ij klmn, ...
            + 10.0d0*derv(3, 0)*r(3)*r(3)                       !ijk lmn, ...

    end if

    if (use_gga) then
      call gga_distribute(s, 2.0d0*v, 0, 0, fmat, gao, gao1)
    else
      call lda_distribute(s,          0, 0, fmat, gao)
    end if

  end subroutine
#endif /* ifdef MOD_AUTODIFF */

end module
