module aor_xc

#ifdef MOD_AUTODIFF
  use dft_cfg
  use dft_density
  use aor_binomial
  use dft_cb_interface
  use xc_fun

  implicit none

  public aor_xc_new
  public aor_xc_new_save_old
  public init_xc_contributions

  private

! interface eval_expansion
!   module procedure eval_n_expansion
!   module procedure eval_z_expansion
! end interface

  type xc_contribution
    real(8) :: prefactor                         = 1.0d0
    integer :: functional_derv_orders(5)         = 0
    integer :: nr_densities(5)                   = 0
    integer :: s_use_matrix(1:2, XC_MAX_ORDER)   = 0
    integer :: p_use_matrix_l(3:5, XC_MAX_ORDER) = 0
    integer :: p_use_matrix_r(3:5, XC_MAX_ORDER) = 0
    integer :: gradient_matrix                   = 0
    integer :: next                              = 0
  end type

  integer, parameter :: max_nr_of_xc_contributions = 1000

  save

  integer               :: nr_xc_ur
  integer               :: nr_xc_vr

  integer               :: istart_xc_ur
  integer               :: istart_xc_vr

  type(xc_contribution) :: xc_ur(max_nr_of_xc_contributions)
  type(xc_contribution) :: xc_vr(max_nr_of_xc_contributions)

  logical :: debug_me = .true.

contains

  subroutine init_xc_contributions(order)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: order
!   ----------------------------------------------------------------------------
    integer             :: i, j, ilast_accepted
!   ----------------------------------------------------------------------------

    nr_xc_ur = 0
    nr_xc_vr = 0

    select case (order)

    case (2)

      i = 0

      i = i + 1
      xc_ur(i)%nr_densities(1) = 1
      xc_ur(i)%s_use_matrix(1, 1) = 1
      xc_ur(i)%next = i + 1        
      nr_xc_ur      = nr_xc_ur + 1

      i = i + 1
      xc_ur(i)%nr_densities(3)   = 1
      xc_ur(i)%p_use_matrix_l(3, 1) = 0
      xc_ur(i)%p_use_matrix_r(3, 1) = 1
      xc_ur(i)%next = i + 1        
      nr_xc_ur      = nr_xc_ur + 1

      xc_ur(i)%next = 0

      i = 0

      i = i + 1
      xc_vr(i)%nr_densities(1) = 1
      xc_vr(i)%s_use_matrix(1, 1) = 1
      xc_vr(i)%gradient_matrix = 0
      xc_vr(i)%next = i + 1        
      nr_xc_vr      = nr_xc_vr + 1

      i = i + 1
      xc_vr(i)%nr_densities(3) = 1
      xc_vr(i)%p_use_matrix_l(3, 1) = 0
      xc_vr(i)%p_use_matrix_r(3, 1) = 1
      xc_vr(i)%gradient_matrix = 0
      xc_vr(i)%next = i + 1        
      nr_xc_vr      = nr_xc_vr + 1

      i = i + 1
      xc_vr(i)%gradient_matrix = 1
      xc_vr(i)%next = i + 1        
      nr_xc_vr      = nr_xc_vr + 1

      xc_vr(i)%next = 0

!   case (3)

!     i = 0

!     i = i + 1
!     xc_ur(i)%prefactor = 2.0d0
!     xc_ur(i)%nr_n                 = 2
!     xc_ur(i)%n_to_matrix(1)       = 1
!     xc_ur(i)%n_to_matrix(2)       = 2
!     xc_ur(i)%next                 = i + 1

!     i = i + 1
!     xc_ur(i)%prefactor = 2.0d0
!     xc_ur(i)%nr_n                 = 1
!     xc_ur(i)%nr_z                 = 1
!     xc_ur(i)%n_to_matrix(1)       = 1
!     xc_ur(i)%z_to_matrix_left(1)  = 0
!     xc_ur(i)%z_to_matrix_right(1) = 2
!     xc_ur(i)%next                 = i + 1

!     i = i + 1
!     xc_ur(i)%prefactor = 2.0d0
!     xc_ur(i)%nr_n                 = 1
!     xc_ur(i)%nr_z                 = 1
!     xc_ur(i)%n_to_matrix(1)       = 2
!     xc_ur(i)%z_to_matrix_left(1)  = 0
!     xc_ur(i)%z_to_matrix_right(1) = 1
!     xc_ur(i)%next                 = i + 1

!     i = i + 1
!     xc_ur(i)%prefactor = 2.0d0
!     xc_ur(i)%nr_z                 = 2
!     xc_ur(i)%z_to_matrix_left(1)  = 0
!     xc_ur(i)%z_to_matrix_right(1) = 1
!     xc_ur(i)%z_to_matrix_left(2)  = 0
!     xc_ur(i)%z_to_matrix_right(2) = 2
!     xc_ur(i)%next                 = i + 1

!     i = i + 1
!     xc_ur(i)%prefactor = 2.0d0
!     xc_ur(i)%nr_z                 = 2
!     xc_ur(i)%z_to_matrix_left(1)  = 1
!     xc_ur(i)%z_to_matrix_right(1) = 2
!     xc_ur(i)%z_to_matrix_left(2)  = 2
!     xc_ur(i)%z_to_matrix_right(2) = 1
!     xc_ur(i)%next                 = i + 1

!     xc_ur(i)%next = 0

    end select

    ilast_accepted = 0
    istart_xc_ur   = 0
    i = 1
    if (nr_xc_ur > 0) then
      do
        xc_ur(i)%functional_derv_orders = (/1, 0, 0, 0, 0/)
        do j = 1, 5
          xc_ur(i)%functional_derv_orders(j) = xc_ur(i)%functional_derv_orders(j) &
                                             + xc_ur(i)%nr_densities(j)
        end do

!       term is accepted
        if (istart_xc_ur == 0) istart_xc_ur = i
        ilast_accepted = i
!       if (debug_me) then
!         write(*, *) 'blub accepted ur', xc_ur(i)%functional_derv_orders
!       end if

!       term is not accepted
        if (ilast_accepted /= 0) then
!         xc_ur(ilast_accepted)%next = xc_ur(i)%next
        end if

        i = xc_ur(i)%next
        if (i == 0) exit
      end do
    end if

    istart_xc_vr = 1
    i = 1
    if (nr_xc_vr > 0) then
      do
        xc_vr(i)%functional_derv_orders = (/0, 0, 1, 0, 0/)
        do j = 1, 5
          xc_vr(i)%functional_derv_orders(j) = xc_vr(i)%functional_derv_orders(j) &
                                             + xc_vr(i)%nr_densities(j)
        end do

!       term is accepted
        if (istart_xc_vr /= 1) istart_xc_vr = i
        ilast_accepted = i

!       term is not accepted
!       xc_vr(ilast_accepted)%next = xc_vr(i)%next

        i = xc_vr(i)%next
        if (i == 0) exit
      end do
    end if

  end subroutine

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
    integer, parameter     :: max_order = 10

    real(8)                :: t
    real(8)                :: ur
    real(8)                :: vr(3)
    real(8)                :: derv(99)

    real(8)                :: n(0:max_order)
    real(8)                :: z(0:max_order, 0:max_order)
    real(8)                :: gn(3, 0:max_order)
    real(8)                :: z_0


    integer                :: i, j, in, iz, iv
    real(8), external      :: ddot
    integer, external      :: taylor_index
!   ----------------------------------------------------------------------------

!   get perturbed densities

    n(0) = r_0

    do i = 1, nr_dmat
      call get_n(n(i), &
                 isymop(i) - 1,             &
                 dmat(1, 1, 1, i),          &
                 buf,                       &
                 gao)
    end do

!   get perturbed densitiy gradients and gradient norms

    if (use_gga) then

      gn(:, 0) = gr_0

      do i = 1, nr_dmat

        call get_gn(gn(1, i),         &
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
          z(i, j) = 2.0d0*ddot(3, gn(1, i), 1, gn(1, j), 1)
        end do
      end do
!     correct wrong prefactor 2 of the unperturbed z
      z(0, 0) = 0.5d0*z(0, 0)

    end if


!   get functional derivatives
!   ==========================

!   if (use_gga) then
      call eval_gga(derv, order, (/r_0, 0.0d0, z(0, 0), 0.0d0, 0.0d0/))
!   else
!   end if


!   get ur contribution
!   ===================

    ur = 0.0d0

    if (nr_xc_ur > 0) then
      i = 1
      do
        t = w*derv(taylor_index(5, xc_ur(i)%functional_derv_orders))
     
        iv = 1
        do j = 1, xc_ur(i)%nr_densities(iv)
          t = t*n(xc_ur(i)%s_use_matrix(iv, j))
        end do
     
        iv = 3
        do j = 1, xc_ur(i)%nr_densities(iv)
          t = t*z(xc_ur(i)%p_use_matrix_l(iv, j), xc_ur(i)%p_use_matrix_r(iv, j))
        end do
     
        ur = ur + t
     
        i = xc_ur(i)%next
        if (i == 0) exit
      end do
    end if


!   get vr contribution
!   ===================

    vr = 0.0d0

    if (nr_xc_vr > 0) then
      i = 1
      do
        t = w*derv(taylor_index(5, xc_vr(i)%functional_derv_orders))
     
        iv = 1
        do j = 1, xc_vr(i)%nr_densities(iv)
          t = t*n(xc_vr(i)%s_use_matrix(iv, j))
        end do
     
        iv = 3
        do j = 1, xc_vr(i)%nr_densities(iv)
          t = t*z(xc_vr(i)%p_use_matrix_l(iv, j), xc_vr(i)%p_use_matrix_r(iv, j))
        end do
     
        vr = vr + t*gn(:, xc_vr(i)%gradient_matrix)
     
        i = xc_vr(i)%next
        if (i == 0) exit
      end do
    end if

    if (use_gga) then
      call gga_distribute(ur, 2.0d0*vr, 0, fmat, gao, gao1)
    else
      call lda_distribute(ur,           0, fmat, gao)
    end if

  end subroutine

  subroutine aor_xc_new_save_old(order,   &
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
    integer, parameter     :: max_order = 10

    real(8)                :: s
    real(8)                :: s_aux(max_order)
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
     write(*, *) 'blub derv me', taylor_index(5, (/1, 0, 1, 0, 0/)), &
       derv_i(taylor_index(5, (/1, 0, 1, 0, 0/))), r_0, z(0, 0)
 
!fixme
!   if (use_gga) then

!     call lda_eval(LDA_SLATER, order, r(0),          derv_i(1, 1))
!     call gga_eval(GGA_BECKE,  order, r(0), z(0, 0), derv_i(1, 2))
!     call gga_eval(GGA_LYP,    order, r(0), z(0, 0), derv_i(1, 3))

!   else

!     call lda_eval(LDA_SLATER, order, r(0), derv_i(1, 1))
!     call lda_eval(LDA_VWN,    order, r(0), derv_i(1, 2))

!   end if

!   multiply by the weight

    derv_i = derv_i*w

!   distribute derivatives into d

!   if (use_gga) then

      do on = 0, order
        do oz = 0, order
          sum_on_oz = on + oz
          if ((sum_on_oz > 0) .and. (sum_on_oz <= order)) then

!           this order zeroes them out
            derv(on, oz) = derv_i(taylor_index(5, (/on, 0, oz, 0, 0/)))

!           if (oz == 0) then
!             derv(on, oz) = derv(on, oz) &
!                          + taylor_index(derv_i(1, 1), 2, (/on, 0/))
!           end if

          end if
        end do
      end do

!   else

!     do on = 1, order
!       derv(on, 0) = taylor_index(derv_i(1, 1), 2, (/on, 0/)) &
!                   + taylor_index(derv_i(1, 2), 2, (/on, 0/))
!     end do

!   end if

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

    s_aux = 0.0d0

    if      (order == 5) then

      s_aux(3) = 3.0d0*r(2)*r(2)

    else if (order == 6) then

      s_aux(4) = 15.0d0*r(1)*r(2)*r(2)
      s_aux(3) = 10.0d0*r(2)*r(3)

    else if (order == 7) then

      s_aux(5) = 45.0d0*r(1)*r(1)*r(2)*r(2)
      s_aux(4) = 60.0d0*r(1)*r(2)*r(3) &
               + 15.0d0*r(2)*r(2)*r(2)
      s_aux(3) = 15.0d0*r(2)*r(4) &
               + 10.0d0*r(3)*r(3)

    else if (order == 8) then

      s_aux(6) = get_binomial(7, 2)*get_binomial(5, 2)*r(1)*r(1)*r(1)*r(2)*r(2)/2.0d0
      s_aux(5) = get_binomial(7, 2)*get_binomial(5, 2)*get_binomial(3, 2)*r(1)*r(2)*r(2)*r(2)/6.0d0 &
               + get_binomial(7, 2)*get_binomial(5, 3)*r(1)*r(1)*r(2)*r(3)
      s_aux(4) = get_binomial(7, 3)*get_binomial(4, 3)*r(1)*r(3)*r(3)/2.0d0 &
               + get_binomial(7, 4)*get_binomial(3, 2)*r(1)*r(2)*r(4)       &
               + get_binomial(7, 2)*get_binomial(5, 2)*r(2)*r(2)*r(3)/2.0d0
      s_aux(3) = get_binomial(7, 4)*r(3)*r(4) &
               + get_binomial(7, 5)*r(2)*r(5)

    end if

!   highest order terms
    s_aux(order) = r(1)**(order - 1)

!   terms that contain r(1)**(n)*r(i)**1
    do i = 2, order - 2
      k = order - 1
      n = k - i
      s_aux(n + 2) = s_aux(n + 2) + get_binomial(k, n)*r(1)**(n)*r(i)
    end do

    do i = 3, order
      s_aux(i) = s_aux(i)*2.0d0**(i - 2)
    end do

!   compensate AUTO_FUNCTIONAL_SCALING
!   do i = 1, order
!     s_aux(i) = s_aux(i)/(0.0001d0**(i))
!   end do

    do i = 1, order
      s = s + s_aux(i)*derv(i, 0)
    end do

    if (use_gga) then
      call gga_distribute(s, 2.0d0*v, 0, fmat, gao, gao1)
    else
      call lda_distribute(s,          0, fmat, gao)
    end if

  end subroutine
#endif /* ifdef MOD_AUTODIFF */

end module
