module dft_ao_new

! * gao, gao1, ... are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in levy-leblond or x2c4 mode some elements will
!   contain garbage but they are not (= should not be) used

  use dft_cb_interface
  use dft_basis_info

  implicit none

  public dft_getsos

  private

  real(8), parameter :: tiny = 1.0d-30

contains

  subroutine dft_getsos(gao,     &
                        gao1,    &
                        cx,      &
                        cy,      &
                        cz,      &
                        buffer,  &
                        ndim,    &
                        nder)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: gao(ndim)
    real(8), intent(inout) :: gao1(ndim, 3)
    real(8), intent(in)    :: cx
    real(8), intent(in)    :: cy
    real(8), intent(in)    :: cz
    integer, intent(in)    :: ndim
    real(8), intent(inout) :: buffer(*)
    integer, intent(in)    :: nder
!   ----------------------------------------------------------------------------

!   get values in ao basis
    call get_ao(gao,  &
                gao1, &
                cx,   &
                cy,   &
                cz,   &
                ndim, &
                nder)

!   without symmetry skip the rest
    if (nr_boson_sym == 1) return

    if (nder == 0) then
      call dzero(buffer, ndim)
!     transform from ao(bu) to so(bs)
      call ao_to_so_lda(gao, &
                        buffer)
      call dcopy(ndim, buffer, 1, gao, 1)
    end if

    if (nder == 1) then
      call dzero(buffer, ndim*4)
!     transform from ao(bu) to so(bs)
      call ao_to_so_gga(gao,                 &
                        gao1(1, 1),          &
                        gao1(1, 2),          &
                        gao1(1, 3),          &
                        buffer(1),           &
                        buffer(ndim   + 1),  &
                        buffer(ndim*2 + 1),  &
                        buffer(ndim*3 + 1))
      call dcopy(ndim,   buffer,           1, gao,  1)
      call dcopy(ndim*3, buffer(ndim + 1), 1, gao1, 1)
    end if

  end subroutine

  subroutine get_ao(gao,  &
                    gao1, &
                    cx,   &
                    cy,   &
                    cz,   &
                    ndim, &
                    nder)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: gao(ndim)
    real(8), intent(inout) :: gao1(ndim, 3)
    real(8), intent(in)    :: cx
    real(8), intent(in)    :: cy
    real(8), intent(in)    :: cz
    integer, intent(in)    :: ndim
    integer, intent(in)    :: nder
!   ----------------------------------------------------------------------------
    integer                :: i, iend, ishell
    real(8)                :: px, py, pz
!   ----------------------------------------------------------------------------

    iend = dft_nr_shells
!   if levy-leblond or x2c
!   don't work on small component shells
!   this will make problems if you calculate the relativistic current density
!   expression using levy-leblond, in this case comment this out
    if (nr_ls_components .eq. 1) iend = dft_nr_large_shells

    do ishell = 1, iend

      px = cx - shift_x(ishell)
      py = cy - shift_y(ishell)
      pz = cz - shift_z(ishell)

      i = spherical_off(ishell) + 1

      if (nder == 0) then
        call get_gao(ishell, &
                     gao(i), &
                     px, py, pz)
      else if (nder == 1) then
        call get_gao1(ishell,     &
                      gao(i),     &
                      gao1(i, 1), &
                      gao1(i, 2), &
                      gao1(i, 3), &
                      px, py, pz)
      end if

    end do

  end subroutine

  subroutine get_gao(ishell, &
                     gao,    &
                     px, py, pz)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8), intent(out)   :: gao(*)
    real(8), intent(in)    :: px, py, pz
!   ----------------------------------------------------------------------------
    integer                :: i, ic, jc
    real(8)                :: p2, f1, f2
!   ----------------------------------------------------------------------------

    p2 = px**2.0d0 + py**2.0d0 + pz**2.0d0

    f1 = 0.0d0
    if (dabs(contraction_coef1(ishell)) > tiny) then
      f1 = contraction_coef1(ishell)*dexp(-primitive_exp1(ishell)*p2)
    end if

    do ic = 1, cartesian_deg(ishell)
      gao(ic) = f1
    end do

    do i = 1, nr_nonzero_contraction_coefs(ishell)
      f2 = contraction_coef2(ishell, i)*dexp(-primitive_exp2(ishell, i)*p2)
      do ic = 1, cartesian_deg(ishell)
        gao(ic) = gao(ic) + f2
      end do
    end do

!   p
    if (cartesian_deg(ishell) == 3) then

      gao(1) = px*gao(1)
      gao(2) = py*gao(2)
      gao(3) = pz*gao(3)

!   d and higher
    else if (cartesian_deg(ishell) > 3) then

      do ic = 1, cartesian_deg(ishell)
        jc = cartesian_off(ishell) + ic
        gao(ic) = (px**l_exp(jc))*(py**m_exp(jc))*(pz**n_exp(jc))*gao(ic)
      end do

    end if

!   if spherical and d or higher then transform
    if (is_spherical(ishell) .and. (cartesian_deg(ishell) > 3)) then
      call dcopy(cartesian_deg(ishell), gao, 1, sbuffer, 1)
      call mxma(dft_csp(dft_ispadr(n_quantum_nr(ishell))), &
                1,                                         &
                spherical_deg(ishell),                     &
                sbuffer,                                   &
                1,                                         &
                cartesian_deg(ishell),                     &
                gao,                                       &
                1,                                         &
                spherical_deg(ishell),                     &
                spherical_deg(ishell),                     &
                cartesian_deg(ishell),                     &
                1)
    end if

  end subroutine

  subroutine get_gao1(ishell, &
                      gao,    &
                      gaox,   &
                      gaoy,   &
                      gaoz,   &
                      px, py, pz)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8), intent(out)   :: gao(*)
    real(8), intent(out)   :: gaox(*)
    real(8), intent(out)   :: gaoy(*)
    real(8), intent(out)   :: gaoz(*)
    real(8), intent(in)    :: px, py, pz
!   ----------------------------------------------------------------------------
    integer                :: i, ic, jc
    real(8)                :: p2, f1, f2
!   ----------------------------------------------------------------------------


      integer :: l, m, n
!     real(8) :: tga, tgx, tgy, tgz, tx, ty, tz
!     real(8) ::                    p0
      real(8) :: gax, gay, gaz ,ga
      real(8) :: tx, ty, tz, ux, uy, uz, p0, p0t
      real(8) :: tp0

    p2 = px**2.0d0 + py**2.0d0 + pz**2.0d0

    f1 = 0.0d0
    if (dabs(contraction_coef1(ishell)) > tiny) then
      f1 = contraction_coef1(ishell)*dexp(-primitive_exp1(ishell)*p2)
    end if

    if (cartesian_deg(ishell) < 4) then

    do ic = 1, cartesian_deg(ishell)
      gao (ic) = f1
      gaox(ic) = primitive_exp1(ishell)*f1
    end do

    do i = 1, nr_nonzero_contraction_coefs(ishell)
      f2 = contraction_coef2(ishell, i)*dexp(-primitive_exp2(ishell, i)*p2)
      do ic = 1, cartesian_deg(ishell)
        gao (ic) = gao (ic) + f2
        gaox(ic) = gaox(ic) + primitive_exp2(ishell, i)*f2
      end do
    end do

    do ic = 1, cartesian_deg(ishell)
!     gaox(ic) = -2.0d0*px*gaox(ic)
!     gaoy(ic) = -2.0d0*py*gaoy(ic)
!     gaoz(ic) = -2.0d0*pz*gaoz(ic)
!     use symmetry
      gaoz(ic) = -2.0d0*pz*gaox(ic)
      gaoy(ic) = -2.0d0*py*gaox(ic)
      gaox(ic) = -2.0d0*px*gaox(ic)
    end do

    end if

    if (cartesian_deg(ishell) == 3) then

    gao (1) = px*gao (1)
    gao (2) = py*gao (2)
    gao (3) = pz*gao (3)
    gaox(1) = px*gaox(1)
    gaox(2) = py*gaox(2)
    gaox(3) = pz*gaox(3)
    gaoy(1) = px*gaoy(1)
    gaoy(2) = py*gaoy(2)
    gaoy(3) = pz*gaoy(3)
    gaoz(1) = px*gaoz(1)
    gaoz(2) = py*gaoz(2)
    gaoz(3) = pz*gaoz(3)

    gaox(1) = gaox(1) + f1
    gaoy(2) = gaoy(2) + f1
    gaoz(3) = gaoz(3) + f1
    do i = 1, nr_nonzero_contraction_coefs(ishell)
      f2 = contraction_coef2(ishell, i)*dexp(-primitive_exp2(ishell, i)*p2)
      gaox(1) = gaox(1) + f2
      gaoy(2) = gaoy(2) + f2
      gaoz(3) = gaoz(3) + f2
    end do

    end if

!   this lmn lookup is expensive so better keep this loop outside
    if (cartesian_deg(ishell) > 3) then

      do ic = 1,cartesian_deg(ishell)
        jc = cartesian_off(ishell) + ic

        l = l_exp(jc)
        m = m_exp(jc)
        n = n_exp(jc)

        tx = 0.0d0
        ty = 0.0d0
        tz = 0.0d0

        if (l > 0) tx = dfloat(l)*(px**(l-1))*(py**m)*(pz**n)
        if (m > 0) ty = dfloat(m)*(px**l)*(py**(m-1))*(pz**n)
        if (n > 0) tz = dfloat(n)*(px**l)*(py**m)*(pz**(n-1))

        p0 = (px**l)*(py**m)*(pz**n)

        p0t = -2.0d0*primitive_exp1(ishell)*p0
        ux  = px*p0t + tx
        uy  = py*p0t + ty
        uz  = pz*p0t + tz

        gao (ic) = p0*f1
        gaox(ic) = ux*f1
        gaoy(ic) = uy*f1
        gaoz(ic) = uz*f1

        do i = 1, nr_nonzero_contraction_coefs(ishell)
          f2 = contraction_coef2(ishell, i)*dexp(-primitive_exp2(ishell, i)*p2)

          p0t = -2.0d0*primitive_exp2(ishell, i)*p0
          ux  = px*p0t + tx
          uy  = py*p0t + ty
          uz  = pz*p0t + tz

          gao (ic) = gao (ic) + p0*f2
          gaox(ic) = gaox(ic) + ux*f2
          gaoy(ic) = gaoy(ic) + uy*f2
          gaoz(ic) = gaoz(ic) + uz*f2

        end do
      end do

    end if

!   if spherical and d or higher then transform
    if (is_spherical(ishell) .and. (cartesian_deg(ishell) > 3)) then
      call dcopy(cartesian_deg(ishell), gao, 1, sbuffer, 1)
      call mxma(dft_csp(dft_ispadr(n_quantum_nr(ishell))), &
                1,                                         &
                spherical_deg(ishell),                     &
                sbuffer,                                   &
                1,                                         &
                cartesian_deg(ishell),                     &
                gao,                                       &
                1,                                         &
                spherical_deg(ishell),                     &
                spherical_deg(ishell),                     &
                cartesian_deg(ishell),                     &
                1)
      call dcopy(cartesian_deg(ishell), gaox, 1, sbuffer, 1)
      call mxma(dft_csp(dft_ispadr(n_quantum_nr(ishell))), &
                1,                                         &
                spherical_deg(ishell),                     &
                sbuffer,                                   &
                1,                                         &
                cartesian_deg(ishell),                     &
                gaox,                                      &
                1,                                         &
                spherical_deg(ishell),                     &
                spherical_deg(ishell),                     &
                cartesian_deg(ishell),                     &
                1)
      call dcopy(cartesian_deg(ishell), gaoy, 1, sbuffer, 1)
      call mxma(dft_csp(dft_ispadr(n_quantum_nr(ishell))), &
                1,                                         &
                spherical_deg(ishell),                     &
                sbuffer,                                   &
                1,                                         &
                cartesian_deg(ishell),                     &
                gaoy,                                      &
                1,                                         &
                spherical_deg(ishell),                     &
                spherical_deg(ishell),                     &
                cartesian_deg(ishell),                     &
                1)
      call dcopy(cartesian_deg(ishell), gaoz, 1, sbuffer, 1)
      call mxma(dft_csp(dft_ispadr(n_quantum_nr(ishell))), &
                1,                                         &
                spherical_deg(ishell),                     &
                sbuffer,                                   &
                1,                                         &
                cartesian_deg(ishell),                     &
                gaoz,                                      &
                1,                                         &
                spherical_deg(ishell),                     &
                spherical_deg(ishell),                     &
                cartesian_deg(ishell),                     &
                1)
    end if

  end subroutine

  subroutine ao_to_so_lda(ao, so)

!   ----------------------------------------------------------------------------
    real(8), intent(in)  :: ao(*)
    real(8), intent(out) :: so(*)
!   ----------------------------------------------------------------------------
    integer              :: iao, iso, ishell, nshell
    integer              :: ic, icount, iso_bs
    integer              :: i, nr_ao, mul, irep, isym, n
    real(8)              :: f
!   ----------------------------------------------------------------------------

    i = 1
    do irep = 0, dft_maxrep
      nr_ao = dft_nbbas(irep, 1) + dft_nbbas(irep, 2)
      icount = 0
      do iso = i, i + nr_ao - 1

        ishell = iand(ishft(dft_ipind(iso), -16), 65535)
        nshell = iand(ishft(dft_ipind(iso),  -8),   255)

        mul = dft_istbao(ishell)
        n   = n_quantum_nr(ishell)

        do isym = 0, dft_maxopr
          if (iand(isym, mul) == 0) then
            f  = dft_pt(iand(isym, ieor(irep, dft_isymao(n, nshell))))
          endif
        end do

        iao = dft_kstrt(ishell) + nshell

        icount = icount + 1
        if (icount > dft_nbbas(irep, 1)) then
!         is small
          ic = 2
        else
!         is large
          ic = 1
        end if

!       shift by difference boson sorted minus boson unsorted
        iso_bs = iso + dft_ibbas(irep, ic) - dft_icos(irep + 1, ic)

        so(iso_bs) = so(iso_bs) + f*ao(iao)

      end do
      i = i + nr_ao
    end do

  end subroutine

  subroutine ao_to_so_gga(ao, aox, aoy, aoz, so, sox, soy, soz)

!   ----------------------------------------------------------------------------
    real(8), intent(in)  :: ao(*)
    real(8), intent(in)  :: aox(*)
    real(8), intent(in)  :: aoy(*)
    real(8), intent(in)  :: aoz(*)
    real(8), intent(out) :: so(*)
    real(8), intent(out) :: sox(*)
    real(8), intent(out) :: soy(*)
    real(8), intent(out) :: soz(*)
!   ----------------------------------------------------------------------------
    integer              :: iao, iso, ishell, nshell
    integer              :: ic, icount, iso_bs
    integer              :: i, nr_ao, mul, irep, isym, n
    real(8)              :: f
!   ----------------------------------------------------------------------------

    i = 1
    do irep = 0, dft_maxrep
      nr_ao = dft_nbbas(irep, 1) + dft_nbbas(irep, 2)
      icount = 0
      do iso = i, i + nr_ao - 1

        ishell = iand(ishft(dft_ipind(iso), -16), 65535)
        nshell = iand(ishft(dft_ipind(iso),  -8),   255)

        mul = dft_istbao(ishell)
        n   = n_quantum_nr(ishell)

        do isym = 0, dft_maxopr
          if (iand(isym, mul) == 0) then
            f  = dft_pt(iand(isym, ieor(irep, dft_isymao(n, nshell))))
          endif
        end do

        iao = dft_kstrt(ishell) + nshell

        icount = icount + 1
        if (icount > dft_nbbas(irep, 1)) then
!         is small
          ic = 2
        else
!         is large
          ic = 1
        end if

!       shift by difference boson sorted minus boson unsorted
        iso_bs = iso + dft_ibbas(irep, ic) - dft_icos(irep + 1, ic)

        so (iso_bs) = so (iso_bs) + f*ao (iao)
        sox(iso_bs) = sox(iso_bs) + f*aox(iao)
        soy(iso_bs) = soy(iso_bs) + f*aoy(iao)
        soz(iso_bs) = soz(iso_bs) + f*aoz(iao)

      end do
      i = i + nr_ao
    end do

  end subroutine

end module
