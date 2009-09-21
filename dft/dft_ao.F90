module dft_ao

! * gao, gao1, ... are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in levy-leblond or x2c4 mode some elements will
!   contain garbage but they are not (= should not be) used

  use dft_basis_info

  implicit none

  public get_so

  private

  real(8), parameter :: tiny = 1.0d-30

contains

  subroutine get_so(cx,     &
                    cy,     &
                    cz,     &
                    ao,     &
                    buffer, &
                    ndim)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: cx
    real(8), intent(in)    :: cy
    real(8), intent(in)    :: cz
    real(8), intent(inout) :: ao(*)
    real(8), intent(inout) :: buffer(*)
    integer, intent(in)    :: ndim
!   ----------------------------------------------------------------------------
    integer                :: nvec, i
!   ----------------------------------------------------------------------------

    call get_ao(cx,     &
                cy,     &
                cz,     &
                ao,     &
                ndim)

!   without symmetry skip the rest
    if (nr_boson_sym == 1) return

    call dzero(buffer, ndim*ndimfactor)
!   transform from ao(bu) to so(bs)
    call ao_to_so(ndim, ndimfactor, ao, buffer)
    call dcopy(ndim*ndimfactor, buffer, 1, ao, 1)

  end subroutine

  subroutine get_ao(cx,   &
                    cy,   &
                    cz,   &
                    ao,   &
                    ndim)
!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: cx
    real(8), intent(in)    :: cy
    real(8), intent(in)    :: cz
    real(8)                :: ao(*)
    integer, intent(in)    :: ndim
!   ----------------------------------------------------------------------------
    integer                :: i, iend, ishell, isym
    real(8)                :: px, py, pz
    integer                :: ofg, ofm
!   ----------------------------------------------------------------------------

    iend = dft_nr_shells
!   if levy-leblond or x2c
!   don't work on small component shells
!   this will make problems if you calculate the relativistic current density
!   expression using levy-leblond, in this case comment this out
    if (nr_ls_components .eq. 1) iend = dft_nr_large_shells

    i = 1
    do ishell = 1, iend
      do isym = 0, dft_maxopr
        if (iand(isym, dft_istbao(ishell)) == 0) then
          px = cx - dft_pt(iand(dft_isymax(1, 1), isym))*dft_cent(ishell, 1, 1)
          py = cy - dft_pt(iand(dft_isymax(2, 1), isym))*dft_cent(ishell, 2, 1)
          pz = cz - dft_pt(iand(dft_isymax(3, 1), isym))*dft_cent(ishell, 3, 1)

!         this is because get_gao1 does both ao and ao1
          if (need_ao_order(0, 0) .and. .not. need_ao_order(1, 0)) then
            call get_gao(ishell,         &
                         ao(i),   &
                         px, py, pz)
          else if (need_ao_order(1, 0)) then
            ! ao_start(1, 0) indicate the place in ao where
            ! geom derivs should start
            ofg = ao_start(1, 0) + i - 1
            call get_gao1(ishell,                   &
                          ao(i),             &
                          ao(ofg),           &
                          ao(ofg + ndim),    &
                          ao(ofg + 2*ndim),  &
                          px, py, pz)
          end if

          if (need_ao_order(0, 1)) then
            ! ao_start(0, 1) indicate the place in ao where
            ! magn dervs should start
            ofm = ao_start(0, 1) + i - 1
            call get_ao_b1(px, py, pz,             &
                           ishell,                 &
                           ao(i),           &
                           ao(ofm),         &
                           ao(ofm + ndim),  &
                           ao(ofm + 2*ndim))
          end if

          i = i + spherical_deg(ishell)

        end if
      end do
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
    real(8)                :: gao(*)
    real(8)                :: gaox(*)
    real(8)                :: gaoy(*)
    real(8)                :: gaoz(*)
    real(8), intent(in)    :: px, py, pz
!   ----------------------------------------------------------------------------
    integer                :: i, ic, jc
    integer                :: l, m, n
    real(8)                :: p2, f1, f2
    real(8)                :: tx, ty, tz, ux, uy, uz
    real(8)                :: p0, p0t
!   ----------------------------------------------------------------------------

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

  subroutine get_ao_b1(px,     &
                       py,     &
                       pz,     &
                       ishell, &
                       ao,     &
                       ao_bx,  &
                       ao_by,  &
                       ao_bz)

!   ----------------------------------------------------------------------------
    real(8), intent(in)  :: px
    real(8), intent(in)  :: py
    real(8), intent(in)  :: pz
    integer, intent(in)  :: ishell
    real(8), intent(in)  :: ao(*)
    real(8), intent(out) :: ao_bx(*)
    real(8), intent(out) :: ao_by(*)
    real(8), intent(out) :: ao_bz(*)
!   ----------------------------------------------------------------------------
    integer              :: ic
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do ic = 1, cartesian_deg(ishell)

      f = 0.5d0*ao(ic)

      ao_bx(ic) = f*px
      ao_by(ic) = f*py
      ao_bz(ic) = f*pz

    end do
    
  end subroutine

  subroutine ao_to_so(ndim, nvec, ao, so)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: ndim
    integer, intent(in)  :: nvec
    real(8), intent(in)  :: ao(ndim, nvec)
    real(8), intent(out) :: so(ndim, nvec)
!   ----------------------------------------------------------------------------
    integer              :: iao, iso, ishell, nshell
    integer              :: ic, icount, iso_bs
    integer              :: i, nr_ao, mul, irep, isym, n, ivec
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

        do isym = 0, dft_maxopr
          if (iand(isym, mul) == 0) then

            f = dft_pt(iand(isym, ieor(irep, dft_isymao(n, nshell))))

            do ivec = 1, nvec
              so(iso_bs, ivec) = so(iso_bs, ivec) + f*ao(iao, ivec)
            end do

            iao = iao + spherical_deg(ishell)

          endif
        end do
      end do
      i = i + nr_ao
    end do

  end subroutine

end module
