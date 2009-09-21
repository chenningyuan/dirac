module dft_ao

! - in order to calculate aos using this module
!   you need to initialize some data first
!   by calling dft_ao_init from each node
! - *after* the module data is initialized you can allocate
!   the ao vector and a buffer, both with the length total_ao_vector_length*mxbllen
! - (de)allocate them only once (each iteration), not in every point
! - aos are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in levy-leblond or x2c4 mode some elements will
!   contain garbage but they are not (= should not be) used

  use memory_allocator
#ifdef VAR_MPI
  use dirac_mpi
#endif

  implicit none

  public get_so
  public dft_ao_init

  save

  integer, public    :: total_ao_vector_length = 0

  integer, parameter :: max_order = 2
  integer, public    :: ao_start     (0:max_order, 0:max_order) = 0
  integer, public    :: ao_end       (0:max_order, 0:max_order) = 0
  logical, public    :: need_ao_order(0:max_order, 0:max_order) = .false.

  private

  integer              :: nr_ao_slices   = 0
  integer              :: nr_ao               = 0
  integer              :: nr_cartesian_ao     = 0
  integer              :: dft_nr_shells       = 0
  integer              :: dft_nr_large_shells = 0
  integer              :: dft_maxrep          = 0
  integer              :: dft_maxopr          = 0
  integer              :: nr_boson_sym        = 0
  integer              :: nr_ls_components    = 0

  real(8), allocatable :: primitive_exp1(:)
  real(8), allocatable :: primitive_exp2(:, :)
  real(8), allocatable :: contraction_coef1(:)
  real(8), allocatable :: contraction_coef2(:, :)
  integer, allocatable :: nr_nonzero_contraction_coefs(:)
  integer, allocatable :: cartesian_deg(:)
  integer, allocatable :: cartesian_off(:)
  integer, allocatable :: spherical_deg(:)
  integer, allocatable :: spherical_off(:)
  integer, allocatable :: n_quantum_nr(:)
  logical, allocatable :: is_spherical(:)
  integer, allocatable :: dft_istbao(:)
  integer, allocatable :: dft_isymao(:, :)
  integer, allocatable :: dft_ipind(:)
  real(8), allocatable :: dft_pt(:)
  integer, allocatable :: dft_kstrt(:)
  real(8), allocatable :: dft_cent(:, :, :)
  integer              :: dft_isymax(3, 2)

  real(8), allocatable :: dft_csp(:)
  integer, allocatable :: dft_ispadr(:)

  integer, allocatable :: dft_icos(:, :)
  integer, allocatable :: dft_ibbas(:, :)
  integer, allocatable :: dft_nbbas(:, :)

  integer, allocatable :: l_exp(:)
  integer, allocatable :: m_exp(:)
  integer, allocatable :: n_exp(:)
  real(8), allocatable :: sbuffer(:)

  real(8), parameter   :: tiny = 1.0d-30

  logical :: dft_basis_info_allocated = .false.
  logical :: dft_basis_info_set       = .false.
#ifdef VAR_MPI
  logical :: dft_basis_info_synced    = .false.
#endif

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "shells.h"
#include "symmet.h"
#include "aovec.h"
#include "primit.h"
#include "sphtrm.h"
#include "dcbbas.h"
#include "pincom.h"
#include "dgroup.h"
#include "dcbham.h"

#ifdef VAR_MPI
#include "infpar.h"
#include "mpif.h"
#endif

contains

  subroutine get_so(cx,     &
                    cy,     &
                    cz,     &
                    ao,     &
                    buffer)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: cx
    real(8), intent(in)    :: cy
    real(8), intent(in)    :: cz
    real(8), intent(inout) :: ao(*)
    real(8), intent(inout) :: buffer(*)
!   ----------------------------------------------------------------------------

    call get_ao(cx,     &
                cy,     &
                cz,     &
                ao)

!   without symmetry skip the rest
    if (nr_boson_sym == 1) return

!   transform from ao(bu) to so(bs)
    call ao_to_so(ao, buffer)
    call dcopy(total_ao_vector_length, buffer, 1, ao, 1)

  end subroutine

  subroutine get_ao(cx,   &
                    cy,   &
                    cz,   &
                    ao)
!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: cx
    real(8), intent(in)    :: cy
    real(8), intent(in)    :: cz
    real(8)                :: ao(*)
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
                          ao(ofg + nr_ao),    &
                          ao(ofg + 2*nr_ao),  &
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
                           ao(ofm + nr_ao),  &
                           ao(ofm + 2*nr_ao))
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

  subroutine ao_to_so(ao, so)

!   ----------------------------------------------------------------------------
    real(8), intent(in)  :: ao(*)
    real(8), intent(out) :: so(*)
!   ----------------------------------------------------------------------------
    integer              :: iao, iso, ishell, nshell
    integer              :: ic, icount, iso_bs
    integer              :: i, nr_ao_irep, mul, irep, isym, n, ivec, off
    real(8)              :: f
!   ----------------------------------------------------------------------------

    call dzero(so, total_ao_vector_length) 

    i = 1
    do irep = 0, dft_maxrep
      nr_ao_irep = dft_nbbas(irep, 1) + dft_nbbas(irep, 2)
      icount = 0
      do iso = i, i + nr_ao_irep - 1

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

            do ivec = 1, nr_ao_slices
              off = nr_ao*(ivec - 1)
              so(off + iso_bs) = so(off + iso_bs) + f*ao(off + iao)
            end do

            iao = iao + spherical_deg(ishell)

          endif
        end do
      end do
      i = i + nr_ao_irep
    end do

  end subroutine

  subroutine dft_ao_init(is_parallel_calculation, &
                         max_ao_g_order,          &
                         max_ao_m_order)

!   ----------------------------------------------------------------------------
    logical, intent(in) :: is_parallel_calculation
    integer, intent(in) :: max_ao_g_order
    integer, intent(in) :: max_ao_m_order
!   ----------------------------------------------------------------------------
    integer             :: g_order
    integer             :: m_order
    integer             :: iprim, i, j, icount, i_c, i_s, contracted_shell_nr, ishell, isym
#ifdef VAR_MPI
    logical             :: i_am_the_master
#endif
!   ----------------------------------------------------------------------------

    if (.not. dft_basis_info_allocated) then

      call alloc(primitive_exp1, mxshel)
      call alloc(primitive_exp2, mxshel, mxcont)
      call alloc(contraction_coef1, mxshel)
      call alloc(contraction_coef2, mxshel, mxcont)
      call alloc(nr_nonzero_contraction_coefs, mxshel)
      call alloc(cartesian_deg, mxshel)
      call alloc(cartesian_off, mxshel)
      call alloc(spherical_deg, mxshel)
      call alloc(spherical_off, mxshel)
      call alloc(n_quantum_nr, mxshel)
      allocate(is_spherical(mxshel))
      call alloc(dft_istbao, mxshel)
      call alloc(dft_isymao, mxqn, mxaqn)
      call alloc(dft_ipind, mxcorb)
      allocate(dft_pt(0:7))
      call alloc(dft_kstrt, mxshel)
      call alloc(dft_cent, mxshel, 3, 2)

      nr_cartesian_ao = 0
      do ishell = 1, mxshel
        nr_cartesian_ao = nr_cartesian_ao + kckt(ishell)
      end do
      dft_nr_shells       = kmax
      dft_nr_large_shells = nlrgsh

      call alloc(dft_csp,    ncsp)
      call alloc(dft_ispadr, mxqn)

      call alloc(dft_icos, 8, 2)
      allocate(dft_ibbas(0:7, 0:2))
      allocate(dft_nbbas(0:7, 0:2))

#ifdef VAR_MPI
      if (is_parallel_calculation) then
        call dirac_mpi_bcast(nr_cartesian_ao)
      end if
#endif

      call alloc(l_exp,   nr_cartesian_ao)
      call alloc(m_exp,   nr_cartesian_ao)
      call alloc(n_exp,   nr_cartesian_ao)
      call alloc(sbuffer, nr_cartesian_ao)

      dft_basis_info_allocated = .true.

    end if

#ifdef VAR_MPI
    i_am_the_master = .false.
    if (mytid == mparid) i_am_the_master = .true.

!   only master sets values
    if (i_am_the_master) then
#endif

    if (.not. dft_basis_info_set) then

      nr_nonzero_contraction_coefs = 0
      cartesian_deg                = 0
      cartesian_off                = 0
      spherical_deg                = 0
      spherical_off                = 0
      n_quantum_nr                 = 0
      is_spherical                 = .false.

      dft_ipind  = ipind
      dft_istbao = istbao
      dft_isymao = isymao
      dft_pt     = pt
      dft_kstrt  = kstrt
      dft_cent   = cent
      dft_isymax = isymax
      nr_boson_sym                = nbsym
      nr_ls_components            = mc
      nr_ao = ntbas(0)

      dft_csp    = csp
      dft_ispadr = ispadr

      dft_icos  = icos
      dft_ibbas = ibbas
      dft_nbbas = nbbas

      dft_maxrep = maxrep
      dft_maxopr = maxopr

      i_c = 1
      i_s = 1

      do ishell = 1, kmax

        n_quantum_nr(ishell)  = nhkt(ishell)

        cartesian_deg(ishell) = kckt(ishell)
        spherical_deg(ishell) = khkt(ishell)

        cartesian_off(ishell) = i_c - 1
        spherical_off(ishell) = i_s - 1

        call lmnval(n_quantum_nr(ishell), cartesian_deg(ishell), l_exp(i_c), m_exp(i_c), n_exp(i_c))

        i_c = i_c + cartesian_deg(ishell)
        i_s = i_s + spherical_deg(ishell)

        is_spherical(ishell)  = sphr(ishell)

        contracted_shell_nr = numcf(ishell)

        i = jstrt(ishell) + 1
        j = jstrt(ishell) + nuco(ishell)

        contraction_coef1(ishell) = priccf(i, contracted_shell_nr)
        primitive_exp1(ishell)    = priexp(i)

        icount = 0
        do iprim = i + 1, j
          if (dabs(priccf(iprim, contracted_shell_nr)) > 1.0d-30) then
            icount = icount + 1
            nr_nonzero_contraction_coefs(ishell) = icount
            contraction_coef2(ishell, icount)    = priccf(iprim, contracted_shell_nr)
            primitive_exp2(ishell, icount)       = priexp(iprim)
          end if
        end do

      end do

      dft_basis_info_set = .true.

    end if

    need_ao_order = .false.
    nr_ao_slices  = 0

    do g_order = 0, max_ao_g_order
      do m_order = 0, max_ao_m_order
        need_ao_order(g_order, m_order) = .true.
        ao_start     (g_order, m_order) = nr_ao_slices*nr_ao + 1
        nr_ao_slices = nr_ao_slices + (3**g_order)*(3**m_order)
        ao_end       (g_order, m_order) = nr_ao_slices*nr_ao
      end do
    end do

    total_ao_vector_length = nr_ao_slices*nr_ao

#ifdef VAR_MPI
    end if !if (i_am_the_master) then
#endif

#ifdef VAR_MPI
    if (is_parallel_calculation .and. .not. dft_basis_info_synced) then

!     is done above
!     call dirac_mpi_bcast(nr_cartesian_ao             )
      call dirac_mpi_bcast(dft_nr_shells               )
      call dirac_mpi_bcast(dft_nr_large_shells         )
      call dirac_mpi_bcast(dft_maxrep                  )
      call dirac_mpi_bcast(dft_maxopr                  )
      call dirac_mpi_bcast(nr_boson_sym                )
      call dirac_mpi_bcast(nr_ls_components            )
      call dirac_mpi_bcast(nr_ao                       )

      call dirac_mpi_bcast(primitive_exp1              )
      call dirac_mpi_bcast(primitive_exp2              )
      call dirac_mpi_bcast(contraction_coef1           )
      call dirac_mpi_bcast(contraction_coef2           )
      call dirac_mpi_bcast(nr_nonzero_contraction_coefs)
      call dirac_mpi_bcast(cartesian_deg               )
      call dirac_mpi_bcast(cartesian_off               )
      call dirac_mpi_bcast(spherical_deg               )
      call dirac_mpi_bcast(spherical_off               )
      call dirac_mpi_bcast(n_quantum_nr                )
      call dirac_mpi_bcast(is_spherical                )
      call dirac_mpi_bcast(dft_istbao                  )
      call dirac_mpi_bcast(dft_isymao                  )
      call dirac_mpi_bcast(dft_ipind                   )
      call dirac_mpi_bcast(dft_pt                      )
      call dirac_mpi_bcast(dft_kstrt                   )
      call dirac_mpi_bcast(dft_cent                    )
      call dirac_mpi_bcast(dft_isymax                  )

      call dirac_mpi_bcast(dft_csp                     )
      call dirac_mpi_bcast(dft_ispadr                  )

      call dirac_mpi_bcast(dft_icos                    )
      call dirac_mpi_bcast(dft_ibbas                   )
      call dirac_mpi_bcast(dft_nbbas                   )

      call dirac_mpi_bcast(l_exp                       )
      call dirac_mpi_bcast(m_exp                       )
      call dirac_mpi_bcast(n_exp                       )

      call dirac_mpi_bcast(nr_ao_slices    )
      call dirac_mpi_bcast(need_ao_order )
      call dirac_mpi_bcast(ao_start      )
      call dirac_mpi_bcast(ao_end        )
      call dirac_mpi_bcast(total_ao_vector_length)

      dft_basis_info_synced = .true.

    end if
#endif

  end subroutine

end module
