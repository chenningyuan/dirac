module ao_eval

!   this is based on DIRAC and Dalton
!   routines written by T. Helgaker, P. Salek, T. Saue, and others

! - in order to calculate aos using this module
!   you need to initialize some data first
!   by calling ao_eval_init from each node
! - ao_eval_init will set symmetry and basis set information
!   that is needed in each point and then mpi_bcast this info to the nodes
! - we try to not recalculate data which does not change from point to point
!   (like symmetry)
! - *after* the module data is initialized you can allocate
!   the ao vector and a buffer, both with the length total_ao_vector_length
! - (de)allocate them only once (each iteration), not in every point
!   this is also the reason why the ao vector and buffer are not allocated here - this would
!   be very inefficient
! - aos are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in levy-leblond or x2c4 mode some elements will
!   contain garbage but they are not (= should not be) used
! - this module is called dft_ao but should be usable also outside dft
!   this is why some information is passed through ao_eval_init and not queried
!   from dft_cfg, to keep this module generic, please dont "use dft_cfg" here

!   some "bad ideas" = ideas which did not improve things
!   so that future generations do not waste time on this:
! - ao to so transformation using dgemv: possible, but much slower than present scheme
! - derivatives of aos by numerical differentiation: again, turned out to be much slower

!   FAQ
!   ===
!
!   1. how is the ao array sorted on output of get_ao?
!      - it is sorted according to ibbas and nbbas
!        a convenient way to see the sorting of (symmetry adapted) aos
!        in the output is to search for 'Generating Lowdin matrix'
!        as an example, with a d2h calculation i got:
!
!                                 Generating Lowdin matrix:
!                                 -------------------------
!
!        L   Ag    * Deleted:          2(Proj:          2, Lindep:          0)
!        L   B1g   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B2g   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B3g   * Deleted:          0(Proj:          0, Lindep:          0)
!        S   B3u   * Deleted:          2(Proj:          2, Lindep:          0)
!        S   B2u   * Deleted:          2(Proj:          2, Lindep:          0)
!        S   B1u   * Deleted:          6(Proj:          6, Lindep:          0)
!        S   Au    * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B3u   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B2u   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B1u   * Deleted:          1(Proj:          1, Lindep:          0)
!        L   Au    * Deleted:          0(Proj:          0, Lindep:          0)
!        S   Ag    * Deleted:          9(Proj:          9, Lindep:          0)
!        S   B1g   * Deleted:          0(Proj:          0, Lindep:          0)
!        S   B2g   * Deleted:          1(Proj:          1, Lindep:          0)
!        S   B3g   * Deleted:          1(Proj:          1, Lindep:          0)
!
!        so the order is:
!        | L Ag | L B1g | L B2g | L B3g | S B3u | ... | S B3g | and then the derivatives


   use basis_info
   use ao_eval_backend
 
   implicit none
 
   public get_ao
   public ao_eval_init

   private

   integer, parameter :: max_order = 2
   integer, public    :: total_ao_vector_length = 0
   integer, public    :: ao_start     (0:max_order, 0:max_order) = 0
   integer, public    :: ao_end       (0:max_order, 0:max_order) = 0
   logical, public    :: need_ao_order(0:max_order, 0:max_order) = .false.
!  standard LDA calculation: nr_ao_slices = 1
!                            ao vector composition = | ao_g0 |
!                            meaning there are no geometric derivatives
!  standard GGA calculation: nr_ao_slices = 4
!                            ao vector composition = | ao_g0 | ao_gx | ao_gy | ao_gz |
!  with LAO there are more slices
!  also with higher geometric derivatives (soon) there are more and more slices
!  (higher derivatives and mixed derivatives)
!  slices start at ao_start and end at ao_end
!  and the total length of the ao vector is total_ao_vector_length
!  buffer has the same length for convenience (one could get away with less for the buffer)
   integer            :: nr_ao_slices = 0
   logical            :: use_ao_screening = .false.

contains

  subroutine get_ao(cx,     &
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

    if (use_ao_screening) call dzero(ao, total_ao_vector_length)

    call get_ao_sub(cx,     &
                cy,     &
                cz,     &
                ao)

!   without symmetry skip the rest
    if (nr_boson_sym == 1) return

!   transform from ao(bu) to so(bs)
    call ao_to_so(ao, buffer)
    call dcopy(total_ao_vector_length, buffer, 1, ao, 1)

  end subroutine

  subroutine get_ao_sub(cx,   &
                    cy,   &
                    cz,   &
                    ao)
!   ----------------------------------------------------------------------------
    real(8), intent(in) :: cx
    real(8), intent(in) :: cy
    real(8), intent(in) :: cz
    real(8)             :: ao(*)
!   ----------------------------------------------------------------------------
    integer             :: i, iend, ishell, isym, ioff, j
    real(8)             :: px, py, pz, p2, f1, f2
    integer             :: ofg, ofm

!#define DEBUG_SHELL_RANGE
#ifdef DEBUG_SHELL_RANGE
    character(nr_shells) :: shell_range
#endif
!   ----------------------------------------------------------------------------

    iend = nr_shells

#ifdef DEBUG_SHELL_RANGE
    do ishell = 1, iend
       shell_range(ishell:ishell) = ' '
    end do
#endif

    i = 1
    do ishell = 1, iend
      do isym = 0, nr_boson_sym - 1
        if (iand(isym, shell_stabilizer(ishell)) == 0) then

          px = cx - bitstring_parity(iand(sym_axes(1, 1), isym))*shell_center_xyz(1, ishell)
          py = cy - bitstring_parity(iand(sym_axes(2, 1), isym))*shell_center_xyz(2, ishell)
          pz = cz - bitstring_parity(iand(sym_axes(3, 1), isym))*shell_center_xyz(3, ishell)

          p2 = px*px + py*py + pz*pz

!         note that in contrast to dalton:
!       - shell_extent contains already the half_cell_diagonal
!       - is squared (for more performance)

          if (                             &
              (.not. use_ao_screening)     &
              .or.                         &
              (p2 <= shell_extent(ishell)) &
             ) then

#ifdef DEBUG_SHELL_RANGE
            shell_range(ishell:ishell) = '*'
#endif

            g_00_00_00_000 = 0.0d0
            if (need_ao_order(1, 0)) then
               g_00_00_00_100 = 0.0d0
               g_00_00_00_010 = 0.0d0
               g_00_00_00_001 = 0.0d0
            end if
            do j = 1, nr_nonzero_contraction_coefs(ishell)
               f1 = contraction_coef(ishell, j)*dexp(-primitive_exp(ishell, j)*p2)
               g_00_00_00_000 = g_00_00_00_000 + f1
               if (need_ao_order(1, 0)) then
                  f2 = -2.0d0*f1*primitive_exp(ishell, j)
                  g_00_00_00_100 = g_00_00_00_100 + px*f2
                  g_00_00_00_010 = g_00_00_00_010 + py*f2
                  g_00_00_00_001 = g_00_00_00_001 + pz*f2
               end if
            end do

!           this is because get_gao1 does both ao and ao1
            if (need_ao_order(0, 0) .and. .not. need_ao_order(1, 0)) then
              call get_gao(ishell,         &
                           ao(i),   &
                           px, py, pz, p2)
            else if (need_ao_order(1, 0)) then
              ! ao_start(1, 0) indicate the place in ao where
              ! geom derivs should start
              ofg = ao_start(1, 0) + i - 1
              call get_gao1(ishell,                   &
                            ao(i),             &
                            ao(ofg),           &
                            ao(ofg + nr_ao),    &
                            ao(ofg + 2*nr_ao),  &
                            px, py, pz, p2)
            end if

            if (need_ao_order(0, 1)) then
              ! ao_start(0, 1) indicate the place in ao where
              ! magn dervs should start
              ofm = ao_start(0, 1) + i - 1
              call get_ao_b1(px, py, pz,         &
                             ishell,                 &
                             ao(i),           &
                             ao(ofm),         &
                             ao(ofm + nr_ao),  &
                             ao(ofm + 2*nr_ao))
              if (need_ao_order(1, 1)) then
                ofm = ao_start(1, 1) + i - 1
                ofg = ao_start(1, 0) + i - 1
                call get_ao_b1(px, py, pz,         &
                               ishell,                 &
                               ao(ofg),           &
                               ao(ofm),         &
                               ao(ofm + nr_ao),  &
                               ao(ofm + 2*nr_ao))
                ofm = ofm + nr_ao*3
                ofg = ofg + nr_ao
                call get_ao_b1(px, py, pz,         &
                               ishell,                 &
                               ao(ofg),           &
                               ao(ofm),         &
                               ao(ofm + nr_ao),  &
                               ao(ofm + 2*nr_ao))
                ofm = ofm + nr_ao*3
                ofg = ofg + nr_ao
                call get_ao_b1(px, py, pz,         &
                               ishell,                 &
                               ao(ofg),           &
                               ao(ofm),         &
                               ao(ofm + nr_ao),  &
                               ao(ofm + 2*nr_ao))
              end if
            end if

!>          if spherical and d or higher then transform
!> \fixme   spherical lao: nothing is done for them now
!> \todo    this should happen outside this routine
            if (shell_is_spherical(ishell) .and. (l_quantum_nr(ishell) > 2)) then
               call transform_c_to_s(ishell, ao(i))
               if (need_ao_order(1, 0)) then
                  ioff = ao_start(1, 0) + i - 1
                  call transform_c_to_s(ishell, ao(ioff          ))
                  call transform_c_to_s(ishell, ao(ioff +   nr_ao))
                  call transform_c_to_s(ishell, ao(ioff + 2*nr_ao))
               end if
            end if

          end if

          i = i + spherical_deg(ishell)

        end if
      end do
    end do

#ifdef DEBUG_SHELL_RANGE
    write(*, *) 'shell range: [', shell_range(1:iend), ']'
#endif

  end subroutine

  subroutine get_gao(ishell, &
                     gao,    &
                     px, py, pz, p2)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8), intent(out)   :: gao(*)
    real(8), intent(in)    :: px, py, pz, p2
!   ----------------------------------------------------------------------------
    integer                :: i, ic, jc
    real(8)                :: f1, f2
!   ----------------------------------------------------------------------------

    if (cartesian_deg(ishell) < 21) then
       call get_cartesian_ao_g0(l_quantum_nr(ishell), &
                                px, py, pz,           &
                                gao)
    else
!     h and higher
      do ic = 1, cartesian_deg(ishell)
        jc = cartesian_off(ishell) + ic
        gao(ic) = (px**l_exp(jc))*(py**m_exp(jc))*(pz**n_exp(jc))*g_00_00_00_000
      end do

    end if

  end subroutine

  subroutine get_gao1(ishell, &
                      gao,    &
                      gaox,   &
                      gaoy,   &
                      gaoz,   &
                      px, py, pz, p2)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8)                :: gao(*)
    real(8)                :: gaox(*)
    real(8)                :: gaoy(*)
    real(8)                :: gaoz(*)
    real(8), intent(in)    :: px, py, pz, p2
!   ----------------------------------------------------------------------------
    integer                :: i, ic, jc
    integer                :: l, m, n
    real(8)                :: f1, f2
    real(8)                :: tx, ty, tz, ux, uy, uz
    real(8)                :: p0, p0t
!   ----------------------------------------------------------------------------

    if (cartesian_deg(ishell) < 21) then
       call get_cartesian_ao_g1(l_quantum_nr(ishell), &
                                px, py, pz,           &
                                gao,                  &
                                gaox,                 &
                                gaoy,                 &
                                gaoz)


    else
!     h and higher

!     this lmn lookup is expensive so better keep this loop outside
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

        gao (ic) = 0.0d0
        gaox(ic) = 0.0d0
        gaoy(ic) = 0.0d0
        gaoz(ic) = 0.0d0
        do i = 1, nr_nonzero_contraction_coefs(ishell)
          f2 = contraction_coef(ishell, i)*dexp(-primitive_exp(ishell, i)*p2)

          p0t = -2.0d0*primitive_exp(ishell, i)*p0
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

!> \brief transform ao to sa-ao (so)
   subroutine ao_to_so(ao, so)

!     --------------------------------------------------------------------------
      real(8), intent(in)  :: ao(total_ao_vector_length)
      real(8), intent(out) :: so(total_ao_vector_length)
!     --------------------------------------------------------------------------
      integer              :: i, n, off
!     --------------------------------------------------------------------------

      so = 0.0d0
      do n = 1, nr_ao_slices
         off = nr_ao*(n - 1)
         do i = 1, nr_so_ao
            so(off + i_so_ao(1, i)) = so(off + i_so_ao(1, i)) &
                                    + ao(off + i_so_ao(2, i))*f_so_ao(i)
         end do
      end do

   end subroutine

!> \brief initialize ao_eval module
   subroutine ao_eval_init(max_ao_g_order, &
                           max_ao_m_order, &
                           use_screening)

!     --------------------------------------------------------------------------
      integer, intent(in) :: max_ao_g_order
      integer, intent(in) :: max_ao_m_order
      logical, intent(in) :: use_screening
!     --------------------------------------------------------------------------
      integer             :: g_order
      integer             :: m_order
!     --------------------------------------------------------------------------

      use_ao_screening = use_screening
      need_ao_order    = .false.
      nr_ao_slices     = 0

      do g_order = 0, max_ao_g_order
         do m_order = 0, max_ao_m_order
            need_ao_order(g_order, m_order) = .true.
            ao_start     (g_order, m_order) = nr_ao_slices*nr_ao + 1
            nr_ao_slices = nr_ao_slices + (3**g_order)*(3**m_order)
            ao_end       (g_order, m_order) = nr_ao_slices*nr_ao
         end do
      end do

!     total_ao_vector_length = nr_ao_slices*nr_ao
!     radovan: the above can lead to "plus one" or "plus few" out-of-bounds
!              error if spherical gaussians are used
!              it is safe to do instead
      total_ao_vector_length = nr_ao_slices*nr_ao_cartesian
!              although it consumes more memory than needed
!              the allocation could be refined

   end subroutine

!> \brief transform cartesian to spherical
   subroutine transform_c_to_s(ishell, ao)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: ishell
      real(8), intent(inout) :: ao(*)
!     --------------------------------------------------------------------------
!> \fixme allocate cartesian to spherical buffer dynamically
      real(8)                :: sbuffer(10000)
!     --------------------------------------------------------------------------

      call dcopy(cartesian_deg(ishell), ao, 1, sbuffer, 1)
      call mxma(cartesian_to_spherical(spherical_address(l_quantum_nr(ishell))), &
                1,                                                               &
                spherical_deg(ishell),                                           &
                sbuffer,                                                         &
                1,                                                               &
                cartesian_deg(ishell),                                           &
                ao,                                                              &
                1,                                                               &
                spherical_deg(ishell),                                           &
                spherical_deg(ishell),                                           &
                cartesian_deg(ishell),                                           &
                1)

   end subroutine

end module
