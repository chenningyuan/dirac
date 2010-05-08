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

  real(8) :: g000_000

  real(8) :: g000_100
  real(8) :: g000_010
  real(8) :: g000_001

  real(8) :: g100_000
  real(8) :: g010_000
  real(8) :: g001_000

  real(8) :: g100_100
  real(8) :: g010_100
  real(8) :: g001_100
  real(8) :: g100_010
  real(8) :: g010_010
  real(8) :: g001_010
  real(8) :: g100_001
  real(8) :: g010_001
  real(8) :: g001_001

  real(8) :: g200_000
  real(8) :: g110_000
  real(8) :: g101_000
  real(8) :: g020_000
  real(8) :: g011_000
  real(8) :: g002_000

  real(8) :: g200_100
  real(8) :: g110_100
  real(8) :: g101_100
  real(8) :: g020_100
  real(8) :: g011_100
  real(8) :: g002_100
  real(8) :: g200_010
  real(8) :: g110_010
  real(8) :: g101_010
  real(8) :: g020_010
  real(8) :: g011_010
  real(8) :: g002_010
  real(8) :: g200_001
  real(8) :: g110_001
  real(8) :: g101_001
  real(8) :: g020_001
  real(8) :: g011_001
  real(8) :: g002_001

  real(8) :: g300_000
  real(8) :: g210_000
  real(8) :: g201_000
  real(8) :: g120_000
  real(8) :: g111_000
  real(8) :: g102_000
  real(8) :: g030_000
  real(8) :: g021_000
  real(8) :: g012_000
  real(8) :: g003_000
  real(8) :: g300_100
  real(8) :: g210_100
  real(8) :: g201_100
  real(8) :: g120_100
  real(8) :: g111_100
  real(8) :: g102_100
  real(8) :: g030_100
  real(8) :: g021_100
  real(8) :: g012_100
  real(8) :: g003_100
  real(8) :: g300_010
  real(8) :: g210_010
  real(8) :: g201_010
  real(8) :: g120_010
  real(8) :: g111_010
  real(8) :: g102_010
  real(8) :: g030_010
  real(8) :: g021_010
  real(8) :: g012_010
  real(8) :: g003_010
  real(8) :: g300_001
  real(8) :: g210_001
  real(8) :: g201_001
  real(8) :: g120_001
  real(8) :: g111_001
  real(8) :: g102_001
  real(8) :: g030_001
  real(8) :: g021_001
  real(8) :: g012_001
  real(8) :: g003_001

  real(8) :: g400_000
  real(8) :: g310_000
  real(8) :: g301_000
  real(8) :: g220_000
  real(8) :: g211_000
  real(8) :: g202_000
  real(8) :: g130_000
  real(8) :: g121_000
  real(8) :: g112_000
  real(8) :: g103_000
  real(8) :: g040_000
  real(8) :: g031_000
  real(8) :: g022_000
  real(8) :: g013_000
  real(8) :: g004_000
  real(8) :: g400_100
  real(8) :: g310_100
  real(8) :: g301_100
  real(8) :: g220_100
  real(8) :: g211_100
  real(8) :: g202_100
  real(8) :: g130_100
  real(8) :: g121_100
  real(8) :: g112_100
  real(8) :: g103_100
  real(8) :: g040_100
  real(8) :: g031_100
  real(8) :: g022_100
  real(8) :: g013_100
  real(8) :: g004_100
  real(8) :: g400_010
  real(8) :: g310_010
  real(8) :: g301_010
  real(8) :: g220_010
  real(8) :: g211_010
  real(8) :: g202_010
  real(8) :: g130_010
  real(8) :: g121_010
  real(8) :: g112_010
  real(8) :: g103_010
  real(8) :: g040_010
  real(8) :: g031_010
  real(8) :: g022_010
  real(8) :: g013_010
  real(8) :: g004_010
  real(8) :: g400_001
  real(8) :: g310_001
  real(8) :: g301_001
  real(8) :: g220_001
  real(8) :: g211_001
  real(8) :: g202_001
  real(8) :: g130_001
  real(8) :: g121_001
  real(8) :: g112_001
  real(8) :: g103_001
  real(8) :: g040_001
  real(8) :: g031_001
  real(8) :: g022_001
  real(8) :: g013_001
  real(8) :: g004_001

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
    integer             :: i, iend, ishell, isym
    real(8)             :: px, py, pz, p2
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

          end if

          i = i + spherical_deg(ishell)

        end if
      end do
    end do

#ifdef DEBUG_SHELL_RANGE
    write(*, *) 'shell range: [', shell_range(1:iend), ']'
#endif

  end subroutine

  subroutine get_cartesian_ao(ishell, &
                              gao,    &
                              px, py, pz, p2)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8), intent(out)   :: gao(*)
    real(8), intent(in)    :: px, py, pz, p2
!   ----------------------------------------------------------------------------
    integer                :: i
!   ----------------------------------------------------------------------------

    g000_000 = 0.0d0
    do i = 1, nr_nonzero_contraction_coefs(ishell)
      g000_000 = g000_000 + contraction_coef(ishell, i)*dexp(-primitive_exp(ishell, i)*p2)
    end do

    if (cartesian_deg(ishell) > 1) then
      g100_000 = g000_000*px
      g010_000 = g000_000*py
      g001_000 = g000_000*pz
    else
!     s
      gao(1)   = g000_000
      return
    end if

    if (cartesian_deg(ishell) > 3) then
      g200_000 = g100_000*px
      g110_000 = g100_000*py
      g101_000 = g100_000*pz
      g020_000 = g010_000*py
      g011_000 = g010_000*pz
      g002_000 = g001_000*pz
    else
!     p
      gao(1)   = g100_000
      gao(2)   = g010_000
      gao(3)   = g001_000
      return
    end if

    if (cartesian_deg(ishell) > 6) then
      g300_000 = g200_000*px
      g210_000 = g200_000*py
      g201_000 = g200_000*pz
      g120_000 = g110_000*py
      g111_000 = g110_000*pz
      g102_000 = g101_000*pz
      g030_000 = g020_000*py
      g021_000 = g020_000*pz
      g012_000 = g011_000*pz
      g003_000 = g002_000*pz
    else
!     d
      gao(1)   = g200_000
      gao(2)   = g110_000
      gao(3)   = g101_000
      gao(4)   = g020_000
      gao(5)   = g011_000
      gao(6)   = g002_000
      return
    end if

    if (cartesian_deg(ishell) > 10) then
      g400_000 = g300_000*px
      g310_000 = g300_000*py
      g301_000 = g300_000*pz
      g220_000 = g210_000*py
      g211_000 = g210_000*pz
      g202_000 = g201_000*pz
      g130_000 = g120_000*py
      g121_000 = g120_000*pz
      g112_000 = g111_000*pz
      g103_000 = g102_000*pz
      g040_000 = g030_000*py
      g031_000 = g030_000*pz
      g022_000 = g021_000*pz
      g013_000 = g012_000*pz
      g004_000 = g003_000*pz
    else
!     f
      gao( 1)  = g300_000
      gao( 2)  = g210_000
      gao( 3)  = g201_000
      gao( 4)  = g120_000
      gao( 5)  = g111_000
      gao( 6)  = g102_000
      gao( 7)  = g030_000
      gao( 8)  = g021_000
      gao( 9)  = g012_000
      gao(10)  = g003_000
      return
    end if

    if (cartesian_deg(ishell) > 15) then
    else
!     g
      gao( 1)  = g400_000
      gao( 2)  = g310_000
      gao( 3)  = g301_000
      gao( 4)  = g220_000
      gao( 5)  = g211_000
      gao( 6)  = g202_000
      gao( 7)  = g130_000
      gao( 8)  = g121_000
      gao( 9)  = g112_000
      gao(10)  = g103_000
      gao(11)  = g040_000
      gao(12)  = g031_000
      gao(13)  = g022_000
      gao(14)  = g013_000
      gao(15)  = g004_000
      return
    end if

  end subroutine

  subroutine get_cartesian_ao1(ishell, &
                              gao,    &
                              gaox,    &
                              gaoy,    &
                              gaoz,    &
                              px, py, pz, p2)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8), intent(out)   :: gao(*)
    real(8), intent(out)   :: gaox(*)
    real(8), intent(out)   :: gaoy(*)
    real(8), intent(out)   :: gaoz(*)
    real(8), intent(in)    :: px, py, pz, p2
!   ----------------------------------------------------------------------------
    integer                :: i
    real(8)                :: f, f2
!   ----------------------------------------------------------------------------

    g000_000 = 0.0d0
    g000_100 = 0.0d0
    g000_010 = 0.0d0
    g000_001 = 0.0d0

    do i = 1, nr_nonzero_contraction_coefs(ishell)
      f = -primitive_exp(ishell, i)
      f2 = contraction_coef(ishell, i)*dexp(f*p2)
      g000_000 = g000_000 +            f2
      g000_100 = g000_100 + px*f*2.0d0*f2
      g000_010 = g000_010 + py*f*2.0d0*f2
      g000_001 = g000_001 + pz*f*2.0d0*f2
    end do

    if (cartesian_deg(ishell) > 1) then
      g100_000 =            px*g000_000
      g010_000 =            py*g000_000
      g001_000 =            pz*g000_000
      g100_100 = g000_000 + px*g000_100
      g010_100 =            py*g000_100
      g001_100 =            pz*g000_100
      g100_010 =            px*g000_010
      g010_010 = g000_000 + py*g000_010
      g001_010 =            pz*g000_010
      g100_001 =            px*g000_001
      g010_001 =            py*g000_001
      g001_001 = g000_000 + pz*g000_001
    else
!     s
      gao(1)  = g000_000
      gaox(1) = g000_100
      gaoy(1) = g000_010
      gaoz(1) = g000_001
      return
    end if

    if (cartesian_deg(ishell) > 3) then
      g200_000 =            px*g100_000
      g110_000 =            py*g100_000
      g101_000 =            pz*g100_000
      g020_000 =            py*g010_000
      g011_000 =            pz*g010_000
      g002_000 =            pz*g001_000
      g200_100 = g100_000 + px*g100_100
      g110_100 =            py*g100_100
      g101_100 =            pz*g100_100
      g020_100 =            py*g010_100
      g011_100 =            pz*g010_100
      g002_100 =            pz*g001_100
      g200_010 =            px*g100_010
      g110_010 = g100_000 + py*g100_010
      g101_010 =            pz*g100_010
      g020_010 = g010_000 + py*g010_010
      g011_010 =            pz*g010_010
      g002_010 =            pz*g001_010
      g200_001 =            px*g100_001
      g110_001 =            py*g100_001
      g101_001 = g100_000 + pz*g100_001
      g020_001 =            py*g010_001
      g011_001 = g010_000 + pz*g010_001
      g002_001 = g001_000 + pz*g001_001
    else
!     p
      gao(1)  = g100_000
      gao(2)  = g010_000
      gao(3)  = g001_000
      gaox(1) = g100_100
      gaox(2) = g010_100
      gaox(3) = g001_100
      gaoy(1) = g100_010
      gaoy(2) = g010_010
      gaoy(3) = g001_010
      gaoz(1) = g100_001
      gaoz(2) = g010_001
      gaoz(3) = g001_001
      return
    end if

    if (cartesian_deg(ishell) > 6) then
      g300_000 =            px*g200_000
      g210_000 =            py*g200_000
      g201_000 =            pz*g200_000
      g120_000 =            py*g110_000
      g111_000 =            pz*g110_000
      g102_000 =            pz*g101_000
      g030_000 =            py*g020_000
      g021_000 =            pz*g020_000
      g012_000 =            pz*g011_000
      g003_000 =            pz*g002_000
      g300_100 = g200_000 + px*g200_100
      g210_100 =            py*g200_100
      g201_100 =            pz*g200_100
      g120_100 =            py*g110_100
      g111_100 =            pz*g110_100
      g102_100 =            pz*g101_100
      g030_100 =            py*g020_100
      g021_100 =            pz*g020_100
      g012_100 =            pz*g011_100
      g003_100 =            pz*g002_100
      g300_010 =            px*g200_010
      g210_010 = g200_000 + py*g200_010
      g201_010 =            pz*g200_010
      g120_010 = g110_000 + py*g110_010
      g111_010 =            pz*g110_010
      g102_010 =            pz*g101_010
      g030_010 = g020_000 + py*g020_010
      g021_010 =            pz*g020_010
      g012_010 =            pz*g011_010
      g003_010 =            pz*g002_010
      g300_001 =            px*g200_001
      g210_001 =            py*g200_001
      g201_001 = g200_000 + pz*g200_001
      g120_001 =            py*g110_001
      g111_001 = g110_000 + pz*g110_001
      g102_001 = g101_000 + pz*g101_001
      g030_001 =            py*g020_001
      g021_001 = g020_000 + pz*g020_001
      g012_001 = g011_000 + pz*g011_001
      g003_001 = g002_000 + pz*g002_001
    else
!     d
      gao(1)  = g200_000
      gao(2)  = g110_000
      gao(3)  = g101_000
      gao(4)  = g020_000
      gao(5)  = g011_000
      gao(6)  = g002_000
      gaox(1) = g200_100
      gaox(2) = g110_100
      gaox(3) = g101_100
      gaox(4) = g020_100
      gaox(5) = g011_100
      gaox(6) = g002_100
      gaoy(1) = g200_010
      gaoy(2) = g110_010
      gaoy(3) = g101_010
      gaoy(4) = g020_010
      gaoy(5) = g011_010
      gaoy(6) = g002_010
      gaoz(1) = g200_001
      gaoz(2) = g110_001
      gaoz(3) = g101_001
      gaoz(4) = g020_001
      gaoz(5) = g011_001
      gaoz(6) = g002_001
      return
    end if

    if (cartesian_deg(ishell) > 10) then
      g400_000 =            px*g300_000
      g310_000 =            py*g300_000
      g301_000 =            pz*g300_000
      g220_000 =            py*g210_000
      g211_000 =            pz*g210_000
      g202_000 =            pz*g201_000
      g130_000 =            py*g120_000
      g121_000 =            pz*g120_000
      g112_000 =            pz*g111_000
      g103_000 =            pz*g102_000
      g040_000 =            py*g030_000
      g031_000 =            pz*g030_000
      g022_000 =            pz*g021_000
      g013_000 =            pz*g012_000
      g004_000 =            pz*g003_000
      g400_100 = g300_000 + px*g300_100
      g310_100 =            py*g300_100
      g301_100 =            pz*g300_100
      g220_100 =            py*g210_100
      g211_100 =            pz*g210_100
      g202_100 =            pz*g201_100
      g130_100 =            py*g120_100
      g121_100 =            pz*g120_100
      g112_100 =            pz*g111_100
      g103_100 =            pz*g102_100
      g040_100 =            py*g030_100
      g031_100 =            pz*g030_100
      g022_100 =            pz*g021_100
      g013_100 =            pz*g012_100
      g004_100 =            pz*g003_100
      g400_010 =            px*g300_010
      g310_010 = g300_000 + py*g300_010
      g301_010 =            pz*g300_010
      g220_010 = g210_000 + py*g210_010
      g211_010 =            pz*g210_010
      g202_010 =            pz*g201_010
      g130_010 = g120_000 + py*g120_010
      g121_010 =            pz*g120_010
      g112_010 =            pz*g111_010
      g103_010 =            pz*g102_010
      g040_010 = g030_000 + py*g030_010
      g031_010 =            pz*g030_010
      g022_010 =            pz*g021_010
      g013_010 =            pz*g012_010
      g004_010 =            pz*g003_010
      g400_001 =            px*g300_001
      g310_001 =            py*g300_001
      g301_001 = g300_000 + pz*g300_001
      g220_001 =            py*g210_001
      g211_001 = g210_000 + pz*g210_001
      g202_001 = g201_000 + pz*g201_001
      g130_001 =            py*g120_001
      g121_001 = g120_000 + pz*g120_001
      g112_001 = g111_000 + pz*g111_001
      g103_001 = g102_000 + pz*g102_001
      g040_001 =            py*g030_001
      g031_001 = g030_000 + pz*g030_001
      g022_001 = g021_000 + pz*g021_001
      g013_001 = g012_000 + pz*g012_001
      g004_001 = g003_000 + pz*g003_001
    else
!     f
      gao( 1)  = g300_000
      gao( 2)  = g210_000
      gao( 3)  = g201_000
      gao( 4)  = g120_000
      gao( 5)  = g111_000
      gao( 6)  = g102_000
      gao( 7)  = g030_000
      gao( 8)  = g021_000
      gao( 9)  = g012_000
      gao(10)  = g003_000
      gaox( 1) = g300_100
      gaox( 2) = g210_100
      gaox( 3) = g201_100
      gaox( 4) = g120_100
      gaox( 5) = g111_100
      gaox( 6) = g102_100
      gaox( 7) = g030_100
      gaox( 8) = g021_100
      gaox( 9) = g012_100
      gaox(10) = g003_100
      gaoy( 1) = g300_010
      gaoy( 2) = g210_010
      gaoy( 3) = g201_010
      gaoy( 4) = g120_010
      gaoy( 5) = g111_010
      gaoy( 6) = g102_010
      gaoy( 7) = g030_010
      gaoy( 8) = g021_010
      gaoy( 9) = g012_010
      gaoy(10) = g003_010
      gaoz( 1) = g300_001
      gaoz( 2) = g210_001
      gaoz( 3) = g201_001
      gaoz( 4) = g120_001
      gaoz( 5) = g111_001
      gaoz( 6) = g102_001
      gaoz( 7) = g030_001
      gaoz( 8) = g021_001
      gaoz( 9) = g012_001
      gaoz(10) = g003_001
      return
    end if

    if (cartesian_deg(ishell) > 15) then
    else
!     g
      gao( 1)  = g400_000
      gao( 2)  = g310_000
      gao( 3)  = g301_000
      gao( 4)  = g220_000
      gao( 5)  = g211_000
      gao( 6)  = g202_000
      gao( 7)  = g130_000
      gao( 8)  = g121_000
      gao( 9)  = g112_000
      gao(10)  = g103_000
      gao(11)  = g040_000
      gao(12)  = g031_000
      gao(13)  = g022_000
      gao(14)  = g013_000
      gao(15)  = g004_000
      gaox( 1) = g400_100
      gaox( 2) = g310_100
      gaox( 3) = g301_100
      gaox( 4) = g220_100
      gaox( 5) = g211_100
      gaox( 6) = g202_100
      gaox( 7) = g130_100
      gaox( 8) = g121_100
      gaox( 9) = g112_100
      gaox(10) = g103_100
      gaox(11) = g040_100
      gaox(12) = g031_100
      gaox(13) = g022_100
      gaox(14) = g013_100
      gaox(15) = g004_100
      gaoy( 1) = g400_010
      gaoy( 2) = g310_010
      gaoy( 3) = g301_010
      gaoy( 4) = g220_010
      gaoy( 5) = g211_010
      gaoy( 6) = g202_010
      gaoy( 7) = g130_010
      gaoy( 8) = g121_010
      gaoy( 9) = g112_010
      gaoy(10) = g103_010
      gaoy(11) = g040_010
      gaoy(12) = g031_010
      gaoy(13) = g022_010
      gaoy(14) = g013_010
      gaoy(15) = g004_010
      gaoz( 1) = g400_001
      gaoz( 2) = g310_001
      gaoz( 3) = g301_001
      gaoz( 4) = g220_001
      gaoz( 5) = g211_001
      gaoz( 6) = g202_001
      gaoz( 7) = g130_001
      gaoz( 8) = g121_001
      gaoz( 9) = g112_001
      gaoz(10) = g103_001
      gaoz(11) = g040_001
      gaoz(12) = g031_001
      gaoz(13) = g022_001
      gaoz(14) = g013_001
      gaoz(15) = g004_001
      return
    end if

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

    g000_000 = 0.0d0
    do i = 1, nr_nonzero_contraction_coefs(ishell)
      g000_000 = g000_000 + contraction_coef(ishell, i)*dexp(-primitive_exp(ishell, i)*p2)
    end do

    if (cartesian_deg(ishell) < 21) then
      call get_cartesian_ao(ishell, &
                            gao,    &
                            px, py, pz, p2)
    else
!     h and higher
      do ic = 1, cartesian_deg(ishell)
        jc = cartesian_off(ishell) + ic
        gao(ic) = (px**l_exp(jc))*(py**m_exp(jc))*(pz**n_exp(jc))*g000_000
      end do

    end if

!   if spherical and d or higher then transform
    if (shell_is_spherical(ishell) .and. (cartesian_deg(ishell) > 3)) then
       call transform_c_to_s(ishell, gao)
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
      call get_cartesian_ao1(ishell, &
                          gao,    &
                          gaox,   &
                          gaoy,   &
                          gaoz,   &
                          px, py, pz, p2)


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

!   if spherical and d or higher then transform
    if (shell_is_spherical(ishell) .and. (cartesian_deg(ishell) > 3)) then
       call transform_c_to_s(ishell, gao)
       call transform_c_to_s(ishell, gaox)
       call transform_c_to_s(ishell, gaoy)
       call transform_c_to_s(ishell, gaoz)
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
