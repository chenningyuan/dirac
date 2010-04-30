module dft_ao

!   this is based on DIRAC and Dalton
!   routines written by T. Helgaker, P. Salek, T. Saue, and others

! - in order to calculate aos using this module
!   you need to initialize some data first
!   by calling dft_ao_init from each node
! - dft_ao_init will set symmetry and basis set information
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
!   this is why some information is passed through dft_ao_init and not queried
!   from dft_cfg, to keep this module generic, please dont "use dft_cfg" here

!   some "bad ideas" = ideas which did not improve things
!   so that future generations do not waste time on this:
! - ao to so transformation using dgemv: possible, but much slower than present scheme
! - derivatives of aos by numerical differentiation: again, turned out to be much slower

!   FAQ
!   ===
!
!   1. how is the ao array sorted on output of get_so?
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
  integer, public      :: nr_shells       = 0
  integer              :: dft_maxrep          = 0
  integer              :: dft_maxopr          = 0
  integer              :: nr_boson_sym        = 0

  real(8), allocatable :: primitive_exp(:, :)
  real(8), allocatable, public :: shell_center(:)
  real(8), allocatable :: contraction_coef(:, :)
  real(8), allocatable :: shell_extent(:)
  integer, allocatable :: nr_nonzero_contraction_coefs(:)
  integer, allocatable, public :: cartesian_deg(:)
  integer, allocatable :: cartesian_off(:)
  integer, allocatable, public :: spherical_deg(:)
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

  integer, parameter   :: max_nr_so_ao = 10000
  integer             , public :: nr_so_ao
  integer, allocatable, public :: i_so_ao(:, :)
  real(8), allocatable, public :: f_so_ao(:)

  logical              :: use_ao_screening = .false.
  logical              :: only_large       = .false.

  real(8), parameter   :: tiny = 1.0d-30

! this is probably eq. 8.2.25 in the pink book
  real(8), parameter   :: g(10) = (/1.0d0, 1.3333d0, 1.6d0, 1.83d0, 2.03d0, 2.22d0, 2.39d0, 2.55d0, 2.70d0, 2.84d0/)

  real(8), parameter   :: cell_size          = 2.25d0
! real(8), parameter   :: half_cell_diagonal = cell_size*0.5d0*sqrt(3.0d0)
! compiler complains about intrinsic in a parameter declaration
  real(8), parameter   :: half_cell_diagonal = cell_size*0.866025404d0

  logical :: dft_basis_info_allocated = .false.
  logical :: dft_basis_info_set       = .false.
#ifdef VAR_MPI
  logical :: dft_basis_info_synced    = .false.
#endif

#ifdef VAR_MPI
#include "infpar.h"
#include "mpif.h"
#endif

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

    if (use_ao_screening) call dzero(ao, total_ao_vector_length)

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
      do isym = 0, dft_maxopr
        if (iand(isym, dft_istbao(ishell)) == 0) then

          px = cx - dft_pt(iand(dft_isymax(1, 1), isym))*dft_cent(ishell, 1, 1)
          py = cy - dft_pt(iand(dft_isymax(2, 1), isym))*dft_cent(ishell, 2, 1)
          pz = cz - dft_pt(iand(dft_isymax(3, 1), isym))*dft_cent(ishell, 3, 1)

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
    real(8)                :: f
!   ----------------------------------------------------------------------------

    g000_000 = 0.0d0
    g000_100 = 0.0d0
    g000_010 = 0.0d0
    g000_001 = 0.0d0

    do i = 1, nr_nonzero_contraction_coefs(ishell)
      f = -primitive_exp(ishell, i)
      g000_000 = g000_000 +            contraction_coef(ishell, i)*dexp(f*p2)
      g000_100 = g000_100 + px*f*2.0d0*contraction_coef(ishell, i)*dexp(f*p2)
      g000_010 = g000_010 + py*f*2.0d0*contraction_coef(ishell, i)*dexp(f*p2)
      g000_001 = g000_001 + pz*f*2.0d0*contraction_coef(ishell, i)*dexp(f*p2)
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
    real(8), intent(in)  :: ao(total_ao_vector_length)
    real(8), intent(out) :: so(total_ao_vector_length)
!   ----------------------------------------------------------------------------
    integer              :: i, n, off
!   ----------------------------------------------------------------------------

    so = 0.0d0

    do n = 1, nr_ao_slices
      off = nr_ao*(n - 1)
      do i = 1, nr_so_ao
        so(off + i_so_ao(1, i)) = so(off + i_so_ao(1, i)) &
                                + ao(off + i_so_ao(2, i))*f_so_ao(i)
      end do
    end do

  end subroutine

  subroutine dft_ao_init(is_parallel_calculation, &
                         max_ao_g_order,          &
                         max_ao_m_order,          &
                         use_screening,           &
                         ao_threshold,            &
                         force_small_ao)

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

!   ----------------------------------------------------------------------------
    logical, intent(in) :: is_parallel_calculation
    integer, intent(in) :: max_ao_g_order
    integer, intent(in) :: max_ao_m_order
    logical, intent(in) :: use_screening
    real(8), intent(in) :: ao_threshold
    logical, intent(in) :: force_small_ao
!   ----------------------------------------------------------------------------
    integer             :: g_order
    integer             :: m_order
    integer             :: iprim, i, j, icount, i_c, i_s, contracted_shell_nr, ishell, isym

    integer              :: iao, iso,         nshell
    integer              :: ic,         iso_bs
    integer              ::    nr_ao_irep, mul, irep,       n, ivec, off
    real(8)              :: f, r
#ifdef VAR_MPI
    logical             :: i_am_the_master
#endif
!   ----------------------------------------------------------------------------

    use_ao_screening = use_screening

    if (.not. dft_basis_info_allocated) then

      call alloc(primitive_exp, mxshel, mxcont)
      call alloc(shell_center, mxshel)
      call alloc(contraction_coef, mxshel, mxcont)
      call alloc(shell_extent, mxshel)
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

!fixme: this is reset later but i need it here for nr_cartesian_ao
      dft_istbao = istbao
      dft_maxopr = maxopr

      nr_cartesian_ao = 0
      do ishell = 1, mxshel
        do isym = 0, dft_maxopr
          if (iand(isym, dft_istbao(ishell)) == 0) then
            nr_cartesian_ao = nr_cartesian_ao + kckt(ishell)
          end if
        end do
      end do


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

      call alloc(i_so_ao, 2, max_nr_so_ao*8)
      call alloc(f_so_ao, max_nr_so_ao*8)
      i_so_ao = 0
      f_so_ao = 0.0d0

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
      if (mc == 1) only_large = .true.
      if (only_large .and. force_small_ao) only_large = .false.
      nr_ao = ntbas(0)

      nr_shells       = kmax
      if (only_large) nr_shells = nlrgsh

      dft_csp    = csp
      dft_ispadr = ispadr

      dft_icos  = icos
      dft_ibbas = ibbas
      dft_nbbas = nbbas
      shell_center = ncent

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

!       get shell extents
        shell_extent(ishell) = 0.0d0

        icount = 0
        do iprim = i, j
          if (dabs(priccf(iprim, contracted_shell_nr)) > 1.0d-30) then
            icount = icount + 1
            nr_nonzero_contraction_coefs(ishell) = icount
            contraction_coef(ishell, icount)    = priccf(iprim, contracted_shell_nr)
            primitive_exp(ishell, icount)       = priexp(iprim)

!           get shell extents
            r = (log(abs(priccf(iprim, contracted_shell_nr))) &
               - log(ao_threshold))/priexp(iprim)

            if (shell_extent(ishell) < r) shell_extent(ishell) = r
          end if
        end do

        shell_extent(ishell) = sqrt(shell_extent(ishell))*g(n_quantum_nr(ishell)) &
                             + half_cell_diagonal

!       better square this once than sqrt p2 in each point
        shell_extent(ishell) = shell_extent(ishell)**2.0d0

      end do

!     ao to so business
      i = 1
      nr_so_ao = 0
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
!           is small
            ic = 2
          else
!           is large
            ic = 1
          end if

!         do not plan to AO-SO transform small component functions
!         in case of zero small-small metric (.LEVY-LEBLOND and .X2C4)
          if (.not. (only_large .and. (ic == 2))) then
      
!           shift by difference boson sorted minus boson unsorted
            iso_bs = iso + dft_ibbas(irep, ic) - dft_icos(irep + 1, ic)
           
            do isym = 0, dft_maxopr
              if (iand(isym, mul) == 0) then
           
                f = dft_pt(iand(isym, ieor(irep, dft_isymao(n, nshell))))
           
                nr_so_ao = nr_so_ao + 1

                if (nr_so_ao > max_nr_so_ao) then
                  call quit('increase max_nr_so_ao')
                end if
           
                i_so_ao(1, nr_so_ao) = iso_bs
                i_so_ao(2, nr_so_ao) = iao
                f_so_ao(nr_so_ao)    = f
           
                iao = iao + spherical_deg(ishell)
           
              endif
            end do

          end if

        end do
        i = i + nr_ao_irep
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

!   total_ao_vector_length = nr_ao_slices*nr_ao
!   radovan: the above can lead to "plus one" or "plus few" out-of-bounds
!            error if spherical gaussians are used
!            it is safe to do instead
    total_ao_vector_length = nr_ao_slices*nr_cartesian_ao
!            although it consumes more memory than needed
!            the allocation could be refined

#ifdef VAR_MPI
    end if !if (i_am_the_master) then
#endif

#ifdef VAR_MPI
    if (is_parallel_calculation .and. .not. dft_basis_info_synced) then

!     is done above
!     call dirac_mpi_bcast(nr_cartesian_ao             )
      call dirac_mpi_bcast(nr_shells               )
      call dirac_mpi_bcast(dft_maxrep                  )
      call dirac_mpi_bcast(dft_maxopr                  )
      call dirac_mpi_bcast(nr_boson_sym                )
      call dirac_mpi_bcast(nr_ao                       )

      call dirac_mpi_bcast(primitive_exp               )
      call dirac_mpi_bcast(shell_center                )
      call dirac_mpi_bcast(contraction_coef            )
      call dirac_mpi_bcast(shell_extent                )
      call dirac_mpi_bcast(nr_nonzero_contraction_coefs)
      call dirac_mpi_bcast(cartesian_deg               )
      call dirac_mpi_bcast(cartesian_off               )
      call dirac_mpi_bcast(spherical_deg               )
      call dirac_mpi_bcast(spherical_off               )
      call dirac_mpi_bcast(n_quantum_nr                )
      call dirac_mpi_bcast(is_spherical                )
      call dirac_mpi_bcast(dft_istbao                  )
      call dirac_mpi_bcast(dft_pt                      )
      call dirac_mpi_bcast(dft_cent                    )
      call dirac_mpi_bcast(dft_isymax                  )

      call dirac_mpi_bcast(dft_csp                     )
      call dirac_mpi_bcast(dft_ispadr                  )

      call dirac_mpi_bcast(l_exp                       )
      call dirac_mpi_bcast(m_exp                       )
      call dirac_mpi_bcast(n_exp                       )

      call dirac_mpi_bcast(nr_ao_slices    )
      call dirac_mpi_bcast(need_ao_order )
      call dirac_mpi_bcast(ao_start      )
      call dirac_mpi_bcast(ao_end        )
      call dirac_mpi_bcast(total_ao_vector_length)
      call dirac_mpi_bcast(only_large)

      call dirac_mpi_bcast(nr_so_ao)
      call dirac_mpi_bcast(i_so_ao)
      call dirac_mpi_bcast(f_so_ao)

      dft_basis_info_synced = .true.

    end if
#endif

  end subroutine

end module
