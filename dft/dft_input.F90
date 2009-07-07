module dft_input

  use dft_cfg
  use character_processing

  implicit none

#ifdef MOD_AUTODIFF
  public parse_functional
#endif
  public parse_dft_input
  public consistency_after_dft_input
  public report_after_dft_input

  private

  real(8), parameter :: tiny = 1.0d-20

contains

  subroutine reset_hfx()

#include "dcbham.h"

    hfxfac = 0.0d0

  end subroutine

  subroutine add_hfx(w)

!   ----------------------------------------------------------------------------
    real(8), intent(in) :: w
!   ----------------------------------------------------------------------------

#include "dcbham.h"

    if (dabs(w) > tiny) then
      hfxfac = hfxfac + w
    end if

  end subroutine

  subroutine add_weight(i, w, is_gga)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: i
    real(8), intent(in) :: w
    logical, intent(in) :: is_gga
!   ----------------------------------------------------------------------------

    if (dabs(w) > tiny) then
      if (is_gga) then
        dft_cfg_gga_functional_weight_array(i) = dft_cfg_gga_functional_weight_array(i) + w
      else
        dft_cfg_lda_functional_weight_array(i) = dft_cfg_lda_functional_weight_array(i) + w
      end if
    end if

  end subroutine

  function get_weight(word, functional)

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: word
    character(*), intent(in)  :: functional
    real(8)                   :: get_weight
!   ----------------------------------------------------------------------------
    integer                   :: i
!   ----------------------------------------------------------------------------

    get_weight = 0.0d0

    if (word_contains(word, '=')) then

      i = index(word, '=')

      if (word(1:i-1) == functional) then
        read(word(i+1:len(word)), *) get_weight
      end if

    else

      if (word == functional) then
        get_weight = 1.0d0
      end if

    end if

  end function

#ifdef MOD_AUTODIFF
  subroutine parse_functional(line)

#include "../autodiff/xc_fortran_interface.h"

!   ----------------------------------------------------------------------------
    character(80)              :: line
!   ----------------------------------------------------------------------------
    character(80), allocatable :: word_array(:)
    integer                    :: i, nr_words
    real(8)                    :: w
!   ----------------------------------------------------------------------------

    dft_cfg_lda_functional_weight_array = 0.0d0
    dft_cfg_gga_functional_weight_array = 0.0d0

    call reset_hfx()

    call lowercase(line)

    nr_words = word_count(line)

    allocate(word_array(nr_words))
    read(line, *) (word_array(i), i = 1, nr_words)

    do i = 1, nr_words


!     first check aliases
!     ===================

!     lda
      w = get_weight(word_array(i), 'lda')
      call add_weight(i=LDA_SLATER, w=w, is_gga=.false.)
      call add_weight(i=LDA_VWN,    w=w, is_gga=.false.)

!     svwn = lda
      w = get_weight(word_array(i), 'svwn')
      call add_weight(i=LDA_SLATER, w=w, is_gga=.false.)
      call add_weight(i=LDA_VWN,    w=w, is_gga=.false.)

!     blyp
      w = get_weight(word_array(i), 'blyp')
      call add_weight(i=LDA_SLATER, w=w, is_gga=.false.)
      call add_weight(i=GGA_BECKE,  w=w, is_gga=.true.)
      call add_weight(i=GGA_LYP,    w=w, is_gga=.true.)

!     b3lyp
      w = get_weight(word_array(i), 'b3lyp')
      call add_weight(i=LDA_SLATER, w=0.80d0*w, is_gga=.false.)
      call add_weight(i=LDA_VWN,    w=0.19d0*w, is_gga=.false.)
      call add_weight(i=GGA_BECKE,  w=0.72d0*w, is_gga=.true.)
      call add_weight(i=GGA_LYP,    w=0.81d0*w, is_gga=.true.)
      call add_hfx(0.20d0*w)

!     pbe
      w = get_weight(word_array(i), 'pbe')
      call add_weight(i=GGA_PBEX,   w=w, is_gga=.true.)
      call add_weight(i=GGA_PBEC,   w=w, is_gga=.true.)


!     now atomic functional acronyms
!     ==============================

      w = get_weight(word_array(i), 'slater')
      call add_weight(i=LDA_SLATER, w=w, is_gga=.false.)

      w = get_weight(word_array(i), 'vwn')
      call add_weight(i=LDA_VWN,    w=w, is_gga=.false.)

      w = get_weight(word_array(i), 'becke')
      call add_weight(i=GGA_BECKE,  w=w, is_gga=.true.)

      w = get_weight(word_array(i), 'lyp')
      call add_weight(i=GGA_LYP,    w=w, is_gga=.true.)

      w = get_weight(word_array(i), 'pbex')
      call add_weight(i=GGA_PBEX,   w=w, is_gga=.true.)

      w = get_weight(word_array(i), 'pbec')
      call add_weight(i=GGA_PBEC,   w=w, is_gga=.true.)

      w = get_weight(word_array(i), 'hfx')
      call add_hfx(w)

    end do

    deallocate(word_array)

    do i = 1, size(dft_cfg_gga_functional_weight_array)
      if (dabs(dft_cfg_gga_functional_weight_array(i)) > tiny) then
        dft_cfg_is_gga_autodiff = .true.
      end if
    end do

  end subroutine
#endif /* ifdef MOD_AUTODIFF */

  function get_lucmd() result(u)

#include "priunit.h"

!   ----------------------------------------------------------------------------
    integer :: u
!   ----------------------------------------------------------------------------

    u = lucmd

  end function

  function get_lupri() result(u)

#include "priunit.h"

!   ----------------------------------------------------------------------------
    integer :: u
!   ----------------------------------------------------------------------------

    u = lupri

  end function

  subroutine set_lscale_dft_gaunt()

#include "dcbham.h"

    lscale_dft_gaunt = .true.

  end subroutine

  subroutine parse_dft_input(word_io)

!   radovan: if you add new keywords then it's
!            nice if the user gets the confirmation
!            that something is happening in the output
!            so write something in report_after_dft_input

!   ----------------------------------------------------------------------------
    character(7), intent(inout) :: word_io
!   ----------------------------------------------------------------------------
    character(7)                :: word
    character(80)               :: line
    character(1)                :: first
    character(4)                :: first_4
    integer                     :: u
!   ----------------------------------------------------------------------------

    if (word_io /= '*DFT   ') return

    u = get_lucmd()

    do while (.true.)

      read (u, '(a7)') word
      first = word(1:1)

      select case (first)
        case ('!')
!         ignore line
        case ('#')
!         ignore line
        case ('*')
          word_io = word
          return
        case ('.')

        select case (word)

          case ('.IMPORT')
            dft_cfg_import_grid = .true.
            read(u, '(a30)') dft_cfg_gridfile

          case ('.EXPORT')
            dft_cfg_export_grid = .true.
            read(u, '(a30)') dft_cfg_gridfile

          case ('.FROMVX')
            call quit('explicit xc potential (.FROMVX) not available')

          case ('.RADINT')
            read(u, *) dft_cfg_radint

          case ('.ANGMIN')
            read(u, *) dft_cfg_angmin

          case ('.ANGINT')
            read(u, *) dft_cfg_angint

          case ('.COARSE')
            dft_cfg_radint = 1.0d-11
            dft_cfg_angint = 35

          case ('.DEBUG')
            dft_cfg_radint = 1.0d-3
            dft_cfg_angint = 10

          case ('.PRINT')
            read(u, *) dft_cfg_print_level

          case ('.INTCHK')
            read(u, *) dft_cfg_integration_check_level

          case ('.NOPRUN')
!           turn off pruning of angular grid
            dft_cfg_no_pruning = .true.

          case ('.ATSIZE')
!           estimate relative atomic sizes for use in the Becke
!           partitioning scheme from atomic contributions
            dft_cfg_estimate_radii = .true.

          case ('.ERRELS')
!           threshold for accepted error in integrated number of electrons
            read(u, *) dft_cfg_accepted_error

          case ('.GRAC  ')
!           gradient regulated asymptotic correction
            dft_cfg_grac = .true.
            read(u, '(a4)') first_4
            if      (first_4 == 'LB94') then
              dft_cfg_asymptote_is_lb94 = .true.
            else if (first_4 == 'LBal') then
              dft_cfg_asymptote_is_lbalpha = .true.
            else if (first_4 == 'LBAl') then
              dft_cfg_asymptote_is_lbalpha = .true.
            else if (first_4 == 'LBAL') then
              dft_cfg_asymptote_is_lbalpha = .true.
            else
              call quit('keyword following .GRAC not recognized')
            end if
            read(u, *) dft_cfg_grac_alpha, &
                       dft_cfg_grac_beta,  &
                       dft_cfg_ac_ip,      &
                       dft_cfg_ac_threshold

          case ('.SAOP  ')
!           statistical averaging of (model) orbital potentials
            dft_cfg_saop = .true.
            read(u, '(a4)') first_4
            if      (first_4 == 'LB94') then
              dft_cfg_asymptote_is_lb94 = .true.
            else if (first_4 == 'LBal') then
              dft_cfg_asymptote_is_lbalpha = .true.
            else if (first_4 == 'LBAl') then
              dft_cfg_asymptote_is_lbalpha = .true.
            else if (first_4 == 'LBAL') then
              dft_cfg_asymptote_is_lbalpha = .true.
            else
              call quit('keyword following .SAOP not recognized')
            end if

          case ('.SAOP! ')
!           "original" SAOP as in JCP 112, 1344 (2000).
            dft_cfg_saop                    = .true.
            dft_cfg_saop_with_response_part = .true.
            dft_cfg_asymptote_is_lbalpha    = .true.
            dft_cfg_alda_hs                 = .true.
            dft_cfg_alda_ha                 = .true.
            call setaldahs(.true.) !activate alda for the h+ part
            call setaldaha(.true.) !activate alda for the h- part

          case ('.SDFT  ')
            dft_cfg_sdft = .true.

          case ('.COLLIN')
            dft_cfg_sdft_collinear = .true.

          case ('.BETA  ')
            dft_cfg_sdft_beta = .true.

          case ('.VRG   ')
            dft_cfg_cdft_vrg = .true.

          case ('.ZIPGRI')
            dft_cfg_zipgrid = .true.

          case ('.GAUTSC')
            call set_lscale_dft_gaunt()

          case ('.ALDA  ')
            dft_cfg_alda_hs = .true.
            dft_cfg_alda_ha = .true.
            call setaldahs(.true.) !activate alda for the h+ part
            call setaldaha(.true.) !activate alda for the h- part

          case ('.ALDA+ ')
!           use alda only for the      hermitian part
            dft_cfg_alda_hs = .true.
            call setaldahs(.true.) !activate alda for the h+ part

          case ('.ALDA- ')
!           use alda only for the anti-hermitian part
            dft_cfg_alda_ha = .true.
            call setaldaha(.true.) !activate alda for the h- part

          case ('.XALDA ')
!           use        (1-f)*S + VWN + f*HFx
!           instead of       S + VWN
            dft_cfg_xalda   = .true.
            dft_cfg_alda_hs = .true.
            dft_cfg_alda_ha = .true.
            call setxalda(.true.)
            call setaldahs(.true.) !activate alda for the h+ part
            call setaldaha(.true.) !activate alda for the h- part

          case ('.XALDA+')
!           use xalda only for the      hermitian part
            dft_cfg_xalda   = .true.
            dft_cfg_alda_hs = .true.
            call setxalda(.true.)
            call setaldahs(.true.) !activate alda for the h+ part

          case ('.XALDA-')
!           use xalda only for the anti-hermitian part
            dft_cfg_xalda   = .true.
            dft_cfg_alda_ha = .true.
            call setxalda(.true.)
            call setaldaha(.true.) !activate alda for the h- part

          case ('.THRESH')
            read(u, *) dft_cfg_dfthr0, &
                       dft_cfg_dfthrl, &
                       dft_cfg_dfthri

          case ('.4CGRID')
!           include small component in the grid generation
!           also for 1- and 2-c dft calculations
            dft_cfg_force_4c_grid = .true.

          case ('.KEEPAO')
            dft_cfg_keep_ao = .true.

          case ('.OLDDFT')
            dft_cfg_old_dft = .true.

          case ('.GRID T')
            read(u, '(a80)') line
!           call dftgridinput(line)

          case default
            call quit('illegal keyword after '//word_io//': '//word)
        end select

        case default
          call quit('illegal keyword after '//word_io//': '//word)
      end select

    end do

  end subroutine

  subroutine consistency_after_dft_input()

    if (dft_cfg_asymptote_is_lb94 .and. dft_cfg_asymptote_is_lbalpha) then
      call quit('LBalpha = LB94 = .true.')
    end if

    if (dft_cfg_sdft_collinear .and. .not. dft_cfg_sdft) then
      call quit('please use .COLLINEAR together with .SDFT')
    end if

  end subroutine

  subroutine report_after_dft_input()

!   ----------------------------------------------------------------------------
    integer :: u
!   ----------------------------------------------------------------------------

    u = get_lupri()

    call header('Kohn-Sham calculation set-up', 2)

!   print level

    write(u, '(a, i5)')      ' * Print level: ', dft_cfg_print_level

!   integration grid

    write(u, '(a)')          ' * Numerical integration grid:'

    if (dft_cfg_import_grid) then

      write(u, '(a)')        '   - integration grid imported from file '//dft_cfg_gridfile

    else

      write(u, '(a)')        '   - radial quadrature according to '// &
                             'R. Lindh, P.-Aa. Malmqvist, and L. Gagliardi'
      write(u, '(a, e12.5)') '     precision of radial quadrature set to: ', dft_cfg_radint
      write(u, '(a, i2)')    '   - angular quadrature using the Lebedev scheme, '// &
                             'exact up to order L = ', dft_cfg_angint

      if (dft_cfg_no_pruning) then
        write(u, '(a)')      '   - pruning of angular grid turned off'
      end if

      if (dft_cfg_integration_check_level > 0) then
        write(u, '(a)')      '   - numerical integration will be tested'
      end if

      write(u, '(a, e12.5)') '   - accepted error in integrated number '// &
                             'of electrons: ', dft_cfg_accepted_error

      if (dft_cfg_estimate_radii) then
        write(u, '(a)')      '   - estimate relative atomic sizes for use in the Becke'
        write(u, '(a)')      '     partitioning scheme from atomic contributions '// &
                             'to the density'
      end if

    end if !if (dft_cfg_import_grid)

    if (dft_cfg_export_grid) then
      write(u, '(a)')        '   - integration grid exported to file'
    end if

!   thresholds

    write(u, '(a)')          ' * DFT thresholds:'
    write(u, '(a, e12.5)')   '   - dfthr0 = ', dft_cfg_dfthr0
    write(u, '(a, e12.5)')   '   - dfthrl = ', dft_cfg_dfthrl
    write(u, '(a, e12.5)')   '   - dfthri = ', dft_cfg_dfthri

!   asymptotic correction

    if (dft_cfg_grac) then
      write(u, '(a)')        ' * Gradient regulated asymptotic correction:'
      write(u, '(a, e12.5)') '   - alpha = ', dft_cfg_grac_alpha
      write(u, '(a, e12.5)') '   - beta  = ', dft_cfg_grac_beta
      write(u, '(a)')        '   - threshold for difference in HOMO eigenvalue'
      write(u, '(a, e12.5)') '     below which asymptotic correction '// &
                             'is switched on: ', dft_cfg_ac_threshold
      write(u, '(a, e12.5)') '   - bulk potential will be shifted by (IP): ', dft_cfg_ac_ip
    end if

    if (dft_cfg_saop) then
      write(u, '(a)')        ' * Statistical averaging of (model) orbital potentials:'
    end if

    if (dft_cfg_saop_with_response_part) then
      write(u, '(a)')        '   - original SAOP as in JCP 112, 1344 (2000) chosen'
      write(u, '(a)')        '   - this activates ALDA xc kernel'
      write(u, '(a)')        '   - the functional under .DFT is expected to be GLLBhole'
    end if

    if (dft_cfg_asymptote_is_lb94) then
      write(u, '(a)')        '   - the asymptotic potential is LB94'
    end if

    if (dft_cfg_asymptote_is_lbalpha) then
      write(u, '(a)')        '   - the asymptotic potential is LBalpha'
    end if

!   sdft

    if (dft_cfg_sdft) then

      write(u, '(a)')        ' * Spin density contribution in response calculations:'

      if (dft_cfg_sdft_collinear) then
        write(u, '(a)')      '   - use the collinear approximation as '// &
                             'a definition of the spin density'
      else
        write(u, '(a)')      '   - use the norm of the spin magnetization vector as '// &
                             'a definition of the spin density (noncollinear definition)'
      end if

    else

      write(u, '(a)')        ' * No spin density contribution in response calculations.'
      write(u, '(a)')        '     (if you want the spin density contribution, use .SDFT)'

    end if !if (sdft) then

!   alda

    if (dft_cfg_alda_hs .or. dft_cfg_alda_ha) then

      write(u, '(a)')        ' * Adiabatic local density approximation (ALDA):'
      write(u, '(a)')        '   - approximate all functional derivatives '// &
                             'beyond the xc potential by SVWN derivatives'

      if (dft_cfg_alda_hs .and. .not. dft_cfg_alda_ha) then
        write(u, '(a)')      '   - use ALDA only for the Hermitian part '// &
                             '(density contribution)'
        write(u, '(a)')      '   - use the proper xc kernel for the anti-Hermitian '// &
                             'part (spin density contribution)'
      end if

      if (dft_cfg_alda_ha .and. .not. dft_cfg_alda_hs) then
        write(u, '(a)')      '   - use ALDA only for the anti-Hermitian part '// &
                             '(spin density contribution)'
        write(u, '(a)')      '   - use the proper xc kernel for the Hermitian '// &
                             'part (density contribution)'
      end if

      if (dft_cfg_xalda) then
        write(u, '(a)')      '   - keep the fraction of exact exchange of the '// &
                             'xc functional in the solution of the response equation'
      else
        write(u, '(a)')      '   - for hybrid functionals exact exchange is '// &
                             'switched off in the solution of the response equation'
        write(u, '(a)')      '     (if you want to keep it, use .XALDA)'
      end if

    end if !if (dft_cfg_alda_hs .or. dft_cfg_alda_ha) then

  end subroutine

end module
