module dft_input

  use dft_cfg
  use character_processing

  implicit none

#ifdef MOD_AUTODIFF
  public parse_functional
#endif
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

    hfxfac = hfxfac + w

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

    use xc_fun

!   radovan: line is made lowercase, then cut into words
!            if word contains "=" the weight is read
!            otherwise it is set to 1.0
!
!            equivalent functionals:
!            lda
!            slater vwn
!            slater=1.0 vwn=1.0
!
!            equivalent functionals:
!            b3lyp
!            slater=0.8 becke=0.72 hfx=0.2 vwn=0.19 lyp=0.81
!            lda slater=-0.2 hfx=0.2 becke=0.72 vwn=-0.81 lyp=0.81

!   ----------------------------------------------------------------------------
    character(80)              :: line
!   ----------------------------------------------------------------------------
    character(80), allocatable :: word_array(:)
    integer                    :: i, nr_words
    real(8)                    :: w
    real(8)                    :: w_array(XC_PARAMS_LEN)
    integer                    :: gga_functionals(4) = (/XC_BECKE_EXCHANGE,  &
                                                         XC_LYP_CORRELATION, &
                                                         XC_PBE_EXCHANGE,    &
                                                         XC_PBE_CORRELATION/)
!   ----------------------------------------------------------------------------

    w_array = 0.0d0

    call reset_hfx()

    line = lowercase(line)

    nr_words = word_count(line)

    allocate(word_array(nr_words))
    read(line, *) (word_array(i), i = 1, nr_words)

    do i = 1, nr_words


!     composite functionals
!     =====================

!     lda
      w = get_weight(word_array(i), 'lda')
      if (dabs(w) > tiny) then
        w_array(XC_SLATER_EXCHANGE) = w_array(XC_SLATER_EXCHANGE) + w
        w_array(XC_VWN_CORRELATION) = w_array(XC_VWN_CORRELATION) + w
      end if

!     svwn = lda
      w = get_weight(word_array(i), 'svwn')
      if (dabs(w) > tiny) then
        w_array(XC_SLATER_EXCHANGE) = w_array(XC_SLATER_EXCHANGE) + w
        w_array(XC_VWN_CORRELATION) = w_array(XC_VWN_CORRELATION) + w
      end if

!     blyp
      w = get_weight(word_array(i), 'blyp')
      if (dabs(w) > tiny) then
        w_array(XC_SLATER_EXCHANGE) = w_array(XC_SLATER_EXCHANGE) + w
        w_array(XC_BECKE_EXCHANGE)  = w_array(XC_BECKE_EXCHANGE)  + w
        w_array(XC_LYP_CORRELATION) = w_array(XC_LYP_CORRELATION) + w
      end if

!     b3lyp
      w = get_weight(word_array(i), 'b3lyp')
      if (dabs(w) > tiny) then
        w_array(XC_SLATER_EXCHANGE) = w_array(XC_SLATER_EXCHANGE) + 0.80d0*w
        w_array(XC_VWN_CORRELATION) = w_array(XC_VWN_CORRELATION) + 0.19d0*w
        w_array(XC_BECKE_EXCHANGE)  = w_array(XC_BECKE_EXCHANGE)  + 0.72d0*w
        w_array(XC_LYP_CORRELATION) = w_array(XC_LYP_CORRELATION) + 0.81d0*w
        call add_hfx(0.20d0*w)
      end if

!     pbe
      w = get_weight(word_array(i), 'pbe')
      if (dabs(w) > tiny) then
        w_array(XC_PBE_EXCHANGE)    = w_array(XC_PBE_EXCHANGE)    + w
        w_array(XC_PBE_CORRELATION) = w_array(XC_PBE_CORRELATION) + w
      end if

!     pbe0
      w = get_weight(word_array(i), 'pbe0')
      if (dabs(w) > tiny) then
        w_array(XC_PBE_EXCHANGE)    = w_array(XC_PBE_EXCHANGE)    + 0.75d0*w
        w_array(XC_PBE_CORRELATION) = w_array(XC_PBE_CORRELATION) +        w
        call add_hfx(0.25d0*w)
      end if


!     "atomic" functionals
!     ====================

      w = get_weight(word_array(i), 'slater')
      if (dabs(w) > tiny) then
        w_array(XC_SLATER_EXCHANGE) = w_array(XC_SLATER_EXCHANGE) + w
      end if

      w = get_weight(word_array(i), 'vwn')
      if (dabs(w) > tiny) then
        w_array(XC_VWN_CORRELATION) = w_array(XC_VWN_CORRELATION) + w
      end if

      w = get_weight(word_array(i), 'becke')
      if (dabs(w) > tiny) then
        w_array(XC_BECKE_EXCHANGE) = w_array(XC_BECKE_EXCHANGE) + w
      end if

      w = get_weight(word_array(i), 'lyp')
      if (dabs(w) > tiny) then
        w_array(XC_LYP_CORRELATION) = w_array(XC_LYP_CORRELATION) + w
      end if

      w = get_weight(word_array(i), 'pbex')
      if (dabs(w) > tiny) then
        w_array(XC_PBE_EXCHANGE) = w_array(XC_PBE_EXCHANGE) + w
      end if

      w = get_weight(word_array(i), 'pbec')
      if (dabs(w) > tiny) then
        w_array(XC_PBE_CORRELATION) = w_array(XC_PBE_CORRELATION) + w
      end if

      w = get_weight(word_array(i), 'hfx')
      if (dabs(w) > tiny) then
        call add_hfx(w)
      end if

    end do

    deallocate(word_array)

    do i = 1, size(gga_functionals)
      if (dabs(w_array(gga_functionals(i))) > tiny) then
        dft_cfg_is_gga_autodiff = .true.
      end if
    end do

    call set_functional(w_array)

  end subroutine
#endif /* ifdef MOD_AUTODIFF */

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
