module input_reader

! radovan: - this module should not initialize anything
!            it is a reader
!          - when we have only modules and no common blocks the initialization
!            problem goes away since variables can be initialized at definition

  use keyword
  use character_processing
  use dirac_cfg

  implicit none

  public read_input
  public move_to_next_star

  private

contains

  subroutine set_lscale_dft_gaunt()

#include "dcbham.h"

    lscale_dft_gaunt = .true.

  end subroutine

  subroutine read_input()

!   radovan: - this routine should not include any common blocks
!              this should be done in called subroutines
!            - keep ifndef traps inside subroutines for better readability
!              of this subroutine

!   ----------------------------------------------------------------------------
    character(kw_length) :: word
    character(kw_length) :: kw_section
!   ----------------------------------------------------------------------------

!   this is to catch keywords that appear before some section starts
    kw_section = '       '

    open(unit_in, file='DIRAC.INP')
    rewind unit_in

    do while (.true.)

      read(unit_in, '(a7)', end=1) word

      select case (word(1:1))

        case ('!', '#')
!         ignore line

        case ('*')
!         it is a section
          kw_section = uppercase(word)

        case default
!         it is a keyword
          select case (kw_section)

!           case ('**METHO', '**WAVE ')
!             call read_input_method(word, kw_section)

            case ('*DFT   ')
              call read_input_dft(word, kw_section)

            case ('*END OF')
              go to 1

            case ('       ')
              call quit('keyword '//word//' without section')

            case default
!             activate as soon as everything is merged to here
!             call quit('section '//kw_section//' not recognized')

          end select

      end select

    end do

1   close(unit_in)

  end subroutine

  subroutine move_to_next_star(word_io)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(inout) :: word_io
!   ----------------------------------------------------------------------------
    character(kw_length)                :: word
!   ----------------------------------------------------------------------------

    do while (.true.)
      read(unit_in, '(a7)') word
      if (word(1:1) == '*') then
        word_io = word
        return
      end if
    end do

  end subroutine

  subroutine read_input_dft(word, kw_section)

    use dft_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    character(80)                    :: line
    character(4)                     :: first_4
!   ----------------------------------------------------------------------------

#ifndef MOD_DFT
    call quit(kw_section//' not available in this version')
#endif

    if (.not. dirac_cfg_dft_calculation) return

    call reset_available_kw_list()

    if (kw_matches(word, '.IMPORT')) then
      dft_cfg_import_grid = .true.
      call kw_read(word, dft_cfg_gridfile)
    end if

    if (kw_matches(word, '.EXPORT')) then
      dft_cfg_export_grid = .true.
      call kw_read(word, dft_cfg_gridfile)
    end if

    if (kw_matches(word, '.FROMVX')) then
      call quit('explicit xc potential (.FROMVX) not available')
    end if

    if (kw_matches(word, '.RADINT')) then
      call kw_read(word, dft_cfg_radint)
    end if

    if (kw_matches(word, '.ANGMIN')) then
      call kw_read(word, dft_cfg_angmin)
    end if

    if (kw_matches(word, '.ANGINT')) then
      call kw_read(word, dft_cfg_angint)
    end if

    if (kw_matches(word, '.COARSE')) then
      dft_cfg_radint = 1.0d-11
      dft_cfg_angint = 35
    end if

    if (kw_matches(word, '.DEBUG ')) then
      dft_cfg_radint = 1.0d-3
      dft_cfg_angint = 10
    end if

    if (kw_matches(word, '.PRINT ')) then
      call kw_read(word, dft_cfg_print_level)
    end if

    if (kw_matches(word, '.INTCHK')) then
      call kw_read(word, dft_cfg_integration_check_level)
    end if

    if (kw_matches(word, '.NOPRUN')) then
!     turn off pruning of angular grid
      dft_cfg_no_pruning = .true.
    end if

    if (kw_matches(word, '.ATSIZE')) then
!     estimate relative atomic sizes for use in the Becke
!     partitioning scheme from atomic contributions
      dft_cfg_estimate_radii = .true.
    end if

    if (kw_matches(word, '.ERRELS')) then
!     threshold for accepted error in integrated number of electrons
      call kw_read(word, dft_cfg_accepted_error)
    end if

    if (kw_matches(word, '.GRAC  ')) then
!     gradient regulated asymptotic correction
      dft_cfg_grac = .true.
      read(unit_in, '(a4)') first_4
      if      (lowercase(first_4) == 'lb94') then
        dft_cfg_asymptote_is_lb94 = .true.
      else if (lowercase(first_4) == 'lbal') then
        dft_cfg_asymptote_is_lbalpha = .true.
      else
        call quit('keyword following .GRAC not recognized')
      end if
      call kw_read(word,               &
                   dft_cfg_grac_alpha, &
                   dft_cfg_grac_beta,  &
                   dft_cfg_ac_ip,      &
                   dft_cfg_ac_threshold)
    end if

    if (kw_matches(word, '.SAOP  ')) then
!     statistical averaging of (model) orbital potentials
      dft_cfg_saop = .true.
      read(unit_in, '(a4)') first_4
      if      (lowercase(first_4) == 'lb94') then
        dft_cfg_asymptote_is_lb94 = .true.
      else if (lowercase(first_4) == 'lbal') then
        dft_cfg_asymptote_is_lbalpha = .true.
      else
        call quit('keyword following .SAOP not recognized')
      end if
    end if

    if (kw_matches(word, '.SAOP! ')) then
!     "original" SAOP as in JCP 112, 1344 (2000).
      dft_cfg_saop                    = .true.
      dft_cfg_saop_with_response_part = .true.
      dft_cfg_asymptote_is_lbalpha    = .true.
      dft_cfg_alda_hs                 = .true.
      dft_cfg_alda_ha                 = .true.
      call setaldahs(.true.) !activate alda for the h+ part
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.SDFT  ')) then
      dft_cfg_sdft = .true.
    end if

    if (kw_matches(word, '.COLLIN')) then
      dft_cfg_sdft_collinear = .true.
    end if

    if (kw_matches(word, '.BETA  ')) then
      dft_cfg_sdft_beta = .true.
    end if

    if (kw_matches(word, '.VRG   ')) then
      dft_cfg_cdft_vrg = .true.
    end if

    if (kw_matches(word, '.ZIPGRI')) then
      dft_cfg_zipgrid = .true.
    end if

    if (kw_matches(word, '.GAUTSC')) then
      call set_lscale_dft_gaunt()
    end if

    if (kw_matches(word, '.ALDA  ')) then
      dft_cfg_alda_hs = .true.
      dft_cfg_alda_ha = .true.
      call setaldahs(.true.) !activate alda for the h+ part
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.ALDA+ ')) then
!     use alda only for the      hermitian part
      dft_cfg_alda_hs = .true.
      call setaldahs(.true.) !activate alda for the h+ part
    end if

    if (kw_matches(word, '.ALDA- ')) then
!     use alda only for the anti-hermitian part
      dft_cfg_alda_ha = .true.
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.XALDA ')) then
!     use        (1-f)*S + VWN + f*HFx
!     instead of       S + VWN
      dft_cfg_xalda   = .true.
      dft_cfg_alda_hs = .true.
      dft_cfg_alda_ha = .true.
      call setxalda(.true.)
      call setaldahs(.true.) !activate alda for the h+ part
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.XALDA+')) then
!     use xalda only for the      hermitian part
      dft_cfg_xalda   = .true.
      dft_cfg_alda_hs = .true.
      call setxalda(.true.)
      call setaldahs(.true.) !activate alda for the h+ part
    end if

    if (kw_matches(word, '.XALDA-')) then
!     use xalda only for the anti-hermitian part
      dft_cfg_xalda   = .true.
      dft_cfg_alda_ha = .true.
      call setxalda(.true.)
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.THRESH')) then
      call kw_read(word,           &
                   dft_cfg_dfthr0, &
                   dft_cfg_dfthrl, &
                   dft_cfg_dfthri)
    end if

    if (kw_matches(word, '.4CGRID')) then
!     include small component in the grid generation
!     also for 1- and 2-c dft calculations
      dft_cfg_force_4c_grid = .true.
    end if

    if (kw_matches(word, '.KEEPAO')) then
      dft_cfg_keep_ao = .true.
    end if

    if (kw_matches(word, '.OLDDFT')) then
      dft_cfg_old_dft = .true.
    end if

    if (kw_matches(word, '.GRID T')) then
      call kw_read(word, line)
!     call dftgridinput(line)
    end if

    call check_whether_kw_found(word, kw_section)

  end subroutine

end module
