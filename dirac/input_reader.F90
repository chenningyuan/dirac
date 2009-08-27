module input_reader

  use character_processing

  implicit none

  public read_input
  public move_to_next_star

  integer, parameter, public :: kw_length = 7
  integer, parameter, public :: unit_in   = 5

  private

  interface kw_read
    module procedure kw_read_c
    module procedure kw_read_i1
    module procedure kw_read_r1
    module procedure kw_read_r2
    module procedure kw_read_r3
    module procedure kw_read_r4
  end interface

  integer, parameter         :: max_nr_kw   = 200
  integer, parameter         :: line_length = 80

  character(kw_length)       :: available_kw_list(max_nr_kw)
  integer                    :: nr_available_kw
  logical                    :: kw_found

contains

  subroutine set_lscale_dft_gaunt()

#include "dcbham.h"

    lscale_dft_gaunt = .true.

  end subroutine

  subroutine read_input()

!   ----------------------------------------------------------------------------
    character(kw_length) :: word
    character(kw_length) :: kw_section
!   ----------------------------------------------------------------------------

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
          kw_section = lowercase(word)

        case default
!         it is a keyword
          select case (kw_section)

!           later these calls can be moved out:
!           call reset_available_kw_list()
!           call check_whether_kw_found(word, kw_section)
!           but only when all sections are in this routine!

!-------------------------------------------------------------------------------
!           dft
            case ('*dft   ')
#ifndef MOD_DFT
              call quit(kw_section//' not available in this version')
#endif
              call reset_available_kw_list()
              call read_input_dft(word)
              call check_whether_kw_found(word, kw_section)
!-------------------------------------------------------------------------------

!           case ('       ')
!             call quit('keyword '//word//' without section')

            case default
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

  subroutine read_input_dft(word)

    use dft_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
!   ----------------------------------------------------------------------------
    character(80)                    :: line
    character(4)                     :: first_4
!   ----------------------------------------------------------------------------

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

  end subroutine

  subroutine kw_read_c(kw_input, c)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    character(*),         intent(out) :: c
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) c
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_i1(kw_input, i)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    integer,              intent(out) :: i
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) i
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r1(kw_input, r)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r2(kw_input, r1, r2)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r1, r2
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r1, r2
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r3(kw_input, r1, r2, r3)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r1, r2, r3
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r1, r2, r3
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r4(kw_input, r1, r2, r3, r4)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r1, r2, r3, r4
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r1, r2, r3, r4
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_error(kw_input)

!   ----------------------------------------------------------------------------
    character(kw_length),   intent(in)  :: kw_input
!   ----------------------------------------------------------------------------
    character(line_length)              :: line
!   ----------------------------------------------------------------------------

    backspace unit_in
    read(unit_in, *) line

    write(*, *) 'error in input line:'
    write(*, *) line
    write(*, *) 'following keyword '//kw_input

    call quit('error in line following keyword '//kw_input)

  end subroutine

  subroutine reset_available_kw_list()

    nr_available_kw = 0
    kw_found        = .false.

  end subroutine

  subroutine check_whether_kw_found(kw_input, kw_section)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: kw_input
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i
!   ----------------------------------------------------------------------------

    if (.not. kw_found) then

      write(*, *) 'illegal keyword '//kw_input//' in section '//kw_section

      write(*, *) 'list of available keywords in section '//kw_section//':'
      do i = 1, nr_available_kw
        write(*, *) available_kw_list(i)
      end do

      call quit('illegal keyword '//kw_input//' in section '//kw_section)

    end if

  end subroutine

  function kw_matches(kw_input, kw_option)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: kw_input
    character(kw_length), intent(in) :: kw_option
    logical                          :: kw_matches
!   ----------------------------------------------------------------------------

    if (lowercase(kw_input) == lowercase(kw_option)) then
      kw_matches = .true.
      kw_found   = .true.
    else
      kw_matches = .false.
      nr_available_kw = nr_available_kw + 1
      available_kw_list(nr_available_kw) = kw_option
    end if

  end function

end module
