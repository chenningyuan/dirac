module aor_input

  use aor_cfg
  use dirac_interface

  implicit none

  public parse_aor_input

  private

contains

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

  subroutine parse_aor_input(word_io)

!   ----------------------------------------------------------------------------
    character(7), intent(inout) :: word_io
!   ----------------------------------------------------------------------------
    character(7)                :: word
    character(80)               :: line
    character(1)                :: first
    integer                     :: u
!   ----------------------------------------------------------------------------

    if (word_io /= '*AORESP') return

    call aor_cfg_init()
    call ask_for_integrals()

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

          case ('.SKIPEE')
            aor_cfg_skip_pp = .true.

          case ('.SKIPEP')
            aor_cfg_skip_pn = .true.

          case ('.THRESH')
            read(u, *) aor_cfg_threshold_response, aor_cfg_threshold_rhs

!         give frequency for a process
          case ('.FREQ_P')
            read(u, *) aor_cfg_freq_process

!         give 4 frequencies
!         \alpha will use only first 2
!         \beta only first 3
          case ('.FREQ  ')
            read(u, *) aor_cfg_freq(1)
            read(u, *) aor_cfg_freq(2)
            read(u, *) aor_cfg_freq(3)
            read(u, *) aor_cfg_freq(4)

!         frequency-dependent \alpha
          case ('.ALPHA ')
            aor_cfg_alpha = .true.

!         frequency-dependent \beta by the n+1 rule
          case ('.BETA  ')
            aor_cfg_beta = .true.

!         frequency-dependent \beta by the 2n+1 rule
          case ('.BETA2N')
            aor_cfg_beta_2np1 = .true.

!         frequency-dependent \gamma
          case ('.GAMMA ')
            aor_cfg_gamma = .true.

!         static \alpha
          case ('.0ALPHA')
            aor_cfg_alpha_lin_static = .true.

!         static \beta
          case ('.0BETA ')
            aor_cfg_beta_lin_static = .true.

!         static \gamma
          case ('.0GAMMA')
            aor_cfg_gamma_lin_static = .true.

          case ('.DC-SHG')
            aor_cfg_gamma_lin_dc_shg = .true.

          case ('.THG   ')
            aor_cfg_gamma_lin_thg = .true.

          case ('.IDRI  ')
            aor_cfg_gamma_lin_idri = .true.

          case ('.DC-KER')
            aor_cfg_gamma_lin_dc_kerr = .true.

          case ('.INTEGR')
            aor_cfg_integrate = .true.

          case ('.VISUAL')
            aor_cfg_visual = .true.

          case ('.TESTXC')
            aor_cfg_test_xc = .true.

          case ('.TESTGE')
            aor_cfg_test_geo = .true.

          case ('.TESTLO')
            aor_cfg_test_london = .true.

!         calculate static derivatives E_{zzz...}
!         of given order
          case ('.ZZZ   ')
            aor_cfg_zzz = .true.
            read(u, *) aor_cfg_zzz_order

          case default
            call quit('illegal keyword after '//word_io//': '//word)
        end select

        case default
          call quit('illegal keyword after '//word_io//': '//word)
      end select

    end do

  end subroutine

end module
