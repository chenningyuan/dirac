module aor_input

  use aor_cfg
  use dirac_interface

  implicit none

  public parse_aor_input

  private

contains

  subroutine parse_aor_input(word_io)

!   ----------------------------------------------------------------------------
    character(7), intent(inout) :: word_io
!   ----------------------------------------------------------------------------
    character(7)                :: word
    character(80)               :: line
    character(1)                :: first
!   ----------------------------------------------------------------------------

    if (word_io /= '*AORESP') return

    do while (.true.)

      read (unit_input, '(a7)') word
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

          case ('.THRESH')
            read(unit_input, *) aor_cfg_threshold_response, aor_cfg_threshold_rhs

!         give 4 frequencies
!         \alpha will use only first 2
!         \beta only first 3
          case ('.FREQ  ')
            read(unit_input, *) aor_cfg_freq(1)
            read(unit_input, *) aor_cfg_freq(2)
            read(unit_input, *) aor_cfg_freq(3)
            read(unit_input, *) aor_cfg_freq(4)

!         give frequency for a process
          case ('.FREQ_P')
            read(unit_input, *) aor_cfg_freq_process

          case ('.SKIPPP', '.SKIPEE')
            aor_cfg_skip_pp = .true.

          case ('.SKIPPN', '.SKIPEP')
            aor_cfg_skip_pn = .true.

!         frequency-dependent \alpha
          case ('.ALPHA ')
            aor_cfg_alpha = .true.

!         static \alpha
          case ('.0ALPHA')
            aor_cfg_alpha_lin_static = .true.

!         frequency-dependent \beta by the n+1 rule
          case ('.BETA  ')
            aor_cfg_beta = .true.

!         frequency-dependent \beta by the 2n+1 rule
          case ('.BETA2N')
            aor_cfg_beta_2np1 = .true.

!         static \beta
          case ('.0BETA ')
            aor_cfg_beta_lin_static = .true.

!         frequency-dependent \gamma
          case ('.GAMMA ')
            aor_cfg_gamma = .true.

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

          case ('.MAGNET')
            aor_cfg_magnetizability = .true.

          case ('.QUADRU')
            aor_cfg_quadrupole = .true.

          case ('.ROA   ')
            aor_cfg_roa = .true.

          case ('.CARS  ')
            aor_cfg_cars = .true.

          case ('.MOLGRD')
            aor_cfg_dg_energy = .true.

          case ('.INTEGR')
            aor_cfg_integrate = .true.

          case ('.VISUAL')
            aor_cfg_visual = .true.

!         calculate static derivatives E_{zzz...}
!         of given order
          case ('.ZZZ   ')
            aor_cfg_zzz = .true.
            read(unit_input, *) aor_cfg_zzz_order

          case default
            call quit('illegal keyword after '//word_io//': '//word)
        end select

        case default
          call quit('illegal keyword after '//word_io//': '//word)
      end select

    end do

  end subroutine

end module
