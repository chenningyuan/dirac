module aor_cfg

  implicit none

  save

  real(8) :: aor_cfg_threshold_response
  real(8) :: aor_cfg_threshold_rhs
  real(8) :: aor_cfg_freq(4)
  real(8) :: aor_cfg_freq_process

  logical :: aor_cfg_skip_pp
  logical :: aor_cfg_skip_pn

  logical :: aor_cfg_alpha
  logical :: aor_cfg_alpha_lin_static
  logical :: aor_cfg_beta
  logical :: aor_cfg_beta_2np1
  logical :: aor_cfg_beta_lin_static
  logical :: aor_cfg_gamma
  logical :: aor_cfg_gamma_lin_static
  logical :: aor_cfg_gamma_lin_dc_shg
  logical :: aor_cfg_gamma_lin_thg
  logical :: aor_cfg_gamma_lin_idri
  logical :: aor_cfg_gamma_lin_dc_kerr

  logical :: aor_cfg_integrate
  logical :: aor_cfg_visual
  logical :: aor_cfg_test_xc
  logical :: aor_cfg_test_geo
  logical :: aor_cfg_test_london

  logical :: aor_cfg_zzz
  integer :: aor_cfg_zzz_order

contains

  subroutine aor_cfg_init()

    aor_cfg_threshold_response = 1.0d-8
    aor_cfg_threshold_rhs      = 1.0d-8
    aor_cfg_freq(4)            = 0.0d0
    aor_cfg_freq_process       = 0.0d0

    aor_cfg_skip_pp            = .false.
    aor_cfg_skip_pn            = .false.

    aor_cfg_alpha              = .false.
    aor_cfg_alpha_lin_static   = .false.
    aor_cfg_beta               = .false.
    aor_cfg_beta_2np1          = .false.
    aor_cfg_beta_lin_static    = .false.
    aor_cfg_gamma              = .false.
    aor_cfg_gamma_lin_static   = .false.
    aor_cfg_gamma_lin_dc_shg   = .false.
    aor_cfg_gamma_lin_thg      = .false.
    aor_cfg_gamma_lin_idri     = .false.
    aor_cfg_gamma_lin_dc_kerr  = .false.

    aor_cfg_integrate          = .false.
    aor_cfg_visual             = .false.
    aor_cfg_test_xc            = .false.
    aor_cfg_test_geo           = .false.
    aor_cfg_test_london        = .false.

    aor_cfg_zzz                = .false.
    aor_cfg_zzz_order          = 0

  end subroutine

end module
