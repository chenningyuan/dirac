module aor_cfg

  implicit none

  save

  real(8), public :: aor_cfg_threshold_response = 1.0d-8
  real(8), public :: aor_cfg_threshold_rhs      = 1.0d-8
  real(8), public :: aor_cfg_freq(4)            = 0.0d0
  real(8), public :: aor_cfg_freq_process       = 0.0d0

  logical, public :: aor_cfg_skip_pp            = .false.
  logical, public :: aor_cfg_skip_pn            = .false.

  logical, public :: aor_cfg_alpha              = .false.
  logical, public :: aor_cfg_alpha_lin_static   = .false.
  logical, public :: aor_cfg_beta               = .false.
  logical, public :: aor_cfg_beta_2np1          = .false.
  logical, public :: aor_cfg_beta_lin_static    = .false.
  logical, public :: aor_cfg_gamma              = .false.
  logical, public :: aor_cfg_gamma_lin_static   = .false.
  logical, public :: aor_cfg_gamma_lin_dc_shg   = .false.
  logical, public :: aor_cfg_gamma_lin_thg      = .false.
  logical, public :: aor_cfg_gamma_lin_idri     = .false.
  logical, public :: aor_cfg_gamma_lin_dc_kerr  = .false.
  logical, public :: aor_cfg_magnetizability    = .false.
  logical, public :: aor_cfg_quadrupole         = .false.
  logical, public :: aor_cfg_roa                = .false.
  logical, public :: aor_cfg_cars               = .false.
  logical, public :: aor_cfg_dg_energy          = .false.

  logical, public :: aor_cfg_integrate          = .false.
  logical, public :: aor_cfg_visual             = .false.

  logical, public :: aor_cfg_zzz                = .false.
  integer, public :: aor_cfg_zzz_order          = 0

  private

end module
