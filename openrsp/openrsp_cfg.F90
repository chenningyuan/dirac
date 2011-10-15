module openrsp_cfg

   implicit none

   save

   real(8), public :: openrsp_cfg_threshold_response = 1.0d-8
   real(8), public :: openrsp_cfg_threshold_rhs      = 1.0d-8

   integer, public              :: openrsp_cfg_nr_freq = 0
   real(8), public, allocatable :: openrsp_cfg_freq(:)
   real(8), public              :: openrsp_cfg_freq_process = 0.0d0

   logical, public :: openrsp_cfg_skip_pp = .false.
   logical, public :: openrsp_cfg_skip_pn = .false.

   logical, public :: openrsp_cfg_skip_1el = .false.

   logical, public :: openrsp_cfg_alpha             = .false.
   logical, public :: openrsp_cfg_alpha_lin_static  = .false.
   logical, public :: openrsp_cfg_beta              = .false.
   logical, public :: openrsp_cfg_beta_2np1         = .false.
   logical, public :: openrsp_cfg_beta_lin_static   = .false.
   logical, public :: openrsp_cfg_gamma             = .false.
   logical, public :: openrsp_cfg_gamma_lin_static  = .false.
   logical, public :: openrsp_cfg_gamma_lin_dc_shg  = .false.
   logical, public :: openrsp_cfg_gamma_lin_thg     = .false.
   logical, public :: openrsp_cfg_gamma_lin_idri    = .false.
   logical, public :: openrsp_cfg_gamma_lin_dc_kerr = .false.
   logical, public :: openrsp_cfg_magnetizability   = .false.
   logical, public :: openrsp_cfg_quadrupole        = .false.
   logical, public :: openrsp_cfg_roa               = .false.
   logical, public :: openrsp_cfg_cars              = .false.
   logical, public :: openrsp_cfg_efgb              = .false.
   logical, public :: openrsp_cfg_nolao             = .false.
   logical, public :: openrsp_cfg_jones             = .false.
   logical, public :: openrsp_cfg_cme               = .false.
   logical, public :: openrsp_cfg_dg_energy         = .false.
   logical, public :: openrsp_cfg_test_delta        = .false.
   logical, public :: openrsp_cfg_pv_beta           = .false.
   logical, public :: openrsp_cfg_pv_gamma          = .false.
   logical, public :: openrsp_cfg_pnc               = .false.

   logical, public :: openrsp_cfg_integrate = .false.
   logical, public :: openrsp_cfg_visual    = .false.

   logical, public :: openrsp_cfg_zzz       = .false.
   integer, public :: openrsp_cfg_zzz_order = 0

   integer, parameter :: max_nr_response_functions = 20
   integer, parameter :: max_order = 5

   public response_function
   type response_function
      integer      :: order
      real(8)      :: w(max_order)
      character(1) :: op(max_order)
   end type

   logical, public :: need_1el_f = .false.
   logical, public :: need_1el_q = .false.
   logical, public :: need_1el_v = .false.
   logical, public :: need_1el_o = .false.
   logical, public :: need_1el_b = .false.
   logical, public :: need_1el_g = .false.

   private

end module
