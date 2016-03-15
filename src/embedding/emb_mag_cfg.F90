module fde_mag_cfg


   implicit none

   public fde_rsp_mag_lao
   public fde_rsp_mag_lao_no_s
   public fde_rsp_mag_lao_import
   public fde_rsp_mag_lao_export
   public fde_lao_direct_active_embpot
   public fde_lao_direct_active_embker
   public fde_lao_reorth_active_embker
   public fde_lao_frozen_embker
   public fde_lao_frozen_embker_coulomb
   public fde_magn_expval_lao
   public fde_lao_magn_expval_no_embpot
   public fde_lao_magn_expval_no_uncoup_embker
   public fde_lao_magn_expval_no_coupl_nonadd_embker


   logical, save :: fde_rsp_mag_lao_no_s          = .false.

!  import/export of perturbed density:
   logical, save :: fde_rsp_mag_lao_import        = .false.
   logical, save :: fde_rsp_mag_lao_export        = .false.

!  these indicate contributions to property gradient (in LAO basis):
   logical, save :: fde_rsp_mag_lao               = .false.
   logical, save :: fde_lao_direct_active_embpot  = .false.
   logical, save :: fde_lao_direct_active_embker  = .false.
   logical, save :: fde_lao_reorth_active_embker  = .false.
   logical, save :: fde_lao_frozen_embker         = .false.
   logical, save :: fde_lao_frozen_embker_coulomb = .false.

!  these indicate contributions to expectation value to magnetizability (in LAO basis):
   logical, save :: fde_magn_expval_lao = .false.
   logical, save :: fde_lao_magn_expval_no_embpot = .false.
   logical, save :: fde_lao_magn_expval_no_uncoup_embker = .false.
   logical, save :: fde_lao_magn_expval_no_coupl_nonadd_embker = .false.


end module fde_mag_cfg
