module fde_mag_cfg


   implicit none

   public fde_rsp_mag_lao
   public fde_rsp_mag_lao_import
   public fde_rsp_mag_lao_export
   public fde_lao_direct_active_embpot
   public fde_lao_direct_active_embker
   public fde_lao_reorth_active_embker
   public fde_lao_frozen_embker
   public fde_rsp_mag_lao_no_s

   logical, save :: fde_rsp_mag_lao               = .false.
   logical, save :: fde_rsp_mag_lao_import        = .false.
   logical, save :: fde_rsp_mag_lao_export        = .false.
   logical, save :: fde_lao_direct_active_embpot  = .false.
   logical, save :: fde_lao_direct_active_embker  = .false.
   logical, save :: fde_lao_reorth_active_embker  = .false.
   logical, save :: fde_lao_frozen_embker         = .false.
   logical, save :: fde_rsp_mag_lao_no_s          = .false.


end module fde_mag_cfg
