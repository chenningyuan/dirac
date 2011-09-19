

module fde_evaluators

   use fde_eval_dirac

   public fde_get_density_matrix
   public fde_calculate_active_n_gn
   public fde_calculate_active_elpot
   public fde_cleanup_density_matrix
   public fde_calculate_embpot 

!   interface fde_get_density_matrix
!      module procedure fde_get_dmat_dirac
!   end interface

!   interface fde_calculate_active_n_gn
!      module procedure fde_get_ngn_dirac
!   end interface

!   interface fde_calculate_active_elpot
!      module procedure fde_get_vc_dirac
!   end interface

!   interface fde_cleanup_density_matrix
!      module procedure fde_free_dmat_dirac
!   end interface

   interface fde_calculate_embpot
      module procedure fde_calculate_vemb_dirac
   end interface

!  interface fde_calculate_embrsp
!     module procedure fde_calculate_lr_dirac
!  end interface

!  interface fde_calculate_embqrp
!     module procedure fde_calculate_qr_dirac
!  end interface

end module fde_evaluators
