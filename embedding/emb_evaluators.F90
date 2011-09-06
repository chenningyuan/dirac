

module fde_evaluators

   use fde_eval_dirac

   interface fde_get_density_matrix
      module procedure fde_get_dmat_dirac
   end interface

   interface fde_calculate_active_n_gn
      module procedure fde_get_ngn_dirac
   end interface

   interface fde_calculate_active_elpot
      module procedure fde_get_vc_dirac
   end interface

   interface fde_cleanup_density_matrix
      module procedure fde_free_dmat_dirac
   end interface

end module fde_evaluators
