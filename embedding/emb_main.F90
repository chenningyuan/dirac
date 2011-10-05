!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module fde_mod

   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
   use fde_exportPotentials
!  use fde_import
   use FDE_ImportPotentials
   use fde_integrate

   use fde_inputDirac
   use fde_eval_dirac

   interface fde_init
      module procedure fdeini
   end interface

   interface fde_input
      module procedure fdeinp
   end interface


#ifdef MOD_UNRELEASED
!  interface fde_cleanup_density_matrix
!     module procedure fde_free_dmat_dirac
!  end interface

   interface fde_calculate_embpot
      module procedure fde_calculate_vemb_dirac
   end interface

!  interface fde_calculate_embrsp
!     module procedure fde_calculate_lr_dirac
!  end interface

!  interface fde_calculate_embqrp
!     module procedure fde_calculate_qr_dirac
!  end interface

#endif

end module fde_mod
