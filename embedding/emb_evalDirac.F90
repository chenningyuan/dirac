module fde_eval_dirac

   use memory_allocator
   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
   use fde_integrate
   use vc_grid
   use density_eval
   use ao_eval
   use interface_ao
   use interface_mo_specific
   use xc_max_block_length
   use vc_grid

   implicit none

   public dirac_data_interface

   
   private 

   integer, save :: nba, nba_orb, nba_aux
   integer, save :: dfcoef_unit
   integer, save :: nz = 4

   logical, save :: is_dirac_initialized = .false.

   contains

      
      
!   ----------------------------------------------------------------------------
      subroutine dirac_data_interface 
!   ----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"

         nba_orb = ntbas(0) 
         nba = nba_orb
         dfcoef_unit = LUCOEF

      end subroutine dirac_data_interface 


end module fde_eval_dirac
