!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 Subroutine HIGH_SECTOR_FSCC_METHODS

  use high_sect_execs

#if defined HSFS_DYN_ALLOC
  call Init_HSFS_Modules ! initialize all modules
#endif

  call launch_hsfs_methods

 End Subroutine HIGH_SECTOR_FSCC_METHODS

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 Subroutine Init_HSFS_Modules
!---------------------------------------------------------------------------------------
!
!   initialize the whole bunch of former static allocations in former common blocks
!  which are now replaced by nice F90 modules with allocatable variables
!
!---------------------------------------------------------------------------------------
use hsfs_parameters

use gen_mod
use sym_mod
use cint_mod
use det_mod
use dat_mod
use fock_mod
use cij1_mod
use cij2_mod
use cij3_mod
use cij4_mod
use cij5_mod
use cij6_mod
use cfiles_mod
use rsp1_mod
use rsp2_mod
use dpd_mod
use ctmp_mod
use crle_mod
use prm_mod
use idv_mod
use ijtemp_mod

 call header ('Output from Init_HSFS_Modules',-1)

 ! very important - assigning dimensions of arrays  for various allocations
 call assign_hsfs_parameters
 
 call init_gen_mod
 call init_sym_mod
 call init_cint_mod
 call init_det_mod
 call init_dat_mod
 call init_fock_mod
 call init_cij1_mod
 call init_cij2_mod
 call init_cij3_mod
 call init_cij4_mod
 call init_cij5_mod
 call init_cij6_mod
 call init_cfiles_mod
 call init_rsp1_mod
 call init_rsp2_mod
 call init_dpd_mod
 call init_ctmp_mod
 call init_crle_mod
 call init_prm_mod
 call init_idv_mod
 call init_ijtemp_mod

End Subroutine Init_HSFS_Modules

Subroutine alloc_stat_check(alloc_label,istat)
#include "priunit.h"
  character(*), intent(in) :: alloc_label
  integer, intent(in) :: istat
  if (istat /= 0) then
    write(LUPRI,*) 'Error in allocation, variable :',alloc_label
    write(LUPRI,*) 'istat=',istat
    call quit('error in F90 dynamic allocation !')
  endif
End Subroutine alloc_stat_check
