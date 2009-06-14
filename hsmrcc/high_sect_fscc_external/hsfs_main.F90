
#if defined STANDALONE_HSFS
 Program HIGH_SECTOR_FSCC_DRIVER
#else
 Subroutine HIGH_SECTOR_FSCC_DRIVER
#endif

  use high_sect_execs

#if defined STANDALONE_HSFS
  call get_line_arguments
#endif

#if defined HSFS_DYN_ALLOC
  call Init_HSFS_Modules 
#endif

  call launch_hsfs_methods

 End

 Subroutine Init_HSFS_Modules
!
! iniialize the whole bunch of former static allocations in common blocks
!
use gen_mod
use cint_mod
!use dat_mod
use fock_mod
!use cij1_mod
!use cij2_mod
!use cij3_mod
!use cij4_mod
!use cij5_mod
!use cij6_mod
 
  call init_gen_mod
  call init_cint_mod
! call init_dat_mod
  call init_fock_mod
! call init_cij1_mod
! call init_cij2_mod
! call init_cij3_mod
! call init_cij4_mod
! call init_cij5_mod
! call init_cij6_mod

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

