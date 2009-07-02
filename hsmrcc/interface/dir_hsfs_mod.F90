!
module DIRAC_HSFS
!
! HSFS parameters read in through the DIRAC input
!

! convergence
integer :: iconv_hsfs=8

! print level
integer :: iprint_hsfs=0

integer, parameter :: luhsfs=18

!! HSFS indexes 
!integer, allocatable :: iord_hsfs(:) 

! ... depends only on integral file (no transformation etc..)
logical :: i_file_based=.false.

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine open_dirac_hsfs_file
!----------------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"

open(luhsfs,file="DIRAC_hsfs",form="unformatted",status="unknown")

write(LUPRI,*) 'DIRAC file was open.'

end subroutine open_dirac_hsfs_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine close_dirac_hsfs_file
!----------------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"

close(luhsfs,status="keep")

write(LUPRI,*) 'DIRAC file was closed.'

end subroutine close_dirac_hsfs_file

end module DIRAC_HSFS
