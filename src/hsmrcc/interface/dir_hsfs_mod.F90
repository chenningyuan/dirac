!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
