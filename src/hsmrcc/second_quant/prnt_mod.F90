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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module control_print
!-----------------------------------------------------------
!
!  Contains printing varible for control printouts
! throughout the code.
!
!-----------------------------------------------------------

!use mystring

integer, private :: iprint = 0

contains

subroutine set_iprint(i)
!-----------------------------------------------------------
!
!-----------------------------------------------------------
integer, intent(in) :: i
call QEnter('set_iprint')
iprint = i
call QExit('set_iprint')
end subroutine set_iprint

function get_iprint() result(i)
call QEnter('get_iprint')
   i = iprint
call QExit('get_iprint')
end function get_iprint

end module control_print
