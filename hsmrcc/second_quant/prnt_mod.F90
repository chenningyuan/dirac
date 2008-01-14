!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module control_print
!-----------------------------------------------------------
!
!  Contains printing varible for control printouts
! throughout the code.
!
!-----------------------------------------------------------

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
