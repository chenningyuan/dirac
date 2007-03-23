!
! Module contains basic variables of the DIRAC memmory array
!
MODULE Memmory

! it's dangerously public !
integer :: KWORK, KFREE, LFREE, LWORK ! basic variables of memmory handling
!integer, private :: KWORK, KFREE, LFREE, LWORK ! basic variables of memmory handling
double precision,  pointer :: PW(:) ! points to the big WORK array

  contains

! print out basic memmory variables
subroutine print_mem_param 
#include "priunit.h"
CALL QENTER('print_mem_param')
write(LUPRI,*) 'Memmory parameters: KWORK=',KWORK,' KFREE=',KFREE,' LFREE=',LFREE, &
 ' LWORK=',LWORK 
if (associated(PW)) then
   write(LUPRI,*) 'Pointer PW is associated with the WORK array.'
else
   write(LUPRI,*) 'Pointer PW is NOT associated with the WORK array.'
endif
CALL QEXIT('print_mem_param')
end subroutine print_mem_param

SUBROUTINE Stick_into_array(WORK,LLWORK)
#include "priunit.h"
! integer, intent(in) :: LLWORK
integer*4, intent(in) :: LLWORK
double precision, intent(in),  target :: WORK(LLWORK)
CALL QENTER('Stick_into_array')
if (associated(PW)) then
  write(LUPRI,*) 'Stick_into_array: PW is already associated - Exit !'
  CALL QUIT('Stick_into_array: PW is already associated !')
else
  write(LUPRI,*) 'Stick_into_array: PW is not associated'
endif
PW => WORK ! PW points now to the big WORK array ... 
LWORK = LLWORK ! actualize module variable
write(LUPRI,*) 'LLWORK=',LLWORK
if (size(PW) /= LLWORK) then
 write(LUPRI,*) 'Stick_into_array:  size(PW)=',size(PW),' while LWORK=',LWORK
 write(LUPRI,*) 'somme error in array sticking !'
 CALL QUIT('Stick_into_array: error in length !')
endif
#include "memint.h"
CALL MEMCHK('Stick_into_array',PW,1)
write(LUPRI,*) 'Stick_into_array: PW sticked into WORK, memmory checked !'
CALL QEXIT('Stick_into_array')
END SUBROUTINE Stick_into_array

SUBROUTINE UnStickArray
#include "priunit.h"
if (associated(PW)) then
  write(LUPRI,*) 'PW is already associated'
else
  write(LUPRI,*) 'PW is NOT associated !'
endif
nullify(PW)
END SUBROUTINE UnStickArray

END MODULE Memmory
