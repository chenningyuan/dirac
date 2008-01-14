Module SecQuantOperator

type, public :: SQOper
 !private
   !integer(kind=4) :: p = 0
   integer :: p = 0
   logical :: is_creation_oper = .false.
   !type(SQOper), pointer :: next => null()
end type SQOper
 
 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_sqo(i,is_creat) result (this)
!----------------------------------------------------------------------------
!
!              Returns second quantized operator
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!integer(kind=4), intent(in) :: i
integer, intent(in) :: i
logical, optional, intent(in) :: is_creat
type(SQOper) :: this

call QEnter('make_sqo')

! write(*,'(2x,a,i3,a,l)')   &
!  'make_sqo: entering i=',i,' entering is_creat=',is_creat
  
 if (i>0) then
  this%p = i
  if ( present(is_creat) ) then
    this%is_creation_oper = is_creat
  else
    ! ... by default ...
    this%is_creation_oper = .FALSE.
  endif
 else
   call Quit('make_sqo: bad entering parameter i !')
 endif
call QExit('make_sqo')
end function make_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_sqo(this)
!----------------------------------------------------------------------------
!
!    print attributes of the entering second quantization operator "this"
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), intent(in) :: this
call QEnter('print_sqo')
write(LUPRI,'(2x,a,i3,$)') ' sec.quant.oper. - index no',this%p
if (this%is_creation_oper) then
  write(LUPRI,*) ',creation oper.'
else 
  write(LUPRI,*) ',anihilation oper.'
endif
call QExit('print_sqo')
end subroutine print_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function is_sqo_filled(this) result(is)
!----------------------------------------------------------------------------
!
!    returns info if the second quantization operator is properly filled
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), intent(in) :: this
logical :: is 
call QEnter('is_sqo_filled')

is = .false.
if (this%p > 0) is = .true.

call QExit('is_sqo_filled')
end function is_sqo_filled

End Module SecQuantOperator

