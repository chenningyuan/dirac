
Module KroneckerDelta

  use Control_print

type, public :: KronDelta
 !private
  !integer(kind=4) :: p=0, q=0
  integer :: p=0, q=0
 ! type(KronDelta), pointer :: next => null()
end type KronDelta

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_kd(self)
!------------------------------------------------------------------
!
!    print-out  kronecker delta 
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), intent(in) :: self
call QEnter('print_kd')
!if (get_iprint() >= 5) then
  write(LUPRI,'(3x,a,2i3)') 'kronecker delta, indexes #:',self%p,self%q
!endif
call QExit('print_kd')
end subroutine print_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_kron_delta(i,j) result(this)
!------------------------------------------------------------------
!
!              create pure Kronecker delta
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer,intent(in) :: i,j
type(KronDelta) :: this
call QEnter('make_kron_delta')

 if (i>0 .and. j>0) then
   this%p = i; this%q = j
 else 
  write(LUPRI,'(2x,a,2i3)') 'entering parameters for KD: p,q=',p,q
  write(LUPRI,'(2x,a)') 'entering parameters for kronecker delta must by > 0 !'
  call Quit('entering parameters for kronecker delta must by > 0 !')
 endif

if (get_iprint() >= 15) then
  write(LUPRI,'(2x,a,$)') 'make_kron_delta: created "KronDelta" type'
  call print_kd(this)
endif

call QExit('make_kron_delta')
end function make_kron_delta

End Module KroneckerDelta

