
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
!type(KronDelta), intent(in) :: self
type(KronDelta), pointer :: self
call QEnter('print_kd')
if (.not.associated(self)) then
  write(LUPRI,'(/,2x,a)')  &
  'print_kd: entering pointer "self" NOT allocated !'
  call Quit('print_kd: entering pointer "self" NOT allocated !')
else
  write(LUPRI,'(3x,a,2i3)') 'kronecker delta, indexes #:',self%p,self%q
endif

call QExit('print_kd')
end subroutine print_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_kron_delta(i,j) result(this)
!------------------------------------------------------------------
!
!    create pure Kronecker delta - allocate it in memmory 
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer,intent(in) :: i,j
!type(KronDelta) :: this
type(KronDelta), pointer :: this
call QEnter('make_kron_delta')

 nullify(this) ; allocate(this, stat=istat) 

if (i>0 .and. j>0) then
  this%p = i; this%q = j
else 
  write(LUPRI,'(/,2x,a,2i3)')  &
  'make_kron_delta: entering parameters for KD: p,q=',p,q
  write(LUPRI,'(2x,a)')  &
  'entering parameters for kronecker delta must be > 0 !'
  call Quit('entering parameters for kronecker delta must be > 0 !')
endif

if (get_iprint() >= 15) then
  write(LUPRI,'(2x,a,$)') &
  'make_kron_delta: created (allocated) "KronDelta" type'
  call print_kd(this)
endif

call QExit('make_kron_delta')
end function make_kron_delta

End Module KroneckerDelta

