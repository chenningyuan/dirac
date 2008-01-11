
Module KroneckerDelta

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
!  print-out  kronecker delta 
!
!------------------------------------------------------------------
type(KronDelta), intent(in) :: self
write(*,'(3x,a,2i3)') 'kronecker delta, indexes #:',self%p,self%q
end subroutine print_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_kron_delta(i,j) result(this)
!------------------------------------------------------------------
!
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!integer(kind=4),intent(in) :: i,j
integer,intent(in) :: i,j
type(KronDelta) :: this

 if (i>0 .and. j>0) then
   this%p = i; this%q = j
 else 
  write(LUPRI,*) 'entering parameters for KD: p,q=',p,q
  write(LUPRI,*) 'entering parameters for kronecker delta must by > 0 !'
  call Quit('entering parameters for kronecker delta must by > 0 !')
 endif
end function make_kron_delta

End Module KroneckerDelta

