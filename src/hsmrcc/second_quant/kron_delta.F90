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

nullify(this) ; allocate(this, stat=ierr) 
if (ierr /= 0) then
  write(LUPRI,'(2x,a)')  &
  'make_kron_delta: memory allocation error; ierr=',ierr
  call Quit('make_kron_delta: memory allocation error !')
endif

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_kd(this)
!------------------------------------------------------------------
!  deallocated Krond.delta
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: this
call QEnter('deallocate_kd')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)')  &
  'deallocate_kd: entering pointer "this" NOT allocated !'
  call Quit('deallocate_kd: entering pointer "this" NOT allocated !')
else

 if ( get_iprint() >= 8 )  then
   write(LUPRI,'(/,2x,a)')   & 
   '*** deallocate_kd: going to deallocate this Kron.delta:'
   call print_kd(this)
 endif

 deallocate(this, stat=ierr)
 if (ierr /= 0) then
    write(LUPRI,'(2x,a)')  &
   'deallocate_kd: memory deallocation error; ierr=',ierr
    call Quit('deallocate_kd: memory deallocation error !')
 endif
 nullify(this)
endif

call QExit('deallocate_kd')
end subroutine deallocate_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_kd(source, dest)
!------------------------------------------------------------------
!
!            copy Kronecker deltas:  source -> dest
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: source, dest
call QEnter('copy_kd')

if (.not.associated(source).or..not.associated(dest)) then
  call Quit('copy_kd: "source" or "dest" NOT allocated !')
endif

dest%p = source%p
dest%q = source%q

if ( get_iprint() >=8 ) then
 write(LUPRI,'(/,2x,a)') '*** copy_kd: source'
 call print_kd(source)
 write(LUPRI,'(/,2x,a)') '*** copy_kd: destination'
 call print_kd(dest)
endif

call QExit('copy_kd')
end subroutine copy_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_KD(kd1,kd2) result(equal)
!------------------------------------------------------------------
!
!          Compare two Kronecker deltas, "kd1", "kd2"
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: kd1,kd2
logical :: equal
call QEnter('equal_KD')
equal = .false.
if (associated(kd1).and.associated(kd2)) then
   if ( (kd1%p==kd2%p).and.(kd1%q==kd2%q)) equal = .true.
endif
call QExit('equal_KD')
end function equal_KD

End Module KroneckerDelta

