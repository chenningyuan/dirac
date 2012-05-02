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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Module SecQuantOperator

use control_print

type, public :: SQOper
 !private
  ! ... index 
   integer :: p = 0
   logical :: is_creation_oper = .false.
!  ... operator type; 1-general, 2-occupied, 3-virtual
   integer :: op_type = 1 ! default - general
end type SQOper

character*8, dimension(3) :: op_label = (/"general ","occupied","virtual "/) 

 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_sqo(i,is_creat,is_occup) result (this)
!----------------------------------------------------------------------------
!
!      Returns pointer "this" to ALLOCATED second quantized operator
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i
logical, optional, intent(in) :: is_creat,is_occup
type(SQOper), pointer :: this 
call QEnter('make_sqo')

if (get_iprint() > 70) then
!if (get_iprint() > 7) then
 write(LUPRI,'(2x,a,i3,a,l1)')   &
 'make_sqo: entering i=',i,' entering is_creat=',is_creat

 if ( present(is_occup) ) then
   write(LUPRI,'(2x,a,i3,a,l1)')   &
   'make_sqo: also entering is_occup=',is_occup
 endif

endif

nullify(this); allocate(this, stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a,i3)')  &
  'make_sqo: memory allocation error; aborting; ierr=',ierr
  call Quit('make_sqo: memory allocation error !')
endif
  
if (i>0) then
  this%p = i
  if ( present(is_creat) ) then
    this%is_creation_oper = is_creat
  else
    ! ... by default ...
    this%is_creation_oper = .FALSE.
  endif

  if ( present(is_occup) ) then
    if (is_occup) then
      this%op_type = 2 ! ... occupied (i,j,k,...)
    else
      this%op_type = 3 ! ... virtual (a,b,c,...) / otherwise general (=1)
    endif
  else
    ! ... by default ...
    this%op_type = 1 ! ... otherwise general
  endif
else
   call Quit('make_sqo: bad entering parameter i !')
endif

call QExit('make_sqo')
end function make_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_sqo(source, dest)
!----------------------------------------------------------------------------
!                copy SQoperators: source -> dest
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), pointer :: source, dest
call QEnter('copy_sqo')

if (.not.associated(source).or..not.associated(dest)) then
  call Quit('copy_sqo: "source" or "dest" NOT associated !')
endif

dest%p = source%p
dest%is_creation_oper = source%is_creation_oper
dest%op_type = source%op_type

if (get_iprint() >= 9) then
   write(LUPRI,'(/,2x,a)') 'copy_sqo: source'
   call print_sqo(source)
   write(LUPRI,'(/,2x,a)') 'copy_sqo: copied'
   call print_sqo(dest)
endif

call QExit('copy_sqo')
end subroutine copy_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_sqo(this)
!----------------------------------------------------------------------------
!
!            DEALLOCATE entering sec.quantized operator, this
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), pointer :: this
call QEnter('deallocate_sqo')

if (get_iprint() >=8 ) then
  write(LUPRI,'(/,2x,a)')  &
  '*** deallocate_sqo: going to deallocate this sqo : '
  call print_sqo(this)
endif

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)')  &
  'deallocate_sqo: entering pointer NOT associated ! Abort !'
  call Quit('deallocate_sqo: entering pointer NOT associated !')
else
  deallocate(this, stat = ierr)
  if (ierr /= 0) then
   write(LUPRI,'(/,2x,a)')  &
   'deallocate_sqo: deallocation error !'
   call Quit('deallocate_sqo: deallocation error !')
  endif
  nullify(this)
endif

call QExit('deallocate_sqo')
end subroutine deallocate_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_sqo(this)
!----------------------------------------------------------------------------
!
!    print attributes of the entering second quantization operator "this"
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!type(SQOper), intent(in) :: this
type(SQOper), pointer :: this
call QEnter('print_sqo')

if (.not.associated(this)) then
  write(LUPRI,'(2x,a)') 'print_sqo: operator NOT allocated !' 
  call Quit('print_sqo: operator NOT allocated !')
else
  write(LUPRI,'(2x,a,i3,$)') ' sec.quant.oper. - index no',this%p
  if (this%is_creation_oper) then
   ! write(LUPRI,*) ',creation oper.; ',op_label(this%op_type)
    write(LUPRI,*) op_label(this%op_type),' creation operator'
  else 
   ! write(LUPRI,*) ',anihilation oper.; ',op_label(this%op_type)
    write(LUPRI,*) op_label(this%op_type),' anihilation operator'
  endif
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
type(SQOper), pointer :: this
logical :: is 
call QEnter('is_sqo_filled')

if (.not.associated(this)) then
  write(LUPRI,'(2x,a)') 'is_sqo_filled: operator NOT allocated ! Abort !'
  call Quit('is_sqo_filled: operator NOT allocated ! Abort !')
else
  is = .false.
  if (this%p > 0) is = .true.
endif

call QExit('is_sqo_filled')
end function is_sqo_filled

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_SQ_operators(oper1,oper2) result (equal)
!----------------------------------------------------------------------------
!
!         Compare two SQ operators - their all three attributes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), pointer :: oper1,oper2
logical :: equal
call QEnter('equal_SQ_operators')

equal = .false.

if ( associated(oper1) .and. associated(oper2) ) then

  if  (  (oper1%p == oper2%p)   .and.   &
      ( (oper1%is_creation_oper .and. oper2%is_creation_oper) .or.    &
        (.not.oper1%is_creation_oper .and. .not.oper2%is_creation_oper) ) .and.  &
       (oper1%op_type == oper2%op_type) ) then

       equal = .true.

 endif

else
  call Quit('equal_SQ_operators: entering operators NOT associated !')
endif

call QExit('equal_SQ_operators')
end function equal_SQ_operators

End Module SecQuantOperator

