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
!  Contains printing varible for control printouts
! throughout the code.
!-----------------------------------------------------------
!use mystring
integer, private :: iprint = 0
contains
subroutine set_iprint(i)
integer, intent(in) :: i
iprint = i
end subroutine set_iprint
function get_iprint() result(i)
   i = iprint
end function get_iprint
end module control_print

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
!    prints out  kronecker delta 
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: self
if (.not.associated(self)) then
  write(LUPRI,'(/,2x,a)')  &
  'print_kd: entering pointer "self" NOT allocated !'
  call Quit('print_kd: entering pointer "self" NOT allocated !')
else
  write(LUPRI,'(3x,a,2i3)') 'kronecker delta, indexes #:',self%p,self%q
endif
end subroutine print_kd
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_kron_delta(i,j) result(this)
!------------------------------------------------------------------
!  creates pure Kronecker delta by allocating it in the memmory 
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer,intent(in) :: i,j
type(KronDelta), pointer :: this
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
end function make_kron_delta

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_kd(this)
!------------------------------------------------------------------
!  deallocated Krond.delta
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: this
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
end subroutine deallocate_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_kd(source, dest)
!------------------------------------------------------------------
!      copy Kronecker deltas:  source -> dest
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: source, dest
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
end subroutine copy_kd

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_KD(kd1,kd2) result(equal)
!------------------------------------------------------------------
!          Compare two Kronecker deltas, "kd1", "kd2"
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(KronDelta), pointer :: kd1,kd2
logical :: equal
equal = .false.
if (associated(kd1).and.associated(kd2)) then
   if ( (kd1%p==kd2%p).and.(kd1%q==kd2%q)) equal = .true.
endif
end function equal_KD
End Module KroneckerDelta

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
!      Returns pointer "this" to ALLOCATED second quantized operator
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i
logical, optional, intent(in) :: is_creat,is_occup
type(SQOper), pointer :: this 
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
end function make_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_sqo(source, dest)
!----------------------------------------------------------------------------
!                copy SQoperators: source -> dest
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), pointer :: source, dest
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
end subroutine copy_sqo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_sqo(this)
!----------------------------------------------------------------------------
!          DEALLOCATE entering sec.quantized operator, "this"
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(SQOper), pointer :: this

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

Module Indexed_Quantities
!#if !defined (VAR_IFORT)
!----------------------------------------------------------------------------
!
!       Module for handling various indexed tems (quantities)
!
!----------------------------------------------------------------------------

use mystring
use control_print

type indexed_term
! assigned unique ID, for instance 1 - Fij, 2 - Tij, 3 - Vpqrs, 4 - Tijab ... 
integer :: id = -1
! ... store the selected name of undexed quantity (like T_ij^ab, F_pq,...)
!character(len=20) :: it_name
type(string), dimension(1) :: it_name
! ... number of indexes (2,4,6...) - must be even number
integer :: nindx 
! ... array for indexes of indexes 
integer, allocatable, dimension(:) :: indxs

end type indexed_term

! ... common name
interface fill_indexes
  module  procedure fill_indexes_2
  module  procedure fill_indexes_4
!  module  procedure fill_indexes_6
!  module  procedure fill_indexes_8
  module  procedure fill_indexes_array
end interface

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_it(id_in, nindx_in, it_name_in ) result(this)
!----------------------------------------------------------------------------
!
!                 Returns ALLOCATED indexed quantity
!
! Entering: id_in - unique ID
!           nindx_in - number of indexes (even number)
!           it_name_in - string to be attached
!
!                  Like  T_ij, Zijab, Xpq, Y_abi ....
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: id_in, nindx_in
character(len=*),optional, intent(in) :: it_name_in
type(indexed_term), pointer :: this

call QEnter('make_it')

if (get_iprint() > 70) then
 write(LUPRI,'(2x,a,i3,a,i3)')   &
  'make_it: entering id=',id_in,' entering nindx=',nindx_in
   if ( present(it_name_in) ) then
     write(LUPRI,'(2x,a,a)')   &
     'make_it: also entering name=',it_name_in
   endif
endif

nullify(this)
allocate(this, stat=ierr)

if (ierr /= 0) then
  write(LUPRI,'(2x,a)') 'make_it: error in "indexed_term" cell allocation !'
  call Quit('make_it: error in in "indexed_term" cell allocation !')
endif

! ... check even number of indexes
if ( MOD(nindx_in,2) /= 0) then
  write(LUPRI,'(2x,a,i3)')  &
  'number of indexes NOT even number !, nindx_in=',nindx_in
  call Quit('make_it: number of indexes NOT even number !')
endif
this%nindx = nindx_in

this%id = id_in

allocate(this%indxs(this%nindx),stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(2x,a)') 'error in allocation of "this%indxs(nindx)" !'
  call Quit('make_it: error in allocation of  "this%indxs(nindx)" !')
endif

this%indxs = -1 

! ... if present string on parameter list, treat it 
if ( present(it_name_in) ) then
   call cancel_string(this%it_name)
   call make_string(this%it_name,it_name_in)
endif

call QExit('make_it')
end function  make_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_it(source, dest)
!----------------------------------------------------------------------------
!
! Copy indexed term : source -> dest, where "source" is completely allocated,
! and "dest" needs allocate some of its attributes.
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: source, dest
call QEnter('copy_it')

if (.not.associated(source).or..not.associated(dest)) then
 call Quit('copy_it: "source" or "dest" NOT associated !')
endif

if (get_iprint() >=8) then
  write(LUPRI,'(/,4x,a)') 'copy_it: going to copy from "source":'
  call print_it(source)
endif

dest%id = source%id
dest%nindx = source%nindx

if (.not.allocated(source%indxs)) then
  call Quit('copy_it: indexes of "source" it NOT allocated ! Wrong !')
else
  allocate(dest%indxs(dest%nindx),stat=ierr)
  if (ierr /= 0) then
    write(LUPRI,'(/,2x,a)') 'copy_it: error in allocation of "this%indxs(nindx)" !'
    call Quit('copy_it: error in allocation of  "this%indxs(nindx)" !')
  endif
  ! ... copy indexes 
  do i = 1, dest%nindx
     dest%indxs(i) = source%indxs(i)
  enddo
endif

! ... now take care of the label if exists
if ( is_string(source%it_name) )  then
 !call cancel_string(dest%it_name)
! ... perhaps it will work ?
 !call make_string(dest%it_name,source%it_name)    
 !call make_string(dest%it_name, get_string(source%it_name) ) 
 call copy_strings(source%it_name, dest%it_name) 
endif

if (get_iprint() >=8) then
  write(LUPRI,'(/,4x,a)') 'copy_it: copied from "source" to "dest":'
  call print_it(dest)
endif

call QExit('copy_it')
end subroutine copy_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine cancel_it(this)
!----------------------------------------------------------------------------
!
!   completely deallocates entering indexed term pointed with "this"
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
call QEnter('cancel_it')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)') 'cancel_it: entering pointer NOT associated !  '
  call Quit('cancel_it: entering pointer NOT associated !')
endif

if (get_iprint() >= 8) then
  write(LUPRI,'(/,2x,a)')   &
  '*** cancel_it: going to deallocate this indexed term :'
  call print_it(this)
endif

if (is_string(this%it_name)) then
 ! ... deallocate the attached string first
 call cancel_string(this%it_name)
endif

if (allocated(this%indxs)) then
  ! ... deallocates indexes
  deallocate(this%indxs,stat=ierr)
  if (ierr /= 0) then
    call Quit('cancel_it: error in deallocation of indexes !')
  endif
else
 ! ... required indexes NOT allocated; return hard exit 
 call print_it(this)
 write(LUPRI,'(2x,a)') 'cancel_it: NO allocated indexes !' & 
                       //'These are required ! Quit!'
 call Quit('cancel_it: NO allocated indexes ! Quit!')
endif

! ... now deallocate the whole cell associated with the entering pointer
deallocate(this,stat=ierr)
if (ierr /= 0) then
 call Quit('cancel_it: error in deallocation of "this" !')
endif
nullify(this)

call QExit('cancel_it')
end subroutine cancel_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes_array(this,indxs_array)
!----------------------------------------------------------------------------
!
!      fills indexed term indexes numbers from array
!
!  entering "indxs_array" contains indexes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indxs_array(:)
call QEnter('fill_indexes_array')

if (.not.associated(this)) then
 call Quit('fill_indexes: pointer NOT associated !')
endif

if (size(indxs_array) /= this%nindx) then
 call Quit('fill_indexes: no indexes !')
endif

do i = 1, this%nindx
   this%indxs(i) = indxs_array(i)
   if (indxs_array(i)<=0) then
     write(LUPRI,'(2x,a,i6)')   &
     'fill_indexes:   something wrong ...indxs_array(i)=', &
     indxs_array(i)
     call Quit('fill_indexes: index <= 0 !')
   endif
enddo

if (get_iprint() > 70) then
  write(LUPRI,'(2x,a,i3)') 'fill_indexes:indexes numbers: ',  &
                           (this%indxs(i),i=1,this%nindx )
endif

call QExit('fill_indexes_array')
end subroutine fill_indexes_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes_2(this,indx1,indx2)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indx1,indx2
call QEnter('fill_indexes_2')

if (.not.associated(this)) then
 call Quit('fill_indexes_2: entering pointer NOT associated !')
endif

if (this%nindx /= 2) then
 call Quit('fill_indexes_2: wrong number of indexes !')
endif

this%indxs(1) = indx1
this%indxs(2) = indx2

call QExit('fill_indexes_2')
end subroutine fill_indexes_2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes_4(this,indx1,indx2,indx3,indx4)
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indx1,indx2,indx3,indx4
call QEnter('fill_indexes_4')

if (.not.associated(this)) then
 call Quit('fill_indexes_4: entering pointer NOT associated !')
endif

if (this%nindx /= 4) then
 write(LUPRI,'(2x,a,i8)') 'fill_indexes_4: wrong "this%nindx" =',this%nindx 
 call Quit('fill_indexes_4: wrong number of indexes !')
endif

this%indxs(1) = indx1
this%indxs(2) = indx2
this%indxs(3) = indx3
this%indxs(4) = indx4

call QExit('fill_indexes_4')
end subroutine fill_indexes_4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_it(this)
!----------------------------------------------------------------------------
!
!        prints allocated scalar indexed quantity - all attributes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
call QEnter('print_it')

if (.not.associated(this)) then
 write(LUPRI,'(/,2x,a,/)')    & 
 'fill_indexes: pointer of "type(indexed_term)" NOT associated !'
 call Quit('fill_indexes: pointer of "type(indexed_term)" NOT associated !')
else
 write(LUPRI,'(3x,a,i1,a,i1,1x,$)')  'indx.term (ID:',this%id,  &
 '); # indexes:',this%nindx
 !write(LUPRI,'(1x,a,6i2,a,$)') ' ::',( this%indxs(i),i=1,this%nindx )
 if (this%nindx == 2) then
   write(LUPRI,'(1x,a,2i2,1x,a,$)') '>>',( this%indxs(i),i=1,this%nindx ),  &
   'label:'
 else if (this%nindx == 4) then
   write(LUPRI,'(1x,a,4i2,1x,a,$)') '>>',( this%indxs(i),i=1,this%nindx ),  &
   'label:'
 else if (this%nindx == 6) then
   write(LUPRI,'(1x,a,6i2,1x,a,$)') '>>',( this%indxs(i),i=1,this%nindx ),  &
   'label:'
 else
   write(LUPRI,'(2x,a,i6)') 'print_it: dead branch; "this%nindx"=',this%nindx
   call Quit('print_it: dead branch for the printout...')
 endif

 if (is_string(this%it_name)) then
   !write(LUPRI,'(3x,a,$)') 'indexed_quantity text description:'
   ! write(LUPRI,'(12x,a,$)') 'label:'
   call print_string(this%it_name)
 endif
 
endif

call QExit('print_it')
end subroutine print_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_IndxTerms(it1,it2) result(equal)
!----------------------------------------------------------------------------
!
!     Compare two indexed terms. Returns ".true." if they are equal.
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: it1, it2
logical :: equal, equal_indxs 
call QEnter('equal_IndxTerms')

equal = .false.

if (associated(it1).and. associated(it2)) then

 if (it1%id==it2%id.and.it1%nindx==it2%nindx) then
   if (allocated(it1%indxs).and.allocated(it2%indxs)) then
     equal_indxs = .true.
     do i = 1, it1%nindx
        if (it1%indxs(i)==it2%indxs(i)) then
            equal_indxs = equal_indxs .and. .true.
        else
            equal_indxs = equal_indxs .and. .false.
        endif
     enddo
     equal = equal_indxs
   endif
 endif

endif

call QExit('equal_IndxTerms')
end function equal_IndxTerms
!#endif
End Module Indexed_Quantities

Module Operators_string
! deactivated for ifort 10.1.015, higher versions working...
!#if !defined (VAR_IFORT)
! ... string of operators ... 
! ... kronecker deltas, creation/anihilation operators, indexed terms
! ...  
use KroneckerDelta
use SecQuantOperator
use Indexed_Quantities
use Control_print
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type string_member
!---------------------------------------------------------------------------------
! ... one string member (creat/anihil. operator, Kronecker delta, indexed term...)
!---------------------------------------------------------------------------------
  type(SQOper),       pointer ::    pSQOper => null()
  type(KronDelta),    pointer :: pKronDelta => null()
  type(indexed_term), pointer :: pIndxTerm  => null()
end type string_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type dll_member
!---------------------------------------------------------------------------------
!
!    ... one member with 2 pointers (to be in the doubly linked list) ...
!
!---------------------------------------------------------------------------------
 type(string_member) :: member
 type(dll_member), pointer :: previous => null(), next => null()
end type dll_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type one_string
!------------------------------------------------------------------
!
!        ... one string of mixed members ("dll_member")           
!
!------------------------------------------------------------------
 real(kind=8) :: factor = 0.0D0  ! default factor ...
! ... head of the DLL of mixed elements
 type(dll_member), pointer :: head => null() 
! ... rear of the DLL of mixed elements
 type(dll_member), pointer :: rear => null()
! ... number of elements in the DLL of "dll_member" type elements
 integer(kind=4) :: n = 0
end type one_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type head_cell
!------------------------------------------------------------------
!
!             ... member of cell of HEADS in DLL
!
!------------------------------------------------------------------
! ... content - pointer to one_string
 type(one_string), pointer :: p_string => null()
! ... pointers to previous, next head cells 
 type(head_cell), pointer :: previous => null()
 type(head_cell), pointer ::     next => null()
end type head_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type heads_list
!------------------------------------------------------------------------
!
!             DLL of heads pointing to strings; 
!
!------------------------------------------------------------------------
! ...  head
 type(head_cell), pointer :: main_head => null()
! ... tail
 type(head_cell), pointer :: main_rear => null()
! ... number of elements in the DLL of "head_cell" type elements
 integer(kind=4) :: n = 0

end type heads_list

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_sqo_cell(i,is_creat) result (this)
!----------------------------------------------------------------------------
!       Returns ALLOCATED second quantized operator cell
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i
logical, optional, intent(in) :: is_creat
type(dll_member), pointer :: this
type(SQOper), pointer :: this_sqo
call QEnter('make_sqo_cell')

this_sqo => make_sqo(i,is_creat)

nullify(this)
allocate(this,stat=ierr)
if (ierr /= 0) then
   write(LUPRI,'(/,2x,a,i4)')  'make_sqo_cell: allocation error ! ierr=',ierr
   call Quit('make_sqo_cell: allocation error !')
endif

this%member%pSQOper  =>  this_sqo

call QExit('make_sqo_cell')
end function make_sqo_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_sqo_cell(this)
!----------------------------------------------------------------------------
!
!      DEALLOCATES sec.quantiz.operator memory cell; returns null
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: this
call QEnter('deallocate_sqo_cell')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)')  &
  'deallocate_sqo_cell: sqo cell pointer NOT associated !'
else

 if (.not.associated(this%member%pSQOper)) then
   write(LUPRI,'(/,2x,a)')  &
   'deallocate_sqo_cell: sqo cell pointer NOT associated !'
 else
  !  deallocate(this%member%pSQOper, stat=ierr)
   ! if (ierr /= 0) then
    !  write(LUPRI,'(/,2x,a,i3)')  &
    !  'deallocate_sqo_cell: error in deallocation of' &
    !  //'"this%member%pSQOper";ierr=',ierr
    !  call Quit('deallocate_sqo_cell: error in deallocation of'&
    !  //'"this%member%pSQOper')
    !endif
    !nullify(this%member%pSQOper)
    call deallocate_sqo(this%member%pSQOper)
 endif

 deallocate(this, stat=ierr)
 if (ierr /= 0) then
   write(LUPRI,'(/,2x,a,i3)')  &
  'deallocate_sqo_cell: error in deallocation of' &
  //'"this";ierr=',ierr
  call Quit('deallocate_sqo_cell: error in deallocation of' &
  //'"this"')
 endif
 nullify(this)
endif

call QExit('deallocate_sqo_cell')
end subroutine deallocate_sqo_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_kd_cell(i,j) result (this)
!----------------------------------------------------------------------------
!
!    Returns cell "this" with the Kronecker delta  "d_ij"
!
!----------------------------------------------------------------------------
use KroneckerDelta
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i,j
type(dll_member), pointer :: this
type(KronDelta), pointer :: this_kd
call QEnter('make_kd_cell')

! ... true Kronecker delta ...
this_kd => make_kron_delta(i,j)

nullify(this)
allocate(this,stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a,i5)')  &
  'make_kd_cell: error in allocation; ierr=  ',ierr
  call Quit('make_kd_cell: error in allocation !')
endif

this%member%pKronDelta => this_kd

! ... check  for older compilers, like pgf90 5.2-4 (opteron64, Tel Aviv)
if (associated(this%member%pSQOper)) then
  if (get_iprint() >= 55) then
    write(LUPRI,'(2x,a)') 'make_kd_cell: pSQOper randomly associated !'
  endif
  this%member%pSQOper => null()
endif

if (associated(this%member%pIndxTerm)) then
  if (get_iprint() >= 55) then
    write(LUPRI,'(2x,a)') 'make_kd_cell: pIndxTerm randomly associated !'
  endif
  this%member%pIndxTerm   =>   null()
endif

if ( get_iprint() >= 35 ) then
  write(LUPRI,'(2x,a,$)') 'make_kd_cell: created Kronecker delta:'
  call print_one_member(this)
 ! call print_one_member(temp)
endif

call QExit('make_kd_cell')
end function make_kd_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_kd_cell(this)
!---------------------------------------------------------------------------
!
!                DEALLOCATE kronecker delta cell
!
!---------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: this
call QEnter('deallocate_kd_cell')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)')  &
  'deallocate_kd_cell: Kron.delta cell pointer NOT associated !'
else

 if (.not.associated(this%member%pKronDelta)) then
   write(LUPRI,'(/,2x,a)')  &
   'deallocate_kd_cell: Kron.delta cell pointer NOT associated !'
 else
   ! deallocate(this%member%pKronDelta, stat=ierr)
   ! if (ierr /= 0) then
    !  write(LUPRI,'(/,2x,a,i3)')  &
    !  'deallocate_kd_cell: error in deallocation of' &
    !  //'"this%member%pKronDelta";ierr=',ierr
    !  call Quit('deallocate_kd_cell: error in deallocation of'&
    !  //'"this%member%pKronDelta')
    !endif
    !nullify(this%member%pKronDelta)
     call deallocate_kd(this%member%pKronDelta)
 endif

 deallocate(this, stat=ierr)
 if (ierr /= 0) then
   write(LUPRI,'(/,2x,a,i3)')  &
  'deallocate_kd_cell: error in deallocation of' &
  //'"this";  ierr=',ierr
  call Quit('deallocate_sqo_cell: error in deallocation of' &
  //'"this"')
 endif
 nullify(this)
endif

call QExit('deallocate_kd_cell')
end subroutine deallocate_kd_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_it_cell(this)
!---------------------------------------------------------------------------
!
!                  DEALLOCATE "indexed_quantity" memory cell 
!
!---------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: this
call QEnter('deallocate_it_cell')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)')  &
  'deallocate_it_cell: IndxTerm cell pointer NOT associated !'
else

 if (.not.associated(this%member%pIndxTerm)) then
   write(LUPRI,'(/,2x,a)')  &
   'deallocate_it_cell: Kron.delta cell pointer NOT associated !'
 else
    call cancel_it(this%member%pIndxTerm)
 endif

 deallocate(this, stat=ierr)
 if (ierr /= 0) then
   write(LUPRI,'(/,2x,a,i3)')  &
  'deallocate_it_cell: error in deallocation of' &
  //'"this";  ierr=',ierr
  call Quit('deallocate_it_cell: error in deallocation of' &
  //'"this"')
 endif
 nullify(this)
endif

call QExit('deallocate_it_cell')
end subroutine deallocate_it_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function create_string(head, factor) result(this)
!---------------------------------------------------------------------------
!
!  constructor, allocates and returns "one_string" cell
!
!  entering "head" - pointer to the first element in string 
!
!  entering "factor" - numerical factor attached to the string (ie +1,-1,...)
!
!---------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer, optional :: head
real(kind=8), intent(in), optional :: factor 
type(one_string), pointer :: this
call QEnter('create_string')

   nullify(this)
   allocate(this, stat = ierr )
   if (ierr /= 0) then
     write(LUPRI,'(/,2x,a,i5)')  &
     'create_string: error in allocation; ierr=',ierr
     call Quit('create_string: error in allocation !')
   endif

   if (present(head)) then
     !if (associated(head))  this%head = head
     this%head => head
   endif

   if (present(factor)) then
     this%factor = factor
   else
     this%factor = 1.0D0
   endif

! ... rather secure zero number of elements and the tail
   this%n = 0
   this%rear => null()  

call QExit ('create_string')
end function create_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
recursive subroutine Insert_SQO_rec(this, head)
!--------------------------------------------------------------------------------------
!
! ... adds second quantized operator to the end of the doubly linked list of elements
!
!--------------------------------------------------------------------------------------
 type(SQOper), intent(inout) :: this 
 type(dll_member), pointer :: head

 if (.not. associated(head)) then
    allocate(head); nullify(head%next) 
    allocate(head%member%pSQOper)
    head%member%pSQOper = make_sqo(this%p,this%is_creation_oper)
    print *,'Insert_SQO: element inserted'; call print_sqo(head%member%pSQOper)
 else
    call Insert_SQO_rec(this,head%next)
 endif 
 end subroutine Insert_SQO_rec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Insert_SQOp_end(this_string,this_oper)
!-------------------------------------------------------------------------------
!
!  insert second quantized operator at the end of the string
!
!  this_string
!
!  this_oper
!
!-------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
type(SQOper), pointer :: this_oper 
type(dll_member), pointer :: temp, temp1

call QEnter('Insert_SQOp_end')

! .. entering string must be allocated - check it !
if (.not.associated(this_string)) then
  write(LUPRI,*) 'Insert_SQOp_end: entering pointer NOT associated ! Abort !'
  call Quit('Insert_SQOp_end: entering pointer NOT associated ! Abort !')
endif

! entering second quantized operator must be alocated - check it !
if (.not.associated(this_oper)) then
  write(LUPRI,'(2x,a)') 'Insert_SQOp_end: entering sec.quant.oper. NOT associated ! Abort !'
  call Quit('Insert_SQOp_end: entering pointer of SQOp NOT associated ! Abort !')
else
! ... operator is associated, now check if it's properly filled 
 if (.not.is_sqo_filled(this_oper)) then
  write(LUPRI,'(2x,a)') 'Insert_SQOp_end: operator not properly filled !'
  call print_sqo(this_oper)
  call Quit('Insert_SQOp_end: operator not properly filled !')
 endif
endif

if (.not.associated(this_string%head)) then

!  ... catch error if there is some
 if (associated(this_string%rear) .or. (this_string%n > 0) ) then
   write(LUPRI,'(/,2x,a)') 'Some error in the empty entering string :'
   write(LUPRI,'(2x,a,l1)')        &
   '...associated(this_string%head)=',associated(this_string%head)
   write(LUPRI,'(2x,a,l1)')        &
   '...associated(this_string%rear)=',associated(this_string%rear)
   write(LUPRI,'(2x,a,i3)') 'number of string elements, this_string%n=',this_string%n
   call Quit('Some error in the empty entering string !')
 endif

! ... the string is empty yet ,ie head => null()  ...
  allocate (this_string%head)
  !this_string%head%member%pSQOper = this_oper
  this_string%head%member%pSQOper => this_oper
  this_string%rear => this_string%head
  this_string%n = 1
else
  ! ... string contains at least one member; add element to the rear
  allocate(temp)  
  temp%previous => this_string%rear
  this_string%rear%next => temp
  temp%next => null()
  this_string%rear => temp
  this_string%rear%member%pSQOper => this_oper ! ... points to the operator
  this_string%n = this_string%n + 1
endif

call QExit('Insert_SQOp_end')
end subroutine Insert_SQOp_end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Insert_KD_end(this_string,this_kd)
!-------------------------------------------------------------------------------
!
!             inserts Kronecker delta at the end of the string
!
!-------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
type(KronDelta), pointer :: this_kd
type(dll_member), pointer :: temp, temp1
call QEnter('Insert_KD_end')

! .. entering string must be allocated - check it !
if (.not.associated(this_string)) then
  write(LUPRI,*) 'Insert_KD_end: entering pointer NOT associated ! Abort !'
  call Quit('Insert_KD_end: entering pointer NOT associated ! Abort !')
endif

! entering second quantized operator must be alocated - check it !
if (.not.associated(this_kd)) then
  write(LUPRI,'(2x,a)') 'Insert_KD_end: entering sec.quant.oper. NOT associated ! Abort !'
  call Quit('Insert_KD_end: entering pointer of SQOp NOT associated ! Abort !')
else
! ... operator is associated, now check if it's properly filled 
! if (.not.is_sqo_filled(this_oper)) then
!  write(LUPRI,'(2x,a)') 'Insert_KD_end: operator not properly filled !'
!  call print_sqo(this_oper)
!  call Quit('Insert_KD_end: operator not properly filled !')
! endif
endif

if (.not.associated(this_string%head)) then

 if (associated(this_string%rear) .or. (this_string%n > 0) ) then
   write(LUPRI,'(2x,a)') 'Some error in the empty entering string !'
   write(LUPRI,'(2x,a,l1)')        &
   '...associated(this_string%rear)=',associated(this_string%rear)
   write(LUPRI,'(2x,a,i3)') 'this_string%n=',this_string%n
   call Quit('Some error in the empty entering string !')
 endif


! ... the string is empty yet ,ie head => null()  ...
  allocate (this_string%head)
  this_string%head%member%pKronDelta => this_kd
  this_string%rear => this_string%head
  this_string%n = 1
else
  ! ... string contains at least one member; add element to the rear
  allocate(temp)  
  temp%previous => this_string%rear
  this_string%rear%next => temp
  temp%next => null()
  this_string%rear => temp
  this_string%rear%member%pKronDelta => this_kd ! ... points to the operator
  this_string%n = this_string%n + 1
endif

call QExit('Insert_KD_end')
end subroutine Insert_KD_end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Insert_IndxTerm_end(this_string,this_it)
!------------------------------------------------------------------------
!
!          Inserts "indexed term" at the end of the string
!
!-------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
type(indexed_term), pointer :: this_it
type(dll_member), pointer :: temp, temp1
call QEnter('Insert_IndxTerm_end')

! .. entering string must be allocated - check it !
if (.not.associated(this_string)) then
  write(LUPRI,*) 'Insert_IndxTerm_end: entering pointer NOT associated ! Abort !'
  call Quit('Insert_IndxTerm_end: entering pointer NOT associated ! Abort !')
endif

! entering second quantized operator must be alocated - check it !
if (.not.associated(this_it)) then
  write(LUPRI,'(2x,a)') 'Insert_IndxTerm_end: entering sec.quant.oper. NOT associated ! Abort !'
  call Quit('Insert_IndxTerm_end: entering pointer of IT NOT associated ! Abort !')
endif

if (.not.associated(this_string%head)) then

 if (associated(this_string%rear) .or. (this_string%n > 0) ) then
   write(LUPRI,'(2x,a)') 'Some error in the empty entering string !'
   write(LUPRI,'(2x,a,l1)')        &
   '...associated(this_string%rear)=',associated(this_string%rear)
   write(LUPRI,'(2x,a,i3)') 'this_string%n=',this_string%n
   call Quit('Some error in the empty entering string !')
 endif

! ... the string is empty yet ,ie head => null()  ...
  allocate (this_string%head)
  this_string%head%member%pIndxTerm => this_it
  this_string%rear => this_string%head
  this_string%n = 1
else

  ! ... string contains at least one member; add element to the rear
  allocate(temp)  
  temp%previous => this_string%rear
  this_string%rear%next => temp
  temp%next => null()
  this_string%rear => temp
  this_string%rear%member%pIndxTerm => this_it ! ... points to the operator
  this_string%n = this_string%n + 1

endif

call QExit('Insert_IndxTerm_end')
end subroutine Insert_IndxTerm_end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 recursive subroutine Insert_KD_rec(this, head)
!-------------------------------------------------------------------------------
!
! ... inserts kron.delta into doubly linked list of mixed elements
!
!--------------------------------------------------------------------------------
 type(KronDelta), intent(inout) :: this 
 type(dll_member), pointer :: head

 if (.not. associated(head)) then
    allocate(head); nullify(head%next) 
    allocate(head%member%pKronDelta)
    head%member%pKronDelta = make_kron_delta(this%p,this%q)
    !head%member%pKronDelta => make_kron_delta(this%p,this%q)
    print *,'Insert_KD: element inserted'; call print_kd(head%member%pKronDelta)
 else
    call Insert_KD_rec(this,head%next)
 endif 
 end subroutine Insert_KD_rec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_all_strings(heads)
!------------------------------------------------------------------------
!
!       print all strings under list "heads"
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!type(heads_list), intent(inout) :: heads
type(heads_list), intent(in) :: heads
type(head_cell), pointer :: temp
call QEnter('print_all_strings')

write(LUPRI,'(/,2x,a)')  &
'***  print_all_strings: printing all strings  ***'

if (associated(heads%main_head)) then

 i = 1
! ... cycle ...
 temp => heads%main_head
 do while (associated(temp))
  write(LUPRI,'(/,5x,a,i3)') '.... *** string #',i
  call print_string_elem(temp%p_string)
  temp => temp%next
  i = i + 1
 enddo
else
  write(LUPRI,'(2x,a)') & 
  'print_all_strings: entering heads%main_head poiner not' &
   //' associated !'
endif

call QExit('print_all_strings')
end subroutine print_all_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_string_elem(p_this_string)
!------------------------------------------------------------------------
!
!   print-out string of elements of the string pointed with "p_this_string"
! together with the factor 
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: p_this_string
call QEnter('print_string_elem')

if (.not.associated(p_this_string)) then
  write(LUPRI,'(/,2x,a)') 'print_string_elem: entering "p_this_string" not associated !' 
  call Quit('print_string_elem: entering "p_this_string" not associated !')
else
  write(LUPRI,'(/,a,i3,a,f12.5)')   & 
  '*** printing  one string of ',p_this_string%n, &
  ' elements and with the attached factor',p_this_string%factor
  call print_members(p_this_string%head)
endif

call QExit('print_string_elem')
end subroutine print_string_elem

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_string_reverse(this_string)
!--------------------------------------------------------------------------
!
!   print-out string of elements of the string pointed with 'this_string'
! together with the factor 
!
!--------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string

call QEnter('print_string_reverse')

if (.not.associated(this_string)) then
  write(LUPRI,'(/,2x,a)') 'print_string_reverse: entering "this_string" not associated !' 
  call Quit('print_string_reverse: entering "this_string" not associated !')
else
  write(LUPRI,'(/,a,f12.5)')   & 
  '*** printing string elements in reverse order; factor:',this_string%factor
  call print_members_reverse(this_string%rear)
endif

call QExit('print_string_reverse')
end subroutine print_string_reverse

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_members(head)
!------------------------------------------------------------------------
!
! prints one string of members (operators, kron.deltas, indexed terms) ...
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: head, temp
call QEnter('print_members')

write(LUPRI,'(a)') '* individual elements in the string:'

if (.not.associated(head)) then
  write(LUPRI,'(2x,a)')   &
 "...string of members (operators, kron.deltas, indexed terms) empty"
else
   temp => head
   do while (associated(temp))
     call print_one_member(temp)
     temp => temp%next
   enddo
endif

call QExit('print_members')
end subroutine print_members

subroutine print_members_reverse(rear)
!------------------------------------------------------------------------
!
!                  print string elements from the rear
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer  :: rear, temp
call QEnter('print_members_reverse')

write(LUPRI,'(a)') & 
'* individual elements of the string in the reverse order:'

if (.not.associated(rear)) then
  write(LUPRI,'(2x,a)')   &
  "...string of members (operators, kron.deltas, indexed terms) empty"
else
   temp => rear 
   do while (associated(temp))
     call print_one_member(temp)
     temp => temp%previous
   enddo
endif

call QExit('print_members_reverse')
end subroutine print_members_reverse

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_one_member(this)
!------------------------------------------------------------------
!
!   print-out one member according to the entering pointer "this"
!
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: this

call QEnter('print_one_member')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)')   & 
  'print_one_member: entering pointer is not associated !'
  !call Quit('print_one_member: entering pointer is not associated !')
else

if (associated(this%member%pSQOper)) then
   call print_sqo(this%member%pSQOper)
else if (associated(this%member%pKronDelta)) then
  ! write(LUPRI,'(2x,a)') 'this%member%pKronDelta !'
   call print_kd(this%member%pKronDelta)
else if (associated(this%member%pIndxTerm)) then
   call print_it(this%member%pIndxTerm)
else
   write(LUPRI,'(/,2x,a)') '*** print_one_member dead branch output ***'
   write(LUPRI,'(2x,a,l1)')   &
   'associated(this%member%pSQOper):',associated(this%member%pSQOper)
   write(LUPRI,'(2x,a,l1)')   &
   'associated(this%member%pKronDelta):',associated(this%member%pKronDelta)
   write(LUPRI,'(2x,a,l1)')   &
   'associated(this%member%pIndxTerm):',associated(this%member%pIndxTerm)
   call Quit('print_one_member: not SQOp, not Kron.delta, not indexed term !')
endif

endif

call QExit('print_one_member')

end subroutine print_one_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_strings(heads)
!------------------------------------------------------------------------
!
!         most important routine - process all strings of elements
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(heads_list), intent(inout) :: heads
type(head_cell), pointer :: temp, temp_prev
!logical :: equal_head_cells
call QEnter('process_strings')

if (get_iprint() >= 5) then
   call Header('*** Output from "process_strings" ***',-1)
 !  write(LUPRI,'(/,2x,a)') 'in process_strings...'
endif

! ... go string after string ...

if (.NOT.associated(heads%main_head)) then
  write(LUPRI,'(/,2x,a)')   &
  'process_strings: heads%main_head NOT associated !'
  call Quit('process_strings: heads%main_head NOT associated !')
endif

i = 1
! ... cycle ... string after string ...
temp => heads%main_head
do while (associated(temp))
  temp_prev => temp ! save the last address
  if (get_iprint() >= 5) then
    write(LUPRI,'(/,2x,a,i5)')  &
    '...going to process string in the list #',i
  endif
  !call process_one_string(temp%p_string)
  call process_one_string(temp)
  temp => temp%next
  i = i + 1
enddo

! ... store (and check) total number of strings

heads%n = i-1
if (heads%n <= 0) then
 write(LUPRI,'(/,2x,a,i6)') 'process_strings: wrong number of strings! i=',i
 call Quit('process_strings: wrong number of strings!')
endif

! ... save the rear ...
heads%main_rear => temp_prev
!if (heads%main_rear .eq. heads%main_head) then
if (equal_head_cells(heads%main_rear,heads%main_head)) then
 if (heads%n  /= 1) then
   write(LUPRI,'(2x,a,i7)')   & 
   'process_strings: something wrong about heads%n=',heads%n
   call Quit('process_strings: something wrong about heads%n !')
 endif
endif

!if (get_iprint() >= 12) then
if (get_iprint() >= 9) then
   write(LUPRI,'(/,2x,a,i5)') '...number of strings:',heads%n
   write(LUPRI,'(/,8x,a)')  'first string in the list... '
   call print_string_elem(heads%main_head%p_string)

   write(LUPRI,'(/,8x,a)')  'last string in the list... '
   !call print_string_elem(temp_prev%p_string)
   call print_string_elem(heads%main_rear%p_string)

!   call Quit('...in process_strings')
endif


if (get_iprint() >= 7) then
   call Header('**** process_strings: printing all strings ****',-1)
   write(LUPRI,'(/,2x,a)') '**** process_strings: printing all strings ****'
   call print_all_strings(heads)
endif

! ... process all strings for removing those giving zero in expect.value
!if ( is_spinors_specif() ) then
 ! write(LUPRI,*) 'spinors specification ! highest occup.',get_upper_occ()
!  call process_exp_val_strings(heads)
!endif

call QExit('process_strings')
end subroutine process_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function is_rearranged_string(this) result(is)
!------------------------------------------------------------------------
! Check whether the string is rearranged: check if  all creation operators 
! are on the left, all anihilation on the right
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
logical :: is, creat, anihil
type(dll_member), pointer :: temp
call QEnter('is_rearranged_string')
is = .false.

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a,/)')   &
  'is_rearranged_string: entering string NOT associated !'
  call Quit('is_rearranged_string: entering string NOT associated !')
endif

if (.not.associated(this%head)) then
! ... empty string ...
  if (this%n /= 0) then
    write(LUPRI,*) 'is_rearranged_string: head=NULL but this%n=',this%n
    call Quit('is_rearranged_string: head=NULL but this%n =/ 0 !!!')
  endif
  is = .true.
else
  if (get_iprint() >= 28) then
    write(LUPRI,'(/,2x,a)') 'is_rearranged_string: elements'
  endif
  temp => this%head; i=0

  anihil = .false.
   creat = .false.

  do while (associated(temp)) 

    if (get_iprint() >= 36) then
      call print_one_member(temp)
    endif

    if  (associated(temp%member%pSQOper)) then
      if (.not.is_sqo_filled(temp%member%pSQOper)) then
         write(LUPRI,'(/,2x,a,/)')  &
         'is_rearranged_string: SQoper NOT filled !'
         call Quit('is_rearranged_string: SQoper NOT filled !')
      endif
      if (get_iprint() >= 36) then
        write(LUPRI,'(2x,a,/)') 'we have sec. quantized operator !'
        !call print_sqo(temp%member%pSQOper)
      endif

      if (temp%member%pSQOper%is_creation_oper) then
         creat = .true.
         if (anihil) then
           is = .false.
           Exit
         endif
      else
         anihil = .true.
         is = .true.
      endif
    endif
    i = i + 1; temp => temp%next
  enddo

!--------------------------------------------------------------------------------
! ... string is rearranged also when there are no creat/anihil. operators,
!    or when there are only creation or only anihilation operators
!--------------------------------------------------------------------------------
  if (      (anihil.and..not.creat).or.  &
       (.not.anihil.and.     creat).or.  & 
       (.not.anihil.and..not.creat) ) is = .true.
endif

call QExit('is_rearranged_string')
end function is_rearranged_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine cancel_strings(heads)
!------------------------------------------------------------------------
! 
!           ...  deallocates all strings in the list
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(heads_list), intent(inout) :: heads
type(head_cell), pointer :: temp, temp_save
call QEnter('cancel_strings')

if (get_iprint() >= 5) then
  call Header('in cancel_strings...',-1)
 ! write(LUPRI,'(/,2x,a)') 'in cancel_strings...'
endif

i = 1
!   ... cycle ... string after string from the head...
temp => heads%main_head
do while (associated(temp))
  if (get_iprint() >= 5) then
     write(LUPRI,'(/,2x,a,i5)')  &
     '...going to deallocate string in the list #',i
  endif
  temp_save => temp
  call deallocate_one_string(temp)
  temp => temp_save%next
  i = i + 1
enddo

call QExit('cancel_strings')
end subroutine cancel_strings
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine deallocate_one_string(this)
!------------------------------------------------------------------------
!
!                  Deallocates one string; returns null
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(head_cell), pointer :: this
type(dll_member), pointer :: temp, temp_save
call QEnter('deallocate_one_string')

if (.not.associated(this)) then
  write(LUPRI,'(/,2x,a)') 'deallocate_one_string: string NOT allocated !'
  call Quit('deallocate_one_string: string NOT allocated !')
endif

if (.not.associated(this%p_string)) then
  write(LUPRI,'(/,2x,a)')  & 
  'deallocate_one_string: string "this%p_string" is NOT allocated !'
  call Quit('deallocate_one_string: string "this%p_string" is NOT allocated !')
endif

if (associated(this%p_string%head)) then
  temp => this%p_string%head
  do while (associated(temp))
    temp_save => temp
   
    if (associated(temp%member%pSQOper)) then
       call deallocate_sqo_cell(temp)
    else if (associated(temp%member%pKronDelta)) then
       call deallocate_kd_cell(temp)
    else if (associated(temp%member%pIndxTerm)) then
       call deallocate_it_cell(temp)
    else 
       call Quit('deallocate_one_string: error branch !')
    endif
    temp => temp_save%next
  enddo

endif

deallocate(this,stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a)') 'deallocate_one_string: error in deallocation !'
  call Quit('deallocate_one_string: error in deallocation of "this" !')
endif
nullify(this)

call QExit('deallocate_one_string')
end subroutine deallocate_one_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_one_string(this)
!------------------------------------------------------------------------
!
!  Process entering string : apply commutator relations
!
!     example:
! ----------------
!
!  P1 i a+ P2 = P1 ( d_ia - a+ i ) P2 = P1 d_ia P2 - P1 a+ i P2
!
!            at the beginning:   P1 i a+ P2
!
! after evaluating commutator:   P1 d_ia P2 - P1 a+ i P2
!
!           where P1, P2 are substrings (collection of elements)
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!type(one_string), pointer :: this
type(head_cell), pointer :: this, this2, this_temp
logical :: found_splitting
!type(dll_member), pointer :: temp
call QEnter('process_one_string')

if (.not.associated(this)) then
 write(LUPRI,'(/,2x,a)')  &
 'process_one_string: entering string for analyse is NOT allocated !'
 call Quit('process_one_string: entering string is NOT allocated !')
endif

if (.not.associated(this%p_string)) then
 write(LUPRI,*) 'process_one_string: this string of elements is NOT allocated !'
 call Quit('process_one_string: this string is NOT allocated !')
endif

if (get_iprint() >= 5) then
 write(LUPRI,'(/,2x,a)')  &
 'process_one_string: processing this ONE entering string...'
 call print_string_elem(this%p_string)
endif

! ... rearrange string - move deltas, indexed terms left
call rearrange_string(this%p_string)

if (get_iprint() >= 5) then
 write(LUPRI,'(/,2x,a)')  &
 'process_one_string: now REARRANGED string going to be processed...'
 call print_string_elem(this%p_string)
endif

! ... initialization ...
found_splitting = .true.
do while ( found_splitting )

  ! ... find first the splitting place 
  itype = -1; ipos = -1
  found_splitting = .false.
  call find_splitting(this%p_string, itype, ipos )

  ! ... allocate one more string, duplicate it and insert behind strings
  if (ipos > 0) then

    found_splitting = .true.

    if ( get_iprint() >= 5 ) then
     write(LUPRI,'(/,2x,a)')  &
     '*** process_one_string - creating next string as splitting place was found ! '
    endif

    ! ... insert the string after the main string, chain it
    allocate(this2)
    this2%p_string => null()

    this2%previous => this

    if (associated(this%next)) then
      this2%next => this%next
      this%next%previous => this2
    else
      nullify( this2%next) ! this2 is the last string
    endif

    this%next => this2

    call duplicate_string(this%p_string, this2%p_string )

    if (get_iprint() >= 7) then
       write(LUPRI,'(/,2x,a)')   &
       'process_one_string: duplicated string:'
       call print_string_elem(this2%p_string)
       if (associated(this2%next)) write(LUPRI,'(2x,a)')  &
       '...next string FOLLOWS ... '
    endif

    call process_2_strings(this%p_string, this2%p_string, ipos)

    ! ... rearrange both updated strings - move deltas, indexed terms left
    call rearrange_string( this%p_string)
    call rearrange_string(this2%p_string)

 endif

enddo

!call quit('process_one_string')

call QExit('process_one_string')
end subroutine process_one_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_2_strings(this_string, copied_string, ipos)
!------------------------------------------------------------------------
!
!  process 2 identical strings - built in commutator relation
! ---------------------------------------------------------------
!
!  P1 i a+ P2 = P1 ( d_ia - a+ i ) P2 = P1 d_ia P2 - P1 a+ i P2
!
!            at the beginning:   P1 i a+ P2  2 times
!
! after evaluating commutator:   P1 d_ia P2 - P1 a+ i P2
!
!             where P1, P2 are bunches of elements 
!
!  ipos - position of i operator in the entering strings
!                              "this_string, copied_string"
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: ipos
type(one_string), pointer ::  this_string, copied_string
!type(dll_member), pointer :: d_ia, temp_i, temp_a
type(dll_member), pointer :: d_ia, temp_i, temp_a
logical, save :: first_entry = .true.
call QEnter('process_2_strings')

if ( get_iprint() >= 5 ) then
  write(LUPRI,'(/,2x,a)')   &
  '*** process_2_strings: going to process 2 subsequent entering strings ***'
  call print_string_elem(this_string)
  call print_string_elem(copied_string)
endif

allocate(temp_i); allocate(temp_a)

temp_i = get_dll_member(this_string,ipos)

if (.not.associated(temp_i%member%pSQOper)) then
   call print_one_member(temp_i)
   call Quit(' process_2_strings: picked element not creat/anihil. operator!')
else if  (.not.(temp_i%member%pSQOper%p > 0) .or.   &
               (temp_i%member%pSQOper%is_creation_oper) ) then
   call print_one_member(temp_i)
   call Quit(' process_2_strings: picked element not creat/anihil. operator!')
endif

temp_a = get_dll_member(this_string,ipos+1)

if (.not.associated(temp_a%member%pSQOper)) then
   call print_one_member(temp_a)
   call Quit(' process_2_strings: picked element not creat/anihil. operator!')
else if(.not.(temp_a%member%pSQOper%p > 0) .or.  &
        .not.(temp_a%member%pSQOper%is_creation_oper) ) then
   call print_one_member(temp_a)
   call Quit(' process_2_strings: picked element not creat/anihil. operator!')
endif

if ( get_iprint() >= 5 ) then
   write(LUPRI,'(/,2x,a,2i3)')  &
   'process_2_strings: picked 2 operators p q+ from the string, positions ', &
   ipos, ipos+1
   call print_one_member(temp_i)
   call print_one_member(temp_a)
endif

write(LUPRI,*) 'd_ia made before: '; call print_string_elem(this_string)
! ... create kronecker delta ...
!d_ia => null()
nullify(d_ia)
allocate(d_ia,stat=ierr)
if (ierr /= 0) then
 write(LUPRI,'(/,2x,a)') 'process_2_strings: error in allocation; ierr=',ierr
 call Quit('process_2_strings: error in allocation!')
endif

d_ia = make_kd_cell(temp_i%member%pSQOper%p, temp_a%member%pSQOper%p)
  write(LUPRI,'(1x,a,$)') 'allocated d_ia:' ;call print_one_member(d_ia)
write(LUPRI,*) 'd_ia made after: '; call print_string_elem(this_string)

! ... adapt first string ...
call remove_element(this_string, ipos)

call insert_element(this_string, ipos, d_ia)
deallocate(d_ia,stat=ierr)
if (ierr /= 0) then
 write(LUPRI,'(/,2x,a)') 'process_2_strings: error in "d_ia" dallocation; ierr=',ierr
 call Quit('process_2_strings: error in "d_ia" deallocation!')
endif

call remove_element(this_string, ipos)

! ... adapt second, copied string ...
call exchange_string_elements(copied_string, ipos, ipos+1) 

! ... change the factor ...
copied_string%factor = - copied_string%factor

if ( get_iprint() >= 5 )  then
 write(LUPRI,'(/,2x,a)')  &
 '*** process_2_strings: adapted 2 subsequent strings'
 call print_string_elem(this_string)
 call print_string_elem(copied_string)
! call quit('...in process_2_strings')
 call FLSHFO(LUPRI)
endif

call QExit('process_2_strings')
end subroutine process_2_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_dll_member(this_string,ipos) result(elem)
!------------------------------------------------------------------------
!
!     returns element from the string at the position 'ipos'
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: ipos
type(dll_member), target :: elem
type(dll_member), pointer :: temp
type(one_string), pointer ::  this_string
call QEnter('get_dll_member')

if ( get_iprint() >= 35 ) then
  write(LUPRI,'(/,2x,a,i3)') '*** get_dll_member *** : wanted element ',ipos
endif

if (.not.associated(this_string)) then
   write(LUPRI,'(/,2x,a)')   &
   'get_dll_member: entering string is NOT associated !'
   call Quit('get_dll_member: entering string is NOT associated !')
else
  if (this_string%n < ipos .or. ipos <= 0) then
    write(LUPRI,'(/,2x,a,i3,a,i3)')   &
    'get_dll_member: number of string elements=',this_string%n,  &
    ' while wrong position in the string, ipos=',ipos
    write(LUPRI,'(2x,a)') 'get_dll_member: control print-out of the whole string'
    call print_string_elem(this_string)
    call Quit('get_dll_member: wrong position of the element in the string !')
  endif
endif

temp => this_string%head; i = 0
do while(associated(temp))
  i = i + 1
  if (i == ipos) then
    elem = temp
    Exit ! leave the cycle when 'ipos'. element is catched
  endif
  temp => temp%next
enddo

if (get_iprint() >= 37) then
  write(LUPRI,'(2x,a,i3,a,$)')  &
  'get_dll_member: picked element from postion ',ipos,' ...'
  call print_one_member(temp) 
endif

call QExit('get_dll_member')
end function get_dll_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine duplicate_string(this_string, copied_string)
!------------------------------------------------------------------------
!
!  duplicate string - allocation of memory for the new one is inside
!
!  on entry: this_string
!
!  on output: copied_string 
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer ::  this_string, copied_string
type(dll_member), pointer :: temp, temp1, temp2
call QEnter('duplicate_string')

if ( get_iprint() >= 25 ) then
  write(LUPRI,'(/,2x,a)') '*** duplicate_string ***'
endif

if (.not.associated(this_string)) then
  write(LUPRI,'(2x,a)') 'entering string NOT allocated ! Exit !'
  call Quit('duplicate_string : entering string NOT allocated !')
else
  if (.not.associated(this_string%head)) then
    write(LUPRI,'(2x,a)') 'entering string  is empty! Exit !'
    call Quit('duplicate_string : entering string is empty !')
  endif 
endif

allocate(copied_string,stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a,i4)')  'duplicate_string: allocation error ! ierr=',ierr
  call Quit('duplicate_string: allocation error !')
endif

temp1 => copied_string%head;  temp2 => temp
temp => this_string%head;  
i = 0
do while (associated(temp))

  i = i + 1

  call insert_member(temp,temp1)

  if (get_iprint() >= 27) then
    write(LUPRI,'(/,2x,a,i3,a,$)')  &
    'duplicate_string: ',i,'. old member element:'
    call print_one_member(temp)
    write(LUPRI,'(2x,a,$)') '              the same new member element:'
    call print_one_member(temp1)
  endif

  if (i == 1) then
    copied_string%head => temp1; temp1%previous => null(); temp1%next => null()
  else if (i > 1) then
    temp2%next => temp1; temp1%previous => temp2
  endif

  temp2 => temp1  ! ... save previous element
  temp => temp%next
enddo

copied_string%rear  =>  temp1
copied_string%n  =  this_string%n
copied_string%factor  =  this_string%factor

! ... after copying element after element print out whole new string
if (get_iprint() >= 25) then
  write(LUPRI,'(/,2x,a)')  &
  '*** duplicate_string: allocated duplicated string:'
  call print_string_elem(copied_string)
endif

!call quit('end in duplicate_string')

call QExit('duplicate_string')
end subroutine duplicate_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine insert_member(member_exist,member_new)
!------------------------------------------------------------------------
! Allocates dll_member "member_new" and copies "member_exist" into it.
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: member_exist, member_new
logical :: is_operator,  is_term 
if (get_iprint() >= 37) then
  write(LUPRI,'(/,2x,a)') 'insert_member: duplicating string elements'
endif
if (.not.associated(member_exist)) then
   write(LUPRI,'(/,2x,a)') 'insert_member: entering element does NOT exist !'
   call Quit('insert_member: entering element does NOT exist !')
endif
allocate(member_new,stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a,i4)')  &
  'insert_member: "member_new" allocation error ! ierr=',ierr
  call Quit('insert_member: "member_new" allocation error !')
endif
! duplicated dll_member
call duplicate_dll_member(member_exist,member_new)
if (get_iprint() >= 37) then
   write(LUPRI,'(/,2x,a,$)') 'insert_member: old member element:'
   call print_one_member(member_exist)
   write(LUPRI,'(2x,a,$)') '           the same new member element:'
   call print_one_member(member_new)
endif
!call quit('insert_member')
end subroutine insert_member
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine duplicate_dll_member(member_exist,member_new)
!------------------------------------------------------------------------
!          copies dll_member:   member_exist -> member_new
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(dll_member), pointer :: member_exist, member_new
if (.not.associated(member_exist).or..not.associated(member_new)) then
 call Quit('duplicate_dll_member: member_exist or member_new NOT associated !')
endif
if (associated(member_exist%member%pSQOper)) then
    allocate(member_new%member%pSQOper,stat=ierr)
    call copy_sqo(member_exist%member%pSQOper,member_new%member%pSQOper)
else if (associated(member_exist%member%pKronDelta)) then
    allocate(member_new%member%pKronDelta,stat=ierr)
    call copy_kd(member_exist%member%pKronDelta,member_new%member%pKronDelta)
else if (associated(member_exist%member%pIndxTerm)) then
    allocate(member_new%member%pIndxTerm,stat=ierr)
    call copy_it(member_exist%member%pIndxTerm,member_new%member%pIndxTerm)
else
  call Quit('duplicate_dll_member: dead branch !')
endif
end subroutine duplicate_dll_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine rearrange_string(this)
!------------------------------------------------------------------------
!
!
!  rearrange string that all indexed terms and kronecker deltas are 
! on the left separated from the rest -  from a bunch 
! of creation/anihilation operators, which cannot be simply exchanged
!
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
type(dll_member), pointer :: temp, temp1
logical :: is_operator,  is_term 
call QEnter('rearrange_string')

if (get_iprint() >= 25) then
  write(LUPRI,'(/,2x,a)') 'going to rearrange the string...'
endif

if (.not.associated(this)) then
 write(LUPRI,'(2x,a)') 'rearrange_string: entering string NOT associated (allocated) !'
 call Quit('rearrange_string: entering string NOT associated (allocated) !')
endif

if (.not.associated(this%head)) then
  if (get_iprint() >= 25) then
    write(LUPRI,'(2x,a)')  & 
    'rearrange_string: entering string is empty ! leaving this routine'
  endif
else

 is_term = .true.

 do while (is_term)

    ! ... go element after element
    ipos = 0
    ipos_term = 0
    temp => this%head
    is_operator = .false. 
    is_term = .false.

    do while (  associated(temp) .and. .not.is_term   )
      ipos = ipos + 1
      if  (associated(temp%member%pSQOper)) is_operator = .true.
      if  (is_operator  .and.   &
          (associated(temp%member%pKronDelta).or.associated(temp%member%pIndxTerm))) then
          ipos_term = ipos; is_term = .true.
          if (get_iprint() >= 28) then
            write(LUPRI,'(2x,a,i3)')  & 
           'rearrange_string: catched element to be moved to the beginnig,'  &
         //' position', ipos_term
          endif
      endif
      temp => temp%next
    enddo

    if (is_term) then

       if (get_iprint() >= 28) then
         write(LUPRI,'(2x,a)') 'after the operator found term !'
         call get_string_element(this,ipos_term,temp1)
         call print_one_member(temp1)
       endif

!    ... move term to the beginning ...
       call move_element_beg(this,ipos_term)

       if ( get_iprint() >= 28 ) then
         write(LUPRI,'(2x,a)') 'term moved to the beginning  !'
         call print_string_elem(this) 
       endif

    endif

  enddo

endif

call QExit('rearrange_string')
end subroutine rearrange_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine move_element_beg(this_string,i)
!------------------------------------------------------------------------
!
!   move element in the string from  position "i" to the very beginning
!
! On input: this_string - pointer to the string
!    i - position of the element to be moved at the beginning
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
type(dll_member), pointer :: temp, temp1, temp2
call QEnter('move_element_beg')

!call check_string(this_string,i)

if (get_iprint() >= 35) then
  write(LUPRI,'(/2x,a,i3)')  &
  'in "move_element_beg:" to be moved element no:',i
endif

!-----------------------------------------------------------
!   ... move only the element which is NOT the first
!-----------------------------------------------------------
if ( i == 1) then

 if (get_iprint() >= 35) then
    write(LUPRI,'(/2x,a,i3)')  &
   'in "move_element_beg:" to be moved element no:',i
 endif

else  ! i>1 

call  get_string_element(this_string,i,temp)

! ... element before ...
temp1 => temp%previous
! ... element after ...
temp2 => temp%next

if (get_iprint() >= 35) then
 write(LUPRI,'(/,2x,a,i3)') '... catched this member no ',i
 call print_one_member(temp)
 write(LUPRI,'(2x,a,$)') '... the element before: '
 call print_one_member(temp1)
 write(LUPRI,'(2x,a,$)') '... the element after: '
 call print_one_member(temp2)
endif

if (get_iprint() >= 35) then
 write(LUPRI,'(/,2x,a)')   &
 '*** move_element_beg: control print-out before reordering the string:'
 call print_string_elem(this_string)
! call print_string_reverse(this_string)
endif

! stop 'here stop ...'

!---------------------------------------------------------------------
! ... move j-th element "temp" to the beginning: redirect pointers
!---------------------------------------------------------------------

temp1%next => temp2

if (associated(temp2)) then
  temp2%previous => temp1
endif

if (get_iprint() >= 105) then
 write(LUPRI,'(/,2x,a)')   &
 '*** move_element_beg: control print-out 1:'
 call print_one_member(temp1%next)
 call print_one_member(temp2%previous)
endif

! stop 'here stop 1...'

!  call print_string(this_string)
!  call print_string_reverse(this_string)

! stop 'here stop 2...'

temp%previous => null() !... first element in the DLL
temp%next => this_string%head

this_string%head%previous =>temp

this_string%head => temp ! ... finally set the head 

endif

if (get_iprint() >= 35) then
 write(LUPRI,'(/,2x,a)')   &
 'move_element_beg: control print-out of reordered string:'
 call print_string_elem(this_string)
endif

call QExit('move_element_beg')
end subroutine move_element_beg

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine remove_element(this_string,i)
!------------------------------------------------------------------------
!  remove i-th element from the entering string
!
!  On input: this_string - pointer to the string
!             i - position
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
type(dll_member), pointer :: temp, temp1, temp2

call get_string_element(this_string,i,temp)

! ... previous element
temp1 => temp%previous
! ... next element
temp2 => temp%next

if (get_iprint() >= 25) then
!if (get_iprint() >= 5) then
  write(LUPRI,'(/2x,a,i3)')  &
  '"remove_element": to be removed element no:',i
  call print_one_member(temp)
  write(LUPRI,'(2x,a,$)') '... the element before: '
  call print_one_member(temp1)
  write(LUPRI,'(2x,a,$)') '... the element after: '
  call print_one_member(temp2)
endif

deallocate(temp, stat = ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a,i5)')  &
  '"remove_element": error in deallocation ! ierr =',ierr
  call Quit('"remove_element": error in deallocation !')
endif

! ... decrease  number of elements ...
this_string%n  =  this_string%n - 1

if (associated(temp1).and.associated(temp2)) then
! ... in the between
  temp1%next => temp2
  temp2%previous => temp1
else if (associated(temp1).and..not.associated(temp2)) then
! ... last element 
   temp1%next => null()
   this_string%rear => temp1
else if (.not.associated(temp1).and.associated(temp2)) then
!... first element
   this_string%head => temp2
   temp2%previous => null()
else if (.not.associated(temp1).and..not.associated(temp2)) then
! ... the only one element in the string
   this_string%head => null()
   this_string%rear => null()
else
  call Quit('remove_element - dead branch !')
endif

if (get_iprint() >= 25) then
!if (get_iprint() >= 5) then
  write(LUPRI,'(/,2x,a)')  &
  'remove_element: printed updated string with removed element:'
  call print_string_elem(this_string)
endif
end subroutine remove_element

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine exchange_string_elements(this_string,i_in,j_in)
!------------------------------------------------------------------------
!
!       Exchange 2 arbitrary elements in the string, i-th and j-th.
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i_in,j_in
type(dll_member), pointer :: temp_i,temp_j,temp_i1,temp_i2,temp_j1,temp_j2
call QEnter('exchange_string_elements')

i = i_in; j = j_in

if (get_iprint() >= 25) then
  write(LUPRI,'(/,2x,a,2i3)') &
  '*** exchange_string_elements: going to exchange 2 string elements ',i,j
endif

if (i == j) then
  write(LUPRI,'(/,2x,a)') 'two identical elements ! Quit !'
  call Quit('two identical elements ! no possible exchanhe !')
else if (i > j) then

  if (get_iprint() >= 25) then
   write(LUPRI,'(/,2x,a)')   &
   'exchange_string_elements: i> j - going to exchange i,j'
  endif

  k = j; j = i; i = k  

endif 

call get_string_element(this_string,i,temp_i)

if (get_iprint() >= 27) then
   write(LUPRI,'(/,2x,a,i3,$)')   &
   'exchange_string_elements: element # ',i
   call print_one_member(temp_i)
endif

call get_string_element(this_string,j,temp_j)

if (get_iprint() >= 27) then
  write(LUPRI,'(/,2x,a,i3,$)')   &
  'exchange_string_elements: element # ',j
  call print_one_member(temp_j)
endif

temp_i1 => temp_i%previous 
temp_i2 => temp_i%next

temp_j1 => temp_j%previous
temp_j2 => temp_j%next

! --------------------------------------------
!         ... redirect pointers ...
! --------------------------------------------
if ( .not. (j == i + 1) ) then
! ... 2 elements, not one after another

temp_i%previous => temp_j1
temp_i%next => temp_j2

temp_j%previous => temp_i1
temp_j%next => temp_i2

if (associated(temp_i1)) then
! ... element before i-th element
 temp_i1%next => temp_j
else
! ... i is the first element, take care of the head 
 this_string%head => temp_j
endif

if (associated(temp_i2)) then
 temp_i2%previous => temp_j
endif

if (associated(temp_j1)) then
 temp_j1%next => temp_i
endif

if (associated(temp_j2)) then
 temp_j2%previous => temp_i
else
! ... j is the last element, take care of the rear 
  this_string%rear => temp_i
endif

else
!-------------------------------------------------
!      ... 2 elements, one after another
!-------------------------------------------------
if (associated(temp_i1)) then
! ... element before i-th element
 temp_i1%next => temp_j
else
! ... i is the first element, take care of the head 
 this_string%head => temp_j
endif

if (associated(temp_j2)) then
 temp_j2%previous => temp_i
else
! ... j is the last element, take care of the rear 
  this_string%rear => temp_i
endif

temp_i%previous => temp_j
temp_i%next => temp_j2

temp_j%previous => temp_i1
temp_j%next => temp_i

endif

if (get_iprint() >= 25) then
   write(LUPRI,'(/,2x,a)')   &
   '*** exchange_string_elements: printing updated string'
   call print_string_elem(this_string)
endif

call QExit('exchange_string_elements')
end subroutine exchange_string_elements

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_string_element(this_string,i,temp)
!------------------------------------------------------------------------
!        returns pointer to the i-th element in the string
!------------------------------------------------------------------------
implicit none
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
integer :: j
type(dll_member), pointer :: temp
call check_string(this_string,i)
j = 1
temp => this_string%head
do while (j < i)
  if (get_iprint() >= 35) then
    write(LUPRI,'(2x,a,i3,a)') 'get_string_element:',j, &
    '. element processed in the string'
    call print_one_member(temp)
  endif
  temp => temp%next
  j = j + 1
enddo
end subroutine get_string_element

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine insert_element(this_string,i_in,this_elem)
!------------------------------------------------------------------------
!
!   insert element into specified position in the string
!
!   on entry:  this_string - entering string
!              i_in        - position in the string
!              this_elem   - entering element
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i_in
type(dll_member), pointer :: this_elem, temp, temp1, temp2
logical :: on_head = .false.

call QEnter('insert_element')

if (get_iprint() >= 25) then
!if (get_iprint() >= 5) then
  write(LUPRI,'(/,2x,a,i3)')  &
  'insert_element: entering position in the string:',i_in
endif 

if (.not.associated(this_string)) then
  write(LUPRI,'(/,2x,a)') 'insert_element: entering "this_string" not associated !' 
  call Quit('insert_element: entering "this_string" not associated !')
else
  if (.not.associated(this_elem)) then
    write(LUPRI,'(/,2x,a)') 'insert_element: entering element "this_elem" not associated !' 
    call Quit('insert_element: entering element "this_elem" not associated !')
  endif
endif

i = i_in

if (i_in > this_string%n) then
  i = this_string%n
  if (get_iprint() >= 27) then
    write(LUPRI,'(2x,a,i3,a,i3)')  &
    'insert_element: entering position ',i_in,' changed to the end, i=',i
  endif 
else if (i_in < 0) then
  i = 1
  if (get_iprint() >= 27) then
    write(LUPRI,'(2x,a,i3,a,i3)')  &
    'insert_element: entering position ',i_in,' changed to the neginning, i=',i
  endif 
endif
 
 on_head = .false. ! pgf90 6.2.4 securing
if (i == 0) then
 on_head = .true.; i = 1
endif

! ... get the pointer to the i-th element of the string
call get_string_element(this_string,i,temp)

if (get_iprint() >= 27) then
  write(LUPRI,'(/,2x,a,i3,a)') '*** insert_element: ',i,'. member cell'
  call print_one_member(temp)
endif

!----------------------------------------------------
! ... ALLOCATES NEW MEMORY POSITION !
!----------------------------------------------------
allocate(temp1, stat = ierr)
if (ierr /= 0) then
  write(LUPRI,'(/,2x,a)') 'insert_element: error in allocation!'
  call Quit('insert_element: error in allocation !')
endif
!temp1 => this_elem

!  ... store entering "this_element" into new fresh memory
temp1 = this_elem 

if (on_head) then
!   ... insert element at the very beginning of the string ...
! -------------------------------------------------------------   
  temp1%next => this_string%head
  temp1%previous => null()

  this_string%head%previous => temp1
  this_string%head => temp1

else
!        ... insert element AFTER the i-th element "temp" 
! -------------------------------------------------------------   
   temp2 => temp%next
! ... 
   temp%next => temp1
   temp1%previous => temp
! ... 
   if (associated(temp2)) then
     temp2%previous => temp1
     temp1%next => temp2
   else
   ! ... temp1 is the last member
     this_string%rear => temp1
   endif
endif

! ... increase counter for number of elements in the string 
this_string%n  =  this_string%n  +  1

if (get_iprint() >= 25) then
!if (get_iprint() >= 5) then
 write(LUPRI,'(/,2x,a)')   &
 '*** insert_element> print the newly updated string:'
 call print_string_elem(this_string)
endif

call QExit('insert_element')
end subroutine insert_element

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine check_string(this_string,i)
!------------------------------------------------------------------------
!
!    checks the string and i-th position of the element
!
!  called when  pick up of i-th element from the string wanted
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
call QEnter('check_string')

if (.not.associated(this_string)) then
 write(LUPRI,'(2x,a)') 'check_string: entering string NOT associated (allocated) !'
 call Quit('check_string: entering string NOT associated (allocated) !')
else
 if (.not.associated(this_string%head)) then
  write(LUPRI,'(2x,a)') 'check_string: entering string is empty !'
  call Quit('check_string: entering string is empty (no member elements) !')
 endif
 if (i<0 .or. i>this_string%n) then
  write(LUPRI,'(2x,a,i3,a,i3)')  'check_string: number of elements in the string:',  &
    this_string%n,' while position to be manipulated:',i
  call Quit('check_string: wrong element position in the string to be manipulated !')
 endif
endif

call QExit('check_string')
end subroutine check_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine find_splitting(this, itype, ipos)
!------------------------------------------------------------------------
!
!  find type and position in the entering string for splitting
!
!  itype , ipos
!
!   itype = 1
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
integer, intent(inout) ::  itype, ipos
type(dll_member), pointer :: temp, temp1
call QEnter('find_splitting')

if (.not.associated(this)) then
 write(LUPRI,'(2x,a)') ': this string is NOT allocated !'
 call Quit(' this string is NOT allocated !')
else
 if (.not.associated(this%head)) then
  write(LUPRI,'(2x,a)') 'find_splitting: this string is empty !'
  call Quit('find_splitting: this string is empty !')
 endif
endif

if (get_iprint() >= 25) then
  write(LUPRI,'(/,2x,a)')  &
  'in find_splitting...searching for splitting position of the entering string'
  call FLSHFO(LUPRI)
endif

if (this%n >= 2) then

  temp => this%head
  !temp1 => temp%next
  if (associated(temp)) temp1 => temp%next
  ipos = 0
  i = 0
  do while ( associated(temp).and.associated(temp1) ) 
    i = i + 1
! ... seeking for a+ i sequence
    if (associated(temp%member%pSQOper).and.associated(temp1%member%pSQOper)) then

      if  (.not.( temp%member%pSQOper%is_creation_oper).and.  &
                (temp1%member%pSQOper%is_creation_oper).and. &
                  temp%member%pSQOper%p /=  temp1%member%pSQOper%p )  then
                    itype = 1; ipos = i
          if ( get_iprint() >= 27 ) then
            write(LUPRI,'(2x,a,i3,a)') 'catched splitting position of p q+, at p:',ipos, &
            ', leaving cycle...'
          endif
         Exit
      endif

    endif
   ! temp => temp%next
     if (associated(temp)) temp => temp%next
   ! temp1 => temp1%next
    if (associated(temp1)) temp1 => temp1%next
  enddo

  if (ipos > 0) then
    if ( get_iprint() >= 25 ) then
      write(LUPRI,'(/,2x,a)')   &
      'string members for splitting: P1 p q+ P2 = P1 d_pq P2 - P1 q+ p P2'
      call print_one_member(temp); call print_one_member(temp1)
    endif
  else
    if ( get_iprint() >= 25 ) then
      write(LUPRI,'(/,2x,a)')   &
      'find_splitting: nothing found for splitting the string !'
    endif
  endif

endif

call QExit('find_splitting')
end subroutine find_splitting

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_head_cells(hc1,hc2) result(equal)
!------------------------------------------------------------------------------
!
!     Check whether two "head_cell" pointers "hc1,hc2" are pointing 
! to the same content - list of strings
!
!------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(head_cell), pointer :: hc1, hc2
logical :: equal
call QEnter('equal_head_cells')

equal = .false.;

if (associated(hc1).and.associated(hc2)) then

  if (associated(hc1%p_string).and. associated(hc2%p_string)) then
     equal = equal_strings(hc1%p_string,hc2%p_string)
  endif

endif

call QExit('equal_head_cells')
end function equal_head_cells

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_strings(str1,str2) result (equal)
!------------------------------------------------------------------------------
!
! Check whether two-strings are equal.  Verify element after element.
!
! Called from: equal_head_cells
!
!------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: str1,str2
type(dll_member), pointer :: temp1, temp2
logical :: equal, equal_elem, all_equal_elem
call QEnter('equal_strings')

equal = .false.
if (associated(str1) .and. associated(str2)) then

 if ( (str1%factor == str2%factor) .and. (str1%n == str2%n)) then
   temp1 => str1%head; temp2 => str2%head
   all_equal_elem = .true.
   do while (associated(temp1).and.associated(temp2))
    equal_elem = equal_elements(temp1%member,temp2%member)
    all_equal_elem = all_equal_elem .and. equal_elem
    temp1 => temp1%next; temp2 => temp2%next
   enddo
  
 endif

endif

call QExit('equal_strings')
end function equal_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function equal_elements(elem1,elem2) result (equal)
!------------------------------------------------------------------------------------
!
!  Compare two entering string members "elem1,elem2". Returns .true. if these
! are equal.
!
!------------------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(string_member), intent(in) :: elem1,elem2
logical :: equal
call QEnter('equal_elements')
equal = .false.

if (associated(elem1%pSQOper).and.associated(elem2%pSQOper)) then
   equal = equal_SQ_operators(elem1%pSQOper,elem2%pSQOper)
else if (associated(elem1%pKronDelta).and.associated(elem2%pKronDelta)) then
   equal = equal_KD(elem1%pKronDelta,elem2%pKronDelta)
else if (associated(elem1%pIndxTerm).and.associated(elem2%pIndxTerm)) then
   equal = equal_IndxTerms(elem1%pIndxTerm,elem2%pIndxTerm)
else
   write(LUPRI,'(2x,a)') 'equal_elements: dead branch !' 
   call Quit('equal_elements: dead branch !')
endif

call QExit('equal_elements')
end function equal_elements

!#endif
End Module Operators_string


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module Evaluate_strings
!#if !defined (VAR_IFORT)
!------------------------------------------------------------------------
!
!     Module for evaluating expect. value of all strings
!
!   <0| strings |0> where |0> is one Slater determinant
!
!------------------------------------------------------------------------

 use Operators_string
 use control_print

  logical, private :: evaluate = .false.

 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_evaluate(log_value)
!------------------------------------------------------------------------
! set the value of "evaluate"
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
logical :: log_value
evaluate = log_value
if (get_iprint() >= 8) then
  write(LUPRI,'(2x,a,l1)') 'set_evaluate: "evaluate" set to:',evaluate 
endif
end subroutine set_evaluate

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function do_evaluate() result(eval)
!------------------------------------------------------------------------
!           returns the proctected value of "evaluate"
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
logical :: eval
  eval = evaluate
end function do_evaluate

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_exp_val_strings(heads)
!------------------------------------------------------------------------
!
!  Process all strings for the expectation value. 
!
!  Called when "evaluate" is .TRUE.
!
!
!
!
!
!
!------------------------------------------------------------------------
!  use Operators_string
#include "implicit.h"
#include "priunit.h"
type(heads_list), intent(inout) :: heads
type(head_cell), pointer :: temp, temp_prev
logical :: do_process_strings = .true., do_process_strings_prev = .false.
call QEnter('process_exp_val_strings')

if (get_iprint() >=7) then
  call Header('*** Output from "process_exp_val_strings" ***',-1)
endif

i=1
do while (do_process_strings.and..not.do_process_strings_prev)
 do_process_strings_prev = do_process_strings
 write(LUPRI,*) 'in evaluate_string_exp_val cycle'
 call evaluate_string_exp_val(heads, do_process_strings)
 i = i + 1
enddo

!temp => heads%main_head; i=1
!do while  (associated(temp))
!  temp_prev => temp 
!  call evaluate_string_exp_val(temp%p_string)
!  temp => temp%next 
!  i = i + 1
!enddo

call quit('in process_exp_val_strings')

call QExit('process_exp_val_strings')
end subroutine process_exp_val_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine evaluate_string_exp_val(heads,process)
!------------------------------------------------------------------------
!                      string must be rearranged 
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(heads_list), intent(inout) :: heads
logical, intent(inout) :: process
call QEnter('evaluate_string_exp_val')

if (get_iprint()>=8) then
   write(LUPRI,'(/,2x,a)') 'evaluate_string_exp_val'
endif

!
! ... search for mixed Kron.deltas, test. presence of sec.quant.op...
!

! VERY COMPLICATED !!!

process = .false.
call QExit('evaluate_string_exp_val')
end subroutine evaluate_string_exp_val

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine change_orbital_index(heads,i_this,i_new,is_occup)
!------------------------------------------------------------------------
!
!   Change the index of general (must be) orbital (spinor) throughout all strings
! 
! 
!                     i_this  -------->   i_new
!                   (general)        (occupied or virtual - when changing
!                                                       creat/anihil. operator)
!
!
!      For instance:    d_ip ----> d_ii
!
!                       X_jq ----> X_ji 
!
!                    p+ q+ jl ---> p+ i+ jl
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i_this,i_new
logical, intent(in) :: is_occup
type(heads_list), intent(inout) :: heads
type(head_cell), pointer :: temp

! ... perhaps check of existence of general i_this in the heads list

! ... go string after string, and then element after element and replace
temp => heads%main_head
do while (associated(temp))

  call change_orbital_index_string(temp%p_string,i_this,i_new,is_occup)

  temp => temp%next ! ... next string in the list

enddo
end subroutine change_orbital_index

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine change_orbital_index_string(this,i_this,i_new,is_occup)
!------------------------------------------------------------------------
!  change orbital(spinor) index in the string of elements;
!  called from "change_orbital_index"
!------------------------------------------------------------------------
#include "priunit.h"
type(one_string), pointer :: this
integer, intent(in) :: i_this,i_new
logical, intent(in) :: is_occup
type(dll_member), pointer :: temp
if (.not.associated(this).or..not.associated(this%head)) then
 call Quit('change_orbital_index_string: string pointer/head NOT associated !')
endif
temp => this%head
do while (associated(temp))
  if      (associated(temp%member%pSQOper)) then
     !   call change_orbital_index_sqo(temp%member%pSQOper,i_this,i_new,is_occup)     
  else if (associated(temp%member%pKronDelta)) then
     !   call change_orbital_index_kd(temp%member%pKronDelta,i_this,i_new,is_occup)     
  else if (associated(temp%member%pIndxTerm)) then
     !   call change_orbital_index_it(temp%member%pIndxTerm,i_this,i_new,is_occup)     
  else
     call Quit('change_orbital_index_string: dead branch !')
  endif
  temp => temp%next
enddo
end subroutine change_orbital_index_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Xevaluate_string_exp_val(this)
!------------------------------------------------------------------------
!
!   evaluate <0| operators_string |0>,
! where |0> is Slater determinant ("Fermi vacuum")
!
!------------------------------------------------------------------------
  use Operators_string
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
logical :: rearranged

if (get_iprint() >= 8) then
  write(LUPRI,'(/,2x,a)') 'evaluate_string_exp_val:'
endif

if (.not.associated(this)) then
  write(LUPRI,'(2x,a)') 'evaluate_string_exp_val: string NOT allocated !'
  call Quit('evaluate_string_exp_val: string NOT allocated !')
endif

if (get_iprint() >= 8) then
  write(LUPRI,'(/,2x,a)') 'this string:'
  call print_string_elem(this)
endif

rearranged = is_rearranged_string(this)
if (rearranged) then
  if (get_iprint() >= 8) then
    write(LUPRI,'(8x,a)') '...this string is rearranged...'
  endif
endif
end subroutine Xevaluate_string_exp_val

end module Evaluate_strings


#if !defined SECQUANT_ROUTINE
Program TEST_SEC_QUANT
#else
subroutine TEST_SEC_QUANT
#endif
 use Operators_string
 use Evaluate_strings
#include "priunit.h"
  type(KronDelta), pointer :: d_ij, d_ab
  type(SQOper), pointer :: i,j,a,b,o1,o2,o3,o4
  type(indexed_term), pointer :: Tia, h_pq, Tijab, Xib
  type(dll_member), pointer :: head => null(), k, d_ik
  type(heads_list) :: HL
  type(head_cell), pointer :: pSH
  type(one_string), pointer ::  pos

!  type(one_string), target ::  os
 ! call set_iprint(10)

  !allocate(d_ij)
  !d_ij => make_kron_delta(1,2); ! call print_kd(d_ij)

  !allocate(d_ab)
 ! d_ab => make_kron_delta(3,4); ! call print_kd(d_ab)

   Tia => make_it(2,2,'Tia')
   !call fill_indexes(Tia,2,3)
   call fill_indexes(Tia,4,3)
   call print_it(Tia)  ; !call quit('here')

   h_pq => make_it(1,2,'h_pq')
   !call fill_indexes(Tia,2,3)
   call fill_indexes(h_pq,1,2)
   call print_it(h_pq)  ; !call quit('here')

  !allocate(i)
  !allocate(o3)
  !i = make_sqo(1,.true.); ! call print_sqo(i)
  !o3 = make_sqo(3); ! call print_sqo(i)

  o4 => make_sqo(4,.false.,.true.)
  o3 => make_sqo(3,.true.,.false.); ! call print_sqo(i)

  o2 => make_sqo(2,.false.); ! call print_sqo(i)
 ! o2 => make_sqo(2,.true.); ! call print_sqo(i)

  !o1 = make_sqo(1,.true.);
  o1 => make_sqo(1,.true.);

!  call print_sqo(a)

  !allocate(pos)
 ! pos = create_string(head,0.50D0)
  pos => create_string(head)

  allocate(pSH)

  pSH%p_string => pos 

 ! pSH%previous => null(); pSH%next => null()

  HL%main_head=>pSH
  pos%head => head 

  call Insert_IndxTerm_end(pos,h_pq)
  call Insert_SQOp_end(pos,o1)
  call Insert_SQOp_end(pos,o2)

  call Insert_IndxTerm_end(pos,Tia)
  call Insert_SQOp_end(pos,o3)
  call Insert_SQOp_end(pos,o4)

 ! call Insert_KD_end(pos,d_ij)
 ! call Insert_KD_end(pos,d_ab)

 ! write(LUPRI,'(/,2x,a)') '*** printing full string:'
 ! call print_string(pos)
 ! call quit('end here...')

  !call print_all_strings(HL)

!  call rearrange_string(pos)
 ! call print_string(pos)

  call set_iprint(10)

 ! call set_upper_occ(2)
  call set_evaluate(.true.)
  call process_strings(HL)

 ! if (do_evaluate()) call process_exp_val_strings(HL)

  call cancel_strings(HL)

 ! call exchange_string_elements(pos,1,3)
 ! call print_string(pos)

 ! call exchange_string_elements(pos,5,3)
 ! call print_string(pos)

 ! call set_iprint(10)

 ! call rearrange_string(pos)
 ! call print_string(pos)

 ! call remove_element(pos,5)
 ! call remove_element(pos,4)
 ! call remove_element(pos,3)
 ! call print_string(pos)

 ! call rearrange_string(pos)
  !call move_element_beg(pos,3)
  !call print_string(pos)

  !call rearrange_string(pos)
  !call print_string(pos)

 ! call move_element_beg(pos,2)
 ! call move_element_beg(pos,5)

!  call remove_element(pos,2)
!  call remove_element(pos,1)
!  call remove_element(pos,1)
!  call remove_element(pos,1)
!  call remove_element(pos,1)

 !call exchange_string_elements(pos,1,3)
 !call exchange_string_elements(pos,2,3)
 ! call exchange_string_elements(pos,1,2)
 !call print_string(pos)

 !allocate(k); allocate(d_ik)
 !k = make_sqo_cell(3); call print_one_member(k)
 !d_ik = make_kd_cell(1,3); call print_one_member(d_ik)
 !call  insert_element(pos,3,k)
 !call  insert_element(pos,1,d_ik)

 !call print_string(pos)

 !call remove_element(pos,1)
 !call remove_element(pos,1)
 !call Insert_SQOp_end(pos,i)
 !call remove_element(pos,1)
 !call remove_element(pos,1)
 !call exchange_string_elements(pos,1,2)
 !call exchange_string_elements(pos,1,2)

 !call print_string(pos)
 !call print_all_strings(HL)

 !call process_strings(HL)
! call process_strings(HL)

!#else
! write(6,*) "Deactivated for ifort compiler version lower than 10.1.015 "
!#endif

#if !defined SECQUANT_ROUTINE
End Program TEST_SEC_QUANT
#else
end subroutine TEST_SEC_QUANT
#endif

