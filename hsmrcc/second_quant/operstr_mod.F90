Module Operators_string
!
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
!
! ... one string member (creat/anihil. operator, Kronecker delta, indexed term...)
!
!---------------------------------------------------------------------------------
    type(SQOper),   pointer :: pSQOper => null()
    type(KronDelta),pointer :: pKronDelta => null()
    type(indexed_term), pointer :: pIndxTerm  => null()
  end type string_member

! ... one member with 2 pointers (to be in the doubly linked list) ...
type dll_member
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
 type(head_cell), pointer :: previous => null(), next => null()
end type head_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type heads_list
!------------------------------------------------------------------------
!
!             DLL of heads of various strings; 
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function create_string(head, factor) result(this)
!---------------------------------------------------------------------------
!
!  constructor, returns "one_string" cell
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
!type(one_string), pointer :: this
type(one_string) :: this

call QEnter('create_string')

!if (associated(head)) then

   !allocate(this)

   if (present(head)) then
     !if (associated(head))  this%head = head
     this%head => head
   endif

   if (present(factor)) then
     this%factor = factor
   else
     this%factor = 1.0D0
   endif

   ! rather secure zero number of elements
   this%n = 0

!else
!  write(*,'(2x,a)') '"one_string" cell not allocated, because ' &
!                  //'entering pointer "head" is not associated !'
!endif

call QExit ('create_string')
end function create_string

! subroutine Insert_KD(head,i,j)
! adds kronecker delta to the end of operator_string list
 !  type(string_member), pointer :: head
 !  integer(kind=4), intent(in) :: i,j
 !  type(string), pointer :: temp
! end subroutine Insert_KD

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
! this_string
! this_oper
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
   write(LUPRI,'(2x,a)') 'Some error in the empty entering string !'
   write(LUPRI,'(2x,a,l)')        &
   '...associated(this_string%rear)=',associated(this_string%rear)
   write(LUPRI,'(2x,a,i3)') 'this_string%n=',this_string%n
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
   write(LUPRI,'(2x,a,l)')        &
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
 recursive subroutine Insert_KD_rec(this, head)
!--------------------------------------------------------------------------------------
!
! ... inserts kron.delta into doubly linked list of mixed elements
!
!--------------------------------------------------------------------------------------
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
!                print all strings under list "heads"
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!type(heads_list), intent(inout) :: heads
type(heads_list), intent(in) :: heads
type(head_cell), pointer :: temp
call QEnter('print_all_strings')

write(LUPRI,*) 'in print_all_strings:'

if (associated(heads%main_head)) then

 i = 1
! ... cycle ...
 temp => heads%main_head
 do while (associated(temp))
  write(LUPRI,'(2x,a,i3)') ' string #',i
  call print_string(temp%p_string)
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
subroutine print_string(head)
!------------------------------------------------------------------------
!
!   print-out string of elements of the string pointed with 'head'
! together with the factor 
!
!------------------------------------------------------------------------
type(one_string), pointer :: head

call QEnter('print_string')
if (.not.associated(head)) then
  write(*,'(/,2x,a)') 'print_string: entering "head" not associated !' 
  call Quit('print_string: entering "head" not associated !')
else
  write(*,'(/,a,f12.5)')   & 
  '*** printing  one string of elements with the factor',head%factor
  call print_members(head%head)
endif

call QExit('print_string')
end subroutine print_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_members(head)
!------------------------------------------------------------------------
!
! prints out one string of members (operators, kron.deltas, indexed terms) ...
!
!------------------------------------------------------------------------
  !type(dll_member), pointer, intent(in) :: head
  type(dll_member), pointer :: head, temp

 write(*,'(a)') '* individual elements in the string:'

 if (.not.associated(head)) then
  write(*,'(2x,a)')   &
 "...string of members (operators, kron.deltas, indexed terms) empty"
 else
   temp => head
   do while (associated(temp))
     call print_one_member(temp)
     temp => temp%next
   enddo
 endif
end subroutine print_members

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine print_one_member(this)
! -----------------------------------------------------------------
!
!     print-out one member according to entering pointer
!
! -----------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
  type(dll_member), pointer :: this

  call QEnter('print_one_member')

  if (.not.associated(this)) then
    write(LUPRI,*) 'print_one_member: entering pointer is not associated !'
    call Quit('print_one_member: entering pointer is not associated !')
  endif

  if (associated(this%member%pSQOper)) then
   call print_sqo(this%member%pSQOper)
  else if (associated(this%member%pKronDelta)) then
   call print_kd(this%member%pKronDelta)
  else
   write(LUPRI,'(2x,a)') '*** print_one_member dead branch output ***'
   write(LUPRI,'(2x,a,l)')   &
   'associated(this%member%pSQOper):',associated(this%member%pSQOper)
   write(LUPRI,'(2x,a,l)')   &
   'associated(this%member%pKronDelta):',associated(this%member%pKronDelta)
   write(LUPRI,'(2x,a,l)')   &
   'associated(this%member%pIndxTerm):',associated(this%member%pIndxTerm)

   call Quit('print_one_member: neither SQOp, nor Kron.delta !')
  endif

 call QExit('print_one_member')

 end subroutine print_one_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_strings(heads)
!------------------------------------------------------------------------
!
!   most important routine - process strings of elements
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(heads_list), intent(inout) :: heads
type(head_cell), pointer :: temp
call QEnter('process_strings')

write(LUPRI,*) 
write(LUPRI,*) 'in process_strings...'

!call print_all_strings(heads)

! ... go string after string ...

if (.NOT.associated(heads%main_head)) then
  write(LUPRI,*) 'heads%main_head NOT associated !'
  call Quit('heads%main_head NOT associated !')
endif

i = 1
! ... cycle ... string after string ...
temp => heads%main_head
do while (associated(temp))
  write(LUPRI,'(2x,a,i3)') ' string #',i
  !call process_one_string(temp%p_string)
  call process_one_string(temp)
  temp => temp%next
  i = i + 1
enddo

call QExit('process_strings')
end subroutine process_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_one_string(this)
!------------------------------------------------------------------------
!
! Process entering string : apply commutator relation
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
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!type(one_string), pointer :: this
type(head_cell), pointer :: this
type(dll_member), pointer :: temp
call QEnter('process_one_string')

if (.not.associated(this)) then
 write(LUPRI,*) 'process_one_string: this string is NOT allocated !'
 call Quit('process_one_string: this string is NOT allocated !')
endif

if (.not.associated(this%p_string)) then
 write(LUPRI,*) 'process_one_string: this string of elements is NOT allocated !'
 call Quit('process_one_string: this string is NOT allocated !')
endif

write(LUPRI,*) 'processing ONE allocated string...'

call print_string(this%p_string)

!... find first splitting place 
call find_splitting(this%p_string)

! .. allocate one more string ...


call QExit('process_one_string')
end subroutine process_one_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine rearrange_string(this)
!------------------------------------------------------------------------
!
!  rearrange string that all indexed terms and kronecker deltas are 
! on the left separated from the rest -  bunch 
! of creation/anihilation operators.
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
type(dll_member), pointer :: temp, temp1
call QEnter('rearrange_string')

write(LUPRI,'(/,2x,a)') 'rearranging the string...'

if (.not.associated(this)) then
 write(LUPRI,'(2x,a)') 'rearrange_string: entering string NOT associated (allocated) !'
 call Quit('rearrange_string: entering string NOT associated (allocated) !')
else
 if (.not.associated(this%head)) then
  write(LUPRI,'(2x,a)') 'rearrange_string: entering string is empty !'
  call Quit('rearrange_string: entering string is empty (no member elements) !')
 endif
endif

! ... go element after element
ipos = 0


call QExit('rearrange_string')
end subroutine rearrange_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine move_element_beg(this_string,i)
!------------------------------------------------------------------------
!
!   move element in string from  position "i" to the very beginning
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
type(dll_member), pointer :: temp, temp1, this_element
call QEnter('move_element_beg')

if (.not.associated(this_string)) then
 write(LUPRI,'(2x,a)') 'move_element_beg: entering string NOT associated (allocated) !'
 call Quit('move_element_beg: entering string NOT associated (allocated) !')
else
 if (.not.associated(this_string%head)) then
  write(LUPRI,'(2x,a)') 'move_element_beg: entering string is empty !'
  call Quit('move_element_beg: entering string is empty (no member elements) !')
 endif
 if (i<0 .or. i>this_string%n) then
  write(LUPRI,'(2x,a,i3,a,i3)')  'number of elements in the string:',  &
   this_string%n, 'while position to be moved at the beginning:',i
  call Quit('move_element_beg: wrong position to be moved ')
 endif
endif

  if (get_iprint() >= 5) then
     write(LUPRI,'(2x,a)') 'in "move_element_beg:" '
  endif

j = 1

 temp => this_string%head
 do while (j < i)
   temp => temp%next
   j = j + 1
   if (get_iprint() >= 5) then
     write(LUPRI,'(2x,i3,a)') j,'. element in the string'
     call print_one_member(temp)
   endif
 enddo

call QExit('move_element_beg')
end subroutine move_element_beg

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine find_splitting(this)
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this

if (.not.associated(this)) then
 write(LUPRI,*) ': this string is NOT allocated !'
 call Quit(' this string is NOT allocated !')
endif

write(LUPRI,*) 'in find_splitting...'

end subroutine find_splitting

End Module Operators_string

