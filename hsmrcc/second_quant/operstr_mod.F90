Module Operators_string
!
! ... string of operators ... 
! ... kronecker deltas, creation/anihilation operators, indexed terms
! ...  

use KroneckerDelta
use SecQuantOperator
use Indexed_Quantities

 ! ... one string member (creat/anihil. operator, Kronecker delta, indexed term...)
  type string_member
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
!
! ... one string of mixed members ("dll_member")
!
 real(kind=8) :: factor = 0.0D0 ! default factor ...
 type(dll_member), pointer :: head => null() ! ... head to the string of members
end type one_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type head_cell
!------------------------------------------------------------------
!
!             ... member of cell of HEADS in DLL
!
!------------------------------------------------------------------
 ! ... content - pointer to string
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
   type(head_cell), pointer :: main_head

end type heads_list

 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function create_string(head, factor) result(this)
!---------------------------------------------------------------------------
!
!  constructor, returns "one_string" cell
!
! entering "head" - pointer to the first element in string 
! entering "factor" - numerical factor attached to the string (ie +1,-1,...)
!
!---------------------------------------------------------------------------
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

 recursive subroutine Insert_SQO(this, head)
! ... adds second quantized operator to the doubly linked list of elements
 type(SQOper), intent(inout) :: this 
 type(dll_member), pointer :: head

 if (.not. associated(head)) then
    allocate(head); nullify(head%next) 
    allocate(head%member%pSQOper)
    head%member%pSQOper = make_sqo(this%p,this%is_creation_oper)
    print *,'Insert_SQO: element inserted'; call print_sqo(head%member%pSQOper)
 else
    call Insert_SQO(this,head%next)
 endif 
 end subroutine Insert_SQO

 recursive subroutine Insert_KD(this, head)
! ... inserts kron.delta into doubly linked list of mixed elements
 type(KronDelta), intent(inout) :: this 
 type(dll_member), pointer :: head

 if (.not. associated(head)) then
    allocate(head); nullify(head%next) 
    allocate(head%member%pKronDelta)
    head%member%pKronDelta = make_kron_delta(this%p,this%q)
    !head%member%pKronDelta => make_kron_delta(this%p,this%q)
    print *,'Insert_KD: element inserted'; call print_kd(head%member%pKronDelta)
 else
    call Insert_KD(this,head%next)
 endif 
 end subroutine Insert_KD
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
! prints out one string of members (operators, kron.deltas, indexed terms) ...
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
! print-out one member
  type(dll_member), pointer :: this

  if (associated(this%member%pSQOper)) then
   call print_sqo(this%member%pSQOper)
  else if (associated(this%member%pKronDelta)) then
   call print_kd(this%member%pKronDelta)
  else
   stop 'neither SQOp, nor Kron.delta !'
  endif

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
! ... cycle ...
temp => heads%main_head
do while (associated(temp))
  write(LUPRI,'(2x,a,i3)') ' string #',i
  call process_one_string(temp%p_string)
  temp => temp%next
  i = i + 1
enddo

call QExit('process_strings')
end subroutine process_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_one_string(this)
!------------------------------------------------------------------------
!
!  basic formula:
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
type(one_string), pointer :: this
type(dll_member), pointer :: temp
call QEnter('process_one_string')

if (.not.associated(this)) then
 write(LUPRI,*) 'process_one_string: this string is NOT allocated !'
 call Quit('process_one_string: this string is NOT allocated !')
endif

if (.not.associated(this%head)) then
 write(LUPRI,*) 'process_one_string: this string of elements is NOT allocated !'
 call Quit('process_one_string: this string is NOT allocated !')
endif

write(LUPRI,*) 'processing ONE allocated string...'

call print_string(this)

!... find first splitting place 

! .. allocate one more string ...


call QExit('process_one_string')
end subroutine process_one_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine rearrange_string(this)
!------------------------------------------------------------------------
!
! rearrange string that all indexed terms and kronecker deltas are 
! on the left separated from bunch of creation/anihilation operators.
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
call QEnter('rearrange_string')

write(LUPRI,*) 'rearrangin string ...'

call QExit('rearrange_string')
end subroutine rearrange_string

End Module Operators_string

