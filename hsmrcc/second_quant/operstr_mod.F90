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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_sqo_cell(i,is_creat) result (this)
!----------------------------------------------------------------------------
!
!                   Returns second quantized operator cell
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i
logical, optional, intent(in) :: is_creat
type(dll_member) :: this
type(SQOper), target :: this_sqo
call QEnter('make_sqo_cell')

 this_sqo = make_sqo(i,is_creat)

 this%member%pSQOper =>  this_sqo

call QExit('make_sqo_cell')
end function make_sqo_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_kd_cell(i,j) result (this)
!----------------------------------------------------------------------------
!
!               Returns cell "this" with the Kronecker delta
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i,j
!type(dll_member) :: this
type(dll_member), target :: this
type(dll_member), pointer :: temp
type(KronDelta), target :: this_kd
call QEnter('make_kd_cell')

! ... true Kronecker delta ...
this_kd = make_kron_delta(i,j)

this%member%pKronDelta => this_kd

! ... check  for  older compilers, like pgf90 5.2-4 (opteron64, Tel Aviv)
if (associated(this%member%pSQOper)) then
  if (get_iprint() >= 55) then
    write(LUPRI,'(2x,a)') 'make_kd_cell: pSQOper randomly associated !'
  endif
  this%member%pSQOper  => null()
endif

if (associated(this%member%pIndxTerm)) then
  if (get_iprint() >= 55) then
    write(LUPRI,'(2x,a)') 'make_kd_cell: pIndxTerm randomly associated !'
  endif
  this%member%pIndxTerm   => null()
endif

!allocate(temp)
temp => this
 
if ( get_iprint() >= 35 ) then
  write(LUPRI,'(2x,a,$)') 'make_kd_cell: created Kronecker delta:'
  !call print_one_member(this)
  call print_one_member(temp)
endif

call QExit('make_kd_cell')
end function make_kd_cell

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

   ! rather secure zero number of elements and the tail
   this%n = 0
   this%rear => null()  
!else
!  write(*,'(2x,a)') '"one_string" cell not allocated, because ' &
!                  //'entering pointer "head" is not associated !'
!endif

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
   write(LUPRI,'(2x,a,l)')        &
   '...associated(this_string%head)=',associated(this_string%head)
   write(LUPRI,'(2x,a,l)')        &
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
subroutine print_string(p_this_string)
!------------------------------------------------------------------------
!
!   print-out string of elements of the string pointed with "p_this_string"
! together with the factor 
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: p_this_string
call QEnter('print_string')

if (.not.associated(p_this_string)) then
  write(LUPRI,'(/,2x,a)') 'print_string: entering "p_this_string" not associated !' 
  call Quit('print_string: entering "p_this_string" not associated !')
else
  write(LUPRI,'(/,a,i3,a,f12.5)')   & 
  '*** printing  one string of ',p_this_string%n, &
  ' elements and with the attached factor',p_this_string%factor
  call print_members(p_this_string%head)
endif

call QExit('print_string')
end subroutine print_string

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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

endif

call QExit('print_one_member')

end subroutine print_one_member

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_strings(heads)
!------------------------------------------------------------------------
!
!      most important routine - process all strings of elements
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(heads_list), intent(inout) :: heads
type(head_cell), pointer :: temp
call QEnter('process_strings')

if (get_iprint() >= 5) then
   write(LUPRI,'(/,2x,a)') 'in process_strings...'
endif

! ... go string after string ...

if (.NOT.associated(heads%main_head)) then
  write(LUPRI,*) 'heads%main_head NOT associated !'
  call Quit('heads%main_head NOT associated !')
endif

i = 1
! ... cycle ... string after string ...
temp => heads%main_head
do while (associated(temp))
  if (get_iprint() >= 5) then
    write(LUPRI,'(/,2x,a,i3)')  &
    '...going to process string #',i
  endif
  !call process_one_string(temp%p_string)
  call process_one_string(temp)
  temp => temp%next
  i = i + 1
enddo

if (get_iprint() >= 7) then
   write(LUPRI,'(/,2x,a)') '**** process_strings: printing all strings ****'
   call print_all_strings(heads)
endif

call QExit('process_strings')
end subroutine process_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_one_string(this)
!------------------------------------------------------------------------
!
! Process entering string : apply commutator relations
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
!             where P1, P2 are bunches of elements 
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
 call print_string(this%p_string)
endif

! ... rearrange string - move deltas, indexed terms left
call rearrange_string(this%p_string)

if (get_iprint() >= 5) then
 write(LUPRI,'(/,2x,a)')  &
 'process_one_string: now REARRANGED string going to be processed...'
 call print_string(this%p_string)
endif

! ... initialization ...
found_splitting = .true.
do while ( found_splitting )

  ! ... find first the splitting place 
  itype = -1; ipos = -1
  found_splitting = .false.
  call find_splitting(this%p_string, itype, ipos )

  ! ... allocate one more string and duplicate the string 
  if (ipos > 0) then

    found_splitting = .true.

    if ( get_iprint() >= 5 ) then
     write(LUPRI,'(/,2x,a)')  &
     '*** process_one_string - creating next string as splitting place was found ! '
    endif

    allocate(this2)
    this2%p_string => null()

    allocate(this_temp)

    this2%previous => this
    this2%next => this_temp

    this_temp%previous => this2

    this_temp => this%next
    this%next => this2

    call duplicate_string(this%p_string, this2%p_string )

    if (get_iprint() >= 7) then
       write(LUPRI,'(/,2x,a)')   &
       'process_one_string: duplicated string:'
       call print_string(this2%p_string)
    endif

    ! ... insert the string after the main string, chain it
    this2%next      => this%next
    this2%previous  => this
    this%next       => this2

    call process_2_strings(this%p_string, this2%p_string, ipos)

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
type(dll_member), pointer :: pd_ia, temp_i, temp_a
type(dll_member), target :: d_ia
call QEnter('process_2_strings')

if ( get_iprint() >= 5 ) then
  write(LUPRI,'(/,2x,a)')   &
  '*** process_2_strings: going to process 2 subsequent entering strings ***'
  call print_string(this_string)
  call print_string(copied_string)
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

write(LUPRI,*) 'd_ia made before: '; call print_string(this_string)
! ... create kronecker delta ...
!d_ia => null()
!allocate(d_ia)
d_ia = make_kd_cell(temp_i%member%pSQOper%p, temp_a%member%pSQOper%p)
pd_ia => d_ia
!call print_one_member(d_ia)
call print_one_member(pd_ia)
write(LUPRI,*) 'd_ia made after: '; call print_string(this_string)

! ... adapt first string ...
call remove_element(this_string, ipos)

!call insert_element(this_string, ipos, d_ia)
call insert_element(this_string, ipos, pd_ia)

call remove_element(this_string, ipos)

! ... adapt second, copied string ...
call exchange_string_elements(copied_string, ipos, ipos+1) 

! ... change the factor ...
copied_string%factor = - copied_string%factor

if ( get_iprint() >= 5 )  then
 write(LUPRI,'(/,2x,a)')  &
 '*** process_2_strings: adapted 2 subsequent strings'
 call print_string(this_string)
 call print_string(copied_string)
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
    call print_string(this_string)
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
! duplicate string - allocation of memory for the new one is inside
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

allocate(copied_string)

temp1 => copied_string%head;  temp2 => temp
temp => this_string%head;  
i = 0
do while (associated(temp))

  i = i + 1

  call insert_member(i,temp,temp1)

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
  call print_string(copied_string)
endif

!call quit('end in duplicate_string')

call QExit('duplicate_string')
end subroutine duplicate_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine insert_member(i,member_exist,member_new)
!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: i
type(dll_member), pointer :: member_exist, member_new
logical :: is_operator,  is_term 
call QEnter('insert_member')

if (get_iprint() >= 37) then
  write(LUPRI,'(/,2x,a)') 'insert_member: duplicating string elements'
endif

if (.not.associated(member_exist)) then
   write(LUPRI,'(/,2x,a)') 'insert_member: entering element does NOT exist !'
   call Quit('insert_member: entering element does NOT exist !')
endif

allocate(member_new)
member_new = member_exist

if (get_iprint() >= 37) then
   write(LUPRI,'(/,2x,a,i3,a,$)') 'insert_member: ',i,'. old member element:'
   call print_one_member(member_exist)
   write(LUPRI,'(2x,a,$)') '           the same new member element:'
   call print_one_member(member_new)
endif

!call quit('insert_member')
call QExit('insert_member')
end subroutine insert_member


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
         call print_string(this) 
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
!             i - position
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
!    ... move only the element which is NOT the first
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
 call print_string(this_string)
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
 call print_string(this_string)
endif

call QExit('move_element_beg')
end subroutine move_element_beg

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine remove_element(this_string,i)
!------------------------------------------------------------------------
!
!  remove i-th element from the entering string
!
!  On input: this_string - pointer to the string
!             i - position
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
type(dll_member), pointer :: temp, temp1, temp2
call QEnter('remove_element')

call get_string_element(this_string,i,temp)

! ... previous element
temp1 => temp%previous
! ... next element
temp2 => temp%next

!if (get_iprint() >= 25) then
if (get_iprint() >= 5) then
  write(LUPRI,'(/2x,a,i3)')  &
  '"remove_element": to be removed element no:',i
  call print_one_member(temp)
  write(LUPRI,'(2x,a,$)') '... the element before: '
  call print_one_member(temp1)
  write(LUPRI,'(2x,a,$)') '... the element after: '
  call print_one_member(temp2)
endif

deallocate(temp)
! ... decreasenumber of elements ...
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

!if (get_iprint() >= 25) then
if (get_iprint() >= 5) then
  write(LUPRI,'(/,2x,a)')  &
  'remove_element: printed updated string with removed element:'
  call print_string(this_string)
endif

call QExit('remove_element')
end subroutine remove_element

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine exchange_string_elements(this_string,i_in,j_in)
!------------------------------------------------------------------------
!
!          Exchange 2 elements in the string, i-th and j-th.
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
   call print_string(this_string)
endif

call QExit('exchange_string_elements')
end subroutine exchange_string_elements

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_string_element(this_string,i,temp)
!------------------------------------------------------------------------
!
!        returns pointer to the i-th element in the string
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i
type(dll_member), pointer :: temp
call QEnter('get_string_element')

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

call QExit('get_string_element')
end subroutine get_string_element

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine insert_element(this_string,i_in,this_elem)
!------------------------------------------------------------------------
!
!   insert element into the string
!
!   on entry:  this_string
!              i_in
!              this_elem
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this_string
integer, intent(in) :: i_in
type(dll_member), pointer :: this_elem, temp, temp1, temp2
logical :: on_head = .false.

call QEnter('insert_element')

!if (get_iprint() >= 25) then
if (get_iprint() >= 5) then
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

if (i == 0) then
 on_head = .true.; i = 1
endif

! ... get the pointer to the i-th element of the string
call get_string_element(this_string,i,temp)

if (get_iprint() >= 27) then
  write(LUPRI,'(/,2x,a,i3,a)') '*** insert_element: ',i,'. member cell'
  call print_one_member(temp)
endif

allocate(temp1)
!temp1 => this_elem
temp1 = this_elem

if (on_head) then
! ... insert element at the very beginning of the string ...
  temp1%next => this_string%head
  temp1%previous => null()

  this_string%head%previous => temp1
  this_string%head => temp1

else
! ... insert element AFTER the i-th element "temp" 
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

!if (get_iprint() >= 25) then
if (get_iprint() >= 5) then
 write(LUPRI,'(/,2x,a)')   &
 '*** insert_element> print the newly updated string:'
 call print_string(this_string)
endif

call QExit('insert_element')
end subroutine insert_element

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine check_string(this_string,i)
!------------------------------------------------------------------------
!
!    checks the string and i-th position of the element
!
!  call when we want pick up i-th element from the string
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
!   find type and position in the entering string for splitting
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

End Module Operators_string

