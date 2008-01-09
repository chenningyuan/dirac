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
 ! ...  
 type(head_cell), pointer :: previous => null(), next => null()
end type head_cell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
type heads_list
!
!   DLL of heads of various strings; 
!
   type(head_cell), pointer :: main_head

end type heads_list

 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function create_string(head, factor) result(this)
!------------------------------------------------------------
!
!    allocates "one_string" cell, returns it
!
! entering "head" - pointer to the first element in string 
! entering "factor" - numerical factor attached to the string (ie +1,-1,...)
!
!------------------------------------------------------------
type(dll_member), pointer, optional :: head
real(kind=8), intent(in), optional :: factor 
type(one_string), pointer :: this
!type(one_string) :: this

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
subroutine print_string(head)
!------------------------------------------------------------------------
!
!   print-out string of elements of the string pointed with 'head'
! together with the factor 
!
!------------------------------------------------------------------------
type(one_string), pointer :: head

if (.not.associated(head)) then
  write(*,'(/,2x,a)') 'print_string: entering "head" not associated !' 
  call abort
else
  write(*,'(/,a,f12.5)')   & 
  '*** printing out one string of elements with the factor',head%factor
  call print_members(head%head)
endif

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
  !call abort
 else
   temp => head
   do while (associated(temp))
     call print_one_member(temp)
     temp => temp%next
   enddo
 endif
 end subroutine print_members

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine process_strings(heads)
!------------------------------------------------------------------------
!
!        most important routine - process strings of elements
!
!------------------------------------------------------------------------
type(heads_list), intent(inout) :: heads




end subroutine process_strings

End Module Operators_string

