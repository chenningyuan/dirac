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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module Evaluate_strings
#if !defined (VAR_IFORT)
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
!
!                      string must be rearranged 
!
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
call QEnter('change_orbital_index')

! ... perhaps check of existence of general i_this in the heads list

! ... go string after string, and then element after element and replace
temp => heads%main_head
do while (associated(temp))

  call change_orbital_index_string(temp%p_string,i_this,i_new,is_occup)

  temp => temp%next ! ... next string in the list

enddo
call QExit('change_orbital_index')
end subroutine change_orbital_index

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine change_orbital_index_string(this,i_this,i_new,is_occup)
!------------------------------------------------------------------------
!
!  change orbital(spinor) index in the string of elements
!
!  called from "change_orbital_index"
!
!------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(one_string), pointer :: this
integer, intent(in) :: i_this,i_new
logical, intent(in) :: is_occup
type(dll_member), pointer :: temp
call QEnter('change_orbital_index_string')

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


call QExit('change_orbital_index_string')
end subroutine  change_orbital_index_string

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
call QEnter('evaluate_string_exp_val')

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

call QExit('evaluate_string_exp_val')
end subroutine Xevaluate_string_exp_val

#endif
end module Evaluate_strings
