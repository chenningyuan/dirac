!
!
!
Program TEST_SEC_QUANT

  use Operators_string

  use control_print

#include "implicit.h"
#include "priunit.h"

  type(KronDelta), pointer :: d_ij, d_ab

  type(SQOper), pointer :: i,j,a,b,o1,o2,o3

  type(dll_member), pointer :: head => null(), k, d_ik

  type(heads_list) :: HL

  type(head_cell), pointer :: pSH

  type(one_string), pointer ::  pos

!  type(one_string), target ::  os

 ! call set_iprint(10)

  allocate(d_ij)
  d_ij = make_kron_delta(1,2); ! call print_kd(d_ij)

  allocate(d_ab)
  d_ab = make_kron_delta(3,4); ! call print_kd(d_ab)

  !allocate(i)
  allocate(o3)
  !i = make_sqo(1,.true.); ! call print_sqo(i)
  o3 = make_sqo(3); ! call print_sqo(i)

  !allocate(j)
  allocate(o2)
  o2 = make_sqo(2); ! call print_sqo(i)

  allocate(o1)
  o1 = make_sqo(1,.true.);
!  call print_sqo(a)

  allocate(pos)
 ! pos = create_string(head,0.50D0)
  pos = create_string(head)

  allocate(pSH)

  pSH%p_string => pos 
 ! pSH%previous => null(); pSH%next => null()

  HL%main_head=>pSH

  pos%head => head 

 ! call Insert_SQOp_end(pos,b)
  call Insert_SQOp_end(pos,o3)
  call Insert_SQOp_end(pos,o2)
  call Insert_SQOp_end(pos,o1)
 ! call Insert_KD_end(pos,d_ij)
 ! call Insert_KD_end(pos,d_ab)

 ! write(LUPRI,'(/,2x,a)') '*** printing full string:'
 ! call print_string(pos)
 ! call quit('end here...')

  !call print_all_strings(HL)

!  call rearrange_string(pos)
 ! call print_string(pos)

  call set_iprint(10)
 call process_strings(HL)


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

End Program TEST_SEC_QUANT
