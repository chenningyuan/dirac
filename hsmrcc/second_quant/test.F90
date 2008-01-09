!
!
!
Program TEST_SEC_QUANT

  use Operators_string
  type(dll_member), pointer :: head => null()
  
  type(KronDelta) :: d_ij, d_ab

  type(SQOper) :: i,j,a,b

  type(heads_list) :: HL

  type(head_cell), pointer :: pSH

  type(one_string), pointer ::  pos

!  type(one_string), target ::  os

 ! d_ij = make_kron_delta(1,2);  call print_kd(d_ij)
 ! d_ab = make_kron_delta(3,4);  call print_kd(d_ab)

  i = make_sqo(1); ! call print_sqo(i)

!  j=make_sqo(2,.true.)

  a = make_sqo(3,.true.);
!  call print_sqo(a)

!  b=make_sqo(4)

!  os = create_string(head,0.50D0)

  allocate(pos)
  pos = create_string(head,0.50D0)

  allocate(pSH)

!  SH%p_string=>os 

   pSH%p_string => pos 

   HL%main_head=>pSH

   call Insert_SQO(i,head)
 pos%head => head 

   call Insert_SQO(a,head)

 !  call Insert_KD(d_ij,head)
 !  call Insert_SQO(b,head)

! call print_string(pos)
  call print_all_strings(HL)

 call process_strings(HL)





End Program TEST_SEC_QUANT

