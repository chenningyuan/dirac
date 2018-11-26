module intermediates

!Module to make intermediates needed in CC and Lambda equations

     implicit none
     private

     complex(8), parameter :: MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0), &
                              ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0)

     interface get_common_blocks
        module procedure get_common_blocks_exatns
        module procedure get_common_blocks_talsh
     end interface get_common_blocks

     interface get_lambda_intermediates
        module procedure get_lambda_intermediates_exatns
        module procedure get_lambda_intermediates_talsh
     end interface get_lambda_intermediates

     public get_common_blocks
     public get_lambda_intermediates

     contains

       subroutine get_common_blocks_exatns(comm_t,int_t)

!       Get common blocks to be used in lambda intermediates

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t
        type(tens_rcrsv_t) :: one_comm_tensor

        integer :: ierr

!        ierr=exatns_tensor_create(one_comm_tensor,"one_comm_tensor",EXA_DATA_KIND_C8)
!        ierr=exatns_tensor_init(one_comm_tensor,ONE)
!
!!       tau
!        ierr=exatns_tensor_init(comm_t%tau)
!        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,b,i,j)*K()",comm_t%tau,comm_t%t2,one_comm_tensor)
!        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,i)*R(b,j)",comm_t%tau,comm_t%t1,comm_t%t1)
!        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,j)*R(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,MINUS_ONE)
!
!!       vv
!        ierr=exatns_tensor_init(comm_t%gvv)
!        ierr=exatns_tensor_contract("G(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%gvv,int_t%oovv,comm_t%tau,ONE_HALF)
!        ierr=exatns_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,int_t%vovv,comm_t%t1)
!        ierr=exatns_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,MINUS_ONE)
!        ierr=exatns_tensor_contract("G(a,c)+=F(a,c)*K()",comm_t%gvv,comm_t%fvv,one_comm_tensor)
!
!!       ov
!        ierr=exatns_tensor_init(comm_t%hov)
!        ierr=exatns_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,int_t%oovv,comm_t%t1)
!        ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",comm_t%hov,comm_t%fov,one_comm_tensor)
!
!!       oo
!        ierr=exatns_tensor_init(comm_t%goo)
!        ierr=exatns_tensor_contract("G(k,i)+=V(k,l,c,d)*T(c,d,i,l)",comm_t%gvv,int_t%oovv,comm_t%tau,ONE_HALF)
!        ierr=exatns_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%gvv,int_t%ooov,comm_t%t1)
!        ierr=exatns_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%gvv,comm_t%fov,comm_t%t1)
!        ierr=exatns_tensor_contract("G(k,i)+=F(k,i)*K()",comm_t%goo,comm_t%foo,one_comm_tensor)
!
!!       vvvv
!        ierr=exatns_tensor_init(comm_t%b_int)
!!       Permute a and b, so term appears twice
!        ierr=exatns_tensor_contract("B(a,b,c,d)+=V(k,a,c,d)*T(b,k)",comm_t%b_int,int_t%vovv,comm_t%t1)
!        ierr=exatns_tensor_contract("B(a,b,c,d)+=V(b,k,c,d)*T(a,k)",comm_t%b_int,int_t%vovv,comm_t%t1)
!        ierr=exatns_tensor_contract("B(a,b,c,d)+=V(a,b,c,d)*K()",comm_t%b_int,int_t%vvvv,one_comm_tensor)
!
!!       vovo
!        ierr=exatns_tensor_init(comm_t%h_int)
!        ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1)                 ! o(2)v(3)
!        ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,MINUS_ONE) ! o(3)v(2)
!        ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,i)*K()",comm_t%h_int,int_t%vovo,one_comm_tensor)
!
!!       oooo
!        ierr=exatns_tensor_init(comm_t%a_int)
!        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,j)*K()",comm_t%a_int,int_t%oooo,one_comm_tensor)
!        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
!        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1,MINUS_ONE) ! o(4)v(1)
!        ierr=exatns_tensor_init(comm_t%d_int)
!        ierr=exatns_tensor_contract("D(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%d_int,int_t%oovv,comm_t%tau,ONE_QUARTER_C)
!
!        ierr=exatns_tensor_destroy(one_comm_tensor)

       end subroutine get_common_blocks_exatns

       subroutine get_common_blocks_talsh(comm_t,int_t)

!       Get common blocks to be used in lambda intermediates

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens
!        use talsh_ao_to_mo, only : tensor_norm2

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t

        integer :: ierr

!       tau
        ierr=talsh_tensor_init(comm_t%tau)
        ierr=talsh_tensor_add("T(a,b,i,j)+=R(a,b,i,j)",comm_t%tau,comm_t%t2)
        ierr=talsh_tensor_contract("T(a,b,i,j)+=R(a,i)*R(b,j)",comm_t%tau,comm_t%t1,comm_t%t1)
        ierr=talsh_tensor_contract("T(a,b,i,j)+=R(a,j)*R(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,scale=MINUS_ONE)

!       vv
        ierr=talsh_tensor_init(comm_t%gvv)
        ierr=talsh_tensor_contract("G(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%gvv,int_t%oovv,comm_t%tau,scale=ONE_HALF)
        ierr=talsh_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,int_t%vovv,comm_t%t1)
        ierr=talsh_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,scale=MINUS_ONE)
        ierr=talsh_tensor_add("G(a,c)+=F(a,c)",comm_t%gvv,comm_t%fvv)

!       ov
        ierr=talsh_tensor_init(comm_t%hov)
        ierr=talsh_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,int_t%oovv,comm_t%t1)
        ierr=talsh_tensor_add("H(k,c)+=F(k,c)",comm_t%hov,comm_t%fov)

!       oo
        ierr=talsh_tensor_init(comm_t%goo)
        ierr=talsh_tensor_contract("G(k,i)+=V(k,l,c,d)*T(c,d,i,l)",comm_t%goo,int_t%oovv,comm_t%tau,scale=ONE_HALF)
        ierr=talsh_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,int_t%ooov,comm_t%t1)
        ierr=talsh_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1)
        ierr=talsh_tensor_add("G(k,i)+=F(k,i)",comm_t%goo,comm_t%foo)

!       vvvv
        ierr=talsh_tensor_init(comm_t%b_int)
!       Permute a and b, so term appears twice
        ierr=talsh_tensor_contract("B(a,b,c,d)+=V(k,a,c,d)*T(b,k)",comm_t%b_int,int_t%vovv,comm_t%t1)
        ierr=talsh_tensor_contract("B(a,b,c,d)+=V(b,k,c,d)*T(a,k)",comm_t%b_int,int_t%vovv,comm_t%t1)
        ierr=talsh_tensor_add("B(a,b,c,d)+=V(a,b,c,d)",comm_t%b_int,int_t%vvvv)

!       vovo
        ierr=talsh_tensor_init(comm_t%h_int)
        ierr=talsh_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1)                 ! o(2)v(3)
        ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
        ierr=talsh_tensor_add("H(a,k,c,i)+=V(a,k,c,i)",comm_t%h_int,int_t%vovo)

!       oooo
        ierr=talsh_tensor_init(comm_t%a_int)
        ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",comm_t%a_int,int_t%oooo)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        ierr=talsh_tensor_init(comm_t%d_int)
        ierr=talsh_tensor_contract("D(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%d_int,int_t%oovv,comm_t%tau, &
                                    scale=ONE_QUARTER_C)

       end subroutine get_common_blocks_talsh

       subroutine get_lambda_intermediates_exatns(comm_t,int_t)

!       Get fixed intermediates of lambda equations

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens,exatns_comm_tens

        type(exatns_comm_tens),intent(inout) :: comm_t
        type(exatns_intg_tens),intent(inout) :: int_t

        integer :: ierr

!       Do stuff

       end subroutine get_lambda_intermediates_exatns

       subroutine get_lambda_intermediates_talsh(comm_t,int_t,lambda_t)

!       Get fixed intermediates of lambda equations

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens
        use talsh_ao_to_mo, only : tensor_norm2

        type(talsh_comm_tens),intent(inout)   :: comm_t
        type(talsh_intg_tens),intent(inout)   :: int_t
        type(talsh_lambda_tens),intent(inout) :: lambda_t
        type(talsh_tens_t)                    :: other1_tensor

        integer :: one_dims(1)

        integer :: ierr

!       Create and initialize scalars that are to be used as tensors in contractions
        one_dims = 1
        ierr=talsh_tensor_construct(other1_tensor,C8,one_dims(1:0),init_val=ONE)

!       Two-body intermediates


!       One-body intermediates
        ierr=talsh_tensor_contract("F(c,a)+=G+(a,c)*K()",lambda_t%fbar_vv,comm_t%gvv,other1_tensor)
        ierr=talsh_tensor_add("F(k,c)+=H(k,c)",lambda_t%fbar_ov,comm_t%hov)
        ierr=talsh_tensor_contract("F(i,k)+=G+(k,i)*K()",lambda_t%fbar_oo,comm_t%goo,other1_tensor)

        !Debug
        print*, "Fbar_oo = ", tensor_norm2(lambda_t%fbar_oo)
        print*, "Fbar_ov = ", tensor_norm2(lambda_t%fbar_ov)
        print*, "Fbar_vv = ", tensor_norm2(lambda_t%fbar_vv)

       end subroutine get_lambda_intermediates_talsh

end module
