module intermediates

!Module to make intermediates needed in CC and Lambda equations

     implicit none
     private

     complex(8), parameter :: ZERO=(0.D0,0.D0), MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0), &
                              ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                              MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0)
    real(8), parameter     :: ONE_QUARTER=0.25D0

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
        use talsh_ao_to_mo, only : tensor_norm2

        type(talsh_comm_tens),intent(inout) :: comm_t
        type(talsh_intg_tens),intent(inout) :: int_t

        integer :: ierr

!       tau ! validated
        ierr=talsh_tensor_init(comm_t%tau)
        ierr=talsh_tensor_add("T(a,b,i,j)+=R(a,b,i,j)",comm_t%tau,comm_t%t2)
        ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1)                  ! o(2)v(2)
        ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,scale=MINUS_ONE)  ! o(2)v(2)

!       vv ! validated
        ierr=talsh_tensor_init(comm_t%hvv)
        ierr=talsh_tensor_add("H(a,c)+=F(a,c)",comm_t%hvv,comm_t%fvv)
        ierr=talsh_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%hvv,int_t%oovv,comm_t%tau,scale=ONE_HALF) ! o(2)v(3)
        ierr=talsh_tensor_init(comm_t%gvv)
        ierr=talsh_tensor_add("G(a,c)+=H(a,c)",comm_t%gvv,comm_t%hvv)
        ierr=talsh_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,int_t%vovv,comm_t%t1)             ! o(1)v(3)
        ierr=talsh_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,scale=MINUS_ONE) ! o(1)v(2)

!       vovo ! more tests
        ierr=talsh_tensor_init(comm_t%h_int)
        ierr=talsh_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1)                 ! o(2)v(3)
        ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
        ierr=talsh_tensor_add("H(a,k,c,i)+=V(a,k,c,i)",comm_t%h_int,int_t%vovo)

!       oo ! validated
        ierr=talsh_tensor_init(comm_t%hoo)
        ierr=talsh_tensor_add("H(k,i)+=F(k,i)",comm_t%hoo,comm_t%foo,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("G(k,i)+=V(k,l,c,d)*T(c,d,l,i)",comm_t%hoo,int_t%oovv,comm_t%tau,scale=ONE_HALF) ! o(3)v(2)
        ierr=talsh_tensor_init(comm_t%goo)
        ierr=talsh_tensor_add("G(k,i)+=H(k,i)",comm_t%goo,comm_t%hoo)
        ierr=talsh_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,int_t%ooov,comm_t%t1,scale=MINUS_ONE)     ! o(3)v(1)
        ierr=talsh_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1)                         ! o(2)v(1)

!       ov ! validated
        ierr=talsh_tensor_init(comm_t%hov)
        ierr=talsh_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,int_t%oovv,comm_t%t1) ! o(2)v(2)
        ierr=talsh_tensor_add("H(k,c)+=F(k,c)",comm_t%hov,comm_t%fov)

!       oooo ! validated
        ierr=talsh_tensor_init(comm_t%a_int)
        ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",comm_t%a_int,int_t%oooo)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%a_int,int_t%oovv,comm_t%tau, &
                                    scale=ONE_HALF)

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

       subroutine get_lambda_intermediates_talsh(comm_t,int_t,lambda_t,nocc,nvir)

!       Get fixed intermediates of lambda equations

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens
        use talsh_ao_to_mo, only : tensor_norm2

        type(talsh_comm_tens),intent(inout)   :: comm_t
        type(talsh_intg_tens),intent(inout)   :: int_t
        type(talsh_lambda_tens),intent(inout) :: lambda_t
        type(talsh_tens_t)                    :: other1_tensor
        type(talsh_tens_t)                    :: tau2_tensor

        integer, intent(in) :: nocc, nvir
        integer             :: one_dims(1)
        integer             :: vvoo_dims(4)

        integer :: ierr

!       Create and initialize scalars that are to be used as tensors in contractions
        one_dims = 1
        ierr=talsh_tensor_construct(other1_tensor,C8,one_dims(1:0),init_val=ONE)

        !Auxiliary tensor
        vvoo_dims(1:2) = nvir * 2
        vvoo_dims(3:4) = nocc * 2
        ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)

!       One-body intermediates
!       vv
        ierr=talsh_tensor_add("F(e,a)+=G(e,a)",lambda_t%fbar_vv,comm_t%gvv)
        ierr=talsh_tensor_contract("F(e,a)+=V(m,e,f,a)*T(f,m)",lambda_t%fbar_vv,int_t%vovv,comm_t%t1,scale=MINUS_TWO)
!       ov
        ierr=talsh_tensor_add("F(m,e)+=H(m,e)",lambda_t%fbar_ov,comm_t%hov)
!       oo
        ierr=talsh_tensor_add("F(i,m)+=G(i,m)",lambda_t%fbar_oo,comm_t%goo,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("F(i,m)+=R(i,e)*T(e,m)",lambda_t%fbar_oo,comm_t%fov,comm_t%t1,scale=TWO)

!       Two-body intermediates
!       vvvv
!       Permute e and f, so term appears twice
        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",lambda_t%w_vvvv,int_t%vovv,comm_t%t1)
        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",lambda_t%w_vvvv,int_t%vovv,comm_t%t1,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",lambda_t%w_vvvv,int_t%oovv,comm_t%tau,scale=ONE_HALF)
        ierr=talsh_tensor_add("W(e,f,a,b)+=V(e,f,a,b)",lambda_t%w_vvvv,int_t%vvvv)

!       vovo bar
        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(b,f,n,j)",lambda_t%wbar_vovo,int_t%oovv,comm_t%t2)
        ierr=talsh_tensor_add("W(b,m,e,j)+=V(b,m,e,j)",lambda_t%wbar_vovo,int_t%vovo)

!       vvvo
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V+(a,m,e,f)*K()",lambda_t%w_vvvo,int_t%vovv,other1_tensor)
!       Need to permute e and f, so term appears twice
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(e,n,a,g)*T(f,g,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(f,n,a,g)*T(g,e,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=F(n,a)*T(e,f,m,n)",lambda_t%w_vvvo,lambda_t%fbar_ov,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(n,o,m,a)*T(e,f,o,n)",lambda_t%w_vvvo,int_t%ooov,comm_t%tau,scale=ONE_HALF)
!       Need to permute e and f, so term appears twice
        ierr=talsh_tensor_contract("W(e,f,a,m)+=R(f,n,a,m)*T(e,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,n,a,m)*T(f,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,f,a,g)*T(g,m)",lambda_t%w_vvvo,lambda_t%w_vvvv,comm_t%t1)

!       vovv
        ierr=talsh_tensor_contract("W(a,m,e,f)+=V(m,n,e,f)*T(a,n)",lambda_t%w_vovv,int_t%oovv,comm_t%t1)
        ierr=talsh_tensor_add("W(a,m,e,f)+=V(a,m,e,f)",lambda_t%w_vovv,int_t%vovv)

!       vovo
        ierr=talsh_tensor_add("W(b,m,e,j)+=H(b,m,e,j)",lambda_t%w_vovo,comm_t%h_int)
        ierr=talsh_tensor_init(tau2_tensor)
        ierr=talsh_tensor_add("T(f,b,j,n)+=R(f,b,j,n)",tau2_tensor,comm_t%t2)
        ierr=talsh_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",tau2_tensor,comm_t%t1,comm_t%t1)
        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",lambda_t%w_vovo,int_t%oovv,tau2_tensor)

!       oooo
        ierr=talsh_tensor_add("W(i,j,m,n)+=A(i,j,m,n)",lambda_t%w_oooo,comm_t%a_int)

!       ovoo
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V+(m,n,i,e)*K()",lambda_t%w_ovoo,int_t%ooov,other1_tensor)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=F(i,f)*T(f,e,m,n)",lambda_t%w_ovoo,lambda_t%fbar_ov,comm_t%t2)
!       Need to permute m and n, so term appears twice
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(i,o,m,f)*T(e,f,n,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(o,i,n,f)*T(e,f,m,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
!       Need to permute m and n, so term appears twice
        ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,n)*T(f,m)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,m)*T(f,n)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=R(o,i,m,n)*T(e,o)",lambda_t%w_ovoo,lambda_t%w_oooo,comm_t%t1)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(e,i,f,g)*T(g,f,m,n)",lambda_t%w_ovoo,int_t%vovv,comm_t%tau,scale=ONE_HALF)

!       ooov
        ierr=talsh_tensor_contract("W(m,n,i,e)+=T(f,i)*V(m,n,f,e)",lambda_t%w_ooov,comm_t%t1,int_t%oovv)
        ierr=talsh_tensor_add("W(m,n,i,e)+=V(m,n,i,e)",lambda_t%w_ooov,int_t%ooov)

        !Debug
        print*, "T1amp      = ", tensor_norm2(comm_t%t1)
        print*, "T2amp      = ", tensor_norm2(comm_t%t2)*ONE_QUARTER
!        print*, "Tau        = ", tensor_norm2(comm_t%tau)*ONE_QUARTER
!        print*, "Fbar_oo    = ", tensor_norm2(lambda_t%fbar_oo)*0.5D0
        print*, "Fbar_ov    = ", tensor_norm2(lambda_t%fbar_ov)
!        print*, "Fbar_vv    = ", tensor_norm2(lambda_t%fbar_vv)*0.5D0
!        print*, "w_vvvv     = ", tensor_norm2(lambda_t%w_vvvv)*ONE_QUARTER
        print*, "wbar_vovo  = ", tensor_norm2(lambda_t%wbar_vovo)
        print*, "w_vvvo     = ", tensor_norm2(lambda_t%w_vvvo)*0.5D0
        print*, "w_vovv     = ", tensor_norm2(lambda_t%w_vovv)*0.5D0
        print*, "w_vovo     = ", tensor_norm2(lambda_t%w_vovo)
        print*, "w_ovoo     = ", tensor_norm2(lambda_t%w_ovoo)*0.5D0
        print*, "w_ooov     = ", tensor_norm2(lambda_t%w_ooov)*0.5D0
        print*, "w_oooo     = ", tensor_norm2(lambda_t%w_oooo)*ONE_QUARTER

        ierr=talsh_tensor_destruct(other1_tensor)
        ierr=talsh_tensor_destruct(tau2_tensor)

       end subroutine get_lambda_intermediates_talsh

end module
