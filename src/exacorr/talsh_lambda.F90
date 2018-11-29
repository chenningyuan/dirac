module talsh_lambda

!This module contains routines for solving the Lambda equations

    use talsh
    use tensor_algebra
    use exacorr_datatypes, only : talsh_comm_tens,talsh_intg_tens,talsh_lambda_tens
    use intermediates
    use talsh_common_routines

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0),MINUS_ONE=(-1.D0,0.D0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    private

    public solve_lambda_equations

    contains

      subroutine solve_lambda_equations(comm_t,int_t,l1_tensor,l2_tensor,nocc,nvir, &
                                            eps_occ,eps_vir,t_target_precision)

!      Routine for solving the lambda equations

       use talsh_ao_to_mo, only       : tensor_norm2

       integer, intent(in) :: nocc, nvir
!      fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
       type(talsh_comm_tens), intent(inout) :: comm_t
       type(talsh_intg_tens), intent(inout) :: int_t
       type(talsh_lambda_tens)              :: lambda_t
!      Mixed type intermediates
       type(talsh_tens_t) :: goo_tensor, gvv_tensor
!      solution tensors
       type(talsh_tens_t) :: l1_tensor, l2_tensor
       type(talsh_tens_t) :: s1_tensor, s2_tensor
!      Aux tensor
       type(talsh_tens_t) :: oo_aux_tensor
!      Tensor dimensions
       integer, dimension(2) :: oo_dims, ov_dims, vv_dims
       integer, dimension(4) :: oooo_dims,ovoo_dims,ooov_dims,vovo_dims,oovv_dims, &
                                vvvo_dims,vovv_dims,vvvv_dims
       integer               :: one_dims(1)
!      Orbital energies for scaling
       real(8), intent(in)   :: eps_occ(:), eps_vir(:)
       integer :: iccsd, nccsd=40
       real(8), intent(in) :: t_target_precision
       real(8) :: t_convergence

       integer :: ierr

       call print_date('Entered lambda solver routine')

!      Construct fixed intermediates
       oo_dims = nocc * 2
       ierr=talsh_tensor_construct(lambda_t%fbar_oo,C8,oo_dims,init_val=ZERO)
       ov_dims(1) = nocc * 2
       ov_dims(2) = nvir * 2
       ierr=talsh_tensor_construct(lambda_t%fbar_ov,C8,ov_dims,init_val=ZERO)
       vv_dims = nvir * 2
       ierr=talsh_tensor_construct(lambda_t%fbar_vv,C8,vv_dims,init_val=ZERO)
       oooo_dims = nocc * 2
       ierr=talsh_tensor_construct(lambda_t%w_oooo,C8,oooo_dims,init_val=ZERO)
       ovoo_dims(1)   = nocc * 2
       ovoo_dims(2)   = nvir * 2
       ovoo_dims(3:4) = nocc * 2
       ierr=talsh_tensor_construct(lambda_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
       ooov_dims(1:3) = nocc * 2
       ooov_dims(4)   = nvir * 2
       ierr=talsh_tensor_construct(lambda_t%w_ooov,C8,ooov_dims,init_val=ZERO)
       vovo_dims(1) = nvir * 2
       vovo_dims(2) = nocc * 2
       vovo_dims(3) = nvir * 2
       vovo_dims(4) = nocc * 2
       ierr=talsh_tensor_construct(lambda_t%w_vovo,C8,vovo_dims,init_val=ZERO)
       ierr=talsh_tensor_construct(lambda_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)
       vvvo_dims(1:3) = nvir * 2
       vvvo_dims(4)   = nocc * 2
       ierr=talsh_tensor_construct(lambda_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
       vovv_dims(1)   = nvir * 2
       vovv_dims(2)   = nocc * 2
       vovv_dims(3:4) = nvir * 2
       ierr=talsh_tensor_construct(lambda_t%w_vovv,C8,vovv_dims,init_val=ZERO)
       vvvv_dims = nvir * 2
       ierr=talsh_tensor_construct(lambda_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)

!      Get fixed intermediates
       call print_date('Calculating fixed intermediates')
       call get_lambda_intermediates(comm_t,int_t,lambda_t,nocc,nvir)
       call print_date('Finished calculating fixed intermediates')

!      Construct solution tensors
       ierr=talsh_tensor_construct(s1_tensor,C8,ov_dims,init_val=ZERO)
       oovv_dims(1:2) = nocc * 2
       oovv_dims(3:4) = nvir * 2
       ierr=talsh_tensor_construct(s2_tensor,C8,oovv_dims,init_val=ZERO)

!      Construct mixed type intermediates
       ierr=talsh_tensor_construct(goo_tensor,C8,oo_dims,init_val=ZERO)
       ierr=talsh_tensor_construct(gvv_tensor,C8,vv_dims,init_val=ZERO)

!      Construct auxiliary tensor
       ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)

       call print_date('Begin iterating')
       do iccsd = 1, nccsd

!       Form mixed intermediates
        ierr=talsh_tensor_init(gvv_tensor)
        ierr=talsh_tensor_contract("G(e,a)+=L(m,n,a,f)*T(f,e,m,n)",gvv_tensor,l2_tensor,comm_t%t2,scale=ONE_HALF)
        ierr=talsh_tensor_init(goo_tensor)
        ierr=talsh_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,comm_t%t2,scale=ONE_HALF)

        !if (debug) then
        if (.true.) then
           print*, "L1  = ", tensor_norm2(l1_tensor)
           print*, "L2  = ", tensor_norm2(l2_tensor)*ONE_QUARTER
           print*, "goo = ", tensor_norm2(goo_tensor)
           print*, "gvv = ", tensor_norm2(gvv_tensor)
        end if

!       Calculate s2 tensor
        ierr=talsh_tensor_init(s2_tensor)
        ierr=talsh_tensor_add("S(i,j,a,b)+=V(i,j,a,b)",s2_tensor,int_t%oovv)

!       Need to permute a and b, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)

        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,scale=ONE_HALF)

!       Need to permute i and j, as well as a and b, so term appears four times
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!       Need to permute a and b, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,gvv_tensor)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,gvv_tensor)

        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,goo_tensor)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,goo_tensor)

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_init(oo_aux_tensor)
        ierr=talsh_tensor_contract("O(i,m)+=L(i,e)*T(e,m)",oo_aux_tensor,l1_tensor,comm_t%t1)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,j,a,b)*O(i,m)",s2_tensor,int_t%oovv,oo_aux_tensor)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,i,b,a)*O(j,m)",s2_tensor,int_t%oovv,oo_aux_tensor)

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*V(e,j,a,b)",s2_tensor,l1_tensor,int_t%vovv)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*V(e,i,b,a)",s2_tensor,l1_tensor,int_t%vovv)

!       Need to permute i and j, as well as a and b, so term appears four times
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)

        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,scale=ONE_HALF)


!       Calculate s1 tensor
        ierr=talsh_tensor_init(s1_tensor)
        ierr=talsh_tensor_add("S(i,a)+=F(i,a)",s1_tensor,lambda_t%fbar_ov)
        ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)
        ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,scale=ONE_HALF)
        ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,gvv_tensor,lambda_t%w_vovv)
        ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,goo_tensor,lambda_t%w_ooov)
        ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,scale=ONE_HALF)


        if (ierr.ne.TALSH_SUCCESS) then
           print*," ccdriver: error in evaluation of lambda equations, code=",ierr
           return
        end if

!       prepare for next iteration
        ierr=talsh_tensor_init(l1_tensor)
        ierr=talsh_tensor_add("L1(i,a)+=S1(i,a)",l1_tensor,s1_tensor)
        ierr=talsh_tensor_init(l2_tensor)
        ierr=talsh_tensor_add("L2(i,j,a,b)+=S2(i,j,a,b)",l2_tensor,s2_tensor)
        call scale_with_denominators (eps_occ,eps_vir,nocc,l1_tensor,l2_tensor)

        call diis (l1_tensor,l2_tensor,t_convergence=t_convergence)

        write (*,'(A,I4,A,ES10.1)') "  Iteration ",iccsd," convergence =",t_convergence
        if (t_convergence .lt. t_target_precision) exit

       end do

       call print_date("Finished iterating")

       if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
       write(*,'(A,ES22.16)'),"  Final L1amp = ", tensor_norm2(l1_tensor)
       write(*,'(A,ES22.16)'),"  Final L2amp = ", tensor_norm2(l2_tensor)*ONE_QUARTER

       ierr=talsh_tensor_destruct(s1_tensor)
       ierr=talsh_tensor_destruct(s2_tensor)
       ierr=talsh_tensor_destruct(goo_tensor)
       ierr=talsh_tensor_destruct(gvv_tensor)
       ierr=talsh_tensor_destruct(oo_aux_tensor)
       ierr=talsh_tensor_destruct(lambda_t%fbar_oo)
       ierr=talsh_tensor_destruct(lambda_t%fbar_ov)
       ierr=talsh_tensor_destruct(lambda_t%fbar_vv)
       ierr=talsh_tensor_destruct(lambda_t%w_oooo)
       ierr=talsh_tensor_destruct(lambda_t%w_ooov)
       ierr=talsh_tensor_destruct(lambda_t%w_ovoo)
       ierr=talsh_tensor_destruct(lambda_t%w_vovo)
       ierr=talsh_tensor_destruct(lambda_t%w_vvvo)
       ierr=talsh_tensor_destruct(lambda_t%w_vovv)
       ierr=talsh_tensor_destruct(lambda_t%w_vvvv)
       ierr=talsh_tensor_destruct(lambda_t%wbar_vovo)

       call print_date('Leaving lambda solver routine')

      end subroutine solve_lambda_equations

end module talsh_lambda
