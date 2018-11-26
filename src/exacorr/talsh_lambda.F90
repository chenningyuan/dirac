module talsh_lambda

!This module contains routines for solving the Lambda equations

    use talsh
    use tensor_algebra
    use exacorr_datatypes, only : talsh_comm_tens,talsh_intg_tens,talsh_lambda_tens
    use intermediates
    use talsh_common_routines

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0)
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
       type(talsh_tens_t) :: xoo_tensor, xvv_tensor
!      solution tensors
       type(talsh_tens_t) :: l1_tensor, l2_tensor
       type(talsh_tens_t) :: s1_tensor, s2_tensor
!      Tensor dimensions
       integer, dimension(2) :: oo_dims, ov_dims, vv_dims
       integer, dimension(4) :: oooo_dims,ovoo_dims,ooov_dims,vovo_dims,oovv_dims, &
                                vvvo_dims,vovv_dims,vvvv_dims
       integer               :: one_dims(1)
!      Orbital energies for scaling
       real(8), intent(in)   :: eps_occ(:), eps_vir(:)
       integer :: iccsd, nccsd=30
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
       call get_lambda_intermediates(comm_t,int_t,lambda_t)
       call print_date('Finished calculating fixed intermediates')

!      Construct solution tensors
       ierr=talsh_tensor_construct(s1_tensor,C8,ov_dims,init_val=ZERO)
       oovv_dims(1:2) = nocc * 2
       oovv_dims(3:4) = nvir * 2
       ierr=talsh_tensor_construct(s2_tensor,C8,oovv_dims,init_val=ZERO)

!      Construct mixed type intermediates
       ierr=talsh_tensor_construct(xoo_tensor,C8,oo_dims,init_val=ZERO)
       ierr=talsh_tensor_construct(xvv_tensor,C8,vv_dims,init_val=ZERO)

       call print_date('Begin iterating')
       do iccsd = 1, nccsd

!       Form mixed intermediates
        ierr=talsh_tensor_init(xvv_tensor)
        ierr=talsh_tensor_contract("X(c,a)+=L(k,l,a,d)*T(d,c,k,l)",xvv_tensor,l2_tensor,comm_t%t2,scale=ONE_HALF)
        ierr=talsh_tensor_init(xoo_tensor)
        ierr=talsh_tensor_contract("X(i,k)+=L(i,l,c,d)*T(c,d,k,l)",xoo_tensor,l2_tensor,comm_t%t2,scale=ONE_HALF)

        !if (debug) then
        if (.true.) then
           print*, "L1  = ", tensor_norm2(l1_tensor)
           print*, "L2  = ", tensor_norm2(l2_tensor)*ONE_QUARTER
           print*, "xoo = ", tensor_norm2(xoo_tensor)
           print*, "xvv = ", tensor_norm2(xvv_tensor)
        end if

!       Calculate s2 tensor
        ierr=talsh_tensor_init(s2_tensor)



!       Calculate s1 tensor
        ierr=talsh_tensor_init(s1_tensor)



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

        write (*,'(A,I4,A,E10.1)') "  Iteration ",iccsd," convergence =",t_convergence
        if (t_convergence .lt. t_target_precision) exit

       end do

       ierr=talsh_tensor_destruct(s1_tensor)
       ierr=talsh_tensor_destruct(s2_tensor)
       ierr=talsh_tensor_destruct(xoo_tensor)
       ierr=talsh_tensor_destruct(xvv_tensor)
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
