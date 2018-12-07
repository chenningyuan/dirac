module talsh_lambda

!This module contains routines for solving the Lambda equations

    use talsh
    use tensor_algebra
    use exacorr_datatypes, only : talsh_comm_tens,talsh_intg_tens,talsh_lambda_tens
    use intermediates

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    private

    public solve_lambda_equations

    interface diis
        module procedure diis
    end interface

    contains

      subroutine solve_lambda_equations(comm_t,int_t,l1_tensor,l2_tensor,nocc,nvir, &
                                            eps_occ,eps_vir,t_target_precision,CCD)

!      Routine for solving the lambda equations

       use talsh_ao_to_mo, only       : tensor_norm2
       use talsh_common_routines, only : scale_with_denominators

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

       logical, intent(in) :: CCD

       integer :: ierr, prints=1

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

!      Obtain intermediates from CCSD ==> DO TIME TEST
       call print_date("Calculating common blocks")
       call get_tau(comm_t,CCD)
       call get_hvv(comm_t,int_t%oovv)
       call get_gvv(comm_t,int_t%vovv,CCD)
       call get_h(comm_t,int_t,CCD)
       call get_hoo(comm_t,int_t%oovv)
       call get_goo(comm_t,int_t%ooov,CCD)
       call get_hov(comm_t,int_t%oovv,CCD)
       call get_a(comm_t,int_t,CCD)
       call print_date("Finished calculating common blocks")

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

!      Re-initialize diis, so it can be used in lambda routine
       call diis (l1_tensor,l2_tensor,mode=1)

       call print_date('Begin iterating')
       do iccsd = 1, nccsd

!       Form mixed intermediates
        ierr=talsh_tensor_init(gvv_tensor)
        ierr=talsh_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,m,n)",gvv_tensor,l2_tensor,comm_t%t2,scale=MINUS_ONE_HALF)
        ierr=talsh_tensor_init(goo_tensor)
        ierr=talsh_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,comm_t%t2,scale=ONE_HALF)

        if (prints .ge. 1) then
           print*, "L1  = ", tensor_norm2(l1_tensor)
           print*, "L2  = ", tensor_norm2(l2_tensor)*ONE_QUARTER
           print*, "goo = ", tensor_norm2(goo_tensor)
           print*, "gvv = ", tensor_norm2(gvv_tensor)
        end if

!       Calculate s1 tensor
        ierr=talsh_tensor_init(s1_tensor)
        ierr=talsh_tensor_add("S(i,a)+=F(i,a)",s1_tensor,lambda_t%fbar_ov)
        if (prints .ge. 2) print*, "S1 init = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)
        if (prints .ge. 2) print*, "S1 after term1 = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,scale=MINUS_ONE)
        if (prints .ge. 2) print*, "S1 after term2 = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,goo_tensor,lambda_t%w_ooov)
        if (prints .ge. 2) print*, "S1 after term3 = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,scale=MINUS_ONE)
        if (prints .ge. 2) print*, "S1 after term4 = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,scale=ONE_HALF)
        if (prints .ge. 2) print*, "S1 after term5 = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,scale=ONE_HALF)
        if (prints .ge. 2) print*, "S1 after term6 = ", tensor_norm2(s1_tensor)
        ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,gvv_tensor,lambda_t%w_vovv)
        if (prints .ge. 2) print*, "S1 after term7 = ", tensor_norm2(s1_tensor)

!       Calculate s2 tensor
        ierr=talsh_tensor_init(s2_tensor)
        ierr=talsh_tensor_add("S(i,j,a,b)+=V(i,j,a,b)",s2_tensor,int_t%oovv)
        if (prints .ge. 2) print*, "S2 after initialization = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute a and b, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
        if (prints .ge. 2) print*, "S2 after term1 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
        if (prints .ge. 2) print*, "S2 after term2 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,scale=ONE_HALF)
        if (prints .ge. 2) print*, "S2 after term3 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute i and j, as well as a and b, so term appears four times
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        !if (prints .ge. 2) print*, "S2 after term4 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute a and b, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,gvv_tensor)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,gvv_tensor)
        if (prints .ge. 2) print*, "S2 after term5 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute a and b, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)
        if (prints .ge. 2) print*, "S2 after term6 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,goo_tensor)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,goo_tensor)
        if (prints .ge. 2) print*, "S2 after term7 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute i and j, so term appears twice
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)
        if (prints .ge. 2) print*, "S2 after term8 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

!       Need to permute i and j, as well as a and b, so term appears four times
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
        if (prints .ge. 2) print*, "S2 after term9 = ", tensor_norm2(s2_tensor)*ONE_QUARTER

        ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,scale=ONE_HALF)
        if (prints .ge. 2) print*, "S2 after term10 = ", tensor_norm2(s2_tensor)*ONE_QUARTER


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

!      Destruct tensors
       call diis (l1_tensor,l2_tensor,mode=2)
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

      subroutine diis (l1_tensor,l2_tensor,t_convergence,mode)

       integer, parameter :: maxdim = 5 ! maximum number of tensors that will be stored 
       integer, parameter :: mindim = 3 ! minimum number of tensors to allow extrapolation
       integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
       real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

       type(talsh_tens_t), intent(inout)  :: l1_tensor, l2_tensor
       integer            :: ierr, i, tens_rank
       integer, save      :: l1_dims(1:2), l2_dims(1:4)
       logical            :: init, clean_up
       integer, save      :: itens, jtens, ntens
       complex(8), save   :: b_diis(maxdim,maxdim)
       complex(8)         :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
       real(8)            :: b_max, b_scale
       type(talsh_tens_t), save :: l1_inp_tensor, l2_inp_tensor
       type(talsh_tens_t), save :: l1_res_tensor(maxdim), l2_res_tensor(maxdim)
       type(talsh_tens_t), save :: l1_err_tensor(maxdim), l2_err_tensor(maxdim)
       integer            :: ipiv(0:maxdim), info  ! needed for ZHESV solver
       real(8)            :: work(maxdim,maxdim)   ! needed for ZHESV solver
       type(talsh_tens_t), save  :: scalar_tensor
       integer:: scalar_dims(1)
       complex(8), pointer :: scalar_tens(:)
       type(C_PTR):: body_p

       init       = .false.
       clean_up   = .false.

       if (present(mode)) then
          if (mode.eq.1) init       = .true.
          if (mode.eq.2) clean_up   = .true.
       end if 

!      define scalar tensor needed to extract the error
       scalar_dims(1) = 1
       ierr=talsh_tensor_construct(scalar_tensor,C8,scalar_dims(1:0),init_val=ZERO)
       ierr=talsh_tensor_get_body_access(scalar_tensor,body_p,C8,0,DEV_HOST)
       call c_f_pointer(body_p,scalar_tens,scalar_dims)

       if (init) then
!         get dimensions of the tensors
          ierr = talsh_tensor_dimensions(l1_tensor,tens_rank,l1_dims)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
          ierr = talsh_tensor_dimensions(l2_tensor,tens_rank,l2_dims)
          if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!         initialize DIIS matrix
          b_diis      = ZERO
!         create storage for input, result and error tensors (in memory for the time being)
          ierr=talsh_tensor_construct(l1_inp_tensor,C8,l1_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(l2_inp_tensor,C8,l2_dims,init_val=ZERO)
          do i = 1, maxdim
             ierr=talsh_tensor_construct(l1_res_tensor(i),C8,l1_dims,init_val=ZERO)
             ierr=talsh_tensor_construct(l1_err_tensor(i),C8,l1_dims,init_val=ZERO)
             ierr=talsh_tensor_construct(l2_res_tensor(i),C8,l2_dims,init_val=ZERO)
             ierr=talsh_tensor_construct(l2_err_tensor(i),C8,l2_dims,init_val=ZERO)
          end do
          itens = 0
          ntens = 0
!         store the tensors in memory as current input tensors, needed to calculate the error vector later on
          ierr=talsh_tensor_init(l1_inp_tensor)
          ierr=talsh_tensor_add("L1O(i,a)+=L1N(i,a)",l1_inp_tensor,l1_tensor)
          ierr=talsh_tensor_init(l2_inp_tensor)
          ierr=talsh_tensor_add("L2O(i,j,a,b)+=L2N(i,j,a,b)",l2_inp_tensor,l2_tensor)
          return ! no extrapolation possible
       end if

       if (clean_up) then
!         de-allocate tensors
          ierr=talsh_tensor_destruct(l1_inp_tensor)
          ierr=talsh_tensor_destruct(l2_inp_tensor)
          do i = 1, maxdim
             ierr=talsh_tensor_destruct(l1_res_tensor(i))
             ierr=talsh_tensor_destruct(l1_err_tensor(i))
             ierr=talsh_tensor_destruct(l2_res_tensor(i))
             ierr=talsh_tensor_destruct(l2_err_tensor(i))
          end do
          return
       end if

!      update counters
       itens = itens + 1                                                                                                                    
       if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
       ntens = ntens + 1                    ! number of tensors stored after returning
       if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!      store the result tensors
       ierr=talsh_tensor_init(l1_res_tensor(itens))
       ierr=talsh_tensor_add("L1R(i,a)+=L1(i,a)",l1_res_tensor(itens),l1_tensor)
       ierr=talsh_tensor_init(l2_res_tensor(itens))
       ierr=talsh_tensor_add("L2R(i,j,a,b)+=L2(i,j,a,b)",l2_res_tensor(itens),l2_tensor)

!      calculate the error tensors by comparing the result of the t-equations with the tensors used as input
       ierr=talsh_tensor_init(l1_err_tensor(itens))
       ierr=talsh_tensor_add("L1E(i,a)+=L1R(i,a)",l1_err_tensor(itens),l1_tensor)
       ierr=talsh_tensor_add("L1E(i,a)+=L1I(i,a)",l1_err_tensor(itens),l1_inp_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_init(l2_err_tensor(itens))
       ierr=talsh_tensor_add("L2E(i,j,a,b)+=L2R(i,j,a,b)",l2_err_tensor(itens),l2_tensor)
       ierr=talsh_tensor_add("L2E(i,j,a,b)+=L2I(i,j,a,b)",l2_err_tensor(itens),l2_inp_tensor,scale=MINUS_ONE)

!      update the b_diis array that holds the scalar products of the error tensors
       ierr=talsh_tensor_init(scalar_tensor)
       do i = 1, ntens
          ierr=talsh_tensor_init(scalar_tensor)
          ierr=talsh_tensor_contract("E()+=L+(i,a)*L(i,a)",scalar_tensor,l1_err_tensor(itens),l1_err_tensor(i))
          ierr=talsh_tensor_contract("E()+=L+(i,j,a,b)*L(i,j,a,b)",scalar_tensor,l2_err_tensor(itens),l2_err_tensor(i))
          b_diis(itens,i) = scalar_tens(1)
          b_diis(i,itens) = dconjg(scalar_tens(1))
       end do

       if (present(t_convergence)) then
          t_convergence = dsqrt(dble(b_diis(itens,itens)))
       end if 

!      transfer to a work array in the way requested by the solver and normalize the matrix
       b_max = maxval (abs(b_diis))
       if (b_max .gt. 1.D-30) then
           b_scale = 1.D0 / b_max
       else
           print*, "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
           return
       end if
       bb = ZERO
       bb(1:ntens,0) = MINUS_ONE
       bb(0,1:ntens) = MINUS_ONE
       bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!      solve the linear system of equations using the LAPACK solver ZHESV
       coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
       coeff(1:ntens) = ZERO 
       call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!      make the extrapolated amplitudes if the space is large enough
       if (ntens.ge.mindim) then
         ierr=talsh_tensor_init(l1_tensor)
         ierr=talsh_tensor_init(l2_tensor)
          do i = 1, ntens
             ierr=talsh_tensor_add("L1(i,a)+=L1O(i,a)",l1_tensor,l1_res_tensor(i),scale=coeff(i))
             ierr=talsh_tensor_add("L2(i,j,a,b)+=L2O(i,j,a,b)",l2_tensor,l2_res_tensor(i),scale=coeff(i))
          end do
       end if

!      the extrapolated amplitudes are the new input tensors
       ierr=talsh_tensor_init(l1_inp_tensor)
       ierr=talsh_tensor_add("L1O(i,a)+=L1(i,a)",l1_inp_tensor,l1_tensor)
       ierr=talsh_tensor_init(l2_inp_tensor)
       ierr=talsh_tensor_add("L2O(i,j,a,b)+=L2(i,j,a,b)",l2_inp_tensor,l2_tensor)

       ierr=talsh_tensor_destruct(scalar_tensor)

      end subroutine diis



end module talsh_lambda
