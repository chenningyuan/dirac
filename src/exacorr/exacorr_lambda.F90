module exacorr_lambda

!This module contains routines for solving the Lambda equations

#if defined (VAR_MPI)

    use exatensor 
    use exacorr_datatypes, only : exatns_comm_tens,exatns_intg_tens,exatns_lambda_tens
    use intermediates

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    private

    public solve_lambda_equations

    interface diis
      module procedure diis
    end interface

    contains

      subroutine solve_lambda_equations(comm_t,int_t,l1_tensor,l2_tensor,occ_space_id,occ_space_root, &
                                        vir_space_id,vir_space_root,CCD,t_target_precision)

!      Routine for solving the lambda equations

       use exacorr_ao_to_mo, only       : print_tensornorm2

!      fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
       type(exatns_comm_tens), intent(inout) :: comm_t
       type(exatns_intg_tens), intent(inout) :: int_t
       type(exatns_lambda_tens)              :: lambda_t
!      Mixed type intermediates
       type(tens_rcrsv_t) :: goo_tensor, gvv_tensor
!      solution tensors
       type(tens_rcrsv_t) :: l1_tensor, l2_tensor
       type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!      Aux tensor
       type(tens_rcrsv_t) :: oo_aux_tensor
       type(tens_rcrsv_t) :: one_l_tensor
!      Tensor dimensions
       integer, intent(in) :: occ_space_id, vir_space_id
       integer(INTL), intent(in) :: occ_space_root, vir_space_root
       integer, dimension(2) :: oo_id, ov_id, vv_id
       integer, dimension(4) :: oooo_id,ovoo_id,ooov_id,vovo_id,oovv_id, &
                                vvvo_id,vovv_id,vvvv_id
       integer(INTL), dimension(2) :: oo_root, ov_root, vv_root
       integer(INTL), dimension(4) :: oooo_root,ovoo_root,ooov_root,vovo_root,oovv_root, &
                                      vvvo_root,vovv_root,vvvv_root
       integer :: iccsd, nccsd=30
       real(8), intent(in) :: t_target_precision
       real(8) :: t_convergence

       logical, intent(in) :: CCD
       integer :: ierr, prints=1

       call print_date('Entered lambda solver routine')

!      Create all ids and roots
       oo_id = (/occ_space_id, occ_space_id/)
       ov_id = (/occ_space_id, vir_space_id/)
       vv_id = (/vir_space_id, vir_space_id/)
       oo_root = (/occ_space_root, occ_space_root/)
       ov_root = (/occ_space_root, vir_space_root/)
       vv_root = (/vir_space_root, vir_space_root/)
       oooo_id = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
       ovoo_id = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
       ooov_id = (/occ_space_id, occ_space_id, occ_space_id, vir_space_id/)
       oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
       vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
       vvvo_id = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
       vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
       vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
       oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
       ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
       ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
       oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
       vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
       vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
       vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
       vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!      Construct fixed intermediates
       ierr=exatns_tensor_create(lambda_t%fbar_oo,"fbar_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%fbar_ov,"fbar_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%fbar_vv,"fbar_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_oooo,"w_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_ovoo,"w_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_ooov,"w_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_vovo,"w_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%wbar_vovo,"wbar_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_vvvo,"w_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_vovv,"w_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(lambda_t%w_vvvv,"w_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)

!      Obtain intermediates from CCSD
       call print_date("Calculating common blocks")
       call get_tau(comm_t,CCD,one_l_tensor)
       call get_hvv(comm_t,int_t%oovv,one_l_tensor)
       call get_gvv(comm_t,int_t%vovv,CCD,one_l_tensor)
       call get_h(comm_t,int_t,CCD,one_l_tensor)
       call get_hoo(comm_t,int_t%oovv,one_l_tensor)
       call get_goo(comm_t,int_t%ooov,CCD,one_l_tensor)
       call get_hov(comm_t,int_t%oovv,CCD,one_l_tensor)
       call get_a(comm_t,int_t,CCD,one_l_tensor)
       call print_date("Finished calculating common blocks")

!      Get fixed intermediates
       call print_date('Calculating fixed intermediates')
       call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD,one_l_tensor)
       call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t,one_l_tensor)
       call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov,one_l_tensor)
       call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_l_tensor)
       call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD,one_l_tensor)
       call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,occ_space_id,occ_space_root, &
                        vir_space_id,vir_space_root,CCD,one_l_tensor)
       call get_w_oooo(lambda_t%w_oooo,comm_t%a_int,one_l_tensor)
       call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_l_tensor)
       call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD,one_l_tensor)
       call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv,one_l_tensor)
       call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD,one_l_tensor)
       call print_date('Finished calculating fixed intermediates')

       !Debug
       print*, "T1amp      = ", print_tensornorm2(comm_t%t1)
       print*, "T2amp      = ", print_tensornorm2(comm_t%t2)*ONE_QUARTER
!       print*, "Tau        = ", print_tensornorm2(comm_t%tau)*ONE_QUARTER
!       print*, "Fbar_oo    = ", print_tensornorm2(lambda_t%fbar_oo)*0.5D0
       print*, "Fbar_ov    = ", print_tensornorm2(lambda_t%fbar_ov)
!       print*, "Fbar_vv    = ", print_tensornorm2(lambda_t%fbar_vv)*0.5D0
!       print*, "w_vvvv     = ", print_tensornorm2(lambda_t%w_vvvv)*ONE_QUARTER
       print*, "wbar_vovo  = ", print_tensornorm2(lambda_t%wbar_vovo)
       print*, "w_vvvo     = ", print_tensornorm2(lambda_t%w_vvvo)*0.5D0
       print*, "w_vovv     = ", print_tensornorm2(lambda_t%w_vovv)*0.5D0
       print*, "w_vovo     = ", print_tensornorm2(lambda_t%w_vovo)
       print*, "w_ovoo     = ", print_tensornorm2(lambda_t%w_ovoo)*0.5D0
       print*, "w_ooov     = ", print_tensornorm2(lambda_t%w_ooov)*0.5D0
       print*, "w_oooo     = ", print_tensornorm2(lambda_t%w_oooo)*ONE_QUARTER

!      Construct solution tensors
       ierr=exatns_tensor_create(s1_tensor,"s1",oovv_id,oovv_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(s2_tensor,"s2",ov_id,ov_root,EXA_DATA_KIND_C8)

!      Construct mixed type intermediates
       ierr=exatns_tensor_create(goo_tensor,"goo_mixed",oo_id,oo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(gvv_tensor,"gvv_mixed",vv_id,vv_root,EXA_DATA_KIND_C8)

!      Construct auxiliary tensor
       ierr=exatns_tensor_create(oo_aux_tensor,"oo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(one_l_tensor,"one_l",EXA_DATA_KIND_C8)
       ierr=exatns_tensor_init(one_l_tensor,ONE)


!      Reinitialize diis, so it can be used in lambda routine
       call diis (l2_tensor,oovv_id,oovv_root,l1_tensor,ov_id,ov_root,mode=1)

       call print_date('Begin iterating')
       do iccsd = 1, nccsd

!       Form mixed intermediates
        ierr=exatns_tensor_init(gvv_tensor)
        ierr=exatns_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,m,n)",gvv_tensor,l2_tensor,comm_t%t2,MINUS_ONE_HALF)
        ierr=exatns_tensor_init(goo_tensor)
        ierr=exatns_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,comm_t%t2,ONE_HALF)

        if (prints .ge. 1) then
           print*, "L1  = ", print_tensornorm2(l1_tensor)
           print*, "L2  = ", print_tensornorm2(l2_tensor)*ONE_QUARTER
           print*, "goo = ", print_tensornorm2(goo_tensor)
           print*, "gvv = ", print_tensornorm2(gvv_tensor)
        end if

!       Calculate s1 tensor
        ierr=exatns_tensor_init(s1_tensor)
        ierr=exatns_tensor_contract("S(i,a)+=F(i,a)*K()",s1_tensor,lambda_t%fbar_ov,one_l_tensor)
        ierr=exatns_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)
        ierr=exatns_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,MINUS_ONE)
        ierr=exatns_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,goo_tensor,lambda_t%w_ooov)
        ierr=exatns_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,MINUS_ONE)
        ierr=exatns_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,ONE_HALF)
        ierr=exatns_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,ONE_HALF)
        ierr=exatns_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,gvv_tensor,lambda_t%w_vovv)

!       Calculate s2 tensor
        ierr=exatns_tensor_init(s2_tensor)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_l_tensor)

!       Need to permute a and b, so term appears twice
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)

!       Need to permute i and j, so term appears twice
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)

        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,ONE_HALF)

!       Need to permute i and j, as well as a and b, so term appears four times
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!       Need to permute a and b, so term appears twice
        ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,gvv_tensor)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,gvv_tensor)

!       Need to permute a and b, so term appears twice
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)

!       Need to permute i and j, so term appears twice
        ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,goo_tensor)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,goo_tensor)

!       Need to permute i and j, so term appears twice
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)

!       Need to permute i and j, as well as a and b, so term appears four times
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)

        ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,ONE_HALF)


        if (ierr.ne.EXA_SUCCESS) then
           print*," ccdriver: error in evaluation of lambda equations, code=",ierr
           return
        end if

!       prepare for next iteration
        ierr=exatns_tensor_init(l1_tensor)
        ierr=exatns_tensor_contract("L1(i,a)+=S1(i,a)*K()",l1_tensor,s1_tensor,one_l_tensor)
        ierr=exatns_tensor_transform(l1_tensor,"Denominate")
        ierr=exatns_tensor_init(l2_tensor)
        ierr=exatns_tensor_contract("L2(i,j,a,b)+=S2(i,j,a,b)*K()",l2_tensor,s2_tensor,one_l_tensor)
        ierr=exatns_tensor_transform(l2_tensor,"Denominate")

        call diis (l2_tensor,l1_tensor=l1_tensor,t_convergence=t_convergence)

        write (*,'(A,I4,A,ES10.1)') "  Iteration ",iccsd," convergence =",t_convergence
        if (t_convergence .lt. t_target_precision) exit

       end do

       call print_date("Finished iterating")

       if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
       write(*,'(A,ES22.16)'),"  Final L1amp = ", print_tensornorm2(l1_tensor)
       write(*,'(A,ES22.16)'),"  Final L2amp = ", print_tensornorm2(l2_tensor)*ONE_QUARTER

!      Destroy diis tensors
       call diis(l2_tensor,l1_tensor=l1_tensor,mode=2)
       ierr=exatns_tensor_destroy(s1_tensor)
       ierr=exatns_tensor_destroy(s2_tensor)
       ierr=exatns_tensor_destroy(goo_tensor)
       ierr=exatns_tensor_destroy(gvv_tensor)
       ierr=exatns_tensor_destroy(oo_aux_tensor)
       ierr=exatns_tensor_destroy(one_l_tensor)
       ierr=exatns_tensor_destroy(lambda_t%fbar_oo)
       ierr=exatns_tensor_destroy(lambda_t%fbar_ov)
       ierr=exatns_tensor_destroy(lambda_t%fbar_vv)
       ierr=exatns_tensor_destroy(lambda_t%w_oooo)
       ierr=exatns_tensor_destroy(lambda_t%w_ooov)
       ierr=exatns_tensor_destroy(lambda_t%w_ovoo)
       ierr=exatns_tensor_destroy(lambda_t%w_vovo)
       ierr=exatns_tensor_destroy(lambda_t%w_vvvo)
       ierr=exatns_tensor_destroy(lambda_t%w_vovv)
       ierr=exatns_tensor_destroy(lambda_t%w_vvvv)
       ierr=exatns_tensor_destroy(lambda_t%wbar_vovo)

       call print_date('Leaving lambda solver routine')

      end subroutine solve_lambda_equations

      subroutine diis (l2_tensor,oovv_id,oovv_root,l1_tensor,ov_id,ov_root,t_convergence,mode)

       integer, parameter             :: maxdim = 5 ! maximum number of tensors that will be stored
       integer, parameter             :: mindim = 3 ! minimum number of tensors to allow extrapolation
       integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
       real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

       type(tens_rcrsv_t), intent(inout),optional :: l1_tensor
       type(tens_rcrsv_t), intent(inout)          :: l2_tensor
       integer                                    :: ierr, i, tens_rank, j
       integer,dimension(2),optional              :: ov_id
       integer,dimension(4),optional              :: oovv_id
       integer(INTL),dimension(2),optional        :: ov_root
       integer(INTL),dimension(4),optional        :: oovv_root
       logical                                    :: init, clean_up
       integer, save                              :: itens, jtens, ntens
       complex(8), save                           :: b_diis(maxdim,maxdim)
       complex(8)                                 :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
       complex(8)                                 :: scalar
       real(8)                                    :: b_max, b_scale
       type(tens_rcrsv_t), save                   :: l1_inp_tensor, l2_inp_tensor
       type(tens_rcrsv_t), save                   :: l1_res_tensor(maxdim), l2_res_tensor(maxdim)
       type(tens_rcrsv_t), save                   :: l1_err_tensor(maxdim), l2_err_tensor(maxdim)
       integer                                    :: ipiv(0:maxdim), info  ! needed for ZHESV solver
       real(8)                                    :: work(maxdim,maxdim)   ! needed for ZHESV solver
       type(tens_rcrsv_t), save                   :: scalar_tensor
       type(tens_rcrsv_t)                         :: one_diis_tensor

       character(len=1) :: vs_str

       init       = .false.
       clean_up   = .false.

       if (present(mode)) then
          if (mode.eq.1) init       = .true.
          if (mode.eq.2) clean_up   = .true.
       end if

       ierr=exatns_tensor_create(one_diis_tensor,"one_diis",EXA_DATA_KIND_C8)
       ierr=exatns_tensor_init(one_diis_tensor,ONE)

       if (init) then
!         get dimensions of the tensors
          tens_rank=l1_tensor%get_rank(ierr)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
          tens_rank=l2_tensor%get_rank(ierr)
          if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!         initialize DIIS matrix
          b_diis      = ZERO
!         create storage for input, result and error tensors (in memory for the time being)
          ierr=exatns_tensor_create(l1_inp_tensor,"l1_inp_tensor",ov_id,ov_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_create(l2_inp_tensor,"l2_inp_tensor",oovv_id,oovv_root,EXA_DATA_KIND_C8)
          do i = 1, maxdim
             write (vs_str, '(I1)') i
             ierr=exatns_tensor_create(l1_res_tensor(i),"l1_res_tensor_" // vs_str,ov_id,ov_root,EXA_DATA_KIND_C8)
             ierr=exatns_tensor_create(l1_err_tensor(i),"l1_err_tensor_" // vs_str,ov_id,ov_root,EXA_DATA_KIND_C8)
             ierr=exatns_tensor_create(l2_res_tensor(i),"l2_res_tensor_" // vs_str,oovv_id,oovv_root,EXA_DATA_KIND_C8)
             ierr=exatns_tensor_create(l2_err_tensor(i),"l2_err_tensor_" // vs_str,oovv_id,oovv_root,EXA_DATA_KIND_C8)
          end do
          if(ierr.ne.EXA_SUCCESS) call quit('failed to create storage for input')
          itens = 0
          ntens = 0
!         store the tensors in memory as current input tensors, needed to calculate the error vector later on
          ierr=exatns_tensor_init(l1_inp_tensor)
          ierr=exatns_tensor_contract("L1O(i,a)+=L1N(i,a)*K()",l1_inp_tensor,l1_tensor,one_diis_tensor)
          ierr=exatns_tensor_init(l2_inp_tensor)
          ierr=exatns_tensor_contract("L2O(i,j,a,b)+=L2N(i,j,a,b)*K()",l2_inp_tensor,l2_tensor,one_diis_tensor)
          ierr=exatns_tensor_create(scalar_tensor,"scalar_tensor",EXA_DATA_KIND_C8)
          return ! no extrapolation possible
       end if

       if (clean_up) then
!         de-allocate tensors
          ierr=exatns_tensor_destroy(l1_inp_tensor)
          ierr=exatns_tensor_destroy(l2_inp_tensor)
          do i = 1, maxdim
              ierr=exatns_tensor_destroy(l1_res_tensor(i))
              ierr=exatns_tensor_destroy(l1_err_tensor(i))
              ierr=exatns_tensor_destroy(l2_res_tensor(i))
              ierr=exatns_tensor_destroy(l2_err_tensor(i))
          end do
          ierr=exatns_tensor_destroy(scalar_tensor)
          return
       end if

!      update counters
       itens = itens + 1
       if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
       ntens = ntens + 1                    ! number of tensors stored after returning
       if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!      store the result tensors
       ierr=exatns_tensor_init(l1_res_tensor(itens))
       ierr=exatns_tensor_contract("L1R(i,a)+=L1(i,a)*K()",l1_res_tensor(itens),l1_tensor,one_diis_tensor)
       ierr=exatns_tensor_init(l2_res_tensor(itens))
       ierr=exatns_tensor_contract("L2R(i,j,a,b)+=L2(i,j,a,b)*K()",l2_res_tensor(itens),l2_tensor,one_diis_tensor)

!      calculate the error tensors by comparing the result of the t-equations with the tensors used as input
       ierr=exatns_tensor_init(l1_err_tensor(itens))
       ierr=exatns_tensor_contract("L1E(i,a)+=L1R(i,a)*K()",l1_err_tensor(itens),l1_tensor,one_diis_tensor)
       ierr=exatns_tensor_contract("L1E(i,a)+=L1I(i,a)*K()",l1_err_tensor(itens),l1_inp_tensor,one_diis_tensor, &
                                                                                                 prefactor=MINUS_ONE)
       ierr=exatns_tensor_init(l2_err_tensor(itens))
       ierr=exatns_tensor_contract("L2E(i,j,a,b)+=L2R(i,j,a,b)*K()",l2_err_tensor(itens),l2_tensor,one_diis_tensor)
       ierr=exatns_tensor_contract("L2E(i,j,a,b)+=L2I(i,j,a,b)*K()",l2_err_tensor(itens),l2_inp_tensor,one_diis_tensor, &
                                                                                              prefactor=MINUS_ONE)

!      update the b_diis array that holds the scalar products of the error tensors
       do i = 1, ntens
          ierr=exatns_tensor_init(scalar_tensor)
          ierr=exatns_tensor_contract("E()+=L+(i,a)*L(i,a)",scalar_tensor,l1_err_tensor(itens),l1_err_tensor(i))
          ierr=exatns_tensor_contract("E()+=L+(i,j,a,b)*L(i,j,a,b)",scalar_tensor,l2_err_tensor(itens),l2_err_tensor(i))
          ierr=exatns_tensor_get_scalar(scalar_tensor,scalar)
          b_diis(itens,i)=scalar
          b_diis(i,itens)=dconjg(scalar)
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
          ierr=exatns_tensor_init(l1_tensor)
          ierr=exatns_tensor_init(l2_tensor)
          do i = 1, ntens
             ierr=exatns_tensor_contract("L1(i,a)+=L1O(i,a)*K()",l1_tensor,l1_res_tensor(i), &
                                                                              one_diis_tensor,prefactor=coeff(i))
             ierr=exatns_tensor_contract("L2(i,j,a,b)+=L2O(i,j,a,b)*K()",l2_tensor,l2_res_tensor(i),one_diis_tensor, &
                                                                                              prefactor=coeff(i))
          end do
       end if

!      the extrapolated amplitudes are the new input tensors
       ierr=exatns_tensor_init(l1_inp_tensor)
       ierr=exatns_tensor_contract("L1O(i,a)+=L1(i,a)*K()",l1_inp_tensor,l1_tensor,one_diis_tensor)
       ierr=exatns_tensor_init(l2_inp_tensor)
       ierr=exatns_tensor_contract("L2O(i,j,a,b)+=L2(i,j,a,b)*K()",l2_inp_tensor,l2_tensor,one_diis_tensor)

       ierr=exatns_tensor_destroy(one_diis_tensor)

      end subroutine diis

#endif

end module exacorr_lambda
