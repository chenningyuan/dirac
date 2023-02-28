module talsh_eom_rspeq_solver
      use talsh_davidson
      use talsh
      use tensor_algebra
      use talsh_common_routines
      use exacorr_datatypes 
      use exacorr_utils
      use talsh_trial_eomccsd

      
      
      implicit none
      complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                         ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                         MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                         MINUS_TWO=(-2.D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0)

      type, public :: talsh_eom_results_tens_t
            character, pointer :: flavor(:)
            type(talsh_davidson_results_tens_t) :: value_vectors_t
      end type talsh_eom_results_tens_t

      private 

      public solve_eom_lr_right
      public solve_eom_lr_left
      public get_left_bvec

      contains

!----------------------------------------------------------------------------------------------------------------

      subroutine solve_eom_lr_right(exa_input, interm_t, int_t, interm_t_prop, &
                                    omega, tprop_tensor)
            use talsh_sigma_eomccsd
            use relcc_cfg
            !XY: here, we use relcc_cfg as done in talsh_drivers_eomccsd to take variable from *CCDIAG


      ! write(*,*) ""
      ! write(*,'(A)') "*******************************************************************************"
      ! write(*,'(A)') "***************************** CC-CI LR Right ******************************"
      ! write(*,'(A)') "*******************************************************************************"
      ! write(*,*) ""


!     input variables
            type(exacc_input), intent(in) :: exa_input
            type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
            type(talsh_intg_tens_t), intent(inout) :: int_t
            type(talsh_rsp_interm_tens_t) :: interm_t_prop
            complex(kind=8), intent(in)               :: omega

            type(talsh_eom_results_tens_t)  :: results_t
            type(talsh_comm_tens_t)       :: tprop_tensor

!     workspace
            type(talsh_davidson_workspace_tens_t) :: workspace_t
            type(talsh_projectors_tens_t) :: projector_t
!     error code
            integer :: ierr
            character(len=2) :: flavor = "EE"
            integer :: nroot                           ! initialize the number of trail vectors

            integer(INTD), dimension(2) :: vo_dims
            integer(INTD), dimension(4) :: vvoo_dims

            nroot = exa_input%N_eom_trail
            
            call talsh_davidson_setup(unit_output=6, &
                        convergence_threshold=relcc_mfd_convergence_threshold, &
                        solve_right=.true.,&
                        solve_left=.false.,&
                        solve_linear_sys=.true.,&
                        symmetric=.false., &
                        use_preconditioner=.false.,&
                        max_iterations=relcc_mfd_max_iterations,&
                        ! max_subspace_size=nroot,&
                        complex_mode=.true.)

            call print_date("Entering CC-CI Right Linear Systems calculation")

            call talsh_davidson_workspace_create(workspace_t, nroot, interm_t%nocc, interm_t%nvir, flavor, 2)
            ! call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
            ! call talsh_eomccsd_ee_trial_create(workspace_t)
            call talsh_eomccsd_LR_proper_create(workspace_t,interm_t_prop%epsilon_prop_vo)
            call talsh_eomccsd_LR_trial_create(workspace_t)
            

!     set pointers to sigma vector subroutines 
            interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ee

            vo_dims(1) = exa_input%nvir
            vo_dims(2) = exa_input%nocc
            vvoo_dims(1:2) = exa_input%nvir
            vvoo_dims(3:4) = exa_input%nocc

            ierr=talsh_tensor_construct(tprop_tensor%t1,C8,vo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(tprop_tensor%t2,C8,vvoo_dims,init_val=ZERO)


            call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, & 
                                    omega, interm_t_prop%epsilon_prop_vo, interm_t_prop%epsilon_prop_vvoo, tprop_tensor)

            call print_date("Finishing CC-CI Right Linear Systems calculation")

            call talsh_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_right


!----------------------------------------------------------------------------------------------------------------
      
      subroutine solve_eom_lr_left(exa_input, interm_t, int_t, interm_t_prop, omega, bvec_tensor, tprop_bar_tensor)
            use talsh_sigma_eomccsd
            use relcc_cfg

            ! write(*,*) ""
            ! write(*,'(A)') "*******************************************************************************"
            ! write(*,'(A)') "***************************** CC-CI LR Left ******************************"
            ! write(*,'(A)') "*******************************************************************************"
            ! write(*,*) ""
      
      
      !     input variables
                  type(exacc_input), intent(in) :: exa_input
                  type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
                  type(talsh_intg_tens_t), intent(inout) :: int_t
                  type(talsh_rsp_interm_tens_t) :: interm_t_prop
                  complex(kind=8),INTENT(IN)              :: omega
      

                  type(talsh_eom_results_tens_t)  :: results_t
                  type(talsh_comm_tens_t)       :: tprop_bar_tensor
                  type(talsh_comm_tens_t)       :: bvec_tensor
      !     workspace
                  type(talsh_davidson_workspace_tens_t) :: workspace_t
                  type(talsh_projectors_tens_t) :: projector_t
      !     error code
                  integer :: ierr
                  integer :: nroots = 1
                  character(len=2) :: flavor = "EE"

                  integer(INTD), dimension(2) :: ov_dims
                  integer(INTD), dimension(4) :: oovv_dims
                 
                  call talsh_davidson_setup(unit_output=6, &
                              convergence_threshold=relcc_mfd_convergence_threshold, &
                              solve_right=.false.,&
                              solve_left=.true.,&
                              solve_linear_sys=.true.,&
                              symmetric=.false., &
                              use_preconditioner=.false.,&
                              max_iterations=relcc_mfd_max_iterations,&
                              complex_mode=.true.)
      
                  call print_date("Entering CC-CI Left Linear Systems calculation")
      
                  call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2, 'L')
                  ! call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, 'L')
                  ! call talsh_eomccsd_ee_trial_create(workspace_t)
                  call talsh_eomccsd_LR_proper_create(workspace_t,bvec_tensor%t1)
                  call talsh_eomccsd_LR_trial_create(workspace_t)
      
      !     set pointers to sigma vector subroutines 
                  interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ee

                  ov_dims(1) = exa_input%nocc
                  ov_dims(2) = exa_input%nvir
                  oovv_dims(1:2) = exa_input%nocc
                  oovv_dims(3:4) = exa_input%nvir

      
                  ierr=talsh_tensor_construct(tprop_bar_tensor%t1,C8,ov_dims,init_val=ZERO)
                  ierr=talsh_tensor_construct(tprop_bar_tensor%t2,C8,oovv_dims,init_val=ZERO)

      
                  call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, &
                                          omega, bvec_tensor%t1, bvec_tensor%t2, tprop_bar_tensor)
                                          
                  call print_date("Finishing CC-CI Left Linear Systems calculation")
      
                  call talsh_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_left

!---------------------------------------------------------------------------------------------------------------------
      
      subroutine get_left_bvec(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,bvec_tensor)
            use tensor_algebra
            use talsh
            use exacorr_datatypes
            use exacorr_utils
            use talsh_common_routines
            use exacorr_global


            type(exacc_input), intent(in) :: exa_input
            type(talsh_intg_tens_t), intent(inout) :: int_t
            type(talsh_intermediates_tens_t),intent(inout) :: interm_t
            type(talsh_rsp_interm_tens_t) :: interm_t_prop
            type(talsh_tens_t), intent(inout)  :: l2_tensor
            type(talsh_tens_t), intent(inout)  :: l1_tensor
            ! type(talsh_comm_tens_t)       :: Tprop_tensor    !for CC-CI, we need right-vectors
            type(talsh_comm_tens_t)       :: bvec_tensor


            type(talsh_tens_t)       :: Aux_tensor1, one_tensor

            integer(INTD), dimension(4)   :: oovv_dims
            integer(INTD), dimension(2)   :: ov_dims
            integer(C_INT)              :: one_dims(1)
            complex(8), pointer         :: aux_tens(:)
            type(C_PTR):: body_p

            integer :: ierr


            one_dims(1) = 1
            oovv_dims(1:2) = exa_input%nocc
            oovv_dims(3:4) = exa_input%nvir
            ov_dims(1) = exa_input%nocc
            ov_dims(2) = exa_input%nvir


            ierr=talsh_tensor_construct(bvec_tensor%t1,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(bvec_tensor%t2,C8,oovv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(Aux_tensor1,C8,one_dims(1:0),init_val=ZERO)



            ierr=talsh_tensor_add('B(i,a)+=Y(i,a)',bvec_tensor%t1,interm_t_prop%eta_prop_ov,scale=ONE) 
            ierr=talsh_tensor_contract("B(i,a)+=l2(i,j,a,b)*E(b,j)",bvec_tensor%t1,l2_tensor,& 
                                                interm_t_prop%epsilon_prop_vo,scale=ONE)  

            ierr=talsh_tensor_add('B(i,j,a,b)+=Y(i,j,a,b)',bvec_tensor%t2,interm_t_prop%eta_prop_oovv,scale=ONE)  

            
!           EOM-CC
            ierr=talsh_tensor_contract("A1()+=l1(i,a)*E(a,i)",Aux_tensor1,l1_tensor,interm_t_prop%epsilon_prop_vo,scale=ONE)
            ierr=talsh_tensor_contract("A1()+=l2(i,j,a,b)*E(a,b,i,j)",Aux_tensor1,l2_tensor,&
                                          interm_t_prop%epsilon_prop_vvoo,scale=ONE_QUARTER_C)
 
            ierr=talsh_tensor_contract("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
            ierr=talsh_tensor_contract("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
                                                Aux_tensor1,scale=MINUS_ONE) 

            ierr=talsh_tensor_destruct(Aux_tensor1)



!--------------------------------------           CC-CI for doubel check    ----------------------------------------------------------
            ! call get_eta(int_t,interm_t,l1_tensor,l2_tensor,interm_t%eta_ov,interm_t%eta_oovv)
            ! ierr=talsh_tensor_contract("A1()+=l1(i,a)*T(a,i)",Aux_tensor1,l1_tensor,Tprop_tensor%t1,scale=ONE)
            ! ierr=talsh_tensor_contract("A1()+=l2(i,j,a,b)*T(a,b,i,j)",Aux_tensor1,l2_tensor,&
            !                               Tprop_tensor%t2,scale=ONE_QUARTER_C)
      

            ! ierr=talsh_tensor_contract("B(i,a)+=n(i,a)*A1()",bvec_tensor%t1,interm_t%eta_ov,Aux_tensor1,scale=ONE) 
            ! ierr=talsh_tensor_contract("B(i,j,a,b)+=n(i,j,a,b)*A1()",bvec_tensor%t2,interm_t%eta_oovv, &
            !                               Aux_tensor1,scale=ONE)  


            ! ierr=talsh_tensor_init(Aux_tensor1)
            ! ! call get_eta(int_t,interm_t,l1_tensor,l2_tensor,interm_t%eta_ov,interm_t%eta_oovv)

            ! ierr=talsh_tensor_contract("A1()+=n1(i,a)*T(a,i)",Aux_tensor1,interm_t%eta_ov,Tprop_tensor%t1,scale=MINUS_ONE)
            ! ierr=talsh_tensor_contract("A1()+=n2(i,j,a,b)*T(a,b,i,j)",Aux_tensor1,interm_t%eta_oovv,&
            !                               Tprop_tensor%t2,scale=MINUS_ONE_QUARTER_C)


            ! ierr=talsh_tensor_get_body_access(Aux_tensor1,body_p,C8,int(0,C_INT),DEV_HOST)
            ! call c_f_pointer(body_p,aux_tens,one_dims)
            ! print *, "Aux for lamda*epsilon:", aux_tens


            ! ierr=talsh_tensor_contract("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
            ! ierr=talsh_tensor_contract("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
            !                               Aux_tensor1,scale=MINUS_ONE)  

            ! ierr=talsh_tensor_destruct(Aux_tensor1)
            ! ierr=talsh_tensor_destruct(one_tensor)
!--------------------------------------           CC-CI for doubel check    ----------------------------------------------------------

      end subroutine get_left_bvec




      
!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_eta(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)
            ! Used in CC-CI left equation to construct eta vector

            !     input variables
                  type(talsh_intg_tens_t), intent(inout) :: int_t
                  type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
                  type(talsh_tens_t),    intent(inout) :: sigma1_tensor
                  type(talsh_tens_t),    intent(inout) :: sigma2_tensor
                  type(talsh_tens_t),    intent(inout) :: l1_tensor
                  type(talsh_tens_t),    intent(inout) :: l2_tensor
                  type(talsh_tens_t) :: one_tensor
            !     auxiliary tensor
                  type(talsh_tens_t) :: goo_tensor, gvv_tensor, oo_aux
                  integer            :: oo_dims(2), vv_dims(2)
                  integer(C_INT)     :: one_dims(1)
            !     error code
                  integer :: ierr
            
                  one_dims(1) = 1
                  oo_dims = (/interm_t%nocc,interm_t%nocc/)
                  vv_dims = (/interm_t%nvir,interm_t%nvir/)
                  ierr=talsh_tensor_construct(goo_tensor,C8,oo_dims,init_val=ZERO)
                  ierr=talsh_tensor_construct(gvv_tensor,C8,vv_dims,init_val=ZERO)
                  ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
                  ierr=talsh_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
                  ierr=talsh_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
            
            !***************
            !     l1_sigma *
            !***************
            
             ierr=talsh_tensor_init(sigma1_tensor)
            
                  
            !-----------------------------------------
            !     term 2: S1(i,a) += L(i,e) * F(e,a) |
            !-----------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",sigma1_tensor,l1_tensor,interm_t%fbar_vv)
            
            !-----------------------------------------
            !     term 3: S1(i,a) -= L(m,a) * F(i,m) |
            !-----------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",sigma1_tensor,l1_tensor,interm_t%fbar_oo,scale=MINUS_ONE)
            
            !-----------------------------------------------------
            !     term 4: S1(i,a) += 1/2 L(i,m,e,f) * W(e,f,a,m) |
            !-----------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",sigma1_tensor,l2_tensor,interm_t%w_vvvo,scale=ONE_HALF)
            
            !---------------------------------------------
            !     term 5: S1(i,a) -= G(f,e) * W(e,i,f,a) |
            !---------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",sigma1_tensor,gvv_tensor,interm_t%w_vovv)
            
            !---------------------------------------------
            !     term 6: S1(i,a) -= G(n,m) * W(m,i,n,a) |
            !---------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",sigma1_tensor,goo_tensor,interm_t%w_ooov)

                        
            !---------------------------------------------
            !     term 7: S1(i,a) += L(m,e) * W(i,e,a,m) |
            !---------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",sigma1_tensor,l1_tensor,interm_t%w_vovo,scale=MINUS_ONE)
            
            
            !-----------------------------------------------------
            !     term 8: S1(i,a) -= 1/2 L(m,n,a,e) * W(i,e,m,n) |
            !-----------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",sigma1_tensor,l2_tensor,interm_t%w_ovoo,scale=ONE_HALF)
            
            !***************
            !     l2_sigma *
            !***************
            
            ierr=talsh_tensor_init(sigma2_tensor)
            
            !--------------------------------------------------------
            !     term 2: S2(i,j,a,b) += P(a,b) L(i,j,a,e) * F(e,b) |
            !--------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
            
            
            !--------------------------------------------------------
            !     term 3: S2(i,j,a,b) -= P(i,j) L(i,m,a,b) * F(j,m) |
            !--------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
            
            !---------------------------------------------------------
            !     term 4: S2(i,j,a,b) += 1/2 L(m,n,a,b) * W(i,j,m,n) |
            !---------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",sigma2_tensor,l2_tensor,interm_t%w_oooo,scale=ONE_HALF)
            
            
            !-------------------------------------------------------------------
            !     term 5: S2(i,j,a,b) += P(i,j) P(a,b) L(i,m,a,e) * W(j,e,b,m) |
            !-------------------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            
            !--------------------------------------------------------
            !     term 6: S2(i,j,a,b) -= P(a,b) V(i,j,a,e) * G(e,b) |
            !--------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",sigma2_tensor,int_t%oovv,gvv_tensor)   !XY
            ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,b,e)*G(e,a)",sigma2_tensor,int_t%oovv,gvv_tensor,scale=MINUS_ONE)   !XY
            
            
            !-------------------------------------------------
            !     term 7: S2(i,j,a,b) -= P(a,b) L(m,a) * W(i,j,m,b) |
            !-------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(i,j,m,b)",sigma2_tensor,l1_tensor,interm_t%w_ooov,scale=MINUS_ONE)  !XY
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",sigma2_tensor,l1_tensor,interm_t%w_ooov)   !XY
            
            !--------------------------------------------------------
            !     term 8: S2(i,j,a,b) -= P(i,j) V(i,m,a,b) * G(j,m) |
            !--------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,i,a,b)*G(j,m)",sigma2_tensor,int_t%oovv,goo_tensor)   !XY
            ierr=talsh_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",sigma2_tensor,int_t%oovv,goo_tensor)   !XY
            
            !---------------------------------------------------------------------
            !     term 9: S2(i,j,a,b) -= P(i,j) V(m,j,a,b) * ( L(i,e) * T(e,m) ) |
            !---------------------------------------------------------------------
            
                  ! if (.not.interm_t%ccd) then
                  ! !   oo_dims = (/interm_t%nocc,interm_t%nocc/)
                  !   ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
                  ! !   ierr=talsh_tensor_contract("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1)
                  !   ierr=talsh_tensor_contract("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1,scale=MINUS_ONE)  !XY
                  !   ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,j,a,b)*U(i,m)",sigma2_tensor,int_t%oovv,oo_aux)
                  !   ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,m,a,b)*U(j,m)",sigma2_tensor,int_t%oovv,oo_aux)
                  !   ierr=talsh_tensor_destruct(oo_aux)
                  ! end if
            
            !---------------------------------------------------------
            !     term 10: S2(i,j,a,b) += P(i,j) L(i,e) * V(e,j,a,b) |
            !---------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",sigma2_tensor,l1_tensor,interm_t%w_vovv)   !XY
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",sigma2_tensor,l1_tensor,interm_t%w_vovv)   !XY
            
            
            !------------------------------------------------------------
            !     term 11: S2(i,j,a,b) += P(i,j) P(a,b) L(i,a) * F(j,b) |
            !------------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
            
            !--------------------------------------------------------
            !     term 12: S2(i,j,a,b) += 1/2 L(i,j,e,f) W(e,f,a,b) |
            !--------------------------------------------------------
            
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",sigma2_tensor,l2_tensor,interm_t%w_vvvv,scale=ONE_HALF)

            
            ierr=talsh_tensor_destruct(goo_tensor)
            ierr=talsh_tensor_destruct(gvv_tensor)
            
            
      end subroutine get_eta

                
end module
