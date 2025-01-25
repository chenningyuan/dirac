module talsh_eom_rspeq_solver
      use talsh_davidson
      use talsh
      use tensor_algebra
      use talsh_common_routines
      use exacorr_datatypes 
      use exacorr_utils
      use talsh_trial_eomccsd
      use talsh_eom

      
      
      implicit none
      complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                         ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                         MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                         MINUS_TWO=(-2.D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0)

      private 

      public solve_eom_lr_right
      public solve_eom_lr_left 
      public get_left_bvec
      public solve_cc_lr_left, get_left_bvec_ccci, get_left_bvec_cccc 

      contains

!----------------------------------------------------------------------------------------------------------------

      subroutine solve_eom_lr_right(exa_input, interm_t, int_t, interm_t_prop, &
                                    omega, tprop_tensor)
            use talsh_sigma_eomccsd
            use mfsolver_dirac_input

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

            nroot = exa_input%N_eom_trial
            
            call talsh_davidson_setup(unit_output=6, &
                        convergence_threshold=mfsolver_options_convergence_threshold, &
                        max_subspace_size=mfsolver_options_max_subspace_size,&
                        verbose=mfsolver_options_verbose,&
                        solve_right=.true.,&
                        solve_left=.false.,&
                        solve_linear_sys=.true.,&
                        symmetric=.false., &
                        use_preconditioner=.true.,&
                        max_iterations=mfsolver_options_max_iterations,&
                        ! max_subspace_size=nroot,&
                        complex_mode=.true.)

            call print_date("Entering CC-CI Right Linear Systems calculation")

            call talsh_davidson_workspace_create(workspace_t, nroot, interm_t%nocc, interm_t%nvir, flavor, 2)
            call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
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


            call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, flavor,& 
                                    omega, interm_t_prop%epsilon_prop_vo, interm_t_prop%epsilon_prop_vvoo, tprop_tensor)

            call print_date("Finishing CC-CI Right Linear Systems calculation")

            call talsh_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_right


!----------------------------------------------------------------------------------------------------------------
      
      subroutine solve_eom_lr_left(exa_input, interm_t, int_t, interm_t_prop, omega, bvec_tensor, tprop_bar_tensor)
            use talsh_sigma_eomccsd
            use mfsolver_dirac_input

            ! write(*,*) ""
            ! write(*,'(A)') "*******************************************************************************"
            ! write(*,'(A)') "***************************** EOM LR Left ******************************"
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
                              convergence_threshold=mfsolver_options_convergence_threshold, &
                              max_subspace_size=mfsolver_options_max_subspace_size,&
                              verbose=mfsolver_options_verbose,&
                              solve_right=.false.,&
                              solve_left=.true.,&
                              solve_linear_sys=.true.,&
                              symmetric=.false., &
                              use_preconditioner=.true.,&
                              max_iterations=mfsolver_options_max_iterations,&
                              complex_mode=.true.)
      
                  call print_date("Entering CC-CI Left Linear Systems calculation")
      
                  call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2, 'L')
                  call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, 'L')
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

      
                  call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, flavor,&
                                          -omega, bvec_tensor%t1, bvec_tensor%t2, tprop_bar_tensor)   !left equation use -omega
                                          
                  call print_date("Finishing CC-CI Left Linear Systems calculation")
      
                  call talsh_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_left

!---------------------------------------------------------------------------------------------------------------------
      
      subroutine solve_cc_lr_left(exa_input, interm_t, int_t, interm_t_prop, omega, bvec_tensor, tprop_bar_tensor)
            use talsh_sigma_eomccsd
            use mfsolver_dirac_input

            ! write(*,*) ""
            ! write(*,'(A)') "*******************************************************************************"
            ! write(*,'(A)') "***************************** CC-CC LR Left ******************************"
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
                              convergence_threshold=mfsolver_options_convergence_threshold, &
                              max_subspace_size=mfsolver_options_max_subspace_size,&
                              verbose=mfsolver_options_verbose,&
                              solve_right=.false.,&
                              solve_left=.true.,&
                              solve_linear_sys=.true.,&
                              symmetric=.false., &
                              use_preconditioner=.true.,&
                              max_iterations=mfsolver_options_max_iterations,&
                              complex_mode=.true.)
      
                  call print_date("Entering CC-CC Left Linear Systems calculation")
      
                  call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2, 'L')
                  call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, 'L')
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

      
                  call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, flavor,&
                                          -omega, bvec_tensor%t1, bvec_tensor%t2, tprop_bar_tensor)   !left equation use -omega
                                          
                  call print_date("Finishing CC-CC Left Linear Systems calculation")
      
                  call talsh_davidson_workspace_destroy(workspace_t)

      end subroutine solve_cc_lr_left

!---------------------------------------------------------------------------------------------------------------------
      
      subroutine get_left_bvec(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,bvec_tensor)

      !XY: Although the final QR result of EOM is same as CC-CI, the response equations and 
      ! the left perturbed ampiltudes tbar are not same. The point is the left bevctor EOM used is not the one of CC-CI 
      ! CC-CI left bvector is depended on the right perturbed ampiltudes. So, here we choose the EOM type and comment out the CC-CI as backup
      ! see more details: J. Chem Phys, 142, p.114109, 2015, & J. Chem Phys, 144, p.024102, 2016

            use tensor_algebra
            use talsh
            use exacorr_datatypes
            use exacorr_utils
            use talsh_common_routines
            use exacorr_global


            type(exacc_input), intent(in) :: exa_input
            type(talsh_intg_tens_t), intent(inout) :: int_t
            type(talsh_intermediates_tens_t),intent(inout) :: interm_t
            type(talsh_rsp_interm_tens_t),intent(inout) :: interm_t_prop
            type(talsh_tens_t), intent(inout)  :: l2_tensor
            type(talsh_tens_t), intent(inout)  :: l1_tensor
            type(talsh_comm_tens_t), intent(inout)       :: bvec_tensor


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


      end subroutine get_left_bvec



!-----------------------------------------------------------------------------------------------------------------------------
      
      subroutine get_left_bvec_ccci(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,Tprop_tensor,bvec_tensor)

            !XY: the ccci left bvector in the left response equations. 
            !    This routine is just for backup since we implemented EOM for quadratic response
      
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
                  type(talsh_comm_tens_t), intent(inout)       :: Tprop_tensor    !for CC-CI, we need right-vectors
                  type(talsh_comm_tens_t), intent(inout)       :: bvec_tensor
      
      
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
      
                  
                  ierr=talsh_tensor_contract("A1()+=l1(i,a)*E(a,i)",Aux_tensor1,l1_tensor,interm_t_prop%epsilon_prop_vo,scale=ONE)
                  ierr=talsh_tensor_contract("A1()+=l2(i,j,a,b)*E(a,b,i,j)",Aux_tensor1,l2_tensor,&
                                                interm_t_prop%epsilon_prop_vvoo,scale=ONE_QUARTER_C)
       
                  ierr=talsh_tensor_contract("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
                  ierr=talsh_tensor_contract("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
                                                      Aux_tensor1,scale=MINUS_ONE) 
      
                  ierr=talsh_tensor_destruct(Aux_tensor1)
      

                  call get_eta(int_t,interm_t,l1_tensor,l2_tensor,interm_t%eta_ov,interm_t%eta_oovv)
                  ierr=talsh_tensor_contract("A1()+=l1(i,a)*T(a,i)",Aux_tensor1,l1_tensor,Tprop_tensor%t1,scale=ONE)
                  ierr=talsh_tensor_contract("A1()+=l2(i,j,a,b)*T(a,b,i,j)",Aux_tensor1,l2_tensor,&
                                                Tprop_tensor%t2,scale=ONE_QUARTER_C)
            
      
                  ierr=talsh_tensor_contract("B(i,a)+=n(i,a)*A1()",bvec_tensor%t1,interm_t%eta_ov,Aux_tensor1,scale=ONE) 
                  ierr=talsh_tensor_contract("B(i,j,a,b)+=n(i,j,a,b)*A1()",bvec_tensor%t2,interm_t%eta_oovv, &
                                                Aux_tensor1,scale=ONE)  
      
      
                  ierr=talsh_tensor_init(Aux_tensor1)
                  ! call get_eta(int_t,interm_t,l1_tensor,l2_tensor,interm_t%eta_ov,interm_t%eta_oovv)
      
                  ierr=talsh_tensor_contract("A1()+=n1(i,a)*T(a,i)",Aux_tensor1,interm_t%eta_ov,Tprop_tensor%t1,scale=MINUS_ONE)
                  ierr=talsh_tensor_contract("A1()+=n2(i,j,a,b)*T(a,b,i,j)",Aux_tensor1,interm_t%eta_oovv,&
                                                Tprop_tensor%t2,scale=MINUS_ONE_QUARTER_C)
      
      
                  ierr=talsh_tensor_get_body_access(Aux_tensor1,body_p,C8,int(0,C_INT),DEV_HOST)
      
      
                  ierr=talsh_tensor_contract("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
                  ierr=talsh_tensor_contract("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
                                                Aux_tensor1,scale=MINUS_ONE)  
      
                  ierr=talsh_tensor_destruct(Aux_tensor1)
                  ierr=talsh_tensor_destruct(one_tensor)
      
            end subroutine get_left_bvec_ccci

!-----------------------------------------------------------------------------------------------------------------------------
      
            subroutine get_left_bvec_cccc(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,Tprop_tensor,bvec_tensor)

                  !The cc-cc left bvector in the left response equations. 
                  !This routine can be used in cc-cc quadratic response function 
                  !and asymmetric cc-cc linear response function
                  !Written by Xiang Yuan Janunary 2024, Lille, France

            
                  use tensor_algebra
                  use talsh
                  use exacorr_datatypes
                  use exacorr_utils
                  use talsh_common_routines
                  use exacorr_global
                  use talsh_eom_LR
            
            
                  type(exacc_input), intent(in) :: exa_input
                  type(talsh_intg_tens_t), intent(inout) :: int_t
                  type(talsh_intermediates_tens_t),intent(inout) :: interm_t
                  type(talsh_rsp_interm_tens_t) :: interm_t_prop
                  type(talsh_tens_t), intent(inout)  :: l2_tensor
                  type(talsh_tens_t), intent(inout)  :: l1_tensor
                  type(talsh_comm_tens_t)       :: Tprop_tensor    !for CC-CI, we need right-vectors
                  type(talsh_comm_tens_t)       :: bvec_tensor
            
            
                  type(talsh_tens_t)       :: Aux_tensor1, one_tensor
                  type(talsh_comm_tens_t)  :: B  !Ft^x
            
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
                  call get_tF(B%t1,B%t2,Tprop_tensor%t1,Tprop_tensor%t2,l1_tensor,l2_tensor,&
                        interm_t,int_t,exa_input%nocc,exa_input%nvir)
            
                  ierr=talsh_tensor_add('B(i,a)+=Y(i,a)',bvec_tensor%t1,interm_t_prop%eta_prop_ov,scale=ONE)             
                  ierr=talsh_tensor_add('B(i,j,a,b)+=Y(i,j,a,b)',bvec_tensor%t2,interm_t_prop%eta_prop_oovv,scale=ONE)  

                  ierr=talsh_tensor_add('B(i,a)+=Ft(i,a)',bvec_tensor%t1,B%t1,scale=ONE)             
                  ierr=talsh_tensor_add('B(i,j,a,b)+=Ft(i,j,a,b)',bvec_tensor%t2,B%t2,scale=ONE) 
            
      end subroutine get_left_bvec_cccc

!------------------------------------------------------------------------------------------------------------------------------------

      subroutine get_eta(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)
            ! Used in CC-CI left equation to construct eta vector. In EOM-QR, we don't need this term.

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
                  integer            :: oo_dims(2), vv_dims(2), ov_dims(2)
                  integer            :: oovv_dims(4)
                  integer(C_INT)     :: one_dims(1)
            !     error code
                  integer :: ierr
            
                  one_dims(1) = 1
                  oo_dims = (/interm_t%nocc,interm_t%nocc/)
                  vv_dims = (/interm_t%nvir,interm_t%nvir/)
                  ov_dims = (/interm_t%nocc,interm_t%nvir/)
                  oovv_dims = (/interm_t%nocc,interm_t%nocc,interm_t%nvir,interm_t%nvir/)
                  ierr=talsh_tensor_construct(goo_tensor,C8,oo_dims,init_val=ZERO)
                  ierr=talsh_tensor_construct(gvv_tensor,C8,vv_dims,init_val=ZERO)
                  ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
                  ierr=talsh_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
                  ierr=talsh_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
            
            !***************
            !     l1_sigma *
            !***************
            ierr=talsh_tensor_construct(sigma2_tensor,C8,oovv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(sigma1_tensor,C8,ov_dims,init_val=ZERO)

            !  ierr=talsh_tensor_init(sigma1_tensor)
            
                  
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
            
            ! ierr=talsh_tensor_init(sigma2_tensor)
            
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

      !-------------------------------------------------------------------------------------------------------------------
!! LH start Ft^x

subroutine get_tF(B1,B2,tX1,tX2,l1_tensor,l2_tensor,interm_t,int_t,nocc,nvir)

      type(talsh_tens_t) :: B1, B2
      type(talsh_tens_t), intent(inout) :: tX1, tX2
      type(talsh_tens_t), intent(inout) :: l1_tensor,l2_tensor
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t 
      type(talsh_intg_tens_t), intent(inout) :: int_t
      !Tensor dimensions
      integer, intent(in) :: nocc, nvir
      integer(INTD), dimension(2) :: oo_dims, vv_dims, ov_dims, vo_dims
      integer(INTD), dimension(4) :: vvoo_dims, ovoo_dims, ovov_dims
      integer(INTD), dimension(4) :: vovo_dims, vovv_dims, oovo_dims
      integer(INTD), dimension(4) :: vooo_dims, ooov_dims, oovv_dims
      integer(INTD), dimension(4) :: oooo_dims, vvvv_dims, ovvo_dims
      integer(C_INT)              :: one_dims(1)
  
      integer(INTD)                   :: ierr
      logical :: debug_local = .false.
      
  
  
      type(talsh_tens_t) :: Aux_oo, Aux_vv, Aux_vo, Aux_ov
      type(talsh_tens_t) :: Aux_oovv, Aux_ooov, Aux_vooo
      type(talsh_tens_t) :: Aux_oooo, Aux_vvvv, Aux_ovov
      type(talsh_tens_t) :: Aux_vovv, Aux_vovo, Aux_oovo
      type(talsh_tens_t) :: Aux_ovvo
      type(talsh_tens_t) :: Aux_one
  
      oo_dims(1:2) = nocc
      vv_dims(1:2) = nvir
      
      ov_dims(1)  = nocc
      ov_dims(2)  = nvir
  
      vo_dims(1)  = nvir
      vo_dims(2)  = nocc
  
      vvvv_dims(1:4) = nvir
      oooo_dims(1:4) = nocc
  
      vooo_dims(1:4) = nocc
      vooo_dims(1) = nvir        
  
      ooov_dims(1:4) =nocc
      ooov_dims(4) = nvir
  
      oovv_dims(1:2) = nocc
      oovv_dims(3:4) = nvir
  
      ovov_dims(1:4) = nocc
      ovov_dims(2)   = nvir
      ovov_dims(4)   = nvir
  
      vovo_dims(1:4) = nvir
      vovo_dims(2) = nocc
      vovo_dims(4) = nocc
  
      vovv_dims(1:4) = nvir
      vovv_dims(2) = nocc
      
      oovo_dims(1:4) = nocc
      oovo_dims(3)   = nvir
  
      ovvo_dims(1:4) = nocc
      ovvo_dims(2:3) = nvir
  
      one_dims(1) = 1
  
      ierr=talsh_tensor_construct(B1,C8,ov_dims,init_val=ZERO)        
      ierr=talsh_tensor_construct(B2,C8,oovv_dims,init_val=ZERO)        
  
      ierr=talsh_tensor_construct(Aux_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_vo,C8,vo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_ov,C8,ov_dims,init_val=ZERO)
      
      ierr=talsh_tensor_construct(Aux_vvvv,C8,vvvv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_oovv,C8,oovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_vooo,C8,vooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_ooov,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_ovov,C8,ovov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_vovv,C8,vovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_oovv,C8,oovv_dims,init_val=ZERO)
  
      ierr=talsh_tensor_construct(Aux_oovo,C8,oovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Aux_ovvo,C8,ovvo_dims,init_val=ZERO)
  
      ierr=talsh_tensor_construct(Aux_one,C8,one_dims(1:0),init_val=ZERO)
  
  
  ! l=0  ; t=1  ; h=-2
      ! U1 = + (\ts ck ) \ts ai \W kica       -> V oovv
      ierr=talsh_tensor_contract("B(k,c)+=T(a,i)*W(k,i,c,a)",B1,tX1,int_t%oovv)
      if (debug_local) print*, "ierr U1", ierr
  
  
  ! l=-1 ; t=1 ; h=-1
      ! c7 = - (\ts ck) \ts ai \lambdas ic \F ka
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,i)+=T(a,i)*F(k,a)",Aux_oo,tX1,interm_t%Fbar_ov)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,i)*L(i,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE)        
      if (debug_local) print*, "ierr c7", ierr
  
      ! c8 = - (\ts ck) \ts ai \lambdas ka \F ic        
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,i)+=L(k,a)*T(a,i)",Aux_oo,l1_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,i)*F(i,c)",B1,Aux_oo,interm_t%Fbar_ov,scale=MINUS_ONE)
      if (debug_local) print*, "ierr c8", ierr
  
      ! g2 = + (\ts ck) \ts ai \lambdas ck \W eica (vovv)
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(e,c)+=T(a,i)*W(e,i,c,a)",Aux_vv,tX1,interm_t%w_vovv)
      if (debug_local) print*, "ierr g2_1st", ierr
      ierr=talsh_tensor_contract("B(k,c)+=A(e,c)*L(k,e)",B1,Aux_vv,l1_tensor)
      if (debug_local) print*, "ierr g2_2nd", ierr
  
      ! k2 = - (\ts ck) \ts ai \lambdas mc \W kima (ooov)
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,m)+=T(a,i)*W(k,i,m,a)",Aux_oo,tX1,interm_t%w_ooov)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,m)*L(m,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE)
      if (debug_local) print*, "ierr k2", ierr
  
      ! k3 = - (\ts ck) \ts ai \lambdas ma \W ikmc (ooov)
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(m,i)+=L(m,a)*T(a,i)",Aux_oo,l1_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,c)+=A(m,i)*W(i,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(m,i)*W(i,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE_HALF)  !XY
      if (debug_local) print*, "ierr k3", ierr
  
      ! g3 = + (\ts ck) \ts ai \lambda ie \W ekac (vovv)
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(a,e)+=L(i,e)*T(a,i)",Aux_vv,l1_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,c)+=A(a,e)*W(e,k,a,c)",B1,Aux_vv,interm_t%w_vovv)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(a,e)*W(e,k,a,c)",B1,Aux_vv,interm_t%w_vovv,scale=ONE_HALF)    !XY
      if (debug_local) print*, "ierr g3", ierr
  
  
  ! l=-1; t=2 ; h=-2
      ! l5 = - one_half (\ts ck) \td abij \lambdas ic \W kjab (oovv) -> V oovv
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,i)+=T(a,b,i,j)*W(k,j,a,b)",Aux_oo,tX2,int_t%oovv)
      ! ierr=talsh_tensor_contract("B(c,k)+=A(k,i)*L(i,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE_HALF)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,i)*L(i,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE_HALF)   !XY
      if (debug_local) print*, "ierr l5_2", ierr
  
      ! l11 = - one_half (\ts ck) \td abij \lambdas ka \W ijcb (oovv) -> V oovv
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(a,c)+=T(a,b,i,j)*W(i,j,c,b)",Aux_vv,tX2,int_t%oovv)
      ! ierr=talsh_tensor_contract("B(c,k)+=A(a,c)*L(k,a)",B1,Aux_vv,l1_tensor,scale=MINUS_ONE_HALF)
      ierr=talsh_tensor_contract("B(k,c)+=A(a,c)*L(k,a)",B1,Aux_vv,l1_tensor,scale=MINUS_ONE_HALF)  !XY
      if (debug_local) print*, "ierr l11", ierr
  
      ! l18 = + (\ts ck) \td abij \lambdas ia \W jkbc (oovv) -> V oovv
      ierr=talsh_tensor_init(Aux_vo)
      ierr=talsh_tensor_contract("A(b,j)+=T(a,b,i,j)*L(i,a)",Aux_vo,tX2,l1_tensor)
      ! ierr=talsh_tensor_contract("B(c,k)+=A(b,j)*W(j,k,b,c)",B1,Aux_vo,int_t%oovv)
      ierr=talsh_tensor_contract("B(k,c)+=A(b,j)*W(j,k,b,c)",B1,Aux_vo,int_t%oovv)   !XY
      if (debug_local) print*, "ierr l18", ierr
  
  
  ! l=-2; t=1 ; h=0
      ! f3 = - (\ts ck) \ts ai \lambdad mkae \W iemc (ovov)
      ierr=talsh_tensor_init(Aux_ooov)
      if (debug_local) print *, "ierr init f3", ierr
      ierr=talsh_tensor_contract("A(m,k,i,e)+=L(m,k,a,e)*T(a,i)",Aux_ooov,l2_tensor,tX1)
      if (debug_local) print *, "ierr get Aux_ooov", ierr
      ierr=talsh_tensor_contract("B(k,c)+=A(m,k,i,e)*W(e,i,c,m)",B1,Aux_ooov,interm_t%w_vovo,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(m,k,i,e)*W(e,i,c,m)",B1,Aux_ooov,interm_t%w_vovo,scale=MINUS_ONE_HALF) !XY
      if (debug_local) print*, "ierr f3", ierr
  
      ! f4 = - (\ts ck) \ts ai \lambdad imec \W ekam (vovo)
      ierr=talsh_tensor_init(Aux_vooo)
      ierr=talsh_tensor_contract("A(e,k,i,m)+=T(a,i)*W(e,k,a,m)",Aux_vooo,tX1,interm_t%w_vovo)
      ierr=talsh_tensor_contract("B(k,c)+=A(e,k,i,m)*L(i,m,e,c)",B1,Aux_vooo,l2_tensor,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(e,k,i,m)*L(i,m,e,c)",B1,Aux_vooo,l2_tensor,scale=MINUS_ONE_HALF)     !XY
      if (debug_local) print*, "ierr f4", ierr
  
      ! i2 = +one_half  (\ts ck) \ts ai \lambdad kief \W efca (vvvv)
      ierr=talsh_tensor_init(Aux_oovv)
      ierr=talsh_tensor_contract("A(k,i,c,a)+=L(k,i,e,f)*W(e,f,c,a)",Aux_oovv,l2_tensor,interm_t%w_vvvv)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,i,c,a)*T(a,i)",B1,Aux_oovv,tX1,scale=ONE_HALF)
      if (debug_local) print*, "ierr i2", ierr
  
      ! n2= + one_half (\ts ck) \ts ai  \lambda mnca \W kimn (oooo)
      ierr=talsh_tensor_init(Aux_oovv)
      ierr=talsh_tensor_contract("A(k,i,c,a)+=L(m,n,c,a)*W(k,i,m,n)",Aux_oovv,l2_tensor,interm_t%w_oooo)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,i,c,a)*T(a,i)",B1,Aux_oovv,tX1,scale=ONE_HALF)
      if (debug_local) print*, "ierr n2", ierr
  
  ! l=-2; t=2 ; h=-1
      ! c9 = -one_half (\ts ck) \td abij \lambdad ijac \F kb (ov)
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(b,c)+=L(i,j,a,c)*T(a,b,i,j)",Aux_vv,l2_tensor,tX2)
      ierr=talsh_tensor_contract("B(k,c)+=A(b,c)*F(k,b)",B1,Aux_vv,interm_t%Fbar_ov,scale=MINUS_ONE_HALF)
      if (debug_local) print*, "ierr c9", ierr
  
      ! c10 = - one_half (\ts ck) \td abij \lambdad ikab \F jc (ov)
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,j)+=L(i,k,a,b)*T(a,b,i,j)",Aux_oo,l2_tensor,tX2)
      ierr=talsh_tensor_contract("B(k,c)+=A(k,j)*F(j,c)",B1,Aux_oo,interm_t%Fbar_ov,scale=MINUS_ONE_HALF)
      if (debug_local) print*, "ierr c10", ierr
  
      ! g6 = - (\ts ck) \td abij \lambdad ikeb \W ejac (vovv)
      ierr=talsh_tensor_init(Aux_vovo)
      ierr=talsh_tensor_contract("A(a,k,e,j)+=L(i,k,e,b)*T(a,b,i,j)",Aux_vovo,l2_tensor,tX2)
      ierr=talsh_tensor_contract("B(k,c)+=A(a,k,e,j)*W(e,j,a,c)",B1,Aux_vovo,interm_t%w_vovv,scale=MINUS_ONE) 
      ! ierr=talsh_tensor_contract("B(k,c)+=A(a,k,e,j)*W(e,j,a,c)",B1,Aux_vovo,interm_t%w_vovv,scale=MINUS_ONE_HALF)  !XY
      if (debug_local) print*, "ierr g6", ierr
  
  
      ! k6 = + (\ts ck) \td abij \lambdad mjac \W ikmb (ooov)
      ierr=talsh_tensor_init(Aux_ovov)
      if (debug_local) print*, "ierr k6_1", ierr
      ! ierr=talsh_tensor_contract("A(m,b,i,c)+=L(i,k,e,b)*T(a,b,i,j)",Aux_ovov,l2_tensor,tX2)
      ierr=talsh_tensor_contract("A(m,b,i,c)+=L(m,j,a,c)*T(a,b,i,j)",Aux_ovov,l2_tensor,tX2)     !XY
      if (debug_local) print*, "ierr k6_2", ierr
      ierr=talsh_tensor_contract("B(k,c)+=A(m,b,i,c)*W(i,k,m,b)",B1,Aux_ovov,interm_t%w_ooov)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(m,b,i,c)*W(i,k,m,b)",B1,Aux_ovov,interm_t%w_ooov,scale=ONE_HALF)  !XY
      if (debug_local) print*, "ierr k6", ierr
  
      ! g7 = - one_quarter (\ts ck) \td abij \lambdad ijec \W ekab (vovv)
      ierr=talsh_tensor_init(Aux_vooo)
      ierr=talsh_tensor_contract("A(e,k,i,j)+=T(a,b,i,j)*W(e,k,a,b)",Aux_vooo,tX2,interm_t%w_vovv)
      ierr=talsh_tensor_contract("B(k,c)+=A(e,k,i,j)*L(i,j,e,c)",B1,Aux_vooo,l2_tensor,scale=MINUS_ONE_QUARTER_C)
      if (debug_local) print*, "ierr g7", ierr
  
      !k7 = + one_quarter (\ts ck) \td abij \lambdad mkab \W ijmc (ooov)
      ierr=talsh_tensor_init(Aux_oooo)
      ierr=talsh_tensor_contract("A(m,k,i,j)+=L(m,k,a,b)*T(a,b,i,j)",Aux_oooo,l2_tensor,tX2)
      ierr=talsh_tensor_contract("B(k,c)+=A(m,k,i,j)*W(i,j,m,c)",B1,Aux_oooo,interm_t%w_ooov,scale=ONE_QUARTER_C)
      if (debug_local) print*, "ierr k7", ierr
  
      !g9 = +one_half (\ts ck) \td abij  \lambdad ijae \W ekbc (vovv)
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(b,e)+=L(i,j,a,e)*T(a,b,i,j)",Aux_vv,l2_tensor,tX2)
      ierr=talsh_tensor_contract("B(k,c)+=A(b,e)*W(e,k,b,c)",B1,Aux_vv,interm_t%w_vovv,scale=ONE_HALF)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(b,e)*W(e,k,b,c)",B1,Aux_vv,interm_t%w_vovv,scale=ONE_QUARTER_C)   !XY
      if (debug_local) print*, "ierr g9", ierr
  
      !k9 = - one_half (\ts ck) \td abij  \lambdad imab \W jkmc (ooov)
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(m,j)+=L(i,m,a,b)*T(a,b,i,j)",Aux_oo,l2_Tensor,tX2)
      ierr=talsh_tensor_contract("B(k,c)+=A(m,j)*W(j,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE_HALF)
      ! ierr=talsh_tensor_contract("B(k,c)+=A(m,j)*W(j,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE_QUARTER_C)  !XY
      if (debug_local) print*, "ierr k9", ierr
  
  
  ! !********************Additional******************************************************************************
  
  ! eq l22 bis :: \B kc = -1/2 \ld mnec \V kifa \td efmn \tXs ai                             
      ierr=talsh_tensor_init(Aux_vv)                                                             
      ierr=talsh_tensor_contract("A(f,c)+=L(m,n,e,c)*T(e,f,m,n)",Aux_vv,l2_tensor,interm_t%t2)     
      if(debug_local) write(*,*) "ierr l22 bis  =",ierr                                     
      ierr=talsh_tensor_init(Aux_oovv)                                                           
      ierr=talsh_tensor_contract("B(k,i,c,a)+=A(f,c)*V(k,i,f,a)",Aux_oovv,Aux_vv,int_t%oovv)     
      if(debug_local) write(*,*) "ierr l22 bis =",ierr                                      
      ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                            B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                
      if(debug_local) write(*,*) "ierr l22 bis=",ierr                                       
                                                                                                 
  ! eq l23 bis :: \B kc = -1/2 \ld mkef \V nica \td efmn \tXs ai                             
      ierr=talsh_tensor_init(Aux_oo)                                                             
      ierr=talsh_tensor_contract("A(k,n)+=L(m,k,e,f)*T(e,f,m,n)",Aux_oo,l2_tensor,interm_t%t2)     
      if(debug_local) write(*,*) "ierr l23 bis  =",ierr                                     
      ierr=talsh_tensor_init(Aux_oovv)                                                           
      ierr=talsh_tensor_contract("B(k,i,c,a)+=A(k,n)*V(n,i,c,a)",Aux_oovv,Aux_oo,int_t%oovv)     
      if(debug_local) write(*,*) "ierr l23 bis =",ierr                                      
      ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                            B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                
      if(debug_local) write(*,*) "ierr l23 bis=",ierr                                       
                                                                                                 
  ! eq l32 bis :: \B kc = -1/2 \ld mnea \V ikfc \td efmn \tXs ai                             
      ierr=talsh_tensor_init(Aux_vv)                                                             
      ierr=talsh_tensor_contract("A(f,a)+=L(m,n,e,a)*T(e,f,m,n)",Aux_vv,l2_tensor,interm_t%t2)     
      if(debug_local) write(*,*) "ierr l32 bis  =",ierr                                     
      ierr=talsh_tensor_init(Aux_oovv)                                                           
      ierr=talsh_tensor_contract("B(k,i,c,a)+=A(f,a)*V(i,k,f,c)",Aux_oovv,Aux_vv,int_t%oovv)     
      if(debug_local) write(*,*) "ierr l32 bis =",ierr                                      
      ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                            B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                
      if(debug_local) write(*,*) "ierr l32 bis=",ierr                                       
                                                                                                 
  ! eq l34 bis :: \B kc = -1/2 \ld mief \V nkac \td efmn \tXs ai                             
      ierr=talsh_tensor_init(Aux_oo)                                                             
      ierr=talsh_tensor_contract("A(i,n)+=L(m,i,e,f)*T(e,f,m,n)",Aux_oo,l2_tensor,interm_t%t2)     
      if(debug_local) write(*,*) "ierr l34 bis  =",ierr                                     
      ierr=talsh_tensor_init(Aux_oovv)                                                           
      ierr=talsh_tensor_contract("B(k,i,c,a)+=A(i,n)*V(n,k,a,c)",Aux_oovv,Aux_oo,int_t%oovv)     
      if(debug_local) write(*,*) "ierr l34 bis =",ierr                                      
      ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                            B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                    
      if(debug_local) write(*,*) "ierr l34 bis=",ierr                                       
                                                                                                 
  ! !************************************************************************************************************
  
  
  ! l=-1; t=1 -> h=-2 
      ! l3  = + \lambdas kc (\tYd || \Pmm cd \Pmm kl)  cdkl  \tXs ai  \W lida (oovv)
      ierr=talsh_tensor_init(Aux_ov)
      ierr=talsh_tensor_contract("A(l,d)+=T(a,i)*V(l,i,d,a)",Aux_ov,tX1,int_t%oovv)
      if (debug_local) print*, "ierr l3_1", ierr
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,d)*L(k,c)",B2,Aux_ov,l1_tensor)
      if (debug_local) print*, "ierr l3_2", ierr
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,d)*L(l,c)",B2,Aux_ov,l1_tensor,scale=MINUS_ONE)
      if (debug_local) print*, "ierr l3_3", ierr
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,c)*L(k,d)",B2,Aux_ov,l1_tensor,scale=MINUS_ONE)
      if (debug_local) print*, "ierr l3_4", ierr
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,c)*L(l,d)",B2,Aux_ov,l1_tensor)
      if (debug_local) print*, "ierr l3", ierr
  
  
      ! l16 = - \lambdas ka ( \tYd cdkl || \Pmm kl | one_half )  \tXs ai \W ilcd 
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,i)+=L(k,a)*T(a,i)",Aux_oo,l1_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,i)*V(i,l,c,d)",B2,Aux_oo,int_t%oovv,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,i)*V(i,k,c,d)",B2,Aux_oo,int_t%oovv)
      if (debug_local) print*, "ierr l16", ierr
  
      ! l17 = - \lambdas ic (\tYd cdkl || \Pmm cd | one_half ) \tXs ai \W klad
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(a,c)+=L(i,c)*T(a,i)",Aux_vv,l1_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,c)*V(k,l,a,d)",B2,Aux_vv,int_t%oovv,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,d)*V(k,l,a,c)",B2,Aux_vv,int_t%oovv)
      if (debug_local) print*, "ierr l17", ierr
  
  
  ! l=-2; t=1 -> h=-1
  
      ! c11 = -  \lambdad kicd ( \tYd cdkl || \Pmm kl | one_half ) \tXs ai \F la (ov)
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(l,i)+=T(a,i)*F(l,a)",Aux_oo,tX1,interm_t%Fbar_ov)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,i)*L(k,i,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,i)*L(l,i,c,d)",B2,Aux_oo,l2_tensor)
      if (debug_local) print*, "ierr c11", ierr
  
      ! c12 = -one_half (\lambdad klca  \ \Pmm cd | ONE_HALF )\tYd cdkl\tXs ai  \F id (ov) 
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(a,d)+=T(a,i)*F(i,d)",Aux_vv,tX1,interm_t%Fbar_ov)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,d)*L(k,l,c,a)",B2,Aux_vv,l2_tensor,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,c)*L(k,l,d,a)",B2,Aux_vv,l2_tensor)
      if (debug_local) print*, "ierr c12", ierr
  
      ! g10 = + (\lambdad klce (\tYd cdkl || \Pmm cd | ONE_HALF)  \tXs ai \tYd cdkl \W eida (vovv)
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(e,d)+=T(a,i)*W(e,i,d,a)",Aux_vv,tX1,interm_t%w_vovv)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(e,d)*L(k,l,c,e)",B2,Aux_vv,l2_tensor)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(e,c)*L(k,l,d,e)",B2,Aux_vv,l2_tensor,scale=MINUS_ONE)
      if (debug_local) print*, "ierr g10", ierr
  
      ! k10= - \lambdad kmcd (\tYd cdkl || \Pmm kl | ONE_HALF) \tXs ai \tYd cdkl\W lima 
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(l,m)+=T(a,i)*W(l,i,m,a)",Aux_oo,tX1,interm_t%w_ooov)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,m)*L(k,m,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,m)*L(l,m,c,d)",B2,Aux_oo,l2_tensor)
      if (debug_local) print*, "ierr k10", ierr
  
      ! k4 = + \lambdad mlca (\tYd cdkl || \Pmm klcd )   \tXs ai  \W kimd (ooov)
      ierr=talsh_tensor_init(Aux_oovo)
      ierr=talsh_tensor_contract("A(m,l,c,i)+=L(m,l,c,a)*T(a,i)",Aux_oovo,l2_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,c,i)*W(k,i,m,d)",B2,Aux_oovo,interm_t%w_ooov)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,c,i)*W(l,i,m,d)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,d,i)*W(k,i,m,c)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,d,i)*W(l,i,m,c)",B2,Aux_oovo,interm_t%w_ooov)
  
      ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,c,i)*W(k,i,m,d)",B2,Aux_oovo,interm_t%w_ooov,scale=ONE_HALF)
      ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,c,i)*W(l,i,m,d)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE_HALF)
      ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,d,i)*W(k,i,m,c)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE_HALF)
      ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,d,i)*W(l,i,m,c)",B2,Aux_oovo,interm_t%w_ooov,scale=ONE_HALF)
      if (debug_local) print*, "ierr k4", ierr
  
      ! k5 = + \lambdad micd ( \tYd cdkl one_quarter )  \tXs ai \W klma (ooov)
      ierr=talsh_tensor_init(Aux_oooo)
      ierr=talsh_tensor_contract("A(k,l,m,i)+=T(a,i)*W(k,l,m,a)",Aux_oooo,tX1,interm_t%w_ooov)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,l,m,i)*L(m,i,c,d)",B2,Aux_oooo,l2_tensor)
      ! ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,l,m,i)*L(m,i,c,d)",B2,Aux_oooo,l2_tensor,scale=ONE_HALF)  !XY
      if (debug_local) print*, "ierr k5", ierr
  
      ! g4 = -  \lambdad kied( \tYd cdkl || \Pmm klcd ) \tXs ai \W elca (vovv)
      ierr=talsh_tensor_init(Aux_vovo)
      ! ierr=talsh_tensor_contract("A(e,l,c,i)+=T(a,i)*W(k,l,m,a)",Aux_vovo,tX1,interm_t%w_ooov)
      ierr=talsh_tensor_contract("A(e,l,c,i)+=T(a,i)*W(e,l,c,a)",Aux_vovo,tX1,interm_t%w_vovv)  !XY
  
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,l,c,i)*L(k,i,e,d)",B2,Aux_vovo,l2_tensor,scale=MINUS_ONE)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,k,c,i)*L(l,i,e,d)",B2,Aux_vovo,l2_tensor)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,l,d,i)*L(k,i,e,c)",B2,Aux_vovo,l2_tensor)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,k,d,i)*L(l,i,e,c)",B2,Aux_vovo,l2_tensor,scale=MINUS_ONE)
      if (debug_local) print*, "ierr g4", ierr
  
      ! g5 = - \lambdad klea (\tYd cdkl one_quarter )  \tXs ai \W eicd (vovv) 
      ierr=talsh_tensor_init(Aux_oovo)
      ierr=talsh_tensor_contract("A(k,l,e,i)+=L(k,l,e,a)*T(a,i)",Aux_oovo,l2_tensor,tX1)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,l,e,i)*W(e,i,c,d)",B2,Aux_oovo,interm_t%w_vovv,scale=MINUS_ONE)
      if (debug_local) print*, "ierr g5", ierr
      
  
  ! l=-2; t=2 -> h=-2
      ! l7 = -one_half \lambdad klca (\tYd cdkl || \Pmm cd one_half || one_quarter ) \tXd abij  \W ijdb (oovv)
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(a,d)+=T(a,b,i,j)*W(i,j,d,b)",Aux_vv,tX2,int_t%oovv)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,d)*L(k,l,c,a)",B2,Aux_vv,l2_tensor,scale=MINUS_ONE_HALF)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,c)*L(k,l,d,a)",B2,Aux_vv,l2_tensor,scale=ONE_HALF)
      if (debug_local) print*, "ierr l7", ierr
  
      ! l9 =  -one_half \lambdad kicd (\tYd cdkl || \Pmm kl one_half|| one_quarter )  \tXd abij  \W ljab (oovv)
      ierr=talsh_tensor_init(Aux_oo) 
      ! ierr=talsh_tensor_contract("A(l,j)+=T(a,b,i,j)*W(l,j,a,b)",Aux_oo,TX2,int_t%oovv)
      ierr=talsh_tensor_contract("A(l,i)+=T(a,b,i,j)*V(l,j,a,b)",Aux_oo,TX2,int_t%oovv) !XY
      if (debug_local) print*, "ierr l9_1", ierr
      ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,j)*L(k,i,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE_HALF)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,i)*L(k,i,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE_HALF) !XY
      if (debug_local) print*, "ierr l9_2", ierr
      ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,j)*L(l,i,c,d)",B2,Aux_oo,l2_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,i)*L(l,i,c,d)",B2,Aux_oo,l2_tensor,scale=ONE_HALF) !XY
      if (debug_local) print*, "ierr l9", ierr
  
      ! l12 = +one_sixteenth \lambdad ijcd (\tYd cdkl || one_quarter((ij),(ab))  ) \tXd abij \W klab (oovv)
      ierr=talsh_tensor_init(Aux_oooo)
      ierr=talsh_tensor_contract("A(k,l,i,j)+=T(a,b,i,j)*V(k,l,a,b)",Aux_oooo,tX2,int_t%oovv)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,l,i,j)*L(i,j,c,d)",B2,Aux_oooo,l2_tensor,scale=ONE_QUARTER_C)
      if (debug_local) print*, "ierr l12", ierr
  
      ! l13 = +one_sixteenth \lambdad klab \tYd cdkl \tXd abij \W ijcd (oovv)
      ierr=talsh_tensor_init(Aux_oooo)
      ierr=talsh_tensor_contract("A(k,l,i,j)+=L(k,l,a,b)*T(a,b,i,j)",Aux_oooo,l2_tensor,tX2)
      ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,l,i,j)*V(i,j,c,d)",B2,Aux_oooo,int_t%oovv,scale=ONE_QUARTER_C)
      if (debug_local) print*, "ierr l13", ierr
  
      ! l20 : +\lambdad kjcb (\tYd cdkl | \Pmm kl \Pmm cd) \tXd abij \W lida (oovv)
      ierr=talsh_tensor_init(Aux_ovvo)
      ierr=talsh_tensor_contract("A(l,b,d,j)+=W(l,i,d,a)*T(a,b,i,j)",Aux_ovvo,int_t%oovv,tX2)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,b,d,j)*L(k,j,c,b)",B2,Aux_ovvo,l2_tensor)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,b,d,j)*L(l,j,c,b)",B2,Aux_ovvo,l2_tensor,scale=MINUS_ONE)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,b,c,j)*L(k,j,d,b)",B2,Aux_ovvo,l2_tensor,scale=MINUS_ONE)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,b,c,j)*L(l,j,d,b)",B2,Aux_ovvo,l2_tensor)
      if (debug_local) print*, "ierr l20", ierr
  
      ! l30 = - one_half \lambdad mnef (\tYd cdkl | \Pmm kl  one_half (cd)  ) \tXd abij \W jlcd
      ierr=talsh_tensor_init(Aux_oo)
      ierr=talsh_tensor_contract("A(k,j)+=L(i,k,a,b)*T(a,b,i,j)",Aux_oo,l2_tensor,tX2)
      ! ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,j)*W(j,l,c,d)",B2,Aux_oo,int_t%oovv,scale=ONE_HALF)
      ! ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,j)*W(j,k,c,d)",B2,Aux_oo,int_t%oovv,scale=MINUS_ONE_HALF)
  
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,j)*W(j,l,c,d)",B2,Aux_oo,int_t%oovv,scale=MINUS_ONE_HALF)  !XY
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,j)*W(j,k,c,d)",B2,Aux_oo,int_t%oovv,scale=ONE_HALF) !XY
      if (debug_local) print*, "ierr l30", ierr
  
      ! l31 = - one_half \lambdad ijac (\tYd cdkl | \Pmm cd one_half (ij)  ) \tXd abij \W klbd
      ierr=talsh_tensor_init(Aux_vv)
      ierr=talsh_tensor_contract("A(b,c)+=L(i,j,a,c)*T(a,b,i,j)",Aux_vv,l2_tensor,tX2)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(b,c)*W(k,l,b,d)",B2,Aux_vv,int_t%oovv,scale=MINUS_ONE_HALF)
      ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(b,d)*W(k,l,b,c)",B2,Aux_vv,int_t%oovv,scale=ONE_HALF)
      if (debug_local) print*, "ierr l31", ierr
  
  
  
      ierr=talsh_tensor_destruct(Aux_oo)
      ierr=talsh_tensor_destruct(Aux_vv)
      ierr=talsh_tensor_destruct(Aux_vo)
      ierr=talsh_tensor_destruct(Aux_ov)
      
      ierr=talsh_tensor_destruct(Aux_vvvv)
      ierr=talsh_tensor_destruct(Aux_oooo)
      ierr=talsh_tensor_destruct(Aux_oovv)
      ierr=talsh_tensor_destruct(Aux_vooo)
      ierr=talsh_tensor_destruct(Aux_ooov)
      ierr=talsh_tensor_destruct(Aux_vovo)
      ierr=talsh_tensor_destruct(Aux_ovov)
      ierr=talsh_tensor_destruct(Aux_vovv)
      ierr=talsh_tensor_destruct(Aux_oovv)
  
      ierr=talsh_tensor_destruct(Aux_oovo)
      ierr=talsh_tensor_destruct(Aux_ovvo)
  
  
      end subroutine get_tF
                
end module
