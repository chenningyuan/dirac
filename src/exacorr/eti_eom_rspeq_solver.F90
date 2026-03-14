module eti_eom_rspeq_solver
      use eti_davidson
      use exacorr_tensor_interface
      use tensor_algebra
      use eti_common_routines
      use exacorr_datatypes 
      use exacorr_cfg
      use exacorr_utils
      use eti_trial_eomccsd
      use eti_eom
      use exacorr_scaling_parameters
      
      implicit none

      private 

      public solve_eom_lr_right
      public solve_eom_lr_left 
      public get_left_bvec
      public solve_cc_lr_left, get_left_bvec_ccci, get_left_bvec_cccc 
      public print_lr_setup
      public solve_1st_order_rspeq_left_right_driver

      contains

!----------------------------------------------------------------------------------------------------------------

      subroutine print_lr_setup(rspeq_side,operator_label,operator,operator_index,freq)
         complex(kind=8), intent(in) :: freq
         integer, intent(in) :: operator_index
         character(1), intent(in) :: rspeq_side 
         character(1), intent(in) :: operator_label
         type(eti_rsp_interm_tens_t), intent(inout) :: operator


         print *, ""
         print *, "Solving "//rspeq_side//"HS linear response equation"
!        when forming e.g 2-photon transition moments, do not print information on indivudial operators 
         if (operator_index.gt.0) then
                print *, "   Operator "//operator_label//" (",operator_index, &
             &      ") : ",operator%property_label
         end if
         print *, "   Frequency  : ",freq
         print *, ""
      end subroutine


! driver to solve both right and left linear response equations, storing the results in Tprop_tensor and Tprop_bar_tensor
      subroutine solve_1st_order_rspeq_left_right_driver(Tprop_tensor, Tprop_bar_tensor, exa_input, &
                    int_t,interm_t,Interm_t_tensor,l1_tensor,l2_tensor, &
                    Bvec_tensor, omega_rhs, omega_lhs, CC_CC, prop_label, prop_index)

         type(eti_comm_tens_t),intent(inout)       :: Tprop_tensor  !right_vector, output
         type(eti_comm_tens_t),intent(inout)       :: Tprop_bar_tensor !left_vector, output

         type(exacc_input), intent(in) :: exa_input
         type(eti_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
         type(eti_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
         type(eti_rsp_interm_tens_t),intent(inout) :: Interm_t_tensor
         type(eti_comm_tens_t),intent(inout)       :: Bvec_tensor
         complex(kind=8),intent(in)                  :: omega_rhs
         complex(kind=8),intent(in), optional        :: omega_lhs
         type(tensor_t), intent(inout)           :: l2_tensor
         type(tensor_t), intent(inout)           :: l1_tensor
         logical, intent(in) :: CC_CC
         character*1, intent(in) :: prop_label
         integer, intent(in) :: prop_index
         complex(kind=8) :: omega_lhs_local

!        if omega_lhs is not provided, we assume it to be -\omega_rhs; this will be ok for linear response
!        without broadening
         if (present(omega_lhs)) then 
            omega_lhs_local = omega_lhs
         else   
            omega_lhs_local = -omega_rhs
         end if

         call print_lr_setup("R",prop_label,Interm_t_tensor,prop_index,omega_rhs)
         call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensor, omega_rhs, Tprop_tensor)

! will use Tprop_tensor and other intermediates to calculate bvec, which is appears as RHS in the in the LHS response equations [tbar(omega_lhs)*A(omega_lhs) = -bvec]
         if (CC_CC) then
            call get_left_bvec_cccc(exa_input,int_t,interm_t,Interm_t_tensor, &
                    l1_tensor,l2_tensor, Tprop_tensor, Bvec_tensor)
         else
            call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensor,l1_tensor,l2_tensor, Bvec_tensor)
         end if

         call print_lr_setup("L",prop_label,Interm_t_tensor,prop_index,omega_lhs_local)
         if (CC_CC) then
!            call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensor, omega_lhs_local, Bvec_tensor, Tprop_bar_tensor)
            call solve_cc_lr_left(exa_input, interm_t, int_t,  omega_lhs_local, Bvec_tensor, Tprop_bar_tensor)
         else
!            call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensor, omega_lhs_local, Bvec_tensor, Tprop_bar_tensor)
            call solve_eom_lr_left(exa_input, interm_t, int_t, omega_lhs_local, Bvec_tensor, Tprop_bar_tensor)
         end if

      end subroutine solve_1st_order_rspeq_left_right_driver

      subroutine solve_eom_lr_right(exa_input, interm_t, int_t, interm_t_prop, &
                                    omega, tprop_tensor)
            use eti_sigma_eomccsd
            use mfsolver_dirac_input

!     input variables
            type(exacc_input), intent(in) :: exa_input
            type(eti_intermediates_tens_t),  intent(inout) :: interm_t
            type(eti_intg_tens_t), intent(inout) :: int_t
            type(eti_rsp_interm_tens_t), intent(inout) :: interm_t_prop
            complex(kind=8), intent(in)               :: omega

            type(eti_eom_results_tens_t)                 :: results_t
            type(eti_comm_tens_t) , intent(inout)      :: tprop_tensor

!     workspace
            type(eti_davidson_workspace_tens_t) :: workspace_t
            type(eti_projectors_tens_t) :: projector_t
!     error code
            integer :: ierr
            character(len=2) :: flavor = "EE"
            integer :: nroot

            integer(INTD), dimension(2) :: vo_dims
            integer(INTD), dimension(4) :: vvoo_dims

            nroot = exa_input%N_eom_trial
            call eti_davidson_setup(unit_output=6, &
                        convergence_threshold=mfsolver_options_convergence_threshold, &
                        max_subspace_size=mfsolver_options_max_subspace_size,&
                        verbose=mfsolver_options_verbose,&
                        solve_right=.true.,&
                        solve_left=.false.,&
                        solve_linear_sys=.true.,&
                        symmetric=.false., &
                        use_preconditioner=.true.,&
                        max_iterations=mfsolver_options_max_iterations,&
                        complex_mode=.true.)

            call print_date("Entering CC-CI Right Linear Systems calculation")

            call eti_davidson_workspace_create(workspace_t, nroot, interm_t%nocc, interm_t%nvir, flavor, 2)
            call eti_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
            call eti_eomccsd_LR_proper_create(workspace_t,interm_t_prop%epsilon_prop_vo)
            call eti_eomccsd_LR_trial_create(workspace_t)

!     set pointers to sigma vector subroutines 
            interm_t%sigma_rhs_p => eti_sigma_eomccsd_rhs_ee

            vo_dims(1) = exa_input%nvir
            vo_dims(2) = exa_input%nocc
            vvoo_dims(1:2) = exa_input%nvir
            vvoo_dims(3:4) = exa_input%nocc

            ierr=construct_tensor_t(tprop_tensor%t1,C8,vo_dims,init_val=ZERO)
            ierr=construct_tensor_t(tprop_tensor%t2,C8,vvoo_dims,init_val=ZERO)


            call eti_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, flavor,& 
                                    omega, interm_t_prop%epsilon_prop_vo, interm_t_prop%epsilon_prop_vvoo, tprop_tensor)

            call print_date("Finishing CC-CI Right Linear Systems calculation")

            call eti_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_right


!----------------------------------------------------------------------------------------------------------------
      
      subroutine solve_eom_lr_left(exa_input, interm_t, int_t, omega, bvec_tensor, tprop_bar_tensor)
            use eti_sigma_eomccsd
            use mfsolver_dirac_input

      !     input variables
            type(exacc_input), intent(in) :: exa_input
            type(eti_intermediates_tens_t),  intent(inout) :: interm_t
            type(eti_intg_tens_t), intent(inout) :: int_t
!            type(eti_rsp_interm_tens_t), intent(inout) :: interm_t_prop
            complex(kind=8),INTENT(IN)              :: omega
      
            type(eti_eom_results_tens_t)  :: results_t
            type(eti_comm_tens_t)       :: tprop_bar_tensor
            type(eti_comm_tens_t)       :: bvec_tensor
      !     workspace
            type(eti_davidson_workspace_tens_t) :: workspace_t
            type(eti_projectors_tens_t) :: projector_t
      !     error code
            integer :: ierr
            integer :: nroots = 1
            character(len=2) :: flavor = "EE"

            integer(INTD), dimension(2) :: ov_dims
            integer(INTD), dimension(4) :: oovv_dims
                 
            call eti_davidson_setup(unit_output=6, &
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
      
            call eti_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2, 'L')
            call eti_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, 'L')
            call eti_eomccsd_LR_proper_create(workspace_t,bvec_tensor%t1)
            call eti_eomccsd_LR_trial_create(workspace_t)
      
      !     set pointers to sigma vector subroutines 
            interm_t%sigma_lhs_p => eti_sigma_eomccsd_lhs_ee

            ov_dims(1) = exa_input%nocc
            ov_dims(2) = exa_input%nvir
            oovv_dims(1:2) = exa_input%nocc
            oovv_dims(3:4) = exa_input%nvir

    
            ierr=construct_tensor_t(tprop_bar_tensor%t1,C8,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(tprop_bar_tensor%t2,C8,oovv_dims,init_val=ZERO)

      
            call eti_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, flavor,&
                                          -omega, bvec_tensor%t1, bvec_tensor%t2, tprop_bar_tensor)   !left equation use -omega
                                          
            call print_date("Finishing CC-CI Left Linear Systems calculation")
      
            call eti_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_left

!---------------------------------------------------------------------------------------------------------------------
      
      subroutine solve_cc_lr_left(exa_input, interm_t, int_t,  omega, bvec_tensor, tprop_bar_tensor)
            use eti_sigma_eomccsd
            use mfsolver_dirac_input

      !     input variables
            type(exacc_input), intent(in) :: exa_input
            type(eti_intermediates_tens_t),  intent(inout) :: interm_t
            type(eti_intg_tens_t), intent(inout) :: int_t
!            type(eti_rsp_interm_tens_t), intent(inout) :: interm_t_prop
            complex(kind=8),INTENT(IN)              :: omega
      

            type(eti_eom_results_tens_t)               :: results_t
            type(eti_comm_tens_t), intent(inout)       :: tprop_bar_tensor
            type(eti_comm_tens_t), intent(inout)       :: bvec_tensor
      !     workspace
            type(eti_davidson_workspace_tens_t) :: workspace_t
            type(eti_projectors_tens_t) :: projector_t
      !     error code
            integer :: ierr
            integer :: nroots = 1
            character(len=2) :: flavor = "EE"

            integer(INTD), dimension(2) :: ov_dims
            integer(INTD), dimension(4) :: oovv_dims
            
            call eti_davidson_setup(unit_output=6, &
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
      
            call eti_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2, 'L')
            call eti_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, 'L')
            ! call eti_eomccsd_ee_trial_create(workspace_t)
            call eti_eomccsd_LR_proper_create(workspace_t,bvec_tensor%t1)
            call eti_eomccsd_LR_trial_create(workspace_t)
      
      !     set pointers to sigma vector subroutines 
            interm_t%sigma_lhs_p => eti_sigma_eomccsd_lhs_ee

            ov_dims(1) = exa_input%nocc
            ov_dims(2) = exa_input%nvir
            oovv_dims(1:2) = exa_input%nocc
            oovv_dims(3:4) = exa_input%nvir

      
            ierr=construct_tensor_t(tprop_bar_tensor%t1,C8,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(tprop_bar_tensor%t2,C8,oovv_dims,init_val=ZERO)

            call eti_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, flavor,&
                                    -omega, bvec_tensor%t1, bvec_tensor%t2, tprop_bar_tensor)   !left equation use -omega
                                    
            call print_date("Finishing CC-CC Left Linear Systems calculation")
      
            call eti_davidson_workspace_destroy(workspace_t)

      end subroutine solve_cc_lr_left

!---------------------------------------------------------------------------------------------------------------------
      
      subroutine get_left_bvec(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,bvec_tensor)

      !XY: Although the final QR result of EOM is same as CC-CI, the response equations and 
      ! the left perturbed ampiltudes tbar are not same. The point is the left bevctor EOM used is not the one of CC-CI 
      ! CC-CI left bvector is depended on the right perturbed ampiltudes. So, here we choose the EOM type and comment out the CC-CI as backup
      ! see more details: J. Chem Phys, 142, p.114109, 2015, & J. Chem Phys, 144, p.024102, 2016

            use tensor_algebra
            use exacorr_datatypes
            use exacorr_utils
            use eti_common_routines
            use exacorr_global


            type(exacc_input), intent(in) :: exa_input
            type(eti_intg_tens_t), intent(inout) :: int_t
            type(eti_intermediates_tens_t),intent(inout) :: interm_t
            type(eti_rsp_interm_tens_t),intent(inout) :: interm_t_prop
            type(tensor_t), intent(inout)  :: l2_tensor
            type(tensor_t), intent(inout)  :: l1_tensor
            type(eti_comm_tens_t), intent(inout)       :: bvec_tensor


            type(tensor_t)       :: Aux_tensor1, one_tensor

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


            ierr=construct_tensor_t(bvec_tensor%t1,C8,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(bvec_tensor%t2,C8,oovv_dims,init_val=ZERO)
            ierr=construct_tensor_t(Aux_tensor1,C8,one_dims(1:0),init_val=ZERO)


            ierr=t_einsum('B(i,a)+=Y(i,a)',bvec_tensor%t1,interm_t_prop%eta_prop_ov,scale=ONE) 
            ierr=t_einsum("B(i,a)+=l2(i,j,a,b)*E(b,j)",bvec_tensor%t1,l2_tensor,& 
                                                interm_t_prop%epsilon_prop_vo,scale=ONE)  

            ierr=t_einsum('B(i,j,a,b)+=Y(i,j,a,b)',bvec_tensor%t2,interm_t_prop%eta_prop_oovv,scale=ONE)  

            
!           EOM-CC
            ierr=t_einsum("A1()+=l1(i,a)*E(a,i)",Aux_tensor1,l1_tensor,interm_t_prop%epsilon_prop_vo,scale=ONE)
            ierr=t_einsum("A1()+=l2(i,j,a,b)*E(a,b,i,j)",Aux_tensor1,l2_tensor,&
                                          interm_t_prop%epsilon_prop_vvoo,scale=ONE_QUARTER)
 
            ierr=t_einsum("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
            ierr=t_einsum("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
                                                Aux_tensor1,scale=MINUS_ONE) 

            ierr=destroy_tensor_t(Aux_tensor1)


      end subroutine get_left_bvec



!-----------------------------------------------------------------------------------------------------------------------------
      
      subroutine get_left_bvec_ccci(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,Tprop_tensor,bvec_tensor)

            !XY: the ccci left bvector in the left response equations. 
            !    This routine is just for backup since we implemented EOM for quadratic response
      
                  use tensor_algebra
                  use exacorr_datatypes
                  use exacorr_utils
                  use eti_common_routines
                  use exacorr_global
      
      
                  type(exacc_input), intent(in) :: exa_input
                  type(eti_intg_tens_t), intent(inout) :: int_t
                  type(eti_intermediates_tens_t),intent(inout) :: interm_t
                  type(eti_rsp_interm_tens_t) :: interm_t_prop
                  type(tensor_t), intent(inout)  :: l2_tensor
                  type(tensor_t), intent(inout)  :: l1_tensor
                  type(eti_comm_tens_t), intent(inout)       :: Tprop_tensor    !for CC-CI, we need right-vectors
                  type(eti_comm_tens_t), intent(inout)       :: bvec_tensor
      
      
                  type(tensor_t)       :: Aux_tensor1, one_tensor
      
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
      
      
                  ierr=construct_tensor_t(bvec_tensor%t1,C8,ov_dims,init_val=ZERO)
                  ierr=construct_tensor_t(bvec_tensor%t2,C8,oovv_dims,init_val=ZERO)
                  ierr=construct_tensor_t(Aux_tensor1,C8,one_dims(1:0),init_val=ZERO)
      
      
      
                  ierr=t_einsum('B(i,a)+=Y(i,a)',bvec_tensor%t1,interm_t_prop%eta_prop_ov,scale=ONE) 
                  ierr=t_einsum("B(i,a)+=l2(i,j,a,b)*E(b,j)",bvec_tensor%t1,l2_tensor,& 
                                                      interm_t_prop%epsilon_prop_vo,scale=ONE)  
      
                  ierr=t_einsum('B(i,j,a,b)+=Y(i,j,a,b)',bvec_tensor%t2,interm_t_prop%eta_prop_oovv,scale=ONE)  
      
                  
                  ierr=t_einsum("A1()+=l1(i,a)*E(a,i)",Aux_tensor1,l1_tensor,interm_t_prop%epsilon_prop_vo,scale=ONE)
                  ierr=t_einsum("A1()+=l2(i,j,a,b)*E(a,b,i,j)",Aux_tensor1,l2_tensor,&
                                                interm_t_prop%epsilon_prop_vvoo,scale=ONE_QUARTER)
       
                  ierr=t_einsum("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
                  ierr=t_einsum("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
                                                      Aux_tensor1,scale=MINUS_ONE) 
      
                  ierr=destroy_tensor_t(Aux_tensor1)
      

                  call get_eta(int_t,interm_t,l1_tensor,l2_tensor,interm_t%eta_ov,interm_t%eta_oovv)
                  ierr=t_einsum("A1()+=l1(i,a)*T(a,i)",Aux_tensor1,l1_tensor,Tprop_tensor%t1,scale=ONE)
                  ierr=t_einsum("A1()+=l2(i,j,a,b)*T(a,b,i,j)",Aux_tensor1,l2_tensor,&
                                                Tprop_tensor%t2,scale=ONE_QUARTER)
            
      
                  ierr=t_einsum("B(i,a)+=n(i,a)*A1()",bvec_tensor%t1,interm_t%eta_ov,Aux_tensor1,scale=ONE) 
                  ierr=t_einsum("B(i,j,a,b)+=n(i,j,a,b)*A1()",bvec_tensor%t2,interm_t%eta_oovv, &
                                                Aux_tensor1,scale=ONE)  
      
      
                  ierr=initialize_tensor_t(Aux_tensor1)
                  ! call get_eta(int_t,interm_t,l1_tensor,l2_tensor,interm_t%eta_ov,interm_t%eta_oovv)
      
                  ierr=t_einsum("A1()+=n1(i,a)*T(a,i)",Aux_tensor1,interm_t%eta_ov,Tprop_tensor%t1,scale=MINUS_ONE)
                  ierr=t_einsum("A1()+=n2(i,j,a,b)*T(a,b,i,j)",Aux_tensor1,interm_t%eta_oovv,&
                                                Tprop_tensor%t2,scale=MINUS_ONE_QUARTER)
      
      
                  ierr=direct_data_pointer(Aux_tensor1,body_p)
      
                  ierr=t_einsum("B(i,a)+=l1(i,a)*A1()",bvec_tensor%t1,l1_tensor,Aux_tensor1,scale=MINUS_ONE) 
                  ierr=t_einsum("B(i,j,a,b)+=l2(i,j,a,b)*A1()",bvec_tensor%t2,l2_tensor, &
                                                Aux_tensor1,scale=MINUS_ONE)  
      
                  ierr=destroy_tensor_t(Aux_tensor1)
                  ierr=destroy_tensor_t(one_tensor)
      
            end subroutine get_left_bvec_ccci

!-----------------------------------------------------------------------------------------------------------------------------
      
            subroutine get_left_bvec_cccc(exa_input,int_t,interm_t,interm_t_prop,l1_tensor,l2_tensor,Tprop_tensor,bvec_tensor)

                  !The cc-cc left bvector in the left response equations. 
                  !This routine can be used in cc-cc quadratic response function 
                  !and asymmetric cc-cc linear response function
                  !Written by Xiang Yuan Janunary 2024, Lille, France

            
                  use tensor_algebra
                  use exacorr_datatypes
                  use exacorr_utils
                  use eti_common_routines
                  use exacorr_global
                  use eti_eom_LR, only : get_tF
            
            
                  type(exacc_input), intent(in) :: exa_input
                  type(eti_intg_tens_t), intent(inout) :: int_t
                  type(eti_intermediates_tens_t),intent(inout) :: interm_t
                  type(eti_rsp_interm_tens_t),intent(inout) :: interm_t_prop
                  type(tensor_t), intent(inout)  :: l2_tensor
                  type(tensor_t), intent(inout)  :: l1_tensor
                  type(eti_comm_tens_t), intent(inout) :: Tprop_tensor    !for CC-CI, we need right-vectors
                  type(eti_comm_tens_t), intent(inout) :: bvec_tensor
            
            
                  type(tensor_t)       :: Aux_tensor1, one_tensor
                  type(eti_comm_tens_t)  :: B  !Ft^x
            
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
            
            
                  ierr=construct_tensor_t(bvec_tensor%t1,C8,ov_dims,init_val=ZERO)
                  ierr=construct_tensor_t(bvec_tensor%t2,C8,oovv_dims,init_val=ZERO)

                  ierr=construct_tensor_t(B%t1,C8,ov_dims,init_val=ZERO)
                  ierr=construct_tensor_t(B%t2,C8,oovv_dims,init_val=ZERO)

                  print *, 'DEBUG, get_left_bvec_cccc : norm bvec #1',&
                            tensor_norm2(bvec_tensor%t1),tensor_norm2(bvec_tensor%t2)
                  call get_tF(B%t1,B%t2,Tprop_tensor%t1,Tprop_tensor%t2,l1_tensor,l2_tensor,&
                        interm_t,int_t,exa_input%nocc,exa_input%nvir)
            
                  ierr=t_einsum('B(i,a)+=Y(i,a)',bvec_tensor%t1,interm_t_prop%eta_prop_ov,scale=ONE)             
                  ierr=t_einsum('B(i,j,a,b)+=Y(i,j,a,b)',bvec_tensor%t2,interm_t_prop%eta_prop_oovv,scale=ONE)  

                  ierr=t_einsum('B(i,a)+=Ft(i,a)',bvec_tensor%t1,B%t1,scale=ONE)             
                  ierr=t_einsum('B(i,j,a,b)+=Ft(i,j,a,b)',bvec_tensor%t2,B%t2,scale=ONE) 
                  print *, 'DEBUG, get_left_bvec_cccc : norm bvec #2',&
                            tensor_norm2(bvec_tensor%t1),tensor_norm2(bvec_tensor%t2)

                  ierr=destroy_tensor_t(B%t1)
                  ierr=destroy_tensor_t(B%t2)

            
      end subroutine get_left_bvec_cccc

!------------------------------------------------------------------------------------------------------------------------------------

      subroutine get_eta(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)
            ! Used in CC-CI left equation to construct eta vector. In EOM-QR, we don't need this term.

            !     input variables
                  type(eti_intg_tens_t), intent(inout) :: int_t
                  type(eti_intermediates_tens_t),  intent(inout) :: interm_t
                  type(tensor_t),    intent(inout) :: sigma1_tensor
                  type(tensor_t),    intent(inout) :: sigma2_tensor
                  type(tensor_t),    intent(inout) :: l1_tensor
                  type(tensor_t),    intent(inout) :: l2_tensor
                  type(tensor_t) :: one_tensor
            !     auxiliary tensor
                  type(tensor_t) :: goo_tensor, gvv_tensor, oo_aux
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
                  ierr=construct_tensor_t(goo_tensor,C8,oo_dims,init_val=ZERO)
                  ierr=construct_tensor_t(gvv_tensor,C8,vv_dims,init_val=ZERO)
                  ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
                  ierr=t_einsum("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
                  ierr=t_einsum("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
            
            !***************
            !     l1_sigma *
            !***************
            ierr=construct_tensor_t(sigma2_tensor,C8,oovv_dims,init_val=ZERO)
            ierr=construct_tensor_t(sigma1_tensor,C8,ov_dims,init_val=ZERO)

            !  ierr=initialize_tensor_t(sigma1_tensor)
            
                  
            !-----------------------------------------
            !     term 2: S1(i,a) += L(i,e) * F(e,a) |
            !-----------------------------------------
            
            ierr=t_einsum("S(i,a)+=L(i,e)*F(e,a)",sigma1_tensor,l1_tensor,interm_t%fbar_vv)
            
            !-----------------------------------------
            !     term 3: S1(i,a) -= L(m,a) * F(i,m) |
            !-----------------------------------------
            
            ierr=t_einsum("S(i,a)+=L(m,a)*F(i,m)",sigma1_tensor,l1_tensor,interm_t%fbar_oo,scale=MINUS_ONE)
            
            !-----------------------------------------------------
            !     term 4: S1(i,a) += 1/2 L(i,m,e,f) * W(e,f,a,m) |
            !-----------------------------------------------------
            
            ierr=t_einsum("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",sigma1_tensor,l2_tensor,interm_t%w_vvvo,scale=ONE_HALF)
            
            !---------------------------------------------
            !     term 5: S1(i,a) -= G(f,e) * W(e,i,f,a) |
            !---------------------------------------------
            
            ierr=t_einsum("S(i,a)+=G(f,e)*W(e,i,a,f)",sigma1_tensor,gvv_tensor,interm_t%w_vovv)
            
            !---------------------------------------------
            !     term 6: S1(i,a) -= G(n,m) * W(m,i,n,a) |
            !---------------------------------------------
            
            ierr=t_einsum("S(i,a)+=G(n,m)*W(i,m,n,a)",sigma1_tensor,goo_tensor,interm_t%w_ooov)

                        
            !---------------------------------------------
            !     term 7: S1(i,a) += L(m,e) * W(i,e,a,m) |
            !---------------------------------------------
            
            ierr=t_einsum("S(i,a)+=L(m,e)*W(e,i,a,m)",sigma1_tensor,l1_tensor,interm_t%w_vovo,scale=MINUS_ONE)
            
            
            !-----------------------------------------------------
            !     term 8: S1(i,a) -= 1/2 L(m,n,a,e) * W(i,e,m,n) |
            !-----------------------------------------------------
            
            ierr=t_einsum("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",sigma1_tensor,l2_tensor,interm_t%w_ovoo,scale=ONE_HALF)
            
            !***************
            !     l2_sigma *
            !***************
            
            ! ierr=initialize_tensor_t(sigma2_tensor)
            
            !--------------------------------------------------------
            !     term 2: S2(i,j,a,b) += P(a,b) L(i,j,a,e) * F(e,b) |
            !--------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
            ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
            
            
            !--------------------------------------------------------
            !     term 3: S2(i,j,a,b) -= P(i,j) L(i,m,a,b) * F(j,m) |
            !--------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
            ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
            
            !---------------------------------------------------------
            !     term 4: S2(i,j,a,b) += 1/2 L(m,n,a,b) * W(i,j,m,n) |
            !---------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",sigma2_tensor,l2_tensor,interm_t%w_oooo,scale=ONE_HALF)
            
            
            !-------------------------------------------------------------------
            !     term 5: S2(i,j,a,b) += P(i,j) P(a,b) L(i,m,a,e) * W(j,e,b,m) |
            !-------------------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            ierr=t_einsum("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            ierr=t_einsum("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
            
            !--------------------------------------------------------
            !     term 6: S2(i,j,a,b) -= P(a,b) V(i,j,a,e) * G(e,b) |
            !--------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",sigma2_tensor,int_t%oovv,gvv_tensor)   !XY
            ierr=t_einsum("S(i,j,a,b)+=V(i,j,b,e)*G(e,a)",sigma2_tensor,int_t%oovv,gvv_tensor,scale=MINUS_ONE)   !XY
            
            
            !-------------------------------------------------
            !     term 7: S2(i,j,a,b) -= P(a,b) L(m,a) * W(i,j,m,b) |
            !-------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(m,a)*W(i,j,m,b)",sigma2_tensor,l1_tensor,interm_t%w_ooov,scale=MINUS_ONE)  !XY
            ierr=t_einsum("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",sigma2_tensor,l1_tensor,interm_t%w_ooov)   !XY
            
            !--------------------------------------------------------
            !     term 8: S2(i,j,a,b) -= P(i,j) V(i,m,a,b) * G(j,m) |
            !--------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=V(m,i,a,b)*G(j,m)",sigma2_tensor,int_t%oovv,goo_tensor)   !XY
            ierr=t_einsum("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",sigma2_tensor,int_t%oovv,goo_tensor)   !XY
            
            !---------------------------------------------------------------------
            !     term 9: S2(i,j,a,b) -= P(i,j) V(m,j,a,b) * ( L(i,e) * T(e,m) ) |
            !---------------------------------------------------------------------
            
                  ! if (.not.interm_t%ccd) then
                  ! !   oo_dims = (/interm_t%nocc,interm_t%nocc/)
                  !   ierr=construct_tensor_t(oo_aux,C8,oo_dims,init_val=ZERO)
                  ! !   ierr=t_einsum("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1)
                  !   ierr=t_einsum("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1,scale=MINUS_ONE)  !XY
                  !   ierr=t_einsum("S(i,j,a,b)+=V(m,j,a,b)*U(i,m)",sigma2_tensor,int_t%oovv,oo_aux)
                  !   ierr=t_einsum("S(i,j,a,b)+=V(i,m,a,b)*U(j,m)",sigma2_tensor,int_t%oovv,oo_aux)
                  !   ierr=destroy_tensor_t(oo_aux)
                  ! end if
            
            !---------------------------------------------------------
            !     term 10: S2(i,j,a,b) += P(i,j) L(i,e) * V(e,j,a,b) |
            !---------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",sigma2_tensor,l1_tensor,interm_t%w_vovv)   !XY
            ierr=t_einsum("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",sigma2_tensor,l1_tensor,interm_t%w_vovv)   !XY
            
            
            !------------------------------------------------------------
            !     term 11: S2(i,j,a,b) += P(i,j) P(a,b) L(i,a) * F(j,b) |
            !------------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(i,a)*F(j,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
            ierr=t_einsum("S(i,j,a,b)+=L(j,a)*F(i,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
            ierr=t_einsum("S(i,j,a,b)+=L(i,b)*F(j,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
            ierr=t_einsum("S(i,j,a,b)+=L(j,b)*F(i,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
            
            !--------------------------------------------------------
            !     term 12: S2(i,j,a,b) += 1/2 L(i,j,e,f) W(e,f,a,b) |
            !--------------------------------------------------------
            
            ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",sigma2_tensor,l2_tensor,interm_t%w_vvvv,scale=ONE_HALF)

            
            ierr=destroy_tensor_t(goo_tensor)
            ierr=destroy_tensor_t(gvv_tensor)
            
            
      end subroutine get_eta
                
end module
