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
                         MINUS_TWO=(-2.D0,0.D0)

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

      subroutine solve_eom_lr_right(exa_input, interm_t, int_t, interm_t_prop, omega, tprop_tensor)
            use talsh_sigma_eomccsd

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
                        solve_right=.true.,&
                        solve_left=.false.,&
                        solve_linear_sys=.true.,&
                        symmetric=.false., &
                        use_preconditioner=.false.,&
                        ! max_subspace_size=nroot,&
                        complex_mode=.true.)

            call print_date("Entering CC-CI Right Linear Systems calculation")

            call talsh_davidson_workspace_create(workspace_t, nroot, interm_t%nocc, interm_t%nvir, flavor, 2)
            call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
            call talsh_eomccsd_ee_trial_create(workspace_t)

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
      
      subroutine solve_eom_lr_left(interm_t, int_t, interm_t_prop, omega, bvec_tensor, tprop_bar_tensor)
            use talsh_sigma_eomccsd

            ! write(*,*) ""
            ! write(*,'(A)') "*******************************************************************************"
            ! write(*,'(A)') "***************************** CC-CI LR Left ******************************"
            ! write(*,'(A)') "*******************************************************************************"
            ! write(*,*) ""
      
      
      !     input variables
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
                 
                  call talsh_davidson_setup(unit_output=6, &
                              solve_right=.false.,&
                              solve_left=.true.,&
                              solve_linear_sys=.true.,&
                              symmetric=.false., &
                              use_preconditioner=.false.,&
                              complex_mode=.true.)
      
                  call print_date("Entering EOM-EE eigenvalue calculation")
      
                  call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2)
                  call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
      !     call talsh_eomccsd_ee_trial_create(workspace_t)
      
      !     set pointers to sigma vector subroutines 
                  interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ee
      
                  call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t, &
                                          omega, bvec_tensor%t1, bvec_tensor%t2)
                                          
                  call print_date("Finishing EOM-EE eigenvalue calculation")
      
                  call talsh_davidson_workspace_destroy(workspace_t)

      end subroutine solve_eom_lr_left

!---------------------------------------------------------------------------------------------------------------------
      
      subroutine get_left_bvec(exa_input,interm_t,interm_t_prop,l1_tensor,l2_tensor,tprop_tensor,bvec_tensor)
            use tensor_algebra
            use talsh
            use exacorr_datatypes
            use exacorr_utils
            use talsh_common_routines
            use exacorr_global


            type(exacc_input), intent(in) :: exa_input
            type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
            type(talsh_rsp_interm_tens_t) :: interm_t_prop
            type(talsh_tens_t), intent(inout)  :: l2_tensor
            type(talsh_tens_t), intent(inout)  :: l1_tensor
            type(talsh_comm_tens_t)       :: tprop_tensor
            type(talsh_comm_tens_t)       :: bvec_tensor


            type(talsh_tens_t)       :: Aux_tensor1, Aux_tensor2

            integer(INTD), dimension(4)   :: oovv_dims
            integer(INTD), dimension(2)   :: ov_dims
            integer(C_INT)              :: one_dims(1)

            integer :: ierr


            one_dims(1) = 1
            oovv_dims(1:2) = exa_input%nocc
            oovv_dims(3:4) = exa_input%nvir
            ov_dims(1) = exa_input%nocc
            ov_dims(2) = exa_input%nvir


            ierr=talsh_tensor_construct(bvec_tensor%t1,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(bvec_tensor%t2,C8,oovv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(Aux_tensor1,C8,one_dims(1:0),init_val=ZERO)
            ierr=talsh_tensor_construct(Aux_tensor2,C8,one_dims(1:0),init_val=ZERO)


            ierr=talsh_tensor_contract("A1()+=l1(i,a)*TX(a,i)",Aux_tensor1,l1_tensor,tprop_tensor%t1,scale=ONE)
            ierr=talsh_tensor_contract("A1()+=l2(i,j,a,b)*TX(a,b,i,j)",Aux_tensor1,l2_tensor,tprop_tensor%t2,scale=ONE)

            ierr=talsh_tensor_contract("A2()+=N(i,a)*TX(a,i)",Aux_tensor2,interm_t%eta_ov,tprop_tensor%t1,scale=ONE)
            ierr=talsh_tensor_contract("A2()+=N(i,j,a,b)*TX(a,b,i,j)",Aux_tensor2,interm_t%eta_oovv,tprop_tensor%t2,scale=ONE)


            ierr=talsh_tensor_add('B(i,a)+=Y(i,a)',bvec_tensor%t1,interm_t_prop%eta_prop_ov,scale=MINUS_ONE)   
            ierr=talsh_tensor_contract("B(i,a)+=l1(i,a)*A2()",bvec_tensor%t1,l1_tensor,Aux_tensor2,scale=ONE)   
            ierr=talsh_tensor_contract("B(i,a)+=A1()*N(i,a)",bvec_tensor%t1,interm_t%eta_ov,Aux_tensor2,scale=ONE)   


            ierr=talsh_tensor_add('B(i,j,a,b)+=Y(i,j,a,b)',bvec_tensor%t2,interm_t_prop%eta_prop_oovv,scale=MINUS_ONE)   
            ierr=talsh_tensor_contract("B(i,j,a,b)+=l2(i,j,a,b)*A2()",bvec_tensor%t2,l2_tensor,Aux_tensor2,scale=ONE)   
            ierr=talsh_tensor_contract("B(i,j,a,b)+=A1()*N(i,j,a,b)",bvec_tensor%t2,interm_t%eta_oovv,Aux_tensor2,scale=ONE) 

            ierr=talsh_tensor_destruct(Aux_tensor1)
            ierr=talsh_tensor_destruct(Aux_tensor2)


      end subroutine get_left_bvec


end module
