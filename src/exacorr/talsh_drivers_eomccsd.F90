module talsh_eom

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates
  use talsh_ao_to_mo, only : tensor_norm2
  use exacorr_global, only : debug_exa
  use exacorr_utils,  only : print_date
  use talsh_davidson
  use talsh_trial_eomccsd

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0),ONE=(1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)
  real(8),    parameter :: ONE_QUARTER=0.25D0

  type, public :: talsh_eom_results_tens_t
     character, pointer :: flavor(:)
     type(talsh_davidson_results_tens_t) :: value_vectors_t
  end type talsh_eom_results_tens_t

  private

  public talsh_eom_eigval_driver

  contains

    subroutine talsh_steom_eigval_driver (results_t,exa_input,eps_occ,eps_vir,interm_t,int_t)

!     input variables
      type(exacc_input),     intent(in   )           :: exa_input
      type(talsh_intg_tens_t), intent(inout)           :: int_t
      real(8),               intent(in)              :: eps_occ(:), eps_vir(:)
      type(talsh_intermediates_tens_t), intent(inout)  :: interm_t
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t

!     local variables
      type(talsh_davidson_workspace_tens_t) :: workspace_t
      type(talsh_projectors_tens_t) :: projector_t

      write(*,*) ""
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A)') "************************** ST Equation of Motion ******************************"
      write(*,'(A)') "*******************************************************************************"
      write(*,*) ""

      call talsh_davidson_setup(unit_output=6, &
                                solve_right=.true.,&
                                solve_left=.false.,&
                                symmetric=.false., &
                                complex_mode=.true.)

      call print_date("Entering the steom_eigval_driver routine")
!     todo : store the results from each of ip/ea runs 
!     we also need to control whether we need to do IP and EA, depending on the options. for now sticking to STEOM-EE, which needs both
      call print_date("Setting up RHS EOM-IP run")
      if (.not.exa_input%ccd) then
!        call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, "IP", 2)
!        call talsh_davidson_create_results(eomip_results_t)
!        call eom_ip(eomip_results_t,left,right,exa_input,workspace_t,projector_t,interm_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
!        call eom_ip(exa_input,workspace_t,projector_t,interm_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
!        call talsh_davidson_workspace_destroy(workspace_t)
         call print_date("Done with  RHS EOM-IP run")

         call print_date("Setting up RHS EOM-EA run")
!        call talsh_davidson_create_results(eomea_results_t)
!        call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, "EA", 2)
!        call eom_ea(exa_input,workspace_t,projector_t,interm_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
!        call talsh_davidson_workspace_destroy(workspace_t)
      else
!        call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, "IP", 2)
!        call eom_ip(exa_input,workspace_t,projector_t,interm_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
!        call talsh_davidson_workspace_destroy(workspace_t)
!        call print_date("Done with  RHS EOM-IP run")

!        call print_date("Setting up RHS EOM-EA run")
!        call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, "EA", 2)
!        call eom_ea(exa_input,workspace_t,projector_t,interm_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
!        call talsh_davidson_workspace_destroy(workspace_t)
      end if
      call print_date("Done with  RHS EOM-EA run")

! build the S coefficients from the EOM-IP/EA vectors that go into (Hbar*S)C in the G matrix
! use stuff in eomip_results_t and eomea_results_t

! build the "g" coefficients 

! build the G matrix over the active space

! diagonalize G and get the STEOM eigenvalues. G should likely be a Fortran array

! store eigenvalues and eigenvectors, to do transition moments etc

      print*, ""
      call print_date("Leaving the steom_eigval_driver routine")
      return


    end subroutine


! i think we should split drivers into:
!
! *_eigval_driver : solves for the eigenvalues and eigenvectors, left and/or right, stores results in derived type keeping eigenvalue/vector results. 
!                      this will use e.g davidson to solve (A - lambda*I)*X = 0 by expanding X into trial vectors B. for this solutions, we can optionally provide projection operators (to do e.g. cvs)
! 
! *_dmat_driver   : forms density matrices (1-body, 2-body) from results from *_eigval_driver, and additionally getting info from ground-state lambda, to prepare density or transition density matrices 
!
! *_expval_driver : uses results from *_eigenvalue_driver and from *_dmat_driver to build density matrices and contract them with property operators
!
! *_lrsp_driver   : solves for modified eigenvalue equations, using *_eigval_driver, for the different perturbations specified, forms the linear response function
!                      this will use e.g davidson to solve the generalized linear system (A - omega*S)*X - C = 0, with omega a frequency (number), S a metric (matrix) and C a vector (=property gradient, etc), 
!                      we should also be able to apply projectors to the (Hbar - omegaS)X part.
!                      for the generalized linear problem, we should store the vectors for each omega, keeping track of which C the solutions are associated to
! *_qrsp_driver   : uses *_linresp_driver to solve the linear response equations involved in the quadratic response, and forms the quadratic response function
!
! and within each use the different variants : EOM-IP/EE/EA/../-CC* where CC* refers to {P-,}EOM-{CC{SD,SDT},MBPT2}, CCn,
    subroutine talsh_eom_eigval_driver (results_t,exa_input,eps_occ,eps_vir,interm_t,int_t)

!     input variables
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
      type(exacc_input),     intent(in   )           :: exa_input
      type(talsh_intg_tens_t), intent(inout)         :: int_t
      type(talsh_intermediates_tens_t), intent(inout):: interm_t
      real(8),               intent(in)              :: eps_occ(:), eps_vir(:)
!     eom tensors


      write(*,*) ""
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A)') "***************************** Equation of Motion ******************************"
      write(*,'(A)') "*******************************************************************************"
      write(*,*) ""

      call talsh_davidson_setup(unit_output=6, &
                                solve_right=.true.,&
                                solve_left=.false.,&
                                symmetric=.false., &
                                complex_mode=.true.)

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      select case(exa_input%eom_flavour)
         case("IP")
            call eom_ccsd_eigenvalue_ip(results_t,exa_input,interm_t,int_t)
         case("EA")
            call eom_ccsd_eigenvalue_ea(results_t,exa_input,interm_t,int_t)
!        case("EE")
!           call eom_ccsd_eigenvalue_ee(results_t,exa_input,interm_t,int_t)
         case default
            call quit ('Error: No valid EOM-flavour chosen')
      end select

      if (exa_input%eom_properties) then  
         call talsh_davidson_setup(unit_output=6, &
                                   solve_right=.false.,&
                                   solve_left=.true.,&
                                   symmetric=.false., &
                                   complex_mode=.true.)
         select case(exa_input%eom_flavour)
            case("IP")
               call eom_ccsd_eigenvalue_ip(results_t,exa_input,interm_t,int_t)
            case("EA")
               call eom_ccsd_eigenvalue_ea(results_t,exa_input,interm_t,int_t)
!           case("EE")
!              call eom_ccsd_eigenvalue_ee(results_t,exa_input,interm_t,int_t)
            case default
               call quit ('Error: No valid EOM-flavour chosen')
         end select
      end if
!     make some noise so that we know we are leaving
      print*, ""
      call print_date("Leaving eom_eigval_driver routine")
      return

    end subroutine 

    subroutine eom_ccsd_eigenvalue_ip(results_t,exa_input,interm_t,int_t)
       use talsh_sigma_eomccsd

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
!     workspace
      type(talsh_davidson_workspace_tens_t) :: workspace_t
      type(talsh_projectors_tens_t) :: projector_t
!     error code
      integer :: ierr

      call print_date("Entering EOM-IP energy calculation")

      call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, exa_input%eom_flavour, 2)
      call talsh_eomccsd_ip_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
      call talsh_eomccsd_ip_trial_create(workspace_t)

!     set pointers to sigma vector subroutines 
      interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ip
      interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ip

      call talsh_davidson_driver(results_t%value_vectors_t,exa_input,int_t,interm_t, workspace_t, projector_t)
      call print_date("Finishing EOM-IP eigenvalue calculations")

      call talsh_davidson_workspace_destroy(workspace_t)

    end subroutine eom_ccsd_eigenvalue_ip


    subroutine eom_ccsd_eigenvalue_ea(results_t,exa_input,interm_t,int_t)
       use talsh_sigma_eomccsd

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
!     workspace
      type(talsh_davidson_workspace_tens_t) :: workspace_t
      type(talsh_projectors_tens_t) :: projector_t
!     error code
      integer :: ierr

      call print_date("Entering EOM-EA energy calculation")

      call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, exa_input%eom_flavour, 2)
!     call talsh_eomccsd_ea_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
!     call talsh_eomccsd_ea_trial_create(workspace_t)

!     set pointers to sigma vector subroutines 
!      interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ea
!      interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ea

      call talsh_davidson_driver(results_t%value_vectors_t,exa_input,int_t,interm_t, workspace_t, projector_t)
      call print_date("Finishing EOM-EA eigenvalue calculations")

      call talsh_davidson_workspace_destroy(workspace_t)

    end subroutine eom_ccsd_eigenvalue_ea



    subroutine eom_ccsd_eigenvalue_ee(exa_input,workspace_t,projector_t,interm_t,t2_tensor,&
     & int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)
       use talsh_sigma_eomccsd, only : sigma_eomccsd_rhs_ee, sigma_eomccsd_lhs_ee

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(talsh_projectors_tens_t), intent(inout) :: projector_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: l2_tensor,t1_tensor
!     tensor dims
      integer :: nocc, nvir
      integer, dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
      integer, dimension(4) :: oooo_dims, ooov_dims, ovoo_dims, oovv_dims, vvoo_dims, vovo_dims, vovv_dims, vvvo_dims, vvvv_dims
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr

      CCD = exa_input%ccd

!************************
!     construct tensors *
!************************

      nocc = exa_input%nocc
      nvir = exa_input%nvir

!     set tensor dims
      oo_dims = (/nocc,nocc/)
      ov_dims = (/nocc,nvir/)
      vo_dims = (/nvir,nocc/)
      vv_dims = (/nvir,nvir/)
      oooo_dims = (/nocc,nocc,nocc,nocc/)
      ooov_dims = (/nocc,nocc,nocc,nvir/)
      ovoo_dims = (/nocc,nvir,nocc,nocc/)
      if (exa_input%eom_properties) oovv_dims = (/nocc,nocc,nvir,nvir/)
      vvoo_dims = (/nvir,nvir,nocc,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      vovv_dims = (/nvir,nocc,nvir,nvir/)
      vvvo_dims = (/nvir,nvir,nvir,nocc/)
      vvvv_dims = (/nvir,nvir,nvir,nvir/)

!     create all fixed intermediates
      ierr=talsh_tensor_construct(interm_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(interm_t%goo,C8,oo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(interm_t%gvv,C8,vv_dims,init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(interm_t%w_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%w_ooov,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%w_vovv,C8,vovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(interm_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

      call get_wbar_vovo(interm_t%wbar_vovo,t2_tensor,int_t)
      if (.not.CCD) then
        call get_fbar_oo(interm_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_fbar_ov(interm_t%fbar_ov,fov_tensor,CCD,int_t%oovv,t1_tensor)
        call get_fbar_vv(interm_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_w_oooo (interm_t%w_oooo,int_t,t2_tensor,t1_tensor)
        call get_w_ooov (interm_t%w_ooov,int_t,t1_tensor)
        call get_w_ovoo (interm_t%w_ovoo,int_t,interm_t%fbar_ov,interm_t%w_oooo,interm_t%wbar_vovo,t2_tensor,t1_tensor)
        call get_w_vovo (interm_t%w_vovo,int_t,t2_tensor,t1_tensor)
        call get_w_vovv (interm_t%w_vovv,int_t,t1_tensor)
        call get_w_vvvv (interm_t%w_vvvv,int_t,t2_tensor,t1_tensor)
        call get_w_vvvo (interm_t%w_vvvo,int_t,interm_t%w_vvvv,interm_t%wbar_vovo,interm_t%fbar_ov,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(interm_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(interm_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(interm_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_oooo (interm_t%w_oooo,int_t,t2_tensor)
        call get_w_ooov (interm_t%w_ooov,int_t)
        call get_w_ovoo (interm_t%w_ovoo,int_t,interm_t%fbar_ov,interm_t%w_oooo,interm_t%wbar_vovo,t2_tensor)
        call get_w_vovo (interm_t%w_vovo,int_t,t2_tensor)
        call get_w_vovv (interm_t%w_vovv,int_t)
        call get_w_vvvv (interm_t%w_vvvv,int_t,t2_tensor)
        call get_w_vvvo (interm_t%w_vvvo,int_t,interm_t%w_vvvv,interm_t%wbar_vovo,interm_t%fbar_ov,t2_tensor)
      end if
      ierr=talsh_tensor_destruct(interm_t%wbar_vovo)

!***********************************
!     calculate right sigma vector *
!***********************************

      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, exa_input%eom_flavour)
      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_intermediates_trial_create(workspace_t, exa_input%eom_flavour)

!     call talsh_davidson_driver(exa_input,int_t,interm_t,t2_tensor, workspace_t, projector_t)
      call print_date("Calculated EOM-EE RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

!     create mixed type intermediates
      if (exa_input%eom_properties) then
        call get_goo_lambda(interm_t%goo,l2_tensor,t2_tensor)
        call get_gvv_lambda(interm_t%gvv,l2_tensor,t2_tensor)
      end if

!     if (exa_input%eom_properties) then
!       if (.not.CCD) then
!         call sigma_left_ee(exa_input,int_t,interm_t,t1_tensor)
!       else
!         call sigma_left_ee(exa_input,int_t,intermediates_t)
!       end if
!       call print_date("Calculated EOM-EE LHS")
!     end if

!**************
!     cleanup *
!**************

      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(interm_t%goo)
        ierr=talsh_tensor_destruct(interm_t%gvv)
      end if
      ierr=talsh_tensor_destruct(interm_t%fbar_oo)
      ierr=talsh_tensor_destruct(interm_t%fbar_ov)
      ierr=talsh_tensor_destruct(interm_t%fbar_vv)
      ierr=talsh_tensor_destruct(interm_t%w_oooo)
      ierr=talsh_tensor_destruct(interm_t%w_ooov)
      ierr=talsh_tensor_destruct(interm_t%w_ovoo)
      ierr=talsh_tensor_destruct(interm_t%w_vovo)
      ierr=talsh_tensor_destruct(interm_t%w_vovv)
      ierr=talsh_tensor_destruct(interm_t%w_vvvo)
      ierr=talsh_tensor_destruct(interm_t%w_vvvv)

    end subroutine eom_ccsd_eigenvalue_ee

end module
