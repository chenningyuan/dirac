module talsh_eom

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates
  use talsh_common_routines, only : tensor_norm2, print_tensor
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

    subroutine talsh_steom_eigval_driver (results_t,eom_config,eps_occ,eps_vir,interm_t,int_t)

!     input variables
      type(eom_config_t),     intent(in)             :: eom_config
      type(talsh_intg_tens_t), intent(inout)           :: int_t
      real(8),               intent(in)              :: eps_occ(:), eps_vir(:)
      type(talsh_intermediates_tens_t), intent(inout)  :: interm_t
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t

!     local variables
      type(talsh_eom_results_tens_t) :: eomip_results_t, eomea_results_t

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
      ! call talsh_davidson_create_results(eomip_results_t)
      ! call talsh_davidson_create_results(eomea_results_t)

      call print_date("Setting up RHS EOM-IP run")
      call eom_ccsd_eigenvalue_ip(eomip_results_t,eom_config%nroots,interm_t,int_t)
      call print_date("Done with  RHS EOM-IP run")

      call print_date("Setting up RHS EOM-EA run")
      call eom_ccsd_eigenvalue_ea(eomea_results_t,eom_config%nroots,interm_t,int_t)
      call print_date("Done with  RHS EOM-EA run")

      if (eom_config%do_expval) then
         call talsh_davidson_setup(unit_output=6, &
                                   solve_right=.false.,&
                                   solve_left=.true.,&
                                   symmetric=.false., &
                                   complex_mode=.true.)
         call print_date("Setting up LHS EOM-IP run")
         call eom_ccsd_eigenvalue_ip(eomip_results_t,eom_config%nroots,interm_t,int_t)
         call print_date("Done with  LHS EOM-IP run")

         call print_date("Setting up LHS EOM-EA run")
         call eom_ccsd_eigenvalue_ea(eomea_results_t,eom_config%nroots,interm_t,int_t)
         call print_date("Done with  LHS EOM-EA run")
      end if

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
    subroutine talsh_eom_eigval_driver (results_t,eom_config,interm_t,int_t)
! aspg: for now, use relcc_cfg to get the infromation on the *CCDIAG variables, like whether or not to do overlap etc
! perhaps these variables have to be moved around. this will definitely not work with a distributed version since the 
! variables need to be broadcast, and that requires a call to relcc_sync_cw(), which is done in relccsd but obviously not in exacorr
      use relcc_cfg

!     input variables
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
      type(eom_config_t),     intent(in)             :: eom_config
      type(talsh_intg_tens_t), intent(inout)         :: int_t
      type(talsh_intermediates_tens_t), intent(inout):: interm_t

      write(*,*) ""
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A)') "***************************** Equation of Motion ******************************"
      write(*,'(A)') "*******************************************************************************"
      write(*,*) ""

! aspg : need to get from hardcoded options to getting them from input (*CCDIAG) on the calls to talsh_davidson_setup
      call talsh_davidson_setup(unit_output=6, &
                                convergence_threshold=relcc_mfd_convergence_threshold, &
                                solve_right=.true.,&
                                solve_left=.false.,&
                                symmetric=.false., &
                                overlap_sorting=relcc_mfd_overlap_sorting,&
                                max_subspace_size=relcc_mfd_max_subspace_size,&
                                max_iterations=relcc_mfd_max_iterations,&
                                verbose=relcc_mfd_verbose,&
                                complex_mode=.true.)

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      ! print *,'debug: eom_input%flavor is__',eom_config%flavor,'___'
      ! print *,'debug: eom_input%flavor is__',eom_config%flavor(1:2),'___'
      select case(eom_config%flavor)
         case("IP")
            call eom_ccsd_eigenvalue_ip(results_t,eom_config%nroots,interm_t,int_t)
         case("EA")
            call eom_ccsd_eigenvalue_ea(results_t,eom_config%nroots,interm_t,int_t)
         case("EE")
            call eom_ccsd_eigenvalue_ee(results_t,eom_config%nroots,interm_t,int_t)
            ! call eom_ccsd_eigenvalue_ee(results_t,eom_config%nroots,interm_t,int_t,'L')  ! test left equation
         case default
            call quit ('Error: No valid EOM-flavour chosen')
      end select

      if (eom_config%do_expval) then  
         call talsh_davidson_setup(unit_output=6, &
                                   solve_right=.false.,&
                                   solve_left=.true.,&
                                   symmetric=.false., &
!                                  overlap_sorting=.true.,&
                                   complex_mode=.true.)
         select case(eom_config%flavor)
            case("IP")
               call eom_ccsd_eigenvalue_ip(results_t,eom_config%nroots,interm_t,int_t)
            case("EA")
               call eom_ccsd_eigenvalue_ea(results_t,eom_config%nroots,interm_t,int_t)
            case("EE")
               call eom_ccsd_eigenvalue_ee(results_t,eom_config%nroots,interm_t,int_t)
            case default
               call quit ('Error: No valid EOM-flavour chosen')
         end select
      end if
!     make some noise so that we know we are leaving
      print*, ""
      call print_date("Leaving eom_eigval_driver routine")
      return

    end subroutine 

    subroutine eom_ccsd_eigenvalue_ip(results_t, nroots, interm_t, int_t)
       use talsh_sigma_eomccsd

!     input variables
      integer, intent(in) :: nroots
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
!     workspace
      type(talsh_davidson_workspace_tens_t) :: workspace_t
      type(talsh_projectors_tens_t) :: projector_t
!     error code
      integer :: ierr
      character(len=2) :: flavor = "IP"

      call print_date("Entering EOM-IP energy calculation")

      call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2)
      call talsh_eomccsd_ip_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
      call talsh_eomccsd_ip_trial_create(workspace_t)
      call talsh_eomccsd_projector_create_kill_repeated_indexes(projector_t,workspace_t)

!     set pointers to sigma vector subroutines 
      interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ip
      interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ip

      call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t)
      call print_date("Finishing EOM-IP eigenvalue calculation")

      call talsh_davidson_workspace_destroy(workspace_t)

    end subroutine eom_ccsd_eigenvalue_ip


    subroutine eom_ccsd_eigenvalue_ea(results_t, nroots, interm_t, int_t)
       use talsh_sigma_eomccsd

!     input variables
      integer, intent(in) :: nroots
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
!     workspace
      type(talsh_davidson_workspace_tens_t) :: workspace_t
      type(talsh_projectors_tens_t) :: projector_t
!     error code
      integer :: ierr
      character(len=2) :: flavor = "EA"

      call print_date("Entering EOM-EA eigenvalue calculation")

      call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2)
      call talsh_eomccsd_ea_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv)
      call talsh_eomccsd_ee_trial_create(workspace_t)

!     set pointers to sigma vector subroutines 
      interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ea
      interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ea

      call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t)
      call print_date("Finishing EOM-EA eigenvalue calculation")

      call talsh_davidson_workspace_destroy(workspace_t)

    end subroutine eom_ccsd_eigenvalue_ea

    subroutine talsh_eomccsd_projector_create_kill_repeated_indexes(projector_t,workspace_t)
       type(talsh_projectors_tens_t), intent(inout) :: projector_t
       type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
       
       logical :: debug = .false.
       integer :: p, q, r, s, ierr, pmax, qmax, rmax, smax
       complex(8), pointer, contiguous:: proj1_tens(:,:), proj2_tens(:,:,:,:)
       type(C_PTR):: proj1_body_p, proj2_body_p
            
      ierr=talsh_tensor_construct(projector_t%p1,C8,workspace_t%tens1_dims(1:workspace_t%tens1_rank),init_val=ONE)
      ierr=talsh_tensor_construct(projector_t%p2,C8,workspace_t%tens2_dims(1:workspace_t%tens2_rank),init_val=ONE)

      projector_t%p1_rank = workspace_t%tens1_rank
      projector_t%p2_rank = workspace_t%tens2_rank

      allocate(projector_t%p1_dims(projector_t%p1_rank))
      allocate(projector_t%p2_dims(projector_t%p2_rank))

      projector_t%p1_dims(1:projector_t%p1_rank) = workspace_t%tens1_dims(1:workspace_t%tens1_rank)
      projector_t%p2_dims(1:projector_t%p2_rank) = workspace_t%tens2_dims(1:workspace_t%tens2_rank)
                  
      pmax  = projector_t%p2_dims(1) 
      qmax  = projector_t%p2_dims(2) 
      rmax  = projector_t%p2_dims(3) 
      smax  = projector_t%p2_dims(4) 

      ierr = talsh_tensor_get_body_access(projector_t%p1,proj1_body_p,C8,0,DEV_HOST)
      ierr = talsh_tensor_get_body_access(projector_t%p2,proj2_body_p,C8,0,DEV_HOST)

      call c_f_pointer(proj1_body_p,proj1_tens,projector_t%p1_dims(1:projector_t%p1_rank))
      call c_f_pointer(proj2_body_p,proj2_tens,projector_t%p2_dims(1:projector_t%p2_rank))

! this will create a projector that will zero out any set of repeated indexes
      if (rmax.eq.smax) then
         do s = 1, smax
            proj2_tens(:,:,s,s) = ZERO
         end do
     end if
     if (pmax.eq.qmax) then
         do p = 1, pmax
            proj2_tens(p,p,:,:) = ZERO
         end do
      end if


      if (debug) then
          write (6,*) 'debug >>> in projector_create_kill_repeated_indexes, p1 >>>'
          call print_tensor(projector_t%p1, 1.0d-18, "p1")
          write (6,*) 'debug <<< in projector_create_kill_repeated_indexes, p1 <<<'
          write (6,*) 'debug >>> in projector_create_kill_repeated_indexes, p2 >>>'
          call print_tensor(projector_t%p2, 1.0d-18, "p2")
          write (6,*) 'debug <<< in projector_create_kill_repeated_indexes, p2 <<<'
      end if


    end subroutine

    subroutine eom_ccsd_eigenvalue_ee(results_t, nroots, interm_t, int_t, left_right)
       use talsh_sigma_eomccsd

!     input variables
      integer, intent(in) :: nroots
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      character*1, OPTIONAL                              :: left_right
      type(talsh_eom_results_tens_t), intent(inout)  :: results_t
!     workspace
      type(talsh_davidson_workspace_tens_t) :: workspace_t
      type(talsh_projectors_tens_t) :: projector_t
!     error code
      integer :: ierr
      character(len=2) :: flavor = "EE"

      call print_date("Entering EOM-EE eigenvalue calculation")

      call talsh_davidson_workspace_create(workspace_t, nroots, interm_t%nocc, interm_t%nvir, flavor, 2, left_right)
      call talsh_eomccsd_ee_Hdiag_create(workspace_t, interm_t%Fbar_oo, interm_t%Fbar_vv, left_right)
      call talsh_eomccsd_ee_trial_create(workspace_t)
      call talsh_eomccsd_projector_create_kill_repeated_indexes(projector_t,workspace_t)

!     set pointers to sigma vector subroutines 
      interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ee
      interm_t%sigma_lhs_p => talsh_sigma_eomccsd_lhs_ee

      call talsh_davidson_driver(results_t%value_vectors_t,int_t,interm_t, workspace_t, projector_t)
      call print_date("Finishing EOM-EE eigenvalue calculation")

      call talsh_davidson_workspace_destroy(workspace_t)

    end subroutine eom_ccsd_eigenvalue_ee

end module
