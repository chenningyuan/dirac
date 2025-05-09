module talsh_QR_TPA
!This module contains the routines needed to calculate CC-CI(EOM-CC) response property


    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global
    use talsh_eom_rspeq_solver
    use talsh_ccrsp_common
    use talsh_eom
    use talsh_eom_LR, only: get_tF
    use talsh_cc_QR, only: get_lXAY, get_Gtt, get_tYFX

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public talsh_QR_TPA_driver

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine talsh_QR_TPA_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,eom_results_t)
!This subroutine drives the coupled-cluster two-photon absorption calculation 
!Written by Xiang Yuan June 2023, Lille, France


    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    ! integer                                     :: rr ! rank of response 
    type(talsh_eom_results_tens_t) :: eom_results_t
    
!   Wave-function tensor needed to construct the response function
    logical          :: CC_CC
    

    if (exa_input%wf_type == 1) then 
        CC_CC = .true.
        call print_date('Activate CC-CC(LR-CC) TPA')
    else 
        CC_CC = .false.
        call print_date('Activate CC-CI(EOM-CC) TPA')
    end if 

    call get_TPA(exa_input,interm_t,int_t,l1_tensor,l2_tensor,eom_results_t,CC_CC)


end subroutine talsh_QR_TPA_driver



!------------------------------------------------------------------------------------------------------------------

subroutine get_TPA(exa_input,interm_t,int_t,l1_tensor,l2_tensor,eom_results_t,CC_CC)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    type(talsh_eom_results_tens_t),intent(inout):: eom_results_t
    logical, intent(in)                         :: CC_CC

!   Wave-function tensor needed to construct the response function
    integer                       :: f_state  ! target TPA excited state
    complex(kind=8)               :: omega_fin
    complex(kind=8),allocatable   :: omega(:), omega_min(:)
    complex(kind=8)               :: qr_tran_left_AB, qr_tran_left_CD
    complex(kind=8)               :: qr_tran_right_AB, qr_tran_right_CD
    complex(kind=8)               :: qr_tpa
    complex(kind=8),allocatable   :: qr_tpas(:,:,:,:,:)
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR
    type(talsh_comm_tens_t), allocatable :: Bvec_tensors(:), Bvec_tensors_min_w(:)
    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:), Interm_t_tensors_B(:)
    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_C(:), Interm_t_tensors_D(:)

    type(talsh_comm_tens_t)       :: FR_f_tensor
    type(talsh_comm_tens_t)       :: barM_f_tensor

    integer                       :: N_prop_A, N_prop_B, N_prop_C, N_prop_D
    integer                       :: N_oper_tpa
    integer                       :: N_omega1, N_omega2
    integer                       :: i, j, k, l, n, m
    integer, allocatable          :: Prop_A(:), Prop_B(:), Prop_C(:), Prop_D(:)

    type(C_PTR) :: body_p
    complex(8), pointer           ::prop_tens(:,:,:,:)
    integer                       :: tens_rank
    integer                       :: ierr
    logical                       :: exchange
    integer                       :: A, B, C, D
    complex(kind=8)               :: biorth_sign_f

! aspg: this is now only one state, generalize to many states, once we deal the issues with biorthogonality for degenerate states
    f_state = exa_input%TPA_state
    omega_fin = eom_results_t%value_vectors_t%values_rhs(f_state,1)

    Prop_A = exa_input%Prop_A
    Prop_B = exa_input%Prop_B
    Prop_C = exa_input%Prop_C
    Prop_D = exa_input%Prop_D

    biorth_sign_f = ONE

    N_oper_tpa = 4
    N_prop_A = size(Prop_A)
    N_prop_B = size(Prop_B)
    N_prop_C = size(Prop_C)
    N_prop_D = size(Prop_D)


    N_omega1 = exa_input%N_omega_qr1

! aspg, 13/02/2024:
!       the allocations below, and the use of temporary omega and omega_min variables 
!       needs to be revised when dealing with more than one frequency. 
!       in the loops below there's a hardcoded access to elements (1) and (2) of 
!       omega() and omega_min() when these are used as arguments to the LR solvers.
!       this means that either we don't need 2*N_omega1 elements, or that the indexing
!       should be changed for the arguments.
!       whatever the case, prior to adding a factor of two in front of N_omega1  was 
!       one of the reasons making the TPA code segfault with heap corruption, since
!       in the loop 2 elements were being access, whereas only 1 was allocated when 
!       N_omega1 = 1...
!   
    allocate(omega(2*N_omega1))
    allocate(omega_min(2*N_omega1))

    omega = ZERO
    omega_min = omega 

    allocate(Interm_t_tensors_A(N_prop_A))
    allocate(Interm_t_tensors_B(N_prop_B))
    allocate(Interm_t_tensors_C(N_prop_C))
    allocate(Interm_t_tensors_D(N_prop_D))


    allocate(Bvec_tensors(N_oper_tpa))
    allocate(Tprop_tensors(N_oper_tpa))
    allocate(Tprop_tensors_min_w(N_oper_tpa))
    allocate(Tprop_bar_tensors(N_oper_tpa))
    allocate(Tprop_bar_tensors_min_w(N_oper_tpa))


    allocate(qr_tpas(N_prop_D,N_prop_C,N_prop_B,N_prop_A,N_omega1))

    call print_date('Entering CC TPA module')
    print *, ' '
    if (CC_CC) then
       print *, ' - Following QR-CC  formalism'
    else
       print *, ' - Following QR-EOM formalism'
    end if
    print *, ' - Final state frequency w_f : ', omega_fin, ' (ground state -> state ',f_state,')'
    print *, ' '

    
    do  i=1, N_prop_A
        m = Prop_A(i)
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_A(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC)
    end do
    
    
    do  i=1, N_prop_B
        m = Prop_B(i)
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_B(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC)
    end do


    do  i=1, N_prop_C
        m = Prop_C(i)
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_C(i),&
                                  int_t, l1_tensor, l2_tensor, CC_CC)
    end do


    do  i=1, N_prop_D
        m = Prop_D(i)
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_D(i),&
                                  int_t, l1_tensor, l2_tensor, CC_CC)
    end do


    call print_date('Constructed intermediates needed in QR-based CC TPA')

! aspg, first step to make this a bit more general, change hardcoded integers to variable
! names (A=perturbation A etc)
    A = 1
    B = 2
    C = 3
    D = 4

! solve excited state multipliers for omega_fin in the case of CC_CC
    if (CC_CC) then
!         solve for auxiliary excited state multipliers targeting omega_fin (w_f)
            call print_date('Solving auxiliary excited state multipliers for w_f')
            !FR_f_tensor will be used as RHS for linear system to get barM_f
            ierr=talsh_tensor_construct(FR_f_tensor%t1,C8,&
                       (/exa_input%nocc,exa_input%nvir/), &
                       init_val=ZERO)
            ierr=talsh_tensor_construct(FR_f_tensor%t2,C8,&
                       (/exa_input%nocc,exa_input%nocc,exa_input%nvir,exa_input%nvir/), &
                       init_val=ZERO)
            call get_tF(FR_f_tensor%t1, &
                        FR_f_tensor%t2, &
                        eom_results_t%value_vectors_t%rhs1(f_state), &
                        eom_results_t%value_vectors_t%rhs2(f_state), &
                        l1_tensor, &
                        l2_tensor, &
                        interm_t,  &
                        int_t,     &
                        exa_input%nocc, &
                        exa_input%nvir)
            print *, 'DEBUG, norm FR_f',tensor_norm2(FR_f_tensor%t1),tensor_norm2(FR_f_tensor%t2)

            call print_lr_setup("L"," ",Interm_t_tensors_A(1),-1,omega_fin)
!            call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensor, &
!                                  omega_fin, FR_f_tensor, barM_f_tensor)
            call solve_cc_lr_left(exa_input, interm_t, int_t, &
                                  omega_fin, FR_f_tensor, barM_f_tensor)

            call print_date('Done with auxiliary excited state multipliers for w_f')
            print *, 'DEBUG, norm barM_f',tensor_norm2(barM_f_tensor%t1),tensor_norm2(barM_f_tensor%t2)

            ierr=talsh_tensor_destruct(FR_f_tensor%t1)
            ierr=talsh_tensor_destruct(FR_f_tensor%t2)
        end if

! now go to loop over frequencies 
    do i=1, N_omega1
        ! omega2 = omega_fin - omega1

        omega(A) = dcmplx( real(0.5d0*omega_fin), aimag(0.5d0*omega_fin))
        omega(B) = dcmplx( real(0.5d0*omega_fin),-aimag(0.5d0*omega_fin))

        omega_min(A) = dcmplx(-real(omega(A)),aimag(omega(A)))
        omega_min(B) = dcmplx(-real(omega(B)),aimag(omega(B)))

        print *, ' '
        print *, ' - Employing the following frequencies in TPA calculation (in atomic units)'
        print *, '   w  (=  w_A) : ', omega(A)      ! +w  + i\gamma
        print *, "   w' (=  w_B) : ", omega(B)      ! +w' - i\gamma
        print *, '  -w  (= -w_A) : ', omega_min(A)  ! -w  + i\gamma 
        print *, "  -w' (= -w_B) : ", omega_min(B)  ! -w' - i\gamma
        print *, ' '

! note that in the solution of linear systems below we solve for perturbed
! amplitudes and multipliers (for pertubation P) in the order:
! (1) t^P(w'), tbar^P(-w') (for either QR-CC and QR-EOMCC)
! (2) t^P(-w'), tbar^P(w') (for QR-CC only)
! instead of
! (1') t^P(w'), tbar^P(w') 
! (2') t^P(-w'), tbar^P(-w') 
! due to the fact that for QR-EOMCC 2-photon transition moments are built using either 
! perturbed amplitudes or petruberd multipliers, whereas for QR-CC left transition moments
! both amplitudes and multipliers for each frequency and perturbation are needed. 

        do n=1, N_prop_A
            ! solve for tA(w'), bartA(-w')
            call solve_left_right_rsp_driver(exa_input,                &
                                            int_t,                     &
                                            interm_t,                  &
                                            Interm_t_tensors_A(n),     &
                                            l1_tensor,l2_tensor,       &
                                            Bvec_tensors(A),           &
                                            omega(B),                  & 
                                            Tprop_tensors(A),          &  ! tA(+w' - i\gamma)
                                            omega_min(B),              &  
                                            Tprop_bar_tensors_min_w(A),&  ! bartA(-w' - i\gamma)
                                            CC_CC)
            if (CC_CC) then
                ! solve for tA(w), bartA(-w) 
                call solve_left_right_rsp_driver(exa_input,            &
                                            int_t,                     &
                                            interm_t,                  &
                                            Interm_t_tensors_A(n),     &
                                            l1_tensor,l2_tensor,       &
                                            Bvec_tensors(A),           &
                                            omega(A),                  & 
                                            Tprop_tensors_min_w(A),    &  ! tA(+w  + i\gamma)
                                            omega_min(A),              &  
                                            Tprop_bar_tensors(A),      &  ! bartA(-w  + i\gamma)
                                            CC_CC)
           end if
            
            do m=1, N_prop_B
                ! solve for tB(w), bartB(-w)
                call solve_left_right_rsp_driver(exa_input,                &
                                                int_t,                     &
                                                interm_t,                  &
                                                Interm_t_tensors_B(m),     &
                                                l1_tensor,                 &
                                                l2_tensor,                 &
                                                Bvec_tensors(B),           &
                                                omega(A),                  &
                                                Tprop_tensors(B),          &  ! tB(+w  + i\gamma)
                                                omega_min(A),              &
                                                Tprop_bar_tensors_min_w(B),&  ! bartB(-w  + i\gamma)
                                                CC_CC)
                if (CC_CC) then
                    ! solve for tB(w'), bartB(-w')
                    call solve_left_right_rsp_driver(exa_input,        &
                                            int_t,                     &
                                            interm_t,                  &
                                            Interm_t_tensors_B(m),     &
                                            l1_tensor,l2_tensor,       &
                                            Bvec_tensors(B),           &
                                            omega(B),                  & 
                                            Tprop_tensors_min_w(B),    &  ! tB(+w' - i\gamma)
                                            omega_min(B),              &  
                                            Tprop_bar_tensors(B),      &  ! bartB(-w' - i\gamma)
                                            CC_CC)
                end if

                do k=1, N_prop_C
                    call solve_left_right_rsp_driver(exa_input,&
                                                     int_t,    &
                                                     interm_t, &
                                                     Interm_t_tensors_C(k),&
                                                     l1_tensor, &
                                                     l2_tensor, &
                                                     Bvec_tensors(C),&
                                                     omega(B), &
                                                     Tprop_tensors(C), &   ! tC(+w' - i\gamma)
                                                     omega_min(B), &
                                                     Tprop_bar_tensors_min_w(C),& ! bartC(-w' - i\gamma)
                                                     CC_CC)
                    if (CC_CC) then
                        ! solve for tB(w'), bartB(-w')
                        call solve_left_right_rsp_driver(exa_input,    &
                                            int_t,                     &
                                            interm_t,                  &
                                            Interm_t_tensors_C(k),     &
                                            l1_tensor,l2_tensor,       &
                                            Bvec_tensors(C),           &
                                            omega(A),                  & 
                                            Tprop_tensors_min_w(C),    &  ! tC(+w + i\gamma)
                                            omega_min(A),              &  
                                            Tprop_bar_tensors(C),      &  ! bartC(-w + i\gamma)
                                            CC_CC)
                    end if
                
                    do l=1, N_prop_D
                        call solve_left_right_rsp_driver(exa_input, &
                                                         int_t, &
                                                         interm_t, &
                                                         Interm_t_tensors_D(l),&
                                                         l1_tensor, &
                                                         l2_tensor, &
                                                         Bvec_tensors(D), &
                                                         omega(A), &
                                                         Tprop_tensors(D), &  ! tD(+w  + i\gamma)
                                                         omega_min(A), &
                                                         Tprop_bar_tensors_min_w(D), & ! bartD(-w  + i\gamma)
                                                         CC_CC)

                        if (CC_CC) then
                            ! solve for tB(w'), bartB(-w')
                            call solve_left_right_rsp_driver(exa_input,&
                                            int_t,                     &
                                            interm_t,                  &
                                            Interm_t_tensors_D(l),     &
                                            l1_tensor,l2_tensor,       &
                                            Bvec_tensors(D),           &
                                            omega(B),                  & 
                                            Tprop_tensors_min_w(D),    &  ! tD(+w' - i\gamma)
                                            omega_min(B),              &  
                                            Tprop_bar_tensors(D),      &  ! bartD(-w' - i\gamma)
                                            CC_CC)
                        end if

                        call print_date('Solved All 1st response equations')


                        call check_biorthognal(eom_results_t,f_state,biorth_sign_f)

                        qr_tran_left_AB  = ZERO
                        qr_tran_left_CD  = ZERO
                        qr_tran_right_AB = ZERO
                        qr_tran_right_CD = ZERO

                        if (CC_CC) then
                            call get_cc_tran_left(exa_input, &
                                                   int_t,     &
                                                   interm_t,  &
                                                   Interm_t_tensors_A(n), & ! X
                                                   Interm_t_tensors_B(m), & ! Y
                                                   l1_tensor, &
                                                   l2_tensor, &
                                                   Tprop_tensors_min_w(B), & ! tX (-w' - i\gamma)
                                                   Tprop_tensors_min_w(A), & ! tY (w + i\gamma)
                                                   Tprop_bar_tensors_min_w(B), & ! tX (-w' - i\gamma)
                                                   Tprop_bar_tensors_min_w(A), & ! tY (w + i\gamma)
                                                   barM_f_tensor, &
                 eom_results_t%value_vectors_t%rhs1(f_state), &
                 eom_results_t%value_vectors_t%rhs2(f_state), &
                                                   qr_tran_left_AB)

                            call get_cc_tran_left(exa_input, &
                                                   int_t,     &
                                                   interm_t,  &
                                                   Interm_t_tensors_C(k), & ! X
                                                   Interm_t_tensors_D(l), & ! Y
                                                   l1_tensor, &
                                                   l2_tensor, &
                                                   Tprop_tensors_min_w(D), & ! tX (-w' - i\gamma)
                                                   Tprop_tensors_min_w(C), & ! tY (w + i\gamma)
                                                   Tprop_bar_tensors_min_w(D), & ! tX (-w' - i\gamma)
                                                   Tprop_bar_tensors_min_w(C), & ! tY (w + i\gamma)
                                                   barM_f_tensor, &
                 eom_results_t%value_vectors_t%rhs1(f_state), &
                 eom_results_t%value_vectors_t%rhs2(f_state), &
                                                   qr_tran_left_CD)

                            call get_cc_tran_right(exa_input, &
                                                   int_t,     &
                                                   interm_t,  &
                                                   Interm_t_tensors_A(n), & ! X
                                                   Interm_t_tensors_B(m), & ! Y
                                                   l1_tensor, &
                                                   l2_tensor, &
                                                   Tprop_tensors(B), & ! tX (-w' - i\gamma)
                                                   Tprop_tensors(A), & ! tY (w + i\gamma)
                 eom_results_t%value_vectors_t%lhs1(f_state), &
                 eom_results_t%value_vectors_t%lhs2(f_state), &
                                                   qr_tran_right_AB)

                            call get_cc_tran_right(exa_input, &
                                                   int_t,     &
                                                   interm_t,  &
                                                   Interm_t_tensors_C(k), & ! X
                                                   Interm_t_tensors_D(l), & ! Y
                                                   l1_tensor, &
                                                   l2_tensor, &
                                                   Tprop_tensors(D), & ! tX (-w' - i\gamma)
                                                   Tprop_tensors(C), & ! tY (w + i\gamma)
                 eom_results_t%value_vectors_t%lhs1(f_state), &
                 eom_results_t%value_vectors_t%lhs2(f_state), &
                                                   qr_tran_right_CD)

                        else
                            call get_eom_tran_left(exa_input, &
                                                  interm_t,   &
                                                  Interm_t_tensors_A(n), &
                                                  Interm_t_tensors_B(m), & 
                                                  Tprop_bar_tensors_min_w(A), &
                                                  Tprop_bar_tensors_min_w(B), &
                                                  l1_tensor, &
                                                  l2_tensor, &
                                                  eom_results_t, &
                                                  qr_tran_left_AB, &
                                                  f_state)

                            call get_eom_tran_left(exa_input, &
                                                  interm_t,   &
                                                  Interm_t_tensors_C(k), &
                                                  Interm_t_tensors_D(l), & 
                                                  Tprop_bar_tensors_min_w(C), &
                                                  Tprop_bar_tensors_min_w(D),&
                                                  l1_tensor, &
                                                  l2_tensor, &
                                                  eom_results_t, &
                                                  qr_tran_left_CD, &
                                                  f_state)

                            call get_eom_tran_right(exa_input, &
                                                   interm_t,&
                                                   Interm_t_tensors_A(n), &
                                                   Interm_t_tensors_B(m), & 
                                                   Tprop_tensors(B), &
                                                   Tprop_tensors(A), &
                                                   l1_tensor, &
                                                   l2_tensor, &
                                                   eom_results_t, &
                                                   qr_tran_right_AB, &
                                                   f_state)

                            call get_eom_tran_right(exa_input, &
                                                   interm_t,&
                                                   Interm_t_tensors_C(k),&
                                                   Interm_t_tensors_D(l),& 
                                                   Tprop_tensors(D), &
                                                   Tprop_tensors(C),&
                                                   l1_tensor, &
                                                   l2_tensor, &
                                                   eom_results_t, &
                                                   qr_tran_right_CD, &
                                                   f_state)
                        end if
! aspg: workaround to deal with sign change due to biorthogonality issues 
                        qr_tran_left_AB = qr_tran_left_AB*biorth_sign_f
                        qr_tran_left_CD = qr_tran_left_CD*biorth_sign_f

                        print *, " - Two-photon transition moments"
                        print *, "   LHS T_AB: ", qr_tran_left_AB
                        print *, "   LHS T_CD: ", qr_tran_left_CD
                        print *, "   RHS T_AB: ", qr_tran_right_AB
                        print *, "   RHS T_CD: ", qr_tran_right_CD

                        call get_scattering_amplitudes (qr_tran_left_AB,  &
                                                        qr_tran_left_CD,  &
                                                        qr_tran_right_AB, &
                                                        qr_tran_right_CD, &
                                                        qr_tpa)
                        qr_tpas(l,k,m,n,i) = qr_tpa
                       
                        call destroy_Tprop_tensors(Tprop_tensors,       &
                                                   Tprop_tensors_min_w, &
                                                   D,                   &
                                                   Bvec_tensors,        &
                                                   Tprop_bar_tensors,   &
                                                   Tprop_bar_tensors_min_w)
                    end do
                    call destroy_Tprop_tensors(Tprop_tensors,       &
                                               Tprop_tensors_min_w, &
                                               C,                   &
                                               Bvec_tensors,        &
                                               Tprop_bar_tensors,   &
                                               Tprop_bar_tensors_min_w)
                end do
                call destroy_Tprop_tensors(Tprop_tensors,       &
                                           Tprop_tensors_min_w, &
                                           B,                   &
                                           Bvec_tensors,        &
                                           Tprop_bar_tensors,   &
                                           Tprop_bar_tensors_min_w)
            end do
            call destroy_Tprop_tensors(Tprop_tensors,       &
                                       Tprop_tensors_min_w, &
                                       A,                   &
                                       Bvec_tensors,        &
                                       Tprop_bar_tensors,   &
                                       Tprop_bar_tensors_min_w)
        end do

    end do


    do i=1, N_prop_A
        call destroy_Interm_t_tensors(Interm_t_tensors_A(i),CC_CC)
    end do

    do i=1, N_prop_B
        call destroy_Interm_t_tensors(Interm_t_tensors_B(i),CC_CC)
    end do
        
    do i=1, N_prop_C
        call destroy_Interm_t_tensors(Interm_t_tensors_C(i),CC_CC)
    end do

    do i=1, N_prop_D
        call destroy_Interm_t_tensors(Interm_t_tensors_D(i),CC_CC)
    end do

    if (CC_CC) then
        write(*,"(A54)"), "*************CC-CC(QR-CC) Two-photon absorption************"
    else
        write(*,"(A54)"), "*************CC-CI(EOM-CC) Two-photon absorption************"
    end if
    
    do n=1, N_omega1
        do m=1, N_prop_A
            do i=1, N_prop_B
                do j=1, N_prop_C
                    do k=1, N_prop_D  
                        write(*,"(F16.8,A,F16.8)"), real(qr_tpas(k,j,i,m,n)), ' \ ', aimag(qr_tpas(k,j,i,m,n))
                    end do
                end do 
            end do
        end do 
    end do
    write(*,"(A54)"), "**********************************************************"  

    DEALLOCATE(omega)
    DEALLOCATE(omega_min)
    DEALLOCATE(Tprop_tensors)
    DEALLOCATE(Tprop_tensors_min_w)
    DEALLOCATE(Tprop_bar_tensors)
    DEALLOCATE(Tprop_bar_tensors_min_w)
    DEALLOCATE(Interm_t_tensors_A)
    DEALLOCATE(Interm_t_tensors_B)
    DEALLOCATE(Interm_t_tensors_C)
    DEALLOCATE(Interm_t_tensors_D)
    DEALLOCATE(Bvec_tensors)
    DEALLOCATE(qr_tpas)

    call print_date('Leave QR-CC TPA') 


end subroutine get_TPA


!---------------------------------------------------------------------------------------------------------------------------------------

subroutine solve_left_right_rsp_driver(exa_input,int_t,interm_t,Interm_t_tensor,l1_tensor,l2_tensor, &
                    Bvec_tensor, omega1, Tprop_tensor, omega2, Tprop_bar_tensor, CC_CC)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_rsp_interm_tens_t),intent(inout) :: Interm_t_tensor
    type(talsh_comm_tens_t),intent(inout)       :: Bvec_tensor
    complex(kind=8),intent(in)                  :: omega1, omega2
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    type(talsh_comm_tens_t),intent(inout)       :: Tprop_tensor  !right_vector
    type(talsh_comm_tens_t),intent(inout)       :: Tprop_bar_tensor !left_vector
    logical, intent(in) :: CC_CC

    call print_lr_setup("R"," ",Interm_t_tensor,-1,omega1)
    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensor, omega1, Tprop_tensor)
    call print_lr_setup("L"," ",Interm_t_tensor,-1,omega2)
    if (CC_CC) then
!        call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensor, omega2, Bvec_tensor, Tprop_bar_tensor)
        call get_left_bvec_cccc(exa_input,int_t,interm_t,Interm_t_tensor, &
                    l1_tensor,l2_tensor, Tprop_tensor, Bvec_tensor)
        call solve_cc_lr_left(exa_input, interm_t, int_t, omega2, Bvec_tensor, Tprop_bar_tensor)
    else
!       call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensor, omega2, Bvec_tensor, Tprop_bar_tensor)
       call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensor,l1_tensor,l2_tensor, Bvec_tensor)
       call solve_eom_lr_left(exa_input, interm_t, int_t, omega2, Bvec_tensor, Tprop_bar_tensor)
    end if


end subroutine solve_left_right_rsp_driver

!------------------------------------------------------------------------------------------------------------------------------------------

subroutine get_cc_tran_right(exa_input, int_t,interm_t, &
                             X_operator, Y_operator, &
                             l1_tensor,l2_tensor, &
                             tX, tY, &
                             Ls_f, Ld_f, &
                             qr_tran_right_XY) !, f_state)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_rsp_interm_tens_t),intent(inout) :: X_operator
    type(talsh_rsp_interm_tens_t),intent(inout) :: Y_operator
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    type(talsh_tens_t), intent(inout)           :: Ls_f ! eigenvalue L1 final state
    type(talsh_tens_t), intent(inout)           :: Ld_f ! eigenvalue L2 final state
    type(talsh_comm_tens_t),intent(inout)       :: tX  ! X perturbed amplitude
    type(talsh_comm_tens_t),intent(inout)       :: tY  ! Y perturbed amplitude
!    type(talsh_eom_results_tens_t),intent(inout):: eom_results_t
    complex(kind=8),intent(inout)               :: qr_tran_right_XY ! 2-photon transition amplitude ^{CC}F_{XY}^{f0}

    type(talsh_comm_tens_t)                     :: LfMN_tensor

    integer(C_INT)     :: one_dims(1)
    integer(INTD), dimension(2)   :: ov_dims
    integer(INTD), dimension(4)   :: oovv_dims
    integer(INTD)  ::  ov_rank, oovv_rank

    integer(C_INT)          :: ierr
    type(talsh_tens_t)      :: r_tensor
    complex(8), pointer     :: r_tens(:)
    type(C_PTR)             :: body_p

    complex(8) :: tran1, tran2, tran
    logical :: debug_local = .true.

!
! calculate RHS 2-photon transition moment from QR for X, Y perturbation
!
! see eqs 5.51-5.57 of Christiansen, Hattig and Jorgensen, IJCQ 1998, doi: 10.1002/(SICI)1097-461X(1998)68:1<1::AID-QUA1>3.0.CO;2-Z
! or for a more straightfoward version, eq 3 in Faber and Coriani, JCTC 2018, doi: 10.1021/acs.jctc.8b01020 
!
! ^{CC}F_{XY}^{f0} = -L_f [  A^X * t^Y(\omega + i\gamma) 
!                          + A^Y * t^X(-\omega^{'} - i\gamma}) 
!                          + B * t^Y(\omega + i\gamma) * t^X(-\omega^{'} - i\gamma})]
!
!
    qr_tran_right_XY = 0

    ierr=talsh_tensor_construct(r_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(r_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,r_tens,one_dims)

    ierr=talsh_tensor_construct(LfMN_tensor%t1,C8,&
                       (/interm_t%nocc,interm_t%nvir/), &
                       init_val=ZERO)
    ierr=talsh_tensor_construct(LfMN_tensor%t2,C8,&
                       (/interm_t%nocc,interm_t%nocc,interm_t%nvir,interm_t%nvir/), &
                       init_val=ZERO)
! 1) form L_f * B * t^Y(\omega + i\gamma) * t^X(-\omega^{'} - i\gamma})
!
! recall we reuse Ft to obtain BtY, with special case (remove U1 vua is_qrcc=.true.) 
! AND at the same time replace l1, l2 by the LHS eigenvector (Ls_f, Ld_f) determined 
! previously

    call get_tF(LfMN_tensor%t1,LfMN_tensor%t2, tY%t1,tY%t2, Ls_f,Ld_f, &
                interm_t,int_t,interm_t%nocc, interm_t%nvir, is_qrcc=.true.)

! contract (L_f *B * tY)^m_e with tX^e_m for singles contribution and
! and      (L_f *B * tY)^m_e with tX^{ef}_{mn} for doubles contribution to 
! 2-photon transition moment 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("Tf0XY()+=LfBtY(m,e)*tX(e,m)",r_tensor,&
                               LfMN_tensor%t1,tX%t1)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("Tf0XY()+=LfBtY(m,n,e,f)*tX(e,f,m,n)",r_tensor,&
                               LfMN_tensor%t2,tX%t2,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(/,2x,A)') 'DEBUG in get_cc_tran_right'
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (1) LfBtXtY    contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    end if

!   add -LfBtXtY to transition moment
    qr_tran_right_XY = qr_tran_right_XY - r_tens(1)
!
! 2) form  L_f * A^Y * t^X
!
! add        (L_f AX)^m_e  to (MAB_tensor)^m_e
! add        (L_f AX)^m_e  to (MAB_tensor)^{mn}_{ef}

    ierr=talsh_tensor_init(LfMN_tensor%t1)
    ierr=talsh_tensor_init(LfMN_tensor%t2)

    call get_lXAY(LfMN_tensor, Y_operator%prop_bar_oo, Y_operator%prop_bar_vv, &
                   Y_operator%prop_bar_ov, interm_t%t1, interm_t%t2, Ls_f,Ld_f, &
                   interm_t%nocc, interm_t%nvir)

    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("Tf0XY()+=LfAY(m,e)*tX(e,m)",r_tensor,&
                               LfMN_tensor%t1,tX%t1)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("Tf0XY()+=LfAY(m,n,e,f)*tX(e,f,m,n)",r_tensor,&
                               LfMN_tensor%t2,tX%t2,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
         write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
               'DEBUG (2) LfAYtX     contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    endif

!   add -LfAYtX to transition moment
    qr_tran_right_XY = qr_tran_right_XY - r_tens(1)
!
! 3) form  L_f * A^X * t^Y
!
! add        (L_f AX)^m_e  to (MAB_tensor)^m_e
! add        (L_f AX)^m_e  to (MAB_tensor)^{mn}_{ef}

    ierr=talsh_tensor_init(LfMN_tensor%t1)
    ierr=talsh_tensor_init(LfMN_tensor%t2)

    call get_lXAY(LfMN_tensor, X_operator%prop_bar_oo, X_operator%prop_bar_vv, &
                   X_operator%prop_bar_ov, interm_t%t1, interm_t%t2, Ls_f,Ld_f, &
                   interm_t%nocc, interm_t%nvir)

    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("Tf0XY()+=LfAX(m,e)*tY(e,m)",r_tensor,&
                               LfMN_tensor%t1,tY%t1)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("Tf0XY()+=LfAX(m,n,e,f)*tY(e,f,m,n)",r_tensor,&
                               LfMN_tensor%t2,tY%t2,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (3) LfAXtY     contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
        write(*,'(2x,A,/)') 'DEBUG exiting get_cc_tran_right'
    endif

!   add -LAXtY to transition moment
    qr_tran_right_XY = qr_tran_right_XY - r_tens(1)
!
    ierr=talsh_tensor_destruct(LfMN_tensor%t1)
    ierr=talsh_tensor_destruct(LfMN_tensor%t2)
    ierr=talsh_tensor_destruct(r_tensor)

end subroutine get_cc_tran_right

!------------------------------------------------------------------------------------------------------------------------------------------

subroutine get_cc_tran_left(exa_input, int_t,interm_t, &
                             X_operator, Y_operator, &
                             l1_tensor,l2_tensor, &
                             tX, tY,  &
                             tbar_X, tbar_Y, &
                             barM_f, &
                             Rs_f, Rd_f, &
                             qr_tran_left_XY)


    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_rsp_interm_tens_t),intent(inout)    :: X_operator
    type(talsh_rsp_interm_tens_t),intent(inout)    :: Y_operator
    type(talsh_tens_t), intent(inout)              :: l2_tensor
    type(talsh_tens_t), intent(inout)              :: l1_tensor
    type(talsh_comm_tens_t),intent(inout)          :: tX  ! X perturbed amplitudes
    type(talsh_comm_tens_t),intent(inout)          :: tY  ! Y perturbed amplitudes
    type(talsh_comm_tens_t),intent(inout)          :: tbar_X ! X perturbed multipliers
    type(talsh_comm_tens_t),intent(inout)          :: tbar_Y ! Y perturbed amplitude
    type(talsh_comm_tens_t),intent(inout)          :: barM_f ! excited state multipliers
    type(talsh_tens_t), intent(inout)              :: Rs_f ! eigenvalue L1 final state
    type(talsh_tens_t), intent(inout)              :: Rd_f ! eigenvalue L2 final state
!   type(talsh_eom_results_tens_t),intent(inout):: eom_results_t
    complex(kind=8),intent(inout)                  :: qr_tran_left_XY ! 2-photon transition amplitude ^{CC}F_{XY}^{0f}

    type(talsh_comm_tens_t)                        :: MNRf_tensor

    integer(C_INT)     :: one_dims(1)
    integer(INTD), dimension(2)   :: ov_dims
    integer(INTD), dimension(4)   :: oovv_dims
    integer(INTD)  ::  ov_rank, oovv_rank

    integer(C_INT)          :: ierr
    type(talsh_tens_t)      :: r_tensor
    complex(8), pointer     :: r_tens(:)
    type(C_PTR)             :: body_p

    complex(8) :: tran1, tran2, tran, tmom_b_XY
    logical :: debug_local = .true.

!
! see eqs 5.51-5.57 of Christiansen, Hattig and Jorgensen, IJCQ 1998, doi: 10.1002/(SICI)1097-461X(1998)68:1<1::AID-QUA1>3.0.CO;2-Z
! or for a more straightfoward version, eq 11 in Faber and Coriani, JCTC 2018, doi: 10.1021/acs.jctc.8b01020 
!
! ^{CC}F_{XY}^{0f} = 
!                  - {
!                       \bar{t}^X(\omega^{'} - i\gamma}) *[ A^Y + B * t^Y(-\omega + i\gamma) ]
!                     + \bar{t}^Y(-\omega + i\gamma)     *[ A^X + B * t^X(\omega^{'} - i\gamma}) ]
!                     + F^X * t^Y(-\omega + i\gamma)
!                     + F^Y * t^X(\omega^{'} - i\gamma})
!                     + G * t^X(\omega^{'} - i\gamma}) * t^Y(-\omega + i\gamma) 
!                    } * R_f       (1)
!
!                  - \bar{M}_f(\omega_f) *
!                    [  A^X * t^Y(-\omega + i\gamma) 
!                     + A^Y * t^X(\omega^{'} - i\gamma}) 
!                     + B * t^Y(-\omega + i\gamma) * t^X(\omega^{'} - i\gamma})]   (2)
!
! remarks for term (2):
! - it needs excited state multipliers as well (eq 13), obtained from solving
!
! \bar{M}_f(\omega_f)(A + \omega_f I) = -FR_f
!
! with \omega_f the frequency for the transition to state f. The \bar{M}_f(\omega_f) tensor should
! be obtained outside of this routine, and pass as an argument
!
! - it can be formed by calling get_cc_tran_right(), provided we use 
! (a) t^X(\omega^{'} - i\gamma) instead of t^X(-\omega^{'} - i\gamma})
! (b) t^Y(-\omega + i\gamma) instead of t^Y(\omega + i\gamma)
!
    qr_tran_left_XY = 0

    ierr=talsh_tensor_construct(r_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(r_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,r_tens,one_dims)

    ierr=talsh_tensor_construct(MNRf_tensor%t1,C8,&
                       (/interm_t%nocc,interm_t%nvir/), &
                       init_val=ZERO)
    ierr=talsh_tensor_construct(MNRf_tensor%t2,C8,&
                       (/interm_t%nocc,interm_t%nocc,interm_t%nvir,interm_t%nvir/), &
                       init_val=ZERO)

! term (1) contributions
!
! 1. form (1.a.1): \bar{t}^X(\omega^{'} - i\gamma}) * A^Y * R_f  
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)
    call get_lXAY(MNRf_tensor, Y_operator%prop_bar_oo, Y_operator%prop_bar_vv, &
                   Y_operator%prop_bar_ov, interm_t%t1, interm_t%t2, &
                   tbar_X%t1, tbar_X%t2, interm_t%nocc, interm_t%nvir)

! contract (tbarX * AY)^m_e with Rf^e_m for singles contribution and
! and      (tbarX * AY)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=tbarXAY(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=tbarXAY(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(/,2x,A)') 'DEBUG in get_cc_tran_left'
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (1) tbarXAYRf  contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    endif

!   add -tbarXAYRf to transition moment
    qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! 2. form (1.a.2): \bar{t}^X(\omega^{'} - i\gamma}) * B * t^Y(-\omega + i\gamma) * R_f 
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)
! recall we reuse Ft to obtain BtY, with special case (remove U1 vua is_qrcc=.true.) AND 
! at the same time  replace l1, l2 by perburtbed l1, l2. which were calculated before (tbar_x)
    call get_tF(MNRf_tensor%t1,MNRf_tensor%t2, tY%t1,tY%t2, tbar_X%t1, tbar_X%t2, &
                interm_t,int_t,interm_t%nocc, interm_t%nvir, is_qrcc=.true.)

! contract (tbarX *B * tY)^m_e with Rf^e_m for singles contribution and
! and      (tbarX *B * tY)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=tbarXBtY(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=tbarXBtY(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (2) tbarXBtYRf contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    end if

!   add -tbarXBtYRf to transition moment
    qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! 3. form (1.b.1): \bar{t}^Y(-\omega + i\gamma) * A^X * R_f
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)
    call get_lXAY(MNRf_tensor, X_operator%prop_bar_oo, X_operator%prop_bar_vv, &
                   X_operator%prop_bar_ov, interm_t%t1, interm_t%t2, &
                   tbar_Y%t1, tbar_Y%t2, interm_t%nocc, interm_t%nvir)
! contract (-tbarY * AX)^m_e with Rf^e_m for singles contribution and
! and      (-tbarY * AX)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYAX(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYAX(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (3) tbarYAXRf  contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    endif

!   add -tbarYAXRf to transition moment
    qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! 4. form (1.b.2): \bar{t}^Y(-\omega + i\gamma) * B * t^X(\omega^{'} - i\gamma}) * R_f
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)
! recall we reuse Ft to obtain BtX, with special case (remove U1 vua is_qrcc=.true.) 
    call get_tF(MNRf_tensor%t1,MNRf_tensor%t2, tX%t1,tX%t2, tbar_Y%t1, tbar_Y%t2, &
                interm_t,int_t,interm_t%nocc, interm_t%nvir, is_qrcc=.true.)
!
! contract (-tbarY *B * tX)^m_e with Rf^e_m for singles contribution and
! and      (-tbarY *B * tX)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYBtX(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYBtX(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (4) tbarYBtXRf contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    end if

!    add -tbarYBtXRf to transition moment
     qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! 5. form (1.c.1): F^X * t^Y(-\omega + i\gamma) * R_f
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)
    call get_tYFX(MNRf_tensor,X_operator%prop_bar_ov,tY%t1,tY%t2,&
                  l1_tensor,l2_tensor,interm_t%nocc,interm_t%nvir)
! contract (FX tY)^m_e with Rf^e_m for singles contribution and
! and      (FX tY)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYBtX(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYBtX(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (5) FXtYRf     contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    end if

!   add -FXtYRf to transition moment
    qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! 6. form (1.c.2): F^Y * t^X(\omega^{'} - i\gamma}) * R_f
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)
    call get_tYFX(MNRf_tensor,Y_operator%prop_bar_ov,tX%t1,tX%t2,&
                  l1_tensor,l2_tensor,interm_t%nocc,interm_t%nvir)
! contract (FX tY)^m_e with Rf^e_m for singles contribution and
! and      (FX tY)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYBtX(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=tbarYBtX(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (6) FYtXRf     contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    end if

!   add -FYtXRf to transition moment
    qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! 7. form (1.d): G * t^X(\omega^{'} - i\gamma}) * t^Y(-\omega + i\gamma) * R_f
!
    ierr=talsh_tensor_init(MNRf_tensor%t1)
    ierr=talsh_tensor_init(MNRf_tensor%t2)

    call get_Gtt(MNRf_tensor,int_t,interm_t%t1,tX%t1,tX%t2,tY%t1,tY%t2,&
                l1_tensor,l2_tensor,interm_t%nocc, interm_t%nvir)
! contract (G * tX * tY)^m_e with Rf^e_m for singles contribution and
! and      (G * tX * tY)^m_e with Rf^{ef}_{mn} for doubles contribution to 2PTM 
    ierr=talsh_tensor_init(r_tensor)
    ierr=talsh_tensor_contract("T0fXY()+=GtXtY(m,e)*Rf(e,m)",r_tensor,&
                               MNRf_tensor%t1,Rs_f)
    if (debug_local) tran1=r_tens(1)
    ierr=talsh_tensor_contract("T0fXY()+=GtXtY(m,n,e,f)*Rf(e,f,m,n)",r_tensor,&
                               MNRf_tensor%t2,Rd_f,scale=ONE_QUARTER_C)
    if (debug_local) then
        tran=r_tens(1)
        tran2=tran-tran1
        write(*,'(2x,A,2E16.8,2x,A,2E16.8,2x,A,2E16.8)') &
              'DEBUG (7) GtXtYRf    contrib. to 2P-Tr.Mom.:',tran,'1:',tran1,'2:',tran2
    end if

!   add -GtXtYRf to transition moment
    qr_tran_left_XY = qr_tran_left_XY - r_tens(1)
!
! form term (2)
!
    tmom_b_XY = 0
 
    if (debug_local) then
        write(*,'(2x,A)') &
        'DEBUG (8-10) calculating M(w_f)[AXtY + AYtX + BtXtY] contrib. to 2P-Tr.Mom.'
    endif

    call get_cc_tran_right(exa_input, int_t,interm_t, X_operator, Y_operator, &
                             l1_tensor,l2_tensor, tX, tY, barM_f%t1,barM_f%t2, tmom_b_XY)

    if (debug_local) then
        write(*,'(2x,A,2E16.8)') &
        'DEBUG (8-10) total -M(w_f)[AXtY + AYtX + BtXtY] contrib. to 2P-Tr.Mom.',tmom_b_XY
        write(*,'(2x,A)') 'DEBUG exiting get_cc_tran_left'
    endif

!   add -M(w_f)[AXtY + AYtX + BtXtY] to transition moment
    qr_tran_left_XY = qr_tran_left_XY - tmom_b_XY
!
    ierr=talsh_tensor_destruct(MNRf_tensor%t1)
    ierr=talsh_tensor_destruct(MNRf_tensor%t2)
    ierr=talsh_tensor_destruct(r_tensor)    

end subroutine get_cc_tran_left

!------------------------------------------------------------------------------------------------------------------------------------------

subroutine get_eom_tran_left(exa_input, interm_t,Interm_t_tensor_X,Interm_t_tensor_Y,& 
                    Tprop_bar_tensor_X, Tprop_bar_tensor_min_Y,l1_tensor, l2_tensor, &
                    eom_results_t, eom_tran_left_XY,f_state)


    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    type(talsh_eom_results_tens_t), intent(inout):: eom_results_t
    type(talsh_comm_tens_t), intent(inout)       ::  Tprop_bar_tensor_X, Tprop_bar_tensor_min_Y !QR
    type(talsh_rsp_interm_tens_t), intent(inout) :: Interm_t_tensor_X,Interm_t_tensor_Y
    complex(kind=8), intent(inout)             :: eom_tran_left_XY
    integer,intent(in)                    :: f_state  !i excited state


    type(talsh_tens_t)         :: aux_1, aux_2, aux_3, aux_4
    type(talsh_tens_t)         :: Aux_ov, Aux_oovv

    integer(C_INT)     :: one_dims(1)
    integer(INTD), dimension(2)   :: ov_dims
    integer(INTD), dimension(4)   :: oovv_dims
    integer(INTD)  ::  ov_rank, oovv_rank

    integer(C_INT)          :: ierr
    type(talsh_tens_t)      :: result_tensor
    complex(8), pointer     :: result_tens(:)
    type(C_PTR)             :: body_p



! XY: EOM-CCSD left transition dipole moment 
    one_dims(1) = 1

    ierr=talsh_tensor_construct(aux_1,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(aux_2,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(aux_3,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(aux_4,C8,one_dims(1:0),init_val=ZERO)


    ov_rank=talsh_tensor_rank(l1_tensor)
    oovv_rank=talsh_tensor_rank(l2_tensor)
    ierr = talsh_tensor_dimensions(l1_tensor,ov_rank,ov_dims)
    ierr = talsh_tensor_dimensions(l2_tensor,oovv_rank,oovv_dims)
    ierr = talsh_tensor_construct(Aux_ov,C8,ov_dims,init_val=ZERO)
    ierr = talsh_tensor_construct(Aux_oovv,C8,oovv_dims,init_val=ZERO)


    ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,result_tens,one_dims)


    

!------------------------------------------------------
! aux_1 =\lambda * \epsilonX                          |
!------------------------------------------------------

!Single:
    ierr=talsh_tensor_contract("Aux1()+=Lamda(i,a)*EpiX(a,i)",& 
        aux_1,l1_tensor,Interm_t_tensor_X%epsilon_prop_vo,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux1()+=Lamda(i,j,a,b)*EpiX(a,b,i,j)",& 
        aux_1,l2_tensor,Interm_t_tensor_X%epsilon_prop_vvoo,scale=ONE_QUARTER_C)


!------------------------------------------------------
! aux_2 =\lambda * \epsilonY                          |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux2()+=Lamda(i,a)*EpiY(a,i)",& 
        aux_2,l1_tensor,Interm_t_tensor_Y%epsilon_prop_vo,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux2()+=Lamda(i,j,a,b)*EpiY(a,b,i,j)",& 
        aux_2,l2_tensor,Interm_t_tensor_Y%epsilon_prop_vvoo,scale=ONE_QUARTER_C)



!------------------------------------------------------
! aux_3 =\lambda * \right_state                       |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux3()+=Lamda(i,a)*RStat(a,i)",& 
        aux_3,l1_tensor,eom_results_t%value_vectors_t%rhs1(f_state),scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux3()+=Lamda(i,j,a,b)*RSta(a,b,i,j)",& 
        aux_3,l2_tensor,eom_results_t%value_vectors_t%rhs2(f_state),scale=ONE_QUARTER_C)



!------------------------------------------------------
! aux4 =+ ^{EOM}t_barY*\epsilonX                      |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux4()+=tbarY(i,a)*EpiX(a,i)",& 
        aux_4,Tprop_bar_tensor_min_Y%t1,Interm_t_tensor_X%epsilon_prop_vo,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux2()+=tbarY(i,j,a,b)*EpiX(a,b,i,j)",& 
        aux_4,Tprop_bar_tensor_min_Y%t2,Interm_t_tensor_X%epsilon_prop_vvoo,scale=ONE_QUARTER_C)



!------------------------------------------------------
! aux4 =+ ^{EOM}t_barX*\epsilonY                      |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux4()+=tbarX(i,a)*EpiY(a,i)",& 
        aux_4,Tprop_bar_tensor_X%t1,Interm_t_tensor_X%epsilon_prop_vo,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux4()+=tbarX(i,j,a,b)*EpiX(a,b,i,j)",& 
        aux_4,Tprop_bar_tensor_X%t2,Interm_t_tensor_X%epsilon_prop_vvoo,scale=ONE_QUARTER_C)




!------------------------------------------------------
! Aux_ov =+ eta_pri_ov=^{EOM}t_barX * ^{EOM}A^{Y}     |
!------------------------------------------------------
    call get_eta_prop_ov_prime(Aux_ov,interm_t%t2, interm_t%t1,&
                Tprop_bar_tensor_X%t2,Tprop_bar_tensor_X%t1,&
                Interm_t_tensor_Y,ov_dims)



!------------------------------------------------------
! Aux_ov =+ eta_pri_ov=^{EOM}t_barY * ^{EOM}A^{X}     |
!------------------------------------------------------
    call get_eta_prop_ov_prime(Aux_ov,interm_t%t2, interm_t%t1,&
                Tprop_bar_tensor_min_Y%t2,Tprop_bar_tensor_min_Y%t1,&
                Interm_t_tensor_X,ov_dims)

    

!------------------------------------------------------
! Aux_ov =- aux_1*^{EOM}t_barY_ov                     |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(i,a)+=tbarY(i,a)*aux1()",& 
        Aux_ov,Tprop_bar_tensor_min_Y%t1,aux_1,scale=MINUS_ONE)



!------------------------------------------------------
! Aux_ov =- aux_2*^{EOM}t_barX_ov                     |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(i,a)+=tbarX(i,a)*aux2()",& 
        Aux_ov,Tprop_bar_tensor_X%t1,aux_2,scale=MINUS_ONE)



!------------------------------------------------------
! Aux_oovv =+ eta_pri_oovv=^{EOM}t_barX * ^{EOM}A^{Y} |
!------------------------------------------------------
    call get_eta_prop_oovv_prime(Aux_oovv, interm_t%t1,&
                Tprop_bar_tensor_X%t2,Tprop_bar_tensor_X%t1,&
                Interm_t_tensor_Y,ov_dims)


!------------------------------------------------------
! Aux_oovv =+ eta_pri_oovv=^{EOM}t_barY * ^{EOM}A^{X} |
!------------------------------------------------------
    call get_eta_prop_oovv_prime(Aux_oovv,interm_t%t1,&
                Tprop_bar_tensor_min_Y%t2,Tprop_bar_tensor_min_Y%t1,&
                Interm_t_tensor_X,ov_dims)



!------------------------------------------------------
! Aux_oovv =- aux_1*^{EOM}t_barY_oovv                 |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(i,j,a,b)+=tbarY(i,j,a,b)*aux1()",& 
        Aux_oovv,Tprop_bar_tensor_min_Y%t2,aux_1,scale=MINUS_ONE)



!------------------------------------------------------
! Aux_oovv =- aux_2*^{EOM}t_barX_oovv                 |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(i,j,a,b)+=tbarX(i,j,a,b)*aux2()",& 
        Aux_oovv,Tprop_bar_tensor_X%t2,aux_2,scale=MINUS_ONE)


!------------------------------------------------------
! Left_TDM =- Aux_ov*right_state_vo                   |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=Aux(i,a)*Rstat(a,i)",& 
        result_tensor,Aux_ov,eom_results_t%value_vectors_t%rhs1(f_state),scale=MINUS_ONE)



!------------------------------------------------------
! Left_TDM =- Aux_oovv*right_state_vvoo               |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=Aux(i,j,a,b)*RSta(a,b,i,j)",& 
        result_tensor,Aux_oovv,eom_results_t%value_vectors_t%rhs2(f_state),scale=MINUS_ONE_QUARTER_C)



!------------------------------------------------------
! Left_TDM =+ aux3*aux4                               |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=aux3(i)*aux4(i)",& 
        result_tensor,aux_3,aux_4,scale=ONE)


    eom_tran_left_XY = result_tens(1)


! Destroy aux tensors
    ierr=talsh_tensor_destruct(aux_1)
    ierr=talsh_tensor_destruct(aux_2)
    ierr=talsh_tensor_destruct(aux_3)   
    ierr=talsh_tensor_destruct(aux_4)
    ierr=talsh_tensor_destruct(Aux_ov)
    ierr=talsh_tensor_destruct(Aux_oovv)
    ierr=talsh_tensor_destruct(result_tensor)
    ! deallocate (result_tens)
end subroutine get_eom_tran_left



!-------------------------------------------------------------------------------------------------------------------
subroutine get_eom_tran_right(exa_input, interm_t,&
                Interm_t_tensor_X,Interm_t_tensor_Y,& 
                Tprop_tensor_Y,Tprop_tensor_min_w_X,&
                l1_tensor, l2_tensor, &
                eom_results_t, eom_tran_right_XY, f_state)



    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    type(talsh_eom_results_tens_t), intent(inout) :: eom_results_t
    type(talsh_comm_tens_t), intent(inout)  ::  Tprop_tensor_Y, Tprop_tensor_min_w_X !QR
    type(talsh_rsp_interm_tens_t), intent(inout)   :: Interm_t_tensor_X,Interm_t_tensor_Y
    complex(kind=8), intent(inout)               :: eom_tran_right_XY
    integer, intent(in)                    :: f_state  !i excited state
            
            
    type(talsh_tens_t)         :: aux_1, aux_2, aux_3, aux_4
    type(talsh_tens_t)         :: Aux_vo, Aux_vvoo

    type(talsh_tens_t)         :: eta_ov_prim,  eta_oovv_prim
            
    integer(C_INT)     :: one_dims(1)
    integer(INTD), dimension(2)   :: vo_dims, ov_dims
    integer(INTD), dimension(4)   :: vvoo_dims, oovv_dims
    integer(INTD)  ::  vo_rank, vvoo_rank
            
    integer(C_INT)          :: ierr
    type(talsh_tens_t)      :: result_tensor
    complex(8), pointer     :: result_tens(:)
    type(C_PTR)             :: body_p
            
            
            
! XY: EOM-CCSD right transition dipole moment 
    one_dims(1) = 1
            
    ierr=talsh_tensor_construct(aux_1,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(aux_2,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(aux_3,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(aux_4,C8,one_dims(1:0),init_val=ZERO)

            
    vo_rank=talsh_tensor_rank(interm_t%t1)
    vvoo_rank=talsh_tensor_rank(interm_t%t2)
    ierr = talsh_tensor_dimensions(interm_t%t1,vo_rank,vo_dims)
    ierr = talsh_tensor_dimensions(interm_t%t2,vvoo_rank,vvoo_dims)

    ov_dims(1) = vo_dims(2)
    ov_dims(2) = vo_dims(1)
    oovv_dims(1:2) = vvoo_dims(3:4)
    oovv_dims(3:4) = vvoo_dims(1:2)

    ierr = talsh_tensor_construct(Aux_vo,C8,vo_dims,init_val=ZERO)
    ierr = talsh_tensor_construct(Aux_vvoo,C8,vvoo_dims,init_val=ZERO)


    ierr = talsh_tensor_construct(eta_ov_prim,C8,ov_dims,init_val=ZERO)
    ierr = talsh_tensor_construct(eta_oovv_prim,C8,oovv_dims,init_val=ZERO)
            
    ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,result_tens,one_dims)


!------------------------------------------------------
! aux_1 =\lambda * \epsilonX                          |
!------------------------------------------------------

!Single:
    ierr=talsh_tensor_contract("Aux1()+=Lamda(i,a)*EpiX(a,i)",& 
        aux_1,l1_tensor,Interm_t_tensor_X%epsilon_prop_vo,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux1()+=Lamda(i,j,a,b)*EpiX(a,b,i,j)",& 
        aux_1,l2_tensor,Interm_t_tensor_X%epsilon_prop_vvoo,scale=ONE_QUARTER_C)


!------------------------------------------------------
! aux_2 =\lambda * \epsilonY                          |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux2()+=Lamda(i,a)*EpiY(a,i)",& 
        aux_2,l1_tensor,Interm_t_tensor_Y%epsilon_prop_vo,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux2()+=Lamda(i,j,a,b)*EpiY(a,b,i,j)",& 
        aux_2,l2_tensor,Interm_t_tensor_Y%epsilon_prop_vvoo,scale=ONE_QUARTER_C)
    


!------------------------------------------------------
! aux_3 =\lambda * \tY                                |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux3()+=Lamda(i,a)*TY(a,i)",& 
        aux_3,l1_tensor,Tprop_tensor_Y%t1,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux3()+=Lamda(i,j,a,b)*TY(a,b,i,j)",& 
        aux_3,l2_tensor,Tprop_tensor_Y%t2,scale=ONE_QUARTER_C)



!------------------------------------------------------
! aux_4 =\lambda * \tX                                |
!------------------------------------------------------
!Single:
    ierr=talsh_tensor_contract("Aux4()+=Lamda(i,a)*TX(a,i)",& 
        aux_4,l1_tensor,Tprop_tensor_min_w_X%t1,scale=ONE)


!Double:
    ierr=talsh_tensor_contract("Aux4()+=Lamda(i,j,a,b)*TX(a,b,i,j)",& 
        aux_4,l2_tensor,Tprop_tensor_min_w_X%t2,scale=ONE_QUARTER_C)



!------------------------------------------------------
! Aux_vo =- aux1 * \tY_vo                             |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,i)+=aux1()*Ty(a,i)",& 
        Aux_vo,aux_1,Tprop_tensor_Y%t1,scale=MINUS_ONE)


!------------------------------------------------------
! Aux_vo =- aux2 * \tX_vo                             |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,i)+=aux2()*Tx(a,i)",& 
        Aux_vo,aux_2,Tprop_tensor_min_w_X%t1,scale=MINUS_ONE)


!------------------------------------------------------
! Aux_vo =- aux3 * \epsilonX_vo                       |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,i)+=aux3()*epsilonX(a,i)",& 
        Aux_vo,aux_3,Interm_t_tensor_X%epsilon_prop_vo,scale=MINUS_ONE)


!------------------------------------------------------
! Aux_vo =- aux4 * \epsilonY_vo                       |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,i)+=aux4()*epsilonY(a,i)",& 
        Aux_vo,aux_4,Interm_t_tensor_Y%epsilon_prop_vo,scale=MINUS_ONE)



!------------------------------------------------------
! Aux_vvoo =- aux1 * \tY_vvoo                         |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,b,i,j)+=aux1()*Ty(a,b,i,j)",& 
        Aux_vvoo,aux_1,Tprop_tensor_Y%t2,scale=MINUS_ONE)


!------------------------------------------------------
! Aux_vvoo =- aux2 * \tX_vvoo                         |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,b,i,j)+=aux2()*Tx(a,b,i,j)",& 
        Aux_vvoo,aux_2,Tprop_tensor_min_w_X%t2,scale=MINUS_ONE)
        
!------------------------------------------------------
! Aux_vvoo =- aux3 * \epsilonX_vvoo                   |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,b,i,j)+=aux3()*epsilonX(a,b,i,j)",& 
        Aux_vvoo,aux_3,Interm_t_tensor_X%epsilon_prop_vvoo,scale=MINUS_ONE)



!------------------------------------------------------
! Aux_vvoo =- aux4 * \epsilonY_vvoo                   |
!------------------------------------------------------
    ierr=talsh_tensor_contract("Aux(a,b,i,j)+=aux4()*epsilonY(a,b,i,j)",& 
        Aux_vvoo,aux_4,Interm_t_tensor_Y%epsilon_prop_vvoo,scale=MINUS_ONE)



!------------------------------------------------------
! eta_ov_prim_X_ov = leftstate * A^{X}               |
!------------------------------------------------------
    ierr=talsh_tensor_init(eta_ov_prim)
    call get_eta_prop_ov_prime(eta_ov_prim,interm_t%t2, interm_t%t1,&
                eom_results_t%value_vectors_t%lhs2(f_state),&
                eom_results_t%value_vectors_t%lhs1(f_state),&
                Interm_t_tensor_X,ov_dims)



!------------------------------------------------------
! eta_ov_prim_X_oovv = leftstate * A^{X}             |
!------------------------------------------------------
    ierr=talsh_tensor_init(eta_oovv_prim)
    call get_eta_prop_oovv_prime(eta_oovv_prim,interm_t%t1,&
                eom_results_t%value_vectors_t%lhs2(f_state),&
                eom_results_t%value_vectors_t%lhs1(f_state),&
                Interm_t_tensor_X,ov_dims)


    ! call print_tensor(eta_ov_prim, 1.0d-10, "eta_ov_prim")

!------------------------------------------------------
! Right_TDM =- eta_ov_prim_X_ov *tY_vo                |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=Eta_pri_X(i,a)*Ty(a,i)",& 
        result_tensor,eta_ov_prim,Tprop_tensor_Y%t1,scale=MINUS_ONE)



!------------------------------------------------------
! Right_TDM =- eta_ov_prim_X_oovv *tY_vvoo             |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=Eta_pri_X(i,j,a,b)*Ty(a,b,i,j)",& 
        result_tensor,eta_oovv_prim,Tprop_tensor_Y%t2,scale=MINUS_ONE_QUARTER_C)



!------------------------------------------------------
! eta_ov_prim_Y_ov = leftstate * A^{Y}               |
!------------------------------------------------------
    ierr=talsh_tensor_init(eta_ov_prim)
    call get_eta_prop_ov_prime(eta_ov_prim,interm_t%t2, interm_t%t1,&
            eom_results_t%value_vectors_t%lhs2(f_state),&
            eom_results_t%value_vectors_t%lhs1(f_state),&
            Interm_t_tensor_Y,ov_dims)



!------------------------------------------------------
! eta_ov_prim_Y_oovv = leftstate * A^{Y}             |
!------------------------------------------------------
    ierr=talsh_tensor_init(eta_oovv_prim)
    call get_eta_prop_oovv_prime(eta_oovv_prim,interm_t%t1,&
            eom_results_t%value_vectors_t%lhs2(f_state),&
            eom_results_t%value_vectors_t%lhs1(f_state),&
            Interm_t_tensor_Y,ov_dims)


!------------------------------------------------------
! Right_TDM =- eta_ov_prim_Y_ov *tX_vo                |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=Eta_pri_Y(i,a)*Tx(a,i)",& 
            result_tensor,eta_ov_prim,Tprop_tensor_min_w_X%t1,scale=MINUS_ONE)


!------------------------------------------------------
! Right_TDM =- eta_ov_prim_Y_oovv *tX_vvoo             |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=Eta_pri_Y(i,j,a,b)*Tx(a,b,i,j)",& 
            result_tensor,eta_oovv_prim,Tprop_tensor_min_w_X%t2,scale=MINUS_ONE_QUARTER_C)


!------------------------------------------------------
! Right_TDM =- leftstate_ov * Aux_vo                  |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=LState(i,a)*Aux(a,i)",& 
            result_tensor,eom_results_t%value_vectors_t%lhs1(f_state),Aux_vo,scale=MINUS_ONE)

        
!------------------------------------------------------
! Right_TDM =- leftstate_oovv * Aux_vvoo              |
!------------------------------------------------------
    ierr=talsh_tensor_contract("TDML()+=LState(i,j,a,b)*Aux(a,b,i,j)",& 
            result_tensor,eom_results_t%value_vectors_t%lhs2(f_state),Aux_vvoo,scale=MINUS_ONE_QUARTER_C)



    eom_tran_right_XY = result_tens(1)

    ! Destroy aux tensors
    ierr=talsh_tensor_destruct(aux_1)
    ierr=talsh_tensor_destruct(aux_2)
    ierr=talsh_tensor_destruct(aux_3)   
    ierr=talsh_tensor_destruct(aux_4)
    ierr=talsh_tensor_destruct(Aux_vo)
    ierr=talsh_tensor_destruct(Aux_vvoo)
    ierr=talsh_tensor_destruct(eta_ov_prim)
    ierr=talsh_tensor_destruct(eta_oovv_prim)
    ierr=talsh_tensor_destruct(result_tensor)
    ! deallocate (result_tens)

end subroutine get_eom_tran_right


!-------------------------------------------------------------------------------------------------------------------

subroutine get_scattering_amplitudes(qr_tran_left_AB, qr_tran_left_CD, &
                qr_tran_right_AB,qr_tran_right_CD,qr_tpa)

    complex(kind=8), intent(in)    :: qr_tran_left_AB, qr_tran_left_CD
    complex(kind=8), intent(in)    :: qr_tran_right_AB, qr_tran_right_CD
    complex(kind=8), intent(out)   :: qr_tpa

    qr_tpa = (qr_tran_left_AB*qr_tran_right_CD)*0.5 
    qr_tpa = qr_tpa + conjg(qr_tran_left_CD*qr_tran_right_AB)*0.5

end subroutine



!--------------------------------------------------------------------------------------------------------------------
subroutine get_eta_prop_ov_prime(eta_prop_ov_tensor, t2_tensor,t1_tensor,&
                                l2_tensor,l1_tensor,Interm_t_tensor,ov_dims)

    type(talsh_tens_t),intent(inout) :: eta_prop_ov_tensor
    ! type(talsh_orb_prop_integ_tens_t) :: orb_int_tensor
    type(talsh_rsp_interm_tens_t),intent(inout) :: Interm_t_tensor
    type(talsh_tens_t),intent(inout) :: t2_tensor,t1_tensor
    type(talsh_tens_t),intent(inout) :: l2_tensor,l1_tensor
    integer(INTD), dimension(2), intent(in)   :: ov_dims


    type(talsh_tens_t)            :: LR_Aux_oo,LR_Aux_vv,LR_Aux_vo
    integer(INTD), dimension(2)   :: vv_dims, oo_dims
    integer :: ierr

    type(talsh_tens_t)            :: Transpose_vo
    integer(INTD), dimension(2)   :: vo_dims
    integer(C_INT)     :: one_dims(1)
    type(talsh_tens_t) :: uncon_tensor
    logical :: debug_local = .false.


    one_dims(1) = 1
    vo_dims(1) = ov_dims(2)
    vo_dims(2) = ov_dims(1)
    oo_dims(1:2) = ov_dims(1)
    vv_dims(1:2) = ov_dims(2)


    ierr=talsh_tensor_construct(LR_Aux_oo,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(LR_Aux_vv,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(LR_Aux_vo,C8,vo_dims,init_val=ZERO)


    ierr=talsh_tensor_construct(Transpose_vo,C8,vo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
    ierr=talsh_tensor_contract('A(a,i)+=Y+(i,a)*K()',Transpose_vo,Interm_t_tensor%prop_bar_ov,one_tensor)
    ierr=talsh_tensor_construct(uncon_tensor,C8,one_dims(1:0),init_val=ZERO)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_tran ', ierr



!XY  \Eta_1 ia=+\lambdad ijab \epison ej
    ierr=talsh_tensor_contract('E(i,a)+=L(i,j,a,e)*Epis(e,j)',eta_prop_ov_tensor,l2_tensor,&
                            Interm_t_tensor%epsilon_prop_vo,scale=ONE)                    
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o1 ', ierr

!LH     \Eta_1 ia=+\Ys ea \lambdas ie
    ierr=talsh_tensor_contract('E(i,a)+=Y(e,a)*L(i,e)',eta_prop_ov_tensor,Interm_t_tensor%prop_bar_vv,l1_tensor,scale=ONE)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o3 ', ierr

!LH     \Eta_1 ia=-\Ys im \lambdas ma
    ierr=talsh_tensor_contract('E(i,a)+=Y(i,m)*L(m,a)',eta_prop_ov_tensor,Interm_t_tensor%prop_bar_oo,l1_tensor,scale=MINUS_ONE)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o4 ', ierr

!LH     \Eta_1=-\Ys ma \ts em \lambdas ie
    ierr=talsh_tensor_init(LR_Aux_oo)       
    ierr=talsh_tensor_contract('H(i,m)+=T(e,m)*L(i,e)',LR_Aux_oo,t1_tensor,l1_tensor,scale=ONE)
    ierr=talsh_tensor_contract('E(i,a)+=H(i,m)*Y(m,a)',eta_prop_ov_tensor,LR_Aux_oo,Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o5 ', ierr

!LH     \Eta_1=-\Ys ie \ts em \llambdas ma
    ierr=talsh_tensor_init(LR_Aux_oo)
    ierr=talsh_tensor_contract('H(i,m)+=Y(i,e)*T(e,m)',Lr_Aux_oo,Interm_t_tensor%prop_bar_ov,t1_tensor,scale=ONE)
    ierr=talsh_tensor_contract('E(i,a)+=H(i,m)*L(m,a)',eta_prop_ov_tensor,LR_Aux_oo,l1_tensor,scale=MINUS_ONE)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o6 ', ierr


!LH     \Eta_1= -one_half \td femn \lambdad mife \Ys na 
    ierr=talsh_tensor_init(LR_Aux_oo)
    ierr=talsh_tensor_contract('H(i,n)+=T(f,e,m,n)*L(m,i,f,e)',LR_Aux_oo,t2_tensor,l2_tensor)
    ierr=talsh_Tensor_contract('E(i,a)+=H(i,n)*Y(n,a)',& 
          eta_prop_ov_tensor,LR_Aux_oo,Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE_HALF)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o28 ', ierr

!LH     \Eta_1=-one_half \td  fenm \lambdad nmfa \Ys ie 
    ierr=talsh_tensor_init(LR_Aux_vv)
    ierr=talsh_tensor_contract('H(e,a)+=T(f,e,n,m)*L(n,m,f,a)',LR_Aux_vv,t2_tensor,l2_tensor)
    ierr=talsh_tensor_contract('E(i,a)+=H(e,a)*Y(i,e)',&
          eta_prop_ov_tensor,LR_Aux_vv,Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE_HALF)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_o29 ', ierr



    ierr=talsh_tensor_destruct(LR_Aux_vv)
    ierr=talsh_tensor_destruct(LR_Aux_vo)
    ierr=talsh_tensor_destruct(LR_Aux_oo)


end subroutine


!-------------------------------------------------------------------------------------------------------------------------------
subroutine get_eta_prop_oovv_prime( eta_prop_oovv_tensor,t1_tensor,l2_tensor,l1_tensor,Interm_t_tensor,ov_dims)

    type(talsh_tens_t), intent(inout) :: eta_prop_oovv_tensor
    type(talsh_rsp_interm_tens_t), intent(inout) :: Interm_t_tensor
    type(talsh_tens_t), intent(inout) :: t1_tensor
    type(talsh_tens_t), intent(inout) :: l2_tensor,l1_tensor
    integer(INTD), dimension(2)   :: ov_dims


    type(talsh_tens_t)            :: LR_Aux_oo,LR_Aux_vv
    integer(INTD), dimension(4)   :: oovv_dims
    integer(INTD), dimension(2)   :: oo_dims, vv_dims
    integer :: ierr
    logical :: debug_local = .false.
    integer(C_INT)     :: one_dims(1)
    type(talsh_tens_t) :: uncon_tensor


    one_dims(1) = 1
    oo_dims(1:2) = ov_dims(1)
    vv_dims(1:2) = ov_dims(2)
    ierr=talsh_tensor_construct(uncon_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_construct(LR_Aux_oo,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(LR_Aux_vv,C8,vv_dims,init_val=ZERO)


!LH     \Eta_2=+ \Pmm ij \Pmm ab \lambdas ia \Ys jb
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(i,a)*Y(j,b)',eta_prop_oovv_tensor,l1_tensor,Interm_t_tensor%prop_bar_ov,scale=ONE)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(j,a)*Y(i,b)',eta_prop_oovv_tensor,l1_tensor,&
                                                        Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(i,b)*Y(j,a)',eta_prop_oovv_tensor,l1_tensor,&
                                                        Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(j,b)*Y(i,a)',eta_prop_oovv_tensor,l1_tensor,Interm_t_tensor%prop_bar_ov,scale=ONE)
    if(debug_local) print *,'LH_talsh eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=+\Pmm ab \lambdad ijae \Ys eb 
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(i,j,a,e)*Y(e,b)', &
             eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_vv,scale=ONE)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(i,j,b,e)*Y(e,a)', & 
             eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_vv,scale=MINUS_ONE)
    if(debug_local) print *,'LH_talsh eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2= - \Pmm ij \lambdad imab \Ys jm 
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(i,m,a,b)*Y(j,m)', &
             eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_oo,scale=MINUS_ONE)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=L(j,m,a,b)*Y(i,m)', &
             eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_oo,scale=ONE)
    if(debug_local) print *,'LH_talsh eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=-\Pmm ab \ts em \lambdad ijae \Ys mb 
    ierr=talsh_tensor_init(LR_Aux_vv)
    ierr=talsh_tensor_contract('H(e,b)+=T(e,m)*Y(m,b)',LR_Aux_vv,t1_tensor,Interm_t_tensor%prop_bar_ov)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=H(e,b)*L(i,j,a,e)',eta_prop_oovv_tensor,LR_Aux_vv,l2_tensor,scale=MINUS_ONE)       
    ierr=talsh_tensor_contract('E(i,j,a,b)+=H(e,a)*L(i,j,b,e)',eta_prop_oovv_tensor,LR_Aux_vv,l2_tensor,scale=ONE)         
    if(debug_local) print *,'LH_talsh eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=-\Pmm ij \ts em \lambdad imab\Ys je 
    ierr=talsh_tensor_init(LR_Aux_oo)
    ierr=talsh_tensor_contract('H(j,m)+=T(e,m)*Y(j,e)',LR_Aux_oo,t1_tensor,Interm_t_tensor%prop_bar_ov)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=H(j,m)*L(i,m,a,b)',eta_prop_oovv_tensor,LR_Aux_oo,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract('E(i,j,a,b)+=H(i,m)*L(j,m,a,b)',eta_prop_oovv_tensor,LR_Aux_oo,l2_tensor,scale=ONE)
    if(debug_local) print *,'LH_talsh eta_prop_oovv_tensor ierr ', ierr


    ierr=talsh_tensor_destruct(LR_Aux_vv)
    ierr=talsh_tensor_destruct(LR_Aux_oo)


end subroutine


!------------------------------------------------------------------------------------------------------------------------------

subroutine check_biorthognal(eom_results_t,f_state, biorth_sign_i)
    type(talsh_eom_results_tens_t), intent(inout)        :: eom_results_t
    integer, intent(in)                    :: f_state  ! target (final) excited state in TPA
    ! aspg, 16/02/2024
    ! the TPA/EOM code needs more work  to ensure  the biorthonormality of L and R eigenvectors.
    ! we have observed that even when L*R is (near) unity, in absolute value, it may be that one of the 
    ! two eigenvectors has a -1 phase.  when that happens, the TPA cross section will have the correct 
    ! magnitude but will have a spurious sign change.
    ! variable  biorth_sign_i will track that and flip sign in case L*R is not positive for a TPA target state.. 
    complex(kind=8), intent(inout)         :: biorth_sign_i 


    integer(C_INT)     :: one_dims(1)

    integer(C_INT)          :: ierr
    type(talsh_tens_t)      :: result_tensor
    complex(8), pointer     :: result_tens(:)
    complex(8)              :: LRnorm
    type(C_PTR)             :: body_p
    integer                 :: i, j, N
    real(kind=8)            :: threshold_small = 1.0d-8
    real(kind=8)            :: LRnorm_re_residual
    logical                 :: LRnorm_re_is_one, LRnorm_re_is_zero, LRnorm_im_is_zero

    one_dims(1) = 1
    ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,result_tens,one_dims)

     N = size(eom_results_t%value_vectors_t%lhs1)
    print *, "Verifying biorthogonality of CC eigenvectors in TPA calculation"
    print *, "  tolerance on deviation from strict biorthogonality (Re and Im):",threshold_small
    print *, "  Number of states:", N 
    print *, "  Target state: ",f_state

    do i=1, N ! Left eigenstates 
        do j=1, N ! Right eigentstates 

           ierr=talsh_tensor_init(result_tensor, ZERO)

           ierr=talsh_tensor_contract("Bio()+=LState(i,a)*RState(a,i)",& 
               result_tensor,eom_results_t%value_vectors_t%lhs1(i),&
               eom_results_t%value_vectors_t%rhs1(j),scale=ONE)

           ierr=talsh_tensor_contract("Bio()+=LState(i,j,a,b)*RState(a,b,i,j)",& 
               result_tensor,eom_results_t%value_vectors_t%lhs2(i),&
               eom_results_t%value_vectors_t%rhs2(j),scale=ONE_QUARTER_C)

           LRnorm = result_tens(1) 

           LRnorm_re_is_zero = .false.
           LRnorm_im_is_zero = .false.
           if (abs(real(LRnorm)).lt.threshold_small) LRnorm_re_is_zero = .true.
           if (abs(aimag(LRnorm)).lt.threshold_small) LRnorm_im_is_zero = .true.
   
           ! L(i)*R(i) = ONE  and L(i)*R(j) = ZERO  when  L and R should be biorthonormal. 
           ! when that is not the case we issue warnings
           if ((i.eq.j)) then
! for i equal to j, and deviation from an acceptably small value for the imaginary part
! is should be flagged. for the real part, we check not only the magnitude but also the sign
! of the real part, that should only differ from 1 by an acceptably small amount. 

              LRnorm_re_is_one = .false.
              if (.not.LRnorm_re_is_zero) then
                 LRnorm_re_residual = (real(LRnorm) - real(ONE))
                 if (abs(LRnorm_re_residual).lt.threshold_small) LRnorm_re_is_one = .true.
              end if
              
              if ((.not.LRnorm_im_is_zero).or.(.not.LRnorm_re_is_one)) then 
                 print *, 'Warning ! LHS state',i,'and RHS state',j,' L*R = ',LRnorm,' when it should be',ONE
                 if ((i.eq.f_state).and.(LRnorm_re_residual.le.(-1.0d0))) then
              ! for the target state, we see whether to flip the sign. 
                    biorth_sign_i = MINUS_ONE
                    print *, '   note: sign will be flipped for target state'
                 end if
              end if
           else 
! for i different than j, any deviation from an acceptably small value for the real or imaginary
! parts represents a problem, so we flag it
              if ((.not.LRnorm_re_is_zero).or.(.not.LRnorm_im_is_zero)) &
                 print *, 'Warning ! LHS state',i,'and RHS state',j,' L*R = ',LRnorm,' when it should be',ZERO
           end if
        end do
    end do

    ierr=talsh_tensor_destruct(result_tensor)

end subroutine 


end module
