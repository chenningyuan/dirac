module talsh_ccrsp
!This module contains the routines needed to calculate CC-CI(EOM-CC) response property


    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global
    use talsh_eom_rspeq_solver
    use talsh_eom_QR
    use talsh_eom_LR
    use talsh_ccrsp_common
    use talsh_cc_QR

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public talsh_CCRSP_driver

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine talsh_CCRSP_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,rr)
!This subroutine drives the response property calculation 
!it now includes linear and quadratic response
!Written by Xiang Yuan Janunary 2022, Lille, France

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    integer                                     :: rr ! rank of response 
    
!   Wave-function tensor needed to construct the response function
    logical          :: CC_CC
    

    if (exa_input%wf_type == 1) then 
        CC_CC = .true.
        call print_date('Activate CC-CC(LR-CC) response')
!       call print_date('Activate     CC (CC-CC) response')
    else 
        CC_CC = .false.
        call print_date('Activate CC-CI(EOM-CC) response')
!       call print_date('Activate EOM-CC (CC-CI) response')
    end if 

    if (rr == 1) then
        if (exa_input%ccrsp_symmetric) then
            call talsh_ccsd_lr_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,CC_CC)
        else
            call talsh_ccsd_lr_asym_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,CC_CC)
        endif
    else if (rr==2) then
        call talsh_ccsd_qr_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,CC_CC)
    else 
        stop 'Unknown option in response calculation. Only LR and QR are implemented.'
    end if


end subroutine talsh_CCRSP_driver


!------------------------------------------------------------------------------------------------------------------

subroutine talsh_ccsd_lr_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,CC_CC)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    logical                                     :: CC_CC
    

!   Wave-function tensor needed to construct the response function

    complex(kind=8)               :: omega1, omega2
    complex(kind=8),allocatable   :: omega(:), omega_min(:)
    complex(kind=8)               :: eom_lr
    complex(kind=8),allocatable   :: eom_lrs(:,:,:)
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR

    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:), Interm_t_tensors_B(:)
    integer                       :: N_prop, N_prop_A, N_prop_B
    integer                       :: N_oper_lr
    integer                       :: N_omega
    integer                       :: i, j, k, n, m
    integer, allocatable          :: Prop_A(:), Prop_B(:)

! 
    type(C_PTR) :: body_p
    complex(8), pointer    ::prop_tens(:,:,:,:)
    integer            :: tens_rank
    integer :: ierr
    logical          :: exchange
    

    N_prop = 3

    call print_date('Do linear response')

    Prop_A = exa_input%Prop_A
    Prop_B = exa_input%Prop_B

    N_oper_lr = 2
    N_omega = exa_input%N_omega_lr
    N_prop_A = size(Prop_A)
    N_prop_B = size(Prop_B)

    allocate(omega(N_oper_lr))
    allocate(omega_min(N_oper_lr))
    allocate(Interm_t_tensors_A(N_prop_A))
    allocate(Interm_t_tensors_B(N_prop_B))
    allocate(Tprop_tensors(N_oper_lr))
    allocate(Tprop_tensors_min_w(N_oper_lr))
    allocate(eom_lrs(N_omega,N_prop_A,N_prop_B))


    do  i=1, N_prop_A
        m = Prop_A(i)
        exchange = .False.
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_A(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC, exchange)
    end do
    
    
    do  i=1, N_prop_B
        m = Prop_B(i)
        exchange = .True.
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_B(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC, exchange)
    end do

    call print_date('Constructed intermediates needed in CC Linear Response')


    do n=1, N_omega
        omega1 =  dcmplx(exa_input%w_lr_re(n), exa_input%w_lr_im(n))
        omega2 =  -omega1
        omega(1) = omega1
        omega(2) = omega2
        print *, "Read frequency1:", omega(1)
        print *, "Read frequency2:", omega(2)

        omega_min(1) = dcmplx(-real(omega(1)),aimag(omega(1)))
        omega_min(2) = dcmplx(-real(omega(2)),aimag(omega(2)))
        print *, "Read frequency min_1:", omega_min(1)
        print *, "Read frequency min_2:", omega_min(2)


        do i=1, N_prop_A
            call print_lr_setup("R","A",i,omega(1))
            call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega(1), Tprop_tensors(1))
            if (omega(1).ne.ZERO) then
                call print_lr_setup("R","A",i,omega_min(1))
                call solve_eom_lr_right(exa_input, interm_t, int_t, &
                                    Interm_t_tensors_A(i), omega_min(1), Tprop_tensors_min_w(1))
            end if

            do j=1, N_prop_B
                call print_lr_setup("R","B",j,omega(2))
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), omega(2), Tprop_tensors(2))
                if (omega(2).ne.ZERO) then
                    call print_lr_setup("R","B",j,omega_min(2))
                    call solve_eom_lr_right(exa_input, interm_t, int_t, &
                                            Interm_t_tensors_B(j), omega_min(2), Tprop_tensors_min_w(2))
                end if

                call print_date('Solved All 1st response right equations')

                if (omega(1).ne.ZERO) then
                    call get_eom_lr(exa_input,interm_t, int_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,& 
                        Tprop_tensors, Tprop_tensors_min_w, i,j,& 
                        l1_tensor,l2_tensor,eom_lr,CC_CC)
                else
                    call print_date('Note: Calculating static linear response.')
                    call get_eom_lr(exa_input,interm_t, int_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,& 
                        Tprop_tensors, Tprop_tensors, i,j,& 
                        l1_tensor,l2_tensor,eom_lr,CC_CC)
                end if

                eom_lrs(n,i,j) = eom_lr
                write(*,"(A4,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), ">> = "
                write(*,"(F16.8,A,F16.8)"), real(eom_lrs(n,i,j)), ' \ ', aimag(eom_lrs(n,i,j))

                call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,2)
            end do

                call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,1)
            end do
        
        end do


        do i=1, N_prop_A
            call destroy_Interm_t_tensors(Interm_t_tensors_A(i),CC_CC)
        end do

        do i=1, N_prop_B
            call destroy_Interm_t_tensors(Interm_t_tensors_B(i),CC_CC)
        end do
        


        if (CC_CC) then 
            write(*,"(A54)"), "*************CC-CC(LR-CC) linear response*************"
        else
            write(*,"(A54)"), "*************CC-CI(EOM-CC) linear response************"
        end if
        
        do n=1, N_omega
            do i=1, N_prop_A
                do j=1, N_prop_B
                    write(*,"(A4,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), ">> = "
                    write(*,"(F16.8,A,F16.8)"), real(eom_lrs(n,i,j)), ' \ ', aimag(eom_lrs(n,i,j))
                end do
            end do
        end do


        DEALLOCATE(omega)
        DEALLOCATE(omega_min)
        DEALLOCATE(Tprop_tensors)
        DEALLOCATE(Tprop_tensors_min_w)
        DEALLOCATE(Interm_t_tensors_A)
        DEALLOCATE(Interm_t_tensors_B)
        DEALLOCATE(eom_lrs)

        write(*,"(A54)"), "**********************************************************"  

        if (CC_CC) then     
            call print_date("Finish CC-CC Linear Response Calculations")
        else
            call print_date("Finish CC-CI Linear Response Calculations")
        end if 


end subroutine talsh_ccsd_lr_driver

!------------------------------------------------------------------------------------------------------------------

subroutine talsh_ccsd_lr_asym_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,CC_CC)

    use talsh_gradient, only : talsh_dm_driver, talsh_dm_set_1rdm_dump_filename

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    logical, intent(in)                         :: CC_CC
    

!   Wave-function tensor needed to construct the response function

    complex(kind=8)               :: omega1, omega2, omega_B, minus_omega_B
    complex(kind=8),allocatable   :: omega(:), omega_min(:)
    complex(kind=8)               :: eom_lr
    complex(kind=8),allocatable   :: eom_lrs(:,:,:)
    ! RHS response vectors 
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:) 
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors_min_w(:) 
    ! LHS response vectors
    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:) 
    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors_min_w(:) 
    ! ??
    type(talsh_comm_tens_t), allocatable :: Bvec_tensors(:)
    type(talsh_comm_tens_t), allocatable :: Bvec_tensors_min_w(:)

    ! intermediates for response
    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:), Interm_t_tensors_B(:)


    integer                       :: N_prop, N_prop_A, N_prop_B
    integer                       :: N_oper_lr
    integer                       :: N_omega
    integer                       :: i, j, k, n, m

    integer                       :: A, B ! indices for operators in <<A;B>>, so A=1, B=2 but here 

    ! mapping to the list of selected property integrals, out of those available (e.g. in MDPROP)
    integer, allocatable          :: Prop_A(:), Prop_B(:)

! 
    type(C_PTR) :: body_p
    complex(8), pointer    ::prop_tens(:,:,:,:)
    integer            :: tens_rank
    integer :: ierr
    character(len=1)   :: char_a, char_b, char_omega
    
    N_prop = 3 ! why this value ?

    call print_date('Do linear response, asymmetric formulation')

    Prop_A = exa_input%Prop_A
    Prop_B = exa_input%Prop_B

    N_oper_lr  = 2 ! will only solve for T^B, so once debugged see if this can be set to 1
    
    A    = 1
    B    = 2 ! operator B will later have index 1 here, since we're only solving for it ?

    N_omega = exa_input%N_omega_lr

    N_prop_A = size(Prop_A)
    N_prop_B = size(Prop_B)

    allocate(Interm_t_tensors_A(N_prop_A))
    allocate(Interm_t_tensors_B(N_prop_B))

    allocate(eom_lrs(N_omega,N_prop_A,N_prop_B))

    allocate(omega(N_oper_lr))
    allocate(omega_min(N_oper_lr))

    allocate(Tprop_tensors(N_oper_lr))
    allocate(Tprop_tensors_min_w(N_oper_lr))
    allocate(Tprop_bar_tensors(N_oper_lr))
    allocate(Tprop_bar_tensors_min_w(N_oper_lr))

    allocate(Bvec_tensors(N_oper_lr))
    allocate(Bvec_tensors_min_w(N_oper_lr))

    do  i=1, N_prop_A
        m = Prop_A(i)
        ! calculates \epsilon^A, \eta^A
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_A(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC)
    end do
    
    do  i=1, N_prop_B
        m = Prop_B(i)
        ! calculates \epsilon^B, \eta^B
        call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_B(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC)
    end do

    call print_date('Constructed intermediates needed in asymmetric CC Linear Response')


    do n=1, N_omega
        omega1 =   dcmplx(exa_input%w_lr_re(n), exa_input%w_lr_im(n))
        omega2 =  -omega1
        omega(A) = omega1
        omega(B) = omega2
        print *, "Read frequency  w_A:", omega(A)
        print *, "Read frequency  w_B:", omega(B)

        omega_min(A) = dcmplx(-real(omega(A)),aimag(omega(A)))
        omega_min(B) = dcmplx(-real(omega(B)),aimag(omega(B)))
        print *, "Read frequency -w_A: ", omega_min(A),A
        print *, "Read frequency -w_B: ", omega_min(B),B

!
! note the asymmetric response only requires the solutiion of equations for perturbation B
!
        do j=1, N_prop_B
            if (CC_CC) then
!
!               solve for +w
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                     omega(B), Tprop_tensors(B))
                call get_left_bvec_cccc(exa_input,int_t,interm_t,Interm_t_tensors_B(j), &
                    l1_tensor,l2_tensor, Tprop_tensors(B), Bvec_tensors(B))
                call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                    omega(B), Bvec_tensors(B), Tprop_bar_tensors(B))
!                    
!               solve for -w
                if (omega(B).ne.ZERO) then 
                    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                         omega_min(B), Tprop_tensors_min_w(B))
                    call get_left_bvec_cccc(exa_input,int_t,interm_t,Interm_t_tensors_B(j), &
                         l1_tensor,l2_tensor, Tprop_tensors_min_w(B), Bvec_tensors_min_w(B))
                    call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                             omega_min(B),Bvec_tensors_min_w(B),Tprop_bar_tensors_min_w(B))
                end if 

            else ! using EOM approximation to LHS equation, consistent with QR
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                     omega(B), Tprop_tensors(B))
                call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_B(j), &
                     l1_tensor,l2_tensor, Bvec_tensors(B))
                call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                     omega(B), Bvec_tensors(B), Tprop_bar_tensors(B))
!               solve for -w
                if (omega(B).ne.ZERO) then 
                    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                         omega_min(B), Tprop_tensors_min_w(B))
                    call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_B(j), &
                         l1_tensor,l2_tensor, Bvec_tensors_min_w(B))
                    call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), &
                         omega_min(B), Bvec_tensors_min_w(B), Tprop_bar_tensors_min_w(B))
                end if 

            end if ! cc-cc or eom; do something for actual cc-ci ?

            call print_date('Solved All 1st response right and left equations for operator B')

            if (exa_input%do_1rdm_lr) then
            ! build/save density matrices here ?
            !
               call print_date('Saving perturbed density matrices for operator B(+w)')
               write (char_a,'(I1)'), i
               write (char_b,'(I1)'), j
               write (char_omega,'(I1)'), n
               call talsh_dm_set_1rdm_dump_filename('CR+'//char_omega//'_'//char_a//'_'//char_b)
               ! call talsh_dm_driver(exa_input, &
               !                   interm_t%t2,l2_tensor, &
               !                   interm_t%t1,l1_tensor, &
               !                   int_t, &
               !                   tX1_tensor=Tprop_tensors(B)%t1, &
               !                   lX1_tensor=Tprop_bar_tensors(B)%t1, &
               !                   tX1_tensor=Tprop_tensors(B)%t2, &
               !                   lX1_tensor=Tprop_bar_tensors(B)%t2)
               call talsh_dm_driver(exa_input, &
                                 interm_t%t2,l2_tensor, &
                                 t1_tensor=interm_t%t1, &
                                 l1_tensor=l1_tensor, &
                                 int_t=int_t, &
                                 tX1_tensor=Tprop_tensors(B)%t1, &
                                 tX2_tensor=Tprop_tensors(B)%t2)
               call talsh_dm_set_1rdm_dump_filename('CL+'//char_omega//'_'//char_a//'_'//char_b)
               call talsh_dm_driver(exa_input, &
                                 interm_t%t2,l2_tensor, &
                                 t1_tensor=interm_t%t1, &
                                 l1_tensor=l1_tensor, &
                                 int_t=int_t, &
                                 lX1_tensor=Tprop_bar_tensors(B)%t1, &
                                 lX2_tensor=Tprop_bar_tensors(B)%t2)
           
                if (omega(B).ne.ZERO) then
                    call print_date('Saving perturbed density matrices for operator B(-w)')
                    ! call talsh_dm_set_1rdm_dump_filename('cc-'//char_omega//'_'//char_a//'_'//char_b)
                    ! call talsh_dm_driver(exa_input, &
                    !              interm_t%t2,l2_tensor, &
                    !              interm_t%t1,l1_tensor, &
                    !              int_t, &
                    !              tX1_tensor=Tprop_tensors_min_w(B)%t1, &
                    !              lX1_tensor=Tprop_bar_tensors_min_w(B)%t1, &
                    !              tX2_tensor=Tprop_tensors_min_w(B)%t2, &
                    !              lX2_tensor=Tprop_bar_tensors_min_w(B)%t2)
                    call talsh_dm_set_1rdm_dump_filename('CR-'//char_omega//'_'//char_a//'_'//char_b)
                    call talsh_dm_driver(exa_input, &
                                 interm_t%t2,l2_tensor, &
                                 t1_tensor=interm_t%t1, &
                                 l1_tensor=l1_tensor, &
                                 int_t=int_t, &
                                 tX1_tensor=Tprop_tensors_min_w(B)%t1, &
                                 tX2_tensor=Tprop_tensors_min_w(B)%t2)
                    call talsh_dm_set_1rdm_dump_filename('CL-'//char_omega//'_'//char_a//'_'//char_b)
                    call talsh_dm_driver(exa_input, &
                                 interm_t%t2,l2_tensor, &
                                 t1_tensor=interm_t%t1, &
                                 l1_tensor=l1_tensor, &
                                 int_t=int_t, &
                                 lX1_tensor=Tprop_bar_tensors_min_w(B)%t1, &
                                 lX2_tensor=Tprop_bar_tensors_min_w(B)%t2)
               end if
            end if

            do i=1, N_prop_A 
                if (omega(B).eq.ZERO) then
                    call print_date('Note: Calculating static linear response.')
                    ! Tprop_*_tensors_min_w not used/initialized
                    call get_asym_lr(exa_input,interm_t, int_t, &
                            Interm_t_tensors_A(i),  & 
                            Tprop_tensors, Tprop_tensors, &
                            Tprop_bar_tensors, Tprop_bar_tensors, &
                            l1_tensor,l2_tensor,eom_lr)
                else
                    call get_asym_lr(exa_input,interm_t, int_t, &
                            Interm_t_tensors_A(i),  & 
                            Tprop_tensors, Tprop_tensors_min_w, &
                            Tprop_bar_tensors, Tprop_bar_tensors_min_w, &
                            l1_tensor,l2_tensor,eom_lr)
                end if 
                
                eom_lrs(n,i,j) = eom_lr
                write(*,"(A4,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), ">> = "
                write(*,"(F16.8,A,F16.8)"), real(eom_lrs(n,i,j)), ' \ ', aimag(eom_lrs(n,i,j))
            end do ! loop over A operators for a given B

! instead of destroy, should save these to file or to results variable for later reuse
                
            call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,2)

            end do ! loop over B operators
    end do ! omega 


    do i=1, N_prop_A
        call destroy_Interm_t_tensors(Interm_t_tensors_A(i),CC_CC)
    end do

    do i=1, N_prop_B
        call destroy_Interm_t_tensors(Interm_t_tensors_B(i),CC_CC)
    end do
    
    if (CC_CC) then 
        write(*,"(A54)"), "*************CC-CC(LR-CC) linear response*************"
    else
        write(*,"(A54)"), "*************CC-CI(EOM-CC) linear response************"
    end if
        
    do n=1, N_omega
        do i=1, N_prop_A
            do j=1, N_prop_B
                write(*,"(A4,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), ">> = "
                write(*,"(F16.8,A,F16.8)"), real(eom_lrs(n,i,j)), ' \ ', aimag(eom_lrs(n,i,j))
            end do
        end do
    end do


    DEALLOCATE(omega)
    DEALLOCATE(omega_min)
    DEALLOCATE(Tprop_tensors)
    DEALLOCATE(Tprop_tensors_min_w)
    DEALLOCATE(Tprop_bar_tensors)
    DEALLOCATE(Tprop_bar_tensors_min_w)
    DEALLOCATE(Interm_t_tensors_A)
    DEALLOCATE(Interm_t_tensors_B)
    DEALLOCATE(Bvec_tensors)
    DEALLOCATE(eom_lrs)

    write(*,"(A54)"), "**********************************************************"  

    if (CC_CC) then     
        call print_date("Finish asymmetric CC-CC Linear Response Calculations")
    else
        call print_date("Finish asymmetric CC-CI Linear Response Calculations")
    end if 


end subroutine talsh_ccsd_lr_asym_driver

!------------------------------------------------------------------------------------------------------------------

subroutine talsh_ccsd_qr_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,CC_CC)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    logical, intent(in)                         :: CC_CC
    

!   Wave-function tensor needed to construct the response function

    complex(kind=8)               :: omega1, omega2, omega3
    complex(kind=8),allocatable   :: omega(:), omega_min(:)
    complex(kind=8)               :: eom_qr
    complex(kind=8),allocatable   :: eom_qrs(:,:,:,:,:)
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR

    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR

    type(talsh_comm_tens_t), allocatable :: Bvec_tensors(:), Bvec_tensors_min_w(:)

    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:), Interm_t_tensors_B(:), Interm_t_tensors_C(:)
    integer                       :: N_prop, N_prop_A, N_prop_B, N_prop_C
    integer                       :: N_oper_lr, N_oper_qr
    integer                       :: N_omega1, N_omega2
    integer                       :: i, j, k, n, m
    integer, allocatable          :: Prop_A(:), Prop_B(:), Prop_C(:)

    type(C_PTR) :: body_p
    complex(8), pointer    ::prop_tens(:,:,:,:)
    integer            :: tens_rank
    integer :: ierr
    logical          :: exchange
    

    N_prop = 3

    call print_date('Do quadratic response')

    N_oper_qr = 3

    Prop_A = exa_input%Prop_A
    Prop_B = exa_input%Prop_B
    Prop_C = exa_input%Prop_C

    N_oper_qr = 3
    N_prop_A = size(Prop_A)
    N_prop_B = size(Prop_B)
    N_prop_C = size(Prop_C)

    allocate(Interm_t_tensors_A(N_prop_A))
    allocate(Interm_t_tensors_B(N_prop_B))
    allocate(Interm_t_tensors_C(N_prop_C))


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


    call print_date('Constructed intermediates needed in CC Quadratic Response')


    call do_cc_qr(exa_input,interm_t,int_t,l1_tensor,l2_tensor,&
                Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C, CC_CC)

end subroutine talsh_ccsd_qr_driver


!----------------------------------------------------------------------------------------------------------------------------------

subroutine do_cc_qr(exa_input,interm_t,int_t,l1_tensor,l2_tensor,&
                    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C, CC_CC)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)              :: l2_tensor
    type(talsh_tens_t), intent(inout)              :: l1_tensor
    type(talsh_rsp_interm_tens_t), allocatable     :: Interm_t_tensors_A(:) ! store exp values for (list of) A operators
    type(talsh_rsp_interm_tens_t), allocatable     :: Interm_t_tensors_B(:) ! store exp values for (list of) B operators
    type(talsh_rsp_interm_tens_t), allocatable     :: Interm_t_tensors_C(:) ! store exp values for (list of) C operators
    logical, intent(in)                            :: CC_CC

!   Wave-function tensor needed to construct the response function

    complex(kind=8)               :: omega1, omega2, omega3
    complex(kind=8),allocatable   :: omega(:), omega_min(:)
    complex(kind=8)               :: cc_qr
    complex(kind=8),allocatable   :: cc_qrs(:,:,:,:,:)
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR
    type(talsh_comm_tens_t), allocatable :: Bvec_tensors(:), Bvec_tensors_min_w(:)

    
    integer                       :: N_prop, N_prop_A, N_prop_B, N_prop_C
    integer                       :: N_oper_lr, N_oper_qr
    integer                       :: N_omega1, N_omega2
    integer                       :: i, j, k, n, m
    integer, allocatable          :: Prop_A(:), Prop_B(:), Prop_C(:)

    type(C_PTR) :: body_p
    complex(8), pointer    ::prop_tens(:,:,:,:)
    integer            :: tens_rank
    integer :: ierr
    logical          :: exchange
    

    N_prop = 3

    call print_date('Do quadratic response')

    N_oper_qr = 3

    Prop_A = exa_input%Prop_A
    Prop_B = exa_input%Prop_B
    Prop_C = exa_input%Prop_C

    N_oper_qr = 3
    N_prop_A = size(Prop_A)
    N_prop_B = size(Prop_B)
    N_prop_C = size(Prop_C)

    N_omega1 = exa_input%N_omega_qr1
    N_omega2 = exa_input%N_omega_qr2


    allocate(omega(N_oper_qr))
    allocate(omega_min(N_oper_qr))

    omega = ZERO 
    omega_min = ZERO 

    allocate(Bvec_tensors(N_oper_qr))
    allocate(Tprop_tensors(N_oper_qr))
    allocate(Tprop_tensors_min_w(N_oper_qr))
    allocate(Tprop_bar_tensors(N_oper_qr))
    allocate(Tprop_bar_tensors_min_w(N_oper_qr))
    allocate(cc_qrs(N_prop_A,N_prop_B,N_prop_C,N_omega2,N_omega1))


    do n=1, N_omega1
        omega1 = dcmplx(exa_input%w_qr1_re(n), exa_input%w_qr1_im(n))
        do m=1, N_omega2
            omega2 = dcmplx(exa_input%w_qr2_re(m), exa_input%w_qr2_im(m)) 
            omega3 = -(omega1+omega2)

            omega(1) = omega3   !omega(1): frequency for A operator
            omega(2) = omega1   !omega(2): frequency for B operator
            omega(3) = omega2   !omega(3): frequency for C operator
            print *, "Read frequency1:", omega(1)
            print *, "Read frequency2:", omega(2)
            print *, "Read frequency3:", omega(3)
      
      
            omega_min(1) = dcmplx(-real(omega(1)),aimag(omega(1)))
            omega_min(2) = dcmplx(-real(omega(2)),aimag(omega(2)))
            omega_min(3) = dcmplx(-real(omega(3)),aimag(omega(3)))
            print *, "Read frequency min_1:", omega_min(1)
            print *, "Read frequency min_2:", omega_min(2)
            print *, "Read frequency min_3:", omega_min(3)


! This loop is to determine which three properties to calculate
        do i=1, N_prop_A

            call print_lr_setup("R","A",i,omega(1))
            call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega(1), Tprop_tensors(1))
            if ((omega(2).ne.ZERO).or.(omega(3).ne.ZERO)) then 
                call print_lr_setup("R","A",i,omega_min(1))
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega_min(1), Tprop_tensors_min_w(1))
            end if 

            if (CC_CC) then
                call get_left_bvec_cccc(exa_input,int_t,interm_t,Interm_t_tensors_A(i), &
                                        l1_tensor,l2_tensor, Tprop_tensors(1), Bvec_tensors(1))
            else
                call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_A(i), &
                                   l1_tensor,l2_tensor, Bvec_tensors(1))
            end if

            call print_lr_setup("L","A",i,omega(1))
            if (CC_CC) then
                call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), &
                                       omega(1), Bvec_tensors(1), Tprop_bar_tensors(1))
            else
                call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), &
                                       omega(1), Bvec_tensors(1), Tprop_bar_tensors(1))
            end if 

            if ((omega(2).ne.ZERO).or.(omega(3).ne.ZERO)) then 
                call print_lr_setup("L","A",i,omega_min(1))
                if (CC_CC) then
                    call solve_cc_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), &
                                      omega_min(1), Bvec_tensors(1), Tprop_bar_tensors_min_w(1))
                else
                     call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), &
                                      omega_min(1), Bvec_tensors(1), Tprop_bar_tensors_min_w(1))
                end if 
            end if

            do j=1, N_prop_B
                
                call print_lr_setup("R","B",j,omega(2))
                call solve_eom_lr_right(exa_input, interm_t, int_t, &
                                        Interm_t_tensors_B(j), &
                                        omega(2), Tprop_tensors(2))
                if ((omega(2).ne.ZERO).or.(omega(3).ne.ZERO)) then 
                    call print_lr_setup("R","B",j,omega_min(2))
                    call solve_eom_lr_right(exa_input, interm_t, int_t, &
                           Interm_t_tensors_B(j), omega_min(2), Tprop_tensors_min_w(2))
                end if 

                if (CC_CC) then
                     call get_left_bvec_cccc(exa_input,int_t,interm_t,Interm_t_tensors_B(j), &
                                        l1_tensor,l2_tensor, Tprop_tensors(2), Bvec_tensors(2))
                else
                     call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_B(j), &
                                        l1_tensor,l2_tensor, Bvec_tensors(2))
                end if 

                call print_lr_setup("L","B",j,omega(2))
                if (CC_CC) then
                    call solve_cc_lr_left(exa_input, interm_t, int_t, &
                                           Interm_t_tensors_B(j), omega(2), &
                                           Bvec_tensors(2), Tprop_bar_tensors(2))
                else
                    call solve_eom_lr_left(exa_input, interm_t, int_t, &
                                           Interm_t_tensors_B(j), omega(2), &
                                           Bvec_tensors(2), Tprop_bar_tensors(2))
                end if

                if ((omega(2).ne.ZERO).or.(omega(3).ne.ZERO)) then 
                    call print_lr_setup("L","B",j,omega_min(2))
                    if (CC_CC) then
                        call solve_cc_lr_left(exa_input, interm_t, int_t, &
                                       Interm_t_tensors_B(j), omega_min(2),&
                                       Bvec_tensors(2), Tprop_bar_tensors_min_w(2))
                    else
                        call solve_eom_lr_left(exa_input, interm_t, int_t, &
                                       Interm_t_tensors_B(j), omega_min(2),&
                                       Bvec_tensors(2), Tprop_bar_tensors_min_w(2))
                    end if
                end if 

                do k=1, N_prop_C

                    call print_lr_setup("R","C",k,omega(3))
                    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_C(k), omega(3), Tprop_tensors(3))
                    if ((omega(2).ne.ZERO).or.(omega(3).ne.ZERO)) then 
                        call print_lr_setup("R","C",k,omega_min(3))
                        call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_C(k), &
                            omega_min(3), Tprop_tensors_min_w(3))
                    end if      

                    if (CC_CC) then
                        call get_left_bvec_cccc(exa_input,int_t,interm_t, &
                                           Interm_t_tensors_C(k), &
                                           l1_tensor,l2_tensor, &
                                           Tprop_tensors(3), Bvec_tensors(3))
                    else
                        call get_left_bvec(exa_input,int_t,interm_t, &
                                           Interm_t_tensors_C(k), &
                                           l1_tensor,l2_tensor, &
                                           Bvec_tensors(3))
                    end if

                    call print_lr_setup("L","C",k,omega(3))
                    if (CC_CC) then
                        call solve_cc_lr_left(exa_input, interm_t, int_t, &
                                               Interm_t_tensors_C(k), omega(3), &
                                               Bvec_tensors(3), Tprop_bar_tensors(3))
                    else
                        call solve_eom_lr_left(exa_input, interm_t, int_t, &
                                               Interm_t_tensors_C(k), omega(3), &
                                               Bvec_tensors(3), Tprop_bar_tensors(3))
                    end if 
                        
                    if ((omega(2).ne.ZERO).or.(omega(3).ne.ZERO)) then 
                        call print_lr_setup("L","C",k,omega_min(3))
                        if (CC_CC) then
                            call solve_cc_lr_left(exa_input, interm_t, int_t, &
                                           Interm_t_tensors_C(k), omega_min(3),&
                                           Bvec_tensors(3), Tprop_bar_tensors_min_w(3))
                        else
                            call solve_eom_lr_left(exa_input, interm_t, int_t, &
                                           Interm_t_tensors_C(k), omega_min(3),&
                                           Bvec_tensors(3), Tprop_bar_tensors_min_w(3))
                        end if
                    end if 

                    call print_date('Solved All 1st response equations') 

                    if ((omega(1).ne.ZERO).or.(omega(2).ne.ZERO)) then
                        if (CC_CC) then
                            call get_cccc_qr(exa_input, int_t, interm_t,&
                                    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,& 
                                    Tprop_tensors, Tprop_bar_tensors,& 
                                    Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                                    i, j, k, l1_tensor, l2_tensor, cc_qr)
                        else
                            call get_eom_qr(exa_input, interm_t,&
                                    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,& 
                                    Tprop_tensors, Tprop_bar_tensors,& 
                                    Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                                    i, j, k, l1_tensor, l2_tensor, cc_qr)
                        end if
                    else
                        call print_date('Note: Calculating static quadratic response.')
                        if (CC_CC) then
                            call get_cccc_qr(exa_input, int_t, interm_t,&
                                    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,& 
                                    Tprop_tensors, Tprop_bar_tensors,& 
                                    Tprop_tensors, Tprop_bar_tensors,&
                                    i, j, k, l1_tensor, l2_tensor, cc_qr)
                        else
                            call get_eom_qr(exa_input, interm_t,&
                                    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,& 
                                    Tprop_tensors, Tprop_bar_tensors,& 
                                    Tprop_tensors, Tprop_bar_tensors,&
                                    i, j, k, l1_tensor, l2_tensor, cc_qr)
                        end if
                end if 

                    cc_qrs(i,j,k,m,n) = cc_qr

                    call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,3,Bvec_tensors,&
                    Tprop_bar_tensors,Tprop_bar_tensors_min_w)

                end do
                call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,2,Bvec_tensors,&
                                            Tprop_bar_tensors,Tprop_bar_tensors_min_w)

            end do
            call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,1,Bvec_tensors,&
                                        Tprop_bar_tensors,Tprop_bar_tensors_min_w)

        end do
    
      end do
    end do 


        do i=1, N_prop_A
            call destroy_Interm_t_tensors(Interm_t_tensors_A(i),.false.)
        end do

        do i=1, N_prop_B
            call destroy_Interm_t_tensors(Interm_t_tensors_B(i),.false.)
        end do
        
        do i=1, N_prop_C
            call destroy_Interm_t_tensors(Interm_t_tensors_C(i),.false.)
        end do


        write(*,"(A54)"), "*************CC-CI(EOM-CC) quadratic response************"
        do n=1, N_omega1
          do m=1, N_omega2
            do i=1, N_prop_A
              do j=1, N_prop_B
                do k=1, N_prop_C
                    write(*,"(A4,I1,A2,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), &
                         ";X", Prop_C(k),">> = "

                    write(*,"(F16.8,A,F16.8)"), real(cc_qrs(i,j,k,m,n)), ' \ ', aimag(cc_qrs(i,j,k,m,n))
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
        DEALLOCATE(Bvec_tensors)
        DEALLOCATE(cc_qrs)


    call print_date('Leave CC response') 


end subroutine do_cc_qr

!------------------------------------------------------------------------------------------------------------------



end module
