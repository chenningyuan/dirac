module talsh_cc_ci
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

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public talsh_CC_CI_driver

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine talsh_CC_CI_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor,rr)
!This subroutine drives the CC-CI response property calculation 
!it now includes linear and quadratic response
!Written by Xiang Yuan Janunary 2022, Lille, France

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    integer                                     :: rr ! rank of response 
    

!   Wave-function tensor needed to construct the response function

    complex(kind=8)               :: omega1, omega2, omega3
    complex(kind=8),allocatable   :: omega(:), omega_min(:)
    complex(kind=8)               :: eom_qr, eom_lr
    complex(kind=8),allocatable   :: eom_qrs(:,:,:), eom_lrs(:,:,:)
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR

    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR

    type(talsh_comm_tens_t), allocatable :: Bvec_tensors(:), Bvec_tensors_min_w(:)
    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:), Interm_t_tensors_B(:), Interm_t_tensors_C(:)
    integer                       :: N_prop, N_prop_A, N_prop_B, N_prop_C
    integer                       :: N_oper_lr, N_oper_qr
    integer                       :: N_omega
    integer                       :: i, j, k, n, m
    integer, allocatable          :: Prop_A(:), Prop_B(:), Prop_C(:)


    integer            :: vvoo_dims(4)
    type(C_PTR) :: body_p
    complex(8), pointer    ::prop_tens(:,:,:,:)
    integer            :: tens_rank
    integer :: ierr
    logical          :: CC_CC
    

    N_prop = 3


    if (exa_input%wf_type == 1) then 
        CC_CC = .true.
        call print_date('Activate CC-CC(LR-CC) response')
    else 
        CC_CC = .false.
        call print_date('Activate CC-CI(EOM-CC) response')
    end if 



    if (rr==1) then 
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
            call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_A(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC)
        end do
    
    
        do  i=1, N_prop_B
            m = Prop_B(i)
            call form_ccresp_intermediates(exa_input, interm_t, m, Interm_t_tensors_B(i),&
                                      int_t, l1_tensor, l2_tensor, CC_CC)
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
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega(1), Tprop_tensors(1))
                call solve_eom_lr_right(exa_input, interm_t, int_t, &
                                    Interm_t_tensors_A(i), omega_min(1), Tprop_tensors_min_w(1))

                do j=1, N_prop_B
                    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), omega(2), Tprop_tensors(2))
                    call solve_eom_lr_right(exa_input, interm_t, int_t, &
                                            Interm_t_tensors_B(j), omega_min(2), Tprop_tensors_min_w(2))

                    call print_date('Solved All 1st response right equations')

                    call get_eom_lr(exa_input,interm_t, int_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,& 
                        Tprop_tensors, Tprop_tensors_min_w, i,j,& 
                        l1_tensor,l2_tensor,eom_lr,CC_CC)

                    eom_lrs(n,i,j) = eom_lr
                    write(*,"(A4,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), ">> = "
                    write(*,"(F16.8,A2,F16.8)"), real(eom_lrs(n,i,j)), '\', aimag(eom_lrs(n,i,j))

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
                    write(*,"(F16.8,A2,F16.8)"), real(eom_lrs(n,i,j)), '\', aimag(eom_lrs(n,i,j))
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


    else 

        call print_date('Do quadratic response')

        N_oper_qr = 3

        Prop_A = exa_input%Prop_A
        Prop_B = exa_input%Prop_B
        Prop_C = exa_input%Prop_C

        N_oper_qr = 3
        N_prop_A = size(Prop_A)
        N_prop_B = size(Prop_B)
        N_prop_C = size(Prop_C)

        allocate(omega(N_oper_qr))
        allocate(omega_min(N_oper_qr))
        allocate(Interm_t_tensors_A(N_prop_A))
        allocate(Interm_t_tensors_B(N_prop_B))
        allocate(Interm_t_tensors_C(N_prop_C))
        allocate(Bvec_tensors(N_oper_qr))
        allocate(Tprop_tensors(N_oper_qr))
        allocate(Tprop_tensors_min_w(N_oper_qr))
        allocate(Tprop_bar_tensors(N_oper_qr))
        allocate(Tprop_bar_tensors_min_w(N_oper_qr))
        allocate(eom_qrs(N_prop_A,N_prop_B,N_prop_C))



        omega1=(0.0,0.0)
        omega2=(-0.0,0.0)
        omega3 = -(omega1+omega2)
        omega(1) = omega1
        omega(2) = omega2
        omega(3) = omega3
        print *, "Read frequency1:", omega(1)
        print *, "Read frequency2:", omega(2)
        print *, "Read frequency3:", omega(3)


        omega_min(1) = dcmplx(-real(omega(1)),aimag(omega(1)))
        omega_min(2) = dcmplx(-real(omega(2)),aimag(omega(2)))
        omega_min(3) = dcmplx(-real(omega(3)),aimag(omega(3)))
        print *, "Read frequency min_1:", omega_min(1)
        print *, "Read frequency min_2:", omega_min(2)
        print *, "Read frequency min_3:", omega_min(3)

        
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

        
! This loop is to determine which three properties to calculate
        do i=1, N_prop_A
            call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_A(i),l1_tensor,l2_tensor, &
                            Bvec_tensors(1))
            call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega(1), Tprop_tensors(1))
            ! call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), -omega(1), Tprop_tensors_min_w(1))
            call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega_min(1), Tprop_tensors_min_w(1))
            call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega(1), &
                                Bvec_tensors(1), Tprop_bar_tensors(1))
            ! call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), -omega(1),&
            !                     Bvec_tensors(1), Tprop_bar_tensors_min_w(1))

            call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_A(i), omega_min(1),&
                                Bvec_tensors(1), Tprop_bar_tensors_min_w(1))

            do j=1, N_prop_B
                call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_B(j),l1_tensor,l2_tensor, &
                                Bvec_tensors(2))
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), omega(2), Tprop_tensors(2))
                ! call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), -omega(2), Tprop_tensors_min_w(2))
                call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_B(j), omega_min(2), Tprop_tensors_min_w(2))
                call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), omega(2), &
                                Bvec_tensors(2), Tprop_bar_tensors(2))
                ! call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), -omega(2),&
                !                 Bvec_tensors(2), Tprop_bar_tensors_min_w(2))
                call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_B(j), omega_min(2),&
                                Bvec_tensors(2), Tprop_bar_tensors_min_w(2))

                do k=1, N_prop_C
                    call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors_C(k),l1_tensor,l2_tensor, &
                                Bvec_tensors(3))
                    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_C(k), omega(3), Tprop_tensors(3))
                    ! call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_C(k), -omega(3), Tprop_tensors_min_w(3))
                    call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors_C(k), omega_min(3), Tprop_tensors_min_w(3))
                    call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_C(k), omega(3), &
                                Bvec_tensors(3), Tprop_bar_tensors(3))
                    ! call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_C(k), -omega(3),&
                    !             Bvec_tensors(3), Tprop_bar_tensors_min_w(3))

                    call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors_C(k), omega_min(3),&
                                Bvec_tensors(3), Tprop_bar_tensors_min_w(3))
                    call print_date('Solved All 1st response equations') 


                    call get_eom_qr(exa_input, interm_t,&
                                    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,& 
                                    Tprop_tensors, Tprop_bar_tensors,& 
                                    Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                                    i, j, k, l1_tensor, l2_tensor, eom_qr)
                    
                    eom_qrs(i,j,k) = eom_qr

                    call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,3,Bvec_tensors,&
                    Tprop_bar_tensors,Tprop_bar_tensors_min_w)

                        !CC-CI left equation
                        ! call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors(m),l1_tensor,l2_tensor, &
                        !                     Tprop_tensors(n),Bvec_tensors(n))
                        ! call get_left_bvec(exa_input,int_t,interm_t,Interm_t_tensors(m),l1_tensor,l2_tensor, &
                        !                     Tprop_tensors_min_w(n),Bvec_tensors_min_w(n))
                        ! call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors(m), omega(n), &
                        !                         Bvec_tensors(n), Tprop_bar_tensors(n))
                        ! call solve_eom_lr_left(exa_input, interm_t, int_t, Interm_t_tensors(m), -omega(n),&
                        !                         Bvec_tensors_min_w(n), Tprop_bar_tensors_min_w(n))
                        ! ierr=talsh_tensor_destruct(Bvec_tensors_min_w(n)%t1)  for CC-CI
                        ! ierr=talsh_tensor_destruct(Bvec_tensors_min_w(n)%t2)

                end do
                call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,2,Bvec_tensors,&
                                            Tprop_bar_tensors,Tprop_bar_tensors_min_w)

            end do
            call destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,1,Bvec_tensors,&
                                        Tprop_bar_tensors,Tprop_bar_tensors_min_w)

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


        write(*,"(A54)"), "*************CC-CI(EOM-CC) quadratic response************"
        do i=1, N_prop_A
            do j=1, N_prop_B
                do k=1, N_prop_C
                    write(*,"(A4,I1,A2,I1,A2,I1,A5)"), "<<X", Prop_A(i), ";X", Prop_B(j), &
                         ";X", Prop_C(k),">> = "

                    write(*,"(F16.8,A2,F16.8)"), real(eom_qrs(i,j,k)), '\', aimag(eom_qrs(i,j,k))
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
        DEALLOCATE(eom_qrs)

    end if

    call print_date('Leave CC response')    

end subroutine talsh_CC_CI_driver


!------------------------------------------------------------------------------------------------------------------


subroutine form_ccresp_intermediates(exa_input, interm_t, i_prop, interm_t_prop1, &
                                   int_t, l1_tensor, l2_tensor, CC_CC)
            use intermediates_eom_rsp
            use intermediates_eom_ccsd

                type(exacc_input), intent(in) :: exa_input
                type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
                integer, intent(in)                   :: i_prop 
                type(talsh_rsp_interm_tens_t) :: interm_t_prop1
                type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
                type(talsh_tens_t), intent(inout)           :: l2_tensor
                type(talsh_tens_t), intent(inout)           :: l1_tensor
                ! CHARACTER, optional           :: WF

                integer(INTD), dimension(4)   :: vvoo_dims, ovoo_dims, vvvo_dims, oovv_dims
                integer(INTD), dimension(2)   :: ov_dims, oo_dims, vv_dims, vo_dims

                type(talsh_orb_prop_integ_tens_t)   :: orb_int_tensor1
                type(talsh_tens_t)                  :: LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo
                integer                   :: nmo1          ! the length of the mo basis
                ! integer,allocatable       :: mo1_list(:)   ! and their indices
                integer                   :: nmo2          ! the length of the mo basis
                ! integer                   :: mo2_list(:)   ! and their indices
                logical :: CC_CC
                
        
                oovv_dims(1:2) = exa_input%nocc
                oovv_dims(3:4) = exa_input%nvir

                vvoo_dims(1:2) = exa_input%nvir
                vvoo_dims(3:4) = exa_input%nocc

                ovoo_dims(1) = exa_input%nocc
                ovoo_dims(2) = exa_input%nvir
                ovoo_dims(3:4) = exa_input%nocc

                vvvo_dims(1:3) = exa_input%nvir
                vvvo_dims(4) = exa_input%nocc

                ov_dims(1) = exa_input%nocc
                ov_dims(2) = exa_input%nvir
                vo_dims(1) = exa_input%nvir
                vo_dims(2) = exa_input%nocc

                oo_dims(1:2) = exa_input%nocc
                vv_dims(1:2) = exa_input%nvir


                ! get property integral in MO basis
                call get_prop_integral(exa_input,orb_int_tensor1,i_prop)

                ! get epsilon_prop  intermediates_eom_rsp.F90

                call construct_aux_tensor(LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

                call get_epsilon_prop_vo(interm_t_prop1%epsilon_prop_vo,orb_int_tensor1,&
                                        interm_t%t2,interm_t%t1,LR_Aux_oo,vo_dims)
                call get_epsilon_prop_vvoo(interm_t_prop1%epsilon_prop_vvoo,orb_int_tensor1,&
                                        interm_t%t2,interm_t%t1,LR_Aux_oo,LR_Aux_vv,vvoo_dims)


                ! get O_prop
                call get_O_vo(interm_t_prop1%O_vo,orb_int_tensor1,vo_dims)

                ! get eta_prop

                if (CC_CC) then
                    call get_eta_prop_ov(interm_t_prop1%eta_prop_ov,orb_int_tensor1,l1_tensor,l2_tensor,&
                                    interm_t%t2,interm_t%t1,LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)
                else 
                    call get_eta_prop_ov(interm_t_prop1%eta_prop_ov,orb_int_tensor1,l1_tensor,l2_tensor,&
                                    interm_t%t2,interm_t%t1,LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)
                    call get_eta_prop_ov_eom(interm_t_prop1%eta_prop_ov_eom,interm_t_prop1%eta_prop_ov,&
                                    interm_t_prop1%epsilon_prop_vo,l2_tensor,ov_dims)
                end if 

                call get_eta_prop_oovv(interm_t_prop1%eta_prop_oovv,orb_int_tensor1,l1_tensor,l2_tensor,&
                                    interm_t%t1,LR_Aux_oo,LR_Aux_vv,oovv_dims)

                ! get prop_bar
                call get_prop_bar_oo(interm_t_prop1%prop_bar_oo,orb_int_tensor1,interm_t%t1,oo_dims)
                call get_prop_bar_vv(interm_t_prop1%prop_bar_vv,orb_int_tensor1,interm_t%t1,vv_dims)
                call get_prop_bar_ov(interm_t_prop1%prop_bar_ov,orb_int_tensor1,ov_dims)


                ! get Q
                ! call get_Q_ovoo(interm_t_prop1%Q_ovoo,interm_t_prop1%prop_bar_ov,interm_t%t2,ovoo_dims)                               
                ! call get_Q_vvvo(interm_t_prop1%Q_vvvo,interm_t_prop1%prop_bar_ov,interm_t%t2,vvvo_dims)

                ! eta is independent on property so that in eom-cc energy intermediates. in intermediates_eom.F90
                ! call get_eta_ov(interm_t%eta_ov,interm_t%fbar_ov,ov_dims)
                ! call get_eta_oovv(interm_t%eta_oovv,int_t%oovv,oovv_dims)

                call destroy_aux_tensor(LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo)


end subroutine form_ccresp_intermediates


!------------------------------------------------------------------------------------------------------------------

subroutine get_prop_integral(exa_input,orb_int_tensor,i_prop)
        type(exacc_input), intent(in)      :: exa_input
        type(talsh_orb_prop_integ_tens_t)  :: orb_int_tensor
        integer                            :: i_prop
        integer                            :: min_occ
        integer                            :: nocc, nvir

        min_occ = minval(exa_input%mo_occ)
        nocc = exa_input%nocc
        nvir = exa_input%nvir

        call get_prop_integral_driver(orb_int_tensor%prop_oo,nocc,exa_input%mo_occ,nocc,exa_input%mo_occ,min_occ,i_prop)
        call get_prop_integral_driver(orb_int_tensor%prop_vv,nvir,exa_input%mo_vir,nvir,exa_input%mo_vir,min_occ,i_prop)
        call get_prop_integral_driver(orb_int_tensor%prop_ov,nocc,exa_input%mo_occ,nvir,exa_input%mo_vir,min_occ,i_prop)
end subroutine 


!------------------------------------------------------------------------------------------------------------------
subroutine get_prop_integral_driver(orb_int_tensor,nmo1,mo1_list,nmo2,mo2_list,min_occ,i_prop)
    ! Read proporty integrals from  MDPROP file 
    type(talsh_tens_t)        :: orb_int_tensor
    
    integer                :: nmo1          ! the length of the mo basis
    integer                :: mo1_list(:)   ! and their indices
    integer                :: nmo2          ! the length of the mo basis
    integer                :: mo2_list(:)   ! and their indices
    integer                :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

    type(C_PTR) :: body_p
    integer(C_INT) :: ierr
    integer                      :: nocc, nvir
    integer                      :: i_prop
    complex(kind=8), pointer, contiguous :: orbital_integral_tens(:,:)


    integer(INTD),dimension(2)   :: mo_dims


    mo_dims(1) = nmo1
    mo_dims(2) = nmo2


    ierr=talsh_tensor_construct(orb_int_tensor,C8,mo_dims,init_val=ZERO)
    ierr=talsh_tensor_get_body_access(orb_int_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,orbital_integral_tens,mo_dims)


    !call add_finitefield(orbital_integral_tens,exa_input,iff,mo1_list,nmo1,mo2_list,nmo2,min_occ)
    call get_ff_mat (orbital_integral_tens,i_prop,mo1_list,nmo1,mo2_list,nmo2,min_occ)
    ! call print_date('Get property integrals')    


end subroutine get_prop_integral_driver
!-------------------------------------------------------------------------------------------------------------------

subroutine construct_aux_tensor(LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

        integer(INTD), dimension(2)   :: ov_dims, oo_dims, vv_dims, vo_dims
        type(talsh_tens_t)            :: LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo
        integer(C_INT) :: ierr

        ierr=talsh_tensor_construct(LR_Aux_ov,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(LR_Aux_vo,C8,vo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(LR_Aux_vv,C8,vv_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(LR_Aux_oo,C8,oo_dims,init_val=ZERO)

end subroutine

!--------------------------------------------------------------------------------------------------------------------

subroutine destroy_aux_tensor(LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo)
    type(talsh_tens_t)            :: LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo
    integer(C_INT) :: ierr

    ierr=talsh_tensor_destruct(LR_Aux_oo)
    ierr=talsh_tensor_destruct(LR_Aux_vv)
    ierr=talsh_tensor_destruct(LR_Aux_ov)
    ierr=talsh_tensor_destruct(LR_Aux_vo)

end subroutine 

!------------------------------------------------------------------------------------------------------------------

subroutine destroy_Tprop_tensors(Tprop_tensors,Tprop_tensors_min_w,i,Bvec_tensors,&
                                Tprop_bar_tensors,Tprop_bar_tensors_min_w)
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
    integer        :: i
    type(talsh_comm_tens_t), allocatable, OPTIONAL :: Bvec_tensors(:)
    type(talsh_comm_tens_t), allocatable, OPTIONAL :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR                        
    
    integer(C_INT) :: ierr
    
    ierr=talsh_tensor_destruct(Tprop_tensors(i)%t1)
    ierr=talsh_tensor_destruct(Tprop_tensors(i)%t2)
    ierr=talsh_tensor_destruct(Tprop_tensors_min_w(i)%t1)
    ierr=talsh_tensor_destruct(Tprop_tensors_min_w(i)%t2)

    if (present(Bvec_tensors)) then
        ierr=talsh_tensor_destruct(Bvec_tensors(i)%t1)
        ierr=talsh_tensor_destruct(Bvec_tensors(i)%t2)
    end if
    if (present(Tprop_bar_tensors)) then 
        ierr=talsh_tensor_destruct(Tprop_bar_tensors(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_bar_tensors(i)%t2)
    end if
    if (present(Tprop_bar_tensors_min_w)) then
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_min_w(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_min_w(i)%t2)
    end if


end subroutine 


!---------------------------------------------------------------------------------------------------------------

subroutine destroy_Interm_t_tensors(interm_t_prop1, CC_CC)
    type(talsh_rsp_interm_tens_t) :: interm_t_prop1
    logical          :: CC_CC
    integer          :: ierr

    ierr=talsh_tensor_destruct(interm_t_prop1%epsilon_prop_vo)
    ierr=talsh_tensor_destruct(interm_t_prop1%epsilon_prop_vvoo)

    ierr=talsh_tensor_destruct(interm_t_prop1%eta_prop_oovv)
    ierr=talsh_tensor_destruct(interm_t_prop1%eta_prop_ov)

    if(.not.CC_CC) then
        ierr=talsh_tensor_destruct(interm_t_prop1%eta_prop_ov_eom)
    end if

    ierr=talsh_tensor_destruct(interm_t_prop1%prop_bar_oo)
    ierr=talsh_tensor_destruct(interm_t_prop1%prop_bar_ov)    
    ierr=talsh_tensor_destruct(interm_t_prop1%prop_bar_vv)
    ierr=talsh_tensor_destruct(interm_t_prop1%O_vo)

end subroutine



end module
