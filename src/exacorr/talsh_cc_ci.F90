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
    !   Written by Xiang Yuan Janunary 2022, Lille, France

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    integer                                     :: rr ! rank of response 
    

!   Wave-function tensor needed to construct the response function

    complex(kind=8)               :: omega1, omega2, omega3
    type(talsh_rsp_interm_tens_t) :: interm_t_prop1, interm_t_prop2, interm_t_prop3  
    type(talsh_comm_tens_t)       :: tprop_tensor1, tprop_tensor2, tprop_tensor3
    type(talsh_comm_tens_t)       :: tprop_bar_tensor1, tprop_bar_tensor2, tprop_bar_tensor3
    type(talsh_comm_tens_t)       :: bvec_tensor1, bvec_tensor2, bvec_tensor3
    complex(kind=8)               :: eom_qr, eom_lr
    type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_minus_fre(:)
    type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:)
    type(talsh_comm_tens_t), allocatable :: Bvec_tensors(:)
    type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors(:)
    integer                       :: N_prop, i, j


    integer            :: vvoo_dims(4)
    type(C_PTR) :: body_p
    complex(8), pointer    ::prop_tens(:,:,:,:)
    integer            :: tens_rank
    integer :: ierr
    
    omega1=(0.0,0.0)
    omega1 = dcmplx(exa_input%w_lr_re, exa_input%w_lr_im)
    omega2=(0.0,0.0)
    omega3=(0.0,0.0)


    N_prop = 3
    allocate(Tprop_tensors(N_prop))
    allocate(Tprop_tensors_minus_fre(N_prop))
    allocate(Tprop_bar_tensors(N_prop))
    allocate(Bvec_tensors(N_prop))
    allocate(Interm_t_tensors(N_prop))

    call print_date('Enter CC-CI response')
    print *, "Read frequency:", omega1

    do  i=1, N_prop
        call form_eomccsd_qr_intermediates(exa_input, interm_t, i, Interm_t_tensors(i),&
                                  int_t, l1_tensor, l2_tensor)
    end do

    ! call print_tensor(Interm_t_tensors(3)%epsilon_prop_vo,1.0d-10,"vo part of epison_prop3_before")
    ! call print_tensor(Interm_t_tensors(3)%epsilon_prop_vvoo,1.0d-10,"T2 part of epison_prop3_before")


    ! ierr=talsh_tensor_dimensions(Interm_t_tensors(3)%epsilon_prop_vvoo,tens_rank,vvoo_dims)
    ! ierr=talsh_tensor_get_body_access(Interm_t_tensors(3)%epsilon_prop_vvoo,body_p,C8,int(0,C_INT),DEV_HOST)
    ! call c_f_pointer(body_p,prop_tens,vvoo_dims)


    ! call print_tensor(Interm_t_tensors(3)%epsilon_prop_vo,1.0d-10,"vo part of epison_prop3_after")
    ! call print_tensor(Interm_t_tensors(3)%epsilon_prop_vvoo,1.0d-10,"T2 part of epison_prop3_after")
    ! print *, "Norm^2 of Interm_t_tensors(3)%epsilon_prop_vo",  tensor_norm2(Interm_t_tensors(3)%epsilon_prop_vo)
    ! print *, "Norm^2 of Interm_t_tensors(3)%epsilon_prop_vo",  tensor_norm2(Interm_t_tensors(3)%epsilon_prop_vvoo)

    ! call print_tensor(Interm_t_tensors(3)%eta_prop_ov,1.0d-10,"vo part of eta_prop3_after")
    ! call print_tensor(Interm_t_tensors(3)%eta_prop_oovv,1.0d-10,"T2 part of eta_prop3_after")
    ! print *, "Norm^2 of Interm_t_tensors(3)%eta_prop_ov",  tensor_norm2(Interm_t_tensors(3)%eta_prop_ov)
    ! print *, "Norm^2 of Interm_t_tensors(3)%eta_prop_oovv",  tensor_norm2(Interm_t_tensors(3)%eta_prop_oovv)

    ! call print_tensor(Interm_t_tensors(3)%eta_prop_ov,1.0d-10,"vo part of eta_prop3_after")
    ! call print_tensor(Interm_t_tensors(3)%eta_prop_oovv,1.0d-10,"T2 part of eta_prop3_after")

    call print_date('Constructed intermediates needed in CC-CI response')


    
    do i=1, N_prop
        call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors(i), omega1, Tprop_tensors(i))
        call solve_eom_lr_right(exa_input, interm_t, int_t, Interm_t_tensors(i), -omega1, Tprop_tensors_minus_fre(i))
    end do
    call print_date('Solved All CC-CI response right equations')


    ! call print_tensor(Tprop_tensors(1)%t1,1.0d-10,"T1 part of T_prop1")
    ! call print_tensor(Tprop_tensors(1)%t2,1.0d-10,"T2 part of T_prop1")


    ! call print_tensor(Tprop_tensors(3)%t1,1.0d-10,"T1 part of T_prop3")
    ! call print_tensor(Tprop_tensors(3)%t2,1.0d-10,"T2 part of T_prop3")
    ! print *, "Norm^2 of Tprop_tensors(3)%t1",  tensor_norm2(Tprop_tensors(3)%t1)
    ! print *, "Norm^2 of Tprop_tensors(3)%t2",  tensor_norm2(Tprop_tensors(3)%t2)


    if (rr==1) then 
        call print_date('Do linear response')
        do i=1, N_prop
            do j=1, N_prop
                call get_eom_lr(exa_input,interm_t, int_t, Interm_t_tensors(i), Interm_t_tensors(j), & 
                Tprop_tensors(i), Tprop_tensors(j), Tprop_tensors_minus_fre(i), Tprop_tensors_minus_fre(j), & 
                    l1_tensor,l2_tensor,eom_lr)
            end do
        end do

    else 
        call print_date('Do quadratic response')
        do i=1, N_prop
            call get_left_bvec(exa_input,interm_t,Interm_t_tensors(i),l1_tensor,l2_tensor,&
                                Tprop_tensors(i),Bvec_tensors(i))
        end do

        do i=1, N_prop
            call solve_eom_lr_left(interm_t, int_t, Interm_t_tensors(i), omega1, Bvec_tensors(i), Tprop_bar_tensors(i))
        end do                         
    
        ! call get_eom_qr(exa_input,interm_t, interm_t_prop1, interm_t_prop2, interm_t_prop3,& 
        !             tprop_tensor1, tprop_tensor2, tprop_tensor3,& 
        !             tprop_bar_tensor1, tprop_bar_tensor2, tprop_bar_tensor3,&
        !             l1_tensor,l2_tensor,eom_qr)
    end if

    call print_date('Leave CC-CI response')    

end subroutine talsh_CC_CI_driver


!------------------------------------------------------------------------------------------------------------------


subroutine form_eomccsd_qr_intermediates(exa_input, interm_t, i_prop, interm_t_prop1, &
                                   int_t, l1_tensor, l2_tensor)
            use intermediates_eom_rsp
            use intermediates_eom_ccsd

                type(exacc_input), intent(in) :: exa_input
                type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
                integer, intent(in)                   :: i_prop 
                type(talsh_rsp_interm_tens_t) :: interm_t_prop1
                type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
                type(talsh_tens_t), intent(inout)           :: l2_tensor
                type(talsh_tens_t), intent(inout)           :: l1_tensor

                integer(INTD), dimension(4)   :: vvoo_dims, ovoo_dims, vvvo_dims, oovv_dims
                integer(INTD), dimension(2)   :: ov_dims, oo_dims, vv_dims, vo_dims

                type(talsh_orb_prop_integ_tens_t)   :: orb_int_tensor1
                type(talsh_tens_t)                  :: LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo
                integer                   :: nmo1          ! the length of the mo basis
                ! integer,allocatable       :: mo1_list(:)   ! and their indices
                integer                   :: nmo2          ! the length of the mo basis
                ! integer                   :: mo2_list(:)   ! and their indices
                
        
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


                ! eta is independent on property so that in eom-cc energy intermediates. in intermediates_eom.F90
                call get_eta_ov(interm_t%eta_ov,interm_t%fbar_ov,ov_dims)
                call get_eta_oovv(interm_t%eta_oovv,int_t%oovv,oovv_dims)

                ! call print_date("get eta_oovv?")

                ! get property integral in MO basis
                call get_prop_integral(exa_input,orb_int_tensor1,i_prop)

                ! call form_eom_qr_prop_inter(interm_t%t1, interm_t%t2, interm_t_prop1, interm_t_prop2, interm_t_prop3, &
                !                 l1_tensor, l2_tensor)

                ! get epsilon_prop  intermediates_eom_rsp.F90

                call construct_aux_tensor(LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

                call get_epsilon_prop_vo(interm_t_prop1%epsilon_prop_vo,orb_int_tensor1,&
                                        interm_t%t2,interm_t%t1,LR_Aux_oo,vo_dims)
                call get_epsilon_prop_vvoo(interm_t_prop1%epsilon_prop_vvoo,orb_int_tensor1,&
                                        interm_t%t2,interm_t%t1,LR_Aux_oo,LR_Aux_vv,vvoo_dims)


                ! call print_tensor(interm_t_prop1%epsilon_prop_vo, 1.0d-10, "Xi1")
                ! call print_tensor(interm_t_prop1%epsilon_prop_vvoo, 1.0d-10, "Xi2")
                
                ! print *, "NORM^2 of episilonX",  tensor_norm2(interm_t_prop1%epsilon_prop_vo) + &
                !                                 tensor_norm2(interm_t_prop1%epsilon_prop_vvoo)


                ! get O_prop
                call get_O_vo(interm_t_prop1%O_vo,orb_int_tensor1,vo_dims)

                ! get eta_prop
                call get_eta_prop_ov(interm_t_prop1%eta_prop_ov,orb_int_tensor1,l1_tensor,l2_tensor,&
                                    interm_t%t2,interm_t%t1,LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)
                call get_eta_prop_oovv(interm_t_prop1%eta_prop_oovv,orb_int_tensor1,l1_tensor,l2_tensor,&
                                    interm_t%t1,LR_Aux_oo,LR_Aux_vv,oovv_dims)
                call get_eta_prop_ov_eom(interm_t_prop1%eta_prop_ov_eom,interm_t_prop1%eta_prop_ov,&
                                    interm_t_prop1%epsilon_prop_vo,l2_tensor,ov_dims)


                ! get prop_bar
                call get_prop_bar_oo(interm_t_prop1%prop_bar_oo,orb_int_tensor1,interm_t%t1,oo_dims)
                call get_prop_bar_vv(interm_t_prop1%prop_bar_vv,orb_int_tensor1,interm_t%t1,vv_dims)
                call get_prop_bar_ov(interm_t_prop1%prop_bar_ov,orb_int_tensor1,ov_dims)


                ! get Q
                call get_Q_ovoo(interm_t_prop1%Q_ovoo,interm_t_prop1%prop_bar_ov,interm_t%t2,ovoo_dims)                               
                call get_Q_vvvo(interm_t_prop1%Q_vvvo,interm_t_prop1%prop_bar_ov,interm_t%t2,vvvo_dims)

                call destroy_aux_tensor(LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo)


end subroutine form_eomccsd_qr_intermediates


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



end module
