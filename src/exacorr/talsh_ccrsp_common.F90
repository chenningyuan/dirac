module talsh_ccrsp_common
!This module contains the routines needed to calculate CC-CI(EOM-CC) response property


    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public form_ccresp_intermediates
    public destroy_Interm_t_tensors
    public destroy_Tprop_tensors
    public destroy_Tprop_tensor
    public destroy_aux_tensor
    public construct_aux_tensor
    public get_prop_integral

    contains


!--------------------------------------------------------------------------------------------------------------------------

subroutine form_ccresp_intermediates(exa_input, interm_t, i_prop, interm_t_prop1, &
                                   int_t, l1_tensor, l2_tensor, CC_CC, exchange)
            use intermediates_eom_rsp
            use intermediates_eom_ccsd

                type(exacc_input), intent(in) :: exa_input
                type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
                integer, intent(in)                   :: i_prop 
                type(talsh_rsp_interm_tens_t), intent(inout):: interm_t_prop1
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
                logical, optional          :: exchange
                
        
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
                call get_prop_integral(exa_input,orb_int_tensor1,i_prop,exchange)

                interm_t_prop1%property_label = orb_int_tensor1%property_label

                ! get epsilon_prop  intermediates_eom_rsp.F90

                call construct_aux_tensor(LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

                call get_epsilon_prop_vo(interm_t_prop1%epsilon_prop_vo,orb_int_tensor1,&
                                        interm_t%t2,interm_t%t1,LR_Aux_oo,vo_dims)
                call get_epsilon_prop_vvoo(interm_t_prop1%epsilon_prop_vvoo,orb_int_tensor1,&
                                        interm_t%t2,interm_t%t1,LR_Aux_oo,LR_Aux_vv,vvoo_dims)



                ! call print_tensor(interm_t_prop1%epsilon_prop_vo, 1.0D-6, 'epsilon_prop_vo')
                ! call print_tensor(interm_t_prop1%epsilon_prop_vvoo, 1.0D-6, 'epsilon_prop_vvoo')


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

subroutine get_prop_integral(exa_input,orb_int_tensor,i_prop,exchange)
        type(exacc_input), intent(in)      :: exa_input
        type(talsh_orb_prop_integ_tens_t), intent(inout)  :: orb_int_tensor
        integer                            :: i_prop
        integer                            :: min_occ
        integer                            :: nocc, nvir
        logical, optional          :: exchange

        min_occ = minval(exa_input%mo_occ)
        nocc = exa_input%nocc
        nvir = exa_input%nvir

        call get_prop_integral_driver(orb_int_tensor%prop_oo,nocc,exa_input%mo_occ,&
                                     nocc,exa_input%mo_occ,min_occ,i_prop,exchange,orb_int_tensor%property_label)
        call get_prop_integral_driver(orb_int_tensor%prop_vv,nvir,exa_input%mo_vir,nvir,exa_input%mo_vir,min_occ,i_prop,exchange)
        call get_prop_integral_driver(orb_int_tensor%prop_ov,nocc,exa_input%mo_occ,nvir,exa_input%mo_vir,min_occ,i_prop,exchange)
        call get_prop_integral_driver(orb_int_tensor%prop_vo,nvir,exa_input%mo_vir,nocc,exa_input%mo_occ,min_occ,i_prop,exchange)

        ! call print_tensor(orb_int_tensor%prop_oo, 1.0D-6, 'orb_int_tensor%prop_oo')
        ! call print_tensor(orb_int_tensor%prop_vv, 1.0D-6, 'orb_int_tensor%prop_vv')
        ! call print_tensor(orb_int_tensor%prop_ov, 1.0D-6, 'orb_int_tensor%prop_ov')
        ! call print_tensor(orb_int_tensor%prop_vo, 1.0D-6, 'orb_int_tensor%prop_vo')

end subroutine 


!------------------------------------------------------------------------------------------------------------------
subroutine get_prop_integral_driver(orb_int_tensor,nmo1,mo1_list,nmo2,mo2_list,min_occ,i_prop,exchange,prop_label)
    ! Read proporty integrals from  MDPROP file 
    type(talsh_tens_t), intent(inout)        :: orb_int_tensor
    character*8, intent(inout), optional :: prop_label
    
    integer                :: nmo1          ! the length of the mo basis
    integer                :: mo1_list(:)   ! and their indices
    integer                :: nmo2          ! the length of the mo basis
    integer                :: mo2_list(:)   ! and their indices
    integer                :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

    type(C_PTR) :: body_p
    integer(C_INT) :: ierr
    integer                      :: nocc, nvir
    integer                      :: i_prop
    logical, optional          :: exchange

    complex(kind=8), pointer, contiguous :: orbital_integral_tens(:,:)

    integer                      :: i,j
    real(8)                      :: prop_real, prop_img, prop

    integer(INTD),dimension(2)   :: mo_dims


    mo_dims(1) = nmo1
    mo_dims(2) = nmo2


    ierr=talsh_tensor_construct(orb_int_tensor,C8,mo_dims,init_val=ZERO)
    ierr=talsh_tensor_get_body_access(orb_int_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,orbital_integral_tens,mo_dims)


    !call add_finitefield(orbital_integral_tens,exa_input,iff,mo1_list,nmo1,mo2_list,nmo2,min_occ)
    if (present(prop_label)) then 
       call get_ff_mat (orbital_integral_tens,i_prop,mo1_list,nmo1,mo2_list,nmo2,min_occ,p_label=prop_label)
    else
       call get_ff_mat (orbital_integral_tens,i_prop,mo1_list,nmo1,mo2_list,nmo2,min_occ)
    end if
    ! if (exchange) then
    ! do j=1, size(orbital_integral_tens,2)
    !     do i=1, size(orbital_integral_tens,1)
    !         prop_real = real(orbital_integral_tens(i,j))
    !         prop_img = aimag(orbital_integral_tens(i,j))
    !         prop = prop_img
    !         prop_img = prop_real
    !         prop_real = prop
    !         orbital_integral_tens(i,j) = dcmplx(prop_real,prop_img)
    !     end do
    ! end do
    ! end if
    ! call print_date('Get property integrals') 

end subroutine get_prop_integral_driver
!-------------------------------------------------------------------------------------------------------------------

subroutine construct_aux_tensor(LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

        integer(INTD), dimension(2)   :: ov_dims, oo_dims, vv_dims, vo_dims
        type(talsh_tens_t), intent(inout)            :: LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo
        integer(C_INT) :: ierr

        ierr=talsh_tensor_construct(LR_Aux_ov,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(LR_Aux_vo,C8,vo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(LR_Aux_vv,C8,vv_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(LR_Aux_oo,C8,oo_dims,init_val=ZERO)

end subroutine

!--------------------------------------------------------------------------------------------------------------------

subroutine destroy_aux_tensor(LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo)
    type(talsh_tens_t), intent(inout)      :: LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo
    integer(C_INT) :: ierr

    ierr=talsh_tensor_destruct(LR_Aux_oo)
    ierr=talsh_tensor_destruct(LR_Aux_vv)
    ierr=talsh_tensor_destruct(LR_Aux_ov)
    ierr=talsh_tensor_destruct(LR_Aux_vo)

end subroutine 

!------------------------------------------------------------------------------------------------------------------

subroutine destroy_Tprop_tensor(Tprop_tensor)
    type(talsh_comm_tens_t), intent(inout) :: Tprop_tensor
    
    integer(C_INT) :: ierr
    
    ierr=talsh_tensor_destruct(Tprop_tensor%t1)
    ierr=talsh_tensor_destruct(Tprop_tensor%t2)

end subroutine

!------------------------------------------------------------------------------------------------------------------

subroutine destroy_Tprop_tensors(Tprop_tensors_wB,    Tprop_tensors_min_wB, &
                                i, &
                                Bvec_tensors, Bvec_tensors_min_w, &
                                Tprop_bar_tensors_wB, Tprop_bar_tensors_min_wB, &
                                Tprop_tensors_wC,     Tprop_tensors_min_wC, &
                                Tprop_bar_tensors_wC, Tprop_bar_tensors_min_wC   )
    type(talsh_comm_tens_t), intent(inout), allocatable :: Tprop_tensors_wB(:)     !LR
    type(talsh_comm_tens_t), intent(inout), allocatable :: Tprop_tensors_min_wB(:) !LR
    integer        :: i
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Bvec_tensors(:)
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Bvec_tensors_min_w(:)
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_wB(:)
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_min_wB(:)
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_tensors_wC(:)  !QR
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_tensors_min_wC(:) !QR
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_wC(:) !QR
    type(talsh_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_min_wC(:) !QR
    
    integer(C_INT) :: ierr
    
    ierr=talsh_tensor_destruct(Tprop_tensors_wB(i)%t1)
    ierr=talsh_tensor_destruct(Tprop_tensors_wB(i)%t2)
    if (present(Tprop_bar_tensors_wB)) then 
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_wB(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_wB(i)%t2)
    end if

    ierr=talsh_tensor_destruct(Tprop_tensors_min_wB(i)%t1)
    ierr=talsh_tensor_destruct(Tprop_tensors_min_wB(i)%t2)
    if (present(Tprop_bar_tensors_min_wB)) then
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_min_wB(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_min_wB(i)%t2)
    end if

    if (present(Bvec_tensors)) then
        ierr=talsh_tensor_destruct(Bvec_tensors(i)%t1)
        ierr=talsh_tensor_destruct(Bvec_tensors(i)%t2)
    end if

    if (present(Bvec_tensors_min_w)) then
        ierr=talsh_tensor_destruct(Bvec_tensors_min_w(i)%t1)
        ierr=talsh_tensor_destruct(Bvec_tensors_min_w(i)%t2)
    end if

    if (present(Tprop_tensors_wC)) then
        ierr=talsh_tensor_destruct(Tprop_tensors_wC(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_tensors_wC(i)%t2)
    end if
    if (present(Tprop_bar_tensors_wC)) then 
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_wC(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_wC(i)%t2)
    end if

    if (present(Tprop_tensors_min_wC)) then
        ierr=talsh_tensor_destruct(Tprop_tensors_min_wC(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_tensors_min_wC(i)%t2)
    end if
    if (present(Tprop_bar_tensors_min_wC)) then
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_min_wC(i)%t1)
        ierr=talsh_tensor_destruct(Tprop_bar_tensors_min_wC(i)%t2)
    end if


end subroutine 


!---------------------------------------------------------------------------------------------------------------

subroutine destroy_Interm_t_tensors(interm_t_prop1, CC_CC)
    type(talsh_rsp_interm_tens_t), intent(inout) :: interm_t_prop1
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
