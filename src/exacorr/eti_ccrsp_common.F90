module eti_ccrsp_common
!This module contains the routines needed to calculate CC-CI(EOM-CC) response property


    use tensor_algebra
    use exacorr_tensor_interface
    use exacorr_datatypes
    use exacorr_scaling_parameters
    use exacorr_cfg
    use exacorr_utils
    use eti_common_routines
    use exacorr_global

    implicit none

    type(tensor_t) :: one_tensor

    private

    public form_ccresp_intermediates
    public destroy_Interm_t_tensors
    public destroy_Tprop_tensors
    public destroy_Tprop_tensor
    public destroy_aux_tensor
    public get_prop_integral

    contains


!--------------------------------------------------------------------------------------------------------------------------

subroutine form_ccresp_intermediates(exa_input, interm_t, i_prop, interm_t_prop1, &
                                   int_t, CC_CC, exchange)
            use intermediates_eom_rsp
            use intermediates_eom_ccsd

                type(exacc_input), intent(in) :: exa_input
                type(eti_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
                integer, intent(in)                   :: i_prop 
                type(eti_rsp_interm_tens_t), intent(inout):: interm_t_prop1
                type(eti_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
                ! CHARACTER, optional           :: WF

                integer(INTD), dimension(4)   :: vvoo_dims, ovoo_dims, vvvo_dims, oovv_dims
                integer(INTD), dimension(2)   :: ov_dims, oo_dims, vv_dims, vo_dims

                type(eti_orb_prop_integ_tens_t)   :: orb_int_tensor1
                type(tensor_t)                  :: LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo
                integer                   :: nmo1          ! the length of the mo basis
                ! integer,allocatable       :: mo1_list(:)   ! and their indices
                integer                   :: nmo2          ! the length of the mo basis
                ! integer                   :: mo2_list(:)   ! and their indices
                logical :: CC_CC
                logical, optional          :: exchange
                type(eti_precision_t) :: eti_precision

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

                call construct_aux_tensor(exa_input%eti_precision,LR_Aux_ov,LR_Aux_vv,&
                                          LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

                print *,'a-a'
                call get_epsilon_prop_vo(interm_t_prop1%epsilon_prop_vo,orb_int_tensor1,&
                                        interm_t%gs%r2,interm_t%gs%r1,LR_Aux_oo,vo_dims)
                print *,'a-b'
                call get_epsilon_prop_vvoo(interm_t_prop1%epsilon_prop_vvoo,orb_int_tensor1,&
                                        interm_t%gs%r2,interm_t%gs%r1,LR_Aux_oo,LR_Aux_vv,vvoo_dims)
                print *,'a-c'



                ! call print_tensor(interm_t_prop1%epsilon_prop_vo, 1.0D-6, 'epsilon_prop_vo')
                ! call print_tensor(interm_t_prop1%epsilon_prop_vvoo, 1.0D-6, 'epsilon_prop_vvoo')


                ! get O_prop
                call get_O_vo(interm_t_prop1%O_vo,orb_int_tensor1,vo_dims)
                print *,'a-d'

                ! get eta_prop

                call get_eta_prop_ov(interm_t_prop1%eta_prop_ov,orb_int_tensor1,&
                                     interm_t%gs%l1, interm_t%gs%l2, &
                                     interm_t%gs%r2, interm_t%gs%r1, &
                                     LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)
                print *,'a-e'
                if (.not.CC_CC) then
                    call get_eta_prop_ov_eom(interm_t_prop1%eta_prop_ov_eom,interm_t_prop1%eta_prop_ov,&
                                    interm_t_prop1%epsilon_prop_vo,interm_t%gs%l2,ov_dims)
                print *,'a-f'
                end if 

                call get_eta_prop_oovv(interm_t_prop1%eta_prop_oovv,orb_int_tensor1,&
                                       interm_t%gs%l1, interm_t%gs%l2, & 
                                       interm_t%gs%r1,LR_Aux_oo,LR_Aux_vv,oovv_dims)
                print *,'a-g'

                ! get prop_bar
                call get_prop_bar_oo(interm_t_prop1%prop_bar_oo,orb_int_tensor1,interm_t%gs%r1,oo_dims)
                print *,'a-h'
                call get_prop_bar_vv(interm_t_prop1%prop_bar_vv,orb_int_tensor1,interm_t%gs%r1,vv_dims)
                print *,'a-i'
                call get_prop_bar_ov(interm_t_prop1%prop_bar_ov,orb_int_tensor1,ov_dims)
                print *,'a-j'


                ! get Q
                ! call get_Q_ovoo(interm_t_prop1%Q_ovoo,interm_t_prop1%prop_bar_ov,interm_t%gs%r2,ovoo_dims)                               
                ! call get_Q_vvvo(interm_t_prop1%Q_vvvo,interm_t_prop1%prop_bar_ov,interm_t%gs%r2,vvvo_dims)

                ! eta is independent on property so that in eom-cc energy intermediates. in intermediates_eom.F90
                ! call get_eta_ov(interm_t%eta_ov,interm_t%fbar_ov,ov_dims)
                ! call get_eta_oovv(interm_t%eta_oovv,int_t%oovv,oovv_dims)

                call destroy_aux_tensor(LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo)
                print *,'a-k'


end subroutine form_ccresp_intermediates


!------------------------------------------------------------------------------------------------------------------

subroutine get_prop_integral(exa_input,orb_int_tensor,i_prop,exchange)
        type(exacc_input), intent(in)      :: exa_input
        type(eti_orb_prop_integ_tens_t), intent(inout)  :: orb_int_tensor
        integer                            :: i_prop
        integer                            :: min_occ
        integer                            :: nocc, nvir
        logical, optional          :: exchange

        min_occ = minval(exa_input%mo_occ)
        nocc = exa_input%nocc
        nvir = exa_input%nvir

        call get_prop_integral_driver(orb_int_tensor%prop_oo,nocc,exa_input%mo_occ,nocc,&
                                      exa_input%mo_occ,min_occ,i_prop,exchange,         &
                                      exa_input%eti_precision,                          &
                                      orb_int_tensor%property_label)
        call get_prop_integral_driver(orb_int_tensor%prop_vv,nvir,exa_input%mo_vir,nvir,&
                                      exa_input%mo_vir,min_occ,i_prop,exchange,         &
                                      exa_input%eti_precision)
        call get_prop_integral_driver(orb_int_tensor%prop_ov,nocc,exa_input%mo_occ,nvir,&
                                      exa_input%mo_vir,min_occ,i_prop,exchange,         &
                                      exa_input%eti_precision)
        call get_prop_integral_driver(orb_int_tensor%prop_vo,nvir,exa_input%mo_vir,nocc,&
                                      exa_input%mo_occ,min_occ,i_prop,exchange,         &
                                      exa_input%eti_precision)

        ! call print_tensor(orb_int_tensor%prop_oo, 1.0D-6, 'orb_int_tensor%prop_oo')
        ! call print_tensor(orb_int_tensor%prop_vv, 1.0D-6, 'orb_int_tensor%prop_vv')
        ! call print_tensor(orb_int_tensor%prop_ov, 1.0D-6, 'orb_int_tensor%prop_ov')
        ! call print_tensor(orb_int_tensor%prop_vo, 1.0D-6, 'orb_int_tensor%prop_vo')

end subroutine 


!------------------------------------------------------------------------------------------------------------------
subroutine get_prop_integral_driver(orb_int_tensor,nmo1,mo1_list,nmo2,mo2_list,min_occ,i_prop,&
                                    exchange,eti_precision,prop_label)
    ! Read proporty integrals from  MDPROP file 
    type(tensor_t), intent(inout)        :: orb_int_tensor
    type(eti_precision_t), intent(in)    :: eti_precision
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

    ierr=construct_tensor_t(orb_int_tensor,eti_precision%storage_t,mo_dims,init_val=ZERO)
    print *,'ierr in get_prop_integral_driver',ierr

    ierr=direct_data_pointer(orb_int_tensor,body_p)
    call c_f_pointer(body_p,orbital_integral_tens,mo_dims)

    print *,'debug orb_int_tensor%datatype',orb_int_tensor%datatype
    print *,'debug eti_precision%storage_t', eti_precision%storage_t

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
    print *,'debug exited get_prop_integral_driver'

end subroutine get_prop_integral_driver
!-------------------------------------------------------------------------------------------------------------------

subroutine construct_aux_tensor(eti_precision,LR_Aux_ov,LR_Aux_vv,LR_Aux_oo,LR_Aux_vo, ov_dims, vo_dims, oo_dims, vv_dims)

        integer(INTD), dimension(2), intent(in)  :: ov_dims, oo_dims, vv_dims, vo_dims
        type(tensor_t), intent(inout)            :: LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo
        integer(C_INT) :: ierr
        type(eti_precision_t), intent(in) :: eti_precision

        ierr=construct_tensor_t(LR_Aux_ov,eti_precision%storage_t,ov_dims,init_val=ZERO)
        ierr=construct_tensor_t(LR_Aux_vo,eti_precision%storage_t,vo_dims,init_val=ZERO)
        ierr=construct_tensor_t(LR_Aux_vv,eti_precision%storage_t,vv_dims,init_val=ZERO)
        ierr=construct_tensor_t(LR_Aux_oo,eti_precision%storage_t,oo_dims,init_val=ZERO)

end subroutine

!--------------------------------------------------------------------------------------------------------------------

subroutine destroy_aux_tensor(LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo)
    type(tensor_t), intent(inout)      :: LR_Aux_ov,LR_Aux_vo,LR_Aux_vv,LR_Aux_oo
    integer(C_INT) :: ierr

    ierr=destroy_tensor_t(LR_Aux_oo)
    ierr=destroy_tensor_t(LR_Aux_vv)
    ierr=destroy_tensor_t(LR_Aux_ov)
    ierr=destroy_tensor_t(LR_Aux_vo)

end subroutine 

!------------------------------------------------------------------------------------------------------------------

subroutine destroy_Tprop_tensor(Tprop_tensor)
    type(eti_comm_tens_t), intent(inout) :: Tprop_tensor
    
    integer(C_INT) :: ierr
    
    ierr=destroy_tensor_t(Tprop_tensor%t1)
    ierr=destroy_tensor_t(Tprop_tensor%t2)

end subroutine

!------------------------------------------------------------------------------------------------------------------

subroutine destroy_Tprop_tensors(Tprop_tensors_wB,    Tprop_tensors_min_wB, &
                                i, &
                                Bvec_tensors, Bvec_tensors_min_w, &
                                Tprop_bar_tensors_wB, Tprop_bar_tensors_min_wB, &
                                Tprop_tensors_wC,     Tprop_tensors_min_wC, &
                                Tprop_bar_tensors_wC, Tprop_bar_tensors_min_wC   )
    type(eti_comm_tens_t), intent(inout), allocatable :: Tprop_tensors_wB(:)     !LR
    type(eti_comm_tens_t), intent(inout), allocatable :: Tprop_tensors_min_wB(:) !LR
    integer        :: i
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Bvec_tensors(:)
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Bvec_tensors_min_w(:)
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_wB(:)
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_min_wB(:)
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_tensors_wC(:)  !QR
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_tensors_min_wC(:) !QR
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_wC(:) !QR
    type(eti_comm_tens_t), intent(inout), allocatable, OPTIONAL :: Tprop_bar_tensors_min_wC(:) !QR
    
    integer(C_INT) :: ierr
    
    ierr=destroy_tensor_t(Tprop_tensors_wB(i)%t1)
    ierr=destroy_tensor_t(Tprop_tensors_wB(i)%t2)
    if (present(Tprop_bar_tensors_wB)) then 
        ierr=destroy_tensor_t(Tprop_bar_tensors_wB(i)%t1)
        ierr=destroy_tensor_t(Tprop_bar_tensors_wB(i)%t2)
    end if

    ierr=destroy_tensor_t(Tprop_tensors_min_wB(i)%t1)
    ierr=destroy_tensor_t(Tprop_tensors_min_wB(i)%t2)
    if (present(Tprop_bar_tensors_min_wB)) then
        ierr=destroy_tensor_t(Tprop_bar_tensors_min_wB(i)%t1)
        ierr=destroy_tensor_t(Tprop_bar_tensors_min_wB(i)%t2)
    end if

    if (present(Bvec_tensors)) then
        ierr=destroy_tensor_t(Bvec_tensors(i)%t1)
        ierr=destroy_tensor_t(Bvec_tensors(i)%t2)
    end if

    if (present(Bvec_tensors_min_w)) then
        ierr=destroy_tensor_t(Bvec_tensors_min_w(i)%t1)
        ierr=destroy_tensor_t(Bvec_tensors_min_w(i)%t2)
    end if

    if (present(Tprop_tensors_wC)) then
        ierr=destroy_tensor_t(Tprop_tensors_wC(i)%t1)
        ierr=destroy_tensor_t(Tprop_tensors_wC(i)%t2)
    end if
    if (present(Tprop_bar_tensors_wC)) then 
        ierr=destroy_tensor_t(Tprop_bar_tensors_wC(i)%t1)
        ierr=destroy_tensor_t(Tprop_bar_tensors_wC(i)%t2)
    end if

    if (present(Tprop_tensors_min_wC)) then
        ierr=destroy_tensor_t(Tprop_tensors_min_wC(i)%t1)
        ierr=destroy_tensor_t(Tprop_tensors_min_wC(i)%t2)
    end if
    if (present(Tprop_bar_tensors_min_wC)) then
        ierr=destroy_tensor_t(Tprop_bar_tensors_min_wC(i)%t1)
        ierr=destroy_tensor_t(Tprop_bar_tensors_min_wC(i)%t2)
    end if


end subroutine 


!---------------------------------------------------------------------------------------------------------------

subroutine destroy_Interm_t_tensors(interm_t_prop1, CC_CC)
    type(eti_rsp_interm_tens_t), intent(inout) :: interm_t_prop1
    logical          :: CC_CC
    integer          :: ierr

    ierr=destroy_tensor_t(interm_t_prop1%epsilon_prop_vo)
    ierr=destroy_tensor_t(interm_t_prop1%epsilon_prop_vvoo)
    ierr=destroy_tensor_t(interm_t_prop1%eta_prop_oovv)
    ierr=destroy_tensor_t(interm_t_prop1%eta_prop_ov)

    if(.not.CC_CC) then
        ierr=destroy_tensor_t(interm_t_prop1%eta_prop_ov_eom)
    end if

    ierr=destroy_tensor_t(interm_t_prop1%prop_bar_oo)
    ierr=destroy_tensor_t(interm_t_prop1%prop_bar_ov)    
    ierr=destroy_tensor_t(interm_t_prop1%prop_bar_vv)
    ierr=destroy_tensor_t(interm_t_prop1%O_vo)

end subroutine



end module
