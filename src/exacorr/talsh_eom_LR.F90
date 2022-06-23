module talsh_eom_LR
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global
    use talsh_eom_rspeq_solver
    ! use talsh_cc_ci

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public get_eom_lr

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine get_eom_lr(exa_input,interm_t,interm_t_prop1, interm_t_prop2, interm_t_prop3, &
                tprop_tensor1, tprop_tensor2, &
                l1_tensor,l2_tensor, eom_lr)
!This subroutine drives the CC-CI linear response property calculation
    !   Written by Xiang Yuan Janunary 2022, Lille, France

            use talsh_sigma_eomccsd

            type(exacc_input), intent(in) :: exa_input
            type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
            type(talsh_tens_t), intent(inout)           :: l2_tensor
            type(talsh_tens_t), intent(inout)           :: l1_tensor
            type(talsh_rsp_interm_tens_t) :: interm_t_prop1, interm_t_prop2, interm_t_prop3  
            type(talsh_comm_tens_t)            :: tprop_tensor1, tprop_tensor2


            integer(INTD), dimension(4)   :: oovv_dims
            integer(INTD), dimension(2)   :: ov_dims
            type(talsh_comm_tens_t)       :: sigma_prop_tensor12, sigma_prop_tensor21
            complex(kind=8)               :: quarsp, quarsp_con

            integer                       :: ierr

            complex(kind=8)               :: eom_lr


            call print_date('Evaluate CC-CI LR')

            oovv_dims(1:2) = exa_input%nocc
            oovv_dims(3:4) = exa_input%nvir
            ov_dims(1) = exa_input%nocc
            ov_dims(2) = exa_input%nvir


            ierr=talsh_tensor_construct(sigma_prop_tensor12%t1,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(sigma_prop_tensor12%t2,C8,oovv_dims,init_val=ZERO)

            ierr=talsh_tensor_construct(sigma_prop_tensor21%t1,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(sigma_prop_tensor21%t2,C8,oovv_dims,init_val=ZERO)


            ! Obtain eta_(nu)^(XY)(wY)
            call sigma_eomccsd_prop(interm_t_prop1, tprop_tensor2, sigma_prop_tensor12)     ! talsh_sigma_eomccsd.F90
            call sigma_eomccsd_prop(interm_t_prop2, tprop_tensor1, sigma_prop_tensor21)


            call get_linrspf(interm_t, interm_t_prop1, interm_t_prop2, interm_t_prop3,&
                    tprop_tensor1, tprop_tensor2, &
                    sigma_prop_tensor12, sigma_prop_tensor21,  &
                    l1_tensor, l2_tensor, quarsp)

            quarsp_con=conjg(quarsp)
            eom_lr=(quarsp+quarsp_con)*0.5

end subroutine get_eom_lr

!------------------------------------------------------------------------------------------------------------------

subroutine get_linrspf(interm_t, interm_t_prop1, interm_t_prop2, interm_t_prop3, &
                tprop_tensor1, tprop_tensor2, &
                sigma_prop_tensor12, sigma_prop_tensor21, &
                l1_tensor, l2_tensor, quarsp)

        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_tens_t), intent(inout)           :: l2_tensor
        type(talsh_tens_t), intent(inout)           :: l1_tensor
        type(talsh_rsp_interm_tens_t) :: interm_t_prop1, interm_t_prop2, interm_t_prop3  
        type(talsh_comm_tens_t)            :: tprop_tensor1, tprop_tensor2
        type(talsh_comm_tens_t)       :: sigma_prop_tensor12, sigma_prop_tensor21
        complex(kind=8)               :: quarsp
        complex(kind=8)               :: tem_quarsp_iter

!--------------------------------------------------------------------------------------
!       << X;Y >>ωY = P(X,Y) [term1 * term2 + term3 + term4]                          |
!--------------------------------------------------------------------------------------
            call get_linrsp_driver(tem_quarsp_iter,tprop_tensor1, tprop_tensor2, &
                    sigma_prop_tensor12,l1_tensor,l2_tensor, &
                    interm_t,interm_t_prop1)
            quarsp=quarsp+tem_quarsp_iter

           tem_quarsp_iter=complex(0.0,0.0)

            call get_linrsp_driver(tem_quarsp_iter,tprop_tensor2, tprop_tensor1, &
                    sigma_prop_tensor21,l1_tensor,l2_tensor, &
                    interm_t,interm_t_prop1)

            quarsp=quarsp-tem_quarsp_iter


end subroutine get_linrspf

!------------------------------------------------------------------------------------------------------------------

subroutine get_linrsp_driver(tem_quarsp_iter,tprop_tensor1, tprop_tensor2, &
                sigma_prop_tensor12,l1_tensor,l2_tensor,&
                interm_t,interm_t_prop1)

    type(talsh_intermediates_tens_t),intent(inout) :: interm_t 
    type(talsh_rsp_interm_tens_t),intent(inout)    :: interm_t_prop1
    type(talsh_tens_t), intent(inout)              :: l2_tensor
    type(talsh_tens_t), intent(inout)              :: l1_tensor
    complex(kind=8)               :: tem_quarsp_iter

    type(talsh_comm_tens_t)            :: tprop_tensor1, tprop_tensor2
        
    type(talsh_comm_tens_t)       :: sigma_prop_tensor12



    !  Point connected with tensors 
    type(talsh_tens_t)          :: result_tensor
    integer(C_INT)              :: ierr
    complex(8), pointer         :: result_tens(:)
    integer                     :: initialize_talsh

    type(C_PTR):: body_p
    integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
    integer(C_INT):: host_arg_max

    integer(C_INT)              :: one_dims(1)

    complex(8)                  :: term1, term2, term3, term4


    one_dims(1) = 1
    ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,result_tens,one_dims)

 
!--------------------------------------------------------------------------------------
!       Term1 -= t_bar0(nu)*tX(nu)   XY: do we need to transpose Tlx(i,a)?            |
!--------------------------------------------------------------------------------------
    ierr=talsh_tensor_contract("E()+=Tl0(a,i)*Tx(a,i)",& 
            result_tensor,l1_tensor,tprop_tensor1%t1,scale=MINUS_ONE)
    term1=result_tens(1)
    ierr=talsh_tensor_init(result_tensor)
    ierr=talsh_tensor_contract("E()+=Tl0(a,b,i,j)*Tx(a,b,i,j)", &
            result_tensor,l2_tensor,tprop_tensor1%t2,scale=MINUS_ONE)
    term1=term1+result_tens(1)


!--------------------------------------------------------------------------------------
!       Term2 += tY(mu)*eta(mu)                                                       |
!--------------------------------------------------------------------------------------
    ierr=talsh_tensor_init(result_tensor)
    ierr=talsh_tensor_contract("E()+=Ty(a,i)*n(i,a)",result_tensor,tprop_tensor2%t1,interm_t%eta_ov)
    term2=result_tens(1)
    ierr=talsh_tensor_init(result_tensor)
    ierr=talsh_tensor_contract("E()+=Ty(a,b,i,j)*n(i,j,a,b)",result_tensor,tprop_tensor2%t2,interm_t%eta_oovv)
    term2=term2+result_tens(1)


!--------------------------------------------------------------------------------------
!       Term3 += tbar0(nu)*SigmaXY(nu)                                                |
!--------------------------------------------------------------------------------------
    ierr=talsh_tensor_init(result_tensor)
    ierr=talsh_tensor_contract("E()+=Tl0(a,i)*Sigxy(a,i)",result_tensor,l1_tensor,sigma_prop_tensor12%t1)
    term3=result_tens(1)
    ierr=talsh_tensor_init(result_tensor)
    ierr=talsh_tensor_contract("E()+=Tl0(a,b,i,j)*Sigxy(a,b,i,j)",result_tensor,l2_tensor,sigma_prop_tensor12%t2)
    term3=term3+result_tens(1)


!--------------------------------------------------------------------------------------
!       Term4 += tY(mu)*OX(mu)                                                       |
!--------------------------------------------------------------------------------------
    ierr=talsh_tensor_init(result_tensor)
    ierr=talsh_tensor_contract("E()+=Ty(a,i)*O(a,i)",result_tensor,tprop_tensor2%t1,interm_t_prop1%O_vo)
    term4=result_tens(1)
    ! ierr=talsh_tensor_init(result_tensor)
    ! ierr=talsh_tensor_contract("E()+=Lambda(a,b,i,j)*Ty(i,j,a,b)",result_tensor,tprop_tensor3%t1,interm_t%O_oovv)
    ! term5=term5+result_tens(1)

    tem_quarsp_iter=term1*term2+term3+term4

end subroutine 

end module
