module talsh_eom_QR
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

    public get_eom_qr

    contains

!-----------------------------------------------------------------------------------------------------------------


subroutine get_eom_qr(exa_input, interm_t, Interm_t_tensors,& 
                        Tprop_tensors, Tprop_bar_tensors,& 
                        Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                        i, j, k, l1_tensor, l2_tensor, eom_qr)
        use talsh_sigma_eomccsd

        type(exacc_input), intent(in) :: exa_input
        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_tens_t), intent(inout)           :: l2_tensor
        type(talsh_tens_t), intent(inout)           :: l1_tensor
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors(:)


        type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
        type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR

        integer                       :: i, j, k
        
        complex(kind=8)               :: quarsp, quarsp_con

        integer                       :: ierr

        complex(kind=8)               :: eom_qr


        call print_date('Evaluate CC-CI QR')



        call get_quarspf(exa_input, interm_t, Interm_t_tensors,&
                        Tprop_tensors,Tprop_bar_tensors,& 
                        i, j, k, l1_tensor, l2_tensor, quarsp)



        call get_quarspf(exa_input, interm_t, Interm_t_tensors,&
                        Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                        i, j, k, l1_tensor, l2_tensor, quarsp_con)

        quarsp_con=conjg(quarsp_con)

        eom_qr=(quarsp+quarsp_con)*0.5

end subroutine get_eom_qr


!------------------------------------------------------------------------------------------------------------------

subroutine get_quarspf(exa_input, interm_t, Interm_t_tensors, &
                            Tprop_tensors,Tprop_bar_tensors,&
                            i, j, k, l1_tensor, l2_tensor, quarsp)         
            type(exacc_input), intent(in) :: exa_input
            type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
            type(talsh_tens_t), intent(inout)           :: l2_tensor
            type(talsh_tens_t), intent(inout)           :: l1_tensor
            type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors(:)
            type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:),Tprop_bar_tensors(:)
            integer                            :: i, j, k

                
            type(talsh_comm_tens_t)       :: sigma_prop_tensor12, sigma_prop_tensor13, sigma_prop_tensor23
            type(talsh_comm_tens_t)       :: sigma_prop_tensor21, sigma_prop_tensor31, sigma_prop_tensor32
            complex(kind=8)               :: quarsp
            complex(kind=8)               :: tem_quarsp_iter


            !XYZ   X is 1; Y is 2; Z is 3
            call construct_AyTz(exa_input,Interm_t_tensors(i),Tprop_tensors(3),sigma_prop_tensor13)
            ! print *, "error1?"

            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(2), Tprop_tensors(3), &
                                Tprop_bar_tensors(1), Tprop_bar_tensors(2), Tprop_bar_tensors(3), &
                                Interm_t_tensors(i), Interm_t_tensors(k), &
                                sigma_prop_tensor13,l1_tensor,l2_tensor)
            ! print *, "error2?"
                                
            quarsp=quarsp+tem_quarsp_iter



            !XZY
            call construct_AyTz(exa_input,Interm_t_tensors(i),Tprop_tensors(2),sigma_prop_tensor12)
            ! print *, "error3?"

            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(3), Tprop_tensors(2), &
                                Tprop_bar_tensors(1), Tprop_bar_tensors(3), Tprop_bar_tensors(2), &
                                Interm_t_tensors(i), Interm_t_tensors(j), &
                                sigma_prop_tensor12,l1_tensor,l2_tensor)
                
            ! print *, "error4?"
            quarsp=quarsp+tem_quarsp_iter


            !YXZ
            call construct_AyTz(exa_input,Interm_t_tensors(j),Tprop_tensors(3),sigma_prop_tensor23)
            ! print *, "error5?"

            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(1), Tprop_tensors(3), &
                                Tprop_bar_tensors(2), Tprop_bar_tensors(1), Tprop_bar_tensors(3), &
                                Interm_t_tensors(j), Interm_t_tensors(k), &
                                sigma_prop_tensor23,l1_tensor,l2_tensor)


            ! print *, "error6?"
            quarsp=quarsp+tem_quarsp_iter


            !YZX
            call construct_AyTz(exa_input,Interm_t_tensors(j),Tprop_tensors(1),sigma_prop_tensor21)
            ! print *, "error7?"

            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(3), Tprop_tensors(1), &
                                Tprop_bar_tensors(2), Tprop_bar_tensors(3), Tprop_bar_tensors(1), &
                                Interm_t_tensors(j), Interm_t_tensors(i), &
                                sigma_prop_tensor21,l1_tensor,l2_tensor)

            ! print *, "error8?"
            quarsp=quarsp+tem_quarsp_iter


            !ZYX
            call construct_AyTz(exa_input,Interm_t_tensors(k),Tprop_tensors(1),sigma_prop_tensor31)
            ! print *, "error9?"

            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(2), Tprop_tensors(1), &
                                Tprop_bar_tensors(3), Tprop_bar_tensors(2), Tprop_bar_tensors(1), &
                                Interm_t_tensors(k), Interm_t_tensors(i), &
                                sigma_prop_tensor31,l1_tensor,l2_tensor)
            ! print *, "error10?"
            quarsp=quarsp+tem_quarsp_iter


            !ZXY
            call construct_AyTz(exa_input,Interm_t_tensors(k),Tprop_tensors(2),sigma_prop_tensor32)
            ! print *, "error11?"

            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(1), Tprop_tensors(2), &
                                Tprop_bar_tensors(3), Tprop_bar_tensors(1), Tprop_bar_tensors(2), &
                                Interm_t_tensors(k), Interm_t_tensors(j), &
                                sigma_prop_tensor32,l1_tensor,l2_tensor)
            
            ! print *, "error12?"
            quarsp=quarsp+tem_quarsp_iter

end subroutine get_quarspf

!------------------------------------------------------------------------------------------------------------------

subroutine get_quarsp_driver(tem_quarsp_iter,tprop_tensor1, tprop_tensor2, &
                            tprop_bar_tensor1, tprop_bar_tensor2, tprop_bar_tensor3, &
                            interm_t_prop1, interm_t_prop2, &
                            sigma_prop_tensor,l1_tensor,l2_tensor)

            type(talsh_rsp_interm_tens_t)    :: interm_t_prop1, interm_t_prop2
            type(talsh_tens_t)               :: l2_tensor
            type(talsh_tens_t)               :: l1_tensor
            type(talsh_comm_tens_t)          :: tprop_tensor1, tprop_tensor2
            type(talsh_comm_tens_t)          :: tprop_bar_tensor1, tprop_bar_tensor2, tprop_bar_tensor3
            type(talsh_comm_tens_t)          :: sigma_prop_tensor
            complex(kind=8)                  :: tem_quarsp_iter

!  Point connected with tensors 
            type(talsh_tens_t)          :: result_tensor
            integer(C_INT)              :: ierr
            complex(8), pointer         :: result_tens(:)
            integer                     :: initialize_talsh

            type(C_PTR):: body_p
            integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
            integer(C_INT):: host_arg_max

            integer(C_INT)              :: one_dims(1)
            integer(INTD), dimension(2)   :: ov_dims
            integer(INTD)  ::  ov_rank
            type(talsh_tens_t) :: Aux_ov

            complex(8)                  :: term1, term2, term3, term4, term5, term6


            one_dims(1) = 1
            ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
            ! print *, "ierr1 quarsp driver: ", ierr
            ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,result_tens,one_dims)
            ! print *, "ierr2 quarsp driver: ", ierr


            ov_rank=talsh_tensor_rank(tprop_bar_tensor1%t1)
            ierr = talsh_tensor_dimensions(tprop_bar_tensor1%t1,ov_rank,ov_dims)
            ierr=talsh_tensor_construct(Aux_ov,C8,ov_dims,init_val=ZERO)

        ! Term1: -Tbar_X*Ty    
            ierr=talsh_tensor_contract("E()+=Tlx(i,a)*Ty(a,i)",& 
                                    result_tensor,tprop_bar_tensor1%t1,tprop_tensor1%t1,scale=MINUS_ONE)
                                    ! print *, "ierr3 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Tlx(i,j,a,b)*Ty(a,b,i,j)", &
                                    result_tensor,tprop_bar_tensor1%t2,tprop_tensor1%t2,scale=MINUS_ONE_QUARTER_C)
                                    ! print *, "ierr4 quarsp driver: ", ierr
            term1=result_tens(1)
            print *, "term1: ", term1


        ! Term2: Lambda*Epison_Z    
            ierr=talsh_tensor_init(result_tensor)
            ! print *, "ierr5 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Epz(a,i)*L(i,a)",result_tensor,interm_t_prop2%epsilon_prop_vo,l1_tensor)
            ! print *, "ierr6 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Epz(a,b,i,j)*L(i,j,a,b)",result_tensor,&
                                        interm_t_prop2%epsilon_prop_vvoo,l2_tensor,scale=ONE_QUARTER_C)
                                        ! print *, "ierr7 quarsp driver: ", ierr
            term2=result_tens(1)
            print *, "term2: ", term2


        ! Term3: Tbar_Y*Sigma_XZ 
            ierr=talsh_tensor_init(result_tensor)
            ! print *, "ierr8 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Tly(i,a)*Sixz(a,i)",result_tensor,tprop_bar_tensor2%t1,sigma_prop_tensor%t1)
            ! print *, "ierr9 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Tly(i,j,a,b)*Sigxz(a,b,i,j)",result_tensor,&
                                        tprop_bar_tensor2%t2,sigma_prop_tensor%t2,scale=ONE_QUARTER_C)
                                        ! print *, "ierr10 quarsp driver: ", ierr
            term3=result_tens(1)
            print *, "term3: ", term3


        ! Term4: -Tbar_Y(D)*Epison_X(S)*TZ(S)
            ierr=talsh_tensor_init(result_tensor)
            ! print *, "ierr11 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("Aux(i,a)+=Tly(i,j,a,b)*EpX(b,j)",Aux_ov,&
                                        tprop_bar_tensor2%t2,interm_t_prop1%epsilon_prop_vo,scale=MINUS_ONE)

            ! print *, "ierr12 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Aux(i,a)*TZ(a,i)",result_tensor,Aux_ov,tprop_tensor2%t1,scale=ONE)
            ! print *, "ierr13 quarsp driver: ", ierr
            term4=result_tens(1)
            print *, "term4: ", term4


        ! Term5: -Lamda*TY
            ierr=talsh_tensor_init(result_tensor)
            ! print *, "ierr14 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=L(i,a)*Tl(a,i)",result_tensor,l1_tensor,tprop_tensor1%t1,scale=MINUS_ONE)
            ! print *, "ierr15 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=L(i,j,a,b)*Tl(a,b,i,j)",result_tensor,&
                                        l2_tensor,tprop_tensor1%t2,scale=MINUS_ONE_QUARTER_C)
                                        ! print *, "ierr16 quarsp driver: ", ierr
            term5=result_tens(1)
            print *, "term5: ", term5


        ! Term6: Tbar_Z*Epison_X
            ierr=talsh_tensor_init(result_tensor)
            ! print *, "ierr17 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Tl(i,a)*Epz(a,i)",result_tensor,&
                                        tprop_bar_tensor3%t1,interm_t_prop1%epsilon_prop_vo,scale=ONE)
                                        ! print *, "ierr18 quarsp driver: ", ierr
            ierr=talsh_tensor_contract("E()+=Tl(i,j,a,b)*Epz(a,b,i,j)",result_tensor,&
                                        tprop_bar_tensor3%t2,interm_t_prop1%epsilon_prop_vvoo,scale=ONE_QUARTER_C)
                                        ! print *, "ierr19 quarsp driver: ", ierr
            term6=result_tens(1)
            print *, "term6: ", term6

            
            tem_quarsp_iter=term1*term2+term3+term4+term5*term6

            ierr=talsh_tensor_destruct(Aux_ov)
            ! print *, "ierr20 quarsp driver: ", ierr

end subroutine 

!------------------------------------------------------------------------------------------------------------------

subroutine construct_AyTz(exa_input,Interm_t_tensor ,tprop_tensor ,sigma_prop_tensor)
    type(exacc_input), intent(in) :: exa_input
    type(talsh_rsp_interm_tens_t) :: Interm_t_tensor
    type(talsh_comm_tens_t)       :: tprop_tensor     
    type(talsh_comm_tens_t)       :: sigma_prop_tensor

    integer(INTD), dimension(4)   :: vvoo_dims
    integer(INTD), dimension(2)   :: vo_dims
    integer                       :: ierr


    vvoo_dims(1:2) = exa_input%nvir
    vvoo_dims(3:4) = exa_input%nocc
    vo_dims(1) = exa_input%nvir
    vo_dims(2) = exa_input%nocc



    ierr=talsh_tensor_construct(sigma_prop_tensor%t1,C8,vo_dims,init_val=ZERO)
    ! print *, "ierr1: ", ierr
    ierr=talsh_tensor_construct(sigma_prop_tensor%t2,C8,vvoo_dims,init_val=ZERO)
    ! print *, "ierr2: ", ierr
    call get_sigma_eomccsd_prop(Interm_t_tensor, tprop_tensor, sigma_prop_tensor)    

end subroutine construct_AyTz

!------------------------------------------------------------------------------------------------------------------

subroutine get_sigma_eomccsd_prop(interm_t_prop1, tprop_tensor2, sigma_prop_tensor12)

    type(talsh_rsp_interm_tens_t) :: interm_t_prop1
    type(talsh_comm_tens_t)            :: tprop_tensor2
    type(talsh_comm_tens_t)       :: sigma_prop_tensor12

    integer :: ierr


!***************
!l1_sigma_prop *
!***************

!--------------------------------------------------------
!     term 1: Sxz(a,i) += Xbar(a,e)* tY(e,i) |
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,i)+=Xbar(a,e)*Tz(e,i)",sigma_prop_tensor12%t1,interm_t_prop1%prop_bar_vv,tprop_tensor2%t1)
    ! print *, "ierr3: ", ierr

!--------------------------------------------------------
!     term 2: Sxz(a,i) += -Xbar(m,i)* tY(a,m) |
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,i)+=Xbar(m,i)*Tz(a,m)", &
                sigma_prop_tensor12%t1,interm_t_prop1%prop_bar_oo,tprop_tensor2%t1,scale=MINUS_ONE)
                ! print *, "ierr4: ", ierr

!--------------------------------------------------------
!     term 3: Sxz(a,i) += Xbar(m,e)* tY(e,a,m,i) |
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,i)+=Xbar(m,e)*Tz(e,a,m,i)", &
                sigma_prop_tensor12%t1,interm_t_prop1%prop_bar_ov,tprop_tensor2%t2)
                ! print *, "ierr5: ", ierr


!***************
!l2_sigma_prop *
!***************

!--------------------------------------------------------
!     term 1: Sxz(a,b,i,j) += -P(ab) Q(m,b,i,j)* tY(a,m) |
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Q(m,b,i,j)*Tz(a,m)", &
                sigma_prop_tensor12%t2,interm_t_prop1%Q_ovoo,tprop_tensor2%t1,scale=MINUS_ONE)
                ! print *, "ierr6: ", ierr
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Q(m,a,i,j)*Tz(b,m)", &
                sigma_prop_tensor12%t2,interm_t_prop1%Q_ovoo,tprop_tensor2%t1,scale=ONE)
                ! print *, "ierr7: ", ierr
    
!--------------------------------------------------------
!     term 2: Sxz(a,b,i,j) += P(ij) Q(a,b,j,e)* tY(e,i)   |
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Q(a,b,e,j)*Tz(e,i)", &
                sigma_prop_tensor12%t2,interm_t_prop1%Q_vvvo,tprop_tensor2%t1,scale=ONE)
                ! print *, "ierr8: ", ierr
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Q(a,b,e,i)*Tz(e,j)", &
                sigma_prop_tensor12%t2,interm_t_prop1%Q_vvvo,tprop_tensor2%t1,scale=MINUS_ONE)
                ! print *, "ierr9: ", ierr

!--------------------------------------------------------
!     term 3: Sxz(a,b,i,j) += P(ab) Xbar(b,e)* tY(a,e,i,j)|
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Xbar(b,e)*Tz(a,e,i,j)", &
                sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_vv,tprop_tensor2%t2,scale=ONE)
                ! print *, "ierr10: ", ierr
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Xbar(a,e)*Tz(b,e,i,j)", &
                sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_vv,tprop_tensor2%t2,scale=MINUS_ONE)
                ! print *, "ierr11: ", ierr

!--------------------------------------------------------
!     term 4: Sxz(a,b,i,j) += -P(ij) Xbar(m,j)* tY(a,b,i,m)|
!--------------------------------------------------------
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Xbar(m,j)*Tz(a,b,i,m)", &
                sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_oo,tprop_tensor2%t2,scale=MINUS_ONE)
                ! print *, "ierr12: ", ierr
    ierr=talsh_tensor_contract("Sxz(a,b,i,j)+=Xbar(m,i)*Tz(a,b,j,m)", &
                sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_oo,tprop_tensor2%t2,scale=ONE)
                ! print *, "ierr13: ", ierr


  end subroutine get_sigma_eomccsd_prop


end module
