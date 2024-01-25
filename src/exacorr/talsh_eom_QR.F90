module talsh_eom_QR
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global
    use talsh_eom_rspeq_solver

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


subroutine get_eom_qr(exa_input, interm_t,&
                        Interm_t_tensors_A,Interm_t_tensors_B, Interm_t_tensors_C,& 
                        Tprop_tensors, Tprop_bar_tensors,& 
                        Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                        i, j, k, l1_tensor, l2_tensor, eom_qr)
   
!   This subroutine drives the EOM-CCSD Quadratic response property calculation
!   Written by Xiang Yuan Janunary 2022, Lille, France

        use talsh_sigma_eomccsd
        type(exacc_input), intent(in) :: exa_input
        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_tens_t), intent(inout)           :: l2_tensor
        type(talsh_tens_t), intent(inout)           :: l1_tensor
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:)
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_B(:)
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_C(:)


        type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
        type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR

        integer                       :: i, j, k
        integer                       :: T
        
        complex(kind=8)               :: quarsp, quarsp_con

        integer                       :: ierr

        complex(kind=8)               :: eom_qr
        real(kind=8)                  :: tem_real, tem_imag

        type(talsh_tens_t)          :: result_tensor
        integer(C_INT)              :: one_dims(1)
        complex(8), pointer         :: result_tens(:)
        type(C_PTR):: body_p


        call print_date('Evaluate CC-CI QR')


        call get_quarspf(exa_input, interm_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,&
                        Tprop_tensors,Tprop_bar_tensors,& 
                        i, j, k, l1_tensor, l2_tensor, quarsp)

        call get_quarspf(exa_input, interm_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,&
                        Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                        i, j, k, l1_tensor, l2_tensor, quarsp_con)

        print *, "quadratic response: +w", quarsp
        print *, "quadratic response: -w", quarsp_con

        T = (-1)**(exa_input%N_Time_anti)

        if (T==-1) then
            print *, "We now flip the real and imaginary parts"
            tem_real = real(quarsp)
            tem_imag = aimag(quarsp)
            quarsp = cmplx(tem_imag,tem_real)

            tem_real = real(quarsp_con)
            tem_imag = aimag(quarsp_con)
            quarsp_con = cmplx(tem_imag,tem_real)

            print *, "New quadratic response: +w", quarsp
            print *, "New quadratic response: -w", quarsp_con
        end if

        quarsp_con=conjg(quarsp_con)
        eom_qr=(quarsp+quarsp_con)*0.5

end subroutine get_eom_qr


!------------------------------------------------------------------------------------------------------------------

subroutine get_quarspf(exa_input, interm_t,&
                        Interm_t_tensors_A, Interm_t_tensors_B, Interm_t_tensors_C,&
                        Tprop_tensors,Tprop_bar_tensors,&
                        i, j, k, l1_tensor, l2_tensor, quarsp)         
            type(exacc_input), intent(in) :: exa_input
            type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
            type(talsh_tens_t), intent(inout)           :: l2_tensor
            type(talsh_tens_t), intent(inout)           :: l1_tensor
            type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:)
            type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_B(:)
            type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_C(:)
            type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:),Tprop_bar_tensors(:)
            integer                            :: i, j, k

                
            type(talsh_comm_tens_t)       :: sigma_prop_tensor12, sigma_prop_tensor13, sigma_prop_tensor23
            type(talsh_comm_tens_t)       :: sigma_prop_tensor21, sigma_prop_tensor31, sigma_prop_tensor32
            complex(kind=8)               :: quarsp
            complex(kind=8)               :: tem_quarsp_iter
            type(talsh_comm_tens_t)       :: bvec_tensor


            !XYZ   X is 1; Y is 2; Z is 3  Ax*Tly

            call construct_TlXAy(exa_input,Interm_t_tensors_A(i),interm_t,Tprop_bar_tensors(2),sigma_prop_tensor12)
            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(2), Tprop_tensors(3), &
                                Tprop_bar_tensors(1), Tprop_bar_tensors(2), Tprop_bar_tensors(3), &
                                Interm_t_tensors_A(i), Interm_t_tensors_C(k), &
                                sigma_prop_tensor12,l1_tensor,l2_tensor)
                                
            quarsp=quarsp+tem_quarsp_iter



            !XZY   Ax*Tlz

            call construct_TlXAy(exa_input,Interm_t_tensors_A(i),interm_t,Tprop_bar_tensors(3),sigma_prop_tensor13)
            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(3), Tprop_tensors(2), &
                                Tprop_bar_tensors(1), Tprop_bar_tensors(3), Tprop_bar_tensors(2), &
                                Interm_t_tensors_A(i), Interm_t_tensors_B(j), &
                                sigma_prop_tensor13,l1_tensor,l2_tensor)

            quarsp=quarsp+tem_quarsp_iter


            !YXZ   Ay*Tlx

            call construct_TlXAy(exa_input,Interm_t_tensors_B(j),interm_t,Tprop_bar_tensors(1),sigma_prop_tensor21)
            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(1), Tprop_tensors(3), &
                                Tprop_bar_tensors(2), Tprop_bar_tensors(1), Tprop_bar_tensors(3), &
                                Interm_t_tensors_B(j), Interm_t_tensors_C(k), &
                                sigma_prop_tensor21,l1_tensor,l2_tensor)

            quarsp=quarsp+tem_quarsp_iter


            !YZX  Ay*Tlz

            call construct_TlXAy(exa_input,Interm_t_tensors_B(j),interm_t,Tprop_bar_tensors(3),sigma_prop_tensor23)
            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(3), Tprop_tensors(1), &
                                Tprop_bar_tensors(2), Tprop_bar_tensors(3), Tprop_bar_tensors(1), &
                                Interm_t_tensors_B(j), Interm_t_tensors_A(i), &
                                sigma_prop_tensor23,l1_tensor,l2_tensor)
            quarsp=quarsp+tem_quarsp_iter


            !ZYX Az*Tly
       
            call construct_TlXAy(exa_input,Interm_t_tensors_C(k),interm_t,Tprop_bar_tensors(2),sigma_prop_tensor32)
            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(2), Tprop_tensors(1), &
                                Tprop_bar_tensors(3), Tprop_bar_tensors(2), Tprop_bar_tensors(1), &
                                Interm_t_tensors_C(k), Interm_t_tensors_A(i), &
                                sigma_prop_tensor32,l1_tensor,l2_tensor)
            quarsp=quarsp+tem_quarsp_iter


            !ZXY Az*Tlx
            
            call construct_TlXAy(exa_input,Interm_t_tensors_C(k),interm_t,Tprop_bar_tensors(1),sigma_prop_tensor31)
            call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(1), Tprop_tensors(2), &
                                Tprop_bar_tensors(3), Tprop_bar_tensors(1), Tprop_bar_tensors(2), &
                                Interm_t_tensors_C(k), Interm_t_tensors_B(j), &
                                sigma_prop_tensor31,l1_tensor,l2_tensor)
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

            type(C_PTR):: body_p, body_p2, body_p3
            integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
            integer(C_INT):: host_arg_max

            integer(C_INT)              :: one_dims(1)
            integer(INTD), dimension(2)   :: ov_dims, vo_dims
            integer(INTD), dimension(4)   :: oovv_dims
            integer(INTD)  ::  ov_rank
            type(talsh_tens_t) :: Aux_ov, Aux_vo

            complex(8)                  :: term1, term2, term3, term4, term5, term6

            complex(8), pointer         :: tbar_tens(:,:), tbar_tens2(:,:,:,:)


            one_dims(1) = 1
            ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
            ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,result_tens,one_dims)


            ov_rank=talsh_tensor_rank(tprop_bar_tensor1%t1)
            ierr = talsh_tensor_dimensions(tprop_bar_tensor1%t1,ov_rank,ov_dims)
            ierr=talsh_tensor_construct(Aux_ov,C8,ov_dims,init_val=ZERO)
            vo_dims(1) = ov_dims(2)
            vo_dims(2) = ov_dims(1)
            ierr=talsh_tensor_construct(Aux_vo,C8,vo_dims,init_val=ZERO)


        ! Term1: -Tbar_X*Ty    
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("E()+=Tlx(i,a)*Ty(a,i)",& 
                                    result_tensor,tprop_bar_tensor1%t1,tprop_tensor1%t1,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("E()+=Tlx(i,j,a,b)*Ty(a,b,i,j)", &
                                    result_tensor,tprop_bar_tensor1%t2,tprop_tensor1%t2,scale=MINUS_ONE_QUARTER_C)
            term1=result_tens(1)


        ! Term2: Lambda*Epison_Z    
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("E()+=Epz(a,i)*L(i,a)",result_tensor,interm_t_prop2%epsilon_prop_vo,l1_tensor)
            ierr=talsh_tensor_contract("E()+=Epz(a,b,i,j)*L(i,j,a,b)",result_tensor,&
                                        interm_t_prop2%epsilon_prop_vvoo,l2_tensor,scale=ONE_QUARTER_C)
            term2=result_tens(1)


!           Term3: eta_primeXY*TZ
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("E()+=Sixy(i,a)*Tz(a,i)",result_tensor,sigma_prop_tensor%t1,tprop_tensor2%t1)
            ierr=talsh_tensor_contract("E()+=Sixy(i,j,a,b)*Tz(a,b,i,j)",result_tensor,&
                                        sigma_prop_tensor%t2,tprop_tensor2%t2,scale=ONE_QUARTER_C) 

            term3=result_tens(1)


        ! Term4: -Lamda*TY
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("E()+=L(i,a)*Tl(a,i)",result_tensor,l1_tensor,tprop_tensor1%t1,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("E()+=L(i,j,a,b)*Tl(a,b,i,j)",result_tensor,&
                                        l2_tensor,tprop_tensor1%t2,scale=MINUS_ONE_QUARTER_C)
            term4=result_tens(1)


        ! Term5: Tbar_Z*Epison_X -------------------------------------------------------------------------------
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("E()+=Tl(i,a)*Epz(a,i)",result_tensor,&
                                        tprop_bar_tensor3%t1,interm_t_prop1%epsilon_prop_vo,scale=ONE)
            

            ierr=talsh_tensor_contract("E()+=Tl(i,j,a,b)*Epz(a,b,i,j)",result_tensor,&
                                        tprop_bar_tensor3%t2,interm_t_prop1%epsilon_prop_vvoo,scale=ONE_QUARTER_C)

            term5=result_tens(1)

            
            tem_quarsp_iter=term1*term2+term3+term4*term5

            ierr=talsh_tensor_destruct(Aux_ov)
            ierr=talsh_tensor_destruct(Aux_vo)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------

subroutine construct_TlXAy(exa_input,Interm_t_tensor ,interm_t,tprop_bar_tensor ,sigma_prop_tensor)
    type(exacc_input), intent(in) :: exa_input
    type(talsh_rsp_interm_tens_t) :: Interm_t_tensor
    type(talsh_comm_tens_t)       :: tprop_bar_tensor     
    type(talsh_comm_tens_t)       :: sigma_prop_tensor
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t

    integer(INTD), dimension(4)   :: oovv_dims
    integer(INTD), dimension(2)   :: ov_dims, vv_dims, oo_dims, vo_dims
    integer(C_INT)     :: one_dims(1)
    type(talsh_tens_t)            :: LR_Aux_oo,LR_Aux_vv,LR_Aux_vo
    integer                       :: ierr


    one_dims(1) = 1
    oovv_dims(3:4) = exa_input%nvir
    oovv_dims(1:2) = exa_input%nocc
    ov_dims(2) = exa_input%nvir
    ov_dims(1) = exa_input%nocc

    vv_dims(1:2) = exa_input%nvir
    oo_dims(1:2) = exa_input%nocc
    vo_dims(1) = exa_input%nvir
    vo_dims(2) = exa_input%nocc

    ierr=talsh_tensor_construct(sigma_prop_tensor%t1,C8,ov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(sigma_prop_tensor%t2,C8,oovv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

    ierr=talsh_tensor_construct(LR_Aux_oo,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(LR_Aux_vo,C8,vo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(LR_Aux_vv,C8,vv_dims,init_val=ZERO)
    

    call get_eta_prop_ov_prime(sigma_prop_tensor%t1,Interm_t_tensor,&
                    tprop_bar_tensor%t1,tprop_bar_tensor%t2,interm_t%t2,interm_t%t1,&
                    LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)
   

    call get_eta_prop_oovv_prime(sigma_prop_tensor%t2,Interm_t_tensor,&
                    tprop_bar_tensor%t1,tprop_bar_tensor%t2,interm_t%t1,&
                    LR_Aux_oo,LR_Aux_vv,oovv_dims)


    ierr=talsh_tensor_destruct(LR_Aux_vv)
    ierr=talsh_tensor_destruct(LR_Aux_oo)
    ierr=talsh_tensor_destruct(LR_Aux_vo)
    ierr=talsh_tensor_destruct(one_tensor)


end subroutine construct_TlXAy

!-----------------------------------------------------------------------------------------------------------
  subroutine get_eta_prop_ov_prime(eta_prop_ov_tensor,Interm_t_tensor,l1_tensor,l2_tensor,&
                            t2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)

    type(talsh_tens_t) :: eta_prop_ov_tensor
    type(talsh_rsp_interm_tens_t) :: Interm_t_tensor
    type(talsh_tens_t) :: t2_tensor,t1_tensor
    type(talsh_tens_t) :: l2_tensor,l1_tensor


    type(talsh_tens_t)            :: LR_Aux_oo,LR_Aux_vv,LR_Aux_vo
    integer(INTD), dimension(2)   :: ov_dims
    integer :: ierr

    type(talsh_tens_t)            :: Transpose_vo
    integer(INTD), dimension(2)   :: vo_dims
    integer(C_INT)     :: one_dims(1)
    type(talsh_tens_t) :: uncon_tensor
    logical :: debug_local = .false.


    one_dims(1) = 1
    vo_dims(1) = ov_dims(2)
    vo_dims(2) = ov_dims(1)
    ierr=talsh_tensor_construct(eta_prop_ov_tensor,C8,ov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Transpose_vo,C8,vo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
    ierr=talsh_tensor_contract('A(a,i)+=Y+(i,a)*K()',Transpose_vo,Interm_t_tensor%prop_bar_ov,one_tensor)
    ierr=talsh_tensor_construct(uncon_tensor,C8,one_dims(1:0),init_val=ZERO)
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr_tran ', ierr



!LH Eat1

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

end subroutine get_eta_prop_ov_prime



!-----------------------------------------------------------------------------------------------------------------
subroutine get_eta_prop_oovv_prime(eta_prop_oovv_tensor,Interm_t_tensor,&
    l1_tensor,l2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,oovv_dims)

    type(talsh_tens_t) :: eta_prop_oovv_tensor
    type(talsh_rsp_interm_tens_t) :: Interm_t_tensor
    type(talsh_tens_t) :: t1_tensor
    type(talsh_tens_t) :: l2_tensor,l1_tensor
    type(talsh_tens_t)            :: LR_Aux_oo,LR_Aux_vv
    integer(INTD), dimension(4)   :: oovv_dims
    integer :: ierr
    logical :: debug_local = .false.
    integer(C_INT)     :: one_dims(1)
    type(talsh_tens_t) :: uncon_tensor


    one_dims(1) = 1
    ierr=talsh_tensor_construct(eta_prop_oovv_tensor,C8,oovv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(uncon_tensor,C8,one_dims(1:0),init_val=ZERO)

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
    if(debug_local) print *,'LH_talsh eta_prop_ov_tensor ierr ', ierr

end subroutine get_eta_prop_oovv_prime

!
! below are routines for the contractions needed to assemble the standand CC QR functions
! code by Loic Halbert, 06/05/2023, to be verified
!
!
    subroutine lXAY(Yoo,Yvv,Yov,comm_t,l1_tensor,l2_tensor,nocc,nvir)

    integer(INTD) :: ierr

    type(talsh_tens_t), intent(inout)    :: Yoo, Yvv, Yov

    type(talsh_comm_tens_t), intent(inout) :: comm_t

    integer,             intent(in   )   :: nocc, nvir

    type(talsh_tens_t),  intent(inout)   :: l1_tensor, l2_tensor

    type(talsh_tens_t)  :: lXAY_me_tensor, lXAY_mnef_tensor
      
    type(talsh_tens_t)  :: oo_aux_tensor, vv_aux_tensor

    integer, dimension(2) :: oo_dims, vv_dims, ov_dims
    integer, dimension(4) :: oovv_dims

    oo_dims=(/nocc,nocc/)
    vv_dims=(/nvir,nvir/)
    ov_dims=(/nocc,nvir/)
    oovv_dims=(/nocc,nocc,nvir,nvir/)

    ierr=talsh_tensor_construct(lXAY_me_tensor,C8,ov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(lXAY_mnef_tensor,C8,oovv_dims,init_val=ZERO)

    ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)

    ierr=talsh_tensor_add("W(m,o)+=Y(m,o)",oo_aux_tensor,Yoo)
    if(ierr.ne.0) write(*,*) "ierr b1-a  =",ierr
    ierr=talsh_tensor_contract("W(m,o)+=Y(m,g)*t(g,o)",oo_aux_tensor,Yov,comm_t%t1)
    if(ierr.ne.0) write(*,*) "ierr d3-a  =",ierr 
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(m,o)*L(o,e)",lXAY_me_tensor,oo_aux_tensor,l1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr b1/d3  =",ierr 

    ierr=talsh_tensor_add("W(g,e)+=Y(g,e)",vv_aux_tensor,Yvv)
    if(ierr.ne.0) write(*,*) "ierr c1-a  =",ierr
    ierr=talsh_tensor_contract("W(g,e)+=t(g,o)*Y(o,e)",vv_aux_tensor,comm_t%t1,Yov)
    if(ierr.ne.0) write(*,*) "ierr d3-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(g,e)*L(m,g)",lXAY_me_tensor,vv_aux_tensor,l1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr c1/d3  =",ierr

    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_contract("W(m,p)+=L(o,m,g,h)*T(g,h,o,p)",oo_aux_tensor,l2_tensor,comm_t%t2)
    if(ierr.ne.0) write(*,*) "ierr d4-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(m,p)*Y(p,e)",lXAY_me_tensor,oo_aux_tensor,Yov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr d4-b  =",ierr

    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("W(h,e)+=L(o,p,g,e)*T(g,h,o,p)",vv_aux_tensor,l2_tensor,comm_t%t2)
    if(ierr.ne.0) write(*,*) "ierr d5-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(h,e)*Y(m,h)",lXAY_me_tensor,vv_aux_tensor,Yov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr d5-b  =",ierr

    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,e)*Y(n,f)",lXAY_mnef_tensor,l1_tensor,Yov)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,f)*Y(n,e)",lXAY_mnef_tensor,l1_tensor,Yov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(n,e)*Y(m,f)",lXAY_mnef_tensor,l1_tensor,Yov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(n,f)*Y(m,e)",lXAY_mnef_tensor,l1_tensor,Yov)
    if(ierr.ne.0) write(*,*) "ierr d2-d  =",ierr

    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_add("W(m,o)+=Y(m,o)",oo_aux_tensor,Yoo)
    if(ierr.ne.0) write(*,*) "ierr b2-a =",ierr
    ierr=talsh_tensor_contract("W(m,o)+=Y(m,g)*t(g,o)",oo_aux_tensor,Yov,comm_t%t1)
    if(ierr.ne.0) write(*,*) "ierr d4-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(o,n,e,f)*W(m,o)",lXAY_mnef_tensor,l2_tensor,oo_aux_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(o,m,e,f)*W(n,o)",lXAY_mnef_tensor,l2_tensor,oo_aux_tensor)
    if(ierr.ne.0) write(*,*) "ierr b2/d4  =",ierr

    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_add("W(g,e)+=Y(g,e)",vv_aux_tensor,Yvv)
    if(ierr.ne.0) write(*,*) "ierr c2-a =",ierr
    ierr=talsh_tensor_contract("W(g,e)+=Y(o,e)*t(g,o)",vv_aux_tensor,Yov,comm_t%t1)
    if(ierr.ne.0) write(*,*) "ierr d5-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,n,g,f)*W(g,e)",lXAY_mnef_tensor,l2_tensor,vv_aux_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,n,g,e)*W(g,f)",lXAY_mnef_tensor,l2_tensor,vv_aux_tensor)
    if(ierr.ne.0) write(*,*) "ierr c2/d5  =",ierr


    ierr=talsh_tensor_destruct(oo_aux_tensor)
    ierr=talsh_tensor_destruct(vv_aux_tensor)

    end subroutine lXAY

    subroutine tYFX(Xov,tYs,tYd,l1_tensor,l2_tensor,nocc,nvir)

    integer(INTD) :: ierr

    type(talsh_tens_t), intent(inout)    :: Xov

    type(talsh_tens_t), intent(inout)    :: tYs,tYd        
    integer,             intent(in   )   :: nocc, nvir

    type(talsh_tens_t),  intent(inout)   :: l1_tensor, l2_tensor

    type(talsh_tens_t)  :: tYFX_me_tensor, tYFX_mnef_tensor

    type(talsh_tens_t) :: oovv_aux_tensor
    type(talsh_tens_t) :: ov_aux_tensor, vv_aux_tensor, oo_aux_tensor

    integer, dimension(2) :: oo_dims, vv_dims, ov_dims
    integer, dimension(4) :: oovv_dims

    oo_dims=(/nocc,nocc/)
    vv_dims=(/nvir,nvir/)
    ov_dims=(/nocc,nvir/)
    oovv_dims=(/nocc,nocc,nvir,nvir/)

    ierr=talsh_tensor_construct(tYFX_me_tensor,C8,ov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(tYFX_mnef_tensor,C8,oovv_dims,init_val=ZERO)

    ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(oovv_aux_tensor,C8,oovv_dims,init_val=ZERO)  

    ierr=talsh_tensor_contract("W(k,m,e,c)+=L(k,e)*X(m,c)",oovv_aux_tensor,l1_tensor,Xov)
    if(ierr.ne.0) write(*,*) "ierr e7-a  =",ierr
    ierr=talsh_tensor_contract("W(k,m,e,c)+=L(m,c)*X(k,e)",oovv_aux_tensor,l1_tensor,Xov)
    if(ierr.ne.0) write(*,*) "ierr e8-a  =",ierr           
    ierr=talsh_tensor_contract("tYFX(m,e)+=W(k,m,e,c)*tY(c,k)",tYFX_me_tensor,oovv_aux_tensor,tYs, scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr e7-e8  =",ierr           

    ierr=talsh_tensor_contract("W(d,e)+=L(k,l,c,e)*tY(c,d,k,l)",vv_aux_tensor,l2_tensor,tYd)
    if(ierr.ne.0) write(*,*) "ierr e9-a  =",ierr
    ierr=talsh_tensor_contract("tYFX(m,e)+=W(d,e)*X(m,d)",tYFX_me_tensor,vv_aux_tensor,Xov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr e9-b  =",ierr           
    ierr=talsh_tensor_contract("W(m,l)+=L(k,m,c,d)*tY(c,d,k,l)",oo_aux_tensor,l2_tensor,tYd)
    if(ierr.ne.0) write(*,*) "ierr e10-a  =",ierr
    ierr=talsh_tensor_contract("tYFX(m,e)+=W(m,l)*X(l,e)",tYFX_me_tensor,oo_aux_tensor,Xov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr e10-b  =",ierr           

    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_contract("W(m,k)+=X(m,c)*tY(c,k)",oo_aux_tensor,Xov,tYs)
    if(ierr.ne.0) write(*,*) "ierr e11-a  =",ierr           
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(m,k)*L(k,n,e,f)",tYFX_mnef_tensor,oo_aux_tensor,L2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(n,k)*L(k,m,e,f)",tYFX_mnef_tensor,oo_aux_tensor,L2_tensor)
    if(ierr.ne.0) write(*,*) "ierr e11-b  =",ierr

    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("W(c,e)+=X(k,e)*tY(c,k)",vv_aux_tensor,Xov,tYs)
    if(ierr.ne.0) write(*,*) "ierr e12-a  =",ierr           
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(c,e)*L(m,n,c,f)",tYFX_mnef_tensor,vv_aux_tensor,L2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(c,f)*L(m,n,c,e)",tYFX_mnef_tensor,vv_aux_tensor,L2_tensor)
    if(ierr.ne.0) write(*,*) "ierr e12-b  =",ierr

    ierr=talsh_tensor_destruct(oovv_aux_tensor)
    ierr=talsh_tensor_destruct(vv_aux_tensor)
    ierr=talsh_tensor_destruct(oo_aux_tensor)
    ierr=talsh_tensor_destruct(ov_aux_tensor)

    end subroutine tYFX

    subroutine Gtt(comm_t,int_t,tXs,tXd,tYs,tYd,l1_tensor,l2_tensor,nocc,nvir)

    integer(INTD) :: ierr

     type(talsh_comm_tens_t), intent(inout) :: comm_t
     type(talsh_intg_tens_t), intent(inout) :: int_t

     type(talsh_tens_t), intent(inout)    :: tXs,tXd,tYs,tYd

     integer,             intent(in   )   :: nocc, nvir

     type(talsh_tens_t),  intent(inout) :: l1_tensor, l2_tensor

     type(talsh_tens_t) :: Gtt_me_tensor,Gtt_mnef_tensor

     type(talsh_tens_t) :: G1_tensor

     type(talsh_tens_t) :: ooov_aux_tensor, ovov_aux_tensor, vvvo_aux_tensor
     type(talsh_tens_t) :: oooo_aux_tensor, vovv_aux_tensor, voov_aux_tensor
     type(talsh_tens_t) :: vooo_aux_tensor, vvoo_aux_tensor

     integer, dimension(2) :: ov_dims

     integer, dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
     integer, dimension(4) :: oooo_dims, vovv_dims, voov_dims 
     integer, dimension(4) :: vooo_dims, vvoo_dims

     integer, dimension(4) :: oovv_dims

     ov_dims   = (/nocc, nvir/)

     ooov_dims = (/nocc,nocc,nocc,nvir/)
     ovov_dims = (/nocc,nvir,nocc,nvir/)
     vvvo_dims = (/nvir,nvir,nvir,nocc/)
     oooo_dims = (/nocc,nocc,nocc,nocc/)
     vovv_dims = (/nvir,nocc,nvir,nvir/)
     voov_dims = (/nvir,nocc,nocc,nvir/)
     vooo_dims = (/nvir,nocc,nocc,nocc/)
     vvoo_dims = (/nvir,nvir,nocc,nocc/)

     oovv_dims = (/nocc,nocc,nvir,nvir/)

     ierr=talsh_tensor_construct(Gtt_me_tensor,C8,ov_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(Gtt_mnef_tensor,C8,oovv_dims,init_val=ZERO)

     ierr=talsh_tensor_construct(G1_tensor,C8,ov_dims,init_val=ZERO)

     ierr=talsh_tensor_construct(ooov_aux_tensor,C8,ooov_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(ovov_aux_tensor,C8,ovov_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vvvo_aux_tensor,C8,vvvo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(oooo_aux_tensor,C8,oooo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vovv_aux_tensor,C8,vovv_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(voov_aux_tensor,C8,voov_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vooo_aux_tensor,C8,vooo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vvoo_aux_tensor,C8,vvoo_dims,init_val=ZERO)
     
     !!!     (Gtt)^{m}_{e} 
     ! g2b/k2b inv k3b/g3b

     ! g2b / k2b
     ! 1.Create the loop V*tXs to get intermediate fov
     call Gtt_loop(G1_tensor,int_t,tXs)
     ! 2. Use the fov to get k2b / g2b
     call Gtt_k2k3g2g3(Gtt_me_tensor,int_t,G1_tensor,tYs,tYd,l1_tensor,l2_tensor,nocc,nvir)

     ! k3b / g3b   : inverse
     ierr=talsh_tensor_init(G1_tensor) 
     ! 1. Create the loop  
     call Gtt_loop(G1_tensor,int_t,tYs)
     ! 2. use fov
     call Gtt_k2k3g2g3(Gtt_me_tensor,int_t,G1_tensor,tXs,tXd,l1_tensor,l2_tensor,nocc,nvir)
       
     ! c7b,c8b : same by inversion \tX et tY
     call Gtt_c7c8(Gtt_me_tensor,tXs,tYs,l1_tensor,int_t,nocc,nvir)
     call Gtt_c7c8(Gtt_me_tensor,tXs,tYs,l1_tensor,int_t,nocc,nvir)

     ! f3b, f3d =inv(f4c,f3e) / f3c, f4e =inv( f4b, f4d) 
     ! f3b + f3d + f4c + f3e
            ! create V + V*t
     ierr=talsh_tensor_init(ooov_aux_tensor) 
     ierr=talsh_tensor_add("W(o1,o2,o3,v1)+=V(o1,o2,o3,v1)",ooov_aux_tensor,int_t%ooov)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d W_a  =",ierr         
     ierr=talsh_tensor_contract("W(m,i,o,c)+=T(g,o)*V(m,i,g,c)",ooov_aux_tensor,comm_t%t1,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d W_b  =",ierr

     ! f3b + f3d
     call f3bf3d(ooov_aux_tensor,Gtt_me_tensor,l2_tensor,tXs,tYs,nocc,nvir)
     ! f4c + f3e 
     call f3bf3d(ooov_aux_tensor,Gtt_me_tensor,l2_tensor,tYs,tXs,nocc,nvir)

     ! now use intermediate V+V*t for n2b+n2c
     ierr=talsh_tensor_init(ovov_aux_tensor)
     ierr=talsh_tensor_contract("W(i,c,n,e)+=W(i,k,n,e)*tY(c,k)",ovov_aux_tensor,ooov_aux_tensor,tYs)
     if(ierr.ne.0) write(*,*) "ierr n2b+n2d _c  =",ierr
            ! contraction with tXs ai
     ierr=talsh_tensor_init(vvvo_aux_tensor)
     ierr=talsh_tensor_contract("W(a,c,n,e)+=W(i,c,n,e)*tX(a,i)",vvvo_aux_tensor,oooo_aux_tensor,tXs)
     if(ierr.ne.0) write(*,*) "ierr n2b _d  =",ierr   
     ierr=talsh_tensor_contract("G(m,e)+=W(a,c,n,e)*L(n,m,a,c)",Gtt_me_tensor,vvvo_aux_tensor,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr n2d+n2d _e  =",ierr   
     
     ! f4b + f4d + f3c + f4e
         ! create V + V*t
     ierr=talsh_tensor_init(vovv_aux_tensor)
     ierr=talsh_tensor_add("W(v1,o1,v2,v3)+=V(v1,o1,v2,v3)",vovv_aux_tensor,int_t%vovv)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_a  =",ierr
     ierr=talsh_tensor_contract("W(v1,o1,v2,v3)+=V(o2,o1,v2,v3)*t(v1,o2)",vovv_aux_tensor,int_t%oovv,comm_t%t1)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_b  =",ierr

     ! f4b +f4b
     call f4bf4d(vovv_aux_tensor,Gtt_me_tensor,l2_tensor,tXs,tYs,nocc,nvir)
     ! f3c + f4e
     call f4bf4d(vovv_aux_tensor,Gtt_me_tensor,l2_tensor,tYs,tXs,nocc,nvir)


     ! now use the intermediates for i2b+i2c
     ierr=talsh_tensor_init(voov_aux_tensor)
     ierr=talsh_tensor_contract("W(f,m,i,c)+=W(f,m,a,c)*tX(a,i)",voov_aux_tensor,vovv_aux_tensor,tXs)
     if(ierr.ne.0) write(*,*) "ierr i2b+i2c _c  =",ierr
          ! contraction with tYs ck
     ierr=talsh_tensor_init(vooo_aux_tensor)
     ierr=talsh_tensor_contract("W(f,m,i,k)+=W(f,m,i,c)*tY(c,k)",vooo_aux_tensor,voov_aux_tensor,tYs)
     if(ierr.ne.0) write(*,*) "ierr i2b+i2c _d  =",ierr
     ierr=talsh_tensor_contract("G(m,e)+=W(f,m,i,k)*L(i,k,f,e)",Gtt_me_tensor,vooo_aux_tensor,l2_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr i2b+i2c _e  =",ierr

     ! c9b, g7b, l29, c10b, k7b, l28  - inversion : c11b, k5b, l21, c12b, k10b , l21 bis
     ! c9b, g7b, l29, c10b, k7b, l28
     call Gtt_c9c10g7k7l28l29(Gtt_me_tensor,int_t,tYs,tXd,l2_tensor,nocc,nvir)
     ! c11b, k5b, l21, c12b, k10b, l21 bis
     call Gtt_c9c10g7k7l28l29(Gtt_me_tensor,int_t,tXs,tYd,l2_tensor,nocc,nvir)

     !!!!!!!!!!!!!!
     !!!! \Gd mnef
     !!!!!!!!!!!!!!

     ! l22bis do inversion for l32bis
     ! l23bis  !! SAME intermediate as l22bis   ! do invesrion for l34bis 
      
     call l22bisl23bis(Gtt_mnef_tensor,l2_tensor,int_t,tXs,tYs,nocc,nvir)         
     call l22bisl23bis(Gtt_mnef_tensor,l2_tensor,int_t,tYs,tXs,nocc,nvir)         


     ! i2d
     ierr=talsh_tensor_init(ooov_aux_tensor)
     ierr=talsh_tensor_contract("W(m,n,k,a)+=V(m,n,c,a)*tY(c,k)",ooov_aux_tensor,int_t%oovv,tYs)
     if(ierr.ne.0) write(*,*) "ierr i2d-a  =",ierr
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,n,k,i)+=W(m,n,k,a)*tX(a,i)",oooo_aux_tensor,ooov_aux_tensor,tXs)
     if(ierr.ne.0) write(*,*) "ierr i2d-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(m,n,k,i)*L(k,i,e,f)",Gtt_mnef_tensor,oooo_aux_tensor,l2_tensor)


     ! n2d
     ierr=talsh_tensor_init(vvoo_aux_tensor)
     ierr=talsh_tensor_contract("W(a,c,i,k)+=tx(a,i)*tY(c,k)",vvoo_aux_tensor,tXs,tYs)
     if(ierr.ne.0) write(*,*) "ierr n2d-a  =",ierr
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,n,k,i)+=W(a,c,i,k)*L(m,n,c,a)",oooo_aux_tensor,vvoo_aux_tensor,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr i2d-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(m,n,k,i)*V(k,i,e,f)",Gtt_mnef_tensor,oooo_aux_tensor,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr i2d-c  =",ierr

     ! l28bis=inv(f4e)
     call l28bis(Gtt_mnef_tensor,l2_tensor,int_t,tXs,tYs,nocc,nvir)
     call l28bis(Gtt_mnef_tensor,l2_tensor,int_t,tYs,tXs,nocc,nvir) 

     ierr=talsh_tensor_destruct(G1_tensor)
     ierr=talsh_tensor_destruct(ooov_aux_tensor)
     ierr=talsh_tensor_destruct(ovov_aux_tensor)
     ierr=talsh_tensor_destruct(vvvo_aux_tensor)
     ierr=talsh_tensor_destruct(oooo_aux_tensor)
     ierr=talsh_tensor_destruct(vovv_aux_tensor)
     ierr=talsh_tensor_destruct(voov_aux_tensor)
     ierr=talsh_tensor_destruct(vooo_aux_tensor)
     ierr=talsh_tensor_destruct(vvoo_aux_tensor)
 
     end subroutine Gtt

     subroutine l28bis(G2_tensor,l2_tensor,int_t,tAs,tBs,nocc,nvir)

     real(8)                        :: print_level=2
     integer(INTD) :: ierr

     type(talsh_intg_tens_t), intent(inout) :: int_t

     type(talsh_tens_t),  intent(inout) :: G2_tensor
     type(talsh_tens_t),  intent(inout) :: l2_tensor
     type(talsh_tens_t),  intent(inout) :: tAs,tBs

     integer,             intent(in   ) :: nocc, nvir

     integer, dimension(4) :: oovo_dims,ovvo_dims

     type(talsh_tens_t) :: oovo_aux_tensor ,ovvo_aux_tensor

     oovo_dims = nocc
     oovo_dims(3) = nvir

     ovvo_dims = nocc
     ovvo_dims(2:3) = nvir

     ierr=talsh_tensor_construct(ovvo_aux_tensor,C8,ovvo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(oovo_aux_tensor,C8,oovo_dims,init_val=ZERO)

     ! l28bis / f4e
     ierr=talsh_tensor_init(oovo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,i,e,k)+=V(m,i,e,c)*tY(c,k)",oovo_aux_tensor,int_t%oovv,tBs)
     if(ierr.ne.0) write(*,*) "ierr l28b/f4e-a  =",ierr
     ierr=talsh_tensor_init(ovvo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,a,e,k)+=W(m,i,e,k)*tX(a,i)",ovvo_aux_tensor,oovo_aux_tensor,tAs)
     if(ierr.ne.0) write(*,*) "ierr l28b/f4e-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(m,a,e,k)*L(n,k,f,a)",G2_tensor,ovvo_aux_tensor,l2_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr l28b/f4e-c  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(n,a,e,k)*L(m,k,f,a)",G2_tensor,ovvo_aux_tensor,l2_tensor)
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(m,a,f,k)*L(n,k,e,a)",G2_tensor,ovvo_aux_tensor,l2_tensor)
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(n,a,f,k)*L(m,k,e,a)",G2_tensor,ovvo_aux_tensor,l2_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr l28b/f4e-d  =",ierr


     ierr=talsh_tensor_destruct(ovvo_aux_tensor)
     ierr=talsh_tensor_destruct(oovo_aux_tensor)

     end subroutine

     subroutine l22bisl23bis(G2_tensor,l2_tensor,int_t,tAs,tBs,nocc,nvir)

     real(8)                        :: print_level=2
     integer(INTD) :: ierr

     type(talsh_intg_tens_t), intent(inout) :: int_t

     type(talsh_tens_t),  intent(inout) :: G2_tensor

     type(talsh_tens_t),  intent(inout) :: l2_tensor
     type(talsh_tens_t),  intent(inout) :: tAs,tBs


     integer,             intent(in   ) :: nocc, nvir

     !     tensor dimensions
     integer, dimension(2) :: ov_dims,vv_dims,oo_dims

     type(talsh_tens_t) :: ov_aux_tensor , vv_aux_tensor, oo_aux_tensor

     ov_dims = nvir
     ov_dims(1) = nocc

     vv_dims = nvir
     oo_dims = nocc

     ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)

     ! l22bis do inversionfor l32bis
     ierr=talsh_tensor_init(ov_aux_tensor)
     ierr=talsh_tensor_contract("F(o,f)+=tX(a,i)*(o,i,f,a)",ov_aux_tensor,tAs,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr l22bis-a  =",ierr
     ierr=talsh_tensor_init(vv_aux_tensor)
     ierr=talsh_tensor_contract("F(c,f)+=F(o,f)*tY(c,o)",vv_aux_tensor,ov_aux_tensor,tBs)
     if(ierr.ne.0) write(*,*) "ierr l22bis-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=L(m,n,e,c)*F(c,f)",G2_tensor,l2_tensor,vv_aux_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr l22bis-c  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=L(m,n,f,c)*F(c,e)",G2_tensor,l2_tensor,vv_aux_tensor)
     if(ierr.ne.0) write(*,*) "ierr l22bis-d  =",ierr


     ! l23bis  !! SAME intermediate as l22bis   ! do invesrion for l34bis 
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("F(n,k)+=F(n,c)*tY(c,k)",oo_aux_tensor,ov_aux_tensor,tBs)
     if(ierr.ne.0) write(*,*) "ierr l23bis-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=L(m,k,e,f)*F(n,k)",G2_tensor,l2_tensor,oo_aux_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr l23bis-c  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=L(n,k,e,f)*F(m,k)",G2_tensor,l2_tensor,oo_aux_tensor)
     if(ierr.ne.0) write(*,*) "ierr l23bis-d  =",ierr
     
     ierr=talsh_tensor_destruct(oo_aux_tensor)
     ierr=talsh_tensor_destruct(vv_aux_tensor)
     ierr=talsh_tensor_destruct(ov_aux_tensor)

     end subroutine


     subroutine f4bf4d(vovv_inter,G1_tensor,l2_tensor,tAs,tBs,nocc,nvir)

     real(8)                        :: print_level=2
     integer(INTD) :: ierr

     type(talsh_tens_t),  intent(inout) :: G1_tensor

     type(talsh_tens_t),  intent(inout) :: l2_tensor
     type(talsh_tens_t),  intent(inout) :: tAs,tBs

     type(talsh_tens_t)                 ::  vovv_inter  !Fortran boude : intent

     integer,             intent(in   ) :: nocc, nvir

     !     tensor dimensions
     integer, dimension(4) :: voov_dims, vvov_dims

     type(talsh_tens_t) :: voov_aux_tensor , vvov_aux_tensor

     voov_dims = nvir
     voov_dims(2:3) = nocc

     vvov_dims = nvir
     vvov_dims(3) = nocc
          
     ierr=talsh_tensor_construct(voov_aux_tensor,C8,voov_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vvov_aux_tensor,C8,vvov_dims,init_val=ZERO)

          ! Contraction  with tXs ai
     ierr=talsh_tensor_init(voov_aux_tensor)
     ierr=talsh_tensor_contract("W(f,k,i,e)+=W(f,k,a,e)*tX(a,i)",voov_aux_tensor,vovv_inter,tAs)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_c  =",ierr
          ! contract with tYs ck
     ierr=talsh_tensor_init(vvov_aux_tensor)
     ierr=talsh_tensor_contract("W(f,c,i,e)+=W(f,k,i,e)*tY(c,k)",vvov_aux_tensor,voov_aux_tensor,tBs)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_d  =",ierr
     ierr=talsh_tensor_contract("G(m,e)+=W(f,c,i,e)*L(i,m,f,c)",G1_tensor,voov_aux_tensor,l2_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_e  =",ierr

     ierr=talsh_tensor_destruct(voov_aux_tensor)
     ierr=talsh_tensor_destruct(vvov_aux_tensor)  

     end 

     subroutine f3bf3d(ooov_inter,G1_tensor,l2_tensor,tAs,tBs,nocc,nvir)

     real(8)                        :: print_level=2
     integer(INTD) :: ierr

     type(talsh_tens_t),  intent(inout) :: G1_tensor

     type(talsh_tens_t),  intent(inout) :: l2_tensor
     type(talsh_tens_t),  intent(inout) :: tAs,tBs

     type(talsh_tens_t)                ::  ooov_inter !Fortran boude : intent


     integer,             intent(in   ) :: nocc, nvir

     !     tensor dimensions
     integer, dimension(4) :: oooo_dims, vooo_dims

     type(talsh_tens_t) :: oooo_aux_tensor , vooo_aux_tensor

     oooo_dims = nocc
     vooo_dims =nocc
     vooo_dims(1) =nvir
          
     ierr=talsh_tensor_construct(oooo_aux_tensor,C8,oooo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vooo_aux_tensor,C8,vooo_dims,init_val=ZERO)

            ! contraction with tYs ck
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("W(i,m,n,k)+=W(i,m,n,c)*tY(c,k)",oooo_aux_tensor,ooov_inter,tBs)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d _c  =",ierr
            ! contraction with tXs ai
     ierr=talsh_tensor_init(vooo_aux_tensor)
     ierr=talsh_tensor_contract("W(a,m,n,k)+=W(i,m,n,k)*tX(a,i)",vooo_aux_tensor,oooo_aux_tensor,tAs)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d _d  =",ierr   
     ierr=talsh_tensor_contract("G(m,e)+=W(a,m,n,k)*L(n,k,a,e)",G1_tensor,vooo_aux_tensor,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d _e  =",ierr   
         
     ierr=talsh_tensor_destruct(oooo_aux_tensor)
     ierr=talsh_tensor_destruct(vooo_aux_tensor)  

     end 




     subroutine Gtt_c9c10g7k7l28l29(G1_tensor,int_t,tBs,tAd,l2_tensor,nocc,nvir)

     real(8)                        :: print_level=2         
     integer(INTD) :: ierr

     type(talsh_intg_tens_t), intent(inout) :: int_t

     type(talsh_tens_t),  intent(inout) :: l2_tensor

     integer,             intent(in   ) :: nocc, nvir

     !     tensor dimensions
     integer, dimension(2) :: oo_dims, vv_dims
     type(talsh_tens_t) :: oo_aux_tensor,vv_aux_tensor

     integer, dimension(4) :: vvvv_dims, oooo_dims!, oovv_dims
     integer, dimension(4) :: vovv_dims, ovvo_dims
     integer, dimension(4) :: ooov_dims, ovvv_dims
     integer, dimension(4) :: oovo_dims

     type(talsh_tens_t), intent(inout)    :: G1_tensor

     type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor !, oovv_aux_tensor
     type(talsh_tens_t) :: vovv_aux_tensor,ovvo_aux_tensor
     type(talsh_tens_t) :: ooov_aux_tensor ,ovvv_aux_tensor
     type(talsh_tens_t) :: oovo_aux_tensor
     type(talsh_tens_t), intent(inout) :: tAd,tBs
   

     vv_dims = nvir
     oo_dims = nocc
     ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
     vvvv_dims = nvir
     oooo_dims = nocc
     ierr=talsh_tensor_construct(oooo_aux_tensor,C8,oooo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vvvv_aux_tensor,C8,vvvv_dims,init_val=ZERO)
     vovv_dims = nvir
     vovv_dims(2) = nocc
     ierr=talsh_tensor_construct(vovv_aux_tensor,C8,vovv_dims,init_val=ZERO)
     ovvo_dims = nocc
     ovvo_dims(2:3) = nvir
     ierr=talsh_tensor_construct(ovvo_aux_tensor,C8,ovvo_dims,init_val=ZERO)
     ooov_dims = nocc
     ooov_dims(4) = nvir
     ierr=talsh_tensor_construct(ooov_aux_tensor,C8,ooov_dims,init_val=ZERO)
     ovvv_dims = nvir
     ovvv_dims(1) = nocc
     ierr=talsh_tensor_construct(ovvv_aux_tensor,C8,ovvv_dims,init_val=ZERO)
     oovo_dims = nocc
     oovo_dims(2) = nvir
     ierr=talsh_tensor_construct(oovo_aux_tensor,C8,oovo_dims,init_val=ZERO)

     ! c9b, g7b , l29
     ierr=talsh_tensor_init(vovv_aux_tensor)
     ierr=talsh_tensor_contract("W(c,m,a,e)+=tY(c,k)*V(k,m,a,e)",vovv_aux_tensor,tBs,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr c9b, g7b, l29-a  =",ierr

     ierr=talsh_tensor_init(vv_aux_tensor)
     ierr=talsh_tensor_contract("F(a,c)+=tX(a,b,i,j)*L(i,j,b,c)",vv_aux_tensor,tAd,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr c9b-a  =",ierr
     ierr=talsh_tensor_contract("G(m,e)+=W(c,m,a,e)*F(a,c)",G1_tensor,vovv_aux_tensor,vv_aux_tensor,scale=MINUS_ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr c9b-b  =",ierr
    
     ierr=talsh_tensor_init(vvvv_aux_tensor)
     ierr=talsh_tensor_contract("W(a,b,c,e)+=tX(a,b,i,j)*L(i,j,c,e)",vvvv_aux_tensor,tAd,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr g7b-a  =",ierr
     ierr=talsh_tensor_contract("G(m,e)+=W(a,b,c,e)*W(c,m,a,b)",G1_tensor,vvvv_aux_tensor,vovv_aux_tensor,scale=ONE_QUARTER_C)
     if(ierr.ne.0) write(*,*) "ierr g7b-b  =",ierr

     ierr=talsh_tensor_init(ovvo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,a,c,i)+=tX(a,b,i,j)*L(j,m,b,c)",ovvo_aux_tensor,tAd,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr l29-a  =",ierr
     ierr=talsh_tensor_contract("G(m,e)+=W(m,a,c,i)*W(c,m,a,e)",G1_tensor,ovvo_aux_tensor,vovv_aux_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr l29-b  =",ierr

     ! c10b ,k7b,l28
     ierr=talsh_tensor_init(ooov_aux_tensor)
     ierr=talsh_tensor_contract("W(i,m,k,e)+=tY(c,k)*V(i,m,c,e)",ooov_aux_tensor,tBs,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr c10b ,k7b,l28-a  =",ierr

     ierr=talsh_tensor_init(oo_aux_tensor) 
     ierr=talsh_tensor_contract("F(k,i)+=tX(a,b,i,j)*L(j,k,b,a)",oo_aux_tensor,tAd,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr c10b-a  =",ierr 
     ierr=talsh_tensor_contract("G(m,e)+=W(i,m,k,e)*F(k,i)",G1_tensor,ooov_aux_tensor,oo_aux_tensor,scale=MINUS_ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr c10b-b  =",ierr

     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("W(k,m,i,j)+=tX(a,b,i,j)*L(k,m,a,b)",oooo_aux_tensor,tAd,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr k7b-a  =",ierr 
     ierr=talsh_tensor_contract("G(m,e)+=W(i,j,k,e)*W(k,m,i,j)",G1_tensor,ooov_aux_tensor,oooo_aux_tensor,scale=ONE_QUARTER_C)
     if(ierr.ne.0) write(*,*) "ierr k7b-b  =",ierr 

     ierr=talsh_tensor_init(ovvo_aux_tensor)
     ierr=talsh_tensor_contract("W(k,a,e,i)+=tX(a,b,i,j)*L(j,k,b,e)",ovvo_aux_tensor,tAd,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr l28-a  =",ierr 
     ierr=talsh_tensor_contract("G(m,e)+=W(i,m,a,k)*W(k,a,e,i)",G1_tensor,ovvo_aux_tensor,oovo_aux_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr l28-b  =",ierr 

     ierr=talsh_tensor_destruct(vv_aux_tensor)
     ierr=talsh_tensor_destruct(oo_aux_tensor)
     ierr=talsh_tensor_destruct(oooo_aux_tensor)
     ierr=talsh_tensor_destruct(vvvv_aux_tensor)
     ierr=talsh_tensor_destruct(vovv_aux_tensor)
     ierr=talsh_tensor_destruct(ovvo_aux_tensor)
     ierr=talsh_tensor_destruct(ooov_aux_tensor)
     ierr=talsh_tensor_destruct(ovvv_aux_tensor)
     ierr=talsh_tensor_destruct(oovo_aux_tensor)

     end subroutine Gtt_c9c10g7k7l28l29


     subroutine Gtt_c7c8(G1_tensor,tAs,tBs,l1_tensor,int_t,nocc,nvir)

     type(talsh_tens_t), intent(inout)    :: G1_tensor
     type(talsh_intg_tens_t), intent(inout) :: int_t
     type(talsh_tens_t), intent(inout)    :: tAs,tBs  
     type(talsh_tens_t)                   :: ooov_aux_tensor,oo_aux_tensor 
     type(talsh_tens_t),  intent(inout) :: l1_tensor 
     integer(INTD) :: ierr
     integer,             intent(in   ) :: nocc, nvir

     ierr=talsh_tensor_construct(ooov_aux_tensor,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
     ierr=talsh_tensor_construct(oo_aux_tensor,C8,(/nocc,nocc/),init_val=ZERO)


     ierr=talsh_tensor_contract("W(k,m,i,e)+=tX(a,i)*V(k,m,a,e)",ooov_aux_tensor,tAs,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr Gtt c7bc8b-a  =",ierr         
     ierr=talsh_tensor_contract("F(i,k)+=L(i,c)*T(c,k)",oo_aux_tensor,l1_tensor,tBs)
     if(ierr.ne.0) write(*,*) "ierr Gtt c7bc8b-b =",ierr         
     ierr=talsh_tensor_contract("G(m,e)+=W(k,m,i,e)*F(i,k)",G1_tensor,ooov_aux_tensor,oo_aux_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr Gtt c7bc8b-c  =",ierr         

     ierr=talsh_tensor_destruct(ooov_aux_tensor)
     ierr=talsh_tensor_destruct(oo_aux_tensor)

     end subroutine Gtt_c7c8


     subroutine Gtt_loop(aux_ov_tensor,int_t,tAs)

     type(talsh_tens_t), intent(inout)    :: aux_ov_tensor
     type(talsh_intg_tens_t), intent(inout) :: int_t
     type(talsh_tens_t)                   :: tAs  
     integer(INTD) :: ierr

     ierr=talsh_tensor_contract("G(m,c)+=V(m,i,c,a)*tAs(a,i)",aux_ov_tensor,int_t%oovv,tAs)
     if(ierr.ne.0) write(*,*) "ierr Gtt_loop =",ierr 

     end subroutine Gtt_loop

     subroutine Gtt_k2k3g2g3(Gtt_me_tensor,int_t,G1_tensor,tBs,tBd,l1_tensor,l2_tensor,nocc,nvir)

     type(talsh_tens_t), intent(inout)    :: G1_tensor
     type(talsh_tens_t), intent(inout)    :: l1_tensor,l2_tensor
     type(talsh_intg_tens_t), intent(inout) :: int_t
     type(talsh_tens_t), intent(inout)    :: tBs, tBd
     integer(INTD)                        :: ierr
     integer,             intent(in   )   :: nocc, nvir         

     integer :: print_level = 2      

     integer, dimension(2) :: oo_dims, vv_dims
     type(talsh_tens_t) :: oo_aux_tensor,vv_aux_tensor

     type(talsh_tens_t), intent(inout) :: Gtt_me_tensor


     oo_dims = nocc
     vv_dims = nvir

     ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
     ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)

     ! g2b / g3b
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("Foo(m,k)+=G(m,c)*tY(c,k)",oo_aux_tensor,G1_tensor,tBs)
     if(ierr.ne.0) write(*,*) "ierr g2b-a  =",ierr
     ierr=talsh_tensor_contract("Gtt(m,e)+=F(m,k)*L(k,e)",Gtt_me_tensor,oo_aux_tensor,l1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr g2b-b  =",ierr
     ! k2b / k3b
     ierr=talsh_tensor_init(oo_aux_tensor) 
     ierr=talsh_tensor_contract("Foo(m,k)+=L(m,c)*tY(c,k)",oo_aux_tensor,l1_tensor,tBs)
     if(ierr.ne.0) write(*,*) "ierr k2b-a  =",ierr
     ierr=talsh_tensor_contract("Gtt(m,e)+=Foo(m,k)*G(k,e)",Gtt_me_tensor,oo_aux_tensor,G1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr k2b-b  =",ierr
      
     ierr=talsh_tensor_destruct(oo_aux_tensor)
     ierr=talsh_tensor_destruct(vv_aux_tensor)


     end subroutine Gtt_k2k3g2g3



end module

