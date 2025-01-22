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
                                ONE_HALF=(0.5D0,0.D0), IMAGINARY=(0.0D0,1.0D0), &
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

        complex(kind=8)               :: phase_factor 

        complex(kind=8)               :: eom_qr
        real(kind=8)                  :: tem_real, tem_imag

        type(talsh_tens_t)          :: result_tensor
        integer(C_INT)              :: one_dims(1)
        complex(8), pointer         :: result_tens(:)
        type(C_PTR):: body_p


        call print_date('Evaluate CC-CI QR')

        quarsp = ZERO
        quarsp_con = ZERO

        call get_quarspf(exa_input, interm_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,&
                        Tprop_tensors,Tprop_bar_tensors,& 
                        i, j, k, l1_tensor, l2_tensor, quarsp)

        call get_quarspf(exa_input, interm_t, &
                        Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,&
                        Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                        i, j, k, l1_tensor, l2_tensor, quarsp_con)

        print *, "Number of time-antisymmetric operators in response function",exa_input%N_Time_anti
        if (exa_input%N_Time_anti.ge.1) &
        print *, "Adding phase to response function components (due to how DIRAC scales time-antisymmetric operators)"

        phase_factor = IMAGINARY**(exa_input%N_Time_anti)

        quarsp     = quarsp*phase_factor
        quarsp_con = quarsp_con*phase_factor

        print *, "quadratic response: +w", quarsp
        print *, "quadratic response: -w", quarsp_con

        quarsp_con=conjg(quarsp_con)
        eom_qr=(quarsp+quarsp_con)*0.5

        print *, "EOM-QR after symmetrization ",eom_qr

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

end module

