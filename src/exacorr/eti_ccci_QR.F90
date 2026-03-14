module eti_ccci_QR
    use tensor_algebra
    use exacorr_tensor_interface
    use exacorr_datatypes
    use exacorr_utils
    use eti_common_routines
    use exacorr_global
    use exacorr_scaling_parameters
    use exacorr_cfg, only: exacc_input 

    implicit none

    type(tensor_t) :: one_tensor

    private

    public get_ccci_qr

    contains


subroutine get_ccci_qr(exa_input, interm_t,&
    Interm_t_tensors_A,Interm_t_tensors_B, Interm_t_tensors_C,& 
    Tprop_tensors, Tprop_bar_tensors,& 
    Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
    i, j, k, l1_tensor, l2_tensor, cc_qr)

!   This subroutine drives the CC-CI Quadratic response property calculation
!   Written by Xiang Yuan Janunary 2022, Lille, France
!   At the moment, there is no call for this code. It's just a backup.

use eti_sigma_eomccsd
type(exacc_input), intent(in) :: exa_input
type(eti_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
type(tensor_t), intent(inout)           :: l2_tensor
type(tensor_t), intent(inout)           :: l1_tensor
type(eti_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:)
type(eti_rsp_interm_tens_t), allocatable :: Interm_t_tensors_B(:)
type(eti_rsp_interm_tens_t), allocatable :: Interm_t_tensors_C(:)


type(eti_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
type(eti_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR

integer                       :: i, j, k

complex(kind=8)               :: quarsp, quarsp_con

integer                       :: ierr, T

complex(kind=8)               :: cc_qr
real(kind=8)                  :: tem_real, tem_imag


! type(eti_comm_tens_t)     :: bvec_tensor
type(tensor_t)          :: result_tensor
integer(C_INT)              :: one_dims(1)
complex(8), pointer         :: result_tens(:)
type(C_PTR):: body_p


call print_date('Evaluate CC-CI QR')


call get_quarspf_ccci(exa_input, interm_t, &
    Interm_t_tensors_A,Interm_t_tensors_B,Interm_t_tensors_C,&
    Tprop_tensors,Tprop_bar_tensors,& 
    i, j, k, l1_tensor, l2_tensor, quarsp)

call get_quarspf_ccci(exa_input, interm_t, &
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
    cc_qr=(quarsp+quarsp_con)*0.5

end subroutine get_ccci_qr
!------------------------------------------------------------------------------------------------------------------


subroutine get_quarspf_ccci(exa_input, interm_t,&
    Interm_t_tensors_A, Interm_t_tensors_B, Interm_t_tensors_C,&
    Tprop_tensors,Tprop_bar_tensors,&
    i, j, k, l1_tensor, l2_tensor, quarsp)         
type(exacc_input), intent(in) :: exa_input
type(eti_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
type(tensor_t), intent(inout)           :: l2_tensor
type(tensor_t), intent(inout)           :: l1_tensor
type(eti_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:)
type(eti_rsp_interm_tens_t), allocatable :: Interm_t_tensors_B(:)
type(eti_rsp_interm_tens_t), allocatable :: Interm_t_tensors_C(:)
type(eti_comm_tens_t), allocatable :: Tprop_tensors(:),Tprop_bar_tensors(:)
integer                            :: i, j, k


type(eti_comm_tens_t)       :: sigma_prop_tensor12, sigma_prop_tensor13, sigma_prop_tensor23
type(eti_comm_tens_t)       :: sigma_prop_tensor21, sigma_prop_tensor31, sigma_prop_tensor32
complex(kind=8)               :: quarsp
complex(kind=8)               :: tem_quarsp_iter
type(eti_comm_tens_t)       :: bvec_tensor


!XYZ   X is 1; Y is 2; Z is 3  Ax*Tly
call construct_AyTz(exa_input,Interm_t_tensors_A(i),Tprop_tensors(3),sigma_prop_tensor13)
call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(2), Tprop_tensors(3), &
                    Tprop_bar_tensors(1), Tprop_bar_tensors(2), Tprop_bar_tensors(3), &
                    Interm_t_tensors_A(i), Interm_t_tensors_C(k), &
                    sigma_prop_tensor13,l1_tensor,l2_tensor, bvec_tensor)
            
quarsp=quarsp+tem_quarsp_iter



!XZY   Ax*Tlz
call construct_AyTz(exa_input,Interm_t_tensors_A(i),Tprop_tensors(2),sigma_prop_tensor12)
call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(3), Tprop_tensors(2), &
                    Tprop_bar_tensors(1), Tprop_bar_tensors(3), Tprop_bar_tensors(2), &
                    Interm_t_tensors_A(i), Interm_t_tensors_B(j), &
                    sigma_prop_tensor12,l1_tensor,l2_tensor, bvec_tensor)

quarsp=quarsp+tem_quarsp_iter


!YXZ   Ay*Tlx
call construct_AyTz(exa_input,Interm_t_tensors_B(j),Tprop_tensors(3),sigma_prop_tensor23)
call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(1), Tprop_tensors(3), &
                    Tprop_bar_tensors(2), Tprop_bar_tensors(1), Tprop_bar_tensors(3), &
                    Interm_t_tensors_B(j), Interm_t_tensors_C(k), &
                    sigma_prop_tensor23,l1_tensor,l2_tensor, bvec_tensor)

quarsp=quarsp+tem_quarsp_iter


!YZX  Ay*Tlz
call construct_AyTz(exa_input,Interm_t_tensors_B(j),Tprop_tensors(1),sigma_prop_tensor21)
call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(3), Tprop_tensors(1), &
                    Tprop_bar_tensors(2), Tprop_bar_tensors(3), Tprop_bar_tensors(1), &
                    Interm_t_tensors_B(j), Interm_t_tensors_A(i), &
                    sigma_prop_tensor21,l1_tensor,l2_tensor, bvec_tensor)


quarsp=quarsp+tem_quarsp_iter


!ZYX Az*Tly
call construct_AyTz(exa_input,Interm_t_tensors_C(k),Tprop_tensors(1),sigma_prop_tensor31)
call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(2), Tprop_tensors(1), &
                    Tprop_bar_tensors(3), Tprop_bar_tensors(2), Tprop_bar_tensors(1), &
                    Interm_t_tensors_C(k), Interm_t_tensors_A(i), &
                    sigma_prop_tensor31,l1_tensor,l2_tensor, bvec_tensor)

quarsp=quarsp+tem_quarsp_iter


!ZXY Az*Tlx
call construct_AyTz(exa_input,Interm_t_tensors_C(k),Tprop_tensors(2),sigma_prop_tensor32)
call get_quarsp_driver(tem_quarsp_iter, Tprop_tensors(1), Tprop_tensors(2), &
                    Tprop_bar_tensors(3), Tprop_bar_tensors(1), Tprop_bar_tensors(2), &
                    Interm_t_tensors_C(k), Interm_t_tensors_B(j), &
                    sigma_prop_tensor32,l1_tensor,l2_tensor, bvec_tensor)

quarsp=quarsp+tem_quarsp_iter

end subroutine get_quarspf_ccci

!------------------------------------------------------------------------------------------------------------------

subroutine get_quarsp_driver(tem_quarsp_iter,tprop_tensor1, tprop_tensor2, &
        tprop_bar_tensor1, tprop_bar_tensor2, tprop_bar_tensor3, &
        interm_t_prop1, interm_t_prop2, &
        sigma_prop_tensor,l1_tensor,l2_tensor,bvec_tensor)

type(eti_rsp_interm_tens_t)    :: interm_t_prop1, interm_t_prop2
type(tensor_t)               :: l2_tensor
type(tensor_t)               :: l1_tensor
type(eti_comm_tens_t)          :: tprop_tensor1, tprop_tensor2
type(eti_comm_tens_t)          :: tprop_bar_tensor1, tprop_bar_tensor2, tprop_bar_tensor3
type(eti_comm_tens_t)          :: sigma_prop_tensor
complex(kind=8)                  :: tem_quarsp_iter
type(eti_comm_tens_t),OPTIONAL          :: bvec_tensor

!  Point connected with tensors 
type(tensor_t)          :: result_tensor
integer(C_INT)              :: ierr
complex(8), pointer         :: result_tens(:)
integer                     :: initialize

type(C_PTR):: body_p, body_p2, body_p3

integer(C_INT)              :: one_dims(1)
integer(INTD), allocatable  :: ov_dims(:)
integer(INTD), dimension(2) :: vo_dims
integer(INTD), dimension(4) :: oovv_dims
integer(INTD)  ::  ov_rank
type(tensor_t) :: Aux_ov, Aux_vo

complex(8)                  :: term1, term2, term3, term4, term5, term6

complex(8), pointer         :: tbar_tens(:,:), tbar_tens2(:,:,:,:)






one_dims(1) = 1
ierr=construct_tensor_t(result_tensor,C8,one_dims(1:0),init_val=ZERO)
ierr=direct_data_pointer(result_tensor,body_p)
call c_f_pointer(body_p,result_tens,one_dims)

ierr = get_nmodes(tprop_bar_tensor1%t1,ov_rank)
ierr = get_extents(tprop_bar_tensor1%t1,ov_dims)

ierr=construct_tensor_t(Aux_ov,C8,ov_dims,init_val=ZERO)
vo_dims(1) = ov_dims(2)
vo_dims(2) = ov_dims(1)
ierr=construct_tensor_t(Aux_vo,C8,vo_dims,init_val=ZERO)


! Asymmetric LR
  ierr=t_einsum("E()+=nx(i,a)*Ty(a,i)",& 
                              result_tensor,bvec_tensor%t1,tprop_tensor1%t1,scale=ONE)

  ierr=t_einsum("E()+=nx(i,j,a,b)*Ty(a,b,i,j)", &
                              result_tensor,bvec_tensor%t2,tprop_tensor1%t2,scale=ONE_QUARTER)

  print *, "n*t of Aym LR:", result_tens(1)
  ierr=t_einsum("E()+=tbar(i,a)*Ep(a,i)",& 
                              result_tensor,tprop_bar_tensor1%t1,interm_t_prop1%epsilon_prop_vo,scale=ONE)

  ierr=t_einsum("E()+=tbar(i,j,a,b)*Ep(a,b,i,j)", &
                              result_tensor,tprop_bar_tensor1%t2,interm_t_prop1%epsilon_prop_vvoo,scale=ONE_QUARTER)

  print *, "Aym LR:", result_tens(1)



! Term1: -Tbar_X*Ty    
ierr=initialize_tensor_t(result_tensor)
ierr=t_einsum("E()+=Tlx(i,a)*Ty(a,i)",& 
                result_tensor,tprop_bar_tensor1%t1,tprop_tensor1%t1,scale=MINUS_ONE)
ierr=t_einsum("E()+=Tlx(i,j,a,b)*Ty(a,b,i,j)", &
                result_tensor,tprop_bar_tensor1%t2,tprop_tensor1%t2,scale=MINUS_ONE_QUARTER)
term1=result_tens(1)
! print *, "term1: ", term1


! Term2: Lambda*Epison_Z    
ierr=initialize_tensor_t(result_tensor)
ierr=t_einsum("E()+=Epz(a,i)*L(i,a)",result_tensor,interm_t_prop2%epsilon_prop_vo,l1_tensor)
ierr=t_einsum("E()+=Epz(a,b,i,j)*L(i,j,a,b)",result_tensor,&
                    interm_t_prop2%epsilon_prop_vvoo,l2_tensor,scale=ONE_QUARTER)
term2=result_tens(1)
! print *, "term2: ", term2


!           Term3: eta_primeXY*TZ
ierr=initialize_tensor_t(result_tensor)
ierr=t_einsum("E()+=Sixy(i,a)*Tz(a,i)",result_tensor,sigma_prop_tensor%t1,tprop_tensor2%t1)
ierr=t_einsum("E()+=Sixy(i,j,a,b)*Tz(a,b,i,j)",result_tensor,&
                    sigma_prop_tensor%t2,tprop_tensor2%t2,scale=ONE_QUARTER) 

term3=result_tens(1)
! print *, "term3: ", term3



! Term5: -Lamda*TY
ierr=initialize_tensor_t(result_tensor)
ierr=t_einsum("E()+=L(i,a)*Tl(a,i)",result_tensor,l1_tensor,tprop_tensor1%t1,scale=MINUS_ONE)
ierr=t_einsum("E()+=L(i,j,a,b)*Tl(a,b,i,j)",result_tensor,&
                    l2_tensor,tprop_tensor1%t2,scale=MINUS_ONE_QUARTER)
term5=result_tens(1)
! print *, "term5: ", term5


!EOM-QR-----! Term6: Tbar_Z*Epison_X -------------------------------------------------------------------------------
ierr=initialize_tensor_t(result_tensor)
ierr=t_einsum("E()+=Tl(i,a)*Epz(a,i)",result_tensor,&
                    tprop_bar_tensor3%t1,interm_t_prop1%epsilon_prop_vo,scale=ONE)


ierr=t_einsum("E()+=Tl(i,j,a,b)*Epz(a,b,i,j)",result_tensor,&
                    tprop_bar_tensor3%t2,interm_t_prop1%epsilon_prop_vvoo,scale=ONE_QUARTER)

term6=result_tens(1)
! print *, "term6: ", term6

!CC-CI-QR-----! Term6: Tbar_Z*Epison_X -------------------------------------------------------------------------------

ierr=initialize_tensor_t(result_tensor)
ierr=t_einsum("E()+=T(a,i)*Epz(a,i)",result_tensor,&
                tprop_tensor2%t1,interm_t_prop1%O_vo,scale=ONE)

term6=result_tens(1)
print *, "term6: ", term6


tem_quarsp_iter=term1*term2+term3+term5*term6

ierr=destroy_tensor_t(Aux_ov)
ierr=destroy_tensor_t(Aux_vo)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------

subroutine construct_TlXAy(exa_input,Interm_t_tensor ,interm_t,tprop_bar_tensor ,sigma_prop_tensor)
type(exacc_input), intent(in) :: exa_input
type(eti_rsp_interm_tens_t) :: Interm_t_tensor
type(eti_comm_tens_t)       :: tprop_bar_tensor     
type(eti_comm_tens_t)       :: sigma_prop_tensor
type(eti_intermediates_tens_t),intent(inout) :: interm_t

integer(INTD), dimension(4)   :: oovv_dims
integer(INTD), dimension(2)   :: ov_dims, vv_dims, oo_dims, vo_dims
integer(C_INT)     :: one_dims(1)
type(tensor_t)            :: LR_Aux_oo,LR_Aux_vv,LR_Aux_vo
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

ierr=construct_tensor_t(sigma_prop_tensor%t1,C8,ov_dims,init_val=ZERO)
ierr=construct_tensor_t(sigma_prop_tensor%t2,C8,oovv_dims,init_val=ZERO)
ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)

ierr=construct_tensor_t(LR_Aux_oo,C8,oo_dims,init_val=ZERO)
ierr=construct_tensor_t(LR_Aux_vo,C8,vo_dims,init_val=ZERO)
ierr=construct_tensor_t(LR_Aux_vv,C8,vv_dims,init_val=ZERO)


! call get_sigma_eomccsd_bar_prop(Interm_t_tensor, tprop_bar_tensor, sigma_prop_tensor)


call get_eta_prop_ov_prime(sigma_prop_tensor%t1,Interm_t_tensor,&
tprop_bar_tensor%t1,tprop_bar_tensor%t2,interm_t%t2,interm_t%t1,&
LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)


call get_eta_prop_oovv_prime(sigma_prop_tensor%t2,Interm_t_tensor,&
tprop_bar_tensor%t1,tprop_bar_tensor%t2,interm_t%t1,&
LR_Aux_oo,LR_Aux_vv,oovv_dims)


ierr=destroy_tensor_t(LR_Aux_vv)
ierr=destroy_tensor_t(LR_Aux_oo)
ierr=destroy_tensor_t(LR_Aux_vo)
ierr=destroy_tensor_t(one_tensor)


end subroutine construct_TlXAy

!-----------------------------------------------------------------------------------------------------------
subroutine get_eta_prop_ov_prime(eta_prop_ov_tensor,Interm_t_tensor,l1_tensor,l2_tensor,&
        t2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)

type(tensor_t) :: eta_prop_ov_tensor
type(eti_rsp_interm_tens_t) :: Interm_t_tensor
type(tensor_t) :: t2_tensor,t1_tensor
type(tensor_t) :: l2_tensor,l1_tensor


type(tensor_t)            :: LR_Aux_oo,LR_Aux_vv,LR_Aux_vo
integer(INTD), dimension(2)   :: ov_dims
integer :: ierr

type(tensor_t)            :: Transpose_vo
integer(INTD), dimension(2)   :: vo_dims
integer(C_INT)     :: one_dims(1)
type(tensor_t) :: uncon_tensor
logical :: debug_local = .false.


one_dims(1) = 1
vo_dims(1) = ov_dims(2)
vo_dims(2) = ov_dims(1)
ierr=construct_tensor_t(eta_prop_ov_tensor,C8,ov_dims,init_val=ZERO)
ierr=construct_tensor_t(Transpose_vo,C8,vo_dims,init_val=ZERO)
ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
ierr=t_einsum('A(a,i)+=Y+(i,a)*K()',Transpose_vo,Interm_t_tensor%prop_bar_ov,one_tensor)
ierr=construct_tensor_t(uncon_tensor,C8,one_dims(1:0),init_val=ZERO)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_tran ', ierr



!LH Eat1
!LH     \Eta_1 ia=+\Ys ia
! ierr=t_einsum('E(i,a)+=Y(i,a)',eta_prop_ov_tensor,Interm_t_tensor%prop_bar_ov,scale=ONE)                    
! if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o1 ', ierr

!XY  \Eta_1 ia=+\lambdad ijab \epison ej
ierr=t_einsum('E(i,a)+=L(i,j,a,e)*Epis(e,j)',eta_prop_ov_tensor,l2_tensor,&
        Interm_t_tensor%epsilon_prop_vo,scale=ONE)                    
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o1 ', ierr

!LH     \Eta_1 ia=+\Ys ea \lambdas ie
ierr=t_einsum('E(i,a)+=Y(e,a)*L(i,e)',eta_prop_ov_tensor,Interm_t_tensor%prop_bar_vv,l1_tensor,scale=ONE)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o3 ', ierr

!LH     \Eta_1 ia=-\Ys im \lambdas ma
ierr=t_einsum('E(i,a)+=Y(i,m)*L(m,a)',eta_prop_ov_tensor,Interm_t_tensor%prop_bar_oo,l1_tensor,scale=MINUS_ONE)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o4 ', ierr

!LH     \Eta_1=-\Ys ma \ts em \lambdas ie
ierr=initialize_tensor_t(LR_Aux_oo)       
ierr=t_einsum('H(i,m)+=T(e,m)*L(i,e)',LR_Aux_oo,t1_tensor,l1_tensor,scale=ONE)
ierr=t_einsum('E(i,a)+=H(i,m)*Y(m,a)',eta_prop_ov_tensor,LR_Aux_oo,Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o5 ', ierr

!LH     \Eta_1=-\Ys ie \ts em \llambdas ma
ierr=initialize_tensor_t(LR_Aux_oo)
ierr=t_einsum('H(i,m)+=Y(i,e)*T(e,m)',Lr_Aux_oo,Interm_t_tensor%prop_bar_ov,t1_tensor,scale=ONE)
ierr=t_einsum('E(i,a)+=H(i,m)*L(m,a)',eta_prop_ov_tensor,LR_Aux_oo,l1_tensor,scale=MINUS_ONE)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o6 ', ierr


!LH     \Eta_1= -one_half \td femn \lambdad mife \Ys na 
ierr=initialize_tensor_t(LR_Aux_oo)
ierr=t_einsum('H(i,n)+=T(f,e,m,n)*L(m,i,f,e)',LR_Aux_oo,t2_tensor,l2_tensor)
ierr=t_einsum('E(i,a)+=H(i,n)*Y(n,a)',& 
eta_prop_ov_tensor,LR_Aux_oo,Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE_HALF)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o28 ', ierr

!LH     \Eta_1=-one_half \td  fenm \lambdad nmfa \Ys ie 
ierr=initialize_tensor_t(LR_Aux_vv)
ierr=t_einsum('H(e,a)+=T(f,e,n,m)*L(n,m,f,a)',LR_Aux_vv,t2_tensor,l2_tensor)
ierr=t_einsum('E(i,a)+=H(e,a)*Y(i,e)',&
eta_prop_ov_tensor,LR_Aux_vv,Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE_HALF)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr_o29 ', ierr

end subroutine get_eta_prop_ov_prime



!-----------------------------------------------------------------------------------------------------------------
subroutine get_eta_prop_oovv_prime(eta_prop_oovv_tensor,Interm_t_tensor,&
l1_tensor,l2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,oovv_dims)

type(tensor_t) :: eta_prop_oovv_tensor
! type(talsh_orb_prop_integ_tens_t) :: orb_int_tensor
type(eti_rsp_interm_tens_t) :: Interm_t_tensor
type(tensor_t) :: t1_tensor
type(tensor_t) :: l2_tensor,l1_tensor
type(tensor_t)            :: LR_Aux_oo,LR_Aux_vv
integer(INTD), dimension(4)   :: oovv_dims
integer :: ierr
logical :: debug_local = .false.
integer(C_INT)     :: one_dims(1)
type(tensor_t) :: uncon_tensor


one_dims(1) = 1
ierr=construct_tensor_t(eta_prop_oovv_tensor,C8,oovv_dims,init_val=ZERO)
ierr=construct_tensor_t(uncon_tensor,C8,one_dims(1:0),init_val=ZERO)

!LH     \Eta_2=+ \Pmm ij \Pmm ab \lambdas ia \Ys jb
ierr=t_einsum('E(i,j,a,b)+=L(i,a)*Y(j,b)',eta_prop_oovv_tensor,l1_tensor,Interm_t_tensor%prop_bar_ov,scale=ONE)
ierr=t_einsum('E(i,j,a,b)+=L(j,a)*Y(i,b)',eta_prop_oovv_tensor,l1_tensor,&
                                    Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE)
ierr=t_einsum('E(i,j,a,b)+=L(i,b)*Y(j,a)',eta_prop_oovv_tensor,l1_tensor,&
                                    Interm_t_tensor%prop_bar_ov,scale=MINUS_ONE)
ierr=t_einsum('E(i,j,a,b)+=L(j,b)*Y(i,a)',eta_prop_oovv_tensor,l1_tensor,Interm_t_tensor%prop_bar_ov,scale=ONE)
if(debug_local) print *,'LH eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=+\Pmm ab \lambdad ijae \Ys eb 
ierr=t_einsum('E(i,j,a,b)+=L(i,j,a,e)*Y(e,b)', &
eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_vv,scale=ONE)
ierr=t_einsum('E(i,j,a,b)+=L(i,j,b,e)*Y(e,a)', & 
eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_vv,scale=MINUS_ONE)
if(debug_local) print *,'LH eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2= - \Pmm ij \lambdad imab \Ys jm 
ierr=t_einsum('E(i,j,a,b)+=L(i,m,a,b)*Y(j,m)', &
eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_oo,scale=MINUS_ONE)
ierr=t_einsum('E(i,j,a,b)+=L(j,m,a,b)*Y(i,m)', &
eta_prop_oovv_tensor,l2_tensor,Interm_t_tensor%prop_bar_oo,scale=ONE)
if(debug_local) print *,'LH eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=-\Pmm ab \ts em \lambdad ijae \Ys mb 
ierr=initialize_tensor_t(LR_Aux_vv)
ierr=t_einsum('H(e,b)+=T(e,m)*Y(m,b)',LR_Aux_vv,t1_tensor,Interm_t_tensor%prop_bar_ov)
ierr=t_einsum('E(i,j,a,b)+=H(e,b)*L(i,j,a,e)',eta_prop_oovv_tensor,LR_Aux_vv,l2_tensor,scale=MINUS_ONE)         
ierr=t_einsum('E(i,j,a,b)+=H(e,a)*L(i,j,b,e)',eta_prop_oovv_tensor,LR_Aux_vv,l2_tensor,scale=ONE)         
if(debug_local) print *,'LH eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=-\Pmm ij \ts em \lambdad imab\Ys je 
ierr=initialize_tensor_t(LR_Aux_oo)
ierr=t_einsum('H(j,m)+=T(e,m)*Y(j,e)',LR_Aux_oo,t1_tensor,Interm_t_tensor%prop_bar_ov)
ierr=t_einsum('E(i,j,a,b)+=H(j,m)*L(i,m,a,b)',eta_prop_oovv_tensor,LR_Aux_oo,l2_tensor,scale=MINUS_ONE)
ierr=t_einsum('E(i,j,a,b)+=H(i,m)*L(j,m,a,b)',eta_prop_oovv_tensor,LR_Aux_oo,l2_tensor,scale=ONE)
if(debug_local) print *,'LH eta_prop_ov_tensor ierr ', ierr

end subroutine get_eta_prop_oovv_prime



!----------------------------------------------Back up for use Hbar to approximate AZ------------------------------
!------------------------------------------------------------------------------------------------------------------

subroutine construct_AyTz(exa_input,Interm_t_tensor,tprop_tensor ,sigma_prop_tensor)
type(exacc_input), intent(in) :: exa_input
type(eti_rsp_interm_tens_t) :: Interm_t_tensor
type(eti_comm_tens_t)       :: tprop_tensor     
type(eti_comm_tens_t)       :: sigma_prop_tensor

integer(INTD), dimension(4)   :: vvoo_dims
integer(INTD), dimension(2)   :: vo_dims
integer                       :: ierr


vvoo_dims(1:2) = exa_input%nvir
vvoo_dims(3:4) = exa_input%nocc
vo_dims(1) = exa_input%nvir
vo_dims(2) = exa_input%nocc



ierr=construct_tensor_t(sigma_prop_tensor%t1,C8,vo_dims,init_val=ZERO)
! print *, "ierr1: ", ierr
ierr=construct_tensor_t(sigma_prop_tensor%t2,C8,vvoo_dims,init_val=ZERO)
! print *, "ierr2: ", ierr
call get_sigma_eomccsd_prop(Interm_t_tensor, tprop_tensor, sigma_prop_tensor)    

end subroutine construct_AyTz


!------------------------------------------------------------------------------------------------------------------

subroutine get_sigma_eomccsd_prop(interm_t_prop1, tprop_tensor2, sigma_prop_tensor12)

type(eti_rsp_interm_tens_t) :: interm_t_prop1
type(eti_comm_tens_t)            :: tprop_tensor2
type(eti_comm_tens_t)       :: sigma_prop_tensor12

integer :: ierr


!***************
!l1_sigma_prop *
!***************

!--------------------------------------------------------
!     term 1: Sxz(a,i) += Xbar(a,e)* tY(e,i) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,i)+=Xbar(a,e)*Tz(e,i)",sigma_prop_tensor12%t1,interm_t_prop1%prop_bar_vv,tprop_tensor2%t1)
! print *, "ierr3: ", ierr

!--------------------------------------------------------
!     term 2: Sxz(a,i) += -Xbar(m,i)* tY(a,m) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,i)+=Xbar(m,i)*Tz(a,m)", &
sigma_prop_tensor12%t1,interm_t_prop1%prop_bar_oo,tprop_tensor2%t1,scale=MINUS_ONE)
! print *, "ierr4: ", ierr

!--------------------------------------------------------
!     term 3: Sxz(a,i) += Xbar(m,e)* tY(e,a,m,i) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,i)+=Xbar(m,e)*Tz(e,a,m,i)", &
sigma_prop_tensor12%t1,interm_t_prop1%prop_bar_ov,tprop_tensor2%t2)
! print *, "ierr5: ", ierr


!***************
!l2_sigma_prop *
!***************

!--------------------------------------------------------
!     term 1: Sxz(a,b,i,j) += -P(ab) Q(m,b,i,j)* tY(a,m) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,b,i,j)+=Q(m,b,i,j)*Tz(a,m)", &
sigma_prop_tensor12%t2,interm_t_prop1%Q_ovoo,tprop_tensor2%t1,scale=MINUS_ONE)
! print *, "ierr6: ", ierr
ierr=t_einsum("Sxz(a,b,i,j)+=Q(m,a,i,j)*Tz(b,m)", &
sigma_prop_tensor12%t2,interm_t_prop1%Q_ovoo,tprop_tensor2%t1,scale=ONE)
! print *, "ierr7: ", ierr

!--------------------------------------------------------
!     term 2: Sxz(a,b,i,j) += P(ij) Q(a,b,j,e)* tY(e,i)   |
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,b,i,j)+=Q(a,b,e,j)*Tz(e,i)", &
sigma_prop_tensor12%t2,interm_t_prop1%Q_vvvo,tprop_tensor2%t1,scale=ONE)
! print *, "ierr8: ", ierr
ierr=t_einsum("Sxz(a,b,i,j)+=Q(a,b,e,i)*Tz(e,j)", &
sigma_prop_tensor12%t2,interm_t_prop1%Q_vvvo,tprop_tensor2%t1,scale=MINUS_ONE)
! print *, "ierr9: ", ierr

!--------------------------------------------------------
!     term 3: Sxz(a,b,i,j) += P(ab) Xbar(b,e)* tY(a,e,i,j)|
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,b,i,j)+=Xbar(b,e)*Tz(a,e,i,j)", &
sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_vv,tprop_tensor2%t2,scale=ONE)
! print *, "ierr10: ", ierr
ierr=t_einsum("Sxz(a,b,i,j)+=Xbar(a,e)*Tz(b,e,i,j)", &
sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_vv,tprop_tensor2%t2,scale=MINUS_ONE)
! print *, "ierr11: ", ierr

!--------------------------------------------------------
!     term 4: Sxz(a,b,i,j) += -P(ij) Xbar(m,j)* tY(a,b,i,m)|
!--------------------------------------------------------
ierr=t_einsum("Sxz(a,b,i,j)+=Xbar(m,j)*Tz(a,b,i,m)", &
sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_oo,tprop_tensor2%t2,scale=MINUS_ONE)
! print *, "ierr12: ", ierr
ierr=t_einsum("Sxz(a,b,i,j)+=Xbar(m,i)*Tz(a,b,j,m)", &
sigma_prop_tensor12%t2,interm_t_prop1%prop_bar_oo,tprop_tensor2%t2,scale=ONE)
! print *, "ierr13: ", ierr


end subroutine get_sigma_eomccsd_prop



!-----------------------------------------------------------------------------------------------------------------------------

subroutine get_sigma_eomccsd_bar_prop(interm_t_prop2, tprop_bar_tensor1, sigma_prop_tensor12)

type(eti_rsp_interm_tens_t) :: interm_t_prop2
type(eti_comm_tens_t)       :: tprop_bar_tensor1
type(eti_comm_tens_t)       :: sigma_prop_tensor12
! type(eti_intermediates_tens_t),intent(inout) :: interm_t

integer :: ierr


!***************
!l1_sigma_prop *
!***************

!--------------------------------------------------------
!     term 1: Sxz(i,a) += Xbar(a,e)* tlY(i,e) |
!--------------------------------------------------------
ierr=t_einsum("Sxy(i,a)+=Xbar(a,e)*Tz(i,e)",sigma_prop_tensor12%t1,interm_t_prop2%prop_bar_vv,tprop_bar_tensor1%t1)

!--------------------------------------------------------
!     term 2: Sxz(i,a) += -Xbar(m,i)* tlY(m,a) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(i,a)+=Xbar(m,i)*Tz(m,a)", &
sigma_prop_tensor12%t1,interm_t_prop2%prop_bar_oo,tprop_bar_tensor1%t1,scale=MINUS_ONE)

!--------------------------------------------------------
!     term 3: Sxz(i,a) += Q(e,f,a,m)* tlY(i,m,e,f) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(i,a)+=Q(e,f,a,m)*Tlz(i,m,e,f)", &
sigma_prop_tensor12%t1,interm_t_prop2%Q_vvvo,tprop_bar_tensor1%t2,scale=ONE_HALF)

!--------------------------------------------------------
!     term 3: Sxz(i,a) += -Q(i,e,m,n)* tlY(m,n,a,e) |
!--------------------------------------------------------
ierr=t_einsum("Sxz(i,a)+=Q(i,e,m,n)*Tlz(m,n,a,e)", &
sigma_prop_tensor12%t1,interm_t_prop2%Q_ovoo,tprop_bar_tensor1%t2,scale=MINUS_ONE_HALF)

!***************
!l2_sigma_prop *
!***************

!--------------------------------------------------------
!     term 2: S2(i,j,a,b) += P(a,b) L(i,j,a,e) * F(e,b) |
!--------------------------------------------------------

ierr=t_einsum("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t2,interm_t_prop2%prop_bar_vv)
ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t2,interm_t_prop2%prop_bar_vv)

!--------------------------------------------------------
!     term 3: S2(i,j,a,b) -= P(i,j) L(i,m,a,b) * F(j,m) |
!--------------------------------------------------------

ierr=t_einsum("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t2,interm_t_prop2%prop_bar_oo)
ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t2,interm_t_prop2%prop_bar_oo)
! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0

!------------------------------------------------------------
!     term 11: S2(i,j,a,b) += P(i,j) P(a,b) L(i,a) * F(j,b) |
!------------------------------------------------------------

ierr=t_einsum("S(i,j,a,b)+=L(i,a)*F(j,b)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t1,interm_t_prop2%prop_bar_ov)
ierr=t_einsum("S(i,j,a,b)+=L(j,a)*F(i,b)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t1,interm_t_prop2%prop_bar_ov,scale=MINUS_ONE)
ierr=t_einsum("S(i,j,a,b)+=L(i,b)*F(j,a)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t1,interm_t_prop2%prop_bar_ov,scale=MINUS_ONE)
ierr=t_einsum("S(i,j,a,b)+=L(j,b)*F(i,a)",sigma_prop_tensor12%t2,&
        tprop_bar_tensor1%t1,interm_t_prop2%prop_bar_ov)
end subroutine get_sigma_eomccsd_bar_prop

end module
