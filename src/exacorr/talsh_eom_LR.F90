module talsh_eom_LR
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), IMAGINARY=(0.0d0,1.0d0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0), &
                                ONE_EIGHTH_C=(0.125D0,0.D0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public get_eom_lr
    public get_asym_lr
    public get_tF

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine get_eom_lr(exa_input,interm_t,int_t,&
                Interm_t_tensors_A,Interm_t_tensors_B,&
                tprop_tensors, tprop_tensors_min_w, i,j,&
                l1_tensor,l2_tensor, eom_lr, CC_CC)
!   This subroutine drives the CC-CI and CC-CC linear response property calculation
!   Written by Xiang Yuan Janunary 2022, Lille, France

            use talsh_sigma_eomccsd

            type(exacc_input), intent(in) :: exa_input
            type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
            type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:), Interm_t_tensors_B(:)
            type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
            type(talsh_intg_tens_t),intent(inout)          :: int_t 
            type(talsh_tens_t), intent(inout)           :: l2_tensor
            type(talsh_tens_t), intent(inout)           :: l1_tensor
            integer                                     :: i,j 


            type(talsh_rsp_interm_tens_t) :: interm_t_prop1, interm_t_prop2
            type(talsh_comm_tens_t)            :: tprop_tensor1, tprop_tensor2
            type(talsh_comm_tens_t)            :: tprop_tensor1_minus_fre, tprop_tensor2_minus_fre


            type(talsh_comm_tens_t)            :: B1, B2   !Ft^x
            type(talsh_comm_tens_t)            :: B1_minus_fre, B2_minus_fre   !Ft^x
            integer            :: nocc, nvir
            integer(INTD), dimension(4)   :: vvoo_dims
            integer(INTD), dimension(2)   :: vo_dims
            complex(kind=8)               :: linrsp, linrsp_con

            integer                       :: ierr

            complex(kind=8)               :: phase_factor

            complex(kind=8)               :: eom_lr
            integer                       :: T
            real(kind=8)                  :: tem_real, tem_imag

            logical :: CC_CC
            logical :: debug_get_lrresp = .false.
            type(talsh_tens_t)    :: total_t_t1, total_t_t2


            eom_lr = dcmplx(0.0,0.0)

            vvoo_dims(3:4) = exa_input%nocc
            vvoo_dims(1:2) = exa_input%nvir
            vo_dims(2) = exa_input%nocc
            vo_dims(1) = exa_input%nvir
            nocc = exa_input%nocc
            nvir = exa_input%nvir


            if (CC_CC) then
                call get_tF(B1%t1,B1%t2,Tprop_tensors(1)%t1,Tprop_tensors(1)%t2,l1_tensor,l2_tensor,&
                        interm_t,int_t,nocc,nvir)
                call get_tF(B2%t1,B2%t2,Tprop_tensors(2)%t1,Tprop_tensors(2)%t2,l1_tensor,l2_tensor,&
                        interm_t,int_t,nocc,nvir)

                if (debug_get_lrresp) then 
                    call print_date("Finish CC-CC tF terms (+w)")
                end if  

                call get_tF(B1_minus_fre%t1,B1_minus_fre%t2,Tprop_tensors_min_w(1)%t1,Tprop_tensors_min_w(1)%t2,& 
                        l1_tensor,l2_tensor,interm_t,int_t,nocc,nvir)
                call get_tF(B2_minus_fre%t1,B2_minus_fre%t2,Tprop_tensors_min_w(2)%t1,Tprop_tensors_min_w(2)%t2,& 
                        l1_tensor,l2_tensor,interm_t,int_t,nocc,nvir)

                if (debug_get_lrresp) then 
                    call print_date("Finish CC-CC tF terms (-w)")
                end if  

                if (debug_get_lrresp) then
                    call print_tensor(B1%t1,1.0d-18,"T1 part of tF")
                    call print_tensor(B1%t2,1.0d-18,"T2 part of tF")
                end if 
            end if 


            if (debug_get_lrresp) then
                if (CC_CC) then 
                    call print_date('Evaluate CC-CC LR (+w)')
                else
                    call print_date('Evaluate CC-CI LR (+w)')
                end if 
            end if 


            call get_linrspf(interm_t, Interm_t_tensors_A(i), Interm_t_tensors_B(j), &
                    Tprop_tensors(1), Tprop_tensors(2), &
                    l1_tensor, l2_tensor, linrsp, B1, B2, CC_CC)

            
            if (debug_get_lrresp) then
                if (CC_CC) then 
                    call print_date('Evaluate CC-CC LR (-w)')
                else
                    call print_date('Evaluate CC-CI LR (-w)')
                end if 
            end if 

            call get_linrspf(interm_t, Interm_t_tensors_A(i), Interm_t_tensors_B(j), &
                    Tprop_tensors_min_w(1), Tprop_tensors_min_w(2),&
                    l1_tensor, l2_tensor, linrsp_con, B1_minus_fre, B2_minus_fre, CC_CC)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASPG+LV: Note concerning the addition of a phase to the response functions below 
!
! The 1-body integrals over time-antisymmetric operators are already in the AO basis 
! of DIRAC scaled by imaginary i  so that they can be treated as time-symmetric (real) 
! operators. Such an auxiliary phase factor takes over to the MO basis in which we solve 
! the CC response equations:
!
! RHS (for linear and quadratic response)
!  (\bar{H} - \omega_y I)t^Y(\omega_y) = - \xi^Y  (RHS, for linear and quadratic response)
! LHS (for quadratic response and asymmetric linear response)  
! \bar{t}^X(\omega_x)(\bar{H} - \omega_x I) = - \eta^X - Ft^X(\omega_x) 
!
! where \xi^Y and \eta^X involve the one-body operators (see 10.1021/acs.jctc.3c00812 and 10.1021/acs.jctc.3c01011).
!
! These auxiliary phase factors are therefore present when assembling the response functions, e.g.
!
!   LR, symmetric formulation
!   <<X; Y>>_{\omega_y} = 1/2 C^{\pm\omega} P(X,Y) [ \eta_Y + 1/2 Ft^X(\omega_x) ] t^Y(omega_y)
!
!   LR, asymmetric formulation 
!   <<X; Y>>_{\omega_y} = 1/2 [\eta^X t^Y(\omega_y) + \bar{t}^Y(\omega_y)\xi^X]
!   
!   QR, symmetric formulation
!   <<X; Y, Z>>_{\omega} = 1/2 C^{\pm\omega} P(X,Y,Z) [[ 1/2 F^Y + 1/6 Gt^Y]t^Y + \bar{t}^X [ A^Y + 1/2 Bt^Y]] t^Z 
!    
!   etc (see aforementioned papers for the QR-EOM expression)
!
!   where
!    - \omega_x = -(\omega_y + \omega_z) or \omega_x = -\omega_y,
!    - P a permutation operation exchanging the role of X and Y,
!    - C a symmetrization operator with respect to simultaneous complex conjugation and inversion 
!      of the signs of the frequencies,
!      C^{\pm\omega} = f^{XY}(\omega_x, \omega_y) + [f^{XY}(-\omega_x, -\omega_y)]^*
!
! and we need to correct for this to obtain the final result.
!
! We can do so by multiplying the terms, for example f^{XY}(\omega_x, \omega_y) and f^{XY}(-\omega_x, -\omega_y) 
! for linear response, by a phase factor G
!
! G = i^{N_a},
!
! where N_a is the number of time-antisymmetric operators among X, Y, Z,...; and this should be made BEFORE 
! the action of C^{\pm\omega} if that is applied.
!
! With this, for example :
! - in the case of the elements of the dipole polarizability tensor ((<<x; y>>)), we have two time-symmetric 
!   operators (x,y) so N_a = 0, there is no change
! - in the case of the optical rotation (<<x; \hat{m}_y>>), we have one time-symmetric (x) and one 
!   time-antisymmetric (the magnetic dipole component \hat{m}_y) operator, so N_a = 1 and we have as result G = i; by 
!   applying it the real and imaginary parts of the response function are interchanged
! - in the case of spin-spin couplings between centers M and N (<<\hat{h}_M;\hat{h}_M>>), we have two 
!   time-antisymmetric operators, so N_a = 2 and G = i^2 = -1, as a result we change the sign of f^{XY}
!
! and so on for the possible combinations in the case of quadratic response functions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            print *, "Number of time-antisymmetric operators in response function",exa_input%N_Time_anti
            if (exa_input%N_Time_anti.ge.1) &
            print *, "Adding phase to response function components (due to how DIRAC scales time-antisymmetric operators)"

            phase_factor = IMAGINARY**(exa_input%N_Time_anti)

            linrsp = linrsp*phase_factor
            linrsp_con = linrsp_con*phase_factor

            print *, "linear response: +w", linrsp
            print *, "linear response: -w", linrsp_con

            linrsp_con=conjg(linrsp_con)
            eom_lr=(linrsp+linrsp_con)*0.5

            print *, "LR after symmetrization ",eom_lr


end subroutine get_eom_lr

!------------------------------------------------------------------------------------------------------------------
subroutine get_asym_lr(exa_input,interm_t, int_t, Interm_t_tensor_A, & 
                       Tprop_tensors, Tprop_tensors_min_w, &
                       Tprop_bar_tensors, Tprop_bar_tensors_min_w, &
                       l1_tensor,l2_tensor,value_lr)
!   This subroutine drives the CC-CI and CC-CC linear response property calculation
!   in the asymmetric formulation
!   Written by Andre Gomes

    use talsh_sigma_eomccsd

    type(exacc_input), intent(in) :: exa_input
    type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !  t1, t2 tensor
    type(talsh_intg_tens_t),intent(inout)       :: int_t 
    type(talsh_tens_t), intent(inout)           :: l2_tensor
    type(talsh_tens_t), intent(inout)           :: l1_tensor
    complex(kind=8), intent(inout)              :: value_lr
    ! \eta and \xi (or epsilon) terms for operators A
    type(talsh_rsp_interm_tens_t), intent(inout) :: Interm_t_tensor_A
    ! RHS response vectors 
    type(talsh_comm_tens_t), allocatable, intent(inout) :: Tprop_tensors(:)
    type(talsh_comm_tens_t), allocatable, intent(inout) :: Tprop_tensors_min_w(:)
    ! LHS response vectors
    type(talsh_comm_tens_t), allocatable, intent(inout) :: Tprop_bar_tensors(:) 
    type(talsh_comm_tens_t), allocatable, intent(inout) :: Tprop_bar_tensors_min_w(:) 

    integer(C_INT)              :: ierr, one_dims(1)
    complex(8), pointer         :: result_tens(:)
    type(C_PTR)                 :: body_p
    type(talsh_tens_t)          :: result_tensor

    complex(kind=8)             :: linrsp, linrsp_con, phase_factor

    logical :: debug_local = .true.
    
    ! create and zero out result_tensor to accumulate +\omega contribution
    one_dims(1) = 1
    ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    if (debug_local) print *, 'ierr 1.a',ierr
    ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    if (debug_local) print *, 'ierr 1.b',ierr
    call c_f_pointer(body_p,result_tens,one_dims)

    ierr=talsh_tensor_init(result_tensor, ZERO)
    if (debug_local) print *, 'ierr 1.c',ierr, result_tens(1)

    ! build \eta^A t^B(+omega)
    
    ierr=talsh_tensor_contract("E()+=EtaA(i,a)*TB(a,i)", result_tensor, &
            Interm_t_tensor_A%eta_prop_ov,Tprop_tensors(2)%t1)
    
    if (debug_local) then
    print *, 'ierr 2.a',ierr,result_tens(1)
    print *, '        eta1 norm',tensor_norm2(Interm_t_tensor_A%eta_prop_ov)
    print *, '        tp1  norm',tensor_norm2(Tprop_tensors(2)%t1)
    end if

    ierr=talsh_tensor_contract("E()+=EtaA(i,j,a,b)*TB(a,b,i,j)", result_tensor, &
            Interm_t_tensor_A%eta_prop_oovv,Tprop_tensors(2)%t2, &
            scale=ONE_QUARTER_C)
   
    if (debug_local) then
    print *, 'ierr 2.b',ierr,result_tens(1)
    print *, '        eta2 norm',tensor_norm2(Interm_t_tensor_A%eta_prop_oovv)
    print *, '        tp2  norm',tensor_norm2(Tprop_tensors(2)%t2)
    end if

    ! build \bar(t)^B(+omega)\epsilon^A 
    
    ierr=talsh_tensor_contract("E()+=barTB(i,a)*XiA(a,i)", result_tensor, &
            Tprop_bar_tensors(2)%t1,Interm_t_tensor_A%epsilon_prop_vo)
    
    if (debug_local) then
    print *, 'ierr 3.a',ierr,result_tens(1)
    print *, '        tbp2 norm',tensor_norm2(Tprop_bar_tensors(2)%t1)
    print *, '        xi1  norm',tensor_norm2(Interm_t_tensor_A%epsilon_prop_vo)
    end if

    ierr=talsh_tensor_contract("E()+=barTB(i,j,a,b)*XiA(a,b,i,j)", result_tensor, &
            Tprop_bar_tensors(2)%t2,Interm_t_tensor_A%epsilon_prop_vvoo, &
            scale=ONE_QUARTER_C)
    
    if (debug_local) then
    print *, 'ierr 3.b',ierr,result_tens(1)
    print *, '        tbp2 norm',tensor_norm2(Tprop_bar_tensors(2)%t2)
    print *, '        xi2  norm',tensor_norm2(Interm_t_tensor_A%epsilon_prop_vvoo)
    end if

    ! save result
    linrsp =  result_tens(1)
    if (debug_local) print *, 'linrsp',linrsp

    ! zero out result_tensor to accumulate -\omega contribution

    ierr=talsh_tensor_init(result_tensor, ZERO)
    
    if (debug_local) print *, 'ierr 4.a',ierr,result_tens(1)

    ! build \eta^A t^B(-omega)

    ierr=talsh_tensor_contract("E()+=EtaA(i,a)*TB(a,i)", result_tensor, &
            Interm_t_tensor_A%eta_prop_ov,Tprop_tensors_min_w(2)%t1)
    
    if (debug_local) print *, 'ierr 5.a',ierr,result_tens(1)

    ierr=talsh_tensor_contract("E()+=EtaA(i,j,a,b)*TB(a,b,i,j)", result_tensor, &
            Interm_t_tensor_A%eta_prop_oovv,Tprop_tensors_min_w(2)%t2, &
            scale=ONE_QUARTER_C)
    
    if (debug_local) print *, 'ierr 5.b',ierr,result_tens(1)

    ! build \bar(t)^B(-omega)\epsilon^A 

    ierr=talsh_tensor_contract("E()+=barTB(i,a)*XiA(a,i)", result_tensor, &
            Tprop_bar_tensors_min_w(2)%t1,Interm_t_tensor_A%epsilon_prop_vo)
    
    if (debug_local) print *, 'ierr 6.a',ierr,result_tens(1)

    ierr=talsh_tensor_contract("E()+=barTB(i,j,a,b)*XiA(a,b,i,j)", result_tensor, &
            Tprop_bar_tensors_min_w(2)%t2,Interm_t_tensor_A%epsilon_prop_vvoo, &
            scale=ONE_QUARTER_C)
    
    if (debug_local) print *, 'ierr 6.b',ierr,result_tens(1)

    ! save result
    linrsp_con = result_tens(1)
    
    if (debug_local) print *, 'linrsp_con',linrsp_con

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASPG+LV: Note concerning the addition of a phase to the response functions below 
!
! The 1-body integrals over time-antisymmetric operators are already in the AO basis 
! of DIRAC scaled by imaginary i  so that they can be treated as time-symmetric (real) 
! operators. Such an auxiliary phase factor takes over to the MO basis in which we solve 
! the CC response equations:
!
! RHS (for linear and quadratic response)
!  (\bar{H} - \omega_y I)t^Y(\omega_y) = - \xi^Y  (RHS, for linear and quadratic response)
! LHS (for quadratic response and asymmetric linear response)  
! \bar{t}^X(\omega_x)(\bar{H} - \omega_x I) = - \eta^X - Ft^X(\omega_x) 
!
! where \xi^Y and \eta^X involve the one-body operators (see 10.1021/acs.jctc.3c00812 and 10.1021/acs.jctc.3c01011).
!
! These auxiliary phase factors are therefore present when assembling the response functions, e.g.
!
!   LR, symmetric formulation
!   <<X; Y>>_{\omega_y} = 1/2 C^{\pm\omega} P(X,Y) [ \eta_Y + 1/2 Ft^X(\omega_x) ] t^Y(omega_y)
!
!   LR, asymmetric formulation 
!   <<X; Y>>_{\omega_y} = 1/2 C^{\pm\omega} [\eta^X t^Y(\omega_y) + \bar{t}^Y(\omega_y)\xi^X]
!   
!   QR, symmetric formulation
!   <<X; Y, Z>>_{\omega} = 1/2 C^{\pm\omega} P(X,Y,Z) [[ 1/2 F^Y + 1/6 Gt^Y]t^Y + \bar{t}^X [ A^Y + 1/2 Bt^Y]] t^Z 
!    
!   etc (see aforementioned papers for the QR-EOM expression)
!
!   where
!    - \omega_x = -(\omega_y + \omega_z) or \omega_x = -\omega_y,
!    - P a permutation operation exchanging the role of X and Y,
!    - C a symmetrization operator with respect to simultaneous complex conjugation and inversion 
!      of the signs of the frequencies,
!      C^{\pm\omega} = f^{XY}(\omega_x, \omega_y) + [f^{XY}(-\omega_x, -\omega_y)]^*
!
! and we need to correct for this to obtain the final result.
!
! We can do so by multiplying the terms, for example f^{XY}(\omega_x, \omega_y) and f^{XY}(-\omega_x, -\omega_y) 
! for linear response, by a phase factor G
!
! G = i^{N_a},
!
! where N_a is the number of time-antisymmetric operators among X, Y, Z,...; and this should be made BEFORE 
! the action of C^{\pm\omega} if that is applied.
!
! With this, for example :
! - in the case of the elements of the dipole polarizability tensor ((<<x; y>>)), we have two time-symmetric 
!   operators (x,y) so N_a = 0, there is no change
! - in the case of the optical rotation (<<x; \hat{m}_y>>), we have one time-symmetric (x) and one 
!   time-antisymmetric (the magnetic dipole component \hat{m}_y) operator, so N_a = 1 and we have as result G = i; by 
!   applying it the real and imaginary parts of the response function are interchanged
! - in the case of spin-spin couplings between centers M and N (<<\hat{h}_M;\hat{h}_M>>), we have two 
!   time-antisymmetric operators, so N_a = 2 and G = i^2 = -1, as a result we change the sign of f^{XY}
!
! and so on for the possible combinations in the case of quadratic response functions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            print *, "Number of time-antisymmetric operators in response function",exa_input%N_Time_anti
            if (exa_input%N_Time_anti.ge.1) &
            print *, "Adding phase to response function components (due to how DIRAC scales time-antisymmetric operators)"

            phase_factor = IMAGINARY**(exa_input%N_Time_anti)

            linrsp = linrsp*phase_factor
            linrsp_con = linrsp_con*phase_factor

            print *, "linear response: +w", linrsp
            print *, "linear response: -w", linrsp_con

            linrsp_con=conjg(linrsp_con)
            value_lr=(linrsp+linrsp_con)*0.5

            print *, "LR after symmetrization ",value_lr

! cleanup
            ierr=talsh_tensor_destruct(result_tensor)
            print *, 'ierr 7.a',ierr


end subroutine get_asym_lr

!------------------------------------------------------------------------------------------------------------------

subroutine get_linrspf(interm_t, interm_t_prop1, interm_t_prop2, &
                tprop_tensor1, tprop_tensor2, &
                l1_tensor, l2_tensor, linrsp, B1,B2,CC_CC)

        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_tens_t), intent(inout)           :: l2_tensor
        type(talsh_tens_t), intent(inout)           :: l1_tensor
        type(talsh_rsp_interm_tens_t) :: interm_t_prop1, interm_t_prop2
        type(talsh_comm_tens_t)            :: tprop_tensor1, tprop_tensor2
        ! type(talsh_comm_tens_t)       :: sigma_prop_tensor12, sigma_prop_tensor21
        complex(kind=8)               :: linrsp
        type(talsh_comm_tens_t),optional   :: B1, B2
        logical :: CC_CC
        complex(kind=8)               :: linrsp_iter

!--------------------------------------------------------------------------------------
!      CC-CC  << X;Y >>ωY = P(X,Y) [1/2tX(FtY) + etaX*TY]                             |
!      CC-CI  << X;Y >>ωY = P(X,Y) [-(lambda*tX)*(lambda*epsilonY) + (eom)etaX*TY]    |                          
!--------------------------------------------------------------------------------------
            call get_linrsp_driver(linrsp,tprop_tensor1, tprop_tensor2, &
                    l1_tensor,l2_tensor, &
                    interm_t,interm_t_prop1,B1,CC_CC)

            call get_linrsp_driver(linrsp_iter,tprop_tensor2, tprop_tensor1, &
                    l1_tensor,l2_tensor, &
                    interm_t,interm_t_prop2,B2,CC_CC)
            
            linrsp=linrsp+linrsp_iter


end subroutine get_linrspf

!------------------------------------------------------------------------------------------------------------------

subroutine get_linrsp_driver(tem_linrsp_iter,tprop_tensor1, tprop_tensor2, &
                l1_tensor,l2_tensor,&
                interm_t,interm_t_prop1,B,CC_CC)

    type(talsh_intermediates_tens_t),intent(inout) :: interm_t 
    type(talsh_rsp_interm_tens_t),intent(inout)    :: interm_t_prop1
    type(talsh_tens_t), intent(inout)              :: l2_tensor
    type(talsh_tens_t), intent(inout)              :: l1_tensor
    complex(kind=8)                                :: tem_linrsp_iter

    type(talsh_comm_tens_t)                        :: tprop_tensor1, tprop_tensor2
        
    ! type(talsh_comm_tens_t)       :: sigma_prop_tensor12
    type(talsh_comm_tens_t),optional   :: B
    logical :: CC_CC
    logical :: debug_local = .false.


    !  Point connected with tensors 
    type(talsh_tens_t)          :: result_tensor
    integer(C_INT)              :: ierr, etaX_rank, l2_ranks, epsilon1_rank
    complex(8), pointer         :: result_tens(:)
    complex(8), pointer         :: l2_tens(:,:,:,:), epsilon1_tens(:,:),etaX1_tens(:,:)
    integer                     :: initialize_talsh

    type(C_PTR):: body_p, body_p_1, body_p_2, body_p_3
    integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
    integer(C_INT):: host_arg_max

    integer(C_INT)              :: one_dims(1), l2_dims(4), epsilon1_dims(2),etaX1_dims(2)

    complex(8)                  :: term1, term2, term3
    integer                     :: i,j,a,c, N_occ, N_vir
    complex(8)                  :: contrac_l2


    type(talsh_tens_t)          :: eta_eom_test1, eta_eom_test2
    integer(INTD)  ::  a_rank
    integer(INTD)  ::  a_dims2(1:2)
    integer(INTD)  ::  a_dims4(1:4)


    one_dims(1) = 1
    ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,result_tens,one_dims)

    ! a_rank=talsh_tensor_rank(l1_tensor)
    ! ierr = talsh_tensor_dimensions(l1_tensor,a_rank,a_dims2)
    ! ierr=talsh_tensor_construct(eta_eom_test1,C8,a_dims2(1:2),init_val=ZERO)

    ! a_rank=talsh_tensor_rank(l2_tensor)
    ! ierr = talsh_tensor_dimensions(l2_tensor,a_rank,a_dims4)
    ! ierr=talsh_tensor_construct(eta_eom_test2,C8,a_dims4(1:4),init_val=ZERO)

 

    if (CC_CC) then
!--------------------------------------------------------------------------------------
!       Term1 = 1/2 BX(nu)*tY(nu)  in CC-CC                                          |
!--------------------------------------------------------------------------------------
        ierr=talsh_tensor_contract("E()+=B(i,a)*Ty(a,i)",& 
            result_tensor,B%t1,tprop_tensor2%t1,scale=ONE_HALF)
        term1=result_tens(1)

        if (debug_local) then
            print *, "Contribution of tFt single part(CC-CC): ", term1
        end if 

        ierr=talsh_tensor_contract("E()+=B(i,j,a,b)*Ty(a,b,i,j)", &
            result_tensor,B%t2,tprop_tensor2%t2,scale=ONE_EIGHTH_C)
        term1=result_tens(1)

        if (debug_local) then
            print *, "Contribution of total tFt(CC-CC): ", term1
        end if 

        term2 = dcmplx(1.0,0.0)

    else 
!--------------------------------------------------------------------------------------
!       Term1 -= t_bar0(nu)*tY(nu)                                                    |
!--------------------------------------------------------------------------------------
 
        ierr=talsh_tensor_init(result_tensor)

        ierr=talsh_tensor_contract("E()+=Ty(a,i)*Tl(i,a)",result_tensor,tprop_tensor2%t1,l1_tensor)

        ierr=talsh_tensor_contract("E()+=Ty(a,b,i,j)*Tl(i,j,a,b)", &
                result_tensor,tprop_tensor2%t2,l2_tensor,scale=ONE_QUARTER_C)

        term1=result_tens(1)

        if (debug_local) then 
            print *, "Contribution of lambda*tY total(CC-CI): ", term1
        end if 

!--------------------------------------------------------------------------------------
!       Term2 += t_bar0(mu)*episilonX(mu)                                             |
!--------------------------------------------------------------------------------------

        ierr=talsh_tensor_init(result_tensor)
        ierr=talsh_tensor_contract("E()+=Tl(i,a)*Ex(a,i)",& 
            result_tensor,l1_tensor,interm_t_prop1%epsilon_prop_vo,scale=MINUS_ONE)
        term2=result_tens(1)

        if (debug_local) then 
            print *, "Contribution of lambda*episilonX single part(CC-CI): ", term2
        end if 

        ierr=talsh_tensor_contract("E()+=Tl(i,j,a,b)*Ex(a,b,i,j)", &
            result_tensor,l2_tensor,interm_t_prop1%epsilon_prop_vvoo,scale=MINUS_ONE_QUARTER_C)
        term2=result_tens(1)


        if (debug_local) then 
            print *, "Contribution of lambda*episilonX total(CC-CI): ", term2
        end if 


    end if

    ierr=talsh_tensor_init(result_tensor)
    if (CC_CC) then 
        ierr=talsh_tensor_contract("E()+=nX(i,a)*Ty(a,i)",result_tensor,interm_t_prop1%eta_prop_ov,tprop_tensor2%t1)
    else 
        ierr=talsh_tensor_contract("E()+=nX(i,a)*Ty(a,i)",result_tensor,interm_t_prop1%eta_prop_ov_eom,tprop_tensor2%t1)
    end if 

    term3=result_tens(1)
    if (debug_local) then
        print *, "single part of etax*tY", term3
    end if 
    ierr=talsh_tensor_contract("E()+=nX(i,j,a,b)*Ty(a,b,i,j)",result_tensor,interm_t_prop1%eta_prop_oovv, &
                            tprop_tensor2%t2,scale=ONE_QUARTER_C)
    term3=result_tens(1)

    if (debug_local) then
        print *, "Total Contribution of etaX*tY: ", term3
    end if 

    tem_linrsp_iter=term1*term2+term3
    ! print *, "total:", tem_linrsp_iter


    ! ierr=talsh_tensor_add("ETAEOM(i,a)+=ETA(i,a)",eta_eom_test1,interm_t_prop1%eta_prop_ov_eom,scale=ONE)
    ! ierr=talsh_tensor_add("ETAEOM(i,a)+=L1(i,a)",eta_eom_test1,l1_tensor,scale=term2)

    ! ierr=talsh_tensor_add("ETAEOM(i,j,a,b)+=ETA(i,j,a,b)",eta_eom_test2,interm_t_prop1%eta_prop_oovv,scale=ONE)
    ! ierr=talsh_tensor_add("ETAEOM(i,j,a,b)+=L2(i,j,a,b)",eta_eom_test2,l2_tensor,scale=term2)

    ! call print_tensor(eta_eom_test1,1.0d-18,"T1 part of ETAEOM")
    ! call print_tensor(eta_eom_test2,1.0d-18,"T2 part of ETAEOM")

    ! ierr=talsh_tensor_destruct(eta_eom_test1)
    ! ierr=talsh_tensor_destruct(eta_eom_test2)
    
    ierr=talsh_tensor_destruct(result_tensor)

end subroutine 



!-------------------------------------------------------------------------------------------------------------------
!! LH start Ft^x

subroutine get_tF(B1,B2,tX1,tX2,l1_tensor,l2_tensor,interm_t,int_t,nocc,nvir, is_qrcc)

    type(talsh_tens_t) :: B1, B2
    type(talsh_tens_t), intent(inout) :: tX1, tX2
    type(talsh_tens_t), intent(inout) :: l1_tensor,l2_tensor
    type(talsh_intermediates_tens_t),  intent(inout) :: interm_t 
    type(talsh_intg_tens_t), intent(inout) :: int_t
    logical, optional, intent(in) :: is_qrcc
    !Tensor dimensions
    integer, intent(in) :: nocc, nvir
    integer(INTD), dimension(2) :: oo_dims, vv_dims, ov_dims, vo_dims
    integer(INTD), dimension(4) :: vvoo_dims, ovoo_dims, ovov_dims
    integer(INTD), dimension(4) :: vovo_dims, vovv_dims, oovo_dims
    integer(INTD), dimension(4) :: vooo_dims, ooov_dims, oovv_dims
    integer(INTD), dimension(4) :: oooo_dims, vvvv_dims, ovvo_dims
    integer(C_INT)              :: one_dims(1)

    integer(INTD)                   :: ierr
    logical :: debug_local = .false.
    logical :: add_U1 


    type(talsh_tens_t) :: Aux_oo, Aux_vv, Aux_vo, Aux_ov
    type(talsh_tens_t) :: Aux_oovv, Aux_ooov, Aux_vooo
    type(talsh_tens_t) :: Aux_oooo, Aux_vvvv, Aux_ovov
    type(talsh_tens_t) :: Aux_vovv, Aux_vovo, Aux_oovo
    type(talsh_tens_t) :: Aux_ovvo
    type(talsh_tens_t) :: Aux_one

    oo_dims(1:2) = nocc
    vv_dims(1:2) = nvir
    
    ov_dims(1)  = nocc
    ov_dims(2)  = nvir

    vo_dims(1)  = nvir
    vo_dims(2)  = nocc

    vvvv_dims(1:4) = nvir
    oooo_dims(1:4) = nocc

    vooo_dims(1:4) = nocc
    vooo_dims(1) = nvir        

    ooov_dims(1:4) =nocc
    ooov_dims(4) = nvir

    oovv_dims(1:2) = nocc
    oovv_dims(3:4) = nvir

    ovov_dims(1:4) = nocc
    ovov_dims(2)   = nvir
    ovov_dims(4)   = nvir

    vovo_dims(1:4) = nvir
    vovo_dims(2) = nocc
    vovo_dims(4) = nocc

    vovv_dims(1:4) = nvir
    vovv_dims(2) = nocc
    
    oovo_dims(1:4) = nocc
    oovo_dims(3)   = nvir

    ovvo_dims(1:4) = nocc
    ovvo_dims(2:3) = nvir

    one_dims(1) = 1

    ierr=talsh_tensor_construct(B1,C8,ov_dims,init_val=ZERO)        
    ierr=talsh_tensor_construct(B2,C8,oovv_dims,init_val=ZERO)        

    ierr=talsh_tensor_construct(Aux_oo,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_vv,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_vo,C8,vo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_ov,C8,ov_dims,init_val=ZERO)
    
    ierr=talsh_tensor_construct(Aux_vvvv,C8,vvvv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_oooo,C8,oooo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_oovv,C8,oovv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_vooo,C8,vooo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_ooov,C8,ooov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_vovo,C8,vovo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_ovov,C8,ovov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_vovv,C8,vovv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_oovv,C8,oovv_dims,init_val=ZERO)

    ierr=talsh_tensor_construct(Aux_oovo,C8,oovo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Aux_ovvo,C8,ovvo_dims,init_val=ZERO)

    ierr=talsh_tensor_construct(Aux_one,C8,one_dims(1:0),init_val=ZERO)

!   The U1 diagram should be included if we are using this routine to
!   build (Ft^x) in linear response. In the case we are using it to 
!   build (t^y F^x) in CC quadratic response, we omit it.
    add_U1 = .true.
    if (present(is_qrcc)) then
       if (is_qrcc) add_U1 = .false.
    end if

    if (add_U1) then
! l=0  ; t=1  ; h=-2
       ! U1 = + (\ts ck ) \ts ai \W kica       -> V oovv
       ierr=talsh_tensor_contract("B(k,c)+=T(a,i)*W(k,i,c,a)",B1,tX1,int_t%oovv)
       if (debug_local) print*, "ierr U1", ierr
    else
       if (debug_local) print*, "not adding U1 to tF"
    end if

! l=-1 ; t=1 ; h=-1
    ! c7 = - (\ts ck) \ts ai \lambdas ic \F ka
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,i)+=T(a,i)*F(k,a)",Aux_oo,tX1,interm_t%Fbar_ov)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,i)*L(i,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE)        
    if (debug_local) print*, "ierr c7", ierr

    ! c8 = - (\ts ck) \ts ai \lambdas ka \F ic        
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,i)+=L(k,a)*T(a,i)",Aux_oo,l1_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,i)*F(i,c)",B1,Aux_oo,interm_t%Fbar_ov,scale=MINUS_ONE)
    if (debug_local) print*, "ierr c8", ierr

    ! g2 = + (\ts ck) \ts ai \lambdas ck \W eica (vovv)
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(e,c)+=T(a,i)*W(e,i,c,a)",Aux_vv,tX1,interm_t%w_vovv)
    if (debug_local) print*, "ierr g2_1st", ierr
    ierr=talsh_tensor_contract("B(k,c)+=A(e,c)*L(k,e)",B1,Aux_vv,l1_tensor)
    if (debug_local) print*, "ierr g2_2nd", ierr

    ! k2 = - (\ts ck) \ts ai \lambdas mc \W kima (ooov)
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,m)+=T(a,i)*W(k,i,m,a)",Aux_oo,tX1,interm_t%w_ooov)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,m)*L(m,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE)
    if (debug_local) print*, "ierr k2", ierr

    ! k3 = - (\ts ck) \ts ai \lambdas ma \W ikmc (ooov)
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(m,i)+=L(m,a)*T(a,i)",Aux_oo,l1_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,c)+=A(m,i)*W(i,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(m,i)*W(i,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE_HALF)  !XY
    if (debug_local) print*, "ierr k3", ierr

    ! g3 = + (\ts ck) \ts ai \lambda ie \W ekac (vovv)
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(a,e)+=L(i,e)*T(a,i)",Aux_vv,l1_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,c)+=A(a,e)*W(e,k,a,c)",B1,Aux_vv,interm_t%w_vovv)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(a,e)*W(e,k,a,c)",B1,Aux_vv,interm_t%w_vovv,scale=ONE_HALF)    !XY
    if (debug_local) print*, "ierr g3", ierr


! l=-1; t=2 ; h=-2
    ! l5 = - one_half (\ts ck) \td abij \lambdas ic \W kjab (oovv) -> V oovv
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,i)+=T(a,b,i,j)*W(k,j,a,b)",Aux_oo,tX2,int_t%oovv)
    ! ierr=talsh_tensor_contract("B(c,k)+=A(k,i)*L(i,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,i)*L(i,c)",B1,Aux_oo,l1_tensor,scale=MINUS_ONE_HALF)   !XY
    if (debug_local) print*, "ierr l5_2", ierr

    ! l11 = - one_half (\ts ck) \td abij \lambdas ka \W ijcb (oovv) -> V oovv
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(a,c)+=T(a,b,i,j)*W(i,j,c,b)",Aux_vv,tX2,int_t%oovv)
    ! ierr=talsh_tensor_contract("B(c,k)+=A(a,c)*L(k,a)",B1,Aux_vv,l1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("B(k,c)+=A(a,c)*L(k,a)",B1,Aux_vv,l1_tensor,scale=MINUS_ONE_HALF)  !XY
    if (debug_local) print*, "ierr l11", ierr

    ! l18 = + (\ts ck) \td abij \lambdas ia \W jkbc (oovv) -> V oovv
    ierr=talsh_tensor_init(Aux_vo)
    ierr=talsh_tensor_contract("A(b,j)+=T(a,b,i,j)*L(i,a)",Aux_vo,tX2,l1_tensor)
    ! ierr=talsh_tensor_contract("B(c,k)+=A(b,j)*W(j,k,b,c)",B1,Aux_vo,int_t%oovv)
    ierr=talsh_tensor_contract("B(k,c)+=A(b,j)*W(j,k,b,c)",B1,Aux_vo,int_t%oovv)   !XY
    if (debug_local) print*, "ierr l18", ierr


! l=-2; t=1 ; h=0
    ! f3 = - (\ts ck) \ts ai \lambdad mkae \W iemc (ovov)
    ierr=talsh_tensor_init(Aux_ooov)
    if (debug_local) print *, "ierr init f3", ierr
    ierr=talsh_tensor_contract("A(m,k,i,e)+=L(m,k,a,e)*T(a,i)",Aux_ooov,l2_tensor,tX1)
    if (debug_local) print *, "ierr get Aux_ooov", ierr
    ierr=talsh_tensor_contract("B(k,c)+=A(m,k,i,e)*W(e,i,c,m)",B1,Aux_ooov,interm_t%w_vovo,scale=MINUS_ONE)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(m,k,i,e)*W(e,i,c,m)",B1,Aux_ooov,interm_t%w_vovo,scale=MINUS_ONE_HALF) !XY
    if (debug_local) print*, "ierr f3", ierr

    ! f4 = - (\ts ck) \ts ai \lambdad imec \W ekam (vovo)
    ierr=talsh_tensor_init(Aux_vooo)
    ierr=talsh_tensor_contract("A(e,k,i,m)+=T(a,i)*W(e,k,a,m)",Aux_vooo,tX1,interm_t%w_vovo)
    ierr=talsh_tensor_contract("B(k,c)+=A(e,k,i,m)*L(i,m,e,c)",B1,Aux_vooo,l2_tensor,scale=MINUS_ONE)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(e,k,i,m)*L(i,m,e,c)",B1,Aux_vooo,l2_tensor,scale=MINUS_ONE_HALF)     !XY
    if (debug_local) print*, "ierr f4", ierr

    ! i2 = +one_half  (\ts ck) \ts ai \lambdad kief \W efca (vvvv)
    ierr=talsh_tensor_init(Aux_oovv)
    ierr=talsh_tensor_contract("A(k,i,c,a)+=L(k,i,e,f)*W(e,f,c,a)",Aux_oovv,l2_tensor,interm_t%w_vvvv)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,i,c,a)*T(a,i)",B1,Aux_oovv,tX1,scale=ONE_HALF)
    if (debug_local) print*, "ierr i2", ierr

    ! n2= + one_half (\ts ck) \ts ai  \lambda mnca \W kimn (oooo)
    ierr=talsh_tensor_init(Aux_oovv)
    ierr=talsh_tensor_contract("A(k,i,c,a)+=L(m,n,c,a)*W(k,i,m,n)",Aux_oovv,l2_tensor,interm_t%w_oooo)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,i,c,a)*T(a,i)",B1,Aux_oovv,tX1,scale=ONE_HALF)
    if (debug_local) print*, "ierr n2", ierr

! l=-2; t=2 ; h=-1
    ! c9 = -one_half (\ts ck) \td abij \lambdad ijac \F kb (ov)
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(b,c)+=L(i,j,a,c)*T(a,b,i,j)",Aux_vv,l2_tensor,tX2)
    ierr=talsh_tensor_contract("B(k,c)+=A(b,c)*F(k,b)",B1,Aux_vv,interm_t%Fbar_ov,scale=MINUS_ONE_HALF)
    if (debug_local) print*, "ierr c9", ierr

    ! c10 = - one_half (\ts ck) \td abij \lambdad ikab \F jc (ov)
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,j)+=L(i,k,a,b)*T(a,b,i,j)",Aux_oo,l2_tensor,tX2)
    ierr=talsh_tensor_contract("B(k,c)+=A(k,j)*F(j,c)",B1,Aux_oo,interm_t%Fbar_ov,scale=MINUS_ONE_HALF)
    if (debug_local) print*, "ierr c10", ierr

    ! g6 = - (\ts ck) \td abij \lambdad ikeb \W ejac (vovv)
    ierr=talsh_tensor_init(Aux_vovo)
    ierr=talsh_tensor_contract("A(a,k,e,j)+=L(i,k,e,b)*T(a,b,i,j)",Aux_vovo,l2_tensor,tX2)
    ierr=talsh_tensor_contract("B(k,c)+=A(a,k,e,j)*W(e,j,a,c)",B1,Aux_vovo,interm_t%w_vovv,scale=MINUS_ONE) 
    ! ierr=talsh_tensor_contract("B(k,c)+=A(a,k,e,j)*W(e,j,a,c)",B1,Aux_vovo,interm_t%w_vovv,scale=MINUS_ONE_HALF)  !XY
    if (debug_local) print*, "ierr g6", ierr


    ! k6 = + (\ts ck) \td abij \lambdad mjac \W ikmb (ooov)
    ierr=talsh_tensor_init(Aux_ovov)
    if (debug_local) print*, "ierr k6_1", ierr
    ! ierr=talsh_tensor_contract("A(m,b,i,c)+=L(i,k,e,b)*T(a,b,i,j)",Aux_ovov,l2_tensor,tX2)
    ierr=talsh_tensor_contract("A(m,b,i,c)+=L(m,j,a,c)*T(a,b,i,j)",Aux_ovov,l2_tensor,tX2)     !XY
    if (debug_local) print*, "ierr k6_2", ierr
    ierr=talsh_tensor_contract("B(k,c)+=A(m,b,i,c)*W(i,k,m,b)",B1,Aux_ovov,interm_t%w_ooov)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(m,b,i,c)*W(i,k,m,b)",B1,Aux_ovov,interm_t%w_ooov,scale=ONE_HALF)  !XY
    if (debug_local) print*, "ierr k6", ierr

    ! g7 = - one_quarter (\ts ck) \td abij \lambdad ijec \W ekab (vovv)
    ierr=talsh_tensor_init(Aux_vooo)
    ierr=talsh_tensor_contract("A(e,k,i,j)+=T(a,b,i,j)*W(e,k,a,b)",Aux_vooo,tX2,interm_t%w_vovv)
    ierr=talsh_tensor_contract("B(k,c)+=A(e,k,i,j)*L(i,j,e,c)",B1,Aux_vooo,l2_tensor,scale=MINUS_ONE_QUARTER_C)
    if (debug_local) print*, "ierr g7", ierr

    !k7 = + one_quarter (\ts ck) \td abij \lambdad mkab \W ijmc (ooov)
    ierr=talsh_tensor_init(Aux_oooo)
    ierr=talsh_tensor_contract("A(m,k,i,j)+=L(m,k,a,b)*T(a,b,i,j)",Aux_oooo,l2_tensor,tX2)
    ierr=talsh_tensor_contract("B(k,c)+=A(m,k,i,j)*W(i,j,m,c)",B1,Aux_oooo,interm_t%w_ooov,scale=ONE_QUARTER_C)
    if (debug_local) print*, "ierr k7", ierr

    !g9 = +one_half (\ts ck) \td abij  \lambdad ijae \W ekbc (vovv)
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(b,e)+=L(i,j,a,e)*T(a,b,i,j)",Aux_vv,l2_tensor,tX2)
    ierr=talsh_tensor_contract("B(k,c)+=A(b,e)*W(e,k,b,c)",B1,Aux_vv,interm_t%w_vovv,scale=ONE_HALF)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(b,e)*W(e,k,b,c)",B1,Aux_vv,interm_t%w_vovv,scale=ONE_QUARTER_C)   !XY
    if (debug_local) print*, "ierr g9", ierr

    !k9 = - one_half (\ts ck) \td abij  \lambdad imab \W jkmc (ooov)
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(m,j)+=L(i,m,a,b)*T(a,b,i,j)",Aux_oo,l2_Tensor,tX2)
    ierr=talsh_tensor_contract("B(k,c)+=A(m,j)*W(j,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE_HALF)
    ! ierr=talsh_tensor_contract("B(k,c)+=A(m,j)*W(j,k,m,c)",B1,Aux_oo,interm_t%w_ooov,scale=MINUS_ONE_QUARTER_C)  !XY
    if (debug_local) print*, "ierr k9", ierr


! !********************Additional******************************************************************************

! eq l22 bis :: \B kc = -1/2 \ld mnec \V kifa \td efmn \tXs ai                             
    ierr=talsh_tensor_init(Aux_vv)                                                             
    ierr=talsh_tensor_contract("A(f,c)+=L(m,n,e,c)*T(e,f,m,n)",Aux_vv,l2_tensor,interm_t%t2)     
    if(debug_local) write(*,*) "ierr l22 bis  =",ierr                                     
    ierr=talsh_tensor_init(Aux_oovv)                                                           
    ierr=talsh_tensor_contract("B(k,i,c,a)+=A(f,c)*V(k,i,f,a)",Aux_oovv,Aux_vv,int_t%oovv)     
    if(debug_local) write(*,*) "ierr l22 bis =",ierr                                      
    ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                          B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                
    if(debug_local) write(*,*) "ierr l22 bis=",ierr                                       
                                                                                               
! eq l23 bis :: \B kc = -1/2 \ld mkef \V nica \td efmn \tXs ai                             
    ierr=talsh_tensor_init(Aux_oo)                                                             
    ierr=talsh_tensor_contract("A(k,n)+=L(m,k,e,f)*T(e,f,m,n)",Aux_oo,l2_tensor,interm_t%t2)     
    if(debug_local) write(*,*) "ierr l23 bis  =",ierr                                     
    ierr=talsh_tensor_init(Aux_oovv)                                                           
    ierr=talsh_tensor_contract("B(k,i,c,a)+=A(k,n)*V(n,i,c,a)",Aux_oovv,Aux_oo,int_t%oovv)     
    if(debug_local) write(*,*) "ierr l23 bis =",ierr                                      
    ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                          B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                
    if(debug_local) write(*,*) "ierr l23 bis=",ierr                                       
                                                                                               
! eq l32 bis :: \B kc = -1/2 \ld mnea \V ikfc \td efmn \tXs ai                             
    ierr=talsh_tensor_init(Aux_vv)                                                             
    ierr=talsh_tensor_contract("A(f,a)+=L(m,n,e,a)*T(e,f,m,n)",Aux_vv,l2_tensor,interm_t%t2)     
    if(debug_local) write(*,*) "ierr l32 bis  =",ierr                                     
    ierr=talsh_tensor_init(Aux_oovv)                                                           
    ierr=talsh_tensor_contract("B(k,i,c,a)+=A(f,a)*V(i,k,f,c)",Aux_oovv,Aux_vv,int_t%oovv)     
    if(debug_local) write(*,*) "ierr l32 bis =",ierr                                      
    ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                          B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                
    if(debug_local) write(*,*) "ierr l32 bis=",ierr                                       
                                                                                               
! eq l34 bis :: \B kc = -1/2 \ld mief \V nkac \td efmn \tXs ai                             
    ierr=talsh_tensor_init(Aux_oo)                                                             
    ierr=talsh_tensor_contract("A(i,n)+=L(m,i,e,f)*T(e,f,m,n)",Aux_oo,l2_tensor,interm_t%t2)     
    if(debug_local) write(*,*) "ierr l34 bis  =",ierr                                     
    ierr=talsh_tensor_init(Aux_oovv)                                                           
    ierr=talsh_tensor_contract("B(k,i,c,a)+=A(i,n)*V(n,k,a,c)",Aux_oovv,Aux_oo,int_t%oovv)     
    if(debug_local) write(*,*) "ierr l34 bis =",ierr                                      
    ierr=talsh_tensor_contract("tF(k,c)+=B(k,i,c,a)*tx(a,i)", &                                
                          B1,Aux_oovv,tX1,scale=MINUS_ONE_HALF)                                    
    if(debug_local) write(*,*) "ierr l34 bis=",ierr                                       
                                                                                               
! !************************************************************************************************************


! l=-1; t=1 -> h=-2 
    ! l3  = + \lambdas kc (\tYd || \Pmm cd \Pmm kl)  cdkl  \tXs ai  \W lida (oovv)
    ierr=talsh_tensor_init(Aux_ov)
    ierr=talsh_tensor_contract("A(l,d)+=T(a,i)*V(l,i,d,a)",Aux_ov,tX1,int_t%oovv)
    if (debug_local) print*, "ierr l3_1", ierr
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,d)*L(k,c)",B2,Aux_ov,l1_tensor)
    if (debug_local) print*, "ierr l3_2", ierr
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,d)*L(l,c)",B2,Aux_ov,l1_tensor,scale=MINUS_ONE)
    if (debug_local) print*, "ierr l3_3", ierr
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,c)*L(k,d)",B2,Aux_ov,l1_tensor,scale=MINUS_ONE)
    if (debug_local) print*, "ierr l3_4", ierr
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,c)*L(l,d)",B2,Aux_ov,l1_tensor)
    if (debug_local) print*, "ierr l3", ierr


    ! l16 = - \lambdas ka ( \tYd cdkl || \Pmm kl | one_half )  \tXs ai \W ilcd 
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,i)+=L(k,a)*T(a,i)",Aux_oo,l1_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,i)*V(i,l,c,d)",B2,Aux_oo,int_t%oovv,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,i)*V(i,k,c,d)",B2,Aux_oo,int_t%oovv)
    if (debug_local) print*, "ierr l16", ierr

    ! l17 = - \lambdas ic (\tYd cdkl || \Pmm cd | one_half ) \tXs ai \W klad
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(a,c)+=L(i,c)*T(a,i)",Aux_vv,l1_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,c)*V(k,l,a,d)",B2,Aux_vv,int_t%oovv,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,d)*V(k,l,a,c)",B2,Aux_vv,int_t%oovv)
    if (debug_local) print*, "ierr l17", ierr


! l=-2; t=1 -> h=-1

    ! c11 = -  \lambdad kicd ( \tYd cdkl || \Pmm kl | one_half ) \tXs ai \F la (ov)
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(l,i)+=T(a,i)*F(l,a)",Aux_oo,tX1,interm_t%Fbar_ov)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,i)*L(k,i,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,i)*L(l,i,c,d)",B2,Aux_oo,l2_tensor)
    if (debug_local) print*, "ierr c11", ierr

    ! c12 = -one_half (\lambdad klca  \ \Pmm cd | ONE_HALF )\tYd cdkl\tXs ai  \F id (ov) 
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(a,d)+=T(a,i)*F(i,d)",Aux_vv,tX1,interm_t%Fbar_ov)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,d)*L(k,l,c,a)",B2,Aux_vv,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,c)*L(k,l,d,a)",B2,Aux_vv,l2_tensor)
    if (debug_local) print*, "ierr c12", ierr

    ! g10 = + (\lambdad klce (\tYd cdkl || \Pmm cd | ONE_HALF)  \tXs ai \tYd cdkl \W eida (vovv)
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(e,d)+=T(a,i)*W(e,i,d,a)",Aux_vv,tX1,interm_t%w_vovv)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(e,d)*L(k,l,c,e)",B2,Aux_vv,l2_tensor)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(e,c)*L(k,l,d,e)",B2,Aux_vv,l2_tensor,scale=MINUS_ONE)
    if (debug_local) print*, "ierr g10", ierr

    ! k10= - \lambdad kmcd (\tYd cdkl || \Pmm kl | ONE_HALF) \tXs ai \tYd cdkl\W lima 
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(l,m)+=T(a,i)*W(l,i,m,a)",Aux_oo,tX1,interm_t%w_ooov)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,m)*L(k,m,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,m)*L(l,m,c,d)",B2,Aux_oo,l2_tensor)
    if (debug_local) print*, "ierr k10", ierr

    ! k4 = + \lambdad mlca (\tYd cdkl || \Pmm klcd )   \tXs ai  \W kimd (ooov)
    ierr=talsh_tensor_init(Aux_oovo)
    ierr=talsh_tensor_contract("A(m,l,c,i)+=L(m,l,c,a)*T(a,i)",Aux_oovo,l2_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,c,i)*W(k,i,m,d)",B2,Aux_oovo,interm_t%w_ooov)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,c,i)*W(l,i,m,d)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,d,i)*W(k,i,m,c)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,d,i)*W(l,i,m,c)",B2,Aux_oovo,interm_t%w_ooov)

    ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,c,i)*W(k,i,m,d)",B2,Aux_oovo,interm_t%w_ooov,scale=ONE_HALF)
    ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,c,i)*W(l,i,m,d)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE_HALF)
    ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,l,d,i)*W(k,i,m,c)",B2,Aux_oovo,interm_t%w_ooov,scale=MINUS_ONE_HALF)
    ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(m,k,d,i)*W(l,i,m,c)",B2,Aux_oovo,interm_t%w_ooov,scale=ONE_HALF)
    if (debug_local) print*, "ierr k4", ierr

    ! k5 = + \lambdad micd ( \tYd cdkl one_quarter )  \tXs ai \W klma (ooov)
    ierr=talsh_tensor_init(Aux_oooo)
    ierr=talsh_tensor_contract("A(k,l,m,i)+=T(a,i)*W(k,l,m,a)",Aux_oooo,tX1,interm_t%w_ooov)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,l,m,i)*L(m,i,c,d)",B2,Aux_oooo,l2_tensor)
    ! ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,l,m,i)*L(m,i,c,d)",B2,Aux_oooo,l2_tensor,scale=ONE_HALF)  !XY
    if (debug_local) print*, "ierr k5", ierr

    ! g4 = -  \lambdad kied( \tYd cdkl || \Pmm klcd ) \tXs ai \W elca (vovv)
    ierr=talsh_tensor_init(Aux_vovo)
    ! ierr=talsh_tensor_contract("A(e,l,c,i)+=T(a,i)*W(k,l,m,a)",Aux_vovo,tX1,interm_t%w_ooov)
    ierr=talsh_tensor_contract("A(e,l,c,i)+=T(a,i)*W(e,l,c,a)",Aux_vovo,tX1,interm_t%w_vovv)  !XY

    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,l,c,i)*L(k,i,e,d)",B2,Aux_vovo,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,k,c,i)*L(l,i,e,d)",B2,Aux_vovo,l2_tensor)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,l,d,i)*L(k,i,e,c)",B2,Aux_vovo,l2_tensor)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(e,k,d,i)*L(l,i,e,c)",B2,Aux_vovo,l2_tensor,scale=MINUS_ONE)
    if (debug_local) print*, "ierr g4", ierr

    ! g5 = - \lambdad klea (\tYd cdkl one_quarter )  \tXs ai \W eicd (vovv) 
    ierr=talsh_tensor_init(Aux_oovo)
    ierr=talsh_tensor_contract("A(k,l,e,i)+=L(k,l,e,a)*T(a,i)",Aux_oovo,l2_tensor,tX1)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,l,e,i)*W(e,i,c,d)",B2,Aux_oovo,interm_t%w_vovv,scale=MINUS_ONE)
    if (debug_local) print*, "ierr g5", ierr
    

! l=-2; t=2 -> h=-2
    ! l7 = -one_half \lambdad klca (\tYd cdkl || \Pmm cd one_half || one_quarter ) \tXd abij  \W ijdb (oovv)
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(a,d)+=T(a,b,i,j)*W(i,j,d,b)",Aux_vv,tX2,int_t%oovv)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,d)*L(k,l,c,a)",B2,Aux_vv,l2_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(a,c)*L(k,l,d,a)",B2,Aux_vv,l2_tensor,scale=ONE_HALF)
    if (debug_local) print*, "ierr l7", ierr

    ! l9 =  -one_half \lambdad kicd (\tYd cdkl || \Pmm kl one_half|| one_quarter )  \tXd abij  \W ljab (oovv)
    ierr=talsh_tensor_init(Aux_oo) 
    ! ierr=talsh_tensor_contract("A(l,j)+=T(a,b,i,j)*W(l,j,a,b)",Aux_oo,TX2,int_t%oovv)
    ierr=talsh_tensor_contract("A(l,i)+=T(a,b,i,j)*V(l,j,a,b)",Aux_oo,TX2,int_t%oovv) !XY
    if (debug_local) print*, "ierr l9_1", ierr
    ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,j)*L(k,i,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(l,i)*L(k,i,c,d)",B2,Aux_oo,l2_tensor,scale=MINUS_ONE_HALF) !XY
    if (debug_local) print*, "ierr l9_2", ierr
    ! ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,j)*L(l,i,c,d)",B2,Aux_oo,l2_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,i)*L(l,i,c,d)",B2,Aux_oo,l2_tensor,scale=ONE_HALF) !XY
    if (debug_local) print*, "ierr l9", ierr

    ! l12 = +one_sixteenth \lambdad ijcd (\tYd cdkl || one_quarter((ij),(ab))  ) \tXd abij \W klab (oovv)
    ierr=talsh_tensor_init(Aux_oooo)
    ierr=talsh_tensor_contract("A(k,l,i,j)+=T(a,b,i,j)*V(k,l,a,b)",Aux_oooo,tX2,int_t%oovv)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,l,i,j)*L(i,j,c,d)",B2,Aux_oooo,l2_tensor,scale=ONE_QUARTER_C)
    if (debug_local) print*, "ierr l12", ierr

    ! l13 = +one_sixteenth \lambdad klab \tYd cdkl \tXd abij \W ijcd (oovv)
    ierr=talsh_tensor_init(Aux_oooo)
    ierr=talsh_tensor_contract("A(k,l,i,j)+=L(k,l,a,b)*T(a,b,i,j)",Aux_oooo,l2_tensor,tX2)
    ierr=talsh_tensor_contract("B(k,l,c,d)+=A(k,l,i,j)*V(i,j,c,d)",B2,Aux_oooo,int_t%oovv,scale=ONE_QUARTER_C)
    if (debug_local) print*, "ierr l13", ierr

    ! l20 : +\lambdad kjcb (\tYd cdkl | \Pmm kl \Pmm cd) \tXd abij \W lida (oovv)
    ierr=talsh_tensor_init(Aux_ovvo)
    ierr=talsh_tensor_contract("A(l,b,d,j)+=W(l,i,d,a)*T(a,b,i,j)",Aux_ovvo,int_t%oovv,tX2)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,b,d,j)*L(k,j,c,b)",B2,Aux_ovvo,l2_tensor)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,b,d,j)*L(l,j,c,b)",B2,Aux_ovvo,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,b,c,j)*L(k,j,d,b)",B2,Aux_ovvo,l2_tensor,scale=MINUS_ONE)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,b,c,j)*L(l,j,d,b)",B2,Aux_ovvo,l2_tensor)
    if (debug_local) print*, "ierr l20", ierr

    ! l30 = - one_half \lambdad mnef (\tYd cdkl | \Pmm kl  one_half (cd)  ) \tXd abij \W jlcd
    ierr=talsh_tensor_init(Aux_oo)
    ierr=talsh_tensor_contract("A(k,j)+=L(i,k,a,b)*T(a,b,i,j)",Aux_oo,l2_tensor,tX2)
    ! ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,j)*W(j,l,c,d)",B2,Aux_oo,int_t%oovv,scale=ONE_HALF)
    ! ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,j)*W(j,k,c,d)",B2,Aux_oo,int_t%oovv,scale=MINUS_ONE_HALF)

    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(k,j)*W(j,l,c,d)",B2,Aux_oo,int_t%oovv,scale=MINUS_ONE_HALF)  !XY
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(l,j)*W(j,k,c,d)",B2,Aux_oo,int_t%oovv,scale=ONE_HALF) !XY
    if (debug_local) print*, "ierr l30", ierr

    ! l31 = - one_half \lambdad ijac (\tYd cdkl | \Pmm cd one_half (ij)  ) \tXd abij \W klbd
    ierr=talsh_tensor_init(Aux_vv)
    ierr=talsh_tensor_contract("A(b,c)+=L(i,j,a,c)*T(a,b,i,j)",Aux_vv,l2_tensor,tX2)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(b,c)*W(k,l,b,d)",B2,Aux_vv,int_t%oovv,scale=MINUS_ONE_HALF)
    ierr=talsh_Tensor_contract("B(k,l,c,d)+=A(b,d)*W(k,l,b,c)",B2,Aux_vv,int_t%oovv,scale=ONE_HALF)
    if (debug_local) print*, "ierr l31", ierr



    ierr=talsh_tensor_destruct(Aux_oo)
    ierr=talsh_tensor_destruct(Aux_vv)
    ierr=talsh_tensor_destruct(Aux_vo)
    ierr=talsh_tensor_destruct(Aux_ov)
    
    ierr=talsh_tensor_destruct(Aux_vvvv)
    ierr=talsh_tensor_destruct(Aux_oooo)
    ierr=talsh_tensor_destruct(Aux_oovv)
    ierr=talsh_tensor_destruct(Aux_vooo)
    ierr=talsh_tensor_destruct(Aux_ooov)
    ierr=talsh_tensor_destruct(Aux_vovo)
    ierr=talsh_tensor_destruct(Aux_ovov)
    ierr=talsh_tensor_destruct(Aux_vovv)
    ierr=talsh_tensor_destruct(Aux_oovv)

    ierr=talsh_tensor_destruct(Aux_oovo)
    ierr=talsh_tensor_destruct(Aux_ovvo)


    end subroutine get_tF


end module
