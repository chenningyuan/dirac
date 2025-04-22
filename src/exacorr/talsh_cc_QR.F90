module talsh_cc_QR
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global
    use talsh_eom_LR, only: get_tF

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), IMAGINARY=(0.0D0,1.0D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.D0,0.D0), &
                                ONE_SIXTH=(0.16666666666666666666D0,0.D0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public get_cccc_qr

    contains


    
!-----------------------------------------------------------------------------------------------------------------

subroutine get_cccc_qr(exa_input, int_t, interm_t,&
                        Interm_t_tensors_A,Interm_t_tensors_B, Interm_t_tensors_C,& 
                        Tprop_tensors, Tprop_bar_tensors,& 
                        Tprop_tensors_min_w, Tprop_bar_tensors_min_w,&
                        i, j, k, l1_tensor, l2_tensor, cccc_qr)
   
!   This subroutine drives the CC-CC Quadratic response property calculation
!   Written by Xiang Yuan Janunary 2022, Lille, France

        use talsh_sigma_eomccsd
        type(exacc_input), intent(in) :: exa_input
        type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_tens_t), intent(inout)           :: l2_tensor
        type(talsh_tens_t), intent(inout)           :: l1_tensor
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_A(:)
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_B(:)
        type(talsh_rsp_interm_tens_t), allocatable :: Interm_t_tensors_C(:)


        type(talsh_comm_tens_t), allocatable :: Tprop_tensors(:), Tprop_tensors_min_w(:) !LR
        type(talsh_comm_tens_t), allocatable :: Tprop_bar_tensors(:), Tprop_bar_tensors_min_w(:) !QR

        integer                       :: i, j, k
        
        complex(kind=8)               :: quarsp, quarsp_con
        complex(kind=8)               :: phase_factor 

        integer                       :: ierr
        integer                       :: T

        complex(kind=8)               :: cccc_qr
        real(kind=8)                  :: tem_real, tem_imag


        ! type(talsh_comm_tens_t)     :: bvec_tensor
        type(talsh_tens_t)          :: result_tensor
        integer(C_INT)              :: one_dims(1)
        complex(8), pointer         :: result_tens(:)
        type(C_PTR):: body_p


        call print_date('Evaluate CC-CC QR')

        call get_cccc_quarspf(exa_input,    &
                        int_t,              &
                        interm_t,           &
                        Interm_t_tensors_A, &
                        Interm_t_tensors_B, &
                        Interm_t_tensors_C, &
                        Tprop_tensors,      &
                        Tprop_bar_tensors,  & 
                        i, j, k,            &
                        l1_tensor,          &
                        l2_tensor,          &
                        quarsp)

        call get_cccc_quarspf(exa_input,    &
                        int_t,              &
                        interm_t,           &
                        Interm_t_tensors_A, &
                        Interm_t_tensors_B, &
                        Interm_t_tensors_C, &
                        Tprop_tensors_min_w,&
                        Tprop_bar_tensors_min_w,& 
                        i, j, k,            &
                        l1_tensor,          &
                        l2_tensor,          &
                        quarsp_con)


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

        quarsp     = quarsp*phase_factor
        quarsp_con = quarsp_con*phase_factor

        print *, "quadratic response: +w", quarsp
        print *, "quadratic response: -w", quarsp_con

        quarsp_con=conjg(quarsp_con)
        cccc_qr=(quarsp+quarsp_con)*0.5

        print *, "EOM-QR after symmetrization ",cccc_qr
end subroutine get_cccc_qr

!------------------------------------------------------------------------------------------------------------------


subroutine get_cccc_quarspf(exa_input, int_t, interm_t,&
            Interm_t_tensors_A, Interm_t_tensors_B, Interm_t_tensors_C,&
            Tprop_tensors,Tprop_bar_tensors,&
            i, j, k, l1_tensor, l2_tensor, quarsp)
                 
        type(exacc_input), intent(in) :: exa_input
        type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_rsp_interm_tens_t), allocatable     :: Interm_t_tensors_A(:)
        type(talsh_rsp_interm_tens_t), allocatable     :: Interm_t_tensors_B(:)
        type(talsh_rsp_interm_tens_t), allocatable     :: Interm_t_tensors_C(:)
        type(talsh_comm_tens_t), allocatable           :: Tprop_tensors(:)
        type(talsh_comm_tens_t), allocatable           :: Tprop_bar_tensors(:)
        type(talsh_tens_t), intent(inout)              :: l2_tensor
        type(talsh_tens_t), intent(inout)              :: l1_tensor
        integer, intent(inout)                         :: i, j, k
        complex(kind=8), intent(inout)                 :: quarsp


        type(talsh_comm_tens_t)       :: MAB_tensor
        complex(kind=8)               :: tem_quarsp_iter
        type(talsh_comm_tens_t)       :: bvec_tensor
        type(talsh_tens_t)            :: result_tensor
        integer(INTD) :: ierr

        integer, dimension(2) :: ov_dims
        integer, dimension(4) :: oovv_dims

        print *, "Constuct CC-CC QR response function"

        ! we get the number of occupied and virtuals to build MAB, from
        ! interm_t%nocc, interm_t%nvir and build the scratch MAB_tensors,

        ! in the following we use the notation M(x,y)t(z), and obtain the contributions from all
        ! permutations,that is:
        ! [<<x,y,z>>]_1 = M(x,y)t(z) + M(z,x)t(y) + M(y,z)t(x) 
        !               + M(y,x)t(z) + M(x,z)t(y) + M(z,y)t(x)
        ! or
        ! [<<x,y,z>>]_2 = M(-x,-y)t(-z) + M(-z,-x)t(-y) + M(-y,-z)t(-x) 
        !               + M(-y,-x)t(-z) + M(-x,-z)t(-y) + M(-z,-y)t(-x)
        !
        ! in each case, M(x,y) will have a component M^{XY}_ai, M(X,Y)_abij, etc since they
        ! will contract with the perturbed amplitudes
        !
        ! depending on the perturbed amplitudes that are passed on to this routine
        ! and at the top-level driver we assemble the final quadratic response function
        ! <<x,y,z>> = 0.5*{[<<x,y,z>>]_1 + ([<<x,y,z>>]_2)^*} 

        ! X is index 1; Y is index 2; and Z is index 3 for amplitudes e.g.
        !    Tprop_tensors(1) == T^x(\omega_z)
        !    Tprop_tensors(2) == T^y(\omega_y)
        !    Tprop_tensors(3) == T^z(\omega_z)

        ov_dims=(/interm_t%nocc,interm_t%nvir/)
        oovv_dims=(/interm_t%nocc,interm_t%nocc,interm_t%nvir,interm_t%nvir/)

        ierr=talsh_tensor_construct(result_tensor,C8,(/1/),init_val=ZERO)
        ierr=talsh_tensor_construct(MAB_tensor%t1,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(MAB_tensor%t2,C8,oovv_dims,init_val=ZERO)

        ! now go for the permutations

        ! {XY,YX}Z combinations

        ! initialize MAB to accumulate M^{XY} 
        ierr = talsh_tensor_init(MAB_tensor%t1, ZERO)
        ierr = talsh_tensor_init(MAB_tensor%t2, ZERO)
        ! calculate M^{XY} 
        call form_cccc_qr_MXY(MAB_tensor,       & ! M^{XY}
                         int_t, interm_t,       &
                         Interm_t_tensors_A(i), & ! prop X(i)
                         Interm_t_tensors_B(j), & ! prop Y(j)
                         Tprop_tensors(1),      & ! tX
                         Tprop_bar_tensors(1),  & ! tbarX
                         Tprop_tensors(2),      & ! tY
                         l1_tensor, l2_tensor)    ! ground-state lambda tensors
        !
        ! contract M^{XY}*t^Z (=Tprop_tensors(3)) and accumulate into scalar
        ierr=talsh_tensor_contract("MXYtZ+=MXY(m,e)*tZ(e,m)",result_tensor,MAB_tensor%t1,Tprop_tensors(3)%t1)
        ierr=talsh_tensor_contract("MXYtZ+=MXY(m,n,e,f)*tZ(e,f,m,n)",result_tensor,MAB_tensor%t2,Tprop_tensors(3)%t2,scale=ONE_QUARTER_C)

        ! initialize MAB to accumulate M^{YX} 
        ierr = talsh_tensor_init(MAB_tensor%t1, ZERO)
        ierr = talsh_tensor_init(MAB_tensor%t2, ZERO)
        ! calculate M^{YX} 
        call form_cccc_qr_MXY(MAB_tensor,       & ! M^{YX}
                         int_t, interm_t,       &
                         Interm_t_tensors_B(j), & ! prop Y(j)
                         Interm_t_tensors_A(i), & ! prop X(i)
                         Tprop_tensors(2),      & ! ty
                         Tprop_bar_tensors(2),  & ! tbarY
                         Tprop_tensors(1),      & ! tX
                         l1_tensor, l2_tensor)    ! ground-state lambda tensors
        !
        ! contract M^{YX} with t^Z (=Tprop_tensors(3)) and accumulate into scalar
        ierr=talsh_tensor_contract("MYXtZ+=MYX(m,e)*tZ(e,m)",result_tensor,MAB_tensor%t1,Tprop_tensors(3)%t1)
        ierr=talsh_tensor_contract("MYXtZ+=MYX(m,n,e,f)*tZ(e,f,m,n)",result_tensor,MAB_tensor%t2,Tprop_tensors(3)%t2,scale=ONE_QUARTER_C)


        ! {YZ,ZY}X combinations

        ! initialize MAB to accumulate M^{YZ} 
        ierr = talsh_tensor_init(MAB_tensor%t1, ZERO)
        ierr = talsh_tensor_init(MAB_tensor%t2, ZERO)
        ! calculate M^{XY} 
        call form_cccc_qr_MXY(MAB_tensor,       & ! M^{YZ}
                         int_t, interm_t,       &
                         Interm_t_tensors_B(j), & ! prop Y(j)
                         Interm_t_tensors_C(k), & ! prop Z(k)
                         Tprop_tensors(2),      & ! tY
                         Tprop_bar_tensors(2),  & ! tbarY
                         Tprop_tensors(3),      & ! tZ
                         l1_tensor, l2_tensor)    ! ground-state lambda tensors
        !
        ! contract M^{YZ}*t^X (=Tprop_tensors(1)) and accumulate into scalar
        ierr=talsh_tensor_contract("MYZtX+=MYZ(m,e)*tX(e,m)",result_tensor,MAB_tensor%t1,Tprop_tensors(1)%t1)
        ierr=talsh_tensor_contract("MYZtX+=MYZ(m,n,e,f)*tX(e,f,m,n)",result_tensor,MAB_tensor%t2,Tprop_tensors(1)%t2,scale=ONE_QUARTER_C)

        ! initialize MAB to accumulate M^{ZY} 
        ierr = talsh_tensor_init(MAB_tensor%t1, ZERO)
        ierr = talsh_tensor_init(MAB_tensor%t2, ZERO)
        ! calculate M^{ZY} 
        call form_cccc_qr_MXY(MAB_tensor,       & ! M^{ZY}
                         int_t, interm_t,       &
                         Interm_t_tensors_C(k), & ! prop Z(k)
                         Interm_t_tensors_B(j), & ! prop Y(j)
                         Tprop_tensors(3),      & ! tZ
                         Tprop_bar_tensors(3),  & ! tbarZ
                         Tprop_tensors(2),      & ! tY
                         l1_tensor, l2_tensor)    ! ground-state lambda tensors
        !
        ! contract M^{YX} with t^Z (=Tprop_tensors(3)) and accumulate into scalar
        ierr=talsh_tensor_contract("MYZtX+=MZY(m,e)*tX(e,m)",result_tensor,MAB_tensor%t1,Tprop_tensors(1)%t1)
        ierr=talsh_tensor_contract("MYZtX+=MZY(m,n,e,f)*tX(e,f,m,n)",result_tensor,MAB_tensor%t2,Tprop_tensors(1)%t2,scale=ONE_QUARTER_C)


        ! {ZX,XZ}Y combinations

        ! initialize MAB to accumulate M^{ZX} 
        ierr = talsh_tensor_init(MAB_tensor%t1, ZERO)
        ierr = talsh_tensor_init(MAB_tensor%t2, ZERO)
        ! calculate M^{XY} 
        call form_cccc_qr_MXY(MAB_tensor,       & ! M^{ZX}
                         int_t, interm_t,       &
                         Interm_t_tensors_C(k), & ! prop Z(k)
                         Interm_t_tensors_A(i), & ! prop X(i)
                         Tprop_tensors(3),      & ! tZ
                         Tprop_bar_tensors(3),  & ! tbarZ
                         Tprop_tensors(1),      & ! tX
                         l1_tensor, l2_tensor)    ! ground-state lambda tensors
        !
        ! contract M^{YZ}*t^X (=Tprop_tensors(1)) and accumulate into scalar
        ierr=talsh_tensor_contract("MZXtY+=MZX(m,e)*tY(e,m)",result_tensor,MAB_tensor%t1,Tprop_tensors(2)%t1)
        ierr=talsh_tensor_contract("MZXtY+=MZX(m,n,e,f)*tY(e,f,m,n)",result_tensor,MAB_tensor%t2,Tprop_tensors(2)%t2,scale=ONE_QUARTER_C)

        ! initialize MAB to accumulate M^{XZ} 
        ierr = talsh_tensor_init(MAB_tensor%t1, ZERO)
        ierr = talsh_tensor_init(MAB_tensor%t2, ZERO)
        ! calculate M^{ZY} 
        call form_cccc_qr_MXY(MAB_tensor,       & ! M^{XZ}
                         int_t, interm_t,       &
                         Interm_t_tensors_A(i), & ! prop X(i)
                         Interm_t_tensors_C(k), & ! prop Z(k)
                         Tprop_tensors(1),      & ! tX
                         Tprop_bar_tensors(1),  & ! tbarX
                         Tprop_tensors(3),      & ! tZ
                         l1_tensor, l2_tensor)    ! ground-state lambda tensors
        !
        ! contract M^{YX} with t^Z (=Tprop_tensors(3)) and accumulate into scalar
        ierr=talsh_tensor_contract("MXZtY+=MXZ(m,e)*tY(e,m)",result_tensor,MAB_tensor%t1,Tprop_tensors(2)%t1)
        ierr=talsh_tensor_contract("MXZtY+=MXZ(m,n,e,f)*tY(e,f,m,n)",result_tensor,MAB_tensor%t2,Tprop_tensors(2)%t2,scale=ONE_QUARTER_C)


        ! destroy the MAB tensors

        ierr=talsh_tensor_destruct(MAB_tensor%t1)
        ierr=talsh_tensor_destruct(MAB_tensor%t2)
        ierr=talsh_tensor_destruct(result_tensor)


end subroutine get_cccc_quarspf

subroutine form_cccc_qr_MXY(MAB_tensor, int_t, interm_t, X_operator, Y_operator, &
                            tx, tbar_x, ty, l1_tensor, l2_tensor) 
        type(talsh_comm_tens_t), intent(inout)         :: MAB_tensor
        type(talsh_intg_tens_t),intent(inout)          :: int_t     !two-electron integral tensor
        type(talsh_intermediates_tens_t),intent(inout) :: interm_t  !t1, t2 tensor
        type(talsh_rsp_interm_tens_t), intent(inout)   :: X_operator
        type(talsh_rsp_interm_tens_t), intent(inout)   :: Y_operator
        type(talsh_comm_tens_t), intent(inout)         :: tx, tbar_x, ty
        type(talsh_tens_t), intent(inout)              :: l1_tensor
        type(talsh_tens_t), intent(inout)              :: l2_tensor

        type(talsh_comm_tens_t) :: tYFX_tensor
        type(talsh_comm_tens_t) :: GtXtY_tensor
        type(talsh_comm_tens_t) :: tbarXBtY_tensor
        type(talsh_comm_tens_t) :: tbarXAY_tensor

        integer, dimension(2) :: ov_dims
        integer, dimension(4) :: oovv_dims
        integer(INTD) :: ierr

! local variables

        ov_dims=(/interm_t%nocc,interm_t%nvir/)
        oovv_dims=(/interm_t%nocc,interm_t%nocc,interm_t%nvir,interm_t%nvir/)

        ierr=talsh_tensor_construct(tYFX_tensor%t1,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(tYFX_tensor%t2,C8,oovv_dims,init_val=ZERO)

        ierr=talsh_tensor_construct(GtXtY_tensor%t1,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(GtXtY_tensor%t2,C8,oovv_dims,init_val=ZERO)

        ierr=talsh_tensor_construct(tbarXBtY_tensor%t1,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(tbarXBtY_tensor%t2,C8,oovv_dims,init_val=ZERO)
      
        ierr=talsh_tensor_construct(tbarXAY_tensor%t1,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(tbarXAY_tensor%t2,C8,oovv_dims,init_val=ZERO)

        ! add (1/2) * (tyFx)^m_e  to (MAB_tensor)^m_e
        ! add (1/2) * (tyFx)^{mn}_{ef}  to (MAB_tensor)^{mn}_{ef}

! add equations 7 and 8 from the paper, 

        call get_tYFX(tYFX_tensor,X_operator%prop_bar_ov,ty%t1,ty%t2,l1_tensor,l2_tensor,interm_t%nocc,interm_t%nvir)

! accumulate intermediate into results with correct scaling due to permutations

        ierr=talsh_tensor_add("Mab(m,e)+=tYFX(m,e)",MAB_tensor%t1,tYFX_tensor%t1,scale=ONE_HALF)
        ierr=talsh_tensor_add("Mab(m,n,e,f)+=tYFX(m,n,e,f)",MAB_tensor%t2,tYFX_tensor%t2,scale=ONE_HALF)

        ! add (1/6) * (Gtxty)^m_e to (MAB_tensor)^m_e
        ! add (1/6) * (Gtxty)^{mn}_{ef} to (MAB_tensor)^{mn}_{ef}

        call get_Gtt(GtXtY_tensor,int_t,interm_t%t1,tx%t1,tx%t2,ty%t1,ty%t2,l1_tensor,l2_tensor,interm_t%nocc, interm_t%nvir)

! accumulate intermediate into results with correct scaling due to permutations

        ierr=talsh_tensor_add("Mab(m,e)+=Gtytx(m,e)",MAB_tensor%t1,GtXtY_tensor%t1,scale=ONE_SIXTH)
        ierr=talsh_tensor_add("Mab(m,n,e,f)+=Gtytx(m,n,e,f)",MAB_tensor%t2,GtXtY_tensor%t2,scale=ONE_SIXTH)

        ! add (1/2) *(lx Bty)^m_e to (MAB_tensor)^m_e
        ! add (1/2) *(lx Bty)^m_e to (MAB_tensor)^{mn}_{ef}

        ! recall we reuse Ft to obtain BtY, with special case (remove U1 vua is_qrcc=.true.) AND at the same time
        ! replace l1, l2 by perburtbed l1, l2. which were calculated before (tbar_x)
        call get_tF(tbarXBtY_tensor%t1,tbarXBtY_tensor%t2, tx%t1,tx%t2, tbar_x%t1,tbar_x%t2, &
                    interm_t,int_t,interm_t%nocc, interm_t%nvir, is_qrcc=.true.)

! accumulate intermediate into results with correct scaling due to permutations

        ierr=talsh_tensor_add("Mab(m,e)+=tbarxBty(m,e)",MAB_tensor%t1,tbarXBtY_tensor%t1,scale=ONE_HALF)
        ierr=talsh_tensor_add("Mab(m,n,e,f)+=tbarxBty(m,n,e,f)",MAB_tensor%t2,tbarXBtY_tensor%t2,scale=ONE_HALF)

        ! add        (lt Ay)^m_e  to (MAB_tensor)^m_e
        ! add        (lt Ay)^m_e  to (MAB_tensor)^{mn}_{ef}

        call get_lXAY(tbarXAY_tensor,Y_operator%prop_bar_oo,Y_operator%prop_bar_vv,Y_operator%prop_bar_ov,&
                & interm_t%t1,interm_t%t2, tbar_x%t1,tbar_x%t2,interm_t%nocc, interm_t%nvir)

! accumulate intermediate into results with correct scaling due to permutations

        ierr=talsh_tensor_add("Mab(m,e)+=tbarxAY(m,e)",MAB_tensor%t1,tbarXAY_tensor%t1)
        ierr=talsh_tensor_add("Mab(m,n,e,f)+=tbarxAY(m,n,e,f)",MAB_tensor%t2,tbarXAY_tensor%t2)

! cleanup of local variables

        ierr=talsh_tensor_destruct(tYFX_tensor%t1)
        ierr=talsh_tensor_destruct(tYFX_tensor%t2)

        ierr=talsh_tensor_destruct(GtXtY_tensor%t1)
        ierr=talsh_tensor_destruct(GtXtY_tensor%t2)

        ierr=talsh_tensor_destruct(tbarXBtY_tensor%t1)
        ierr=talsh_tensor_destruct(tbarXBtY_tensor%t2)
      
        ierr=talsh_tensor_destruct(tbarXAY_tensor%t1)
        ierr=talsh_tensor_destruct(tbarXAY_tensor%t2)

end subroutine form_cccc_qr_MXY

!--------------------------------------------------------------------------------------------------------------
!
! below are routines for the contractions needed to assemble the standand CC QR functions
! code by Loic Halbert, 06/05/2023, to be verified
!
!
  subroutine get_lXAY(tbarXAY,Yoo,Yvv,Yov,t1_tensor,t2_tensor,l1_tensor,l2_tensor,nocc,nvir)
    type(talsh_comm_tens_t), intent(inout) :: tbarXAY ! result
    type(talsh_tens_t), intent(inout)      :: Yoo, Yvv, Yov ! blocks of the property integral tensor
    type(talsh_tens_t), intent(inout)      :: t1_tensor, t2_tensor ! ground state amplitudes
    type(talsh_tens_t), intent(inout)      :: l1_tensor, l2_tensor ! pertubed X lambdas
    integer,             intent(in   )     :: nocc, nvir

    integer, dimension(2) :: oo_dims, vv_dims
    type(talsh_tens_t)  :: oo_aux_tensor, vv_aux_tensor
    integer(INTD) :: ierr
 
    oo_dims=(/nocc,nocc/)
    vv_dims=(/nvir,nvir/)

    ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)

    ierr=talsh_tensor_add("W(m,o)+=Y(m,o)",oo_aux_tensor,Yoo)
    if(ierr.ne.0) write(*,*) "ierr b1-a  =",ierr
    ierr=talsh_tensor_contract("W(m,o)+=Y(m,g)*t(g,o)",oo_aux_tensor,Yov,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr d3-a  =",ierr 
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(m,o)*L(o,e)",tbarXAY%t1,oo_aux_tensor,l1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr b1/d3  =",ierr 

    ierr=talsh_tensor_add("W(g,e)+=Y(g,e)",vv_aux_tensor,Yvv)
    if(ierr.ne.0) write(*,*) "ierr c1-a  =",ierr
    ierr=talsh_tensor_contract("W(g,e)+=t(g,o)*Y(o,e)",vv_aux_tensor,t1_tensor,Yov)
    if(ierr.ne.0) write(*,*) "ierr d3-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(g,e)*L(m,g)",tbarXAY%t1,vv_aux_tensor,l1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr c1/d3  =",ierr

    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_contract("W(m,p)+=L(o,m,g,h)*T(g,h,o,p)",oo_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr d4-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(m,p)*Y(p,e)",tbarXAY%t1,oo_aux_tensor,Yov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr d4-b  =",ierr

    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("W(h,e)+=L(o,p,g,e)*T(g,h,o,p)",vv_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr d5-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,e)+=W(h,e)*Y(m,h)",tbarXAY%t1,vv_aux_tensor,Yov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr d5-b  =",ierr

    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,e)*Y(n,f)",tbarXAY%t2,l1_tensor,Yov)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,f)*Y(n,e)",tbarXAY%t2,l1_tensor,Yov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(n,e)*Y(m,f)",tbarXAY%t2,l1_tensor,Yov,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(n,f)*Y(m,e)",tbarXAY%t2,l1_tensor,Yov)
    if(ierr.ne.0) write(*,*) "ierr d2-d  =",ierr

    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_add("W(m,o)+=Y(m,o)",oo_aux_tensor,Yoo)
    if(ierr.ne.0) write(*,*) "ierr b2-a =",ierr
    ierr=talsh_tensor_contract("W(m,o)+=Y(m,g)*t(g,o)",oo_aux_tensor,Yov,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr d4-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(o,n,e,f)*W(m,o)",tbarXAY%t2,l2_tensor,oo_aux_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(o,m,e,f)*W(n,o)",tbarXAY%t2,l2_tensor,oo_aux_tensor)
    if(ierr.ne.0) write(*,*) "ierr b2/d4  =",ierr

    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_add("W(g,e)+=Y(g,e)",vv_aux_tensor,Yvv)
    if(ierr.ne.0) write(*,*) "ierr c2-a =",ierr
    ierr=talsh_tensor_contract("W(g,e)+=Y(o,e)*t(g,o)",vv_aux_tensor,Yov,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr d5-a  =",ierr
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,n,g,f)*W(g,e)",tbarXAY%t2,l2_tensor,vv_aux_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("lXAY(m,n,e,f)+=L(m,n,g,e)*W(g,f)",tbarXAY%t2,l2_tensor,vv_aux_tensor)
    if(ierr.ne.0) write(*,*) "ierr c2/d5  =",ierr

    ierr=talsh_tensor_destruct(oo_aux_tensor)
    ierr=talsh_tensor_destruct(vv_aux_tensor)

    end subroutine get_lXAY

    subroutine get_tYFX(tYFX,Xov,tYs,tYd,l1_tensor,l2_tensor,nocc,nvir)
       ! result, intermediate for QR
       type(talsh_comm_tens_t), intent(inout) :: tYFX
       ! ov block of property integrals
       type(talsh_tens_t), intent(inout)      :: Xov 
       ! y-perturbed t amplitudes
       type(talsh_tens_t), intent(inout)      :: tYs,tYd 
       ! ground-state lambdas
       type(talsh_tens_t),  intent(inout)   :: l1_tensor, l2_tensor 
       integer,             intent(in   )     :: nocc, nvir

       type(talsh_tens_t) :: oovv_aux_tensor 
       type(talsh_tens_t) :: ov_aux_tensor, vv_aux_tensor, oo_aux_tensor

       integer, dimension(2) :: oo_dims, vv_dims, ov_dims
       integer, dimension(4) :: oovv_dims
       integer(INTD) :: ierr


    oo_dims=(/nocc,nocc/)
    vv_dims=(/nvir,nvir/)
    ov_dims=(/nocc,nvir/)
    oovv_dims=(/nocc,nocc,nvir,nvir/)

    ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(oovv_aux_tensor,C8,oovv_dims,init_val=ZERO)  

    ierr=talsh_tensor_contract("W(k,m,e,c)+=L(k,e)*X(m,c)",oovv_aux_tensor,l1_tensor,Xov)
    if(ierr.ne.0) write(*,*) "ierr e7-a  =",ierr
    ierr=talsh_tensor_contract("W(k,m,e,c)+=L(m,c)*X(k,e)",oovv_aux_tensor,l1_tensor,Xov)
    if(ierr.ne.0) write(*,*) "ierr e8-a  =",ierr           
    ierr=talsh_tensor_contract("tYFX(m,e)+=W(k,m,e,c)*tY(c,k)",tYFX%t1,oovv_aux_tensor,tYs, scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr e7-e8  =",ierr           

    ierr=talsh_tensor_contract("W(d,e)+=L(k,l,c,e)*tY(c,d,k,l)",vv_aux_tensor,l2_tensor,tYd)
    if(ierr.ne.0) write(*,*) "ierr e9-a  =",ierr
    ierr=talsh_tensor_contract("tYFX(m,e)+=W(d,e)*X(m,d)",tYFX%t1,vv_aux_tensor,Xov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr e9-b  =",ierr           
    ierr=talsh_tensor_contract("W(m,l)+=L(k,m,c,d)*tY(c,d,k,l)",oo_aux_tensor,l2_tensor,tYd)
    if(ierr.ne.0) write(*,*) "ierr e10-a  =",ierr
    ierr=talsh_tensor_contract("tYFX(m,e)+=W(m,l)*X(l,e)",tYFX%t1,oo_aux_tensor,Xov,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr e10-b  =",ierr           

    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_contract("W(m,k)+=X(m,c)*tY(c,k)",oo_aux_tensor,Xov,tYs)
    if(ierr.ne.0) write(*,*) "ierr e11-a  =",ierr           
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(m,k)*L(k,n,e,f)",tYFX%t2,oo_aux_tensor,L2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(n,k)*L(k,m,e,f)",tYFX%t2,oo_aux_tensor,L2_tensor)
    if(ierr.ne.0) write(*,*) "ierr e11-b  =",ierr

    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("W(c,e)+=X(k,e)*tY(c,k)",vv_aux_tensor,Xov,tYs)
    if(ierr.ne.0) write(*,*) "ierr e12-a  =",ierr           
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(c,e)*L(m,n,c,f)",tYFX%t2,vv_aux_tensor,L2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("tYFX(m,n,e,f)+=W(c,f)*L(m,n,c,e)",tYFX%t2,vv_aux_tensor,L2_tensor)
    if(ierr.ne.0) write(*,*) "ierr e12-b  =",ierr

    ierr=talsh_tensor_destruct(oovv_aux_tensor)
    ierr=talsh_tensor_destruct(vv_aux_tensor)
    ierr=talsh_tensor_destruct(oo_aux_tensor)
    ierr=talsh_tensor_destruct(ov_aux_tensor)

    end subroutine get_tYFX

    subroutine get_Gtt(Gtt_tensor,int_t,t1_tensor,tXs,tXd,tYs,tYd,l1_tensor,l2_tensor,nocc,nvir)
     type(talsh_comm_tens_t), intent(inout) :: Gtt_tensor 
     type(talsh_intg_tens_t), intent(inout) :: int_t
     type(talsh_tens_t),  intent(inout) :: t1_tensor !ground-state  t1 amplitudes
     type(talsh_tens_t),  intent(inout) :: l1_tensor, l2_tensor !ground-state lambdas
     type(talsh_tens_t), intent(inout)    :: tXs,tXd,tYs,tYd ! perturbed Tx, Ty
     integer,             intent(in   )   :: nocc, nvir


     integer(INTD) :: ierr
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
     call Gtt_k2k3g2g3(Gtt_tensor%t1,int_t,G1_tensor,tYs,tYd,l1_tensor,l2_tensor,nocc,nvir)

     ! k3b / g3b   : inverse
     ierr=talsh_tensor_init(G1_tensor) 
     ! 1. Create the loop  
     call Gtt_loop(G1_tensor,int_t,tYs)
     ! 2. use fov
     call Gtt_k2k3g2g3(Gtt_tensor%t1,int_t,G1_tensor,tXs,tXd,l1_tensor,l2_tensor,nocc,nvir)
       
     ! c7b,c8b : same by inversion \tX et tY
     call Gtt_c7c8(Gtt_tensor%t1,tXs,tYs,l1_tensor,int_t,nocc,nvir)
     call Gtt_c7c8(Gtt_tensor%t1,tXs,tYs,l1_tensor,int_t,nocc,nvir)

     ! f3b, f3d =inv(f4c,f3e) / f3c, f4e =inv( f4b, f4d) 
     ! f3b + f3d + f4c + f3e
            ! create V + V*t
     ierr=talsh_tensor_init(ooov_aux_tensor) 
     ierr=talsh_tensor_add("W(o1,o2,o3,v1)+=V(o1,o2,o3,v1)",ooov_aux_tensor,int_t%ooov)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d W_a  =",ierr         
     ierr=talsh_tensor_contract("W(m,i,o,c)+=T(g,o)*V(m,i,g,c)",ooov_aux_tensor,t1_tensor,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr f3b + f3d W_b  =",ierr

     ! f3b + f3d
     call f3bf3d(ooov_aux_tensor,Gtt_tensor%t1,l2_tensor,tXs,tYs,nocc,nvir)
     ! f4c + f3e 
     call f3bf3d(ooov_aux_tensor,Gtt_tensor%t1,l2_tensor,tYs,tXs,nocc,nvir)

     ! now use intermediate V+V*t for n2b+n2c
     ierr=talsh_tensor_init(ovov_aux_tensor)
     ierr=talsh_tensor_contract("W(i,c,n,e)+=W(i,k,n,e)*tY(c,k)",ovov_aux_tensor,ooov_aux_tensor,tYs)
     if(ierr.ne.0) write(*,*) "ierr n2b+n2d _c  =",ierr
            ! contraction with tXs ai
     ierr=talsh_tensor_init(vvvo_aux_tensor)
     ierr=talsh_tensor_contract("W(a,c,n,e)+=W(i,c,n,e)*tX(a,i)",vvvo_aux_tensor,oooo_aux_tensor,tXs)
     if(ierr.ne.0) write(*,*) "ierr n2b _d  =",ierr   
     ierr=talsh_tensor_contract("G(m,e)+=W(a,c,n,e)*L(n,m,a,c)",Gtt_tensor%t1,vvvo_aux_tensor,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr n2d+n2d _e  =",ierr   
     
     ! f4b + f4d + f3c + f4e
         ! create V + V*t
     ierr=talsh_tensor_init(vovv_aux_tensor)
     ierr=talsh_tensor_add("W(v1,o1,v2,v3)+=V(v1,o1,v2,v3)",vovv_aux_tensor,int_t%vovv)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_a  =",ierr
     ierr=talsh_tensor_contract("W(v1,o1,v2,v3)+=V(o2,o1,v2,v3)*t(v1,o2)",vovv_aux_tensor,int_t%oovv,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr f4b + f4d W_b  =",ierr

     ! f4b +f4b
     call f4bf4d(vovv_aux_tensor,Gtt_tensor%t1,l2_tensor,tXs,tYs,nocc,nvir)
     ! f3c + f4e
     call f4bf4d(vovv_aux_tensor,Gtt_tensor%t1,l2_tensor,tYs,tXs,nocc,nvir)


     ! now use the intermediates for i2b+i2c
     ierr=talsh_tensor_init(voov_aux_tensor)
     ierr=talsh_tensor_contract("W(f,m,i,c)+=W(f,m,a,c)*tX(a,i)",voov_aux_tensor,vovv_aux_tensor,tXs)
     if(ierr.ne.0) write(*,*) "ierr i2b+i2c _c  =",ierr
          ! contraction with tYs ck
     ierr=talsh_tensor_init(vooo_aux_tensor)
     ierr=talsh_tensor_contract("W(f,m,i,k)+=W(f,m,i,c)*tY(c,k)",vooo_aux_tensor,voov_aux_tensor,tYs)
     if(ierr.ne.0) write(*,*) "ierr i2b+i2c _d  =",ierr
     ierr=talsh_tensor_contract("G(m,e)+=W(f,m,i,k)*L(i,k,f,e)",Gtt_tensor%t1,vooo_aux_tensor,l2_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr i2b+i2c _e  =",ierr

     ! c9b, g7b, l29, c10b, k7b, l28  - inversion : c11b, k5b, l21, c12b, k10b , l21 bis
     ! c9b, g7b, l29, c10b, k7b, l28
     call Gtt_c9c10g7k7l28l29(Gtt_tensor%t1,int_t,tYs,tXd,l2_tensor,nocc,nvir)
     ! c11b, k5b, l21, c12b, k10b, l21 bis
     call Gtt_c9c10g7k7l28l29(Gtt_tensor%t1,int_t,tXs,tYd,l2_tensor,nocc,nvir)

     !!!!!!!!!!!!!!
     !!!! \Gd mnef
     !!!!!!!!!!!!!!

     ! l22bis do inversion for l32bis
     ! l23bis  !! SAME intermediate as l22bis   ! do invesrion for l34bis 
      
     call l22bisl23bis(Gtt_tensor%t2,l2_tensor,int_t,tXs,tYs,nocc,nvir)         
     call l22bisl23bis(Gtt_tensor%t2,l2_tensor,int_t,tYs,tXs,nocc,nvir)         


     ! i2d
     ierr=talsh_tensor_init(ooov_aux_tensor)
     ierr=talsh_tensor_contract("W(m,n,k,a)+=V(m,n,c,a)*tY(c,k)",ooov_aux_tensor,int_t%oovv,tYs)
     if(ierr.ne.0) write(*,*) "ierr i2d-a  =",ierr
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,n,k,i)+=W(m,n,k,a)*tX(a,i)",oooo_aux_tensor,ooov_aux_tensor,tXs)
     if(ierr.ne.0) write(*,*) "ierr i2d-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(m,n,k,i)*L(k,i,e,f)",Gtt_tensor%t2,oooo_aux_tensor,l2_tensor)


     ! n2d
     ierr=talsh_tensor_init(vvoo_aux_tensor)
     ierr=talsh_tensor_contract("W(a,c,i,k)+=tx(a,i)*tY(c,k)",vvoo_aux_tensor,tXs,tYs)
     if(ierr.ne.0) write(*,*) "ierr n2d-a  =",ierr
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("W(m,n,k,i)+=W(a,c,i,k)*L(m,n,c,a)",oooo_aux_tensor,vvoo_aux_tensor,l2_tensor)
     if(ierr.ne.0) write(*,*) "ierr i2d-b  =",ierr
     ierr=talsh_tensor_contract("G(m,n,e,f)+=W(m,n,k,i)*V(k,i,e,f)",Gtt_tensor%t2,oooo_aux_tensor,int_t%oovv)
     if(ierr.ne.0) write(*,*) "ierr i2d-c  =",ierr

     ! l28bis=inv(f4e)
     call l28bis(Gtt_tensor%t2,l2_tensor,int_t,tXs,tYs,nocc,nvir)
     call l28bis(Gtt_tensor%t2,l2_tensor,int_t,tYs,tXs,nocc,nvir) 

     ierr=talsh_tensor_destruct(G1_tensor)
     ierr=talsh_tensor_destruct(ooov_aux_tensor)
     ierr=talsh_tensor_destruct(ovov_aux_tensor)
     ierr=talsh_tensor_destruct(vvvo_aux_tensor)
     ierr=talsh_tensor_destruct(oooo_aux_tensor)
     ierr=talsh_tensor_destruct(vovv_aux_tensor)
     ierr=talsh_tensor_destruct(voov_aux_tensor)
     ierr=talsh_tensor_destruct(vooo_aux_tensor)
     ierr=talsh_tensor_destruct(vvoo_aux_tensor)
 
     end subroutine get_Gtt

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
