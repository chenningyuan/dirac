module talsh_sigma_eomccsd

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use talsh_common_routines, only : tensor_norm2, print_tensor
  use exacorr_utils,  only : print_date

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0),ONE=(1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)
  real(8),    parameter :: ONE_QUARTER=0.25D0

  private

  public talsh_sigma_eomccsd_rhs_ip
  public talsh_sigma_eomccsd_lhs_ip
  public talsh_sigma_eomccsd_rhs_ea
  public talsh_sigma_eomccsd_lhs_ea
  public talsh_sigma_eomccsd_rhs_ee
  public talsh_sigma_eomccsd_lhs_ee

  interface sigma_eomccsd_rhs_ip
     module procedure talsh_sigma_eomccsd_rhs_ip
  end interface

  interface sigma_eomccsd_lhs_ip
     module procedure talsh_sigma_eomccsd_lhs_ip
  end interface

  public sigma_eomccsd_rhs_ip
  public sigma_eomccsd_lhs_ip

  interface sigma_eomccsd_rhs_ea
     module procedure talsh_sigma_eomccsd_rhs_ea
  end interface

  interface sigma_eomccsd_lhs_ea
     module procedure talsh_sigma_eomccsd_lhs_ea
  end interface

  public sigma_eomccsd_rhs_ea
  public sigma_eomccsd_lhs_ea

  interface sigma_eomccsd_rhs_ee
     module procedure talsh_sigma_eomccsd_rhs_ee
  end interface

  interface sigma_eomccsd_lhs_ee
     module procedure talsh_sigma_eomccsd_lhs_ee
  end interface

  public sigma_eomccsd_rhs_ee
  public sigma_eomccsd_lhs_ee
  

  contains

    subroutine talsh_sigma_eomccsd_rhs_ip(int_t, interm_t, r1_tensor, r2_tensor, sigma1_tensor, sigma2_tensor)

!     input variables
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: r1_tensor
      type(talsh_tens_t),    intent(inout) :: r2_tensor
      type(talsh_tens_t),    intent(inout) :: sigma1_tensor
      type(talsh_tens_t),    intent(inout) :: sigma2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: v_aux
      
!     error code
      integer :: ierr
      real(kind=8) :: norm_s, norm_s_prev
      logical :: debug=.true.

!     test tensor
      type(talsh_tens_t) :: test_tensor


      ierr=talsh_tensor_construct(test_tensor,C8,(/interm_t%nvir,1,interm_t%nocc,interm_t%nocc/),init_val=ZERO)



!***************
!     r1_sigma *
!***************

      norm_s_prev = 0.0d0
!
! below X represents a continnum orbital, the size of that dimension will be 1 when defining the tensors
!--------------------------------------
!     term 1: S1(X,i) -= F(m,i) * R1(X,m) |
!--------------------------------------

      ierr=talsh_tensor_contract("S(X,i)+=F(m,i)*R(X,m)",sigma1_tensor,interm_t%fbar_oo,r1_tensor,scale=MINUS_ONE)
      norm_s = tensor_norm2(sigma1_tensor)
      if (debug) print*, "norm contraction 1 (accumulated, this):",norm_s, (norm_s - norm_s_prev) 
      norm_s_prev = norm_s

!------------------------------------------
!     term 2: S1(X,i) += F(m,e) * R2(e,X,i,m) |
!------------------------------------------

      ierr=talsh_tensor_contract("S(X,i)+=F(m,e)*R(e,X,i,m)",sigma1_tensor,interm_t%fbar_ov,r2_tensor)
      norm_s = tensor_norm2(sigma1_tensor)
      if (debug) print*, "norm contraction 2 (accumulated, this):",norm_s, (norm_s - norm_s_prev) 
      norm_s_prev = norm_s


!--------------------------------------------------
!     term 3: S1(X,i) += 1/2 W(m,n,i,e) * R2(e,X,m,n) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(X,i)+=W(m,n,i,e)*R(e,X,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=ONE_HALF)
      norm_s = tensor_norm2(sigma1_tensor)
      if (debug) print*, "norm contraction 3 (accumulated, this):",norm_s, (norm_s - norm_s_prev) 
      norm_s_prev = norm_s

!***************
!     r2_sigma *
!***************

!---------------------------------------------
!     term 1: S2(a,X,i,j) -= W(m,a,i,j) * R(X,m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,X,i,j)+=W(m,a,j,i)*R(X,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE)

      !Test
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(m,a,j,i)*R(X,m)",test_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 4 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm4_test")
      ! if (debug) call print_tensor(sigma2_tensor, 1.0d-10, "diagm4_sigma2")
      norm_s_prev = norm_s

!---------------------------------------------------------
!     term 2: S2(a,X,i,j) -= P(i,j) F(m,i) * R2(a,X,m,j) |
!---------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=F(m,i)*r(a,X,j,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)!Used
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=F(m,j)*r(a,X,m,i)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)


      ierr=talsh_tensor_contract("S(a,X,i,j)+=F(m,i)*r(a,X,m,j)",sigma2_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)                !XY
      ierr=talsh_tensor_contract("S(a,X,i,j)+=F(m,j)*r(a,X,m,i)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)!XY

      ! Test
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=F(m,i)*r(a,X,j,m)",test_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=F(m,j)*r(a,X,i,m)",test_tensor,interm_t%fbar_oo,r2_tensor)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 5 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm5")
      norm_s_prev = norm_s


!------------------------------------------------------
!     term 3: S2(a,X,i,j) += W(m,n,i,j) * R2(a,X,m,n) |
!------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(m,n,i,j)*r(a,X,m,n)",sigma2_tensor,interm_t%w_oooo,r2_tensor)

      ierr=talsh_tensor_contract("S(a,X,i,j)+=W(m,n,i,j)*r(a,X,m,n)",sigma2_tensor,interm_t%w_oooo,r2_tensor,scale=ONE_HALF) !XY
      !Test
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(m,n,i,j)*r(a,X,m,n)",test_tensor,interm_t%w_oooo,r2_tensor,scale=ONE_HALF)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 6 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm6")
      norm_s_prev = norm_s


!-------------------------------------------------------------
!     term 4: S2(a,X,i,j) -= P(i,j) W(m,a,i,e) * R2(e,X,m,j) |
!-------------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(a,m,e,i)*R(e,X,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)!Used
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(a,m,e,j)*R(e,X,m,i)",sigma2_tensor,interm_t%w_vovo,r2_tensor)


      ierr=talsh_tensor_contract("S(a,X,i,j)+=W(a,m,e,i)*R(e,X,m,j)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)!XY
      ierr=talsh_tensor_contract("S(a,X,i,j)+=W(a,m,e,j)*R(e,X,m,i)",sigma2_tensor,interm_t%w_vovo,r2_tensor)!XY


      ! Test
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(a,m,e,i)*R(e,X,j,m)",test_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=W(a,m,e,j)*R(e,X,i,m)",test_tensor,interm_t%w_vovo,r2_tensor)


      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 7 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm7")
      ! if (debug) call print_tensor(interm_t%w_vovo, 1.0d-10, "w_ovov")
      norm_s_prev = norm_s

!--------------------------------------------------
!     term 5: S2(a,X,i,j) += F(a,e) * R2(e,X,i,j) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a,X,i,j)+=F(a,e)*R(e,X,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)

      ! Test
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=F(a,e)*R(e,X,i,j)",test_tensor,interm_t%fbar_vv,r2_tensor)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 8 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm8")
      norm_s_prev = norm_s


!---------------------------------------------------------------------------
!     term 6: S2(a,X,i,j) += 1/2 ( V(m,n,e,f) * R2(f,X,m,n) ) * T(e,a,i,j) |
!---------------------------------------------------------------------------

      ierr=talsh_tensor_construct(v_aux,C8,(/1,interm_t%nvir/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(X,e)+=V(m,n,e,f)*R(f,X,m,n)",v_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=U(X,e)*T(e,a,i,j)",sigma2_tensor,v_aux,interm_t%t2)

      ierr=talsh_tensor_contract("S(a,X,i,j)+=U(X,e)*T(e,a,i,j)",sigma2_tensor,v_aux,interm_t%t2,scale=MINUS_ONE)  !XY

      !Test
      ! ierr=talsh_tensor_contract("S(a,X,i,j)+=U(X,e)*T(e,a,i,j)",test_tensor,v_aux,interm_t%t2)

      ierr=talsh_tensor_destruct(v_aux)

      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 9 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm9")
      norm_s_prev = norm_s


      ! debug
      if (debug) then
        print*, "norm r1_sigma", tensor_norm2(sigma1_tensor)
        print*, "norm r2_sigma", tensor_norm2(sigma2_tensor),tensor_norm2(sigma2_tensor)*ONE_QUARTER
        print*, "norm r_sigma ", (tensor_norm2(sigma1_tensor) + tensor_norm2(sigma2_tensor))
      end if

    end subroutine talsh_sigma_eomccsd_rhs_ip


    subroutine talsh_sigma_eomccsd_lhs_ip(int_t,interm_t,l1_tensor, l2_tensor, sigma1_tensor, sigma2_tensor)

!     input variables
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: l1_tensor
      type(talsh_tens_t),    intent(inout) :: l2_tensor
      type(talsh_tens_t),    intent(inout) :: sigma1_tensor
      type(talsh_tens_t),    intent(inout) :: sigma2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: v_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(sigma1_tensor)

!-----------------------------------------
!     term 1: S1(i,X) -= F(i,m) * L(m,X) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(i,X)+=F(i,m)*L(m,X)",sigma1_tensor,interm_t%fbar_oo,l1_tensor,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 2: S1(i,X) += 1/2 W(i,e,m,n) * L(m,n,e,X) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,X)+=W(i,e,m,n)*L(m,n,e,X)",sigma1_tensor,interm_t%w_ovoo,l2_tensor,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(sigma2_tensor)

!---------------------------------------------
!     term 1: S2(i,j,a,X) -= W(i,j,m,a) * L(m,X) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,X)+=W(j,i,m,a)*L(m,X)",sigma2_tensor,interm_t%w_ooov,l1_tensor)

!---------------------------------------------
!     term 2: S2(i,j,a,X) += F(e,a) * L(i,j,e,X) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=F(e,a)*L(i,j,e)",sigma2_tensor,interm_t%fbar_vv,l2_tensor)

!----------------------------------------------------
!     term 3: S2(i,j,a,X) -= P(i,j) F(j,m) * L(i,m,a,X) |
!----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,X)+=F(j,m)*L(m,i,a,X)",sigma2_tensor,interm_t%fbar_oo,l2_tensor)
      ierr=talsh_tensor_contract("S(i,j,a,X)+=F(i,m)*L(j,m,a,X)",sigma2_tensor,interm_t%fbar_oo,l2_tensor)

!--------------------------------------------------------
!     term 4: S2(i,j,a,X) += P(i,j) W(j,e,a,m) * L(i,m,e,X) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,X)+=W(e,j,a,m)*L(m,i,e,X)",sigma2_tensor,interm_t%w_vovo,l2_tensor)
      ierr=talsh_tensor_contract("S(i,j,a,X)+=W(e,i,a,m)*L(j,m,e,X)",sigma2_tensor,interm_t%w_vovo,l2_tensor)

!-----------------------------------------------------
!     term 5: S2(i,j,a,X) += 1/2 W(i,j,m,n) * L(m,n,a,X) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,X)+=W(i,j,m,n)*L(m,n,a,X)",sigma2_tensor,interm_t%w_oooo,l2_tensor,scale=ONE_HALF)

!----------------------------------------------------------------------
!     term 6: S2(i,j,a,X) += 1/2 V(i,j,e,a) * ( L(m,n,f,X) * T(e,f,m,n) ) |
!----------------------------------------------------------------------

      ierr=talsh_tensor_construct(v_aux,C8,(/interm_t%nvir, 1/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(e,X)+=L(m,n,f,X)*T(e,f,m,n)",v_aux,l2_tensor,interm_t%t2,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(i,j,a,X)+=V(i,j,e,a)*U(e,X)",sigma2_tensor,int_t%oovv,v_aux)
      ierr=talsh_tensor_destruct(v_aux)

!------------------------------------------------
!     term 7: S2(i,j,a,X) += P(i,j) L(i,X) * F(j,a) |
!------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,X)+=L(i,X)*F(j,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
      ierr=talsh_tensor_contract("S(i,j,a,X)+=L(j,X)*F(i,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)

    end subroutine talsh_sigma_eomccsd_lhs_ip

    subroutine talsh_sigma_eomccsd_rhs_ea(int_t, interm_t, r1_tensor, r2_tensor, sigma1_tensor, sigma2_tensor)
!     input variables
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: sigma1_tensor
      type(talsh_tens_t),    intent(inout) :: sigma2_tensor
      type(talsh_tens_t),    intent(inout) :: r1_tensor
      type(talsh_tens_t),    intent(inout) :: r2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: o_aux
!     error code
      integer :: ierr

!***************
!     r1_sigma *
!***************

      ierr=talsh_tensor_init(sigma1_tensor)

!--------------------------------------
!     term 1: S1(a,X) += F(e,a) * R1(e,X) |
!--------------------------------------

      ierr=talsh_tensor_contract("S(a,X)+=F(e,a)*R(e,X)",sigma1_tensor,interm_t%fbar_vv,r1_tensor)

!------------------------------------------
!     term 2: S1(a,X) += F(m,e) * R2(e,a,m,X) |
!------------------------------------------

      ierr=talsh_tensor_contract("S(a,X)+=F(m,e)*R(e,a,m,X)",sigma1_tensor,interm_t%fbar_ov,r2_tensor)

!--------------------------------------------------
!     term 3: S1(a,X) += 1/2 W(a,m,f,e) * R2(e,f,m,X) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a,X)+=W(a,m,f,e)*R(e,f,m,X)",sigma1_tensor,interm_t%w_vovv,r2_tensor,scale=ONE_HALF)

!***************
!     r2_sigma *
!***************

      ierr=talsh_tensor_init(sigma2_tensor)

!---------------------------------------------
!     term 1: S2(a,b,i,X) += W(a,b,e,i) * R(e,X) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,X)+=W(a,b,e,i)*R(e,X)",sigma2_tensor,interm_t%w_vvvo,r1_tensor)

!-----------------------------------------------------
!     term 2: S2(a,b,i,X) += P(a,b) F(a,e) * R2(e,b,i,X) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,X)+=F(a,e)*r(e,b,i,X)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)
      ierr=talsh_tensor_contract("S(a,b,i,X)+=F(b,e)*r(a,e,i,X)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)

!--------------------------------------------------
!     term 3: S2(a,b,i,X) += W(a,b,e,f) * R2(e,f,i,X) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,X)+=W(a,b,e,f)*r(e,f,i,X)",sigma2_tensor,interm_t%w_vvvv,r2_tensor,scale=ONE_HALF)

!---------------------------------------------------------
!     term 4: S2(a,b,i,X) -= P(a,b) W(a,m,e,i) * R2(e,b,m,X) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,X)+=W(a,m,e,i)*R(b,e,m,X)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ierr=talsh_tensor_contract("S(a,b,i,X)+=W(b,m,e,i)*R(e,a,m,X)",sigma2_tensor,interm_t%w_vovo,r2_tensor)

!----------------------------------------------
!     term 5: S2(a,b,i,X) -= F(m,i) * R2(a,b,m,X) |
!----------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,X)+=F(m,i)*R(b,a,m,X)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)

!-----------------------------------------------------------------------
!     term 6: S2(a,b,i,X) += 1/2 ( V(m,n,e,f) * R2(e,f,n,X) ) * T(a,b,m,i) |
!-----------------------------------------------------------------------

      ierr=talsh_tensor_construct(o_aux,C8,(/interm_t%nocc, 1/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(m,X)+=V(m,n,e,f)*R(e,f,n,X)",o_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(a,b,i,X)+=U(m,X)*T(a,b,m,i)",sigma2_tensor,o_aux,interm_t%t2)
      ierr=talsh_tensor_destruct(o_aux)

      ! debug
      if (.true.) then
        print*, "r1_sigma", tensor_norm2(sigma1_tensor)
        print*, "r2_sigma", tensor_norm2(sigma2_tensor)*ONE_QUARTER
      end if

    end subroutine talsh_sigma_eomccsd_rhs_ea

    subroutine talsh_sigma_eomccsd_lhs_ea(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)

!     input variables
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: sigma1_tensor
      type(talsh_tens_t),    intent(inout) :: sigma2_tensor
      type(talsh_tens_t),    intent(inout) :: l1_tensor
      type(talsh_tens_t),    intent(inout) :: l2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: o_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(sigma1_tensor)

!-------------------------------------
!     term 1: S1(X,a) += F(e,a) * L(X,e) |
!-------------------------------------

      ierr=talsh_tensor_contract("S(X,a)+=F(e,a)*L(X,e)",sigma1_tensor,interm_t%fbar_vv,l1_tensor)

!-------------------------------------------------
!     term 2: S1(X,a) -= 1/2 W(e,f,a,m) * L(m,X,f,e) |
!-------------------------------------------------

      ierr=talsh_tensor_contract("S(X,a)+=W(e,f,a,m)*L(m,X,e,f)",sigma1_tensor,interm_t%w_vvvo,l2_tensor,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(sigma2_tensor)

!---------------------------------------------
!     term 1: S2(i,X,a,b) += W(e,i,b,a) * L(X,e) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,X,a,b)+=W(e,i,b,a)*L(X,e)",sigma2_tensor,interm_t%w_vovv,l1_tensor)

!----------------------------------------------------
!     term 2: S2(i,X,a,b) += P(a,b) F(e,a) * L(i,X,e,b) |
!----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,X,a,b)+=F(e,a)*L(i,X,e,b)",sigma2_tensor,interm_t%fbar_vv,l2_tensor)
      ierr=talsh_tensor_contract("S(i,X,a,b)+=F(e,b)*L(i,X,a,e)",sigma2_tensor,interm_t%fbar_vv,l2_tensor)

!---------------------------------------------
!     term 3: S2(i,X,a,b) -= F(i,m) * L(m,X,a,b) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,X,a,b)+=F(i,m)*L(m,X,b,a)",sigma2_tensor,interm_t%fbar_oo,l2_tensor)

!-----------------------------------------------------
!     term 4: S2(i,X,a,b) += 1/2 W(e,f,a,b) * L(i,X,e,f) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,X,a,b)+=W(e,f,a,b)*L(i,X,e,f)",sigma2_tensor,interm_t%w_vvvv,l2_tensor,scale=ONE_HALF)

!--------------------------------------------------------
!     term 5: S2(i,X,a,b) -= P(a,b) W(e,i,b,m) * L(X,m,a,e) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,X,a,b)+=W(e,i,b,m)*L(X,m,e,a)",sigma2_tensor,interm_t%w_vovo,l2_tensor)
      ierr=talsh_tensor_contract("S(i,X,a,b)+=W(e,i,a,m)*L(X,m,b,e)",sigma2_tensor,interm_t%w_vovo,l2_tensor)

!----------------------------------------------------------------------
!     term 6: S2(i,X,a,b) -= 1/2 ( L(X,o,e,f) * T(e,f,m,o) ) * V(m,i,b,a) |
!----------------------------------------------------------------------

      ierr=talsh_tensor_construct(o_aux,C8,(/1, interm_t%nocc/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(X,m)+=L(X,o,e,f)*T(e,f,o,m)",o_aux,l2_tensor,interm_t%t2,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(i,X,a,b)+=U(X,m)*V(m,i,b,a)",sigma2_tensor,o_aux,int_t%oovv)
      ierr=talsh_tensor_destruct(o_aux)

!------------------------------------------------
!     term 7: S2(i,X,a,b) += P(a,b) F(i,a) * L(X,b) |
!------------------------------------------------

      ierr=talsh_tensor_contract("S(i,X,a,b)+=F(i,a)*L(X,b)",sigma2_tensor,interm_t%fbar_ov,l1_tensor)
      ierr=talsh_tensor_contract("S(i,X,a,b)+=F(i,b)*L(X,a)",sigma2_tensor,interm_t%fbar_ov,l1_tensor,scale=MINUS_ONE)

    end subroutine talsh_sigma_eomccsd_lhs_ea

    subroutine talsh_sigma_eomccsd_rhs_ee(int_t,interm_t,r1_tensor,r2_tensor,sigma1_tensor,sigma2_tensor)

!     input variables
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: sigma1_tensor
      type(talsh_tens_t),    intent(inout) :: sigma2_tensor
      type(talsh_tens_t),    intent(inout) :: r1_tensor
      type(talsh_tens_t),    intent(inout) :: r2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: oo_aux, vv_aux
      integer            :: oo_dims(2), vv_dims(2)
!     error code
      integer :: ierr

      ! type(talsh_tens_t) :: test_tensor

      oo_dims = (/interm_t%nocc,interm_t%nocc/)
      vv_dims = (/interm_t%nvir,interm_t%nvir/)

      ! ierr=talsh_tensor_construct(test_tensor,C8,(/interm_t%nvir,interm_t%nvir,interm_t%nocc,interm_t%nocc/),init_val=ZERO)

!***************
!     r1_sigma *
!***************

      ierr=talsh_tensor_init(sigma1_tensor)

!-----------------------------------------
!     term 1: S1(a,i) += F(a,e) * R(e,i) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=F(a,e)*R(e,i)",sigma1_tensor,interm_t%fbar_vv,r1_tensor)

!-----------------------------------------
!     term 2: S1(a,i) -= F(m,i) * R(a,m) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=F(m,i)*R(a,m)",sigma1_tensor,interm_t%fbar_oo,r1_tensor,scale=MINUS_ONE)

!---------------------------------------------
!     term 3: S1(a,i) += F(m,e) * R(e,a,m,i) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=F(m,e)*R(e,a,m,i)",sigma1_tensor,interm_t%fbar_ov,r2_tensor)

!---------------------------------------------
!     term 4: S1(a,i) += W(m,a,e,i) * R(e,m) |
!---------------------------------------------
      ! Flip the index since we only have w_vovo
      ierr=talsh_tensor_contract("S(a,i)+=W(a,m,e,i)*R(e,m)",sigma1_tensor,interm_t%w_vovo,r1_tensor,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 5: S1(a,i) += 1/2 W(a,m,e,f) * R(e,f,i,m) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=W(a,m,e,f)*R(e,f,i,m)",sigma1_tensor,interm_t%w_vovv,r2_tensor,scale=ONE_HALF)

!-----------------------------------------------------
!     term 6: S1(a,i) += 1/2 W(m,n,i,e) * R(e,a,m,n) |
!-----------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,i)+=W(m,n,i,e)*R(a,e,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=ONE_HALF)

      ierr=talsh_tensor_contract("S(a,i)+=W(m,n,i,e)*R(e,a,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=ONE_HALF) !XY

!***************
!     r2_sigma *
!***************

      ierr=talsh_tensor_init(sigma2_tensor)

!--------------------------------------------------------
!     term 1: S2(a,b,i,j) -= P(a,b) W(m,b,i,j) * R(a,m) |
!--------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,b,j,i)*R(a,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,a,i,j)*R(b,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor)

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,b,i,j)*R(a,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE)     !XY
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,a,i,j)*R(b,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor)     !XY

      ! call print_tensor(sigma2_tensor, 1.0d-10, "diagm7")


!--------------------------------------------------------
!     term 2: S2(a,b,i,j) += P(i,j) W(a,b,e,j) * R(e,i) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,b,e,j)*R(e,i)",sigma2_tensor,interm_t%w_vvvo,r1_tensor)  !Used
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,a,e,i)*R(e,j)",sigma2_tensor,interm_t%w_vvvo,r1_tensor)

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,b,e,i)*R(e,j)",sigma2_tensor,interm_t%w_vvvo,r1_tensor,scale=MINUS_ONE)   !XY
      call print_tensor(interm_t%w_vvvo, 1.0d-10, "w_vvvo")
      call print_tensor(sigma2_tensor, 1.0d-10, "diagm8")

!-------------------------------------------------------------------------
!     term 3: S2(a,b,i,j) += P(a,b) ( W(b,m,f,e) * R(e,m) ) * T(a,f,i,j) |
!-------------------------------------------------------------------------

      ierr=talsh_tensor_construct(vv_aux,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(b,f)+=W(b,m,f,e)*R(e,m)",vv_aux,interm_t%w_vovv,r1_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(b,f)*T2(a,f,i,j)",sigma2_tensor,vv_aux,interm_t%t2)     !Used
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(a,f)*T2(f,b,i,j)",sigma2_tensor,vv_aux,interm_t%t2)

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(a,f)*T2(b,f,i,j)",sigma2_tensor,vv_aux,interm_t%t2,scale=MINUS_ONE)      !XY
      ierr=talsh_tensor_destruct(vv_aux)

!-------------------------------------------------------------------------
!     term 4: S2(a,b,i,j) -= P(i,j) ( W(n,m,j,e) * R(e,m) ) * T(a,b,i,n) |
!-------------------------------------------------------------------------

      ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(n,j)+=W(n,m,j,e)*R(e,m)",oo_aux,interm_t%w_ooov,r1_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,j)*T(a,b,i,n)",sigma2_tensor,oo_aux,interm_t%t2)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,i)*T(a,b,n,j)",sigma2_tensor,oo_aux,interm_t%t2)


      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,j)*T(a,b,i,n)",sigma2_tensor,oo_aux,interm_t%t2,scale=MINUS_ONE)        !XY
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,i)*T(a,b,j,n)",sigma2_tensor,oo_aux,interm_t%t2)        !XY
      ierr=talsh_tensor_destruct(oo_aux)

!--------------------------------------------------------
!     term 5: S2(a,b,i,j) += P(a,b) F(b,e) * R(a,e,i,j) |
!--------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(b,e)*R(a,e,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)   !Used
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(a,e)*R(e,b,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(a,e)*R(b,e,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor,scale=MINUS_ONE)   !XY


!--------------------------------------------------------
!     term 6: S2(a,b,i,j) -= P(i,j) F(m,j) * R(a,b,i,m) |
!--------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,j)*R(a,b,m,i)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,i)*R(a,b,j,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)


      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,j)*R(a,b,i,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,i)*R(a,b,j,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)
!---------------------------------------------------------
!     term 7: S2(a,b,i,j) += 1/2 W(m,n,i,j) * R(a,b,m,n) |
!---------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,n,i,j)*R(a,b,m,n)",sigma2_tensor,interm_t%w_oooo,r2_tensor,scale=ONE_HALF)

!-------------------------------------------------------------------
!     term 8: S2(a,b,i,j) += P(i,j) P(a,b) W(m,b,e,j) * R(a,e,i,m) |
!-------------------------------------------------------------------


      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,m,e,j)*R(a,e,m,i)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,m,e,j)*R(b,e,i,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,m,e,i)*R(a,e,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,m,e,i)*R(b,e,m,j)",sigma2_tensor,interm_t%w_vovo,r2_tensor)


      ! Filp the index since we only have w_vovo, XY
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,m,e,j)*R(a,e,i,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,m,e,j)*R(b,e,i,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,m,e,i)*R(a,e,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,m,e,i)*R(b,e,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)

!---------------------------------------------------------------------------------
!     term 9: S2(a,b,i,j) -= 1/2 P(a,b) ( V(n,m,f,e) * R(e,a,m,n) ) * T(f,b,i,j) |
!---------------------------------------------------------------------------------

      ierr=talsh_tensor_construct(vv_aux,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(a,f)+=V(n,m,f,e)*R(e,a,m,n)",vv_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(a,f)*T(f,b,j,i)",sigma2_tensor,vv_aux,interm_t%t2)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(b,f)*T(f,a,i,j)",sigma2_tensor,vv_aux,interm_t%t2)


      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(a,f)*T(f,b,i,j)",sigma2_tensor,vv_aux,interm_t%t2,scale=MINUS_ONE)   !XY
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(b,f)*T(f,a,i,j)",sigma2_tensor,vv_aux,interm_t%t2)                   !XY
      ierr=talsh_tensor_destruct(vv_aux)

!---------------------------------------------------------------------------------
!     term 10: S2(a,b,i,j) -= 1/2 P(i,j) ( V(n,m,f,e) * R(f,e,i,m) ) * T(b,a,j,n) |
!---------------------------------------------------------------------------------

      ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(n,i)+=V(n,m,f,e)*R(f,e,i,m)",oo_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,i)*T(b,a,n,j)",sigma2_tensor,oo_aux,interm_t%t2)
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,j)*T(b,a,i,n)",sigma2_tensor,oo_aux,interm_t%t2)


      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,i)*T(b,a,n,j)",sigma2_tensor,oo_aux,interm_t%t2)                   !XY
      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,j)*T(b,a,n,i)",sigma2_tensor,oo_aux,interm_t%t2,scale=MINUS_ONE)   !XY
      ierr=talsh_tensor_destruct(oo_aux)

!----------------------------------------------------------
!     term 11: S2(a,b,i,j) += 1/2 W(a,b,e,f) * R(e,f,i,j) |
!----------------------------------------------------------

      ! ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,b,e,f)*R(e,f,i,j)",sigma2_tensor,interm_t%w_vvvv,r2_tensor,scale=ONE_HALF)

    end subroutine talsh_sigma_eomccsd_rhs_ee

    subroutine talsh_sigma_eomccsd_lhs_ee(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)

!     input variables
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_tens_t),    intent(inout) :: sigma1_tensor
      type(talsh_tens_t),    intent(inout) :: sigma2_tensor
      type(talsh_tens_t),    intent(inout) :: l1_tensor
      type(talsh_tens_t),    intent(inout) :: l2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: oo_aux
      integer            :: oo_dims(2)
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(sigma1_tensor)

!--------------------------------
!     term 1: S1(i,a) += F(i,a) |
!--------------------------------

      ierr=talsh_tensor_add("S(i,a)+=F(i,a)",sigma1_tensor,interm_t%fbar_ov)

!-----------------------------------------
!     term 2: S1(i,a) += L(i,e) * F(e,a) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",sigma1_tensor,l1_tensor,interm_t%fbar_vv)

!-----------------------------------------
!     term 3: S1(i,a) -= L(m,a) * F(i,m) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",sigma1_tensor,l1_tensor,interm_t%fbar_oo,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 4: S1(i,a) += 1/2 L(i,m,e,f) * W(e,f,a,m) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",sigma1_tensor,l2_tensor,interm_t%w_vvvo,scale=ONE_HALF)

!---------------------------------------------
!     term 5: S1(i,a) -= G(f,e) * W(e,i,f,a) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",sigma1_tensor,interm_t%gvv,interm_t%w_vovv)

!---------------------------------------------
!     term 6: S1(i,a) -= G(n,m) * W(m,i,n,a) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",sigma1_tensor,interm_t%goo,interm_t%w_ooov)

!---------------------------------------------
!     term 7: S1(i,a) += L(m,e) * W(i,e,a,m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",sigma1_tensor,l1_tensor,interm_t%w_vovo,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 8: S1(i,a) -= 1/2 L(m,n,a,e) * W(i,e,m,n) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",sigma1_tensor,l2_tensor,interm_t%w_ovoo,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(sigma2_tensor)

!----------------------------------------
!     term 1: S2(i,j,a,b) += V(i,j,a,b) |
!----------------------------------------

      ierr=talsh_tensor_add("S(i,j,a,b)+=V(i,j,a,b)",sigma2_tensor,int_t%oovv)

!--------------------------------------------------------
!     term 2: S2(i,j,a,b) += P(a,b) L(i,j,a,e) * F(e,b) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)

!--------------------------------------------------------
!     term 3: S2(i,j,a,b) -= P(i,j) L(i,m,a,b) * F(j,m) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)

!---------------------------------------------------------
!     term 4: S2(i,j,a,b) += 1/2 L(m,n,a,b) * W(i,j,m,n) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",sigma2_tensor,l2_tensor,interm_t%w_oooo,scale=ONE_HALF)

!-------------------------------------------------------------------
!     term 5: S2(i,j,a,b) += P(i,j) P(a,b) L(i,m,a,e) * W(j,e,b,m) |
!-------------------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)

!--------------------------------------------------------
!     term 6: S2(i,j,a,b) += P(a,b) V(i,j,a,e) * G(e,b) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",sigma2_tensor,int_t%oovv,interm_t%gvv)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",sigma2_tensor,int_t%oovv,interm_t%gvv)

!-------------------------------------------------
!     term 7: S2(i,j,a,b) -= L(m,a) * W(i,j,m,b) |
!-------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",sigma2_tensor,l1_tensor,interm_t%w_ooov)

!--------------------------------------------------------
!     term 8: S2(i,j,a,b) -= P(i,j) V(i,m,a,b) * G(j,m) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,i,a,b)*G(j,m)",sigma2_tensor,int_t%oovv,interm_t%goo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",sigma2_tensor,int_t%oovv,interm_t%goo)

!---------------------------------------------------------------------
!     term 9: S2(i,j,a,b) += P(i,j) V(m,j,a,b) * ( L(i,e) * T(e,m) ) |
!---------------------------------------------------------------------

      if (.not.interm_t%ccd) then
        oo_dims = (/interm_t%nocc,interm_t%nocc/)
        ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
        ierr=talsh_tensor_contract("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,j,a,b)*U(i,m)",sigma2_tensor,int_t%oovv,oo_aux)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,m,a,b)*U(j,m)",sigma2_tensor,int_t%oovv,oo_aux)
        ierr=talsh_tensor_destruct(oo_aux)
      end if

!---------------------------------------------------------
!     term 10: S2(i,j,a,b) += P(i,j) L(i,e) * V(e,j,a,b) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*V(e,j,a,b)",sigma2_tensor,l1_tensor,int_t%vovv)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*V(e,i,b,a)",sigma2_tensor,l1_tensor,int_t%vovv)

!------------------------------------------------------------
!     term 11: S2(i,j,a,b) += P(i,j) P(a,b) L(i,a) * F(j,b) |
!------------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)

!--------------------------------------------------------
!     term 12: S2(i,j,a,b) += 1/2 L(i,j,e,f) W(e,f,a,b) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",sigma2_tensor,l2_tensor,interm_t%w_vvvv,scale=ONE_HALF)

    end subroutine talsh_sigma_eomccsd_lhs_ee

end module talsh_sigma_eomccsd
