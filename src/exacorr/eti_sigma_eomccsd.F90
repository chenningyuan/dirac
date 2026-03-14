module eti_sigma_eomccsd

  use exacorr_tensor_interface
  use tensor_algebra
  use exacorr_datatypes
  use exacorr_scaling_parameters
  use eti_common_routines, only : tensor_norm2, print_tensor
 
  implicit none

  private

  public eti_sigma_eomccsd_rhs_ip
  public eti_sigma_eomccsd_lhs_ip
  public eti_sigma_eomccsd_rhs_ea
  public eti_sigma_eomccsd_lhs_ea
  public eti_sigma_eomccsd_rhs_ee
  public eti_sigma_eomccsd_lhs_ee

  interface sigma_eomccsd_rhs_ip
     module procedure eti_sigma_eomccsd_rhs_ip
  end interface

  interface sigma_eomccsd_lhs_ip
     module procedure eti_sigma_eomccsd_lhs_ip
  end interface

  public sigma_eomccsd_rhs_ip
  public sigma_eomccsd_lhs_ip

  interface sigma_eomccsd_rhs_ea
     module procedure eti_sigma_eomccsd_rhs_ea
  end interface

  interface sigma_eomccsd_lhs_ea
     module procedure eti_sigma_eomccsd_lhs_ea
  end interface

  public sigma_eomccsd_rhs_ea
  public sigma_eomccsd_lhs_ea

  interface sigma_eomccsd_rhs_ee
     module procedure eti_sigma_eomccsd_rhs_ee
  end interface

  interface sigma_eomccsd_lhs_ee
     module procedure eti_sigma_eomccsd_lhs_ee
  end interface


  public sigma_eomccsd_rhs_ee
  public sigma_eomccsd_lhs_ee
  

  contains

    subroutine eti_sigma_eomccsd_rhs_ip(int_t, interm_t, r1_tensor, r2_tensor, sigma1_tensor, sigma2_tensor)

!     input variables
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(tensor_t),    intent(inout) :: r1_tensor
      type(tensor_t),    intent(inout) :: r2_tensor
      type(tensor_t),    intent(inout) :: sigma1_tensor
      type(tensor_t),    intent(inout) :: sigma2_tensor
!     auxiliary tensor
      type(tensor_t) :: v_aux
      
!     error code
      integer :: ierr
      real(kind=8) :: norm_s, norm_s_prev
      logical :: debug=.false.

!     test tensor
      type(tensor_t) :: test_tensor, test1_tensor


      ierr=construct_tensor_t(test_tensor,C8,(/interm_t%nvir,1,interm_t%nocc,interm_t%nocc/),init_val=ZERO)
      ierr=construct_tensor_t(test1_tensor,C8,(/1,interm_t%nocc/),init_val=ZERO)



!***************
!     r1_sigma *
!***************

      norm_s_prev = 0.0d0
!
! below X represents a continnum orbital, the size of that dimension will be 1 when defining the tensors
!--------------------------------------
!     term 1: S1(X,i) -= F(m,i) * R1(X,m) |
!--------------------------------------

      ierr=t_einsum("S(X,i)+=F(m,i)*R(X,m)",sigma1_tensor,interm_t%fbar_oo,r1_tensor,scale=MINUS_ONE)
      norm_s = tensor_norm2(sigma1_tensor)
      if (debug) print*, "norm contraction 1 (accumulated, this):",norm_s, (norm_s - norm_s_prev) 
      norm_s_prev = norm_s

!------------------------------------------
!     term 2: S1(X,i) += F(m,e) * R2(e,X,i,m) |
!------------------------------------------

      ierr=t_einsum("S(X,i)+=F(m,e)*R(e,X,i,m)",sigma1_tensor,interm_t%fbar_ov,r2_tensor)   !diam2

      !Test
      ! ierr=t_einsum("S(X,i)+=F(m,e)*R(e,X,i,m)",test1_tensor,interm_t%fbar_ov,r2_tensor)


      norm_s = tensor_norm2(sigma1_tensor)
      if (debug) print*, "norm contraction 2 (accumulated, this):",norm_s, (norm_s - norm_s_prev) 
      ! call print_tensor(interm_t%fbar_ov, 1.0d-10, "fbar_ov")

      ! call print_tensor(test1_tensor, 1.0d-10, "diagm2_sigma1")
      norm_s_prev = norm_s


!--------------------------------------------------
!     term 3: S1(X,i) += 1/2 W(m,n,i,e) * R2(e,X,m,n) |
!--------------------------------------------------

      ! ierr=t_einsum("S(X,i)+=W(m,n,i,e)*R(e,X,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=ONE_HALF) !diam3

      ierr=t_einsum("S(X,i)+=W(m,n,i,e)*R(e,X,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=MINUS_ONE_HALF) !  XYdiam3

      norm_s = tensor_norm2(sigma1_tensor)

      !Test
      ! ierr=t_einsum("S(X,i)+=W(m,n,i,e)*R(e,X,m,n)",test1_tensor,interm_t%w_ooov,r2_tensor,scale=MINUS_ONE)
      if (debug) print*, "norm contraction 3 (accumulated, this):",norm_s, (norm_s - norm_s_prev) 
      ! call print_tensor(test1_tensor, 1.0d-10, "diagm3_sigma1")
      norm_s_prev = norm_s

!***************
!     r2_sigma *
!***************

!---------------------------------------------
!     term 1: S2(a,X,i,j) -= W(m,a,i,j) * R(X,m) |
!---------------------------------------------

      ! ierr=t_einsum("S(a,X,i,j)+=W(m,a,j,i)*R(X,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE) !diagm4


      ierr=t_einsum("S(a,X,i,j)+=W(m,a,i,j)*R(X,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE) !diagm4 XY
      ! !Test
      ! ierr=t_einsum("S(a,X,i,j)+=W(m,a,i,j)*R(X,m)",test_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE_HALF)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 4 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm4_test")
      ! if (debug) call print_tensor(sigma2_tensor, 1.0d-10, "diagm4_sigma2")
      norm_s_prev = norm_s

!---------------------------------------------------------
!     term 2: S2(a,X,i,j) -= P(i,j) F(m,i) * R2(a,X,m,j) |
!---------------------------------------------------------

      ! ierr=t_einsum("S(a,X,i,j)+=F(m,i)*r(a,X,j,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)!Used
      ! ierr=t_einsum("S(a,X,i,j)+=F(m,j)*r(a,X,m,i)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)


      ierr=t_einsum("S(a,X,i,j)+=F(m,i)*r(a,X,m,j)",sigma2_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)                !XY !diam5
      ierr=t_einsum("S(a,X,i,j)+=F(m,j)*r(a,X,m,i)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)!XY !diam5 

      ! Test
      ! ierr=t_einsum("S(a,X,i,j)+=F(m,i)*r(a,X,m,j)",test_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)
      ! ierr=t_einsum("S(a,X,i,j)+=F(m,j)*r(a,X,m,i)",test_tensor,interm_t%fbar_oo,r2_tensor)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 5 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm5")
      norm_s_prev = norm_s


!------------------------------------------------------
!     term 3: S2(a,X,i,j) += W(m,n,i,j) * R2(a,X,m,n) |
!------------------------------------------------------

      ! ierr=t_einsum("S(a,X,i,j)+=W(m,n,i,j)*r(a,X,m,n)",sigma2_tensor,interm_t%w_oooo,r2_tensor)

      ierr=t_einsum("S(a,X,i,j)+=W(m,n,i,j)*r(a,X,m,n)",sigma2_tensor,interm_t%w_oooo,r2_tensor,scale=ONE_HALF) !XY!diam6
      !Test
      ! ierr=t_einsum("S(a,X,i,j)+=W(m,n,i,j)*r(a,X,m,n)",test_tensor,interm_t%w_oooo,r2_tensor,scale=ONE_HALF)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 6 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm6")
      norm_s_prev = norm_s


!-------------------------------------------------------------
!     term 4: S2(a,X,i,j) -= P(i,j) W(m,a,i,e) * R2(e,X,m,j) |
!-------------------------------------------------------------

      ! ierr=t_einsum("S(a,X,i,j)+=W(a,m,e,i)*R(e,X,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)!Used
      ! ierr=t_einsum("S(a,X,i,j)+=W(a,m,e,j)*R(e,X,m,i)",sigma2_tensor,interm_t%w_vovo,r2_tensor)


      ierr=t_einsum("S(a,X,i,j)+=W(a,m,e,i)*R(e,X,m,j)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)!XY !diam7
      ierr=t_einsum("S(a,X,i,j)+=W(a,m,e,j)*R(e,X,m,i)",sigma2_tensor,interm_t%w_vovo,r2_tensor)!XY !diam7


      ! Test
      ! ierr=t_einsum("S(a,X,i,j)+=W(a,m,e,i)*R(e,X,j,m)",test_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)
      ! ierr=t_einsum("S(a,X,i,j)+=W(a,m,e,j)*R(e,X,i,m)",test_tensor,interm_t%w_vovo,r2_tensor)


      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 7 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm7")
      ! if (debug) call print_tensor(interm_t%w_vovo, 1.0d-10, "w_ovov")
      norm_s_prev = norm_s

!--------------------------------------------------
!     term 5: S2(a,X,i,j) += F(a,e) * R2(e,X,i,j) |
!--------------------------------------------------

      ierr=t_einsum("S(a,X,i,j)+=F(a,e)*R(e,X,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor) !diam8

      ! Test
      ! ierr=t_einsum("S(a,X,i,j)+=F(a,e)*R(e,X,i,j)",test_tensor,interm_t%fbar_vv,r2_tensor)
      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 8 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm8")
      norm_s_prev = norm_s


!---------------------------------------------------------------------------
!     term 6: S2(a,X,i,j) += 1/2 ( V(m,n,e,f) * R2(f,X,m,n) ) * T(e,a,i,j) |
!---------------------------------------------------------------------------

      ierr=construct_tensor_t(v_aux,C8,(/1,interm_t%nvir/),init_val=ZERO)
      ierr=t_einsum("U(X,e)+=V(m,n,e,f)*R(f,X,m,n)",v_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ! ierr=t_einsum("S(a,X,i,j)+=U(X,e)*T(e,a,i,j)",sigma2_tensor,v_aux,interm_t%t2)

      ierr=t_einsum("S(a,X,i,j)+=U(X,e)*T(e,a,i,j)",sigma2_tensor,v_aux,interm_t%t2,scale=MINUS_ONE)  !XY !diam9

      !Test
      ! ierr=t_einsum("S(a,X,i,j)+=U(X,e)*T(e,a,i,j)",test_tensor,v_aux,interm_t%t2)

      ierr=destroy_tensor_t(v_aux)

      norm_s = tensor_norm2(sigma2_tensor)
      if (debug) print*, "norm contraction 9 (accumulated, this):",norm_s, (norm_s - norm_s_prev)
      ! if (debug) call print_tensor(test_tensor, 1.0d-10, "diagm9")
      norm_s_prev = norm_s


      ! debug
      if (debug) then
        print*, "norm r1_sigma", tensor_norm2(sigma1_tensor)
        print*, "norm r2_sigma", tensor_norm2(sigma2_tensor),tensor_norm2(sigma2_tensor)*ONE_QUARTER_R
        print*, "norm r_sigma ", (tensor_norm2(sigma1_tensor) + tensor_norm2(sigma2_tensor))
      end if

    end subroutine eti_sigma_eomccsd_rhs_ip


    subroutine eti_sigma_eomccsd_lhs_ip(int_t,interm_t,l1_tensor, l2_tensor, sigma1_tensor, sigma2_tensor)

!     input variables
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(tensor_t),    intent(inout) :: l1_tensor
      type(tensor_t),    intent(inout) :: l2_tensor
      type(tensor_t),    intent(inout) :: sigma1_tensor
      type(tensor_t),    intent(inout) :: sigma2_tensor
!     auxiliary tensor
      type(tensor_t) :: v_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=initialize_tensor_t(sigma1_tensor)

!-----------------------------------------
!     term 1: S1(i,X) -= F(i,m) * L(m,X) |
!-----------------------------------------

      ierr=t_einsum("S(i,X)+=F(i,m)*L(m,X)",sigma1_tensor,interm_t%fbar_oo,l1_tensor,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 2: S1(i,X) += 1/2 W(i,e,m,n) * L(m,n,e,X) |
!-----------------------------------------------------

      ierr=t_einsum("S(i,X)+=W(i,e,m,n)*L(m,n,e,X)",sigma1_tensor,interm_t%w_ovoo,l2_tensor,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=initialize_tensor_t(sigma2_tensor)

!---------------------------------------------
!     term 1: S2(i,j,a,X) -= W(i,j,m,a) * L(m,X) |
!---------------------------------------------

      ierr=t_einsum("S(i,j,a,X)+=W(j,i,m,a)*L(m,X)",sigma2_tensor,interm_t%w_ooov,l1_tensor)

!---------------------------------------------
!     term 2: S2(i,j,a,X) += F(e,a) * L(i,j,e,X) |
!---------------------------------------------

      ierr=t_einsum("S(i,j,a)+=F(e,a)*L(i,j,e)",sigma2_tensor,interm_t%fbar_vv,l2_tensor)

!----------------------------------------------------
!     term 3: S2(i,j,a,X) -= P(i,j) F(j,m) * L(i,m,a,X) |
!----------------------------------------------------

      ierr=t_einsum("S(i,j,a,X)+=F(j,m)*L(m,i,a,X)",sigma2_tensor,interm_t%fbar_oo,l2_tensor)
      ierr=t_einsum("S(i,j,a,X)+=F(i,m)*L(j,m,a,X)",sigma2_tensor,interm_t%fbar_oo,l2_tensor)

!--------------------------------------------------------
!     term 4: S2(i,j,a,X) += P(i,j) W(j,e,a,m) * L(i,m,e,X) |
!--------------------------------------------------------

      ierr=t_einsum("S(i,j,a,X)+=W(e,j,a,m)*L(m,i,e,X)",sigma2_tensor,interm_t%w_vovo,l2_tensor)
      ierr=t_einsum("S(i,j,a,X)+=W(e,i,a,m)*L(j,m,e,X)",sigma2_tensor,interm_t%w_vovo,l2_tensor)

!-----------------------------------------------------
!     term 5: S2(i,j,a,X) += 1/2 W(i,j,m,n) * L(m,n,a,X) |
!-----------------------------------------------------

      ierr=t_einsum("S(i,j,a,X)+=W(i,j,m,n)*L(m,n,a,X)",sigma2_tensor,interm_t%w_oooo,l2_tensor,scale=ONE_HALF)

!----------------------------------------------------------------------
!     term 6: S2(i,j,a,X) += 1/2 V(i,j,e,a) * ( L(m,n,f,X) * T(e,f,m,n) ) |
!----------------------------------------------------------------------

      ierr=construct_tensor_t(v_aux,C8,(/interm_t%nvir, 1/),init_val=ZERO)
      ierr=t_einsum("U(e,X)+=L(m,n,f,X)*T(e,f,m,n)",v_aux,l2_tensor,interm_t%t2,scale=ONE_HALF)
      ierr=t_einsum("S(i,j,a,X)+=V(i,j,e,a)*U(e,X)",sigma2_tensor,int_t%oovv,v_aux)
      ierr=destroy_tensor_t(v_aux)

!------------------------------------------------
!     term 7: S2(i,j,a,X) += P(i,j) L(i,X) * F(j,a) |
!------------------------------------------------

      ierr=t_einsum("S(i,j,a,X)+=L(i,X)*F(j,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
      ierr=t_einsum("S(i,j,a,X)+=L(j,X)*F(i,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)

    end subroutine eti_sigma_eomccsd_lhs_ip

    subroutine eti_sigma_eomccsd_rhs_ea(int_t, interm_t, r1_tensor, r2_tensor, sigma1_tensor, sigma2_tensor)
!     input variables
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(tensor_t),    intent(inout) :: sigma1_tensor
      type(tensor_t),    intent(inout) :: sigma2_tensor
      type(tensor_t),    intent(inout) :: r1_tensor
      type(tensor_t),    intent(inout) :: r2_tensor
!     auxiliary tensor
      type(tensor_t) :: o_aux
!     error code
      integer :: ierr

!***************
!     r1_sigma *
!***************

      ierr=initialize_tensor_t(sigma1_tensor)

!--------------------------------------
!     term 1: S1(a,X) += F(a,e) * R1(e,X) |
!--------------------------------------

      ! ierr=t_einsum("S(a,X)+=F(e,a)*R(e,X)",sigma1_tensor,interm_t%fbar_vv,r1_tensor)
      ierr=t_einsum("S(a,X)+=F(a,e)*R(e,X)",sigma1_tensor,interm_t%fbar_vv,r1_tensor)  !XY
      ! print*, "r1_sigma", tensor_norm2(sigma1_tensor)

!------------------------------------------
!     term 2: S1(a,X) += F(m,e) * R2(e,a,m,X) |
!------------------------------------------

      ierr=t_einsum("S(a,X)+=F(m,e)*R(e,a,m,X)",sigma1_tensor,interm_t%fbar_ov,r2_tensor)
      ! print*, "r1_sigma", tensor_norm2(sigma1_tensor)

!--------------------------------------------------
!     term 3: S1(a,X) += 1/2 W(a,m,f,e) * R2(e,f,m,X) |
!--------------------------------------------------

      ierr=t_einsum("S(a,X)+=W(a,m,f,e)*R(e,f,m,X)",sigma1_tensor,interm_t%w_vovv,r2_tensor,scale=ONE_HALF)
      ! print*, "r1_sigma", tensor_norm2(sigma1_tensor)

!***************
!     r2_sigma *
!***************

      ierr=initialize_tensor_t(sigma2_tensor)

!---------------------------------------------
!     term 1: S2(a,b,i,X) -= W(a,b,e,i) * R(e,X) |
!---------------------------------------------

      ! ierr=t_einsum("S(a,b,i,X)+=W(a,b,e,i)*R(e,X)",sigma2_tensor,interm_t%w_vvvo,r1_tensor)
      ierr=t_einsum("S(a,b,i,X)+=W(a,b,e,i)*R(e,X)",sigma2_tensor,interm_t%w_vvvo,r1_tensor,&
                                    scale=MINUS_ONE) !XY
      ! print*, "r2_sigma", tensor_norm2(sigma2_tensor)*ONE_HALF


!-----------------------------------------------------
!     term 2: S2(a,b,i,X) += P(a,b) F(a,e) * R2(e,b,i,X) |
!-----------------------------------------------------

      ierr=t_einsum("S(a,b,i,X)+=F(a,e)*r(e,b,i,X)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)
      ierr=t_einsum("S(a,b,i,X)+=F(b,e)*r(a,e,i,X)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)

!--------------------------------------------------
!     term 3: S2(a,b,i,X) += W(a,b,e,f) * R2(e,f,i,X) |
!--------------------------------------------------

      ierr=t_einsum("S(a,b,i,X)+=W(a,b,e,f)*r(e,f,i,X)",sigma2_tensor,interm_t%w_vvvv,r2_tensor,scale=ONE_HALF)

!---------------------------------------------------------
!     term 4: S2(a,b,i,X) -= P(a,b) W(a,m,e,i) * R2(e,b,m,X) |
!---------------------------------------------------------

      ierr=t_einsum("S(a,b,i,X)+=W(a,m,e,i)*R(b,e,m,X)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ierr=t_einsum("S(a,b,i,X)+=W(b,m,e,i)*R(e,a,m,X)",sigma2_tensor,interm_t%w_vovo,r2_tensor)

!----------------------------------------------
!     term 5: S2(a,b,i,X) -= F(m,i) * R2(a,b,m,X) |
!----------------------------------------------

      ierr=t_einsum("S(a,b,i,X)+=F(m,i)*R(b,a,m,X)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)

!-----------------------------------------------------------------------
!     term 6: S2(a,b,i,X) -= 1/2 ( V(m,n,e,f) * R2(e,f,n,X) ) * T(a,b,m,i) |
!-----------------------------------------------------------------------

      ierr=construct_tensor_t(o_aux,C8,(/interm_t%nocc, 1/),init_val=ZERO)
      ierr=t_einsum("U(m,X)+=V(m,n,e,f)*R(e,f,n,X)",o_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ! ierr=t_einsum("S(a,b,i,X)+=U(m,X)*T(a,b,m,i)",sigma2_tensor,o_aux,interm_t%t2)
      ierr=t_einsum("S(a,b,i,X)+=U(m,X)*T(a,b,m,i)",sigma2_tensor,o_aux,interm_t%t2,scale=MINUS_ONE) !XY
      ierr=destroy_tensor_t(o_aux)

      ! debug
      ! if (.true.) then
      !   print*, "r1_sigma", tensor_norm2(sigma1_tensor)
      !   print*, "r2_sigma", tensor_norm2(sigma2_tensor)*ONE_HALF
      ! end if

    end subroutine eti_sigma_eomccsd_rhs_ea

    subroutine eti_sigma_eomccsd_lhs_ea(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)

!     input variables
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(tensor_t),    intent(inout) :: sigma1_tensor
      type(tensor_t),    intent(inout) :: sigma2_tensor
      type(tensor_t),    intent(inout) :: l1_tensor
      type(tensor_t),    intent(inout) :: l2_tensor
!     auxiliary tensor
      type(tensor_t) :: o_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=initialize_tensor_t(sigma1_tensor)

!-------------------------------------
!     term 1: S1(X,a) += F(e,a) * L(X,e) |
!-------------------------------------

      ierr=t_einsum("S(X,a)+=F(e,a)*L(X,e)",sigma1_tensor,interm_t%fbar_vv,l1_tensor)

!-------------------------------------------------
!     term 2: S1(X,a) -= 1/2 W(e,f,a,m) * L(m,X,f,e) |
!-------------------------------------------------

      ierr=t_einsum("S(X,a)+=W(e,f,a,m)*L(m,X,e,f)",sigma1_tensor,interm_t%w_vvvo,l2_tensor,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=initialize_tensor_t(sigma2_tensor)

!---------------------------------------------
!     term 1: S2(i,X,a,b) += W(e,i,b,a) * L(X,e) |
!---------------------------------------------

      ierr=t_einsum("S(i,X,a,b)+=W(e,i,b,a)*L(X,e)",sigma2_tensor,interm_t%w_vovv,l1_tensor)

!----------------------------------------------------
!     term 2: S2(i,X,a,b) += P(a,b) F(e,a) * L(i,X,e,b) |
!----------------------------------------------------

      ierr=t_einsum("S(i,X,a,b)+=F(e,a)*L(i,X,e,b)",sigma2_tensor,interm_t%fbar_vv,l2_tensor)
      ierr=t_einsum("S(i,X,a,b)+=F(e,b)*L(i,X,a,e)",sigma2_tensor,interm_t%fbar_vv,l2_tensor)

!---------------------------------------------
!     term 3: S2(i,X,a,b) -= F(i,m) * L(m,X,a,b) |
!---------------------------------------------

      ierr=t_einsum("S(i,X,a,b)+=F(i,m)*L(m,X,b,a)",sigma2_tensor,interm_t%fbar_oo,l2_tensor)

!-----------------------------------------------------
!     term 4: S2(i,X,a,b) += 1/2 W(e,f,a,b) * L(i,X,e,f) |
!-----------------------------------------------------

      ierr=t_einsum("S(i,X,a,b)+=W(e,f,a,b)*L(i,X,e,f)",sigma2_tensor,interm_t%w_vvvv,l2_tensor,scale=ONE_HALF)

!--------------------------------------------------------
!     term 5: S2(i,X,a,b) -= P(a,b) W(e,i,b,m) * L(X,m,a,e) |
!--------------------------------------------------------

      ierr=t_einsum("S(i,X,a,b)+=W(e,i,b,m)*L(X,m,e,a)",sigma2_tensor,interm_t%w_vovo,l2_tensor)
      ierr=t_einsum("S(i,X,a,b)+=W(e,i,a,m)*L(X,m,b,e)",sigma2_tensor,interm_t%w_vovo,l2_tensor)

!----------------------------------------------------------------------
!     term 6: S2(i,X,a,b) -= 1/2 ( L(X,o,e,f) * T(e,f,m,o) ) * V(m,i,b,a) |
!----------------------------------------------------------------------

      ierr=construct_tensor_t(o_aux,C8,(/1, interm_t%nocc/),init_val=ZERO)
      ierr=t_einsum("U(X,m)+=L(X,o,e,f)*T(e,f,o,m)",o_aux,l2_tensor,interm_t%t2,scale=ONE_HALF)
      ierr=t_einsum("S(i,X,a,b)+=U(X,m)*V(m,i,b,a)",sigma2_tensor,o_aux,int_t%oovv)
      ierr=destroy_tensor_t(o_aux)

!------------------------------------------------
!     term 7: S2(i,X,a,b) += P(a,b) F(i,a) * L(X,b) |
!------------------------------------------------

      ierr=t_einsum("S(i,X,a,b)+=F(i,a)*L(X,b)",sigma2_tensor,interm_t%fbar_ov,l1_tensor)
      ierr=t_einsum("S(i,X,a,b)+=F(i,b)*L(X,a)",sigma2_tensor,interm_t%fbar_ov,l1_tensor,scale=MINUS_ONE)

    end subroutine eti_sigma_eomccsd_lhs_ea

    subroutine eti_sigma_eomccsd_rhs_ee(int_t,interm_t,r1_tensor,r2_tensor,sigma1_tensor,sigma2_tensor)

!     input variables
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(tensor_t),    intent(inout) :: sigma1_tensor
      type(tensor_t),    intent(inout) :: sigma2_tensor
      type(tensor_t),    intent(inout) :: r1_tensor
      type(tensor_t),    intent(inout) :: r2_tensor
!     auxiliary tensor
      type(tensor_t) :: oo_aux, vv_aux
      integer            :: oo_dims(2), vv_dims(2)
!     error code
      integer :: ierr

      ! type(tensor_t) :: test_tensor

      oo_dims = (/interm_t%nocc,interm_t%nocc/)
      vv_dims = (/interm_t%nvir,interm_t%nvir/)

      ! ierr=construct_tensor_t(test_tensor,C8,(/interm_t%nvir,interm_t%nvir,interm_t%nocc,interm_t%nocc/),init_val=ZERO)

!***************
!     r1_sigma *
!***************

      ierr=initialize_tensor_t(sigma1_tensor)

!-----------------------------------------
!     term 1: S1(a,i) += F(a,e) * R(e,i) |
!-----------------------------------------

      ierr=t_einsum("S(a,i)+=F(a,e)*R(e,i)",sigma1_tensor,interm_t%fbar_vv,r1_tensor)

!-----------------------------------------
!     term 2: S1(a,i) -= F(m,i) * R(a,m) |
!-----------------------------------------

      ierr=t_einsum("S(a,i)+=F(m,i)*R(a,m)",sigma1_tensor,interm_t%fbar_oo,r1_tensor,scale=MINUS_ONE)

!---------------------------------------------
!     term 3: S1(a,i) += F(m,e) * R(e,a,m,i) |
!---------------------------------------------

      ierr=t_einsum("S(a,i)+=F(m,e)*R(e,a,m,i)",sigma1_tensor,interm_t%fbar_ov,r2_tensor)

!---------------------------------------------
!     term 4: S1(a,i) += W(m,a,e,i) * R(e,m) |
!---------------------------------------------
      ! Flip the index since we only have w_vovo
      ierr=t_einsum("S(a,i)+=W(a,m,e,i)*R(e,m)",sigma1_tensor,interm_t%w_vovo,r1_tensor,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 5: S1(a,i) += 1/2 W(a,m,e,f) * R(e,f,i,m) |
!-----------------------------------------------------

      ierr=t_einsum("S(a,i)+=W(a,m,e,f)*R(e,f,i,m)",sigma1_tensor,interm_t%w_vovv,r2_tensor,scale=ONE_HALF)
      
!-----------------------------------------------------
!     term 6: S1(a,i) += 1/2 W(m,n,i,e) * R(e,a,m,n) |
!-----------------------------------------------------

      ! ierr=t_einsum("S(a,i)+=W(m,n,i,e)*R(a,e,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=ONE_HALF)

      ierr=t_einsum("S(a,i)+=W(m,n,i,e)*R(e,a,m,n)",sigma1_tensor,interm_t%w_ooov,r2_tensor,scale=ONE_HALF) !XY

!***************
!     r2_sigma *
!***************

      ierr=initialize_tensor_t(sigma2_tensor)

!--------------------------------------------------------
!     term 1: S2(a,b,i,j) -= P(a,b) W(m,b,i,j) * R(a,m) |
!--------------------------------------------------------

      ! ierr=t_einsum("S(a,b,i,j)+=W(m,b,j,i)*R(a,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor)
      ! ierr=t_einsum("S(a,b,i,j)+=W(m,a,i,j)*R(b,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor)

      ierr=t_einsum("S(a,b,i,j)+=W(m,b,i,j)*R(a,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor,scale=MINUS_ONE)     !XY
      ierr=t_einsum("S(a,b,i,j)+=W(m,a,i,j)*R(b,m)",sigma2_tensor,interm_t%w_ovoo,r1_tensor)     !XY

      ! call print_tensor(sigma2_tensor, 1.0d-10, "diagm7")


!--------------------------------------------------------
!     term 2: S2(a,b,i,j) += P(i,j) W(a,b,e,j) * R(e,i) |
!--------------------------------------------------------

      ierr=t_einsum("S(a,b,i,j)+=W(a,b,e,j)*R(e,i)",sigma2_tensor,interm_t%w_vvvo,r1_tensor)  !Used
      ! ierr=t_einsum("S(a,b,i,j)+=W(b,a,e,i)*R(e,j)",sigma2_tensor,interm_t%w_vvvo,r1_tensor)

      ierr=t_einsum("S(a,b,i,j)+=W(a,b,e,i)*R(e,j)",sigma2_tensor,interm_t%w_vvvo,r1_tensor,scale=MINUS_ONE)   !XY
      ! call print_tensor(interm_t%w_vvvo, 1.0d-10, "w_vvvo")
      ! call print_tensor(sigma2_tensor, 1.0d-10, "diagm8")

!-------------------------------------------------------------------------
!     term 3: S2(a,b,i,j) += P(a,b) ( W(b,m,f,e) * R(e,m) ) * T(a,f,i,j) |
!-------------------------------------------------------------------------

      ierr=construct_tensor_t(vv_aux,C8,vv_dims,init_val=ZERO)
      ierr=t_einsum("U(b,f)+=W(b,m,f,e)*R(e,m)",vv_aux,interm_t%w_vovv,r1_tensor)
      ierr=t_einsum("S(a,b,i,j)+=U(b,f)*T2(a,f,i,j)",sigma2_tensor,vv_aux,interm_t%t2)     !Used
      ! ierr=t_einsum("S(a,b,i,j)+=U(a,f)*T2(f,b,i,j)",sigma2_tensor,vv_aux,interm_t%t2)

      ierr=t_einsum("S(a,b,i,j)+=U(a,f)*T2(b,f,i,j)",sigma2_tensor,vv_aux,interm_t%t2,scale=MINUS_ONE)      !XY
      ierr=destroy_tensor_t(vv_aux)
      ! call print_tensor(interm_t%w_vovv, 1.0d-10, "w_vovv")
      ! call print_tensor(sigma2_tensor, 1.0d-10, "diagm9")


!-------------------------------------------------------------------------
!     term 4: S2(a,b,i,j) -= P(i,j) ( W(n,m,j,e) * R(e,m) ) * T(a,b,i,n) |
!-------------------------------------------------------------------------

      ierr=construct_tensor_t(oo_aux,C8,oo_dims,init_val=ZERO)
      ierr=t_einsum("U(n,j)+=W(n,m,j,e)*R(e,m)",oo_aux,interm_t%w_ooov,r1_tensor)
      ! ierr=t_einsum("S(a,b,i,j)+=U(n,j)*T(a,b,i,n)",sigma2_tensor,oo_aux,interm_t%t2)
      ! ierr=t_einsum("S(a,b,i,j)+=U(n,i)*T(a,b,n,j)",sigma2_tensor,oo_aux,interm_t%t2)


      ierr=t_einsum("S(a,b,i,j)+=U(n,j)*T(a,b,i,n)",sigma2_tensor,oo_aux,interm_t%t2,scale=MINUS_ONE)        !XY
      ierr=t_einsum("S(a,b,i,j)+=U(n,i)*T(a,b,j,n)",sigma2_tensor,oo_aux,interm_t%t2)        !XY
      ierr=destroy_tensor_t(oo_aux)

!--------------------------------------------------------
!     term 5: S2(a,b,i,j) += P(a,b) F(b,e) * R(a,e,i,j) |
!--------------------------------------------------------

      ierr=t_einsum("S(a,b,i,j)+=F(b,e)*R(a,e,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)   !Used
      ! ierr=t_einsum("S(a,b,i,j)+=F(a,e)*R(e,b,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor)

      ierr=t_einsum("S(a,b,i,j)+=F(a,e)*R(b,e,i,j)",sigma2_tensor,interm_t%fbar_vv,r2_tensor,scale=MINUS_ONE)   !XY


!--------------------------------------------------------
!     term 6: S2(a,b,i,j) -= P(i,j) F(m,j) * R(a,b,i,m) |
!--------------------------------------------------------

      ! ierr=t_einsum("S(a,b,i,j)+=F(m,j)*R(a,b,m,i)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)
      ! ierr=t_einsum("S(a,b,i,j)+=F(m,i)*R(a,b,j,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)


      ierr=t_einsum("S(a,b,i,j)+=F(m,j)*R(a,b,i,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor,scale=MINUS_ONE)   !XY
      ierr=t_einsum("S(a,b,i,j)+=F(m,i)*R(a,b,j,m)",sigma2_tensor,interm_t%fbar_oo,r2_tensor)            !XY
!---------------------------------------------------------
!     term 7: S2(a,b,i,j) += 1/2 W(m,n,i,j) * R(a,b,m,n) |
!---------------------------------------------------------

      ierr=t_einsum("S(a,b,i,j)+=W(m,n,i,j)*R(a,b,m,n)",sigma2_tensor,interm_t%w_oooo,r2_tensor,scale=ONE_HALF)

!-------------------------------------------------------------------
!     term 8: S2(a,b,i,j) += P(i,j) P(a,b) W(m,b,e,j) * R(a,e,i,m) |
!-------------------------------------------------------------------


      ! ierr=t_einsum("S(a,b,i,j)+=W(b,m,e,j)*R(a,e,m,i)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=t_einsum("S(a,b,i,j)+=W(a,m,e,j)*R(b,e,i,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=t_einsum("S(a,b,i,j)+=W(b,m,e,i)*R(a,e,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ! ierr=t_einsum("S(a,b,i,j)+=W(a,m,e,i)*R(b,e,m,j)",sigma2_tensor,interm_t%w_vovo,r2_tensor)


      ! Filp the index since we only have w_vovo, XY
      ierr=t_einsum("S(a,b,i,j)+=W(b,m,e,j)*R(a,e,i,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)
      ierr=t_einsum("S(a,b,i,j)+=W(a,m,e,j)*R(b,e,i,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ierr=t_einsum("S(a,b,i,j)+=W(b,m,e,i)*R(a,e,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor)
      ierr=t_einsum("S(a,b,i,j)+=W(a,m,e,i)*R(b,e,j,m)",sigma2_tensor,interm_t%w_vovo,r2_tensor,scale=MINUS_ONE)

!---------------------------------------------------------------------------------
!     term 9: S2(a,b,i,j) -= 1/2 P(a,b) ( V(n,m,f,e) * R(e,a,m,n) ) * T(f,b,i,j) |
!---------------------------------------------------------------------------------

      ierr=construct_tensor_t(vv_aux,C8,vv_dims,init_val=ZERO)
      ierr=t_einsum("U(a,f)+=V(n,m,f,e)*R(e,a,m,n)",vv_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)
      ! ierr=t_einsum("S(a,b,i,j)+=U(a,f)*T(f,b,j,i)",sigma2_tensor,vv_aux,interm_t%t2)
      ! ierr=t_einsum("S(a,b,i,j)+=U(b,f)*T(f,a,i,j)",sigma2_tensor,vv_aux,interm_t%t2)


      ierr=t_einsum("S(a,b,i,j)+=U(a,f)*T(f,b,i,j)",sigma2_tensor,vv_aux,interm_t%t2,scale=MINUS_ONE)   !XY
      ierr=t_einsum("S(a,b,i,j)+=U(b,f)*T(f,a,i,j)",sigma2_tensor,vv_aux,interm_t%t2)                   !XY
      ierr=destroy_tensor_t(vv_aux)

!---------------------------------------------------------------------------------
!     term 10: S2(a,b,i,j) -= 1/2 P(i,j) ( V(n,m,f,e) * R(f,e,i,m) ) * T(b,a,j,n) |
!---------------------------------------------------------------------------------

      ierr=construct_tensor_t(oo_aux,C8,oo_dims,init_val=ZERO)
      ierr=t_einsum("U(n,i)+=V(n,m,f,e)*R(f,e,i,m)",oo_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)      !Used 
      ! ierr=t_einsum("U(n,j)+=V(n,m,f,e)*R(f,e,j,m)",oo_aux,int_t%oovv,r2_tensor,scale=ONE_HALF)    !XY
      ierr=t_einsum("S(a,b,i,j)+=U(n,i)*T(b,a,n,j)",sigma2_tensor,oo_aux,interm_t%t2)                !Used 
      ierr=t_einsum("S(a,b,i,j)+=U(n,j)*T(b,a,i,n)",sigma2_tensor,oo_aux,interm_t%t2)                !Used


      ! ierr=t_einsum("S(a,b,i,j)+=U(n,j)*T(b,a,i,n)",sigma2_tensor,oo_aux,interm_t%t2,scale=MINUS_ONE)  !XY
      ! ierr=t_einsum("S(a,b,i,j)+=U(n,i)*T(b,a,j,n)",sigma2_tensor,oo_aux,interm_t%t2)   !XY
      ierr=destroy_tensor_t(oo_aux)

!----------------------------------------------------------
!     term 11: S2(a,b,i,j) += 1/2 W(a,b,e,f) * R(e,f,i,j) |
!----------------------------------------------------------

      ierr=t_einsum("S(a,b,i,j)+=W(a,b,e,f)*R(e,f,i,j)",sigma2_tensor,interm_t%w_vvvv,r2_tensor,scale=ONE_HALF)

    end subroutine eti_sigma_eomccsd_rhs_ee

    subroutine eti_sigma_eomccsd_lhs_ee(int_t,interm_t,l1_tensor,l2_tensor,sigma1_tensor,sigma2_tensor)

!     input variables
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(tensor_t),    intent(inout) :: sigma1_tensor
      type(tensor_t),    intent(inout) :: sigma2_tensor
      type(tensor_t),    intent(inout) :: l1_tensor
      type(tensor_t),    intent(inout) :: l2_tensor
      type(tensor_t) :: one_tensor
!     auxiliary tensor
      type(tensor_t) :: goo_tensor, gvv_tensor, oo_aux
      integer            :: oo_dims(2), vv_dims(2)
      integer(C_INT)     :: one_dims(1)
!     error code
      integer :: ierr

      one_dims(1) = 1
      oo_dims = (/interm_t%nocc,interm_t%nocc/)
      vv_dims = (/interm_t%nvir,interm_t%nvir/)
      ierr=construct_tensor_t(goo_tensor,C8,oo_dims,init_val=ZERO)
      ierr=construct_tensor_t(gvv_tensor,C8,vv_dims,init_val=ZERO)
      ierr=t_einsum("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)
      ierr=t_einsum("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,interm_t%t2,scale=ONE_HALF)

!***************
!     l1_sigma *
!***************

      ierr=initialize_tensor_t(sigma1_tensor)

!--------------------------------
!     term 1: S1(i,a) += F(i,a) |
!--------------------------------

      ! ierr=t_einsum("S(i,a)+=F(i,a)",sigma1_tensor,interm_t%fbar_ov)
      ! call print_tensor(sigma1_tensor, 1.0d-10, "sigma1")
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma1_tensor)
    
      
!-----------------------------------------
!     term 2: S1(i,a) += L(i,e) * F(e,a) |
!-----------------------------------------

      ierr=t_einsum("S(i,a)+=L(i,e)*F(e,a)",sigma1_tensor,l1_tensor,interm_t%fbar_vv)
      ! call print_tensor(interm_t%fbar_vv, 1.0d-10, "fbar_vv")
      ! call print_tensor(sigma1_tensor, 1.0d-10, "sigma1")
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma1_tensor)

!-----------------------------------------
!     term 3: S1(i,a) -= L(m,a) * F(i,m) |
!-----------------------------------------

      ierr=t_einsum("S(i,a)+=L(m,a)*F(i,m)",sigma1_tensor,l1_tensor,interm_t%fbar_oo,scale=MINUS_ONE)
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma1_tensor)

!-----------------------------------------------------
!     term 4: S1(i,a) += 1/2 L(i,m,e,f) * W(e,f,a,m) |
!-----------------------------------------------------

      ierr=t_einsum("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",sigma1_tensor,l2_tensor,interm_t%w_vvvo,scale=ONE_HALF)

!---------------------------------------------
!     term 5: S1(i,a) -= G(f,e) * W(e,i,f,a) |
!---------------------------------------------

      ! ierr=t_einsum("S(i,a)+=G(f,e)*W(e,i,a,f)",sigma1_tensor,interm_t%gvv,interm_t%w_vovv)
      ierr=t_einsum("S(i,a)+=G(f,e)*W(e,i,a,f)",sigma1_tensor,gvv_tensor,interm_t%w_vovv)

!---------------------------------------------
!     term 6: S1(i,a) -= G(n,m) * W(m,i,n,a) |
!---------------------------------------------

      ! ierr=t_einsum("S(i,a)+=G(n,m)*W(i,m,n,a)",sigma1_tensor,interm_t%goo,interm_t%w_ooov)
      ierr=t_einsum("S(i,a)+=G(n,m)*W(i,m,n,a)",sigma1_tensor,goo_tensor,interm_t%w_ooov)
      ! call print_tensor(sigma1_tensor, 1.0d-10, "sigma1")
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma1_tensor)

      

!---------------------------------------------
!     term 7: S1(i,a) += L(m,e) * W(i,e,a,m) |
!---------------------------------------------

      ierr=t_einsum("S(i,a)+=L(m,e)*W(e,i,a,m)",sigma1_tensor,l1_tensor,interm_t%w_vovo,scale=MINUS_ONE)
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma1_tensor)


!-----------------------------------------------------
!     term 8: S1(i,a) -= 1/2 L(m,n,a,e) * W(i,e,m,n) |
!-----------------------------------------------------

      ierr=t_einsum("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",sigma1_tensor,l2_tensor,interm_t%w_ovoo,scale=ONE_HALF)
      ! print*, "norm contraction sigma1 (accumulated, this):", tensor_norm2(sigma1_tensor)

!***************
!     l2_sigma *
!***************

      ierr=initialize_tensor_t(sigma2_tensor)

!----------------------------------------
!     term 1: S2(i,j,a,b) += V(i,j,a,b) |
!----------------------------------------

      ! ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,b)",sigma2_tensor,int_t%oovv)

!--------------------------------------------------------
!     term 2: S2(i,j,a,b) += P(a,b) L(i,j,a,e) * F(e,b) |
!--------------------------------------------------------

      ierr=t_einsum("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
      ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",sigma2_tensor,l2_tensor,interm_t%fbar_vv)
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0


!--------------------------------------------------------
!     term 3: S2(i,j,a,b) -= P(i,j) L(i,m,a,b) * F(j,m) |
!--------------------------------------------------------

      ierr=t_einsum("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
      ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",sigma2_tensor,l2_tensor,interm_t%fbar_oo)
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0

!---------------------------------------------------------
!     term 4: S2(i,j,a,b) += 1/2 L(m,n,a,b) * W(i,j,m,n) |
!---------------------------------------------------------

      ierr=t_einsum("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",sigma2_tensor,l2_tensor,interm_t%w_oooo,scale=ONE_HALF)
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0


!-------------------------------------------------------------------
!     term 5: S2(i,j,a,b) += P(i,j) P(a,b) L(i,m,a,e) * W(j,e,b,m) |
!-------------------------------------------------------------------

      ierr=t_einsum("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ierr=t_einsum("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ierr=t_einsum("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",sigma2_tensor,l2_tensor,interm_t%w_vovo)
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0

!--------------------------------------------------------
!     term 6: S2(i,j,a,b) += P(a,b) V(i,j,a,e) * G(e,b) |
!--------------------------------------------------------

      ! ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",sigma2_tensor,int_t%oovv,interm_t%gvv)
      ! ierr=t_einsum("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",sigma2_tensor,int_t%oovv,interm_t%gvv)

      ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",sigma2_tensor,int_t%oovv,gvv_tensor)   !XY
      ierr=t_einsum("S(i,j,a,b)+=V(i,j,b,e)*G(e,a)",sigma2_tensor,int_t%oovv,gvv_tensor,scale=MINUS_ONE)   !XY
      ! print*, "norm contraction 1 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0


!-------------------------------------------------
!     term 7: S2(i,j,a,b) -= P(a,b) L(m,a) * W(i,j,m,b) |
!-------------------------------------------------

      ! ierr=t_einsum("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",sigma2_tensor,l1_tensor,interm_t%w_ooov)
      ierr=t_einsum("S(i,j,a,b)+=L(m,a)*W(i,j,m,b)",sigma2_tensor,l1_tensor,interm_t%w_ooov,scale=MINUS_ONE)  !XY
      ierr=t_einsum("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",sigma2_tensor,l1_tensor,interm_t%w_ooov)   !XY

!--------------------------------------------------------
!     term 8: S2(i,j,a,b) -= P(i,j) V(i,m,a,b) * G(j,m) |
!--------------------------------------------------------

      ! ierr=t_einsum("S(i,j,a,b)+=V(m,i,a,b)*G(j,m)",sigma2_tensor,int_t%oovv,interm_t%goo)
      ! ierr=t_einsum("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",sigma2_tensor,int_t%oovv,interm_t%goo)

      ierr=t_einsum("S(i,j,a,b)+=V(m,i,a,b)*G(j,m)",sigma2_tensor,int_t%oovv,goo_tensor)   !XY
      ierr=t_einsum("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",sigma2_tensor,int_t%oovv,goo_tensor)   !XY

!---------------------------------------------------------------------
!     term 9: S2(i,j,a,b) -= P(i,j) V(m,j,a,b) * ( L(i,e) * T(e,m) ) |
!---------------------------------------------------------------------

      ! if (.not.interm_t%ccd) then
      ! !   oo_dims = (/interm_t%nocc,interm_t%nocc/)
      !   ierr=construct_tensor_t(oo_aux,C8,oo_dims,init_val=ZERO)
      ! !   ierr=t_einsum("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1)
      !   ierr=t_einsum("U(i,m)+=L(i,e)*T(e,m)",oo_aux,l1_tensor,interm_t%t1,scale=MINUS_ONE)  !XY
      !   ierr=t_einsum("S(i,j,a,b)+=V(m,j,a,b)*U(i,m)",sigma2_tensor,int_t%oovv,oo_aux)
      !   ierr=t_einsum("S(i,j,a,b)+=V(i,m,a,b)*U(j,m)",sigma2_tensor,int_t%oovv,oo_aux)
      !   ierr=destroy_tensor_t(oo_aux)
      ! end if

!---------------------------------------------------------
!     term 10: S2(i,j,a,b) += P(i,j) L(i,e) * V(e,j,a,b) |
!---------------------------------------------------------

      ! ierr=t_einsum("S(i,j,a,b)+=L(i,e)*V(e,j,a,b)",sigma2_tensor,l1_tensor,int_t%vovv)
      ! ierr=t_einsum("S(i,j,a,b)+=L(j,e)*V(e,i,b,a)",sigma2_tensor,l1_tensor,int_t%vovv)

      ierr=t_einsum("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",sigma2_tensor,l1_tensor,interm_t%w_vovv)   !XY
      ierr=t_einsum("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",sigma2_tensor,l1_tensor,interm_t%w_vovv)   !XY


!------------------------------------------------------------
!     term 11: S2(i,j,a,b) += P(i,j) P(a,b) L(i,a) * F(j,b) |
!------------------------------------------------------------

      ierr=t_einsum("S(i,j,a,b)+=L(i,a)*F(j,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)
      ierr=t_einsum("S(i,j,a,b)+=L(j,a)*F(i,b)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
      ierr=t_einsum("S(i,j,a,b)+=L(i,b)*F(j,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov,scale=MINUS_ONE)
      ierr=t_einsum("S(i,j,a,b)+=L(j,b)*F(i,a)",sigma2_tensor,l1_tensor,interm_t%fbar_ov)

!--------------------------------------------------------
!     term 12: S2(i,j,a,b) += 1/2 L(i,j,e,f) W(e,f,a,b) |
!--------------------------------------------------------

      ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",sigma2_tensor,l2_tensor,interm_t%w_vvvv,scale=ONE_HALF)
      ! print *, "new sigma vector:" 
      ! call print_tensor(L1_tensor, 1.0d-10, "L1")
      ! call print_tensor(l2_tensor, 1.0d-10, "L2")
      ! call print_tensor(sigma1_tensor, 1.0d-10, "sigma1")
      ! call print_tensor(sigma2_tensor, 1.0d-10, "sigma2")
      ! print*, "norm contraction sigma2 (accumulated, this):", tensor_norm2(sigma2_tensor)/4.0

      ierr=destroy_tensor_t(goo_tensor)
      ierr=destroy_tensor_t(gvv_tensor)


    end subroutine eti_sigma_eomccsd_lhs_ee


end module eti_sigma_eomccsd
