module talsh_eom

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates
  use talsh_ao_to_mo, only : tensor_norm2
  use exacorr_global, only : debug_exa

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0),ONE=(1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)
  real(8),    parameter :: ONE_QUARTER=0.25D0

  private

  public talsh_eom_driver

  contains

    subroutine talsh_eom_driver (exa_input,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)

!     input variables
      type(exacc_input),     intent(in   )           :: exa_input
      type(talsh_tens_t),    intent(inout)           :: t2_tensor
      type(talsh_intg_tens), intent(inout)           :: int_t
      type(talsh_tens_t),    intent(inout)           :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor, l2_tensor
!     eom tensors
      type(talsh_eom_tens) :: eom_t

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      if (.not.exa_input%ccd) then
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EA")
            ! call eom_ea(exa_input,eom_t,t2_tensor,int_t,t1_tensor)
          case("EE")
            ! call eom_ee(exa_input,eom_t,t2_tensor,int_t,l2_tensor,t1_tensor)
          case default
            call quit ('Error: No valid EOM-flavour chosen')
        end select
      else
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EA")
            ! call eom_ea(exa_input,eom_t,t2_tensor,int_t)
          case("EE")
            ! call eom_ee(exa_input,eom_t,t2_tensor,int_t,l2_tensor)
        case default
          call quit ('Error: No valid EOM-flavour chosen')
        end select
      end if

!     make some noise so that we know we are leaving
      print*, " Leaving eom_driver routine"
      return

    end subroutine talsh_eom_driver

    subroutine eom_ip(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!     tensor dims
      integer :: nocc, nvir
      integer, dimension(2) :: oo_dims, ov_dims, vv_dims
      integer, dimension(3) :: voo_dims, oov_dims
      integer, dimension(4) :: oooo_dims, ooov_dims, ovoo_dims, vovo_dims
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr

      CCD = exa_input%ccd

!************************
!     construct tensors *
!************************

      nocc = exa_input%nocc
      nvir = exa_input%nvir

!     set tensor dims
      oo_dims = (/nocc,nocc/)
      ov_dims = (/nocc,nvir/)
      vv_dims = (/nvir,nvir/)
      voo_dims = (/nvir,nocc,nocc/)
      oooo_dims = (/nocc,nocc,nocc,nocc/)
      ooov_dims = (/nocc,nocc,nocc,nvir/)
      ovoo_dims = (/nocc,nvir,nocc,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      if (exa_input%eom_properties) oov_dims = (/nocc,nocc,nvir/)

!     create all fixed intermediates
      ierr=talsh_tensor_construct(eom_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ooov,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

      call get_wbar_vovo(eom_t%wbar_vovo,t2_tensor,int_t)
      if (.not.CCD) then
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv,t1_tensor)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor,t1_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t,t1_tensor)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor,t1_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor)
      end if

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,(/nocc/),init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,(/nocc/),init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,(/nvir/),init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,(/nvir/),init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%r2_coeff,C8,voo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r2_sigma,C8,voo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l2_coeff,C8,oov_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l2_sigma,C8,oov_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      call sigma_right_ip(exa_input,int_t,eom_t,t2_tensor)

!************************************************
!     calculate left eigenvector for properties *
!************************************************

      if (exa_input%eom_properties) call sigma_left_ip(exa_input,int_t,eom_t,t2_tensor)

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(eom_t%r1_coeff)
      ierr=talsh_tensor_destruct(eom_t%r1_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l1_coeff)
        ierr=talsh_tensor_destruct(eom_t%l1_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%r2_coeff)
      ierr=talsh_tensor_destruct(eom_t%r2_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l2_coeff)
        ierr=talsh_tensor_destruct(eom_t%l2_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%fbar_oo)
      ierr=talsh_tensor_destruct(eom_t%fbar_ov)
      ierr=talsh_tensor_destruct(eom_t%fbar_vv)
      ierr=talsh_tensor_destruct(eom_t%w_oooo)
      ierr=talsh_tensor_destruct(eom_t%w_ooov)
      ierr=talsh_tensor_destruct(eom_t%w_ovoo)
      ierr=talsh_tensor_destruct(eom_t%w_vovo)
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

    end subroutine eom_ip

    subroutine sigma_right_ip(exa_input,int_t,eom_t,t2_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: v_aux
!     error code
      integer :: ierr

!***************
!     r1_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%r1_sigma)

!--------------------------------------
!     term 1: S1(i) -= F(m,i) * R1(m) |
!--------------------------------------

      ierr=talsh_tensor_contract("S(i)+=F(m,i)*R(m)",eom_t%r1_sigma,eom_t%fbar_oo,eom_t%r1_coeff,scale=MINUS_ONE)

!------------------------------------------
!     term 2: S1(i) += F(m,e) * R2(e,i,m) |
!------------------------------------------

      ierr=talsh_tensor_contract("S(i)+=F(m,e)*R(e,i,m)",eom_t%r1_sigma,eom_t%fbar_ov,eom_t%r2_coeff)

!--------------------------------------------------
!     term 3: S1(i) += 1/2 W(m,n,i,e) * R2(e,m,n) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(i)+=W(m,n,i,e)*R(e,m,n)",eom_t%r1_sigma,eom_t%w_ooov,eom_t%r2_coeff,scale=ONE_HALF)

!***************
!     r2_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%r2_sigma)

!---------------------------------------------
!     term 1: S2(a,i,j) -= W(m,a,i,j) * R(m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=W(m,a,j,i)*R(m)",eom_t%r2_sigma,eom_t%w_ovoo,eom_t%r1_coeff)

!-----------------------------------------------------
!     term 2: S2(a,i,j) -= P(i,j) F(m,i) * R2(a,m,j) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=F(m,i)*r(a,j,m)",eom_t%r2_sigma,eom_t%fbar_oo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,i,j)+=F(m,j)*r(a,m,i)",eom_t%r2_sigma,eom_t%fbar_oo,eom_t%r2_coeff)

!--------------------------------------------------
!     term 3: S2(a,i,j) += W(m,n,i,j) * R2(a,m,n) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=W(m,n,i,j)*r(a,m,n)",eom_t%r2_sigma,eom_t%w_oooo,eom_t%r2_coeff)

!---------------------------------------------------------
!     term 4: S2(a,i,j) -= P(i,j) W(m,a,i,e) * R2(e,m,j) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=W(a,m,e,i)*R(e,j,m)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,i,j)+=W(a,m,e,j)*R(e,m,i)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)

!----------------------------------------------
!     term 5: S2(a,i,j) += F(a,e) * R2(e,i,j) |
!----------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=F(a,e)*R(e,i,j)",eom_t%r2_sigma,eom_t%fbar_vv,eom_t%r2_coeff)

!-----------------------------------------------------------------------
!     term 6: S2(a,i,j) += 1/2 ( V(m,n,e,f) * R2(f,m,n) ) * T(e,a,i,j) |
!-----------------------------------------------------------------------

      ierr=talsh_tensor_construct(v_aux,C8,(/exa_input%nvir/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(e)+=V(m,n,e,f)*R(f,m,n)",v_aux,int_t%oovv,eom_t%r2_coeff,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(a,i,j)+=U(e)*T(e,a,i,j)",eom_t%r2_sigma,v_aux,t2_tensor)
      ierr=talsh_tensor_destruct(v_aux)

      ! debug
      if (debug_exa) then
        print*, "r1_sigma", tensor_norm2(eom_t%r1_sigma)
        print*, "r2_sigma", tensor_norm2(eom_t%r2_sigma)*ONE_QUARTER
      end if

    end subroutine sigma_right_ip

    subroutine sigma_left_ip(exa_input,int_t,eom_t,t2_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: v_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%l1_sigma)

!-------------------------------------
!     term 1: S1(i) -= F(i,m) * L(m) |
!-------------------------------------

      ierr=talsh_tensor_contract("S(i)+=F(i,m)*L(m)",eom_t%l1_sigma,eom_t%fbar_oo,eom_t%l1_coeff,scale=MINUS_ONE)

!-------------------------------------------------
!     term 2: S1(i) += 1/2 W(i,e,m,n) * L(m,n,e) |
!-------------------------------------------------

      ierr=talsh_tensor_contract("S(i)+=W(i,e,m,n)*L(m,n,e)",eom_t%l1_sigma,eom_t%w_ovoo,eom_t%l2_coeff,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%l2_sigma)

!---------------------------------------------
!     term 1: S2(i,j,a) -= W(i,j,m,a) * L(m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=W(j,i,m,a)*L(m)",eom_t%l2_sigma,eom_t%w_ooov,eom_t%l1_coeff)

!---------------------------------------------
!     term 2: S2(i,j,a) += F(e,a) * L(i,j,e) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=F(e,a)*L(i,j,e)",eom_t%l2_sigma,eom_t%fbar_vv,eom_t%l2_coeff)

!----------------------------------------------------
!     term 3: S2(i,j,a) -= P(i,j) F(j,m) * L(i,m,a) |
!----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=F(j,m)*L(m,i,a)",eom_t%l2_sigma,eom_t%fbar_oo,eom_t%l2_coeff)
      ierr=talsh_tensor_contract("S(i,j,a)+=F(i,m)*L(j,m,a)",eom_t%l2_sigma,eom_t%fbar_oo,eom_t%l2_coeff)

!--------------------------------------------------------
!     term 4: S2(i,j,a) += P(i,j) W(j,e,a,m) * L(i,m,e) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=W(e,j,a,m)*L(m,i,e)",eom_t%l2_sigma,eom_t%w_vovo,eom_t%l2_coeff)
      ierr=talsh_tensor_contract("S(i,j,a)+=W(e,i,a,m)*L(j,m,e)",eom_t%l2_sigma,eom_t%w_vovo,eom_t%l2_coeff)

!-----------------------------------------------------
!     term 5: S2(i,j,a) += 1/2 W(i,j,m,n) * L(m,n,a) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=W(i,j,m,n)*L(m,n,a)",eom_t%l2_sigma,eom_t%w_oooo,eom_t%l2_coeff,scale=ONE_HALF)

!----------------------------------------------------------------------
!     term 6: S2(i,j,a) += 1/2 V(i,j,e,a) * ( L(m,n,f) * T(e,f,m,n) ) |
!----------------------------------------------------------------------

      ierr=talsh_tensor_construct(v_aux,C8,(/exa_input%nvir/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(e)+=L(m,n,f)*T(e,f,m,n)",v_aux,eom_t%l2_coeff,t2_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(i,j,a)+=V(i,j,e,a)*U(e)",eom_t%l2_sigma,int_t%oovv,v_aux)
      ierr=talsh_tensor_destruct(v_aux)

!------------------------------------------------
!     term 7: S2(i,j,a) += P(i,j) L(i) * F(j,a) |
!------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=L(i)*F(j,a)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%fbar_ov)
      ierr=talsh_tensor_contract("S(i,j,a)+=L(j)*F(i,a)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%fbar_ov,scale=MINUS_ONE)

    end subroutine sigma_left_ip

end module talsh_eom
