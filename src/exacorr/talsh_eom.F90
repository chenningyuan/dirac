module talsh_eom

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates

  implicit none

  complex(8), parameter :: MINUS_ONE=(-1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)

  private

  public talsh_eom_driver

  contains

    subroutine talsh_eom_driver (exa_input,t2_tensor,int_t,t1_tensor,l2_tensor)

!     input variables
      type(exacc_input),     intent(in   )           :: exa_input
      type(talsh_tens_t),    intent(inout)           :: t2_tensor
      type(talsh_intg_tens), intent(inout)           :: int_t
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor, l2_tensor
!     tensor dims
      ! integer,dimension(2)       :: oo_id, ov_id, vv_id
      ! integer,dimension(4)       :: oooo_id, ovoo_id, ooov_id, vovo_id, vvvo_id, vovv_id, vvvv_id
!     error code
      integer :: ierr

!     make some noise so that we know we are here
      call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
      print*, " Starting EOM-CC calculation with talsh"

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      select case(exa_input%eom_flavour)
        if (.not.exa_input%ccd)
          case("IP")
            call eom_ip(exa_input,t2_tensor,int_t,t1_tensor)
          case("EA")
            ! call eom_ea(exa_input,t2_tensor,int_t,t1_tensor)
          case("EE")
            ! call eom_ee(exa_input,t2_tensor,int_t,t1_tensor,l2_tensor)
        else
          case("IP")
            call eom_ip(exa_input,t2_tensor,int_t)
          case("EA")
            ! call eom_ea(exa_input,t2_tensor,int_t)
          case("EE")
            ! call eom_ee(exa_input,t2_tensor,int_t)
        end if
        case default
          call quit ('Error: No valid EOM-flavour chosen')
      end select

!     make some noise so that we know we are leaving
      print*, " Leaving eom_driver routine"
      return

    end subroutine exacorr_eom_driver

    subroutine eom_ip(exa_input,t2_tensor,int_t,t1_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!     eom tensors
      type(talsh_tens_t) :: r1_coeff, r2_coeff, l1_coeff, l2_coeff
      type(talsh_tens_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!     fixed intermediates
      type(exatns_lambda_tens) :: lambda_t    
!     tensor dims
      integer, dimension(3) :: voo_dims, oov_dims
!     CCD switch
      logical :: CCD

      CCD = exa_input%ccd

!************************
!     construct tensors *
!************************

!     set tensor dims
      voo_dims = (exa_input%nvir,exa_input%nocc,exa_input%nocc)
      if (exa_input%eom_properties) oov_dims = (exa_input%nocc,exa_input%nocc,exa_input%nvir)

!     create all coefficient and sigma tensors
      if (.not.CCD) then
        ierr=talsh_tensor_construct(r1_coeff,C8,(/exa_input%nocc/),init_val=ZERO)
        ierr=talsh_tensor_construct(r1_sigma,C8,(/exa_input%nocc/),init_val=ZERO)
        if (exa_input%eom_properties) then
          ierr=talsh_tensor_construct(l1_coeff,C8,(/exa_input%nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(l1_sigma,C8,(/exa_input%nvir/),init_val=ZERO)
        end if
      end if
      ierr=talsh_tensor_construct(r2_coeff,C8,voo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(r2_sigma,C8,voo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(l2_coeff,C8,oov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(l2_sigma,C8,oov_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      if (.not.CCD)
        call sigma_right_ip(exa_input,int_t,t2_tensor,r2_coeff,r2_sigma,t1_tensor,r1_coeff,r1_sigma)
      else
        call sigma_right_ip(exa_input,int_t,t2_tensor,r2_coeff,r2_sigma)
      end if

!************************************************
!     calculate left eigenvector for properties *
!************************************************

      if (exa_input%eom_properties) then
        if (.not.CCD) then
          call sigma_left_ip(exa_input,int_t,t2_tensor,l2_coeff,l2_sigma,t1_tensor,l1_coeff,l1_sigma)
        else
          call sigma_left_ip(exa_input,int_t,t2_tensor,l2_coeff,l2_sigma)
        end if
      end if

!**************
!     cleanup *
!**************

      if (.not.CCD) then
        ierr=talsh_tensor_destruct(r1_coeff)
        ierr=talsh_tensor_destruct(r1_sigma)
        if (exa_input%eom_properties) then
          ierr=talsh_tensor_destruct(l1_coeff)
          ierr=talsh_tensor_destruct(l1_sigma)
        end if
      end if
      ierr=talsh_tensor_destruct(r2_coeff)
      ierr=talsh_tensor_destruct(r2_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(l2_coeff)
        ierr=talsh_tensor_destruct(l2_sigma)
      end if

    end subroutine eom_ip

    subroutine sigma_right_ip(exa_input,int_t,t2_tensor,r2_coeff,r2_sigma,t1_tensor,r1_coeff,r1_sigma)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor, r2_coeff, r2_sigma
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor,r1_coeff, r1_sigma
!     auxiliary tensor
      type(talsh_tens_t) :: v_aux
!     error code
      integer :: ierr

!***************
!     r1_sigma *
!***************

      ierr=talsh_tensor_init(r1_sigma,init_val=ZERO)

!--------------------------------------
!     term 1: S1(i) -= F(m,i) * R1(m) |
!--------------------------------------

      ierr=talsh_tensor_contract("S(i)+=F(m,i)*R(m)",r1_sigma,lambda_t%fbar_oo,r1_coeff,scale=MINUS_ONE)

!------------------------------------------
!     term 2: S1(i) += F(m,e) * R2(e,i,m) |
!------------------------------------------

      ierr=talsh_tensor_contract("S(i)+=F(m,e)*R(e,i,m)",r1_sigma,lambda_t%fbar_ov,scale=r2_coef)

!--------------------------------------------------
!     term 3: S1(i) += 1/2 W(m,n,i,e) * R2(e,m,n) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(i)+=W(m,n,i,e)*R(e,m,n)",r1_sigma,lambda_t%w_ooov,r2_coeff,scale=ONE_HALF)

!***************
!     r2_sigma *
!***************

      ierr=talsh_tensor_init(r2_sigma,init_val=ZERO)

!---------------------------------------------
!     term 1: S2(a,i,j) -= W(m,a,i,j) * R(m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=W(m,a,j,i)*R(m)",r2_sigma,lambda_t%w_ovoo,r1_coef)

!-----------------------------------------------------
!     term 2: S2(a,i,j) -= P(i,j) F(m,i) * R2(a,m,j) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=F(m,i)*r(a,j,m)",r2_sigma,lambda_t%fbar_oo,r2_coef)
      ierr=talsh_tensor_contract("S(a,i,j)+=F(m,j)*r(a,m,i)",r2_sigma,lambda_t%fbar_oo,r2_coef)

!--------------------------------------------------
!     term 3: S2(a,i,j) += W(m,n,i,j) * R2(a,m,n) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=W(m,n,i,j)*r(a,m,n)",r2_sigma,lambda_t%w_oooo,r2_coeff)

!---------------------------------------------------------
!     term 4: S2(a,i,j) -= P(i,j) W(m,a,i,e) * R2(e,m,j) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i,j)+=W(a,m,e,i)*R(e,j,m)",r2_sigma,lambda_t%w_vovo,r2_coeff)
      ierr=talsh_tensor_contract("S(a,i,j)+=W(a,m,e,j)*R(e,m,i)",r2_sigma,lambda_t%w_vovo,r2_coeff)

!----------------------------------------------
!     term 5: S2(a,i,j) += F(a,e) * R2(e,i,j) |
!----------------------------------------------

      ierr=talsh_tensor_conract("S(a,i,j)+=F(a,e)*R(e,i,j)",r2_sigma,lambda_t%fbar_vv,r2_coeff)

!-----------------------------------------------------------------------
!     term 6: S2(a,i,j) += 1/2 ( W(m,n,e,f) * R2(f,m,n) ) * T(e,a,i,j) |
!-----------------------------------------------------------------------

      ierr=talsh_tensor_construct(v_aux,C8,(/nvir/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(e)+=W(m,n,e,f)*R(f,m,n)",v_aux,int_t%oovv,r2_coeff,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(a,i,j)+=U(e)*T(e,a,i,j)",r2_sigma,v_aux,t2_tensor)
      ierr=talsh_tensor_destruct(v_aux)

    end subroutine sigma_right_ip

    subroutine sigma_left_ip(exa_input,int_t,t2_tensor,l2_coeff,l2_sigma,t1_tensor,l1_coeff,l1_sigma)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor, l2_coeff, l2_sigma
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor,l1_coeff, l1_sigma
!     auxiliary tensor
      type(talsh_tens_t) :: v_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(l1_sigma,init_val=ZERO)

!-------------------------------------
!     term 1: S1(i) -= F(i,m) * l(m) |
!-------------------------------------

      ierr=talsh_tensor_contract("S(i)+=F(i,m)*L(m)",l1_sigma,lambda_t%fbar_oo,l1_coeff,scale=MINUS_ONE)

!-------------------------------------------------
!     term 2: S1(i) += 1/2 W(i,e,m,n) * l(m,n,e) |
!-------------------------------------------------

      ierr=talsh_tensor_contract("S(i)+=W(i,e,m,n)*L(m,n,e)",l1_sigma,lambda_t%w_ovoo,l2_coeff,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(l2_sigma,init_val=ZERO)

!---------------------------------------------
!     term 1: S2(i,j,a) -= W(i,j,m,a) * l(m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=W(j,i,m,a)*L(m)",l2_sigma,lambda_t%w_ooov,l1_coeff)

!---------------------------------------------
!     term 2: S2(i,j,a) += F(e,a) * l(i,j,e) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=F(e,a)*L(i,j,e)",l2_sigma,lambda_t%fbar_vv,l2_coeff)

!----------------------------------------------------
!     term 3: S2(i,j,a) -= P(i,j) F(j,m) * l(i,m,a) |
!----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=F(j,m)*L(m,i,a)",l2_sigma,lambda_t%fbar_oo,l2_coeff)
      ierr=talsh_tensor_contract("S(i,j,a)+=F(i,m)*L(j,m,a)",l2_sigma,lambda_t%fbar_oo,l2_coeff)

!--------------------------------------------------------
!     term 4: S2(i,j,a) += P(i,j) W(j,e,a,m) * l(i,m,e) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=W(e,j,a,m)*L(m,i,e)",l2_sigma,lambda_t%w_vovo,l2_coeff)
      ierr=talsh_tensor_contract("S(i,j,a)+=W(e,i,a,m)*L(j,m,e)",l2_sigma,lambda_t%w_vovo,l2_coeff)

!-----------------------------------------------------
!     term 5: S2(i,j,a) += 1/2 W(i,j,m,n) * l(m,n,a) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=W(i,j,m,n)*L(m,n,a)",l2_sigma,lambda_t%w_oooo,l2_coeff,scale=ONE_HALF)

!----------------------------------------------------------------------
!     term 6: S2(i,j,a) += 1/2 V(i,j,e,a) * ( l(m,n,f) * t(e,f,m,n) ) |
!----------------------------------------------------------------------

      ierr=talsh_tensor_construct(v_aux,C8,(/nvir/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(e)+=L(m,n,f)*T(e,f,m,n)",v_aux,l2_coeff,comm_t%t2,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(i,j,a)+=V(i,j,e,a)*U(e)",l2_sigma,int_t%oovv,v_aux)
      ierr=talsh_tensor_destruct(v_aux)

!------------------------------------------------
!     term 7: S2(i,j,a) += P(i,j) l(i) * F(j,a) |
!------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a)+=L(i)*F(j,a)",l2_sigma,l1_coeff,lambda_t%fbar_ov)
      ierr=talsh_tensor_contract("S(i,j,a)+=L(j)*F(i,a)",l2_sigma,l1_coeff,lambda_t%fbar_ov,scale=MINUS_ONE)

    end subroutine sigma_left_ip

end module talsh_eom
