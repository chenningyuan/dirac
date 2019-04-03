module talsh_eom

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates
  use talsh_ao_to_mo, only : tensor_norm2
  use exacorr_global, only : debug_exa
  use exacorr_utils,  only : print_date

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

      write(*,*) ""
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A)') "***************************** Equation of Motion ******************************"
      write(*,'(A)') "*******************************************************************************"
      write(*,*) ""

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      if (.not.exa_input%ccd) then
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EA")
            call eom_ea(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EE")
            if (.not.exa_input%eom_properties) then
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor=t1_tensor)
            else
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)
            end if
          case default
            call quit ('Error: No valid EOM-flavour chosen')
        end select
      else
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EA")
            call eom_ea(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EE")
            if (.not.exa_input%eom_properties) then
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
            else
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor)
            end if
        case default
          call quit ('Error: No valid EOM-flavour chosen')
        end select
      end if

!     make some noise so that we know we are leaving
      print*, ""
      call print_date("Leaving eom_driver routine")
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
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,(/nocc/),init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,(/nocc/),init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,(/nocc/),init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,(/nocc/),init_val=ZERO)
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
      call print_date("Calculated EOM-IP RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

      if (exa_input%eom_properties) then
        call sigma_left_ip(exa_input,int_t,eom_t,t2_tensor)
        call print_date("Calculated EOM-IP LHS")
      end if

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

    end subroutine eom_ip

    subroutine eom_ea(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)

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
      integer, dimension(3) :: vvo_dims, ovv_dims
      integer, dimension(4) :: vovo_dims, vovv_dims, vvvo_dims, vvvv_dims
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
      vvo_dims = (/nvir,nvir,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      vovv_dims = (/nvir,nocc,nvir,nvir/)
      vvvo_dims = (/nvir,nvir,nvir,nocc/)
      vvvv_dims = (/nvir,nvir,nvir,nvir/)
      if (exa_input%eom_properties) ovv_dims = (/nocc,nvir,nvir/)

!     create all fixed intermediates
      ierr=talsh_tensor_construct(eom_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovv,C8,vovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

      call get_wbar_vovo(eom_t%wbar_vovo,t2_tensor,int_t)
      if (.not.CCD) then
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv,t1_tensor)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor,t1_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t,t1_tensor)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor,t1_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor)
      end if
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,(/nvir/),init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,(/nvir/),init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,(/nvir/),init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,(/nvir/),init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%r2_coeff,C8,vvo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r2_sigma,C8,vvo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l2_coeff,C8,ovv_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l2_sigma,C8,ovv_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      call sigma_right_ea(exa_input,int_t,eom_t,t2_tensor)
      call print_date("Calculated EOM-EA RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

      if (exa_input%eom_properties) then
        call sigma_left_ea(exa_input,int_t,eom_t,t2_tensor)
        call print_date("Calculated EOM-EA LHS")
      end if

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
      ierr=talsh_tensor_destruct(eom_t%w_vovo)
      ierr=talsh_tensor_destruct(eom_t%w_vovv)
      ierr=talsh_tensor_destruct(eom_t%w_vvvo)
      ierr=talsh_tensor_destruct(eom_t%w_vvvv)

    end subroutine eom_ea

    subroutine eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: l2_tensor,t1_tensor
!     tensor dims
      integer :: nocc, nvir
      integer, dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
      integer, dimension(4) :: oooo_dims, ooov_dims, ovoo_dims, oovv_dims, vvoo_dims, vovo_dims, vovv_dims, vvvo_dims, vvvv_dims
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
      vo_dims = (/nvir,nocc/)
      vv_dims = (/nvir,nvir/)
      oooo_dims = (/nocc,nocc,nocc,nocc/)
      ooov_dims = (/nocc,nocc,nocc,nvir/)
      ovoo_dims = (/nocc,nvir,nocc,nocc/)
      if (exa_input%eom_properties) oovv_dims = (/nocc,nocc,nvir,nvir/)
      vvoo_dims = (/nvir,nvir,nocc,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      vovv_dims = (/nvir,nocc,nvir,nvir/)
      vvvo_dims = (/nvir,nvir,nvir,nocc/)
      vvvv_dims = (/nvir,nvir,nvir,nvir/)

!     create all fixed intermediates
      ierr=talsh_tensor_construct(eom_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%goo,C8,oo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(eom_t%gvv,C8,vv_dims,init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%w_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ooov,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovv,C8,vovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
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
        call get_w_vovv (eom_t%w_vovv,int_t,t1_tensor)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor,t1_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor)
      end if
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,vo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,vo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,ov_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,ov_dims,init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%r2_coeff,C8,vvoo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r2_sigma,C8,vvoo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l2_coeff,C8,oovv_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l2_sigma,C8,oovv_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      call sigma_right_ee(exa_input,int_t,eom_t,t2_tensor)
      call print_date("Calculated EOM-EE RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

!     create mixed type intermediates
      if (exa_input%eom_properties) then
        call get_goo_lambda(eom_t%goo,l2_tensor,t2_tensor)
        call get_gvv_lambda(eom_t%gvv,l2_tensor,t2_tensor)
      end if

      if (exa_input%eom_properties) then
        if (.not.CCD) then
          call sigma_left_ee(exa_input,int_t,eom_t,t1_tensor)
        else
          call sigma_left_ee(exa_input,int_t,eom_t)
        end if
        call print_date("Calculated EOM-EE LHS")
      end if

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
        ierr=talsh_tensor_destruct(eom_t%goo)
        ierr=talsh_tensor_destruct(eom_t%gvv)
      end if
      ierr=talsh_tensor_destruct(eom_t%fbar_oo)
      ierr=talsh_tensor_destruct(eom_t%fbar_ov)
      ierr=talsh_tensor_destruct(eom_t%fbar_vv)
      ierr=talsh_tensor_destruct(eom_t%w_oooo)
      ierr=talsh_tensor_destruct(eom_t%w_ooov)
      ierr=talsh_tensor_destruct(eom_t%w_ovoo)
      ierr=talsh_tensor_destruct(eom_t%w_vovo)
      ierr=talsh_tensor_destruct(eom_t%w_vovv)
      ierr=talsh_tensor_destruct(eom_t%w_vvvo)
      ierr=talsh_tensor_destruct(eom_t%w_vvvv)

    end subroutine eom_ee

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

    subroutine sigma_right_ea(exa_input,int_t,eom_t,t2_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: o_aux
!     error code
      integer :: ierr

!***************
!     r1_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%r1_sigma)

!--------------------------------------
!     term 1: S1(a) += F(e,a) * R1(e) |
!--------------------------------------

      ierr=talsh_tensor_contract("S(a)+=F(e,a)*R(e)",eom_t%r1_sigma,eom_t%fbar_vv,eom_t%r1_coeff)

!------------------------------------------
!     term 2: S1(a) += F(m,e) * R2(e,a,m) |
!------------------------------------------

      ierr=talsh_tensor_contract("S(a)+=F(m,e)*R(e,a,m)",eom_t%r1_sigma,eom_t%fbar_ov,eom_t%r2_coeff)

!--------------------------------------------------
!     term 3: S1(a) += 1/2 W(a,m,f,e) * R2(e,f,m) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a)+=W(a,m,f,e)*R(e,f,m)",eom_t%r1_sigma,eom_t%w_vovv,eom_t%r2_coeff,scale=ONE_HALF)

!***************
!     r2_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%r2_sigma)

!---------------------------------------------
!     term 1: S2(a,b,i) += W(a,b,e,i) * R(e) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i)+=W(a,b,e,i)*R(m)",eom_t%r2_sigma,eom_t%w_vvvo,eom_t%r1_coeff)

!-----------------------------------------------------
!     term 2: S2(a,b,i) += P(a,b) F(a,e) * R2(e,b,i) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i)+=F(a,e)*r(e,b,i)",eom_t%r2_sigma,eom_t%fbar_vv,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i)+=F(b,e)*r(a,e,i)",eom_t%r2_sigma,eom_t%fbar_vv,eom_t%r2_coeff)

!--------------------------------------------------
!     term 3: S2(a,b,i) += W(a,b,e,f) * R2(e,f,i) |
!--------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i)+=W(a,b,e,f)*r(e,f,i)",eom_t%r2_sigma,eom_t%w_vvvv,eom_t%r2_coeff,scale=ONE_HALF)

!---------------------------------------------------------
!     term 4: S2(a,b,i) -= P(a,b) W(a,m,e,i) * R2(e,b,m) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i)+=W(a,m,e,i)*R(b,e,m)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i)+=W(b,m,e,i)*R(e,a,m)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)

!----------------------------------------------
!     term 5: S2(a,b,i) -= F(m,i) * R2(a,b,m) |
!----------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i)+=F(m,i)*R(b,a,m)",eom_t%r2_sigma,eom_t%fbar_oo,eom_t%r2_coeff)

!-----------------------------------------------------------------------
!     term 6: S2(a,b,i) += 1/2 ( V(m,n,e,f) * R2(e,f,n) ) * T(a,b,m,i) |
!-----------------------------------------------------------------------

      ierr=talsh_tensor_construct(o_aux,C8,(/exa_input%nocc/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(m)+=V(m,n,e,f)*R(e,f,n)",o_aux,int_t%oovv,eom_t%r2_coeff,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(a,b,i)+=U(m)*T(a,b,m,i)",eom_t%r2_sigma,o_aux,t2_tensor)
      ierr=talsh_tensor_destruct(o_aux)

      ! debug
      if (debug_exa) then
        print*, "r1_sigma", tensor_norm2(eom_t%r1_sigma)
        print*, "r2_sigma", tensor_norm2(eom_t%r2_sigma)*ONE_QUARTER
      end if

    end subroutine sigma_right_ea

    subroutine sigma_left_ea(exa_input,int_t,eom_t,t2_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: o_aux
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%l1_sigma)

!-------------------------------------
!     term 1: S1(a) += F(e,a) * L(e) |
!-------------------------------------

      ierr=talsh_tensor_contract("S(a)+=F(e,a)*L(e)",eom_t%l1_sigma,eom_t%fbar_vv,eom_t%l1_coeff)

!-------------------------------------------------
!     term 2: S1(a) -= 1/2 W(e,f,a,m) * L(m,f,e) |
!-------------------------------------------------

      ierr=talsh_tensor_contract("S(a)+=W(e,f,a,m)*L(m,e,f)",eom_t%l1_sigma,eom_t%w_vvvo,eom_t%l2_coeff,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%l2_sigma)

!---------------------------------------------
!     term 1: S2(i,a,b) += W(e,i,b,a) * L(e) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a,b)+=W(e,i,b,a)*L(e)",eom_t%l2_sigma,eom_t%w_vovv,eom_t%l1_coeff)

!----------------------------------------------------
!     term 2: S2(i,a,b) += P(a,b) F(e,a) * L(i,e,b) |
!----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a,b)+=F(e,a)*L(i,e,b)",eom_t%l2_sigma,eom_t%fbar_vv,eom_t%l2_coeff)
      ierr=talsh_tensor_contract("S(i,a,b)+=F(e,b)*L(i,a,e)",eom_t%l2_sigma,eom_t%fbar_vv,eom_t%l2_coeff)

!---------------------------------------------
!     term 3: S2(i,a,b) -= F(i,m) * L(m,a,b) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a,b)+=F(i,m)*L(m,b,a)",eom_t%l2_sigma,eom_t%fbar_oo,eom_t%l2_coeff)

!-----------------------------------------------------
!     term 4: S2(i,a,b) += 1/2 W(e,f,a,b) * L(i,e,f) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a,b)+=W(e,f,a,b)*L(i,e,f)",eom_t%l2_sigma,eom_t%w_vvvv,eom_t%l2_coeff,scale=ONE_HALF)

!--------------------------------------------------------
!     term 5: S2(i,a,b) -= P(a,b) W(e,i,b,m) * L(m,a,e) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a,b)+=W(e,i,b,m)*L(m,e,a)",eom_t%l2_sigma,eom_t%w_vovo,eom_t%l2_coeff)
      ierr=talsh_tensor_contract("S(i,a,b)+=W(e,i,a,m)*L(m,b,e)",eom_t%l2_sigma,eom_t%w_vovo,eom_t%l2_coeff)

!----------------------------------------------------------------------
!     term 6: S2(i,a,b) -= 1/2 ( L(o,e,f) * T(e,f,m,o) ) * V(m,i,b,a) |
!----------------------------------------------------------------------

      ierr=talsh_tensor_construct(o_aux,C8,(/exa_input%nocc/),init_val=ZERO)
      ierr=talsh_tensor_contract("U(m)+=L(o,e,f)*T(e,f,o,m)",o_aux,eom_t%l2_coeff,t2_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(i,a,b)+=U(m)*V(m,i,b,a)",eom_t%l2_sigma,o_aux,int_t%oovv)
      ierr=talsh_tensor_destruct(o_aux)

!------------------------------------------------
!     term 7: S2(i,a,b) += P(a,b) F(i,a) * L(b) |
!------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a,b)+=F(i,a)*L(b)",eom_t%l2_sigma,eom_t%fbar_ov,eom_t%l1_coeff)
      ierr=talsh_tensor_contract("S(i,a,b)+=F(i,b)*L(a)",eom_t%l2_sigma,eom_t%fbar_ov,eom_t%l1_coeff,scale=MINUS_ONE)

    end subroutine sigma_left_ea

    subroutine sigma_right_ee(exa_input,int_t,eom_t,t2_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: oo_aux, vv_aux
      integer            :: oo_dims(2), vv_dims(2)
!     error code
      integer :: ierr

      oo_dims = (/exa_input%nocc,exa_input%nocc/)
      vv_dims = (/exa_input%nvir,exa_input%nvir/)

!***************
!     r1_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%r1_sigma)

!-----------------------------------------
!     term 1: S1(a,i) += F(a,e) * R(e,i) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=F(a,e)*R(e,i)",eom_t%r1_sigma,eom_t%fbar_vv,eom_t%r1_coeff)

!-----------------------------------------
!     term 2: S1(a,i) -= F(m,i) * R(a,m) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=F(m,i)*R(a,m)",eom_t%r1_sigma,eom_t%fbar_oo,eom_t%r1_coeff,scale=MINUS_ONE)

!---------------------------------------------
!     term 3: S1(a,i) += F(m,e) * R(e,a,m,i) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=F(m,e)*R(e,a,m,i)",eom_t%r1_sigma,eom_t%fbar_ov,eom_t%r2_coeff)

!---------------------------------------------
!     term 4: S1(a,i) += W(m,a,e,i) * R(e,m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=W(a,m,e,i)*R(e,m)",eom_t%r1_sigma,eom_t%w_vovo,eom_t%r1_coeff,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 5: S1(a,i) += 1/2 W(a,m,e,f) * R(e,f,i,m) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=W(a,m,e,f)*R(e,f,i,m)",eom_t%r1_sigma,eom_t%w_vovv,eom_t%r2_coeff,scale=ONE_HALF)

!-----------------------------------------------------
!     term 6: S1(a,i) -= 1/2 W(m,n,i,e) * R(e,a,m,n) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(a,i)+=W(m,n,i,e)*R(a,e,m,n)",eom_t%r1_sigma,eom_t%w_ooov,eom_t%r2_coeff,scale=ONE_HALF)

!***************
!     r2_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%r2_sigma)

!--------------------------------------------------------
!     term 1: S2(a,b,i,j) -= P(a,b) W(m,b,i,j) * R(a,m) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,b,j,i)*R(a,m)",eom_t%r2_sigma,eom_t%w_ovoo,eom_t%r1_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,a,i,j)*R(b,m)",eom_t%r2_sigma,eom_t%w_ovoo,eom_t%r1_coeff)

!--------------------------------------------------------
!     term 2: S2(a,b,i,j) += P(i,j) W(a,b,e,j) * R(e,i) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,b,e,j)*R(e,i)",eom_t%r2_sigma,eom_t%w_vvvo,eom_t%r1_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,a,e,i)*R(e,j)",eom_t%r2_sigma,eom_t%w_vvvo,eom_t%r1_coeff)

!-------------------------------------------------------------------------
!     term 3: S2(a,b,i,j) += P(a,b) ( W(b,m,f,e) * R(e,m) ) * T(a,f,i,j) |
!-------------------------------------------------------------------------

      ierr=talsh_tensor_construct(vv_aux,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(b,f)+=W(b,m,f,e)*R(e,m)",vv_aux,eom_t%w_vovv,eom_t%r1_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(b,f)*T2(a,f,i,j)",eom_t%r2_sigma,vv_aux,t2_tensor)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(a,f)*T2(f,b,i,j)",eom_t%r2_sigma,vv_aux,t2_tensor)
      ierr=talsh_tensor_destruct(vv_aux)

!-------------------------------------------------------------------------
!     term 4: S2(a,b,i,j) += P(i,j) ( W(n,m,j,e) * R(e,m) ) * T(a,b,i,n) |
!-------------------------------------------------------------------------

      ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(n,j)+=W(n,m,j,e)*R(e,m)",oo_aux,eom_t%w_ooov,eom_t%r1_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,j)*T(a,b,i,n)",eom_t%r2_sigma,oo_aux,t2_tensor)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,i)*T(a,b,n,j)",eom_t%r2_sigma,oo_aux,t2_tensor)
      ierr=talsh_tensor_destruct(oo_aux)

!--------------------------------------------------------
!     term 5: S2(a,b,i,j) += P(a,b) F(b,e) * R(a,e,i,j) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=F(b,e)*R(a,e,i,j)",eom_t%r2_sigma,eom_t%fbar_vv,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=F(a,e)*R(e,b,i,j)",eom_t%r2_sigma,eom_t%fbar_vv,eom_t%r2_coeff)

!--------------------------------------------------------
!     term 6: S2(a,b,i,j) -= P(i,j) F(m,j) * R(a,b,i,m) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,j)*R(a,b,m,i)",eom_t%r2_sigma,eom_t%fbar_oo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,i)*R(a,b,j,m)",eom_t%r2_sigma,eom_t%fbar_oo,eom_t%r2_coeff)

!---------------------------------------------------------
!     term 7: S2(a,b,i,j) += 1/2 W(m,n,i,j) * R(a,b,m,n) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(m,n,i,j)*R(a,b,m,n)",eom_t%r2_sigma,eom_t%w_oooo,eom_t%r2_coeff,scale=ONE_HALF)

!-------------------------------------------------------------------
!     term 8: S2(a,b,i,j) += P(i,j) P(a,b) W(m,b,e,j) * R(a,e,i,m) |
!-------------------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,m,e,j)*R(a,e,m,i)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,m,e,j)*R(b,e,i,m)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(b,m,e,i)*R(a,e,j,m)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,m,e,i)*R(b,e,m,j)",eom_t%r2_sigma,eom_t%w_vovo,eom_t%r2_coeff)

!---------------------------------------------------------------------------------
!     term 9: S2(a,b,i,j) -= 1/2 P(a,b) ( V(n,m,f,e) * R(e,a,m,n) ) * T(f,b,i,j) |
!---------------------------------------------------------------------------------

      ierr=talsh_tensor_construct(vv_aux,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(a,f)+=V(n,m,f,e)*R(e,a,m,n)",vv_aux,int_t%oovv,eom_t%r2_coeff,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(a,f)*T(f,b,j,i)",eom_t%r2_sigma,vv_aux,t2_tensor)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(b,f)*T(f,a,i,j)",eom_t%r2_sigma,vv_aux,t2_tensor)
      ierr=talsh_tensor_destruct(vv_aux)

!---------------------------------------------------------------------------------
!     term 10: S2(a,b,i,j) -= 1/2 P(i,j) ( V(n,m,f,e) * R(f,e,i,m) ) * T(b,a,j,n) |
!---------------------------------------------------------------------------------

      ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_contract("U(n,i)+=V(n,m,f,e)*R(f,e,i,m)",oo_aux,int_t%oovv,eom_t%r2_coeff,scale=ONE_HALF)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,i)*T(b,a,n,j)",eom_t%r2_sigma,oo_aux,t2_tensor)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=U(n,j)*T(b,a,i,n)",eom_t%r2_sigma,oo_aux,t2_tensor)
      ierr=talsh_tensor_destruct(oo_aux)

!----------------------------------------------------------
!     term 11: S2(a,b,i,j) += 1/2 W(a,b,e,f) * R(e,f,i,j) |
!----------------------------------------------------------

      ierr=talsh_tensor_contract("S(a,b,i,j)+=W(a,b,e,f)*R(e,f,i,j)",eom_t%r2_sigma,eom_t%w_vvvv,eom_t%r2_coeff,scale=ONE_HALF)

    end subroutine sigma_right_ee

    subroutine sigma_left_ee(exa_input,int_t,eom_t,t1_tensor)

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!     auxiliary tensor
      type(talsh_tens_t) :: oo_aux
      integer            :: oo_dims(2)
!     error code
      integer :: ierr

!***************
!     l1_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%l1_sigma)

!--------------------------------
!     term 1: S1(i,a) += F(i,a) |
!--------------------------------

      ierr=talsh_tensor_add("S(i,a)+=F(i,a)",eom_t%l1_sigma,eom_t%fbar_ov)

!-----------------------------------------
!     term 2: S1(i,a) += L(i,e) * F(e,a) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",eom_t%l1_sigma,eom_t%l1_coeff,eom_t%fbar_vv)

!-----------------------------------------
!     term 3: S1(i,a) -= L(m,a) * F(i,m) |
!-----------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",eom_t%l1_sigma,eom_t%l1_coeff,eom_t%fbar_oo,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 4: S1(i,a) += 1/2 L(i,m,e,f) * W(e,f,a,m) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",eom_t%l1_sigma,eom_t%l2_coeff,eom_t%w_vvvo,scale=ONE_HALF)

!---------------------------------------------
!     term 5: S1(i,a) -= G(f,e) * W(e,i,f,a) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",eom_t%l1_sigma,eom_t%gvv,eom_t%w_vovv)

!---------------------------------------------
!     term 6: S1(i,a) -= G(n,m) * W(m,i,n,a) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",eom_t%l1_sigma,eom_t%goo,eom_t%w_ooov)

!---------------------------------------------
!     term 7: S1(i,a) += L(m,e) * W(i,e,a,m) |
!---------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",eom_t%l1_sigma,eom_t%l1_coeff,eom_t%w_vovo,scale=MINUS_ONE)

!-----------------------------------------------------
!     term 8: S1(i,a) -= 1/2 L(m,n,a,e) * W(i,e,m,n) |
!-----------------------------------------------------

      ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",eom_t%l1_sigma,eom_t%l2_coeff,eom_t%w_ovoo,scale=ONE_HALF)

!***************
!     l2_sigma *
!***************

      ierr=talsh_tensor_init(eom_t%l2_sigma)

!----------------------------------------
!     term 1: S2(i,j,a,b) += V(i,j,a,b) |
!----------------------------------------

      ierr=talsh_tensor_add("S(i,j,a,b)+=V(i,j,a,b)",eom_t%l2_sigma,int_t%oovv)

!--------------------------------------------------------
!     term 2: S2(i,j,a,b) += P(a,b) L(i,j,a,e) * F(e,b) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%fbar_vv)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%fbar_vv)

!--------------------------------------------------------
!     term 3: S2(i,j,a,b) -= P(i,j) L(i,m,a,b) * F(j,m) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,i,a,b)*F(j,m)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%fbar_oo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%fbar_oo)

!---------------------------------------------------------
!     term 4: S2(i,j,a,b) += 1/2 L(m,n,a,b) * W(i,j,m,n) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%w_oooo,scale=ONE_HALF)

!-------------------------------------------------------------------
!     term 5: S2(i,j,a,b) += P(i,j) P(a,b) L(i,m,a,e) * W(j,e,b,m) |
!-------------------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%w_vovo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%w_vovo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%w_vovo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%w_vovo)

!--------------------------------------------------------
!     term 6: S2(i,j,a,b) += P(a,b) V(i,j,a,e) * G(e,b) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",eom_t%l2_sigma,int_t%oovv,eom_t%gvv)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",eom_t%l2_sigma,int_t%oovv,eom_t%gvv)

!-------------------------------------------------
!     term 7: S2(i,j,a,b) -= L(m,a) * W(i,j,m,b) |
!-------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%w_ooov)

!--------------------------------------------------------
!     term 8: S2(i,j,a,b) -= P(i,j) V(i,m,a,b) * G(j,m) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,i,a,b)*G(j,m)",eom_t%l2_sigma,int_t%oovv,eom_t%goo)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",eom_t%l2_sigma,int_t%oovv,eom_t%goo)

!---------------------------------------------------------------------
!     term 9: S2(i,j,a,b) += P(i,j) V(m,j,a,b) * ( L(i,e) * T(e,m) ) |
!---------------------------------------------------------------------

      if (.not.exa_input%ccd) then
        oo_dims = (/exa_input%nocc,exa_input%nocc/)
        ierr=talsh_tensor_construct(oo_aux,C8,oo_dims,init_val=ZERO)
        ierr=talsh_tensor_contract("U(i,m)+=L(i,e)*T(e,m)",oo_aux,eom_t%l1_coeff,t1_tensor)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(m,j,a,b)*U(i,m)",eom_t%l2_sigma,int_t%oovv,oo_aux)
        ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,m,a,b)*U(j,m)",eom_t%l2_sigma,int_t%oovv,oo_aux)
        ierr=talsh_tensor_destruct(oo_aux)
      end if

!---------------------------------------------------------
!     term 10: S2(i,j,a,b) += P(i,j) L(i,e) * V(e,j,a,b) |
!---------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*V(e,j,a,b)",eom_t%l2_sigma,eom_t%l1_coeff,int_t%vovv)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*V(e,i,b,a)",eom_t%l2_sigma,eom_t%l1_coeff,int_t%vovv)

!------------------------------------------------------------
!     term 11: S2(i,j,a,b) += P(i,j) P(a,b) L(i,a) * F(j,b) |
!------------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%fbar_ov)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%fbar_ov,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%fbar_ov,scale=MINUS_ONE)
      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",eom_t%l2_sigma,eom_t%l1_coeff,eom_t%fbar_ov)

!--------------------------------------------------------
!     term 12: S2(i,j,a,b) += 1/2 L(i,j,e,f) W(e,f,a,b) |
!--------------------------------------------------------

      ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",eom_t%l2_sigma,eom_t%l2_coeff,eom_t%w_vvvv,scale=ONE_HALF)

    end subroutine sigma_left_ee

end module talsh_eom