module intermediates_eom_ccsd

!Module to make intermediates needed in CC and Lambda equations

     implicit none
     private

     complex(8), parameter :: ZERO=(0.D0,0.D0), MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0), &
                              ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                              MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0)
    real(8), parameter     :: ONE_QUARTER=0.25D0

     interface get_tau
        module procedure get_tau_talsh
     end interface get_tau

     interface get_fbar_oo
        module procedure get_fbar_oo_talsh_dt
        module procedure get_fbar_oo_talsh_no_dt
     end interface get_fbar_oo

     interface get_fbar_ov
        module procedure get_fbar_ov_talsh
        module procedure get_fbar_ov_talsh_scratch
     end interface get_fbar_ov

     interface get_fbar_vv
        module procedure get_fbar_vv_talsh
        module procedure get_fbar_vv_talsh_scratch
     end interface get_fbar_vv

     interface get_w_oooo
        module procedure get_w_oooo_talsh
        module procedure get_w_oooo_talsh_scratch
     end interface get_w_oooo

     interface get_w_ooov
        module procedure get_w_ooov_talsh
     end interface get_w_ooov

     interface get_w_ovoo
        module procedure get_w_ovoo_talsh
        module procedure get_w_ovoo_talsh_scratch
     end interface get_w_ovoo

     interface get_w_vovo
        module procedure get_w_vovo_talsh
        module procedure get_w_vovo_talsh_scratch
     end interface get_w_vovo

     interface get_wbar_vovo
        module procedure get_wbar_vovo_talsh
     end interface get_wbar_vovo

     interface get_w_vvvo
        module procedure get_w_vvvo_talsh_dt
        module procedure get_w_vvvo_talsh_no_dt
     end interface get_w_vvvo

     interface get_w_vovv
        module procedure get_w_vovv_talsh
     end interface get_w_vovv

     interface get_w_vvvv
        module procedure get_w_vvvv_talsh_dt
        module procedure get_w_vvvv_talsh_no_dt
     end interface get_w_vvvv

     interface get_goo_lambda
        module procedure get_goo_lambda_talsh
     end interface get_goo_lambda

     interface get_gvv_lambda
        module procedure get_gvv_lambda_talsh
     end interface get_gvv_lambda

     interface form_ccsd_intermediates
        module procedure form_ccsd_intermediates_talsh
     end interface form_ccsd_intermediates

     interface destroy_ccsd_intermediates
        module procedure destroy_ccsd_intermediates_talsh
     end interface destroy_ccsd_intermediates

     interface get_eta_ov
        module procedure get_eta_ov_talsh
     end interface get_eta_ov
  
     interface get_eta_oovv
        module procedure get_eta_oovv_talsh
     end interface get_eta_oovv


!    top-level driver for forming intermediates with CCSD wavefunctions
     public form_ccsd_intermediates
     public destroy_ccsd_intermediates
     public get_eta_ov, get_eta_oovv

     contains

       subroutine get_tau_talsh(tau_tensor,t2_tensor,t1_tensor)

        use talsh

        type(talsh_tens_t), intent(inout) :: tau_tensor,t2_tensor
        type(talsh_tens_t), intent(inout), optional :: t1_tensor
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

!       tau ! validated
        ierr=talsh_tensor_init(tau_tensor)
        ierr=talsh_tensor_add("T(a,b,i,j)+=R(a,b,i,j)",tau_tensor,t2_tensor)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",tau_tensor,t1_tensor,t1_tensor)                  ! o(2)v(2)
            ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",tau_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)  ! o(2)v(2)
        end if

       end subroutine get_tau_talsh

       subroutine get_hoo_talsh(comm_t,oovv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       oo ! validated
        ierr=talsh_tensor_init(comm_t%hoo)
        ierr=talsh_tensor_add("H(k,i)+=F(k,i)",comm_t%hoo,comm_t%foo,scale=MINUS_ONE)
        ierr=talsh_tensor_contract("H(k,i)+=V(k,l,c,d)*T(c,d,l,i)",comm_t%hoo,oovv_tensor,comm_t%tau,scale=ONE_HALF) ! o(3)v(2)

       end subroutine get_hoo_talsh

       subroutine get_hov_talsh_dt(comm_t,oovv_tensor,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_comm_tens_t), intent(inout) :: comm_t
        type(talsh_tens_t),    intent(inout) :: oovv_tensor
        logical,               intent(in   ) :: CCD
        integer :: ierr

!       ov ! validated
        ierr=talsh_tensor_init(comm_t%hov)
        if (.not.CCD) ierr=talsh_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,oovv_tensor,comm_t%t1) ! o(2)v(2)
        ierr=talsh_tensor_add("H(k,c)+=F(k,c)",comm_t%hov,comm_t%fov)

       end subroutine get_hov_talsh_dt

       subroutine get_hov_talsh_no_dt(hov_tensor,oovv_tensor,fov_tensor,t1_tensor)

        use talsh

        type(talsh_tens_t),intent(inout) :: hov_tensor
        type(talsh_tens_t), intent(inout) :: oovv_tensor, fov_tensor
        type(talsh_tens_t), intent(inout), optional :: t1_tensor
        logical :: CCD=.false.
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

!       ov ! validated
        ierr=talsh_tensor_init(hov_tensor)
        if (.not.CCD) ierr=talsh_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",hov_tensor,oovv_tensor,t1_tensor) ! o(2)v(2)
        ierr=talsh_tensor_add("H(k,c)+=F(k,c)",hov_tensor,fov_tensor)

       end subroutine get_hov_talsh_no_dt

       subroutine get_hvv_talsh(comm_t,oovv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       vv ! validated
        ierr=talsh_tensor_init(comm_t%hvv)
        ierr=talsh_tensor_add("H(a,c)+=F(a,c)",comm_t%hvv,comm_t%fvv)
        ierr=talsh_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%hvv,oovv_tensor,comm_t%tau,scale=ONE_HALF) ! o(2)v(3)

       end subroutine get_hvv_talsh

       subroutine get_goo_talsh(comm_t,ooov_tensor,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: ooov_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       oo ! validated
        ierr=talsh_tensor_init(comm_t%goo)
        ierr=talsh_tensor_add("G(k,i)+=H(k,i)",comm_t%goo,comm_t%hoo)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,ooov_tensor,comm_t%t1,scale=MINUS_ONE)     ! o(3)v(1)
            ierr=talsh_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1)                         ! o(2)v(1)
        end if

       end subroutine get_goo_talsh

       subroutine get_gvv_talsh(comm_t,vovv_tensor,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: vovv_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       vv ! validated
        ierr=talsh_tensor_init(comm_t%gvv)
        ierr=talsh_tensor_add("G(a,c)+=H(a,c)",comm_t%gvv,comm_t%hvv)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,vovv_tensor,comm_t%t1)             ! o(1)v(3)
            ierr=talsh_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,scale=MINUS_ONE) ! o(1)v(2)
        end if

       end subroutine get_gvv_talsh

       subroutine get_a_talsh(comm_t,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t, talsh_comm_tens_t

        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

!       oooo ! validated
        ierr=talsh_tensor_init(comm_t%a_int)
        ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",comm_t%a_int,int_t%oooo)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
            ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        end if
        ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%a_int,int_t%oovv,comm_t%tau, &
                                    scale=ONE_HALF)

       end subroutine get_a_talsh

       subroutine get_h_talsh(comm_t,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t, talsh_comm_tens_t

        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

!       vovo ! more tests
        ierr=talsh_tensor_init(comm_t%h_int)
        if (.not.CCD) then
            ierr=talsh_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1)                 ! o(2)v(3)
            ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
        end if
        ierr=talsh_tensor_add("H(a,k,c,i)+=V(a,k,c,i)",comm_t%h_int,int_t%vovo)

       end subroutine get_h_talsh

       subroutine get_fbar_oo_talsh_dt(fbar_oo,comm_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_tens_t),   intent(inout)   :: fbar_oo
        type(talsh_comm_tens_t),intent(inout)   :: comm_t
        logical, intent(in) :: CCD
        integer :: ierr

        ierr=talsh_tensor_add("F(i,m)+=G(i,m)",fbar_oo,comm_t%goo,scale=MINUS_ONE)
        if (.not.CCD) ierr=talsh_tensor_contract("F(i,m)+=R(i,e)*T(e,m)",fbar_oo,comm_t%fov,comm_t%t1,scale=TWO)

       end subroutine get_fbar_oo_talsh_dt

       subroutine get_fbar_oo_talsh_no_dt(fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: fbar_oo
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediate
        type(talsh_tens_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=talsh_tensor_init(fbar_oo)
        ierr=talsh_tensor_add("F(i,m)+=U(i,m)",fbar_oo,foo_tensor)
        if (.not.CCD) then
          ierr=talsh_tensor_contract("F(i,m)+=U(i,e)*T(e,m)",fbar_oo,fov_tensor,t1_tensor)
          ierr=talsh_tensor_contract("F(i,m)+=V(i,n,m,e)*T(e,n)",fbar_oo,int_t%ooov,t1_tensor)
        end if

        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=talsh_tensor_contract("F(i,m)+=V(i,n,e,f)*T(e,f,m,n)",fbar_oo,int_t%oovv,tau_tensor,ONE_HALF)

        ierr=talsh_tensor_destruct(tau_tensor)

       end subroutine get_fbar_oo_talsh_no_dt

       subroutine get_fbar_ov_talsh(fbar_ov,hov_tensor)

        use talsh

        type(talsh_tens_t),intent(inout)   :: fbar_ov
        type(talsh_tens_t),intent(inout)   :: hov_tensor
        integer :: ierr

        ierr=talsh_tensor_add("F(m,e)+=H(m,e)",fbar_ov,hov_tensor)

       end subroutine get_fbar_ov_talsh

       subroutine get_fbar_ov_talsh_scratch(fbar_ov,fov_tensor,CCD,oovv_tensor,t1_tensor)

        use talsh

!       input variables
        type(talsh_tens_t), intent(inout) :: fbar_ov
        type(talsh_tens_t), intent(inout) :: fov_tensor
        type(talsh_tens_t), intent(inout), optional :: oovv_tensor, t1_tensor
        logical,            intent(in   ) :: CCD
!       error code
        integer :: ierr

        ierr=talsh_tensor_init(fbar_ov)
        ierr=talsh_tensor_add("F(m,e)+=U(m,e)",fbar_ov,fov_tensor)
        if (.not.CCD) ierr=talsh_tensor_contract("F(m,e)+=V(m,n,e,f)*T(f,n)",fbar_ov,oovv_tensor,t1_tensor)

       end subroutine get_fbar_ov_talsh_scratch

       subroutine get_fbar_vv_talsh(fbar_vv,gvv_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_comm_tens_t

        type(talsh_tens_t),intent(inout)   :: fbar_vv
        type(talsh_tens_t),intent(inout)   :: gvv_tensor
        integer :: ierr

        ierr=talsh_tensor_add("F(e,a)+=G(e,a)",fbar_vv,gvv_tensor)

       end subroutine get_fbar_vv_talsh

       subroutine get_fbar_vv_talsh_scratch(fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: fbar_vv
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: fvv_tensor, fov_tensor, t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediate
        type(talsh_tens_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=talsh_tensor_init(fbar_vv)
        ierr=talsh_tensor_add("F(e,a)+=U(e,a)",fbar_vv,fvv_tensor)
        if (.not.CCD) then
          ierr=talsh_tensor_contract("F(e,a)+=U(m,a)*T(e,m)",fbar_vv,fov_tensor,t1_tensor,scale=MINUS_ONE)
          ierr=talsh_tensor_contract("F(e,a)+=V(e,m,a,f)*T(f,m)",fbar_vv,int_t%vovv,t1_tensor)
        end if

        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=talsh_tensor_contract("F(e,a)+=V(m,n,a,f)*T(f,e,m,n)",fbar_vv,int_t%oovv,tau_tensor,ONE_HALF)

        ierr=talsh_tensor_destruct(tau_tensor)

       end subroutine get_fbar_vv_talsh_scratch

       subroutine get_w_oooo_talsh(w_oooo,a_int)

        use talsh

        type(talsh_tens_t),intent(inout) :: w_oooo
        type(talsh_tens_t),intent(inout) :: a_int
        integer :: ierr

        ierr=talsh_tensor_add("W(i,j,m,n)+=A(i,j,m,n)",w_oooo,a_int)

       end subroutine get_w_oooo_talsh

       subroutine get_w_oooo_talsh_scratch(w_oooo,int_t,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: w_oooo
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediate
        type(talsh_tens_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=talsh_tensor_add("W(i,j,m,n)+=V(i,j,m,n)",w_oooo,int_t%oooo)
        if (.not.CCD) then
          ierr=talsh_tensor_contract("W(i,j,m,n)+=V(j,i,n,e)*T(e,m)",w_oooo,int_t%ooov,t1_tensor)
          ierr=talsh_tensor_contract("W(i,j,m,n)+=V(i,j,m,e)*T(e,n)",w_oooo,int_t%ooov,t1_tensor)
        end if

        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=talsh_tensor_contract("W(i,j,m,n)+=V(i,j,e,f)*T(e,f,m,n)",w_oooo,int_t%oovv,tau_tensor,ONE_HALF)

        ierr=talsh_tensor_destruct(tau_tensor)

       end subroutine get_w_oooo_talsh_scratch

       subroutine get_w_ooov_talsh(w_ooov,int_t,t1_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: w_ooov
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=talsh_tensor_add("W(m,n,i,e)+=V(m,n,i,e)",w_ooov,int_t%ooov)
        if (.not.CCD) ierr=talsh_tensor_contract("W(m,n,i,e)+=T(f,i)*V(m,n,f,e)",w_ooov,t1_tensor,int_t%oovv)
       
       end subroutine get_w_ooov_talsh

       subroutine get_w_ovoo_talsh(lambda_t,comm_t,int_t,CCD,one_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t, talsh_comm_tens_t, talsh_lambda_tens_t

        type(talsh_lambda_tens_t),intent(inout) :: lambda_t
        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(talsh_tens_t) :: one_tensor
        integer :: ierr

        ierr=talsh_tensor_contract("W(i,e,m,n)+=V+(m,n,i,e)*K()",lambda_t%w_ovoo,int_t%ooov,one_tensor)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=F(i,f)*T(f,e,m,n)",lambda_t%w_ovoo,lambda_t%fbar_ov,comm_t%t2)
!       Need to permute m and n, so term appears twice
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(i,o,m,f)*T(e,f,n,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(o,i,n,f)*T(e,f,m,o)",lambda_t%w_ovoo,int_t%ooov,comm_t%t2)
        if (.not.CCD) then 
!        Need to permute m and n, so term appears twice
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,n)*T(f,m)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,m)*T(f,n)",lambda_t%w_ovoo,lambda_t%wbar_vovo,comm_t%t1)
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(o,i,m,n)*T(e,o)",lambda_t%w_ovoo,lambda_t%w_oooo,comm_t%t1)
        end if
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(e,i,f,g)*T(g,f,m,n)",lambda_t%w_ovoo,int_t%vovv,comm_t%tau,scale=ONE_HALF)

       end subroutine get_w_ovoo_talsh

       subroutine get_w_ovoo_talsh_scratch(w_ovoo,int_t,fbar_ov,w_oooo,wbar_vovo,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t

        type(talsh_tens_t),    intent(inout) :: w_ovoo
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: fbar_ov, w_oooo, wbar_vovo, t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(talsh_tens_t) :: tau_tensor, one_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4), one_dims(1)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V+(m,n,i,e)*K()",w_ovoo,int_t%ooov,one_tensor)
        ierr=talsh_tensor_destruct(one_tensor)

        ierr=talsh_tensor_contract("W(i,e,m,n)+=F(i,f)*T(f,e,m,n)",w_ovoo,fbar_ov,t2_tensor)
!       Need to permute m and n, so term appears twice
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(i,o,m,f)*T(e,f,n,o)",w_ovoo,int_t%ooov,t2_tensor)
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(o,i,n,f)*T(e,f,m,o)",w_ovoo,int_t%ooov,t2_tensor)
        if (.not.CCD) then 
!        Need to permute m and n, so term appears twice
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,n)*T(f,m)",w_ovoo,wbar_vovo,t1_tensor,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(e,i,f,m)*T(f,n)",w_ovoo,wbar_vovo,t1_tensor)
         ierr=talsh_tensor_contract("W(i,e,m,n)+=R(o,i,m,n)*T(e,o)",w_ovoo,w_oooo,t1_tensor)
        end if

        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if
        ierr=talsh_tensor_contract("W(i,e,m,n)+=V(e,i,f,g)*T(g,f,m,n)",w_ovoo,int_t%vovv,tau_tensor,scale=ONE_HALF)
        ierr=talsh_tensor_destruct(tau_tensor)

       end subroutine get_w_ovoo_talsh_scratch

       subroutine get_w_vovo_talsh(w_vovo,comm_t,oovv_tensor,nocc,nvir,CCD)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_comm_tens_t
        use talsh_common_routines, only       : tensor_norm2

        type(talsh_tens_t),intent(inout) :: w_vovo
        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_tens_t),intent(inout) :: oovv_tensor
        integer,intent(in) :: nocc, nvir
        logical, intent(in) :: CCD
        type(talsh_tens_t) :: tau2_tensor
        integer(INTD) :: vvoo_dims(4)
        integer(INTD) :: ierr

!       Auxiliary tensor
        vvoo_dims(1:2) = nvir
        vvoo_dims(3:4) = nocc
        ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)

        print*, "debug: w_vovo tensor norm before contractions (get_w_vovo_talsh)", tensor_norm2(w_vovo)
        ierr=talsh_tensor_add("W(b,m,e,j)+=H(b,m,e,j)",w_vovo,comm_t%h_int)
        ierr=talsh_tensor_init(tau2_tensor)
        ierr=talsh_tensor_add("T(f,b,j,n)+=R(f,b,j,n)",tau2_tensor,comm_t%t2)
        if (.not.CCD) ierr=talsh_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",tau2_tensor,comm_t%t1,comm_t%t1)
        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,tau2_tensor)

        ierr=talsh_tensor_destruct(tau2_tensor)
        print *, "debug: w_vovo tensor norm after contractions (get_w_vovo_talsh)", tensor_norm2(w_vovo)

       end subroutine get_w_vovo_talsh

       subroutine get_w_vovo_talsh_scratch(w_vovo,int_t,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t
        use talsh_common_routines, only       : tensor_norm2

!       input variables
        type(talsh_tens_t),    intent(inout) :: w_vovo
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(talsh_tens_t) :: tau2_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        print*, "debug: w_vovo tensor norm before contractions (get_w_vovo_talsh_scratch)", tensor_norm2(w_vovo)
        ierr=talsh_tensor_add("W(b,m,e,j)+=V(b,m,e,j)",w_vovo,int_t%vovo)
        if (.not.CCD) then
          ierr=talsh_tensor_contract("W(b,m,e,j)+=V(b,m,e,f)*T(f,j)",w_vovo,int_t%vovv,t1_tensor)
          ierr=talsh_tensor_contract("W(b,m,e,j)+=V(n,m,j,e)*T(b,n)",w_vovo,int_t%ooov,t1_tensor)
        end if
 
        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
        ierr=talsh_tensor_add("T(f,b,j,n)+=R(f,b,j,n)",tau2_tensor,t2_tensor)
        print*, "debug: t1 tensor norm (get_w_vovo_talsh_scratch)", tensor_norm2(t1_tensor)
        print*, "debug: t2 tensor norm (get_w_vovo_talsh_scratch)", tensor_norm2(t2_tensor)
        if (.not.CCD) ierr=talsh_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",tau2_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)
        print*, "debug: tau2 tensor norm (get_w_vovo_talsh_scratch)", tensor_norm2(tau2_tensor)
        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,int_t%oovv,tau2_tensor)

        ierr=talsh_tensor_destruct(tau2_tensor)
        print*, "debug: w_vovo tensor norm before contractions (get_w_vovo_talsh_scratch)", tensor_norm2(w_vovo)

       end subroutine get_w_vovo_talsh_scratch

       subroutine get_wbar_vovo_talsh(wbar_vovo,t2_tensor,int_t)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t

        type(talsh_tens_t),intent(inout) :: wbar_vovo
        type(talsh_tens_t),intent(inout) :: t2_tensor
        type(talsh_intg_tens_t),intent(inout) :: int_t
        integer :: ierr

        ierr=talsh_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(b,f,n,j)",wbar_vovo,int_t%oovv,t2_tensor)
        ierr=talsh_tensor_add("W(b,m,e,j)+=V(b,m,e,j)",wbar_vovo,int_t%vovo)

       end subroutine get_wbar_vovo_talsh

       subroutine get_w_vvvo_talsh_dt(lambda_t,comm_t,int_t,CCD,one_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t, talsh_comm_tens_t, talsh_lambda_tens_t

        type(talsh_lambda_tens_t),intent(inout) :: lambda_t
        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(talsh_tens_t) :: one_tensor
        integer :: ierr

        ierr=talsh_tensor_contract("W(e,f,a,m)+=V+(a,m,e,f)*K()",lambda_t%w_vvvo,int_t%vovv,one_tensor)
!       Need to permute e and f, so term appears twice
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(e,n,a,g)*T(f,g,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(f,n,a,g)*T(g,e,m,n)",lambda_t%w_vvvo,int_t%vovv,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=F(n,a)*T(e,f,m,n)",lambda_t%w_vvvo,lambda_t%fbar_ov,comm_t%t2)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(n,o,m,a)*T(e,f,o,n)",lambda_t%w_vvvo,int_t%ooov,comm_t%tau,scale=ONE_HALF)
        if (.not.CCD) then
!        Need to permute e and f, so term appears twice
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(f,n,a,m)*T(e,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1)
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,n,a,m)*T(f,n)",lambda_t%w_vvvo,lambda_t%wbar_vovo,comm_t%t1,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,f,a,g)*T(g,m)",lambda_t%w_vvvo,lambda_t%w_vvvv,comm_t%t1)
        end if

       end subroutine get_w_vvvo_talsh_dt

       subroutine get_w_vvvo_talsh_no_dt(w_vvvo,int_t,w_vvvv,wbar_vovo,fbar_ov,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: w_vvvo
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: w_vvvv, wbar_vovo, fbar_ov, t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(talsh_tens_t) :: tau_tensor, one_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4), one_dims(1)
!       local variables
        logical :: CCD = .false.
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V+(a,m,e,f)*K()",w_vvvo,int_t%vovv,one_tensor)
        ierr=talsh_tensor_destruct(one_tensor)

!       Need to permute e and f, so term appears twice
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(e,n,a,g)*T(f,g,m,n)",w_vvvo,int_t%vovv,t2_tensor)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(f,n,a,g)*T(g,e,m,n)",w_vvvo,int_t%vovv,t2_tensor)
        ierr=talsh_tensor_contract("W(e,f,a,m)+=F(n,a)*T(e,f,m,n)",w_vvvo,fbar_ov,t2_tensor)

        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if
        ierr=talsh_tensor_contract("W(e,f,a,m)+=V(n,o,m,a)*T(e,f,o,n)",w_vvvo,int_t%ooov,tau_tensor,scale=ONE_HALF)
        ierr=talsh_tensor_destruct(tau_tensor)

        if (.not.CCD) then
!        Need to permute e and f, so term appears twice
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(f,n,a,m)*T(e,n)",w_vvvo,wbar_vovo,t1_tensor)
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,n,a,m)*T(f,n)",w_vvvo,wbar_vovo,t1_tensor,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("W(e,f,a,m)+=R(e,f,a,g)*T(g,m)",w_vvvo,w_vvvv,t1_tensor)
        end if

       end subroutine get_w_vvvo_talsh_no_dt

       subroutine get_w_vovv_talsh(w_vovv,int_t,t1_tensor)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: w_vovv
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        if (.not.CCD) ierr=talsh_tensor_contract("W(a,m,e,f)+=V(m,n,e,f)*T(a,n)",w_vovv,int_t%oovv,t1_tensor)
        ierr=talsh_tensor_add("W(a,m,e,f)+=V(a,m,e,f)",w_vovv,int_t%vovv)

       end subroutine get_w_vovv_talsh

       subroutine get_w_vvvv_talsh_dt(w_vvvv,comm_t,int_t,CCD)

        use talsh
        use exacorr_datatypes, only : talsh_intg_tens_t, talsh_comm_tens_t

        type(talsh_tens_t),intent(inout) :: w_vvvv
        type(talsh_comm_tens_t),intent(inout) :: comm_t
        type(talsh_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

        if (.not.CCD) then
!        Permute e and f, so term appears twice
         ierr=talsh_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv,int_t%vovv,comm_t%t1)
         ierr=talsh_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv,int_t%vovv,comm_t%t1,scale=MINUS_ONE)
        end if
        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv,int_t%oovv,comm_t%tau,scale=ONE_HALF)
        ierr=talsh_tensor_add("W(e,f,a,b)+=V(e,f,a,b)",w_vvvv,int_t%vvvv)

       end subroutine get_w_vvvv_talsh_dt

       subroutine get_w_vvvv_talsh_no_dt(w_vvvv,int_t,t2_tensor,t1_tensor)

        use talsh
        use tensor_algebra
        use exacorr_datatypes, only : talsh_intg_tens_t

!       input variables
        type(talsh_tens_t),    intent(inout) :: w_vvvv
        type(talsh_intg_tens_t), intent(inout) :: int_t
        type(talsh_tens_t),    intent(inout) :: t2_tensor
        type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(talsh_tens_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       local variables
        logical :: CCD = .false.
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        if (.not.CCD) then
!        Permute e and f, so term appears twice
         ierr=talsh_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv,int_t%vovv,t1_tensor)
         ierr=talsh_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv,int_t%vovv,t1_tensor,scale=MINUS_ONE)
        end if

        ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr=talsh_tensor_construct(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=talsh_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv,int_t%oovv,tau_tensor,scale=ONE_HALF)
        ierr=talsh_tensor_add("W(e,f,a,b)+=V(e,f,a,b)",w_vvvv,int_t%vvvv)

        ierr=talsh_tensor_destruct(tau_tensor)

       end subroutine get_w_vvvv_talsh_no_dt

       subroutine get_goo_lambda_talsh(goo_tensor,l2_tensor,t2_tensor)

        use talsh

        type(talsh_tens_t), intent(inout) :: goo_tensor
        type(talsh_tens_t), intent(inout) :: l2_tensor, t2_tensor
        integer :: ierr

        ierr=talsh_tensor_init(goo_tensor)
        ierr=talsh_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,t2_tensor,scale=ONE_HALF)

       end subroutine get_goo_lambda_talsh

       subroutine get_gvv_lambda_talsh(gvv_tensor,l2_tensor,t2_tensor)

        use talsh

        type(talsh_tens_t), intent(inout) :: gvv_tensor
        type(talsh_tens_t), intent(inout) :: l2_tensor, t2_tensor
        integer :: ierr

        ierr=talsh_tensor_init(gvv_tensor)
        ierr=talsh_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,n,m)",gvv_tensor,l2_tensor,t2_tensor,scale=ONE_HALF)

       end subroutine get_gvv_lambda_talsh

       subroutine form_ccsd_intermediates_talsh(interm_t, CCD, nocc, nvir, comm_t, int_t, l2_tensor)
          use tensor_algebra
          use talsh
          use exacorr_datatypes
          use exacorr_utils
          use talsh_common_routines

! input variables
          type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
! CCD switch
          logical, intent(in) :: CCD
          integer, intent(in) :: nocc, nvir
          type(talsh_comm_tens_t),    intent(inout) :: comm_t
          type(talsh_intg_tens_t), intent(inout) :: int_t
          type(talsh_tens_t),    intent(inout), optional :: l2_tensor 
! tensor dims
         integer, dimension(2) :: oo_dims, ov_dims, vv_dims, vo_dims
         integer, dimension(4) :: oooo_dims, ooov_dims, ovoo_dims, vovo_dims, vvoo_dims, oovv_dims
         integer, dimension(4) :: vvvv_dims, vovv_dims, vvvo_dims
!     error code
         integer :: ierr

! set tensor dims
         oo_dims = (/nocc,nocc/)
         ov_dims = (/nocc,nvir/)
         vv_dims = (/nvir,nvir/)
         vo_dims = (/nvir,nocc/)

         oooo_dims = (/nocc,nocc,nocc,nocc/)
         ooov_dims = (/nocc,nocc,nocc,nvir/)
         ovoo_dims = (/nocc,nvir,nocc,nocc/)
         vovo_dims = (/nvir,nocc,nvir,nocc/)
         vvoo_dims = (/nvir,nvir,nocc,nocc/)

         vvvv_dims = (/nvir,nvir,nvir,nvir/)
         vvvo_dims = (/nvir,nvir,nvir,nocc/)
         vovv_dims = (/nvir,nocc,nvir,nvir/)

         interm_t%nocc = nocc
         interm_t%nvir = nvir
         interm_t%ccd  = ccd

         print*, "foo_tensor norm ", tensor_norm2(comm_t%foo)
         print*, "fov_tensor norm ", tensor_norm2(comm_t%fov)
         print*, "fvv_tensor norm ", tensor_norm2(comm_t%fvv)

!  copy t1 and t2 from comm_t
         ierr=talsh_tensor_construct(interm_t%t1,C8,vo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%t2,C8,vvoo_dims,init_val=ZERO)

         ierr=talsh_tensor_add("IntmT(a,i)+=T(a,i)",interm_t%t1,comm_t%t1)
         ierr=talsh_tensor_add("IntmT(a,b,i,j)+=T(a,b,i,j)",interm_t%t2,comm_t%t2)
         print*, "copying t2 tensor, norm ", tensor_norm2(comm_t%t2),tensor_norm2(interm_t%t2)

!  create all fixed intermediates
         ierr=talsh_tensor_construct(interm_t%fbar_oo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%fbar_ov,C8,ov_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%fbar_vv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%w_oooo,C8,oooo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%w_ooov,C8,ooov_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%w_vovo,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

         ierr=talsh_tensor_construct(interm_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(interm_t%w_vovv,C8,vovv_dims,init_val=ZERO)

         call get_wbar_vovo(interm_t%wbar_vovo,interm_t%t2,int_t)
         if (.not.CCD) then
            call get_fbar_oo(interm_t%fbar_oo,int_t,comm_t%foo,comm_t%fov,interm_t%t2,interm_t%t1)
            call get_fbar_ov(interm_t%fbar_ov,comm_t%fov,CCD,int_t%oovv,interm_t%t1)
            call get_fbar_vv(interm_t%fbar_vv,int_t,comm_t%fvv,comm_t%fov,interm_t%t2,interm_t%t1)
            call get_w_oooo (interm_t%w_oooo,int_t,interm_t%t2,interm_t%t1)
            call get_w_ooov (interm_t%w_ooov,int_t,interm_t%t1)
            call get_w_ovoo (interm_t%w_ovoo,int_t,interm_t%fbar_ov,interm_t%w_oooo,interm_t%wbar_vovo,interm_t%t2,interm_t%t1)
            call get_w_vovo (interm_t%w_vovo,int_t,interm_t%t2,interm_t%t1)
            call get_w_vvvv (interm_t%w_vvvv,int_t,interm_t%t2,interm_t%t1)
            call get_w_vvvo (interm_t%w_vvvo,int_t,interm_t%w_vvvv,interm_t%wbar_vovo,interm_t%fbar_ov,interm_t%t2,interm_t%t1)
            call get_w_vovv (interm_t%w_vovv,int_t,interm_t%t1)
         else
            call get_fbar_oo(interm_t%fbar_oo,int_t,comm_t%foo,comm_t%fov,interm_t%t2)
            call get_fbar_ov(interm_t%fbar_ov,comm_t%fov,CCD,int_t%oovv)
            call get_fbar_vv(interm_t%fbar_vv,int_t,comm_t%fvv,comm_t%fov,interm_t%t2)
            call get_w_oooo (interm_t%w_oooo,int_t,interm_t%t2)
            call get_w_ooov (interm_t%w_ooov,int_t)
            call get_w_ovoo (interm_t%w_ovoo,int_t,interm_t%fbar_ov,interm_t%w_oooo,interm_t%wbar_vovo,interm_t%t2)
            call get_w_vovo (interm_t%w_vovo,int_t,interm_t%t2)
         end if


!  if information is present for the lagrange multipliers, construct the intermediates
!  required for lhs/lambda equations
         if (present(l2_tensor)) then
            oovv_dims = (/nocc,nocc,nvir,nvir/)
            ierr=talsh_tensor_construct(interm_t%goo,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(interm_t%gvv,C8,vv_dims,init_val=ZERO)
            call get_goo_lambda(interm_t%goo,l2_tensor,interm_t%t2)
            call get_gvv_lambda(interm_t%gvv,l2_tensor,interm_t%t2)
            print*, "g_oo      tensor norm", tensor_norm2(interm_t%goo)
            print*, "g_vv      tensor norm", tensor_norm2(interm_t%gvv)
         end if

         print*, "t1        tensor norm", tensor_norm2(interm_t%t1)
         print*, "t2        tensor norm", tensor_norm2(interm_t%t2)
         print*, "fbar_oo   tensor norm", tensor_norm2(interm_t%fbar_oo)
         print*, "fbar_ov   tensor norm", tensor_norm2(interm_t%fbar_ov)
         print*, "fbar_vv   tensor norm", tensor_norm2(interm_t%fbar_vv)
         print*, "w_oooo    tensor norm", tensor_norm2(interm_t%w_oooo)
         print*, "w_ooov    tensor norm", tensor_norm2(interm_t%w_ooov)
         print*, "w_ovoo    tensor norm", tensor_norm2(interm_t%w_ovoo)
         print*, "w_vovo    tensor norm", tensor_norm2(interm_t%w_vovo)
         print*, "wbar_vovo tensor norm", tensor_norm2(interm_t%wbar_vovo)
         print*, "w_vvvv    tensor norm", tensor_norm2(interm_t%w_vvvv)
         print*, "w_vvvo    tensor norm", tensor_norm2(interm_t%w_vvvo)
         print*, "w_vovv    tensor norm", tensor_norm2(interm_t%w_vovv)
         

         ierr=talsh_tensor_destruct(interm_t%wbar_vovo)

      end subroutine

      subroutine destroy_ccsd_intermediates_talsh(interm_t)
         use tensor_algebra
         use talsh
         use exacorr_datatypes
         use exacorr_utils
         use talsh_common_routines

! input variables
         type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
         integer :: ierr

         print*, "destroying ccsd intermediates for talsh"
         ierr=talsh_tensor_destruct(interm_t%goo)
         ierr=talsh_tensor_destruct(interm_t%gvv)
         ierr=talsh_tensor_destruct(interm_t%fbar_oo)
         ierr=talsh_tensor_destruct(interm_t%fbar_ov)
         ierr=talsh_tensor_destruct(interm_t%fbar_vv)
         ierr=talsh_tensor_destruct(interm_t%w_oooo)
         ierr=talsh_tensor_destruct(interm_t%w_ooov)
         ierr=talsh_tensor_destruct(interm_t%w_ovoo)
         ierr=talsh_tensor_destruct(interm_t%w_vovo)
         ierr=talsh_tensor_destruct(interm_t%w_vovv)
         ierr=talsh_tensor_destruct(interm_t%w_vvvo)
         ierr=talsh_tensor_destruct(interm_t%w_vvvv)
         print*, "done destroying ccsd intermediates"

      end subroutine


      subroutine get_eta_ov_talsh(eta_ov,fbar_ov,ov_dims)
        use talsh
        use tensor_algebra

        type(talsh_tens_t),intent(inout)   :: eta_ov
        type(talsh_tens_t),intent(inout)   :: fbar_ov
        integer :: ierr
        integer(INTD), dimension(2)   :: ov_dims

        ierr=talsh_tensor_construct(eta_ov,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_add("N(i,a)+=F(i,a)",eta_ov,fbar_ov)
      end subroutine


      subroutine get_eta_oovv_talsh(eta_oovv,int_oovv,oovv_dims)
        use talsh
        use tensor_algebra

        type(talsh_tens_t),intent(inout)   :: eta_oovv
        type(talsh_tens_t),intent(inout)   :: int_oovv
        integer :: ierr
        integer(INTD), dimension(4)   :: oovv_dims

        ierr=talsh_tensor_construct(eta_oovv,C8,oovv_dims,init_val=ZERO)
        ierr=talsh_tensor_add("N(i,j,a,b)+=V(i,j,a,b)",eta_oovv,int_oovv)
      end subroutine get_eta_oovv_talsh

end module
