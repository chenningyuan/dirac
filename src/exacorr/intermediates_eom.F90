module intermediates_eom_ccsd

!Module to make intermediates needed in CC and Lambda equations

     use exacorr_tensor_interface
     use tensor_algebra
     use exacorr_datatypes
     use exacorr_utils, only: get_free_fileunit
     use eti_common_routines
     use intermediates, only: get_wbar_vovo, get_goo_lambda, get_gvv_lambda

     implicit none
     private

     complex(8), parameter :: ZERO=(0.D0,0.D0), MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0), &
                              ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                              MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0), MINUS_ONE_HALF=(-0.5D0,0.D0)
    real(8), parameter     :: ONE_QUARTER=0.25D0

     interface get_tau
        module procedure get_tau_eti_no_comm_t
     end interface get_tau

     interface get_fbar_oo
        module procedure get_fbar_oo_eti_no_dt
     end interface get_fbar_oo

     interface get_fbar_ov
        module procedure get_fbar_ov_eti_scratch
     end interface get_fbar_ov

     interface get_fbar_vv
        module procedure get_fbar_vv_eti_scratch
     end interface get_fbar_vv

     interface get_w_oooo
        module procedure get_w_oooo_eti_scratch
     end interface get_w_oooo

     interface get_w_ooov
        module procedure get_w_ooov_talsh
     end interface get_w_ooov

     interface get_w_ovoo
        module procedure get_w_ovoo_eti_scratch
     end interface get_w_ovoo

     interface get_w_vovo
        module procedure get_w_vovo_eti_scratch
     end interface get_w_vovo

     interface get_w_vvvo
        module procedure get_w_vvvo_eti_no_dt
     end interface get_w_vvvo

     interface get_w_vovv
        module procedure get_w_vovv_eti_no_logical
     end interface get_w_vovv

     interface get_w_vvvv
        module procedure get_w_vvvv_eti_no_dt
     end interface get_w_vvvv

     interface form_ccsd_intermediates
        module procedure form_ccsd_intermediates_talsh
     end interface form_ccsd_intermediates

     interface destroy_ccsd_intermediates
        module procedure destroy_ccsd_intermediates_talsh
     end interface destroy_ccsd_intermediates

     interface get_eta_ov
        module procedure get_eta_ov_eti
     end interface get_eta_ov
  
     interface get_eta_oovv
        module procedure get_eta_oovv_talsh
     end interface get_eta_oovv


!    top-level driver for forming intermediates with CCSD wavefunctions
     public form_ccsd_intermediates
     public destroy_ccsd_intermediates
     public get_eta_ov, get_eta_oovv

     contains

       subroutine get_tau_eti_no_comm_t(tau_tensor,t2_tensor,t1_tensor)
        type(tensor_t), intent(inout) :: tau_tensor,t2_tensor
        type(tensor_t), intent(inout), optional :: t1_tensor
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

!       tau ! validated
        ierr=initialize_tensor_t(tau_tensor)
        ierr=t_einsum("T(a,b,i,j)+=R(a,b,i,j)",tau_tensor,t2_tensor)
        if (.not.CCD) then
            ierr=t_einsum("T(a,b,i,j)+=T(a,i)*T(b,j)",tau_tensor,t1_tensor,t1_tensor)                  ! o(2)v(2)
            ierr=t_einsum("T(a,b,i,j)+=T(a,j)*T(b,i)",tau_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)  ! o(2)v(2)
        end if

       end subroutine get_tau_eti_no_comm_t

       subroutine get_hov_eti_no_comm_t(hov_tensor,oovv_tensor,fov_tensor,t1_tensor)
        type(tensor_t),intent(inout) :: hov_tensor
        type(tensor_t), intent(inout) :: oovv_tensor, fov_tensor
        type(tensor_t), intent(inout), optional :: t1_tensor
        logical :: CCD=.false.
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

!       ov ! validated
        ierr=initialize_tensor_t(hov_tensor)
        if (.not.CCD) ierr=t_einsum("H(k,c)+=V(k,l,c,d)*T(d,l)",hov_tensor,oovv_tensor,t1_tensor) ! o(2)v(2)
        ierr=t_einsum("H(k,c)+=F(k,c)",hov_tensor,fov_tensor)

       end subroutine get_hov_eti_no_comm_t

! this has no analog in intermediates.F90 and not only not uses comm_t but also codes a different 
! intermediate Fbar_oo (which follows Shee et al. JCP 148, 174113, 2018 but has been revised/corrected if needed)
       subroutine get_fbar_oo_eti_no_dt(fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: fbar_oo
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: foo_tensor, fov_tensor, t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediate
        type(tensor_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        integer                :: f_oo_rel
        real(8),allocatable,dimension(:) :: Fbar_mi
        complex(8), pointer    ::fbaroo_tens(:,:), fbaroo_tens1(:,:)
        integer(INTD),dimension(2)   :: oo_dims
        integer :: i,j,n,m
        type(C_PTR) :: body_p
        complex(8) :: q

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=initialize_tensor_t(fbar_oo)
        ierr=t_einsum("F(i,m)+=U(i,m)",fbar_oo,foo_tensor)
        if (.not.CCD) then
          ierr=t_einsum("F(i,m)+=U(i,e)*T(e,m)",fbar_oo,fov_tensor,t1_tensor)
          ierr=t_einsum("F(i,m)+=V(i,n,m,e)*T(e,n)",fbar_oo,int_t%ooov,t1_tensor)
        end if

!       ierr=talsh_tensor_dimensions(t2_tensor,tens_rank,vvoo_dims)
        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=t_einsum("F(i,m)+=V(i,n,e,f)*T(e,f,m,n)",fbar_oo,int_t%oovv,tau_tensor,scale=ONE_HALF)
        ierr=destroy_tensor_t(tau_tensor)

       end subroutine get_fbar_oo_eti_no_dt

! This is not the same as get_fbar_ov_talsh in intermediates.F90, as the second contraction (oovv*t1) is only present here.
       subroutine get_fbar_ov_eti_scratch(fbar_ov,fov_tensor,CCD,oovv_tensor,t1_tensor)
!       input variables
        type(tensor_t), intent(inout) :: fbar_ov
        type(tensor_t), intent(inout) :: fov_tensor
        type(tensor_t), intent(inout), optional :: oovv_tensor, t1_tensor
        logical,            intent(in   ) :: CCD
!       error code
        integer :: ierr

        ierr=initialize_tensor_t(fbar_ov)
        ierr=t_einsum("F(m,e)+=U(m,e)",fbar_ov,fov_tensor)
        if (.not.CCD) ierr=t_einsum("F(m,e)+=V(m,n,e,f)*T(f,n)",fbar_ov,oovv_tensor,t1_tensor)

       end subroutine get_fbar_ov_eti_scratch
! this is not the same as get_fbar_vv_talsh in intermediates.F90 since it explicitly constructs Fbar_vv from
! fock matrix, integral and t tensors, in line with Shee et al JCP 2018 mentioned above (with potential corrections)
       subroutine get_fbar_vv_eti_scratch(fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: fbar_vv
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: fvv_tensor, fov_tensor, t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediate
        type(tensor_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr


        integer                :: f_vv_rel
        real(8),allocatable,dimension(:) :: Fbar_ae
        complex(8), pointer    ::fbarvv_tens(:,:), fbarvv_tens1(:,:)
        integer(INTD),dimension(2)   :: vv_dims
        integer :: i,j,n,m
        integer :: i_remain, j_remain, i_quotient, j_quotient, i_pri, j_pri
        type(C_PTR) :: body_p
        complex(8) :: q

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=initialize_tensor_t(fbar_vv)
        ierr=t_einsum("F(e,a)+=U(e,a)",fbar_vv,fvv_tensor)
        if (.not.CCD) then
          ierr=t_einsum("F(e,a)+=U(m,a)*T(e,m)",fbar_vv,fov_tensor,t1_tensor,scale=MINUS_ONE)
          ierr=t_einsum("F(e,a)+=V(e,m,a,f)*T(f,m)",fbar_vv,int_t%vovv,t1_tensor)
        end if

        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=t_einsum("F(e,a)+=V(m,n,a,f)*T(f,e,m,n)",fbar_vv,int_t%oovv,tau_tensor,scale=ONE_HALF)
        ierr=destroy_tensor_t(tau_tensor)

       end subroutine get_fbar_vv_eti_scratch

! this is not the same as get_w_oooo_talsh in intermediates.F90 since it does not 
! use comm_t (so acually recalculates tau internally), and also the first contractions
! sort of look like those coming from get_a_talsh but the sign isn't the same. 
       subroutine get_w_oooo_eti_scratch(w_oooo,int_t,t2_tensor,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: w_oooo
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediate
        type(tensor_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=t_einsum("W(i,j,m,n)+=V(i,j,m,n)",w_oooo,int_t%oooo)
        if (.not.CCD) then
          ierr=t_einsum("W(i,j,m,n)+=V(j,i,n,e)*T(e,m)",w_oooo,int_t%ooov,t1_tensor)
          ierr=t_einsum("W(i,j,m,n)+=V(i,j,m,e)*T(e,n)",w_oooo,int_t%ooov,t1_tensor)
        end if

        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=t_einsum("W(i,j,m,n)+=V(i,j,e,f)*T(e,f,m,n)",w_oooo,int_t%oovv,tau_tensor,scale=ONE_HALF)

        ierr=destroy_tensor_t(tau_tensor)

       end subroutine get_w_oooo_eti_scratch

! this is almost the same as get_w_ooov_talsh: equations are the same
! but making t1_tensor as optional, the CCD is no longer necessary
! as an argument.
       subroutine get_w_ooov_talsh(w_ooov,int_t,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: w_ooov
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=t_einsum("W(m,n,i,e)+=V(m,n,i,e)",w_ooov,int_t%ooov)
        if (.not.CCD) ierr=t_einsum("W(m,n,i,e)+=T(f,i)*V(m,n,f,e)",w_ooov,t1_tensor,int_t%oovv)
       
       end subroutine get_w_ooov_talsh

! this is similar to get_w_ovoo_talsh in intermediates.F90 but it dues not use comm_t or passes CCD
! as argument, and actually explicitly creates and initializes one_tensor to one (which we probably don't need
! as one can do a copy). one does not use lambda_t either (rather storing results in w_ovoo, and explicitly
! using fbar_ov etc) but equations are the same.
       subroutine get_w_ovoo_eti_scratch(w_ovoo,int_t,fbar_ov,w_oooo,wbar_vovo,t2_tensor,t1_tensor)
        type(tensor_t),    intent(inout) :: w_ovoo
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: fbar_ov, w_oooo, wbar_vovo, t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(tensor_t) :: tau_tensor, one_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4), one_dims(1)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        one_dims(1) = 1
        ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
        ierr=t_einsum("W(i,e,m,n)+=V+(m,n,i,e)*K()",w_ovoo,int_t%ooov,one_tensor)
        ierr=destroy_tensor_t(one_tensor)

        ierr=t_einsum("W(i,e,m,n)+=F(i,f)*T(f,e,m,n)",w_ovoo,fbar_ov,t2_tensor)
!       Need to permute m and n, so term appears twice
        ierr=t_einsum("W(i,e,m,n)+=V(i,o,m,f)*T(e,f,n,o)",w_ovoo,int_t%ooov,t2_tensor)
        ierr=t_einsum("W(i,e,m,n)+=V(o,i,n,f)*T(e,f,m,o)",w_ovoo,int_t%ooov,t2_tensor)
        if (.not.CCD) then 
!        Need to permute m and n, so term appears twice
         ierr=t_einsum("W(i,e,m,n)+=R(e,i,f,n)*T(f,m)",w_ovoo,wbar_vovo,t1_tensor,scale=MINUS_ONE)
         ierr=t_einsum("W(i,e,m,n)+=R(e,i,f,m)*T(f,n)",w_ovoo,wbar_vovo,t1_tensor)
         ierr=t_einsum("W(i,e,m,n)+=R(o,i,m,n)*T(e,o)",w_ovoo,w_oooo,t1_tensor)
        end if

        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if
        ierr=t_einsum("W(i,e,m,n)+=V(e,i,f,g)*T(g,f,m,n)",w_ovoo,int_t%vovv,tau_tensor,scale=ONE_HALF)
        ierr=destroy_tensor_t(tau_tensor)

       end subroutine get_w_ovoo_eti_scratch


! this is similar to get_w_vovo_talsh but without using comm_t, or passing on the logical CCD
       subroutine get_w_vovo_eti_scratch(w_vovo,int_t,t2_tensor,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: w_vovo
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(tensor_t) :: tau2_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        ierr=t_einsum("W(b,m,e,j)+=V(b,m,e,j)",w_vovo,int_t%vovo)
        if (.not.CCD) then
          ierr=t_einsum("W(b,m,e,j)+=V(b,m,e,f)*T(f,j)",w_vovo,int_t%vovv,t1_tensor)
          ierr=t_einsum("W(b,m,e,j)+=V(n,m,j,e)*T(b,n)",w_vovo,int_t%ooov,t1_tensor)
        end if
 
        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
        ierr=t_einsum("T(f,b,j,n)+=R(f,b,j,n)",tau2_tensor,t2_tensor)
        if (.not.CCD) ierr=t_einsum("T(f,b,j,n)+=T(f,j)*T(b,n)",tau2_tensor,t1_tensor,t1_tensor)   !XY
        ierr=t_einsum("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,int_t%oovv,tau2_tensor)
        ierr=destroy_tensor_t(tau2_tensor)

       end subroutine get_w_vovo_eti_scratch

! this is quite similar to get_w_vvvo_talsh but without using comm_t or lambda_t explicitly, and that tau is
! recalculated internally. 
! as mentioned before  we can perhaps get rid of one tensor by just using add or copy
       subroutine get_w_vvvo_eti_no_dt(w_vvvo,int_t,w_vvvv,wbar_vovo,fbar_ov,t2_tensor,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: w_vvvo
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: w_vvvv, wbar_vovo, fbar_ov, t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(tensor_t) :: tau_tensor, one_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4), one_dims(1)
!       local variables
        logical :: CCD = .false.
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        one_dims(1) = 1
        ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
        ierr=t_einsum("W(e,f,a,m)+=V+(a,m,e,f)*K()",w_vvvo,int_t%vovv,one_tensor)
        ierr=destroy_tensor_t(one_tensor)

!       Need to permute e and f, so term appears twice
        ierr=t_einsum("W(e,f,a,m)+=V(e,n,a,g)*T(f,g,m,n)",w_vvvo,int_t%vovv,t2_tensor)
        ierr=t_einsum("W(e,f,a,m)+=V(f,n,a,g)*T(g,e,m,n)",w_vvvo,int_t%vovv,t2_tensor)
        ierr=t_einsum("W(e,f,a,m)+=F(n,a)*T(e,f,m,n)",w_vvvo,fbar_ov,t2_tensor)

        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if
        ierr=t_einsum("W(e,f,a,m)+=V(n,o,m,a)*T(e,f,o,n)",w_vvvo,int_t%ooov,tau_tensor,scale=ONE_HALF)
        ierr=destroy_tensor_t(tau_tensor)

        if (.not.CCD) then
!        Need to permute e and f, so term appears twice
         ierr=t_einsum("W(e,f,a,m)+=R(f,n,a,m)*T(e,n)",w_vvvo,wbar_vovo,t1_tensor)
         ierr=t_einsum("W(e,f,a,m)+=R(e,n,a,m)*T(f,n)",w_vvvo,wbar_vovo,t1_tensor,scale=MINUS_ONE)
         ierr=t_einsum("W(e,f,a,m)+=R(e,f,a,g)*T(g,m)",w_vvvo,w_vvvv,t1_tensor)
        end if

       end subroutine get_w_vvvo_eti_no_dt

! this is nearly the same as get_w_vovv_talsh in intermediates.F90, apart
! from removing the logical CCD from the argument list by making t1 optional
       subroutine get_w_vovv_eti_no_logical(w_vovv,int_t,t1_tensor)
!       input variables
        type(tensor_t),    intent(inout) :: w_vovv
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       CCD switch
        logical :: CCD = .false.
!       error code
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        if (.not.CCD) ierr=t_einsum("W(a,m,e,f)+=V(m,n,e,f)*T(a,n)",w_vovv,int_t%oovv,t1_tensor)
        ierr=t_einsum("W(a,m,e,f)+=V(a,m,e,f)",w_vovv,int_t%vovv)

       end subroutine get_w_vovv_eti_no_logical

! this is nearly the same as get_w_vvvv_talsh in intermediates.F90 with respect to the actual 
! contractions, but not using comm_t (so tau is recalculated) and removing the CCD logical from the argument list
       subroutine get_w_vvvv_eti_no_dt(w_vvvv,int_t,t2_tensor,t1_tensor)

!       input variables
        type(tensor_t),    intent(inout) :: w_vvvv
        type(eti_intg_tens_t), intent(inout) :: int_t
        type(tensor_t),    intent(inout) :: t2_tensor
        type(tensor_t),    intent(inout), optional :: t1_tensor
!       local intermediates
        type(tensor_t) :: tau_tensor
        integer            :: tens_rank
        integer            :: vvoo_dims(4)
!       local variables
        logical :: CCD = .false.
        integer :: ierr

        if (.not.present(t1_tensor)) CCD = .true.

        if (.not.CCD) then
!        Permute e and f, so term appears twice
         ierr=t_einsum("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv,int_t%vovv,t1_tensor)
         ierr=t_einsum("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv,int_t%vovv,t1_tensor,scale=MINUS_ONE)
        end if

        ierr = get_nmodes(t2_tensor,tens_rank)
        ierr = get_extents_static(t2_tensor,vvoo_dims)

        ierr=construct_tensor_t(tau_tensor,C8,vvoo_dims,init_val=ZERO)
        if (.not.CCD) then
          call get_tau(tau_tensor,t2_tensor,t1_tensor)
        else
          call get_tau(tau_tensor,t2_tensor)
        end if

        ierr=t_einsum("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv,int_t%oovv,tau_tensor,scale=ONE_HALF)
        ierr=t_einsum("W(e,f,a,b)+=V(e,f,a,b)",w_vvvv,int_t%vvvv)

        ierr=destroy_tensor_t(tau_tensor)

       end subroutine get_w_vvvv_eti_no_dt

       subroutine form_ccsd_intermediates_talsh(interm_t, CCD, nocc, nvir, comm_t, int_t, l2_tensor)
! input variables
          type(eti_intermediates_tens_t),  intent(inout) :: interm_t
! CCD switch
          logical, intent(in) :: CCD
          integer, intent(in) :: nocc, nvir
          type(eti_comm_tens_t),    intent(inout) :: comm_t
          type(eti_intg_tens_t), intent(inout) :: int_t
          type(tensor_t),    intent(inout), optional :: l2_tensor 
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
         ierr=construct_tensor_t(interm_t%t1,C8,vo_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%t2,C8,vvoo_dims,init_val=ZERO)

         ierr=t_einsum("IntmT(a,i)+=T(a,i)",interm_t%t1,comm_t%t1)
         ierr=t_einsum("IntmT(a,b,i,j)+=T(a,b,i,j)",interm_t%t2,comm_t%t2)
         print*, "copying t2 tensor, norm ", tensor_norm2(comm_t%t2),tensor_norm2(interm_t%t2)

!  create all fixed intermediates
         ierr=construct_tensor_t(interm_t%fbar_oo,C8,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%fbar_ov,C8,ov_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%fbar_vv,C8,vv_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%w_oooo,C8,oooo_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%w_ooov,C8,ooov_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%w_vovo,C8,vovo_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

         ierr=construct_tensor_t(interm_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
         ierr=construct_tensor_t(interm_t%w_vovv,C8,vovv_dims,init_val=ZERO)

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
            ierr=construct_tensor_t(interm_t%goo,C8,oo_dims,init_val=ZERO)
            ierr=construct_tensor_t(interm_t%gvv,C8,vv_dims,init_val=ZERO)
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
         

         ierr=destroy_tensor_t(interm_t%wbar_vovo)

      end subroutine

      subroutine destroy_ccsd_intermediates_talsh(interm_t)
! input variables
         type(eti_intermediates_tens_t),  intent(inout) :: interm_t
         integer :: ierr

         print*, "destroying ccsd intermediates for talsh"
         ierr=destroy_tensor_t(interm_t%goo)
         ierr=destroy_tensor_t(interm_t%gvv)
         ierr=destroy_tensor_t(interm_t%fbar_oo)
         ierr=destroy_tensor_t(interm_t%fbar_ov)
         ierr=destroy_tensor_t(interm_t%fbar_vv)
         ierr=destroy_tensor_t(interm_t%w_oooo)
         ierr=destroy_tensor_t(interm_t%w_ooov)
         ierr=destroy_tensor_t(interm_t%w_ovoo)
         ierr=destroy_tensor_t(interm_t%w_vovo)
         ierr=destroy_tensor_t(interm_t%w_vovv)
         ierr=destroy_tensor_t(interm_t%w_vvvo)
         ierr=destroy_tensor_t(interm_t%w_vvvv)
         print*, "done destroying ccsd intermediates"

      end subroutine


      subroutine get_eta_ov_eti(eta_ov,fbar_ov,ov_dims)
        type(tensor_t),intent(inout)   :: eta_ov
        type(tensor_t),intent(inout)   :: fbar_ov
        integer :: ierr
        integer(INTD), dimension(2)   :: ov_dims

        ierr=construct_tensor_t(eta_ov,C8,ov_dims,init_val=ZERO)
        ierr=t_einsum("N(i,a)+=F(i,a)",eta_ov,fbar_ov)
      end subroutine


      subroutine get_eta_oovv_talsh(eta_oovv,int_oovv,oovv_dims)
        type(tensor_t),intent(inout)   :: eta_oovv
        type(tensor_t),intent(inout)   :: int_oovv
        integer :: ierr
        integer(INTD), dimension(4)   :: oovv_dims

        ierr=construct_tensor_t(eta_oovv,C8,oovv_dims,init_val=ZERO)
        ierr=t_einsum("N(i,j,a,b)+=V(i,j,a,b)",eta_oovv,int_oovv)
      end subroutine get_eta_oovv_talsh

end module
