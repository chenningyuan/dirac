module talsh_gradient

!This module contains routines for solving the Lambda equations

    use talsh
    use tensor_algebra

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE_QUARTER=(-0.25D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0)

    private

    public get_dm_mo

    contains

      subroutine get_dm_mo(t1_tensor,t2_tensor,l1_tensor,l2_tensor,nocc,nvir)

!      Routine for obtaining the density matrices in the MO basis

       use talsh_ao_to_mo, only       : tensor_norm2

       integer, intent(in) :: nocc, nvir
!      solution tensors
       type(talsh_tens_t),intent(inout) :: t1_tensor, t2_tensor
       type(talsh_tens_t),intent(inout) :: l1_tensor, l2_tensor
!      gamma tensors
       type(talsh_tens_t) :: gprime_oo, gprime_ov, gprime_vo, gprime_vv
       type(talsh_tens_t) :: gsymm_oo, gsymm_vv
!      Aux tensor
       type(talsh_tens_t) :: one_tensor
       type(talsh_tens_t) :: tau2_tensor
       type(talsh_tens_t) :: ooov_aux_tensor,oo_aux_tensor
!      Tensor dimensions
       integer :: one_dims(1)=1
       integer, dimension(2) :: oo_dims,ov_dims,vo_dims,vv_dims
       integer, dimension(4) :: vvoo_dims,ooov_dims

       integer :: ierr

       call print_date('Entered DM routine')

!      Construct gamma prime tensors
       oo_dims = nocc * 2
       ierr=talsh_tensor_construct(gprime_oo,C8,oo_dims,init_val=ZERO)
       ierr=talsh_tensor_construct(gsymm_oo,C8,oo_dims,init_val=ZERO)
       ov_dims(1) = nocc * 2
       ov_dims(2) = nvir * 2
       ierr=talsh_tensor_construct(gprime_ov,C8,ov_dims,init_val=ZERO)
       vo_dims(1) = nvir * 2
       vo_dims(2) = nocc * 2
       ierr=talsh_tensor_construct(gprime_vo,C8,vo_dims,init_val=ZERO)
       vv_dims = nvir * 2
       ierr=talsh_tensor_construct(gprime_vv,C8,vv_dims,init_val=ZERO)
       ierr=talsh_tensor_construct(gsymm_vv,C8,vv_dims,init_val=ZERO)

!      Construct auxiliary tensor
       ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
       vvoo_dims(1:2) = nvir * 2
       vvoo_dims(3:4) = nocc * 2
       ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
       ooov_dims(1:3) = nocc * 2
       ooov_dims(4)   = nvir * 2
       ierr=talsh_tensor_construct(ooov_aux_tensor,C8,ooov_dims,init_val=ZERO)
       ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)

!      oo
       ierr=talsh_tensor_contract("G(j,i)+=T(e,f,i,m)*L(m,j,e,f)",gprime_oo,t2_tensor,l2_tensor,scale=ONE_HALF)
       ierr=talsh_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gprime_oo,t1_tensor,l1_tensor,scale=MINUS_ONE)

!      vv
       ierr=talsh_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gprime_vv,t2_tensor,l2_tensor,scale=ONE_HALF)
       ierr=talsh_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gprime_vv,t1_tensor,l1_tensor)

!      vo
       ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gprime_vo,t1_tensor)
       ierr=talsh_tensor_add("T(a,e,i,m)+=R(a,e,i,m)",tau2_tensor,t2_tensor)
       ierr=talsh_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gprime_vo,l1_tensor,tau2_tensor)
       ierr=talsh_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gprime_vo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE_QUARTER)
       ierr=talsh_tensor_contract("X(m,n,i,f)+=L(m,n,e,f)*T(e,i)",ooov_aux_tensor,l2_tensor,t1_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(m,n,i,f)*T(a,f,m,n)",gprime_vo,ooov_aux_tensor,t2_tensor,scale=MINUS_ONE_QUARTER)

!      ov
       ierr=talsh_tensor_add("G(i,a)+=L(i,a)",gprime_ov,l1_tensor)

!      vv sym
       ierr=talsh_tensor_add("G(b,a)+=X(b,a)",gsymm_vv,gprime_vv,scale=ONE_HALF)
       ierr=talsh_tensor_contract("G(b,a)+=X(a,b)*K()",gsymm_vv,gprime_vv,one_tensor,scale=ONE_HALF)

!      oo sym
       ierr=talsh_tensor_add("G(j,i)+=X(j,i)",gsymm_oo,gprime_oo,scale=ONE_HALF)
       ierr=talsh_tensor_contract("G(j,i)+=X(i,j)*K()",gsymm_oo,gprime_oo,one_tensor,scale=ONE_HALF)
       
       if (ierr.ne.TALSH_SUCCESS) then
          print*," ccdriver: error in evaluation of density matrices, code=",ierr
          return
       end if

       call print_date("Calculated density matrices in MO basis")
       if (.true.) then
          print*, "gprime_oo  = ", tensor_norm2(gprime_oo)
          print*, "gprime_ov  = ", tensor_norm2(gprime_ov)
          print*, "gprime_vo  = ", tensor_norm2(gprime_vo)
          print*, "gprime_vv  = ", tensor_norm2(gprime_vv)
          print*, "gsymm_oo   = ", tensor_norm2(gsymm_oo)
          print*, "gsymm_vv   = ", tensor_norm2(gsymm_vv)
       end if

       ierr=talsh_tensor_destruct(gprime_oo)
       ierr=talsh_tensor_destruct(gprime_ov)
       ierr=talsh_tensor_destruct(gprime_vo)
       ierr=talsh_tensor_destruct(gprime_vv)
       ierr=talsh_tensor_destruct(tau2_tensor)
       ierr=talsh_tensor_destruct(ooov_aux_tensor)
       ierr=talsh_tensor_destruct(oo_aux_tensor)
       ierr=talsh_tensor_destruct(one_tensor)

      end subroutine get_dm_mo 

end module talsh_gradient
