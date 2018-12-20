module exacorr_gradient

!This module contains routines for solving the Lambda equations

    use exatensor

#if defined (VAR_MPI)

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE_QUARTER=(-0.25D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0)

    private

    public get_dm_mo

    contains

      subroutine get_dm_mo(t1_tensor,t2_tensor,l1_tensor,l2_tensor,dims)

!      Routine for obtaining the density matrices in the MO basis

       use exacorr_ao_to_mo, only     : print_tensornorm2
       use exacorr_datatypes, only    : space_dims

       type(space_dims),intent(in) :: dims
       !integer, intent(in) :: nocc, nvir
!      solution tensors
       type(tens_rcrsv_t),intent(inout) :: t1_tensor, t2_tensor
       type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!      gamma tensors
       type(tens_rcrsv_t) :: gprime_oo, gprime_ov, gprime_vo, gprime_vv
       type(tens_rcrsv_t) :: gsymm_oo, gsymm_vv
!      Aux tensor
       type(tens_rcrsv_t) :: one_dm_tensor
       type(tens_rcrsv_t) :: tau2_tensor
       type(tens_rcrsv_t) :: ooov_aux_tensor,oo_aux_tensor
!      Tensor dimensions
       integer, dimension(2) :: oo_id,ov_id,vo_id,vv_id
       integer, dimension(4) :: vvoo_id,ooov_id
       integer(INTL), dimension(2) :: oo_root,ov_root,vo_root,vv_root
       integer(INTL), dimension(4) :: vvoo_root,ooov_root

       integer :: ierr

       call print_date('Entered DM routine')

       oo_id = (/dims%occ_space_id, dims%occ_space_id/)
       ov_id = (/dims%occ_space_id, dims%vir_space_id/)
       vo_id = (/dims%vir_space_id, dims%occ_space_id/)
       vv_id = (/dims%vir_space_id, dims%vir_space_id/)
       oo_root = (/dims%occ_space_root, dims%occ_space_root/)
       ov_root = (/dims%occ_space_root, dims%vir_space_root/)
       vo_root = (/dims%vir_space_root, dims%occ_space_root/)
       vv_root = (/dims%vir_space_root, dims%vir_space_root/)
       ooov_id = (/dims%occ_space_id, dims%occ_space_id, dims%occ_space_id, dims%vir_space_id/)
       vvoo_id = (/dims%vir_space_id, dims%vir_space_id, dims%occ_space_id, dims%occ_space_id/)
       ooov_root = (/dims%occ_space_root, dims%occ_space_root, dims%occ_space_root, dims%vir_space_root/)
       vvoo_root = (/dims%vir_space_root, dims%vir_space_root, dims%occ_space_root, dims%occ_space_root/)

!      Construct gamma prime tensors
       ierr=exatns_tensor_create(gprime_oo,"gprime_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(gprime_ov,"gprime_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(gprime_vo,"gprime_vo",vo_id,vo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(gprime_vv,"gprime_vo",vv_id,vv_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(gsymm_oo,"gsymm_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(gsymm_vv,"gsymm_vv",vv_id,vv_root,EXA_DATA_KIND_C8)

!      Construct auxiliary tensor
       ierr=exatns_tensor_create(one_dm_tensor,"one_dm",EXA_DATA_KIND_C8)
       ierr=exatns_tensor_init(one_dm_tensor,ONE)
       ierr=exatns_tensor_create(tau2_tensor,"tau2_dm",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(ooov_aux_tensor,"ooov_aux",ooov_id,ooov_root,EXA_DATA_KIND_C8)
       ierr=exatns_tensor_create(oo_aux_tensor,"oo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)

!      oo
       ierr=exatns_tensor_contract("G(j,i)+=T(e,f,i,m)*L(m,j,e,f)",gprime_oo,t2_tensor,l2_tensor,ONE_HALF)
       ierr=exatns_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gprime_oo,t1_tensor,l1_tensor,MINUS_ONE)

!      vv
       ierr=exatns_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gprime_vv,t2_tensor,l2_tensor,ONE_HALF)
       ierr=exatns_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gprime_vv,t1_tensor,l1_tensor)

!      vo
       ierr=exatns_tensor_contract("G(a,i)+=T(a,i)*K()",gprime_vo,t1_tensor,one_dm_tensor)
       ierr=exatns_tensor_contract("T(a,e,i,m)+=R(a,e,i,m)*K()",tau2_tensor,t2_tensor,one_dm_tensor)
       ierr=exatns_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,MINUS_ONE)
       ierr=exatns_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gprime_vo,l1_tensor,tau2_tensor)
       ierr=exatns_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
       ierr=exatns_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gprime_vo,oo_aux_tensor,t1_tensor,MINUS_ONE_QUARTER)
       ierr=exatns_tensor_contract("X(m,n,i,f)+=L(m,n,e,f)*T(e,i)",ooov_aux_tensor,l2_tensor,t1_tensor)
       ierr=exatns_tensor_contract("G(a,i)+=X(m,n,i,f)*T(a,f,m,n)",gprime_vo,ooov_aux_tensor,t2_tensor,MINUS_ONE_QUARTER)

!      ov
       ierr=exatns_tensor_contract("G(i,a)+=L(i,a)*K()",gprime_ov,l1_tensor,one_dm_tensor)

!      vv sym
       ierr=exatns_tensor_contract("G(b,a)+=X(b,a)*K()",gsymm_vv,gprime_vv,one_dm_tensor,ONE_HALF)
       ierr=exatns_tensor_contract("G(b,a)+=X(a,b)*K()",gsymm_vv,gprime_vv,one_dm_tensor,ONE_HALF)

!      oo sym
       ierr=exatns_tensor_contract("G(j,i)+=X(j,i)*K()",gsymm_oo,gprime_oo,one_dm_tensor,ONE_HALF)
       ierr=exatns_tensor_contract("G(j,i)+=X(i,j)*K()",gsymm_oo,gprime_oo,one_dm_tensor,ONE_HALF)
       
       if (ierr.ne.EXA_SUCCESS) then
          print*," ccdriver: error in evaluation of density matrices, code=",ierr
          return
       end if

       call print_date("Calculated density matrices in MO basis")
       if (.true.) then
          print*, "gprime_oo  = ", print_tensornorm2(gprime_oo)
          print*, "gprime_ov  = ", print_tensornorm2(gprime_ov)
          print*, "gprime_vo  = ", print_tensornorm2(gprime_vo)
          print*, "gprime_vv  = ", print_tensornorm2(gprime_vv)
          print*, "gsymm_oo   = ", print_tensornorm2(gsymm_oo)
          print*, "gsymm_vv   = ", print_tensornorm2(gsymm_vv)
       end if

       ierr=exatns_tensor_destroy(gprime_oo)
       ierr=exatns_tensor_destroy(gprime_ov)
       ierr=exatns_tensor_destroy(gprime_vo)
       ierr=exatns_tensor_destroy(gprime_vv)
       ierr=exatns_tensor_destroy(gsymm_oo)
       ierr=exatns_tensor_destroy(gsymm_vv)
       ierr=exatns_tensor_destroy(tau2_tensor)
       ierr=exatns_tensor_destroy(ooov_aux_tensor)
       ierr=exatns_tensor_destroy(oo_aux_tensor)
       ierr=exatns_tensor_destroy(one_dm_tensor)

      end subroutine get_dm_mo 

#endif

end module exacorr_gradient
