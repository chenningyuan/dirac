module exacorr_gradient

! This module contains routines for solving the Lambda equations

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))

  use exatensor
  use exacorr_ao_to_mo,  only : print_tensornorm2
  use exacorr_datatypes
  use exacorr_utils,     only : print_date
  use, intrinsic:: ISO_C_BINDING

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), ONE_QUARTER=(0.25D0,0.D0), &
                           MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0)

  private

  public exacorr_dm_driver

  contains

    subroutine exacorr_dm_driver(exa_input,dims,t2_tensor,l2_tensor,t1_tensor,l1_tensor)

      use talsh_gradient

!     Routine for obtaining the density matrices in the MO basis

!     input variables
      type(exacc_input), intent(in  )  :: exa_input
      type(space_dims),  intent(in   ) :: dims
      type(tens_rcrsv_t),intent(inout) :: t2_tensor, l2_tensor
      type(tens_rcrsv_t),intent(inout), optional :: t1_tensor, l1_tensor
!     density matrices type, output
      type(exatns_dm_tens) :: gamma

      type(tens_rcrsv_t)               :: result_tensor
      complex(8)                     :: result_val

      type(exatns_intg_tens)   :: int_t !! LH intent ? 

!     tensor dimensions
      integer(INTD), dimension(2) :: oo_id,  ov_id,  vo_id,  vv_id
      integer(INTL), dimension(2) :: oo_root,ov_root,vo_root,vv_root

      integer(INTD), dimension(4) :: oovv_id,   vovv_id,   ooov_id
      integer(INTL), dimension(4) :: oovv_root, vovv_root, ooov_root
      integer(INTD), dimension(4) :: vvvv_id,   oooo_id,   ovov_id
      integer(INTL), dimension(4) :: vvvv_root, oooo_root, ovov_root
      integer(INTD), dimension(4) :: vvvo_id,   ovoo_id,   vvoo_id
      integer(INTL), dimension(4) :: vvvo_root, ovoo_root, vvoo_root

      logical :: Add_Energy=.True.

!     CCD switch
      logical :: CCD
!     error code
      integer(INTD) :: ierr

      CCD = exa_input%CCD

!**********************
!     set ids + roots *
!**********************

      call t2_tensor%get_space_ids(vvoo_id,vvoo_root,ierr)
      vv_id   = vvoo_id(1:2)
      vv_root = vvoo_root(1:2)
      oo_id   = vvoo_id(3:4)
      oo_root = vvoo_root(3:4)
      if (.not.CCD) then
        vo_id   = vvoo_id(2:3)
        vo_root = vvoo_root(2:3)
        ov_id(1) = vvoo_id(3)
        ov_id(2) = vvoo_id(1)
        ov_root(1) = vvoo_root(3)
        ov_root(2) = vvoo_root(1)
      end if

!************************************
!     construct gamma prime tensors *
!************************************

      ierr=exatns_tensor_create(gamma%prime_oo,"gprime_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(gamma%prime_oo,ZERO)
      if (.not.CCD) then
        ierr=exatns_tensor_create(gamma%prime_ov,"gprime_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(gamma%symm_ov,"gsymm_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(gamma%prime_vo,"gprime_vo",vo_id,vo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(gamma%symm_vo,"gsymm_vo",vo_id,vo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(gamma%prime_ov,ZERO)
        ierr=exatns_tensor_init(gamma%symm_ov, ZERO)
        ierr=exatns_tensor_init(gamma%prime_vo,ZERO)
        ierr=exatns_tensor_init(gamma%symm_vo, ZERO)
      end if
      ierr=exatns_tensor_create(gamma%prime_vv,"gprime_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_oo,"gsymm_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_vv,"gsymm_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(gamma%prime_vv,ZERO)
      ierr=exatns_tensor_init(gamma%symm_oo,ZERO)
      ierr=exatns_tensor_init(gamma%symm_vv,ZERO)

     if (exa_input%do_2rdm_gs) then
!    Construct double

      oovv_id   = (/dims%occ_space_id,   dims%occ_space_id,   dims%vir_space_id,   dims%vir_space_id  /)
      oovv_root = (/dims%occ_space_root, dims%occ_space_root, dims%vir_space_root, dims%vir_space_root/)
      vovv_id   = (/dims%vir_space_id,   dims%occ_space_id,   dims%vir_space_id,   dims%vir_space_id  /)
      vovv_root = (/dims%vir_space_root, dims%occ_space_root, dims%vir_space_root, dims%vir_space_root/)
      ooov_id   = (/dims%occ_space_id,   dims%occ_space_id,   dims%occ_space_id,   dims%vir_space_id  /)
      ooov_root = (/dims%occ_space_root, dims%occ_space_root, dims%occ_space_root, dims%vir_space_root/)
      vvvv_id   = (/dims%vir_space_id,   dims%vir_space_id,   dims%vir_space_id,   dims%vir_space_id  /)
      vvvv_root = (/dims%vir_space_root, dims%vir_space_root, dims%vir_space_root, dims%vir_space_root/)
      oooo_id   = (/dims%occ_space_id,   dims%occ_space_id,   dims%occ_space_id,   dims%occ_space_id  /)
      oooo_root = (/dims%occ_space_root, dims%occ_space_root, dims%occ_space_root, dims%occ_space_root/)
      ovov_id   = (/dims%occ_space_id,   dims%vir_space_id,   dims%occ_space_id,   dims%vir_space_id  /)
      ovov_root = (/dims%occ_space_root, dims%vir_space_root, dims%occ_space_root, dims%vir_space_root/)
      vvvo_id   = (/dims%vir_space_id,   dims%vir_space_id,   dims%vir_space_id,   dims%occ_space_id  /)
      vvvo_root = (/dims%vir_space_root, dims%vir_space_root, dims%vir_space_root, dims%occ_space_root/)
      ovoo_id   = (/dims%occ_space_id,   dims%vir_space_id,   dims%occ_space_id,   dims%occ_space_id  /)
      ovoo_root = (/dims%occ_space_root, dims%vir_space_root, dims%occ_space_root, dims%occ_space_root/)
      vvoo_id   = (/dims%vir_space_id,   dims%vir_space_id,   dims%occ_space_id,   dims%occ_space_id  /)
      vvoo_root = (/dims%vir_space_root, dims%vir_space_root, dims%occ_space_root, dims%occ_space_root/)

      ierr=exatns_tensor_create(gamma%prime_oovv,"gprime_oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_oovv,"gsymm_oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_vovv,"gprime_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_vovv,"gsymm_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_ooov,"gprime_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_ooov,"gsymm_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_vvvv,"gprime_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_vvvv,"gsymm_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_oooo,"gprime_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_oooo,"gsymm_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_ovov,"gprime_ovov",ovov_id,ovov_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_ovov,"gsymm_ovov",ovov_id,ovov_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_vvvo,"gprime_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_vvvo,"gsymm_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_ovoo,"gprime_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_ovoo,"gsymm_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%prime_vvoo,"gprime_vvoo",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
      ierr=exatns_tensor_create(gamma%symm_vvoo,"gsymm_vvoo",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
     end if 

      ierr=exatns_tensor_init(gamma%prime_oovv,ZERO)
      ierr=exatns_tensor_init(gamma%prime_vovv,ZERO)
      ierr=exatns_tensor_init(gamma%prime_ooov,ZERO)
      ierr=exatns_tensor_init(gamma%prime_vvvv,ZERO)
      ierr=exatns_tensor_init(gamma%prime_oooo,ZERO)
      ierr=exatns_tensor_init(gamma%prime_ovov,ZERO)
      ierr=exatns_tensor_init(gamma%prime_vvvo,ZERO)
      ierr=exatns_tensor_init(gamma%prime_ovoo,ZERO)
      ierr=exatns_tensor_init(gamma%prime_vvoo,ZERO)


!**************************************
!     calculate DM blocks in MO basis *
!**************************************

      if (.not.CCD) then
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,exa_input%print_level,t1_tensor,l1_tensor)
      else
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,exa_input%print_level)
      end if
      call print_date("Calculated density matrices in MO basis")

      if (exa_input%do_2rdm_gs) then
      call calc_dmdouble_mo(gamma,t2_tensor,l2_tensor,dims,CCD,exa_input%print_level,t1_tensor,l1_tensor)

!***********************
! ASPG, LH : Check density matrices CCSD with the calculation of the energy
!***********************

      print *, " "
      print *, "          *****************************************"
      print *, "          *   LH - 072022                         *"
      print *, "          * PROBLEM int_t does not include        *"
      print *, "          * foo, fov, fvo and no intent !!!       *"
      print *, "          *****************************************"
      print *, " "



      print *, " "
      print *, "          *****************************************"
      print *, "          * ASP. Gomes, L. Halbert  07-2022       *"
      print *, "          * Construction of CCSD density Matrices *"
      print *, "          * Verification with Energy calculation  *"
      print *, "          *****************************************"
      print *, " "
!  Symmetry used
      ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(i,a)*F+(i,a)",result_tensor,gamma%symm_ov,int_t%fov)
      ierr=exatns_tensor_get_scalar(result_tensor,result_val)
      if (exa_input%print_level.ge.2) print*, "Final ov symm : ",result_val

      if (Add_Energy .eqv. .False.) ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(a,i)*F(i,a)",result_tensor,gamma%symm_vo,int_t%fov)
      if (exa_input%print_level.ge.2) print*, "Final vo symm : ",result_val

      if (Add_Energy .eqv. .False.) ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(a,b)*F(b,a)",result_tensor,gamma%symm_vv,int_t%fvv)
      if (exa_input%print_level.ge.2) print*, "Final vv symm : ",result_val

      if (Add_Energy .eqv. .False.) ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(i,j)*F(j,i)",result_tensor,gamma%symm_oo,int_t%foo)
      if (exa_input%print_level.ge.2) print*, "Final oo symm : ",result_val

      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(a,b,i,j)*V(i,j,a,b)",result_tensor,gamma%symm_vvoo,int_t%oovv,ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final vvoo symm : ",result_val
      
      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(i,j,a,b)*V+(i,j,a,b)",result_tensor,gamma%symm_oovv,int_t%oovv,ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final oovv symm : ",result_val
      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(i,a,j,k)*V(j,k,i,a)",result_tensor,gamma%symm_ovoo,int_t%ooov,ONE_HALF)
      if (exa_input%print_level.ge.2) print*, "Final ovoo symm : ",result_val

      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(a,b,c,i)*V(c,i,a,b)",result_tensor,gamma%symm_vvvo,int_t%vovv,ONE_HALF)
      if (exa_input%print_level.ge.2) print*, "Final vvvo symm : ",result_val

      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
!      ierr=exatns_tensor_contract("E()+=G(a,i,b,c)*V(b,c,a,i)",result_tensor,gamma%prime_vovv,int_t%vvvo,ONE_HALF)
      ierr=exatns_tensor_contract("E()+=G(a,i,b,c)*V+(a,i,b,c)",result_tensor,gamma%symm_vovv,int_t%vovv,ONE_HALF)
      if (exa_input%print_level.ge.2) print*, "Final vovv symm : ",result_val

      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(i,a,j,b)*V(b,j,a,i)",result_tensor,gamma%symm_ovov,int_t%vovo)
      if (exa_input%print_level.ge.2) print*, "Final ovov symm : ",result_val

      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(a,b,c,d)*V(c,d,a,b)",result_tensor,gamma%symm_vvvv,int_t%vvvv,ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final vvvv symm : ",result_val

      if (Add_Energy .eqv. .False.)       ierr=exatns_tensor_init(result_tensor)
      ierr=exatns_tensor_contract("E()+=G(i,j,k,l)*V(k,l,i,j)",result_tensor,gamma%symm_oooo,int_t%oooo,ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final oooo symm : ",result_val

      print *, " "
      print *, " "
      if (Add_Energy .eqv. .True.) print *, " Energy by CCSD Density Matrix ", result_val
      if (Add_Energy .eqv. .False.) print *, " Energy by CCSD Density Matrix is the sum of all the idividuale terms "
      print *, " "
      print *, " "
      print *, " ************** "
     end if




!***********************
!     write DM to file *
!***********************

      call store_dm_mo(exa_input,gamma)

!     cleanup
      ierr=exatns_tensor_destroy(gamma%prime_oo)
      if (.not.CCD) then
        ierr=exatns_tensor_destroy(gamma%prime_ov)
        ierr=exatns_tensor_destroy(gamma%symm_ov)
        ierr=exatns_tensor_destroy(gamma%prime_vo)
        ierr=exatns_tensor_destroy(gamma%symm_vo)
      end if
      ierr=exatns_tensor_destroy(gamma%prime_vv)
      ierr=exatns_tensor_destroy(gamma%symm_oo)
      ierr=exatns_tensor_destroy(gamma%symm_vv)
       
      if (ierr.ne.EXA_SUCCESS) then
        print*," ccdriver: error in evaluation of density matrices, code=",ierr
        return
      end if

    end subroutine exacorr_dm_driver

    subroutine calc_dm_mo(gamma,t2_tensor,l2_tensor,print_level,t1_tensor,l1_tensor)

!     input variables
      type(exatns_dm_tens), intent(inout) :: gamma
      type(tens_rcrsv_t),   intent(inout) :: t2_tensor, l2_tensor
      type(tens_rcrsv_t),   intent(inout), optional :: t1_tensor, l1_tensor
      integer,              intent(in   ) :: print_level
!     auxiliary tensors
      type(tens_rcrsv_t) :: one_dm_tensor
      type(tens_rcrsv_t) :: tau2_tensor
      type(tens_rcrsv_t) :: vv_aux_tensor,oo_aux_tensor
!     tensor dimensions
      integer(INTD), dimension(2) :: oo_id,   vv_id
      integer(INTL), dimension(2) :: oo_root, vv_root
      integer(INTD), dimension(4) :: vvoo_id
      integer(INTL), dimension(4) :: vvoo_root
!     error code
      integer(INTD) :: ierr
!     CCD switch
      logical       :: CCD = .false.

       if (.not. present(t1_tensor)) CCD = .true.
       call t2_tensor%get_space_ids(vvoo_id,vvoo_root,ierr)
!*********************************
!     construct auxiliary tensor *
!*********************************

      ierr=exatns_tensor_create(one_dm_tensor,"one_dm",EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(one_dm_tensor,ONE)
      if (.not.CCD) then
        vv_id   = vvoo_id(1:2)
        vv_root = vvoo_root(1:2)
        oo_id   = vvoo_id(3:4)
        oo_root = vvoo_root(3:4)

        ierr=exatns_tensor_create(tau2_tensor,"tau2_dm",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(vv_aux_tensor,"vv_aux",vv_id,vv_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_create(oo_aux_tensor,"oo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
        ierr=exatns_tensor_init(tau2_tensor,ZERO)
        ierr=exatns_tensor_init(vv_aux_tensor,ZERO)
        ierr=exatns_tensor_init(oo_aux_tensor,ZERO)
      end if

!********************
!     calculate DMs *
!********************

!----------------------------------------------------------------
!     oo: g'(j,i) = - 1/2 t(ef,im) * l(jm,ef) - t(e,i) * l(j,e) |
!----------------------------------------------------------------
      ierr=exatns_tensor_init(gamma%prime_oo,ZERO)
      ierr=exatns_tensor_contract("G(j,i)+=T(e,f,i,m)*L(m,j,e,f)",gamma%prime_oo,t2_tensor,l2_tensor,ONE_HALF)
      if (.not.CCD) ierr=exatns_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gamma%prime_oo,t1_tensor,l1_tensor,MINUS_ONE)

!--------------------------------------------------------------
!     vv: g'(a,b) = 1/2 t(ae,mn) * l(mn,be) + t(a,m) * l(m,b) |
!--------------------------------------------------------------
      ierr=exatns_tensor_init(gamma%prime_vv,ZERO)
      ierr=exatns_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gamma%prime_vv,t2_tensor,l2_tensor,ONE_HALF)
      if (.not.CCD) ierr=exatns_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gamma%prime_vv,t1_tensor,l1_tensor)

!     no off-diagonal terms in CCD
      if (.not.CCD) then
!----------------------------
!     vo: g'(a,i) += t(a,i) |
!----------------------------
        ierr=exatns_tensor_init(gamma%prime_vo,ZERO)
        ierr=exatns_tensor_contract("G(a,i)+=T(a,i)*K()",gamma%prime_vo,t1_tensor,one_dm_tensor)

!-----------------------------------------------------------
!     vo: g'(a,i) += l(m,e) * (t(ae,im) - t(e,i) * t(a,m)) |
!-----------------------------------------------------------
        ierr=exatns_tensor_init(tau2_tensor,ZERO)
        ierr=exatns_tensor_contract("T(a,e,i,m)+=R(a,e,i,m)*K()",tau2_tensor,t2_tensor,one_dm_tensor)
        ierr=exatns_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,MINUS_ONE)
        ierr=exatns_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,tau2_tensor)

!----------------------------------------------------------------------------
!     vo: g'(a,i) -= 1/2 l(mn,ef) * (t(ef,in) * t(a,m) + t(e,i) * t(af,mn)) |
!----------------------------------------------------------------------------
        ierr=exatns_tensor_init(oo_aux_tensor,ZERO) 
        ierr=exatns_tensor_init(vv_aux_tensor,ZERO) 
        ierr=exatns_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
        ierr=exatns_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,oo_aux_tensor,t1_tensor,MINUS_ONE_HALF)
        ierr=exatns_tensor_contract("X(a,e)+=L(m,n,e,f)*T(a,f,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
        ierr=exatns_tensor_contract("G(a,i)+=X(a,e)*T(e,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,MINUS_ONE_HALF)

!---------------------------
!     ov: g'(i,a) = l(i,a) |
!---------------------------
        ierr=exatns_tensor_init(gamma%prime_ov,ZERO)
        ierr=exatns_tensor_contract("G(i,a)+=L(i,a)*K()",gamma%prime_ov,l1_tensor,one_dm_tensor)

!--------------------------------------------------
!       ov symm: g(i,a) = 1/2 (g'(i,a) + g'(a,i)) |
!--------------------------------------------------
        ierr=exatns_tensor_init(gamma%symm_ov,ZERO)
        ierr=exatns_tensor_contract("G(i,a)+=X(i,a)*K()",gamma%symm_ov,gamma%prime_ov,one_dm_tensor,ONE_HALF)
        ierr=exatns_tensor_contract("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_dm_tensor,ONE_HALF)

!--------------------------------------------------
!       vo symm: g(a,i) = symmetrized c.c. g(i,a)
!--------------------------------------------------
         ierr=exatns_tensor_init(gamma%symm_vo,ZERO)
         ierr=exatns_tensor_contract("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%symm_ov,one_dm_tensor)

      end if

!------------------------------------------------
!     vv symm: g(b,a) = 1/2 (g'(b,a) + g'(a,b)) |
!------------------------------------------------
      ierr=exatns_tensor_init(gamma%symm_vv,ZERO)
      ierr=exatns_tensor_contract("G(b,a)+=X(b,a)*K()",gamma%symm_vv,gamma%prime_vv,one_dm_tensor,ONE_HALF)
      ierr=exatns_tensor_contract("G(b,a)+=X+(a,b)*K()",gamma%symm_vv,gamma%prime_vv,one_dm_tensor,ONE_HALF)

!------------------------------------------------
!     oo symm: g(j,i) = 1/2 (g'(j,i) + g'(i,j)) |
!------------------------------------------------
      ierr=exatns_tensor_init(gamma%symm_oo,ZERO)
      ierr=exatns_tensor_contract("G(j,i)+=X(j,i)*K()",gamma%symm_oo,gamma%prime_oo,one_dm_tensor,ONE_HALF)
      ierr=exatns_tensor_contract("G(j,i)+=X+(i,j)*K()",gamma%symm_oo,gamma%prime_oo,one_dm_tensor,ONE_HALF)

!     debug_exa
      if (print_level.gt.8) then
        print*, "gprime_oo  = ", print_tensornorm2(gamma%prime_oo)
        if (.not.CCD) then
          print*, "gprime_ov  = ", print_tensornorm2(gamma%prime_ov)
          print*, "gsymm_ov   = ", print_tensornorm2(gamma%symm_ov)
          print*, "gprime_vo  = ", print_tensornorm2(gamma%prime_vo)
          print*, "gsymm_vo   = ", print_tensornorm2(gamma%symm_vo)
        end if
        print*, "gprime_vv  = ", print_tensornorm2(gamma%prime_vv)
        print*, "gsymm_oo   = ", print_tensornorm2(gamma%symm_oo)
        print*, "gsymm_vv   = ", print_tensornorm2(gamma%symm_vv)
      end if

!**************
!     cleanup *
!**************

      ierr=exatns_tensor_destroy(one_dm_tensor)
      if (.not.CCD) then
        ierr=exatns_tensor_destroy(tau2_tensor)
        ierr=exatns_tensor_destroy(vv_aux_tensor)
        ierr=exatns_tensor_destroy(oo_aux_tensor)
      end if

    end subroutine calc_dm_mo

    subroutine store_dm_mo(exa_input,gamma)

      use talsh
      use talsh_gradient

!     input variables
      type(exacc_input),    intent(in   ) :: exa_input
      type(exatns_dm_tens), intent(inout) :: gamma
!     talsh type to get local copy
      type(talsh_dm_tens)          :: talsh_gamma
!     error code
      integer :: ierr
!     variables for talsh, density matrix needs to be smaller than buffer size 
      integer(C_SIZE_T)               :: buf_size=1024_8*1024_8*1024_8 
      integer(C_INT)                  :: host_arg_max
!
      character(len=4)                :: wavefun
      wavefun = "CCSD"

      ! Convert to TALSH format and call the TALSH writer
      buf_size=exa_input%talsh_buff*buf_size
      ierr=talsh_init(buf_size,host_arg_max)
      call talsh_dm_initialize_dmtens(exa_input,talsh_gamma)

!     oo + vv
      ierr=exatns_tensor_get_slice(gamma%symm_oo,talsh_gamma%symm_oo)
      ierr=exatns_tensor_get_slice(gamma%symm_vv,talsh_gamma%symm_vv)
!     vo + ov
      if (.not.exa_input%CCD) then
         ierr=exatns_tensor_get_slice(gamma%symm_vo,talsh_gamma%symm_vo)
         ierr=exatns_tensor_get_slice(gamma%symm_ov,talsh_gamma%symm_ov)
      end if

      call talsh_dm_store_mo (exa_input,talsh_gamma,wavefun)
      call talsh_kudm_store_mo(exa_input,talsh_gamma,wavefun)

      ! Clean up and switch off talsh
      call talsh_dm_cleanup_dmtens(exa_input,talsh_gamma)
      ierr = talsh_shutdown()

    end subroutine store_dm_mo

    subroutine calc_dmdouble_mo(gamma,t2_tensor,l2_tensor,dims,CCD,print_level,t1_tensor,l1_tensor)

!     input variables
      type(exatns_dm_tens), intent(inout) :: gamma
      type(tens_rcrsv_t),   intent(inout) :: t2_tensor, l2_tensor
      type(tens_rcrsv_t),   intent(inout), optional :: t1_tensor, l1_tensor
      type(space_dims),     intent(in   ) :: dims
      logical,              intent(in   ) :: CCD
      integer,              intent(in   ) :: print_level
!     auxiliary tensors
      type(tens_rcrsv_t) :: one_dm_tensor
      type(tens_rcrsv_t) :: tau2_tensor
      type(tens_rcrsv_t) :: vv_aux_tensor,oo_aux_tensor

      type(tens_rcrsv_t) :: Hinter_vovv
      type(tens_rcrsv_t) :: Hinter_ooov, Hinter_ovvo
      type(tens_rcrsv_t) :: Hinter_oooo, Hinter_vvoo



!     tensor dimensions
      integer(INTD), dimension(4) :: vvoo_id,   ooov_id
      integer(INTL), dimension(4) :: vvoo_root, ooov_root
      integer(INTD), dimension(4) :: vovv_id,   ovvo_id
      integer(INTL), dimension(4) :: vovv_root, ovvo_root
      integer(INTD), dimension(4) :: oooo_id
      integer(INTL), dimension(4) :: oooo_root



!     error code
      integer(INTD) :: ierr
!*********************************
!     construct auxiliary tensor *
!*********************************

      ierr=exatns_tensor_create(one_dm_tensor,"one_dm",EXA_DATA_KIND_C8)
      ierr=exatns_tensor_init(one_dm_tensor,ONE)

      vvoo_id   = (/dims%vir_space_id,   dims%vir_space_id,   dims%occ_space_id,   dims%occ_space_id  /)
      vvoo_root = (/dims%vir_space_root, dims%vir_space_root, dims%occ_space_root, dims%occ_space_root/)
      ooov_id   = (/dims%occ_space_id,   dims%occ_space_id,   dims%occ_space_id,   dims%vir_space_id  /)
      ooov_root = (/dims%occ_space_root, dims%occ_space_root, dims%occ_space_root, dims%vir_space_root/)
      vovv_id   = (/dims%vir_space_id,   dims%occ_space_id,   dims%vir_space_id,   dims%vir_space_id  /)
      vovv_root = (/dims%vir_space_root, dims%occ_space_root, dims%vir_space_root, dims%vir_space_root/)
      ovvo_id   = (/dims%occ_space_id,   dims%vir_space_id,   dims%vir_space_id,   dims%occ_space_id  /)
      ovvo_root = (/dims%occ_space_root, dims%vir_space_root, dims%vir_space_root, dims%occ_space_root/)
      oooo_id   = (/dims%occ_space_id,   dims%occ_space_id,   dims%occ_space_id,   dims%occ_space_id  /)
      oooo_root = (/dims%occ_space_root, dims%occ_space_root, dims%occ_space_root, dims%occ_space_root/)


      ierr=exatns_tensor_create(Hinter_vovv,"Hinter_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
!      ierr=exatns_tensor_init(Hinter_vovv,ZERO)
      ierr=exatns_tensor_create(Hinter_ooov,"Hinter_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
!      ierr=exatns_tensor_init(Hinter_ooov,ZERO)
      ierr=exatns_tensor_create(Hinter_ovvo,"Hinter_ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8)
!      ierr=exatns_tensor_init(Hinter_ovvo,ZERO)
      ierr=exatns_tensor_create(Hinter_oooo,"Hinter_ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8)
!      ierr=exatns_tensor_init(Hinter_oooo,ZERO)
      ierr=exatns_tensor_create(Hinter_vvoo,"Hinter_ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8)
!      ierr=exatns_tensor_init(Hinter_vvoo,ZERO)


      !! start gamma%prime construction 

      ! gamma oovv 
      ! j1 
!      ierr=talsh_tensor_add("G(i,j,a,b)+=L(i,j,a,b)",gamma%prime_oovv,l2_tensor) 
      ierr=exatns_tensor_contract("G(i,j,a,b)+=L(i,j,a,b)*K()",gamma%prime_oovv,l2_tensor,one_dm_tensor)
      if (ierr.ne.0)  print*, " j1, ierr  ", ierr

      ! gamma vovv
      ! h2
      ierr=exatns_tensor_contract("G(a,i,b,c)+=L(m,i,b,c)*T(a,m)",gamma%prime_vovv,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " h2, ierr  ", ierr

      !gamma ooov
      ! e2
      ierr=exatns_tensor_contract("G(i,j,k,a)+=L(i,j,e,a)*T(e,k)",gamma%prime_ooov,l2_tensor,t1_tensor,MINUS_ONE)
      if (ierr.ne.0) print*, " e2, ierr  ", ierr

      ! gamma vvvv
      ! i3
      ierr=exatns_tensor_contract("G(a,b,c,d)+=L(m,n,c,d)*T(a,b,m,n)",gamma%prime_vvvv,l2_tensor,t2_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " i3, ierr  ", ierr
      ! i4
      ierr=exatns_tensor_init(Hinter_vovv)
      ierr=exatns_tensor_contract("H(a,n,c,d)+=L(m,n,c,d)*T(a,m)",Hinter_vovv,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " i4-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,c,d)+=H(a,n,c,d)*T(b,n)",gamma%prime_vvvv,Hinter_vovv,t1_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " i4-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,c,d)+=H(b,n,c,d)*T(a,n)",gamma%prime_vvvv,Hinter_vovv,t1_tensor,MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " i4-c, ierr  ", ierr

      ! gamma oooo
      ! n3
      ierr=exatns_tensor_contract("G(i,j,k,l)+=L(i,j,e,f)*T(e,f,k,l)",gamma%prime_oooo,l2_tensor,t2_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " n3, ierr  ", ierr
      ! n4
      ierr=exatns_tensor_init(Hinter_ooov)
      ierr=exatns_tensor_contract("H(i,j,k,f)+=L(i,j,e,f)*T(e,k)",Hinter_ooov,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " n4-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(i,j,k,l)+=H(i,j,k,f)*T(f,l)",gamma%prime_oooo,Hinter_ooov,t1_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " n4-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(i,j,k,l)+=H(i,j,l,f)*T(f,k)",gamma%prime_oooo,Hinter_ooov,t1_tensor,MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " n4-c, ierr  ", ierr

      ! Start ovov f7 intermediate for g16 and k16
      ! f7
      ierr=exatns_tensor_contract("G(i,a,j,b)+=L(i,m,e,b)*T(e,a,j,m)",gamma%prime_ovov,l2_tensor,t2_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " f7, ierr  ", ierr
      ! gamma vvvo
      ! g16 sign of f7 and inversion
      ierr=exatns_tensor_contract("G(a,b,c,i)+=H(j,b,i,c)*T(a,j)",gamma%prime_vvvo,gamma%prime_ovov,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " g16-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,c,i)+=H(j,a,i,c)*T(b,j)",gamma%prime_vvvo,gamma%prime_ovov,t1_tensor)
      if (ierr.ne.0)  print*, " g16-b, ierr  ", ierr

      ! g15 + g18
      ierr=exatns_tensor_contract("G(a,b,c,i)+=H(a,c)*T(b,i)",gamma%prime_vvvo,gamma%prime_vv,t1_tensor)
      if (ierr.ne.0)  print*, "g15+g18-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,c,i)+=H(b,c)*T(a,i)",gamma%prime_vvvo,gamma%prime_vv,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, "g15+g18-b, ierr  ", ierr

      ! g17 + g19
      ierr=exatns_tensor_contract("G(a,b,c,i)+=H(a,b,c,d)*T(d,i)",gamma%prime_vvvo,gamma%prime_vvvv,t1_tensor,MINUS_ONE)

      ! g14
      ierr=exatns_tensor_contract("G(a,b,c,i)+=L(m,c)*T(a,b,m,i)",gamma%prime_vvvo,l1_tensor,t2_tensor)
      if (ierr.ne.0)  print*, " g14, ierr  ", ierr

      ! gamma ovoo
      ! k16  sign of f7
      ierr=exatns_tensor_contract("G(i,a,j,k)+=H(i,a,j,e)*T(e,k)",gamma%prime_ovoo,gamma%prime_ovov,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " k16-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(i,a,j,k)+=H(i,a,k,e)*T(e,j)",gamma%prime_ovoo,gamma%prime_ovov,t1_tensor)
      if (ierr.ne.0)  print*, " k16-b, ierr  ", ierr

      ! k15+k18
      ierr=exatns_tensor_contract("G(i,a,j,k)+=H(i,j)*T(a,k)",gamma%prime_ovoo,gamma%prime_oo,t1_tensor)
      if (ierr.ne.0)  print*, " k15+k18-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(i,a,j,k)+=H(i,k)*T(a,j)",gamma%prime_ovoo,gamma%prime_oo,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " k15+k18-b, ierr  ", ierr

      ! k17+k19
      ierr=exatns_tensor_contract("G(i,a,j,k)+=H(i,l,j,k)*T(a,l)",gamma%prime_ovoo,gamma%prime_oooo,t1_tensor)
      if (ierr.ne.0)  print*, " k19-d, ierr  ", ierr

      ! k14
      ierr=exatns_tensor_contract("G(i,a,j,k)+=L(i,e)*T(e,a,j,k)",gamma%prime_ovoo,l1_tensor,t2_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " k14, ierr  ", ierr
      ! End of gamma ovov as f7 was earlier construct
      ! f6
      ierr=exatns_tensor_contract("G(i,a,j,b)+=L(i,b)*T(a,j)",gamma%prime_ovov,l1_tensor,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " f6, ierr  ", ierr
      ! f8
      ierr=exatns_tensor_init(Hinter_ooov)
      ierr=exatns_tensor_contract("H(i,m,j,b)+=L(i,m,e,b)*T(e,j)",Hinter_ooov,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " f8-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(i,a,j,b)+=H(i,m,j,b)*T(a,m)",gamma%prime_ovov,Hinter_ooov,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " f8-b, ierr  ", ierr

      ! gamma vvoo
      ! l1
      ierr=exatns_tensor_contract("G(a,b,i,j)+=T(a,b,i,j)*K()",gamma%prime_vvoo,t2_tensor,one_dm_tensor)
      if (ierr.ne.0)  print*, " l1-a, ierr  ", ierr

      !m5+m9
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,i)*T(a,b,n,j)",gamma%prime_vvoo,gamma%prime_oo,t2_tensor)
      if (ierr.ne.0)  print*, " m9-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,j)*T(a,b,n,i)",gamma%prime_vvoo,gamma%prime_oo,t2_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " m9-c, ierr  ", ierr

      !m7+m11
      if (ierr.ne.0)  print*, " m11-c, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(a,e)*T(e,b,i,j)",gamma%prime_vvoo,gamma%prime_vv,t2_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " m11-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(b,e)*T(e,a,i,j)",gamma%prime_vvoo,gamma%prime_vv,t2_tensor)
      if (ierr.ne.0)  print*, " m11-c, ierr  ", ierr
      ! m20 
      ierr=exatns_tensor_init(Hinter_ovvo)
      ierr=exatns_tensor_contract("H(n,a,f,i)+=L(m,n,e,f)*T(a,e,i,m)",Hinter_ovvo,l2_tensor,t2_tensor)
      if (ierr.ne.0)  print*, " m20-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,a,f,i)*T(b,f,j,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " m20-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,a,f,j)*T(b,f,i,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " m20-c, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,b,f,i)*T(a,f,j,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " m20-d, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,b,f,j)*T(a,f,i,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " m20-e, ierr  ", ierr
      ! m21
      ierr=exatns_tensor_contract("A(f,b,j,n)+=T(f,j)*T(b,n)",Hinter_vvoo,t1_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " m21-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,a,f,i)*A(f,b,j,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo,MINUS_ONE)
      if (ierr.ne.0)  print*, " m21-c, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,a,f,j)*A(f,b,i,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo)
      if (ierr.ne.0)  print*, " m21-d, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,b,f,i)*A(f,a,j,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo)
      if (ierr.ne.0)  print*, " m21-e, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(n,b,f,j)*A(f,a,i,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo,MINUS_ONE)
      if (ierr.ne.0)  print*, " m21-f, ierr  ", ierr

      ! m3+m4+m6+m8 = Pij Pab \ts bj \gammas vo -l2
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(a,i)*T(b,j)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor)
      if (ierr.ne.0)  print*, " m8-c, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(a,j)*T(b,i)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " m8-d, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(b,i)*T(a,j)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, " m8-e, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(b,j)*T(a,i)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor)
      if (ierr.ne.0)  print*, " m8-f, ierr  ", ierr

      ! MINUS l2
      ierr=exatns_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " l2-a, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=T(a,j)*T(b,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " l2-b, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=T(b,i)*T(a,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " l2-c, ierr  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=T(b,j)*T(a,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " l2-d, ierr  ", ierr

      !m12+m14+m15+m24
      ierr=exatns_tensor_init(Hinter_vvoo)
      ierr=exatns_tensor_contract("H(a,b,i,j)+=T(a,b,i,j)*K()",Hinter_vvoo,t2_tensor,one_dm_tensor)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -a, ierr  ", ierr
      ierr=exatns_tensor_contract("H(a,b,i,j)+=T(a,i)*T(b,j)",Hinter_vvoo,t1_tensor,t1_tensor)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -b, ierr  ", ierr
      ierr=exatns_tensor_contract("H(a,b,i,j)+=T(a,j)*T(b,i)",Hinter_vvoo,t1_tensor,t1_tensor,MINUS_ONE)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -c, ierr  ", ierr
      ierr=exatns_tensor_init(Hinter_oooo)
      ierr=exatns_tensor_contract("H(m,n,i,j)+=L(m,n,e,f)*T(e,f,i,j)",Hinter_oooo,l2_tensor,Hinter_vvoo)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -d, ierr  ", ierr

      ierr=exatns_tensor_contract("G(a,b,i,j)+=H(m,n,i,j)*A(a,b,m,n)",&
     &                            gamma%prime_vvoo,Hinter_oooo,Hinter_vvoo,ONE_QUARTER)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -e, ierr  ", ierr

      !! gamma%prime construction

      !! Start Symetrisation
      ! gamma%symm_ooov 
      ierr=exatns_tensor_contract("G(i,j,k,a)+=G(i,j,k,a)*K()",&
     &                            gamma%symm_ooov,gamma%prime_ooov,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm ooov 1-a  ", ierr
      ierr=exatns_tensor_contract("G(i,j,k,a)+=G+(k,a,i,j)*K()",&
     &                            gamma%symm_ooov,gamma%prime_ovoo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm ooov 1-b  ", ierr

      ! gamma%symm_ovoo 
      ierr=exatns_tensor_contract("G(i,a,k,l)+=G(i,a,k,l)*K()",&
     &                            gamma%symm_ovoo,gamma%prime_ovoo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovoo 2-a  ", ierr
      ierr=exatns_tensor_contract("G(i,a,k,l)+=G+(k,l,i,a)*K()",&
     &                            gamma%symm_ovoo,gamma%prime_ooov,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovoo 2-b  ", ierr

      ! gamma%symm_ovov 
      ierr=exatns_tensor_contract("G(i,a,j,b)+=G(i,a,j,b)*K()",&
     &                            gamma%symm_ovov,gamma%prime_ovov,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovov 3-a  ", ierr
      ierr=exatns_tensor_contract("G(i,a,j,b)+=G+(j,b,i,a)*K()",&
     &                            gamma%symm_ovov,gamma%prime_ovov,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovov 3-b  ", ierr

      ! gamma%symm_vvvo 
      ierr=exatns_tensor_contract("G(a,b,c,i)+=G(a,b,c,i)*K()",&
     &                            gamma%symm_vvvo,gamma%prime_vvvo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvo 4-a  ", ierr
      ierr=exatns_tensor_contract("G(a,b,c,i)+=G+(c,i,a,b)*K()",&
     &                            gamma%symm_vvvo,gamma%prime_vovv,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvo 4-b  ", ierr

      ! gamma%symm_vovv
      ierr=exatns_tensor_contract("G(a,i,b,c)+=G(a,i,b,c)*K()",&
     &                            gamma%symm_vovv,gamma%prime_vovv,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vovv 5-a  ", ierr
      ierr=exatns_tensor_contract("G(a,i,b,c)+=G+(b,c,a,i)*K()",&
     &                            gamma%symm_vovv,gamma%prime_vvvo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vovv 5-b  ", ierr

      ! gamma%symm_vvvv
      ierr=exatns_tensor_contract("G(a,b,c,d)+=G(a,b,c,d)*K()",&
     &                            gamma%symm_vvvv,gamma%prime_vvvv,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvv 6-a  ", ierr
      ierr=exatns_tensor_contract("G(a,b,c,d)+=G+(c,d,a,b)*K()",&
     &                            gamma%symm_vvvv,gamma%prime_vvvv,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvv 6-b  ", ierr

      ! gamma%symm_oooo
      ierr=exatns_tensor_contract("G(i,j,k,l)+=G(i,j,k,l)*K()",&
     &                            gamma%symm_oooo,gamma%prime_oooo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm oooo 7-a  ", ierr
      ierr=exatns_tensor_contract("G(i,j,k,l)+=G+(k,l,i,j)*K()",&
     &                            gamma%symm_oooo,gamma%prime_oooo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm oooo 7-b  ", ierr

      ! gamma%symm_oovv 
      ierr=exatns_tensor_contract("G(i,j,a,b)+=G(i,j,a,b)*K()",&
     &                            gamma%symm_oovv,gamma%prime_oovv,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm oovv 8-a  ", ierr
      ierr=exatns_tensor_contract("G(i,j,a,b)+=G+(a,b,i,j)*K()",&
     &                            gamma%symm_oovv,gamma%prime_vvoo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm oovv 8-b  ", ierr

      ! gamma%symm_vvoo 
      ierr=exatns_tensor_contract("G(a,b,i,j)+=G(a,b,i,j)*K()",&
     &                            gamma%symm_vvoo,gamma%prime_vvoo,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvoo 9-a  ", ierr
      ierr=exatns_tensor_contract("G(a,b,i,j)+=G+(i,j,a,b)*K()",&
     &                            gamma%symm_vvoo,gamma%prime_oovv,one_dm_tensor,ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvoo 9-b  ", ierr


      ierr=exatns_tensor_destroy(Hinter_vovv)
      ierr=exatns_tensor_destroy(Hinter_ooov)
      ierr=exatns_tensor_destroy(Hinter_ovvo)
      ierr=exatns_tensor_destroy(Hinter_oooo)
      ierr=exatns_tensor_destroy(Hinter_vvoo)

   end subroutine calc_dmdouble_mo




#endif

end module exacorr_gradient
