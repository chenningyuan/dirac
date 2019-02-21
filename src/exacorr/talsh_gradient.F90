module talsh_gradient

! This module contains routines to calculate density matrices
! See Shee, A.; Visscher, L.; Saue, T. J. Chem. Phys. 2016, 145, 184107.
! Written by Stan Papadopoulos, February 2019

  use talsh
  use tensor_algebra
  use talsh_ao_to_mo, only    : tensor_norm2
  use exacorr_datatypes, only : talsh_dm_tens, exacc_input

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE_QUARTER=(-0.25D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0)

  private

  public talsh_dm_driver

  contains

    subroutine talsh_dm_driver(exa_input,t1_tensor,t2_tensor,l1_tensor,l2_tensor)

!     Routine for obtaining the density matrices in the MO basis

      use exacorr_global, only : get_nao
      
!     input variables
      type(exacc_input), intent(in   ) :: exa_input
      type(talsh_tens_t),intent(inout) :: t1_tensor, t2_tensor ! solution tensors
      type(talsh_tens_t),intent(inout) :: l1_tensor, l2_tensor ! solution tensors
!     density matrices type, output
      type(talsh_dm_tens) :: gamma
      type(talsh_tens_t)  :: dm_ao
!     tensor dimensions
      integer :: nao
      integer, dimension(2) :: oo_dims,ov_dims,vo_dims,vv_dims
      integer, dimension(2) :: dm_ao_dims
!     error code
      integer :: ierr

      call print_date('Entered DM routine')

!************************************
!     construct gamma prime tensors *
!************************************
 
      oo_dims = exa_input%nkr_occ
      ierr=talsh_tensor_construct(gamma%prime_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(gamma%symm_oo,C8,oo_dims,init_val=ZERO)
      ov_dims(1) = exa_input%nkr_occ
      ov_dims(2) = exa_input%nkr_vir
      ierr=talsh_tensor_construct(gamma%prime_ov,C8,ov_dims,init_val=ZERO)
      vo_dims(1) = exa_input%nkr_vir
      vo_dims(2) = exa_input%nkr_occ
      ierr=talsh_tensor_construct(gamma%prime_vo,C8,vo_dims,init_val=ZERO)
      vv_dims = exa_input%nkr_vir
      ierr=talsh_tensor_construct(gamma%prime_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(gamma%symm_vv,C8,vv_dims,init_val=ZERO)

!**************************************
!     calculate DM blocks in MO basis *
!**************************************

      call calc_dm_mo(gamma,t1_tensor,t2_tensor,l1_tensor,l2_tensor,exa_input%nkr_occ,exa_input%nkr_vir)
      call print_date("Calculated density matrices in MO basis")

!******************************
!     MO to AO transformation *
!******************************

!     construct tensor
      nao        = get_nao()
      dm_ao_dims = nao
      ierr=talsh_tensor_construct(dm_ao,C8,dm_ao_dims,init_val=ZERO)

      call mo_to_ao(exa_input,gamma,dm_ao)
      call print_date("MO to AO transformation completed")

!     cleanup
      ierr=talsh_tensor_destruct(gamma%prime_oo)
      ierr=talsh_tensor_destruct(gamma%prime_ov)
      ierr=talsh_tensor_destruct(gamma%prime_vo)
      ierr=talsh_tensor_destruct(gamma%prime_vv)
      ierr=talsh_tensor_destruct(gamma%symm_oo)
      ierr=talsh_tensor_destruct(gamma%symm_vv)

!********************
!     write to disk *
!********************

      call write_to_something(dm_ao,nao)
       
      if (ierr.ne.TALSH_SUCCESS) then
        print*," ccdriver: error in evaluation of density matrices, code=",ierr
        return
      end if

      ierr=talsh_tensor_destruct(dm_ao)

    end subroutine talsh_dm_driver

    subroutine calc_dm_mo(gamma,t1_tensor,t2_tensor,l1_tensor,l2_tensor,nocc,nvir)

      use exacorr_global, only : debug_exa

!     input variables
      type(talsh_dm_tens), intent(inout)  :: gamma
      type(talsh_tens_t),  intent(inout)  :: t1_tensor, t2_tensor
      type(talsh_tens_t),  intent(inout)  :: l1_tensor, l2_tensor
      integer,             intent(in)     :: nocc, nvir
!     auxiliary tensors
      type(talsh_tens_t) :: one_tensor
      type(talsh_tens_t) :: tau2_tensor
      type(talsh_tens_t) :: ooov_aux_tensor,oo_aux_tensor
!     tensor dimensions
      integer :: one_dims(1)=1
      integer, dimension(2) :: oo_dims
      integer, dimension(4) :: vvoo_dims, ooov_dims
!     error code
      integer :: ierr

!*********************************
!     construct auxiliary tensor *
!*********************************
 
      ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
      vvoo_dims(1:2) = nvir
      vvoo_dims(3:4) = nocc
      ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
      ooov_dims(1:3) = nocc
      ooov_dims(4)   = nvir
      ierr=talsh_tensor_construct(ooov_aux_tensor,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)

!********************
!     calculate DMs *
!*********************     

!----------------------------------------------------------------
!     oo: g'(j,i) = - 1/2 t(ef,im) * l(jm,ef) - t(e,i) * l(j,e) |
!----------------------------------------------------------------
      ierr=talsh_tensor_contract("G(j,i)+=T(e,f,i,m)*L(m,j,e,f)",gamma%prime_oo,t2_tensor,l2_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gamma%prime_oo,t1_tensor,l1_tensor,scale=MINUS_ONE)

!--------------------------------------------------------------
!     vv: g'(a,b) = 1/2 t(ae,mn) * l(mn,be) + t(a,m) * l(m,b) |
!--------------------------------------------------------------
      ierr=talsh_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gamma%prime_vv,t2_tensor,l2_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gamma%prime_vv,t1_tensor,l1_tensor)

!----------------------------
!     vo: g'(a,i) += t(a,i) |
!----------------------------
      ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)

!-----------------------------------------------------------
!     vo: g'(a,i) += l(m,e) * (t(ae,im) - t(e,i) * t(a,m)) |
!-----------------------------------------------------------
      ierr=talsh_tensor_add("T(a,e,i,m)+=R(a,e,i,m)",tau2_tensor,t2_tensor)
      ierr=talsh_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)

!----------------------------------------------------------------------------
!     vo: g'(a,i) += 1/4 l(mn,ef) * (t(ef,in) * t(a,m) + t(e,i) * t(af,mn)) |
!----------------------------------------------------------------------------
      ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,tau2_tensor)
      ierr=talsh_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
      ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE_QUARTER)
      ierr=talsh_tensor_contract("X(m,n,i,f)+=L(m,n,e,f)*T(e,i)",ooov_aux_tensor,l2_tensor,t1_tensor)
      ierr=talsh_tensor_contract("G(a,i)+=X(m,n,i,f)*T(a,f,m,n)",gamma%prime_vo,ooov_aux_tensor,t2_tensor,scale=MINUS_ONE_QUARTER)

!---------------------------
!     ov: g'(i,a) = l(i,a) |
!---------------------------
      ierr=talsh_tensor_add("G(i,a)+=L(i,a)",gamma%prime_ov,l1_tensor)

!------------------------------------------------
!     vv symm: g(b,a) = 1/2 (g'(b,a) + g'(a,b)) |
!------------------------------------------------
      ierr=talsh_tensor_add("G(b,a)+=X(b,a)",gamma%symm_vv,gamma%prime_vv,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(b,a)+=X(a,b)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)

!------------------------------------------------
!     oo symm: g(j,i) = 1/2 (g'(j,i) + g'(i,j)) |
!------------------------------------------------
      ierr=talsh_tensor_add("G(j,i)+=X(j,i)",gamma%symm_oo,gamma%prime_oo,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(j,i)+=X(i,j)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)

!     debug
      if (debug_exa) then
        print*, "gprime_oo  = ", tensor_norm2(gamma%prime_oo)
        print*, "gprime_ov  = ", tensor_norm2(gamma%prime_ov)
        print*, "gprime_vo  = ", tensor_norm2(gamma%prime_vo)
        print*, "gprime_vv  = ", tensor_norm2(gamma%prime_vv)
        print*, "gsymm_oo   = ", tensor_norm2(gamma%symm_oo)
        print*, "gsymm_vv   = ", tensor_norm2(gamma%symm_vv)
      end if

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(tau2_tensor)
      ierr=talsh_tensor_destruct(ooov_aux_tensor)
      ierr=talsh_tensor_destruct(oo_aux_tensor)
      ierr=talsh_tensor_destruct(one_tensor)

    end subroutine calc_dm_mo

    subroutine mo_to_ao(exa_input,gamma,dm_ao)

      use exacorr_mo
      use exacorr_global, only : debug_exa, get_nao, get_mo_coefficients

!     input variables
      type(exacc_input),   intent(in   )  :: exa_input
      type(talsh_dm_tens), intent(inout)  :: gamma
      type(talsh_tens_t),  intent(inout)  :: dm_ao
!     store mo_coefficients
      type(cmo)                       :: cspinor
      complex(8), pointer, contiguous :: coeff_fort(:,:)
      type(talsh_tens_t)              :: coeff_a, coeff_b
!     tensor dimensions
      integer :: nao, nocc, nvir
!     DMs
      type(talsh_tens_t)  :: dm_mo
      type(talsh_tens_t)  :: dm_ht
!     DM dims
      integer :: dm_mo_dims(2), dm_ht_dims(2),coeff_dims(2)
!     body pointers
      type(C_PTR)                     :: body_p1, body_p2
      complex(8), pointer, contiguous :: dm_mo_fort(:,:), dm_block(:,:)
!     complete mo list
      integer, allocatable :: mo_list(:)
!     error code
      integer :: ierr

!***************
!     set dims *
!***************
 
      nao           = get_nao()
      nocc          = exa_input%nkr_occ
      nvir          = exa_input%nkr_vir
      dm_mo_dims    = nocc + nvir
      dm_ht_dims(1) = nocc + nvir
      dm_ht_dims(2) = nao
      coeff_dims(1) = nao
      coeff_dims(2) = nocc + nvir

!************************
!     construct tensors *
!************************

      ierr=talsh_tensor_construct(dm_mo,C8,dm_mo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(dm_ht,C8,dm_ht_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(coeff_a,C8,coeff_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(coeff_b,C8,coeff_dims,init_val=ZERO)

!***************************************
!     copy DM blocks into one large DM *
!***************************************

!     -----------
!     | OO | OV |
!     -----------
!     | VO | VV |
!     -----------

!     body acces to final DM
      ierr=talsh_tensor_get_body_access(dm_mo,body_p1,C8,0,DEV_HOST)
      call c_f_pointer(body_p1,dm_mo_fort,dm_mo_dims) ! to use <dm_mo_fort> as a regular Fortran 2d array

!     oo
      ierr=talsh_tensor_get_body_access(gamma%symm_oo,body_p2,C8,0,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array

      dm_mo_fort(1:nocc,1:nocc) = dm_block

!     vo
      ierr=talsh_tensor_get_body_access(gamma%prime_vo,body_p2,C8,0,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nvir,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array

      dm_mo_fort(1+nocc:nocc+nvir,1:nocc) = dm_block

!     ov
      ierr=talsh_tensor_get_body_access(gamma%prime_ov,body_p2,C8,0,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array

      dm_mo_fort(1:nocc,1+nocc:nocc+nvir) = dm_block

!     vv
      ierr=talsh_tensor_get_body_access(gamma%symm_vv,body_p2,C8,0,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nvir,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array

      dm_mo_fort(1+nocc:nocc+nvir,1+nocc:nocc+nvir) = dm_block

!**************************
!     get MO coefficients *
!**************************

      allocate(mo_list(nocc+nvir))
      mo_list(1:nocc          ) = exa_input%mokr_occ
      mo_list(1+nocc:nocc+nvir) = exa_input%mokr_vir
      call get_mo_coefficients(cspinor,mo_list,nocc+nvir)
      deallocate(mo_list)

!     body acces to coeff_a
      ierr=talsh_tensor_get_body_access(coeff_a,body_p1,C8,0,DEV_HOST)
      call c_f_pointer(body_p1,coeff_fort,coeff_dims) ! to use <coeff_fort> as a regular Fortran 2d array
      coeff_fort(:,:) = cspinor%coeff(:,:,1)

!     body acces to coeff_b
      ierr=talsh_tensor_get_body_access(coeff_b,body_p1,C8,0,DEV_HOST)
      call c_f_pointer(body_p1,coeff_fort,coeff_dims) ! to use <coeff_fort> as a regular Fortran 2d array
      coeff_fort(:,:) = cspinor%coeff(:,:,2)

      call dealloc_mo(cspinor)

!******************************
!     MO to AO transformation *
!******************************

!-----------------------------------------------------
!     ht1 alpha: DM_ht(p,s) += DM_MO(p,q) * Ca+(q,s) |
!-----------------------------------------------------
      ierr=talsh_tensor_contract("H(p,s)+=D(p,q)*C+(s,q)",dm_ht,dm_mo,coeff_a)

!----------------------------------------------------
!     ht2 alpha: DM_AO(r,s) += Ca(r,p) * DM_ht(p,s) |
!----------------------------------------------------
      ierr=talsh_tensor_contract("D(r,s)+=C(r,p)*H(p,s)",dm_ao,coeff_a,dm_ht)

!-----------------------------------------------------
!     ht1 beta: DM_ht(p,s) += DM_MO(p,q) * Cb+(q,s) |
!-----------------------------------------------------

!     first reset dm_ht
      ierr=talsh_tensor_init(dm_ht,ZERO)
      ierr=talsh_tensor_contract("H(p,s)+=D(p,q)*C+(s,q)",dm_ht,dm_mo,coeff_b)

!----------------------------------------------------
!     ht2 beta: DM_AO(r,s) += Cb(r,p) * DM_ht(p,s) |
!----------------------------------------------------
      ierr=talsh_tensor_contract("D(r,s)+=C(r,p)*H(p,s)",dm_ao,coeff_b,dm_ht)

!     debug
      print*, "DM_AO = ", tensor_norm2(dm_ao)

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(dm_mo)
      ierr=talsh_tensor_destruct(dm_ht)
      ierr=talsh_tensor_destruct(coeff_a)
      ierr=talsh_tensor_destruct(coeff_b)

    end subroutine mo_to_ao

    subroutine write_to_something(dm_ao,nao)

!     input variables
      type(talsh_tens_t), intent(inout) :: dm_ao
      integer,            intent(in   ) :: nao
!     file variable
      integer :: file_unit
!     body pointers
      type(C_PTR)                     :: body_p
      complex(8), pointer, contiguous :: dm_ao_fort(:,:)
!     error code
      integer :: ierr

!     openn file to read
      call get_free_fileunit(file_unit)
      open(file_unit,file="DM_AO",status="new",form="unformatted",action="write")

!     dimension info
      write(file_unit) nao

!     body acces to  DM_AO
      ierr=talsh_tensor_get_body_access(dm_ao,body_p,C8,0,DEV_HOST)
      call c_f_pointer(body_p,dm_ao_fort,(/nao,nao/)) ! to use <dm_ao_fort> as a regular Fortran 2d array

!     write fortran array of DM_AO
      write(file_unit) dm_ao_fort

      close(file_unit)

    end subroutine write_to_something

end module talsh_gradient
