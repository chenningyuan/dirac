module talsh_gradient

! This module contains routines to calculate density matrices
! See Shee, A.; Visscher, L.; Saue, T. J. Chem. Phys. 2016, 145, 184107.
! Written by Stan Papadopoulos, February 2019

  use talsh
  use tensor_algebra
  use talsh_ao_to_mo, only    : tensor_norm2
  use exacorr_datatypes, only : talsh_dm_tens, exacc_input
  use exacorr_utils, only     : print_date
  use, intrinsic:: ISO_C_BINDING

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0), MINUS_ONE_HALF=(-0.5D0,0.D0), &
                             MINUS_ONE_QUARTER=(-0.25D0,0.D0), ONE_QUARTER=(0.25D0,0.D0), &
                             MINUS_ONE_EIGHT=(-0.125D0,0.D0)

  private

  public talsh_dm_driver

  contains

    subroutine talsh_dm_driver(exa_input,t2_tensor,l2_tensor,t1_tensor,l1_tensor)

!     Routine for obtaining the density matrices in the MO basis

!     input variables
      type(exacc_input), intent(in   ) :: exa_input
      type(talsh_tens_t),intent(inout) :: t2_tensor, l2_tensor ! solution tensors
      type(talsh_tens_t),intent(inout), optional :: t1_tensor, l1_tensor ! solution tensors
!     density matrices type, output
      type(talsh_dm_tens) :: gamma
!     tensor dimensions
      integer(INTD), dimension(2) :: oo_dims,ov_dims,vo_dims,vv_dims
!     CCD switch
      logical :: CCD
!     error code
      integer(INTD) :: ierr

      CCD = exa_input%CCD

!************************************
!     construct gamma prime tensors *
!************************************
 
      oo_dims = exa_input%nocc
      ierr=talsh_tensor_construct(gamma%prime_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(gamma%symm_oo,C8,oo_dims,init_val=ZERO)
      if (.not.CCD) then
        ov_dims(1) = exa_input%nocc
        ov_dims(2) = exa_input%nvir
        ierr=talsh_tensor_construct(gamma%prime_ov,C8,ov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(gamma%symm_ov,C8,ov_dims,init_val=ZERO)
        vo_dims(1) = exa_input%nvir
        vo_dims(2) = exa_input%nocc
        ierr=talsh_tensor_construct(gamma%prime_vo,C8,vo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(gamma%symm_vo,C8,vo_dims,init_val=ZERO)
      end if
      vv_dims = exa_input%nvir
      ierr=talsh_tensor_construct(gamma%prime_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(gamma%symm_vv,C8,vv_dims,init_val=ZERO)

!**************************************
!     calculate DM blocks in MO basis *
!**************************************

      if (.not.CCD) then
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,exa_input%print_level,t1_tensor,l1_tensor)
      else
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,exa_input%print_level)
      end if
      call print_date("Calculated density matrices in MO basis")

!***********************
!     write DM to file *
!***********************

      call store_dm_mo(exa_input,gamma)

!     cleanup
      ierr=talsh_tensor_destruct(gamma%prime_oo)
      if (.not.CCD) then
        ierr=talsh_tensor_destruct(gamma%prime_ov)
        ierr=talsh_tensor_destruct(gamma%symm_ov)
        ierr=talsh_tensor_destruct(gamma%prime_vo)
        ierr=talsh_tensor_destruct(gamma%symm_vo)
      end if
      ierr=talsh_tensor_destruct(gamma%prime_vv)
      ierr=talsh_tensor_destruct(gamma%symm_oo)
      ierr=talsh_tensor_destruct(gamma%symm_vv)

      if (ierr.ne.TALSH_SUCCESS) then
        print*," ccdriver: error in evaluation of density matrices, code=",ierr
        return
      end if


    end subroutine talsh_dm_driver

    subroutine calc_dm_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,print_level,t1_tensor,l1_tensor)

      use talsh_ao_to_mo, only       : tensor_norm2
!     input variables
      type(talsh_dm_tens), intent(inout) :: gamma
      type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
      type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor
      integer,             intent(in   ) :: nocc, nvir, print_level
!     auxiliary tensors
      type(talsh_tens_t) :: one_tensor
      type(talsh_tens_t) :: tau2_tensor
      type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor
!     tensor dimensions
      integer(C_INT)              :: one_dims(1)=1
      integer(INTD), dimension(2) :: oo_dims, vv_dims
      integer(INTD), dimension(4) :: vvoo_dims
!     CCD switch
      logical :: CCD = .false.
!     error code
      integer :: ierr

      if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

!*********************************
!     construct auxiliary tensor *
!*********************************
 
      ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
      if (.not.CCD) then
        vvoo_dims(1:2) = nvir
        vvoo_dims(3:4) = nocc
        ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
        vv_dims = nvir
        ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
        oo_dims = nocc
        ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
      end if

      print *,'T1',tensor_norm2(t1_tensor)
      print *,'T2',tensor_norm2(t2_tensor)*0.25d0
      print *,'L1',tensor_norm2(l1_tensor)
      print *,'L2',tensor_norm2(l2_tensor)*0.25d0

!********************
!     calculate DMs *
!********************

!----------------------------------------------------------------
!     oo: g'(j,i) = - 1/2 t(ef,im) * l(jm,ef) - t(e,i) * l(j,e) |
!----------------------------------------------------------------
      ierr=talsh_tensor_init(gamma%prime_oo,ZERO)
      ierr=talsh_tensor_contract("G(j,i)+=T(e,f,i,m)*L(j,m,e,f)",gamma%prime_oo,t2_tensor,l2_tensor,scale=MINUS_ONE_HALF)
      print *,'G(j,i) a',tensor_norm2(gamma%prime_oo)
      if (.not.CCD) ierr=talsh_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gamma%prime_oo,t1_tensor,l1_tensor,scale=MINUS_ONE)
      print *,'G(j,i) b',tensor_norm2(gamma%prime_oo)

!--------------------------------------------------------------
!     vv: g'(a,b) = 1/2 t(ae,mn) * l(mn,be) + t(a,m) * l(m,b) |
!--------------------------------------------------------------
      ierr=talsh_tensor_init(gamma%prime_vv,ZERO)
      ierr=talsh_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gamma%prime_vv,t2_tensor,l2_tensor,scale=ONE_HALF)
      print *,'G(a,b) a',tensor_norm2(gamma%prime_vv)
      if (.not.CCD) ierr=talsh_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gamma%prime_vv,t1_tensor,l1_tensor)
      print *,'G(a,b) b',tensor_norm2(gamma%prime_vv)

!     no off-diagonal terms in CCD
      if (.not.CCD) then
!----------------------------
!     vo: g'(a,i) += t(a,i) |
!----------------------------
        ierr=talsh_tensor_init(gamma%prime_vo,ZERO)
        ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)
        print *,'G(a,i) a',tensor_norm2(gamma%prime_vo)

!-----------------------------------------------------------
!     vo: g'(a,i) += l(m,e) * (t(ae,im) - t(e,i) * t(a,m)) |
!-----------------------------------------------------------
        ierr=talsh_tensor_init(tau2_tensor,ZERO)
        ierr=talsh_tensor_add("T(a,e,i,m)+=R(a,e,i,m)",tau2_tensor,t2_tensor)
        print *,'G(a,i) tau2 a',tensor_norm2(tau2_tensor)
        ierr=talsh_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)
        print *,'G(a,i) tau2 b',tensor_norm2(tau2_tensor)
        ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,tau2_tensor)
        print *,'G(a,i) b',tensor_norm2(gamma%prime_vo)

!----------------------------------------------------------------------------
!     vo: g'(a,i) -= 1/2 l(mn,ef) * (t(ef,in) * t(a,m) + t(e,i) * t(af,mn)) |
!----------------------------------------------------------------------------
        ierr=talsh_tensor_init(oo_aux_tensor,ZERO)
        ierr=talsh_tensor_init(vv_aux_tensor,ZERO)
        ierr=talsh_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
        print *,'G(a,i) aux oo',tensor_norm2(oo_aux_tensor)
        ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF) 
        print *,'G(a,i) term c',tensor_norm2(gamma%prime_vo)
        ierr=talsh_tensor_contract("X(a,e)+=L(m,n,e,f)*T(a,f,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
        print *,'G(a,i) aux vv',tensor_norm2(vv_aux_tensor)
        ierr=talsh_tensor_contract("G(a,i)+=X(a,e)*T(e,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
        print *,'G(a,i) term d',tensor_norm2(gamma%prime_vo)

!---------------------------
!     ov: g'(i,a) = l(i,a) |
!---------------------------
        ierr=talsh_tensor_init(gamma%prime_ov,ZERO)
        ierr=talsh_tensor_add("G(i,a)+=L(i,a)",gamma%prime_ov,l1_tensor)
        print *,'G(i,a) term a',tensor_norm2(gamma%prime_ov)

!--------------------------------------------------
!       ov symm: g(i,a) = 1/2 (g'(i,a) + g'(a,i)) |
!--------------------------------------------------
        ierr=talsh_tensor_init(gamma%symm_ov,ZERO)
        ierr=talsh_tensor_contract("G(i,a)+=X(i,a)*K()",gamma%symm_ov,gamma%prime_ov,one_tensor,scale=ONE_HALF)
        ierr=talsh_tensor_contract("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_tensor,scale=ONE_HALF)
        print *,'G(i,a) symmetrized',tensor_norm2(gamma%symm_ov)
        

!--------------------------------------------------
!       vo symm: g(a,i) = symmetrized c.c. g(i,a)
!--------------------------------------------------
         ierr=talsh_tensor_init(gamma%symm_vo,ZERO)
         ierr=talsh_tensor_contract("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%symm_ov,one_tensor)
         print *,'G(a,i) symmetrized',tensor_norm2(gamma%symm_vo)

      end if

!------------------------------------------------
!     vv symm: g(b,a) = 1/2 (g'(b,a) + g'(a,b)) |
!------------------------------------------------
      ierr=talsh_tensor_init(gamma%symm_vv,ZERO)
      ierr=talsh_tensor_contract("G(b,a)+=X(b,a)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(b,a)+=X+(a,b)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)
      print *,'G(b,a) symmetrized',tensor_norm2(gamma%symm_vv)

!------------------------------------------------
!     oo symm: g(j,i) = 1/2 (g'(j,i) + g'(i,j)) |
!------------------------------------------------
      ierr=talsh_tensor_init(gamma%symm_oo,ZERO)
      ierr=talsh_tensor_contract("G(j,i)+=X(j,i)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(j,i)+=X+(i,j)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)
      print *,'G(j,i) symmetrized',tensor_norm2(gamma%symm_oo)

!     debug
!     if (print_level.gt.8) then
        if (.not.CCD) then
          print*, "gprime_ov  = ", tensor_norm2(gamma%prime_ov)
          print*, "gprime_vo  = ", tensor_norm2(gamma%prime_vo)
        end if
        print*, "gprime_vv  = ", tensor_norm2(gamma%prime_vv)
        print*, "gprime_oo  = ", tensor_norm2(gamma%prime_oo)

        if (.not.CCD) then
          print*, "gsymm_ov   = ", tensor_norm2(gamma%symm_ov)
          print*, "gsymm_vo   = ", tensor_norm2(gamma%symm_vo)
        end if
        print*, "gsymm_vv   = ", tensor_norm2(gamma%symm_vv)
        print*, "gsymm_oo   = ", tensor_norm2(gamma%symm_oo)
!     end if

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(one_tensor)
      if (.not.CCD) then
        ierr=talsh_tensor_destruct(tau2_tensor)
        ierr=talsh_tensor_destruct(vv_aux_tensor)
        ierr=talsh_tensor_destruct(oo_aux_tensor)
      end if

    end subroutine calc_dm_mo

    subroutine store_dm_mo (exa_input,gamma)

      use exacorr_global, only : get_nsolutions

!     input variables
      type(exacc_input),   intent(in   )  :: exa_input
      type(talsh_dm_tens), intent(inout)  :: gamma
!     matrix dimensions
      integer :: nocc, nvir, nkr_occ, nkr_vir, nesh
!     body pointers
      type(C_PTR)                     :: body_p2
      complex(8), pointer, contiguous :: dm_block(:,:)
      complex(8), allocatable         :: dm_mo_fort(:,:)
      real(8),allocatable             :: dm_diracformat(:,:,:)
!     complete mo list
      integer, allocatable :: mo_list(:)
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr
!     loop variables
      integer :: i, j, norbt

      CCD = exa_input%CCD

!***************
!     set dims *
!***************
 
      nocc    = exa_input%nocc
      nvir    = exa_input%nvir
      nkr_occ = exa_input%nkr_occ
      nkr_vir = exa_input%nkr_vir

!***************************************
!     copy DM blocks into one large DM *
!***************************************

      allocate (dm_mo_fort(nocc+nvir,nocc+nvir))
      dm_mo_fort = ZERO

!     -----------
!     | OO | OV |
!     -----------
!     | VO | VV |
!     -----------

!     oo
      ierr=talsh_tensor_get_body_access(gamma%symm_oo,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_mo_fort(1:nocc,1:nocc) = dm_block

!     vo + ov
      if (.not.CCD) then
        ierr=talsh_tensor_get_body_access(gamma%symm_vo,body_p2,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p2,dm_block,(/nvir,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_mo_fort(1+nocc:nocc+nvir,1:nocc) = dm_block

        ierr=talsh_tensor_get_body_access(gamma%symm_ov,body_p2,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p2,dm_block,(/nocc,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_mo_fort(1:nocc,1+nocc:nocc+nvir) = dm_block
      end if

!     vv
      ierr=talsh_tensor_get_body_access(gamma%symm_vv,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nvir,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_mo_fort(1+nocc:nocc+nvir,1+nocc:nocc+nvir) = dm_block

!*****************************************************************************************************************************
!     Determine where this should go in the full list of MOs (NB: DIRAC assumes restricted densities and uses Kramers pairs) *
!*****************************************************************************************************************************

      norbt = nkr_occ + nkr_vir
      allocate(mo_list(norbt))
      mo_list(1:nkr_occ)                 = exa_input%mokr_occ
      mo_list(1+nkr_occ:nkr_occ+nkr_vir) = exa_input%mokr_vir

      nesh = get_nsolutions() / 2
      allocate (dm_diracformat(nesh,nesh,4))
      dm_diracformat = real(ZERO)
      do j = 1, norbt
        do i = 1, norbt
          ! We copy the unbar-unbar part of the density matrix, would be good to test whether the Kramers-restricted assumption holds
          ! Or even to average over the two density matrices and enforce a Kramers-restricted density
          dm_diracformat(mo_list(i),mo_list(j),1) = real(dm_mo_fort(2*i-1,2*j-1)) 
          dm_diracformat(mo_list(i),mo_list(j),2) = imag(dm_mo_fort(2*i-1,2*j-1))
        end do
      end do

      call store_cc_density ('CCSD', dm_diracformat )

      deallocate (dm_diracformat)
      deallocate (mo_list)
      deallocate (dm_mo_fort)

    end subroutine store_dm_mo

end module talsh_gradient
