module talsh_dm_eomccsd

! This module contains routines to calculate density matrices for EOMCC-SD
! Loic

  use talsh
  use tensor_algebra
  use talsh_ao_to_mo, only    : tensor_norm2
  use exacorr_datatypes, only : talsh_dm_tens, exacc_input
  use exacorr_utils, only     : print_date

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0), MINUS_ONE_HALF=(-0.5D0,0.D0), &
                             ONE_QUARTER=(0.25D0,0.D0), MINUS_ONE_QUARTER=(-0.25D0,0.D0) !Loic_parameter
                             
  private

  public talsh_dm_driver

  contains

    subroutine talsh_dm_eom_driver(exa_input,t2_tensor,r2_tensor,l2_tensor,t1_tensor,r1_tensor,l1_tensor)

!     Routine for obtaining the density matrices in the MO basis

!     input variables
      type(exacc_input), intent(in   ) :: exa_input
      type(talsh_tens_t),intent(inout) :: t2_tensor, l2_tensor ! solution tensors
      type(talsh_tens_t),intent(inout), optional :: t1_tensor, l1_tensor ! solution tensors
      type(talsh_tens_t),intent(inout), optional :: r1_tensor, r2_tensor ! solution tensors  ! Loic_


!     density matrices type, output
      type(talsh_dm_tens) :: gamma
!     tensor dimensions
      integer, dimension(2) :: oo_dims,ov_dims,vo_dims,vv_dims
      integer, dimension(4) :: oooo_dims,vvvv_dims,oovv_dims,vvoo_dims !Loic_
      integer, dimension(4) :: vvvo_dims,ooov_dims,vooo_dims,ovvv_dims !Loic_
      integer, dimension(4) :: vovo_dims                               !Loic_

!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr

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

!************************************
!   construct gamma second tensors  *
!************************************


!**************************************
!     calculate DM blocks in MO basis *
!**************************************

      if (.not.CCD) then
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,t1_tensor,l1_tensor)
      else
        call calc_dm_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir)
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

    subroutine calc_dm_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)

!     input variables
      type(talsh_dm_tens), intent(inout) :: gamma
      type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
      type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

      type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor !Loic_

      integer,             intent(in   ) :: nocc, nvir
!     auxiliary tensors
      type(talsh_tens_t) :: one_tensor
      type(talsh_tens_t) :: tau2_tensor
      type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

      ! Loic_tensor 
      type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
      type(talsh_tens_t) :: MATone_oo_tensor, MATone_vv_tensor    
      type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

      type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor ! Loic !!!! Eviter la redondance
      type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor, oovv_aux_tensor
      type(talsh_tens_t) :: vvoo_aux_vvoo_tensor, vvvo_aux_tensor, ooov_aux_tensor
      type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor, vovo_aux_tensor

!     tensor dimensions
      integer :: one_dims(1)=1
      integer, dimension(2) :: oo_dims, vv_dims
      integer, dimension(4) :: vvoo_dims
      integer, dimension(2) :: oo_dims,ov_dims,vo_dims,vv_dims
      integer, dimension(4) :: oooo_dims,vvvv_dims,oovv_dims,vvoo_dims !Loic_
      integer, dimension(4) :: vvvo_dims,ooov_dims,vooo_dims,ovvv_dims !Loic_
      integer, dimension(4) :: vovo_dims                               !Loic_

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


!********************
! Intermediates     *
!********************

        ! Loic_Taille

        vv_dims = nvir
        ierr=talsh_tensor_construct(MATone_vv_tensor,C8,vv_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(MATtwo_vv_tensor,C8,vv_dims,init_val=ZERO)
        oo_dims = nocc
        ierr=talsh_tensor_construct(MATone_oo_tensor,C8,oo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(MATtwo_oo_tensor,C8,oo_dims,init_val=ZERO)


        vvoo_dims(1:2) = nvir
        vvoo_dims(3:4) = nocc
        ierr=talsh_tensor_construct(tauUP_tensor,C8,vvoo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(tauDown_tensor,C8,vvoo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(tauNo_tensor,C8,vvoo_dims,init_val=ZERO)

        !  Loic_Tenseurs_Intermediaires


       !  MATone_oo
       ierr=talsh_tensor_contract("MATone(i,j)+=R(e,j)*L(j,e)",MATone_oo_tensor,r1_tensor,l1_tensor)
       ierr=talsh_tensor_contract("MATone(i,j)+=R(e,f,j,m)*L(i,m,e,f)",MATone_oo_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
       ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,j)*L(j,e)",MATtwo_oo_tensor,t1_tensor,l1_tensor)
       ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,f,j,m)*L(i,m,e,f)",MATtwo_oo_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
       ierr=talsh_tensor_contract("MATone(a,b)+=R(a,m)*L(m,b)",MATone_vv_tensor,r1_tensor,l1_tensor)
       ierr=talsh_tensor_contract("MATone(a,b)+=R(e,a,m,n)*L(m,n,e,b)",MATone_vv_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
       ierr=talsh_tensor_contract("MATtwo(a,b)+=T(a,m)*L(m,b)",MATtwo_vv_tensor,t1_tensor,l1_tensor)
       ierr=talsh_tensor_contract("MATtwo(a,b)+=T(e,a,m,n)*L(m,n,e,b)",MATtwo_vv_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)

       !  Tau

       ierr=talsh_tensor_add("TauNo(a,b,i,j)+=T(a,b,i,j)",TauNo_tensor,t2_tensor)
       ierr=talsh_tensor_contract("TauNo(a,b,i,j)+=T(a,i)*T(b,j)",TauNo_tensor,t1_tensor,t1_tensor)

       ierr=talsh_tensor_add("TauUp(a,b,i,j)+=T(a,b,i,j)",TauUp_tensor,t2_tensor)
       ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(a,i)*T(b,j)",TauUp_tensor,t1_tensor,t1_tensor)
       ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(b,i)*T(a,j)",TauUp_tensor,t1_tensor,t1_tensor)

       ierr=talsh_tensor_add("TauDown(a,b,i,j)+=T(a,b,i,j)",TauDown_tensor,t2_tensor)
       ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,i)*T(b,j)",TauDown_tensor,t1_tensor,t1_tensor)
       ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,j)*T(b,i)",TauDown_tensor,t1_tensor,t1_tensor)

!********************
!     calculate DMs *
!********************


! rho ij
       ierr=talsh_tensor_add("G(i,j)+=X(i,j)",gamma%prime_oo,MATone_oo_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_contract("X(i,e)+=L(i,m,e,f)*R(f,m)",ov_aux_tensor,l2_tensor,r1_tensor)
       ierr=talsh_tensor_contract("G(i,j)+=X(i,e)*T(e,j)",gamma%prime_oo,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_add("G(i,j)+=X(i,j)",gamma%prime_oo,MATtwo_oo_tensor,scale=MINUS_ONE)

!
!   TAKE CARE no \kronecker_{ij}
!


! rho ab
       ierr=talsh_tensor_add("G(a,b)+=X(a,b)",gamma%prime_vv,MATone_vv_tensor)
       ierr=talsh_tensor_contract("X(n,b)+=L(m,n,e,b)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
       ierr=talsh_tensor_contract("G(a,b)+=X(n,b)*T(a,n)",gamma%prime_vv,ov_aux_tensor,t1_tensor)
       ierr=talsh_tensor_add("G(a,b)+=X(a,b)",gamma%prime_vv,MATtwo_vv_tensor)

! rho ai
       ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,t2_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,MATtwo_oo_tensor,t1_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_contract("X(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t1_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*T(f,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,scale=MINUS_HALF_ONE)
       ierr=talsh_tensor_contract("G(a,i)+=X(n,i)*T(a,n)",gamma%prime_vo,MATone_oo_tensor,t1_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*T(f,i)",gamma%prime_vo,MATone_vv_tensor,t1_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_contract("X(n,f)+=L(m,n,e,f)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(n,f)*TauNo(f,a,i,n)",gamma%prime_vo,vv_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*R(a,e,i,m)",gamma%prime_vo,l1_tensor,r2_tensor)
       ierr=talsh_tensor_contract("X(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*R(f,i)",gamma%prime_vo,vv_aux_tensor,r1_tensor,scale=MINUS_HALF_ONE)
       ierr=talsh_tensor_contract("X(n,i)+=L(m,n,e,f)*T(e,f,m,i)",oo_aux_tensor,l2_tensor,t2_tensor)
       ierr=talsh_tensor_contract("G(a,i)+=X(n,i)*R(a,n)",gamma%prime_vo,oo_aux_tensor,r1_tensor,scale=MINUS_HALF_ONE)
! rho_ia

       ierr=talsh_tensor_add("G(i,a)+=L(i,a)",gamma%prime_ov,l1_tensor)
       ierr=talsh_tensor_contract("G(i,a)+=L(i,m,a,e)*R(e,m)",gamma%prime_ov,l2_tensor,r1_tensor)


!  Symmetrization


       ierr=talsh_tensor_add("G(a,i)+=X(a,i)",gamma%symm_ov,gamma%prime_vo,scale=ONE_HALF)
       ierr=talsh_tensor_add("G(a,i)+=X+(i,a)*K()",gamma%symm_ov,gamma%prime_vo,one_tensor,scale=ONE_HALF)
       ierr=talsh_tensor_add("G(i,a)+=X(i,a)",gamma%symm_vo,gamma%prime_ov,scale=ONE_HALF)
       ierr=talsh_tensor_add("G(i,a)+=X+(a,i)*K()",gamma%symm_vo,gamma%prime_ov,one_tensor,scale=ONE_HALF)


!     debug
      if (.true.) then
        print*, "gprime_oo  = ", tensor_norm2(gamma%prime_oo)
        if (.not.CCD) then
          print*, "gprime_ov  = ", tensor_norm2(gamma%prime_ov)
          print*, "gsymm_ov   = ", tensor_norm2(gamma%symm_ov)
          print*, "gprime_vo  = ", tensor_norm2(gamma%prime_vo)
          print*, "gsymm_vo   = ", tensor_norm2(gamma%symm_vo)
        end if
        print*, "gprime_vv  = ", tensor_norm2(gamma%prime_vv)
        print*, "gsymm_oo   = ", tensor_norm2(gamma%symm_oo)
        print*, "gsymm_vv   = ", tensor_norm2(gamma%symm_vv)
      end if

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
      ierr=talsh_tensor_get_body_access(gamma%symm_oo,body_p2,C8,0,DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_mo_fort(1:nocc,1:nocc) = dm_block

!     vo + ov
      if (.not.CCD) then
        ierr=talsh_tensor_get_body_access(gamma%symm_vo,body_p2,C8,0,DEV_HOST)
        call c_f_pointer(body_p2,dm_block,(/nvir,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_mo_fort(1+nocc:nocc+nvir,1:nocc) = dm_block

        ierr=talsh_tensor_get_body_access(gamma%symm_ov,body_p2,C8,0,DEV_HOST)
        call c_f_pointer(body_p2,dm_block,(/nocc,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_mo_fort(1:nocc,1+nocc:nocc+nvir) = dm_block
      end if

!     vv
      ierr=talsh_tensor_get_body_access(gamma%symm_vv,body_p2,C8,0,DEV_HOST)
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
