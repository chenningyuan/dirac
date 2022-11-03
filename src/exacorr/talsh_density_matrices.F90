module talsh_gradient

! This module contains routines to calculate density matrices
! See Shee, A.; Visscher, L.; Saue, T. J. Chem. Phys. 2016, 145, 184107 (ground-state 1RDM)
! and Halbert, L and Gomes, ASP, in preparation (ground-state 2RDM, to come EOMCC 1RDMs) 
! Written by Stan Papadopoulos, February 2019
! Refactored by Andre Gomes and Loic Halbert summer 2022, to move tensor operations elsewhere so that this file contains common/general manipulation operations,

  use talsh
  use tensor_algebra
  use talsh_common_routines, only : tensor_norm2, print_tensor_elements
  use exacorr_datatypes, only : talsh_dm_tens, exacc_input
  use exacorr_utils, only     : print_date
  use talsh_1rdm_2rdm_ccsd_kernels
  use, intrinsic:: ISO_C_BINDING

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0), MINUS_ONE_HALF=(-0.5D0,0.D0), &
                             MINUS_ONE_QUARTER=(-0.25D0,0.D0), ONE_QUARTER=(0.25D0,0.D0), &
                             MINUS_ONE_EIGHT=(-0.125D0,0.D0)

  private

  public talsh_dm_driver
  public talsh_dm_1rdm_form_mo_fort_from_tensors
  public talsh_dm_1rdm_form_ao_diracformat_from_dm_mo
  public talsh_dm_initialize_dmtens
  public talsh_dm_cleanup_dmtens
  public talsh_dm_store_mo

  contains

    subroutine talsh_dm_energy_from_1rdm_2rdm(recalculated_energy,exa_input,gamma,int_t)
      use exacorr_datatypes

      complex(8), intent (inout) :: recalculated_energy
      type(exacc_input), intent(in   ) :: exa_input
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_dm_tens), intent(inout) :: gamma
      type(C_PTR)         :: res_body_p
      type(talsh_tens_t)  :: res_tensor
      integer(INTD)       :: res_dims(1)
      complex(8), pointer :: res_tens(:)
      integer(INTD) :: ierr

!     setup for calculate the wavefun ground-state energy from the 1RDM and 2RDM
!     this calculations will always be done as a cross-check, if the generation of the ground-state 2RDM is requested
      res_dims(1) = 1
      ierr=talsh_tensor_construct(res_tensor,C8,res_dims(1:0),init_val=ZERO)
      ierr=talsh_tensor_get_body_access(res_tensor,res_body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(res_body_p,res_tens,res_dims)

!  Symmetry used
      ierr=talsh_tensor_init(res_tensor)
      ierr=talsh_tensor_contract("E()+=G(i,a)*F+(i,a)",res_tensor,gamma%symm_ov,int_t%fov)
      if (exa_input%print_level.ge.2) print*, "Final ov symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(a,i)*F(i,a)",res_tensor,gamma%symm_vo,int_t%fov)
      if (exa_input%print_level.ge.2) print*, "Final vo symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(a,b)*F(b,a)",res_tensor,gamma%symm_vv,int_t%fvv)
      if (exa_input%print_level.ge.2) print*, "Final vv symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(i,j)*F(j,i)",res_tensor,gamma%symm_oo,int_t%foo)
      if (exa_input%print_level.ge.2) print*, "Final oo symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(a,b,i,j)*V(i,j,a,b)",res_tensor,gamma%symm_vvoo,int_t%oovv,scale=ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final vvoo symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(i,j,a,b)*V+(i,j,a,b)",res_tensor,gamma%symm_oovv,int_t%oovv,scale=ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final oovv symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(i,j,k,a)*V+(i,j,k,a)",res_tensor,gamma%symm_ooov,int_t%ooov,scale=ONE_HALF)
      if (exa_input%print_level.ge.2) print*, "Final ooov symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(i,a,j,k)*V(j,k,i,a)",res_tensor,gamma%symm_ovoo,int_t%ooov,scale=ONE_HALF)
      if (exa_input%print_level.ge.2) print*, "Final ovoo symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(a,b,c,i)*V(c,i,a,b)",res_tensor,gamma%symm_vvvo,int_t%vovv,scale=ONE_HALF)
      if (exa_input%print_level.ge.2) print*, "Final vvvo symm : ",res_tens(1)

!      ierr=talsh_tensor_contract("E()+=G(a,i,b,c)*V(b,c,a,i)",result_tensor,gamma%prime_vovv,int_t%vvvo,scale=ONE_QUARTER)
      ierr=talsh_tensor_contract("E()+=G(a,i,b,c)*V+(a,i,b,c)",res_tensor,gamma%symm_vovv,int_t%vovv,scale=ONE_HALF)

      ierr=talsh_tensor_contract("E()+=G(i,a,j,b)*V(b,j,a,i)",res_tensor,gamma%symm_ovov,int_t%vovo)
      if (exa_input%print_level.ge.2) print*, "Final ovov symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(a,b,c,d)*V(c,d,a,b)",res_tensor,gamma%symm_vvvv,int_t%vvvv,scale=ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final vvvv symm : ",res_tens(1)

      ierr=talsh_tensor_contract("E()+=G(i,j,k,l)*V(k,l,i,j)",res_tensor,gamma%symm_oooo,int_t%oooo,scale=ONE_QUARTER)
      if (exa_input%print_level.ge.2) print*, "Final oooo symm : ",res_tens(1)

      recalculated_energy = res_tens(1) 

      ierr=talsh_tensor_destruct(res_tensor)
    
    end subroutine

    subroutine talsh_dm_driver(exa_input,t2_tensor,l2_tensor,t1_tensor,l1_tensor, int_t)
      use exacorr_datatypes

!     Routine for obtaining the density matrices in the MO basis

!     input variables
      type(exacc_input), intent(in   ) :: exa_input
      type(talsh_tens_t),intent(inout) :: t2_tensor, l2_tensor ! solution tensors
      type(talsh_tens_t),intent(inout), optional :: t1_tensor, l1_tensor ! solution tensors
      type(talsh_intg_tens), intent(inout), optional :: int_t
!     density matrices type, output
      type(talsh_dm_tens) :: gamma
!     CCD switch
      logical :: CCD
      character(len=4) :: wavefun

! aspg: temporary solution for handling result_tensor (need to move the double-check on the energies outside this routine)
      type(C_PTR)         :: res_body_p
!      type(talsh_tens_t)  :: res_tensor
!      integer(INTD)       :: res_dims(1)
      complex(8), pointer :: res_tens(:)
      complex(8) :: recalculated_energy
! aspg: temporary solution for handling result_tensor (need to move the double-check on the energies outside this routine)

!     error code
      integer(INTD) :: ierr

      CCD = exa_input%CCD
      wavefun = "CCSD"

!      if (exa_input%do_2rdm_gs) then
!!     setup for calculate the wavefun ground-state energy from the 1RDM and 2RDM
!!     this calculations will always be done as a cross-check, if the generation of the ground-state 2RDM is requested
!         res_dims(1) = 1
!         ierr=talsh_tensor_construct(res_tensor,C8,res_dims(1:0),init_val=ZERO)
!         ierr=talsh_tensor_get_body_access(res_tensor,res_body_p,C8,int(0,C_INT),DEV_HOST)
!         call c_f_pointer(res_body_p,res_tens,res_dims)
!      end if
 
!************************************
!     construct gamma prime tensors *
!************************************
      call talsh_dm_initialize_dmtens(exa_input,gamma)

!*****************************************
!     calculate 1-RDM blocks in MO basis *
!*****************************************
      if (exa_input%do_1rdm_gs) then
         if (.not.CCD) then
           call talsh_calc_dm_1rdm_ccsd_GS_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,&
                exa_input%print_level,t1_tensor,l1_tensor)
         else
           call talsh_calc_dm_1rdm_ccsd_GS_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,&
                exa_input%print_level)
         end if
         call print_date("Calculated 1-RDM density matrices in MO basis")
      end if

!*****************************************
!     calculate 2-RDM blocks in MO basis *
!*****************************************
      if (exa_input%do_2rdm_gs) then
         if (.not.CCD) then
           call talsh_calc_dm_2rdm_ccsd_GS_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,&
                exa_input%print_level,t1_tensor,l1_tensor)
         else
           call talsh_calc_dm_2rdm_ccsd_GS_mo(gamma,t2_tensor,l2_tensor,exa_input%nocc,exa_input%nvir,&
                exa_input%print_level)
         end if

         call print_date("Calculated 2-RDM density matrices in MO basis")

         call talsh_dm_energy_from_1rdm_2rdm(recalculated_energy,exa_input,gamma,int_t)
         write (*,'(A)') ""
         write (*,'(A,F20.15)') "  Final "//wavefun//" energy (from 1RDM+2RDM):",real(recalculated_energy)
         write (*,'(A)') ""
         
      end if

!***********************
!     write DM to file *
!***********************
      call talsh_dm_store_mo(exa_input,gamma,wavefun)
      call talsh_kudm_store_mo(exa_input,gamma,wavefun) !also storing KU format!
!     cleanup
      call talsh_dm_cleanup_dmtens(exa_input,gamma)

    end subroutine talsh_dm_driver

    subroutine talsh_dm_cleanup_dmtens(exa_input,dmat)
      type(exacc_input),   intent(in   ) :: exa_input
      type(talsh_dm_tens), intent(inout) :: dmat
!     error code
      integer(INTD) :: ierr

      if (exa_input%do_1rdm_gs) then
         ierr=talsh_tensor_destruct(dmat%prime_oo)
         if (.not.exa_input%ccd) then
           ierr=talsh_tensor_destruct(dmat%prime_ov)
           ierr=talsh_tensor_destruct(dmat%symm_ov)
           ierr=talsh_tensor_destruct(dmat%prime_vo)
           ierr=talsh_tensor_destruct(dmat%symm_vo)
         end if
         ierr=talsh_tensor_destruct(dmat%prime_vv)
         ierr=talsh_tensor_destruct(dmat%symm_oo)
         ierr=talsh_tensor_destruct(dmat%symm_vv)
         call print_date("Destroyed 1RDM storage")
      end if

      if (exa_input%do_2rdm_gs) then
         ierr=talsh_tensor_destruct(dmat%prime_oooo)
         ierr=talsh_tensor_destruct(dmat%prime_vooo)
         ierr=talsh_tensor_destruct(dmat%prime_ooov)
         ierr=talsh_tensor_destruct(dmat%prime_ovoo)
         ierr=talsh_tensor_destruct(dmat%prime_oovv)
         ierr=talsh_tensor_destruct(dmat%prime_vvoo)
         ierr=talsh_tensor_destruct(dmat%prime_ovvo)
         ierr=talsh_tensor_destruct(dmat%prime_voov)
         ierr=talsh_tensor_destruct(dmat%prime_ovov)
         ierr=talsh_tensor_destruct(dmat%prime_vvvo)
         ierr=talsh_tensor_destruct(dmat%prime_vovv)
         ierr=talsh_tensor_destruct(dmat%prime_vvvv)

         ierr=talsh_tensor_destruct(dmat%symm_oooo)
         ierr=talsh_tensor_destruct(dmat%symm_vooo)
         ierr=talsh_tensor_destruct(dmat%symm_ooov)
         ierr=talsh_tensor_destruct(dmat%symm_ovoo)
         ierr=talsh_tensor_destruct(dmat%symm_oovv)
         ierr=talsh_tensor_destruct(dmat%symm_vvoo)
         ierr=talsh_tensor_destruct(dmat%symm_ovvo)
         ierr=talsh_tensor_destruct(dmat%symm_voov)
         ierr=talsh_tensor_destruct(dmat%symm_ovov)
         ierr=talsh_tensor_destruct(dmat%symm_vvvo)
         ierr=talsh_tensor_destruct(dmat%symm_vovv)
         ierr=talsh_tensor_destruct(dmat%symm_vvvv)
         call print_date("Destroyed 2RDM storage")
      end if

      if (ierr.ne.TALSH_SUCCESS) then
        print*," ccdriver: error in evaluation of density matrices, code=",ierr
        return
      end if
    end subroutine

    subroutine talsh_dm_initialize_dmtens(exa_input,dmat)
!     Routine for obtaining the density matrices in the MO basis
!     input variables
      type(exacc_input), intent(in   ) :: exa_input
!     density matrices type, output
      type(talsh_dm_tens), intent(inout) :: dmat 
!     tensor dimensions
      integer(INTD), dimension(2) :: oo_dims,ov_dims,vo_dims,vv_dims
!     CCD switch
      logical :: CCD
!     error code
      integer(INTD) :: ierr
      integer :: nocc, nvir

      nocc = exa_input%nocc
      nvir = exa_input%nvir

!*********************************************
!     construct gamma prime and symm tensors *
!*********************************************

       if (exa_input%do_1rdm_gs) then
          oo_dims = nocc
          ierr=talsh_tensor_construct(dmat%prime_oo,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_oo,C8,oo_dims,init_val=ZERO)
          if (.not.exa_input%ccd) then
            ov_dims(1) = nocc
            ov_dims(2) = nvir
            ierr=talsh_tensor_construct(dmat%prime_ov,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(dmat%symm_ov,C8,ov_dims,init_val=ZERO)
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=talsh_tensor_construct(dmat%prime_vo,C8,vo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(dmat%symm_vo,C8,vo_dims,init_val=ZERO)
          end if
          vv_dims = nvir
          ierr=talsh_tensor_construct(dmat%prime_vv,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_vv,C8,vv_dims,init_val=ZERO)
          call print_date("Created 1RDM storage")
       end if

       if (exa_input%do_2rdm_gs) then
          ierr=talsh_tensor_construct(dmat%prime_oooo,C8,(/nocc,nocc,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_vooo,C8,(/nvir,nocc,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_ooov,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_ovoo,C8,(/nocc,nvir,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_oovv,C8,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_vvoo,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_ovvo,C8,(/nocc,nvir,nvir,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_voov,C8,(/nvir,nocc,nocc,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_ovov,C8,(/nocc,nvir,nocc,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_vvvo,C8,(/nvir,nvir,nvir,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_vovv,C8,(/nvir,nocc,nvir,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%prime_vvvv,C8,(/nvir,nvir,nvir,nvir/),init_val=ZERO)

          ierr=talsh_tensor_construct(dmat%symm_oooo, C8,(/nocc,nocc,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_vooo, C8,(/nvir,nocc,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_ooov, C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_ovoo, C8,(/nocc,nvir,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_oovv, C8,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_vvoo, C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_ovvo, C8,(/nocc,nvir,nvir,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_voov, C8,(/nvir,nocc,nocc,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_ovov, C8,(/nocc,nvir,nocc,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_vvvo, C8,(/nvir,nvir,nvir,nocc/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_vovv, C8,(/nvir,nocc,nvir,nvir/),init_val=ZERO)
          ierr=talsh_tensor_construct(dmat%symm_vvvv, C8,(/nvir,nvir,nvir,nvir/),init_val=ZERO)
          call print_date("Created 2RDM storage")
       end if

       if (ierr.ne.TALSH_SUCCESS) then
         print*," ccdriver: error in evaluation of density matrices, code=",ierr
         return
       end if

    end subroutine

    subroutine talsh_dm_store_mo (exa_input,gamma,wavefun)

      use exacorr_global, only : get_nsolutions

!     input variables
      type(exacc_input),   intent(in   )  :: exa_input
      type(talsh_dm_tens), intent(inout)  :: gamma
      character*(*), intent(in)           :: wavefun
!     matrix dimensions
      integer :: nocc, nvir, nkr_occ, nkr_vir, nesh
!     body pointers
      type(C_PTR)                     :: body_p2
      complex(8), pointer, contiguous :: dm_block(:,:)
      complex(8), allocatable         :: dm_1rdm_mo_fort(:,:)
      real(8),allocatable             :: dm_1rdm_diracformat(:,:,:)
      complex(8), allocatable         :: dm_2rdm_mo_fort(:,:,:,:)
      real(8),allocatable             :: dm_2rdm_diracformat(:,:,:)
!     complete mo list
      integer, allocatable :: mo_list(:)
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr
!     loop variables
      integer :: i, j, norbt
      integer :: i_unbar, j_unbar, j_bar
! debug
      logical :: debug=.false.
      integer :: k

      CCD = exa_input%CCD

!***************
!     set dims *
!***************
 
      nocc    = exa_input%nocc
      nvir    = exa_input%nvir
      nkr_occ = exa_input%nkr_occ
      nkr_vir = exa_input%nkr_vir

! we always setup a 1RDM
! set up storage for 1rdm in MO and AO bases, in fortran storage
      allocate (dm_1rdm_mo_fort(nocc+nvir,nocc+nvir))
      dm_1rdm_mo_fort = ZERO

      norbt = nkr_occ + nkr_vir
      allocate(mo_list(norbt))
      mo_list(1:nkr_occ)                 = exa_input%mokr_occ
      mo_list(1+nkr_occ:nkr_occ+nkr_vir) = exa_input%mokr_vir

      nesh = get_nsolutions() / 2
      allocate (dm_1rdm_diracformat(nesh,nesh,4))
      dm_1rdm_diracformat = real(ZERO)

! form MO basis density matrix, then AO basis density matrix, and finally save AO basis density matrix to file
      call talsh_dm_1rdm_form_mo_fort_from_tensors(dm_1rdm_mo_fort, gamma, nocc, nvir, CCD)
      call talsh_dm_1rdm_form_ao_diracformat_from_dm_mo(dm_1rdm_diracformat, dm_1rdm_mo_fort, mo_list, norbt)
      call store_cc_density (wavefun, dm_1rdm_diracformat, nesh, 4)

      deallocate (dm_1rdm_diracformat)
      deallocate (dm_1rdm_mo_fort)

! in case the 2RDM is requested, we store it in fortran storage as well
      if (exa_input%do_2rdm_gs) then
         allocate (dm_2rdm_mo_fort(nocc+nvir,nocc+nvir,nocc+nvir,nocc+nvir))
         dm_2rdm_mo_fort = ZERO

! form MO basis density matrix, then AO basis density matrix, and finally save AO basis density matrix to file
         call talsh_dm_2rdm_form_mo_fort_from_tensors(dm_2rdm_mo_fort, gamma, nocc, nvir, CCD)
! aspg: TODO is prepare transformation to AO matrix. also, one needs to define the format for that
! aspg: while we do noy yet have a definitive writing to AO basis, we dump the MO basis 2RDM to file using 
!       either a FCDUMP-style format (p q r s Re(Gamma(p,q,r,s) Im(Gamma(p,q,r,s)) or put the whole matrix
!       in a single record and write that to file in binary format 
         call dump_dm_2rdm_to_file_as_binary(4242,"CCSD2RDM",dm_2rdm_mo_fort)
      end if

      deallocate (mo_list)

    end subroutine talsh_dm_store_mo

    subroutine talsh_dm_1rdm_form_mo_fort_from_tensors(dm_1rdm_mo_fort,gamma,nocc,nvir,CCD)
!     input variables
      type(talsh_dm_tens), intent(inout)     :: gamma 
      complex(8), allocatable, intent(inout) :: dm_1rdm_mo_fort(:,:)
!     matrix dimensions
      integer, intent(in) :: nocc, nvir
      logical, intent(in) :: CCD 

!     body pointers
      type(C_PTR)                     :: body_p2
      complex(8), pointer, contiguous :: dm_block(:,:)
!     error code
      integer :: ierr
!     loop variables
      integer :: i, j
! debug
      logical :: debug=.false.
      integer :: k

!***************************************
!     copy DM blocks into one large DM *
!***************************************
!     -----------
!     | OO | OV |
!     -----------
!     | VO | VV |
!     -----------

!     oo
      ierr=talsh_tensor_get_body_access(gamma%symm_oo,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_1rdm_mo_fort(1:nocc,1:nocc) = dm_block

!     vo + ov
      if (.not.CCD) then
        ierr=talsh_tensor_get_body_access(gamma%symm_vo,body_p2,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p2,dm_block,(/nvir,nocc/)) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_1rdm_mo_fort(1+nocc:nocc+nvir,1:nocc) = dm_block

        ierr=talsh_tensor_get_body_access(gamma%symm_ov,body_p2,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p2,dm_block,(/nocc,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
        dm_1rdm_mo_fort(1:nocc,1+nocc:nocc+nvir) = dm_block
      end if

!     vv
      ierr=talsh_tensor_get_body_access(gamma%symm_vv,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nvir,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
      dm_1rdm_mo_fort(1+nocc:nocc+nvir,1+nocc:nocc+nvir) = dm_block

      if (debug) then
        print *,'<<<<< DM, complete, in MO basis and fortran storage'
        k = 1
        do i=1,(nocc+nvir)
           do j=1,(nocc+nvir)
              write(*,'(3X,I4,2x,2F28.10,2x,4I4)') k, dm_1rdm_mo_fort(i,j),i,j
              k = k + 1
           end do
         end do
         print *,'>>>>> DM, complete, in MO basis and fortran storage'
      end if

    end subroutine

    subroutine talsh_dm_2rdm_form_mo_fort_from_tensors(dm_2rdm_mo_fort,gamma,nocc,nvir,CCD)
!     input variables
      type(talsh_dm_tens), intent(inout)     :: gamma
      complex(8), allocatable, intent(inout) :: dm_2rdm_mo_fort(:,:,:,:)
!     matrix dimensions
      integer, intent(in) :: nocc, nvir
      logical, intent(in) :: CCD

!     body pointers
      type(C_PTR)                     :: body_p2
      complex(8), pointer, contiguous :: dm_block(:,:,:,:)
!     error code
      integer :: ierr
!     loop variables
      integer :: i, j
! debug
      logical :: debug=.false.
      integer :: k

!***************************************
!     copy DM blocks into one large DM *
!***************************************
!     -----------------------------
!     | OOOO | OOOV | OOVO | OOVV | 
!     | OVOO | OVOV | OVVO | OVVV |
!     | VOOO | VOOV | VOVO | VOVV |
!     | VVOO | VVOV | VVVO | VVVV |
!     -----------------------------

!     oooo
      ierr=talsh_tensor_get_body_access(gamma%symm_oooo,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nocc,nocc,nocc/)) ! to use <dm_mo_block> as a regular Fortran 4d array
      dm_2rdm_mo_fort(1:nocc,1:nocc,1:nocc,1:nocc) = dm_block

!     ooov
      ierr=talsh_tensor_get_body_access(gamma%symm_ooov,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nocc,nocc,nocc,nvir/)) ! to use <dm_mo_block> as a regular Fortran 4d array
      dm_2rdm_mo_fort(1:nocc,1:nocc,1:nocc,1:nvir) = dm_block

!     vooo
      ierr=talsh_tensor_get_body_access(gamma%symm_vooo,body_p2,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p2,dm_block,(/nvir,nocc,nocc,nocc/)) ! to use <dm_mo_block> as a regular Fortran 4d array
      dm_2rdm_mo_fort(1:nvir,1:nocc,1:nocc,1:nocc) = dm_block


   end subroutine

   subroutine talsh_dm_1rdm_form_ao_diracformat_from_dm_mo(dm_1rdm_diracformat, dm_1rdm_mo_fort, mo_list, norbt)
      complex(8), allocatable, intent(in) :: dm_1rdm_mo_fort(:,:)
      integer, allocatable, intent(in) :: mo_list(:)
      real(8),allocatable, intent(inout)  :: dm_1rdm_diracformat(:,:,:)
      integer, intent(in) :: norbt
!     loop variables
      integer :: i, j
      integer :: i_unbar, j_unbar, j_bar

!*****************************************************************************************************************************
!     Determine where an element of the MO density matrix  should go in the full list of MOs 
!     (NB: DIRAC assumes restricted densities and uses Kramers pairs) 
!*****************************************************************************************************************************

      do j = 1, norbt
        j_unbar = 2*j-1
        j_bar   = 2*j
        do i = 1, norbt
          ! We copy the unbar-unbar part of the density matrix, as well as the unbar-bar part
          ! This follows from eq. 27 of Shee, Visscher, Saue JCP 145 184107, 2016
          ! it would be good to test whether the Kramers-restricted assumption holds
          ! Or even to average over the two density matrices and enforce a Kramers-restricted density
          i_unbar = 2*i-1
          dm_1rdm_diracformat(mo_list(i),mo_list(j),1) = real(dm_1rdm_mo_fort(i_unbar,j_unbar))
          dm_1rdm_diracformat(mo_list(i),mo_list(j),2) = imag(dm_1rdm_mo_fort(i_unbar,j_unbar))
          dm_1rdm_diracformat(mo_list(i),mo_list(j),3) = real(dm_1rdm_mo_fort(i_unbar,j_bar))
          dm_1rdm_diracformat(mo_list(i),mo_list(j),4) = imag(dm_1rdm_mo_fort(i_unbar,j_bar))
        end do
      end do
   end subroutine

!ROUTINES FOR STORING CCDENS IN KUFORMAT ON EXACC FILE
    subroutine talsh_kudm_store_mo (exa_input,gamma,wavefun)

      use exacorr_global, only : get_nsolutions

!     input variables
      type(exacc_input),   intent(in   )  :: exa_input
      type(talsh_dm_tens), intent(inout)  :: gamma
      character*(*), intent(in)           :: wavefun
!     matrix dimensions
      integer :: nocc, nvir, nesh
!     body pointers
      type(C_PTR)                     :: body_p2
      complex(8), pointer, contiguous :: dm_block(:,:)
      complex(8), allocatable         :: dm_1rdm_mo_fort(:,:)
      real(8),allocatable             :: dm_1rdm_storeformat(:,:,:)
!     complete mo list
      integer, allocatable :: mo_list(:)
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr
!     loop variables
      integer :: i, j, norbt
! debug
      logical :: debug=.false.
      integer :: k

      CCD = exa_input%CCD

!***************
!     set dims *
!***************

      nocc    = exa_input%nocc
      nvir    = exa_input%nvir

! set up storage for density matrices in MO and AO bases, in fortran storage
      allocate (dm_1rdm_mo_fort(nocc+nvir,nocc+nvir))
      dm_1rdm_mo_fort = ZERO

      norbt = nocc + nvir
      allocate(mo_list(norbt))
      mo_list(1:nocc) = exa_input%mo_occ
      mo_list(nocc+1:norbt) = exa_input%mo_vir

      nesh = get_nsolutions()
      allocate (dm_1rdm_storeformat(nesh,nesh,2))
      dm_1rdm_storeformat = real(ZERO)

! form MO basis density matrix
      call talsh_dm_1rdm_form_mo_fort_from_tensors(dm_1rdm_mo_fort, gamma, nocc, nvir, CCD)
! form AO basis density matrix
      call talsh_kudm_1rdm_form_ao_storeformat_from_dm_mo(dm_1rdm_storeformat, dm_1rdm_mo_fort, mo_list, norbt)
! save AOm basis density matrix to file
      call talsh_store_kucc_density (wavefun, dm_1rdm_storeformat, nesh, mo_list, nocc)

      deallocate (dm_1rdm_storeformat)
      deallocate (mo_list)
      deallocate (dm_1rdm_mo_fort)

    end subroutine talsh_kudm_store_mo

    subroutine talsh_kudm_1rdm_form_ao_storeformat_from_dm_mo(dm_1rdm_storeformat,dm_1rdm_mo_fort, mo_list, norbt)
      complex(8), allocatable, intent(in) :: dm_1rdm_mo_fort(:,:)
      integer, allocatable, intent(in) :: mo_list(:)
      real(8),allocatable, intent(inout)  :: dm_1rdm_storeformat(:,:,:)
      integer, intent(in) :: norbt
!     loop variables
      integer :: i, j

!*****************************************************************************************************************************
!     Determine where an element of the MO density matrix  should go in the full
!     list of MOs
!*****************************************************************************************************************************

      do j = 1, norbt
        do i = 1, norbt
          ! Extracting real and imaginary parts of KU-DM
          dm_1rdm_storeformat(mo_list(i),mo_list(j),1) = real(dm_1rdm_mo_fort(i,j))
          dm_1rdm_storeformat(mo_list(i),mo_list(j),2) = dimag(dm_1rdm_mo_fort(i,j))
        end do
      end do
    end subroutine

    subroutine talsh_store_kucc_density(wavefun, dm_1rdm_storeformat, nmo, mo_list, nocc)
      use labeled_storage
      !use exacorr_datatypes
      implicit none

      character(len=4) :: wavefun
      real(8), allocatable, intent(inout) :: dm_1rdm_storeformat(:,:,:)
      integer, intent(in)                 :: nmo, nocc
      integer, allocatable, intent(in) :: mo_list(:)
      type(file_info_t), save             :: exacc_file     
      
      call print_date("Writing CCDENS data to EXACC file")
 
      if (exacc_file%status == -1) then
         write(*,*) " Opening HDF5 EXACC file"
         exacc_file%type = 2
         exacc_file%name = "EXACC.h5"
         exacc_file%status = 0
      end if

      call lab_write (exacc_file,'nmo',nmo)
      call lab_write (exacc_file,'CCDENS_RE',rdata=reshape(transpose(dm_1rdm_storeformat(1:nmo,1:nmo,1)), shape=(/nmo**2/) ) )
      call lab_write (exacc_file,'CCDENS_IM',rdata=reshape(transpose(dm_1rdm_storeformat(1:nmo,1:nmo,2)), shape=(/nmo**2/) ) )
      
    end subroutine

    subroutine dump_dm_1rdm_to_file_as_binary(fileno,filename,matrix)
        integer, intent(in) :: fileno
        character(len=*), intent(in) :: filename
        complex(kind=8), intent(in) :: matrix(:,:)
        integer :: irec

        inquire (iolength = irec) matrix
        open (unit=fileno, file=filename, form="unformatted", access="direct", recl=irec)
        write(unit=fileno, rec=1) matrix
        close(unit=fileno)
        print *,''
        print *,' Wrote 1RDM(no+nv,no+nv) to binary file ',filename
        print *,''
    end subroutine

    subroutine dump_dm_2rdm_to_file_as_binary(fileno,filename,matrix)
        integer, intent(in) :: fileno
        character(len=*), intent(in) :: filename
        complex(kind=8), intent(in) :: matrix(:,:,:,:)
        integer :: irec

        inquire (iolength = irec) matrix
        open (unit=fileno, file=filename, form="unformatted", access="direct", recl=irec)
        write(unit=fileno, rec=1) matrix
        close(unit=fileno)
        print *,''
        print *,' Wrote 2RDM(no+nv,no+nv,no+nv,no+nv) to binary file ',filename
        print *,''
    end subroutine

end module talsh_gradient
