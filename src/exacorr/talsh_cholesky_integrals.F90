module talsh_cholesky_integrals
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils
        use talsh_common_routines
        use exacorr_global
        use talsh_ao_to_mo
        use talsh_cholesky_routines

        use, intrinsic:: ISO_C_BINDING

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER=(0.25D0,0.D0), MINUS_ONE_QUARTER=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
        real(8), parameter    :: ONE_QUARTER_R=0.25D0, ONE_HALF_R=0.5D0

        private

        public get_CC_integrals_chol
        public get_cholesky_vec

       contains

       subroutine get_CC_integrals_chol (nocc,nvir,mo_occ,mo_vir,int_t,t_cholesky,print_level)

!        Written by Johann Pototschnig, summer 2019

         implicit none

         integer, intent(in)                  :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)                  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         type(talsh_intg_tens),intent(inout)  :: int_t ! arget integral tensors to be filled and given back to the caller
         real(8), intent(in)                  :: t_cholesky ! cholesky threshold
         integer, intent(in)                  :: print_level

         integer(C_INT)                :: ierr
         integer                       :: nao
         integer(C_INT), dimension(4)  :: oooo_dims,ooov_dims,oovv_dims,ovov_dims,ovvo_dims
         integer(C_INT), dimension(4)  :: vovo_dims,vovv_dims,vvov_dims,vvvv_dims
         type(talsh_tens_t)            :: ovvo_tensor

!        cholesky
         type(talsh_tens_t)              :: chol_tensor, temp_tensor
         type(talsh_tens_t)              :: uoo_tensor, uov_tensor, uvv_tensor
         integer, dimension(3)           :: chol_dims, dims3
         integer(C_INT)                  :: rank
         integer, allocatable         :: mo_list(:)

         type(talsh_tens_t) :: one_chol
         integer(C_INT)     :: scalar_dims(1)
         scalar_dims(1) = 1
         ierr=talsh_tensor_construct(one_chol,C8,scalar_dims(1:0),init_val=ONE)
         
!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get cholesky-vectors and copy them in to a talsh tensor
         call print_date('-Start- Cholesky Decomposition')
         if (print_level.gt.1) write(*,*) "using complete tensor - LV refactored"
         call decompose_cholesky_talsh (chol_tensor, t_cholesky, print_level)
         ierr = talsh_tensor_dimensions(chol_tensor,rank,chol_dims)
         call print_date('-End- Cholesky Decomposition')

!       transform occupied
        dims3=chol_dims
        dims3(1)=nocc
        ierr=talsh_tensor_construct(temp_tensor,C8,dims3,init_val=ZERO)
        dims3(1:2)=nocc
        ierr=talsh_tensor_construct(uoo_tensor,C8,dims3,init_val=ZERO)
        dims3(2)=nvir
        ierr=talsh_tensor_construct(uov_tensor,C8,dims3,init_val=ZERO)

        call ao2mo_ind(chol_tensor,temp_tensor,nocc,mo_occ,11)
        call ao2mo_ind(temp_tensor,uoo_tensor,nocc,mo_occ,21)
        call ao2mo_ind(temp_tensor,uov_tensor,nvir,mo_vir,21)
        ierr=talsh_tensor_init(temp_tensor)
        call ao2mo_ind(chol_tensor,temp_tensor,nocc,mo_occ,12)
        call ao2mo_ind(temp_tensor,uoo_tensor,nocc,mo_occ,22)
        call ao2mo_ind(temp_tensor,uov_tensor,nvir,mo_vir,22)

        if (print_level.gt.2) call print_date(' uoo/uov done ')
        if(print_level > 1) write(*,*) "uoo=",tensor_norm2(uoo_tensor)
        if(print_level > 1) write(*,*) "uov=",tensor_norm2(uov_tensor)

        ierr=talsh_tensor_destruct(temp_tensor)

        dims3=chol_dims
        dims3(1)=nvir
        ierr=talsh_tensor_construct(temp_tensor,C8,dims3,init_val=ZERO)
        dims3(1:2)=nvir
        ierr=talsh_tensor_construct(uvv_tensor,C8,dims3,init_val=ZERO)

        call ao2mo_ind(chol_tensor,temp_tensor,nvir,mo_vir,11)
        call ao2mo_ind(temp_tensor,uvv_tensor,nvir,mo_vir,21)
        ierr=talsh_tensor_init(temp_tensor)
        call ao2mo_ind(chol_tensor,temp_tensor,nvir,mo_vir,12)
        call ao2mo_ind(temp_tensor,uvv_tensor,nvir,mo_vir,22)

        if (print_level.gt.2) call print_date(' uvv done ')
        if(print_level > 1) write(*,*) "uvv=",tensor_norm2(uvv_tensor)

        ierr=talsh_tensor_destruct(temp_tensor)
        ierr=talsh_tensor_destruct(chol_tensor)


        if (print_level.gt.1) call print_date(' ao2mo transformation done ')

!        define dims
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         oooo_dims = nocc
         ooov_dims = nocc
         ooov_dims(4) = nvir
         oovv_dims= nocc
         oovv_dims(3:4) = nvir
         ovov_dims= nocc
         ovov_dims(2) = nvir
         ovov_dims(4) = nvir
         ovvo_dims = nocc
         ovvo_dims(2:3) = nvir
         vovo_dims = nvir
         vovo_dims(2) = nocc
         vovo_dims(4) = nocc
         vovv_dims = nvir
         vovv_dims(2) = nocc
         vvov_dims = nvir
         vvov_dims(3) = nocc
         vvvv_dims = nvir

!        Get antisymmetric oooo tensor
         ierr=talsh_tensor_construct(int_t%oooo,C8,oooo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(temp_tensor,C8,oooo_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L(r,s,u)",temp_tensor,uoo_tensor,uoo_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,int_t%oooo,12)
         ierr=talsh_tensor_destruct(temp_tensor)
         if (print_level.gt.2) call print_date(' oooo done ')

!        Get ooov tensor
         ierr=talsh_tensor_construct(int_t%ooov,C8,ooov_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(temp_tensor,C8,ooov_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L(r,s,u)",temp_tensor,uoo_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,int_t%ooov,12)
         ierr=talsh_tensor_destruct(temp_tensor)
         if (print_level.gt.2) call print_date(' ooov done ')

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
         ierr=talsh_tensor_construct(int_t%vovo,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(temp_tensor,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
         ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L(r,s,u)",temp_tensor,uvv_tensor,uoo_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,int_t%vovo,0)
         ierr=talsh_tensor_destruct(temp_tensor)
         ierr=talsh_tensor_construct(temp_tensor,C8,ovvo_dims,init_val=ZERO)
         !CAREFUL: nasty TALSH bug encountered when using the more logical construct in the next line
         !My (LV) suspicion is that TALSH ignores complex conjugations when placed on the second tensor of the contraction.
         !ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L+(s,r,u)",temp_tensor,uov_tensor,uov_tensor)
         !The following line should be equivalent to the one above but gives a different (and correct) result.
         ierr=talsh_tensor_contract("G(p,q,r,s)+=L+(s,r,u)*K(p,q,u)",temp_tensor,uov_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,ovvo_tensor,0)
         ierr=talsh_tensor_destruct(temp_tensor)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,one_chol,scale=MINUS_ONE)
         ierr=talsh_tensor_destruct(ovvo_tensor)
         if (print_level.gt.2) call print_date(' vovo done ')

!        Clean-up
         ierr=talsh_tensor_destruct(uoo_tensor)

!        Get antisymmetric oovv tensor
         ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(temp_tensor,C8,ovov_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L(r,s,u)",temp_tensor,uov_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,int_t%oovv,12)
         ierr=talsh_tensor_destruct(temp_tensor)
         if (print_level.gt.2) call print_date(' oovv done ')

!        Get vovv tensor
         ierr=talsh_tensor_construct(int_t%vovv,C8,vovv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(temp_tensor,C8,vvov_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L(r,s,u)",temp_tensor,uvv_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,int_t%vovv,34)
         ierr=talsh_tensor_destruct(temp_tensor)
         if (print_level.gt.2) call print_date(' vvov done ')

!        Clean-up
         ierr=talsh_tensor_destruct(uov_tensor)

!        Get vvvv tensor
         ierr=talsh_tensor_construct(int_t%vvvv,C8,vvvv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(temp_tensor,C8,vvvv_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("G(p,q,r,s)+=K(p,q,u)*L(r,s,u)",temp_tensor,uvv_tensor,uvv_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(temp_tensor,int_t%vvvv,12)
         ierr=talsh_tensor_destruct(temp_tensor)
         if (print_level.gt.2) call print_date(' vvvv done ')

!        Clean-up
         ierr=talsh_tensor_destruct(uvv_tensor)

         ierr=talsh_tensor_destruct(one_chol)

       end subroutine get_CC_integrals_chol

   subroutine get_cholesky_vec (nocc,nvir,mo_occ,mo_vir,uov_tensor,uvo_tensor,t_cholesky,print_level)

!        Written by Johann Pototschnig, summer 2019

         implicit none

         integer, intent(in)                  :: nocc,nvir            ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)                  :: mo_occ(:), mo_vir(:) ! the list of occupied and virtual orbitals
         real(8), intent(in)                  :: t_cholesky           ! cholesky threshold
         integer, intent(in)                  :: print_level
         type(talsh_tens_t), intent(out)      :: uov_tensor, uvo_tensor

         integer(C_INT)                  :: ierr
         integer                         :: nao
         type(talsh_tens_t)              :: chol_tensor, temp_tensor
         integer, dimension(3)           :: chol_dims, dims3
         integer(C_INT)                  :: rank
         !integer, allocatable            :: mo_list(:)

         if (print_level.gt.1) call print_date('-Start- cholesky computation')

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get cholesky-vectors and copy them in to a talsh tensor
         call print_date('-Start- Cholesky Decomposition')
         if (print_level.gt.1) write(*,*) "using complete tensor"
         call decompose_cholesky_talsh (chol_tensor, t_cholesky, print_level)
         ierr = talsh_tensor_dimensions(chol_tensor,rank,chol_dims)
         call print_date('-End- Cholesky Decomposition')

!       transform occupied
        dims3=chol_dims
        dims3(1)=nocc
        ierr=talsh_tensor_construct(temp_tensor,C8,dims3,init_val=ZERO)
        dims3(2)=nvir
        ierr=talsh_tensor_construct(uov_tensor,C8,dims3,init_val=ZERO)

        call ao2mo_ind(chol_tensor,temp_tensor,nocc,mo_occ,11)
        call ao2mo_ind(temp_tensor,uov_tensor,nvir,mo_vir,21)
        ierr=talsh_tensor_init(temp_tensor)
        call ao2mo_ind(chol_tensor,temp_tensor,nocc,mo_occ,12)
        call ao2mo_ind(temp_tensor,uov_tensor,nvir,mo_vir,22)

        if (print_level.gt.2) call print_date(' uov done ')
        if(print_level > 1) write(*,*) "uov=",tensor_norm2(uov_tensor)

        ierr=talsh_tensor_destruct(temp_tensor)

        dims3=chol_dims
        dims3(1)=nvir
        ierr=talsh_tensor_construct(temp_tensor,C8,dims3,init_val=ZERO)
        dims3(2)=nocc
        ierr=talsh_tensor_construct(uvo_tensor,C8,dims3,init_val=ZERO)

        call ao2mo_ind(chol_tensor,temp_tensor,nvir,mo_vir,11)
        call ao2mo_ind(temp_tensor,uvo_tensor,nocc,mo_occ,21)
        ierr=talsh_tensor_init(temp_tensor)
        call ao2mo_ind(chol_tensor,temp_tensor,nvir,mo_vir,12)
        call ao2mo_ind(temp_tensor,uvo_tensor,nocc,mo_occ,22)

        if (print_level.gt.2) call print_date(' uvo done ')
        if(print_level > 1) write(*,*) "uvo=",tensor_norm2(uvo_tensor)

        !clean up
        ierr=talsh_tensor_destruct(temp_tensor)
        ierr=talsh_tensor_destruct(chol_tensor)

        if (print_level.gt.1) call print_date('-End- cholesky computation')

   end subroutine get_cholesky_vec

end module talsh_cholesky_integrals
