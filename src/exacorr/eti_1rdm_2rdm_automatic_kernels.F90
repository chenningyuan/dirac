module eti_1rdm_2rdm_automatic_kernels

! Automatic code by G. Fabbro, J. Brandejs, T. Saue, J. Phys. Chem. A 2025, 129, 7360–7372

  use, intrinsic:: ISO_C_BINDING

  use tensor_algebra
  use tapp_top
  use exacorr_tensor_interface
  use exacorr_datatypes
  use exacorr_utils
  use intermediates_cc
  use eti_common_routines
  use exacorr_scaling_parameters
  use exacorr_global

  implicit none

  private
   
  public eti_automatic_dm_1rdm

  contains

    subroutine eti_automatic_dm_1rdm(comm_t,gamma ,nocc,nvir,CCD, &
                        print_level,iterative_level,iterative_scheme,auto_scheme)

!     input variables
      type(eti_dm_tens_t), intent(inout) :: gamma
      integer,             intent(in   ) :: nocc, nvir, print_level
!     auxiliary tensors
      type(tensor_t) :: one_tensor
      type(eti_comm_tens_t), intent(in) :: comm_t
      type(tensor_t)               :: temp_tens
!     tensor dimensions
      integer(C_INT)              :: one_dims(1)=1
!     CCD switch
      logical, intent(in)  :: CCD
!     error code
      integer :: ierr
      integer, intent(in) :: iterative_level,iterative_scheme,auto_scheme
      integer(INTD), dimension(2) :: ov_dims
      type(eti_precision_t) :: eti_precision

      eti_precision%storage_t = comm_t%t2%datatype

!         ov_dims = (/nocc, nvir/)

        ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ONE)    
        if (.not.CCD) then 
         ierr=initialize_tensor_t(gamma%prime_ov,ZERO)
         ierr=initialize_tensor_t(gamma%prime_vo,ZERO)
        endif 

        ierr=initialize_tensor_t(gamma%prime_oo,ZERO)
        ierr=initialize_tensor_t(gamma%prime_vv,ZERO)

       if ((iterative_level.eq.2).and.(iterative_scheme.eq.21)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            !
            else if (auto_scheme.eq.2) then
            ! tenpi
            call eti_tenpi_1rdms_ccsd(nocc,nvir,&
                                        comm_t%l2,&
                                        comm_t%l1,&
                                        comm_t%t1,&
                                        comm_t%t2,&
                                        gamma%prime_oo,&
                                        gamma%prime_ov,&
                                        gamma%prime_vo,&
                                        gamma%prime_vv)
           else if (auto_scheme.eq.3) then 
           call tapp_tenpi_1rdms_ccsd(nocc,nvir,&
                                        comm_t%l2,&
                                        comm_t%l1,&
                                        comm_t%t1,&
                                        comm_t%t2,&
                                        gamma%prime_oo,&
                                        gamma%prime_ov,&
                                        gamma%prime_vo,&
                                        gamma%prime_vv)
            end if
          else if ((iterative_level.eq.2).and.(iterative_scheme.eq.20)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            !
          else if (auto_scheme.eq.2) then
          call eti_tenpi_1rdms_ccd(nocc, nvir, &
                                    comm_t%l2, &
                                    comm_t%t2, &
                                    gamma%prime_oo, &
                                    gamma%prime_vv)
          else if (auto_scheme.eq.3) then 
          call tapp_tenpi_1rdms_ccd(nocc, nvir, &
                                    comm_t%l2, &         
                                    comm_t%t2, &           
                                    gamma%prime_oo, &           
                                    gamma%prime_vv)
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.30)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            ! 
            else if (auto_scheme.eq.2) then
            ! tenpi
              call eti_tenpi_1rdms_ccsdt(nocc,nvir,&
                                                comm_t%l3,&
                                                comm_t%l2,&
                                                comm_t%l1,&
                                                comm_t%t1,&
                                                comm_t%t2,&
                                                comm_t%t3,&
                                                gamma%prime_oo,&
                                                gamma%prime_ov,&
                                                gamma%prime_vo,&
                                                gamma%prime_vv)
           else if (auto_scheme.eq.3) then
           call tapp_tenpi_1rdms_ccsdt(nocc,nvir,&
                                                comm_t%l3,&
                                                comm_t%l2,&
                                                comm_t%l1,&
                                                comm_t%t1,&
                                                comm_t%t2,&
                                                comm_t%t3,&
                                                gamma%prime_oo,&
                                                gamma%prime_ov,&
                                                gamma%prime_vo,&
                                                gamma%prime_vv)
           end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.31)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            ! 
        end if 
          else if ((iterative_level.eq.4).and.(iterative_scheme.eq.40)) then
           if (auto_scheme.eq.1) then 
            ! mb-codegen
            !
           else if (auto_scheme.eq.2) then
            ! tenpi
              call eti_tenpi_1rdms_ccsdtq(nocc,nvir,&
                                                comm_t%l4,&
                                                comm_t%l3,&
                                                comm_t%l2,&
                                                comm_t%l1,&
                                                comm_t%t1,&
                                                comm_t%t2,&
                                                comm_t%t3,&
                                                comm_t%t4,&
                                                gamma%prime_oo,&
                                                gamma%prime_ov,&
                                                gamma%prime_vo,&
                                                gamma%prime_vv) 
           else if (auto_scheme.eq.3) then 
             call tapp_tenpi_1rdms_ccsdtq(nocc,nvir,&
                                                comm_t%l4,&
                                                comm_t%l3,&
                                                comm_t%l2,&
                                                comm_t%l1,&
                                                comm_t%t1,&
                                                comm_t%t2,&
                                                comm_t%t3,&
                                                comm_t%t4,&
                                                gamma%prime_oo,&
                                                gamma%prime_ov,&
                                                gamma%prime_vo,&
                                                gamma%prime_vv)
           end if
       end if

         !     write(*,*) tensor_norm2(gamma%prime_oo)
         !     write(*,*) tensor_norm2(gamma%prime_ov)
         !     write(*,*) tensor_norm2(gamma%prime_vo)
         !     write(*,*) tensor_norm2(gamma%prime_vv)
        if(.not.CCD) then
        !--------------------------------------------------
        !       ov symm: g(i,a) = 1/2 (g'(i,a) + g'(a,i)) |
        !--------------------------------------------------
        ierr=initialize_tensor_t(gamma%symm_ov,ZERO)
        ierr=t_einsum("G(i,a)+=X(i,a)*K()",gamma%symm_ov,gamma%prime_ov,one_tensor,scale=ONE_HALF)
        ierr=t_einsum("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_tensor,scale=ONE_HALF)
        
        !--------------------------------------------------
        !       vo symm: g(a,i) = symmetrized c.c. g(i,a)
        !--------------------------------------------------
         ierr=initialize_tensor_t(gamma%symm_vo,ZERO)
         ierr=t_einsum("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%symm_ov,one_tensor)
        end if
        !------------------------------------------------
        !     vv symm: g(b,a) = 1/2 (g'(b,a) + g'(a,b)) |
        !------------------------------------------------
         ierr=initialize_tensor_t(gamma%symm_vv,ZERO)
         ierr=t_einsum("G(b,a)+=X(b,a)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)
         ierr=t_einsum("G(b,a)+=X+(a,b)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)
        !------------------------------------------------
        !     oo symm: g(j,i) = 1/2 (g'(j,i) + g'(i,j)) |
        !------------------------------------------------
         ierr=initialize_tensor_t(gamma%symm_oo,ZERO)
         ierr=t_einsum("G(j,i)+=X(j,i)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)
         ierr=t_einsum("G(j,i)+=X+(i,j)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)

         ! write(*,*) 'Symmetrized tensors OO', tensor_norm2(gamma%symm_oo)
         ! write(*,*) 'Symmetrized tensors OV', tensor_norm2(gamma%symm_ov)
         ! write(*,*) 'Symmetrized tensors VO', tensor_norm2(gamma%symm_vo)
         ! write(*,*) 'Symmetrized tensors VV', tensor_norm2(gamma%symm_vv)

    end subroutine eti_automatic_dm_1rdm


end module 


