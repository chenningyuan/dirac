!This module contains the routines needed to calculate ae coupled cluster wave function
       module exacorr_cc


        implicit none
        private
        public cc_driver

       contains
!-------------------------------------------------------------------------------------------
        subroutine cc_driver (nocc,nvirt,mo_occ,mo_virt)

!This subroutine drives the ccd energy calculation
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
 
!        Written by Lucas Visscher, winter 2016/2017
 
         use tensor_algebra
         use talsh
         use exacorr_ao_to_mo
         use dirac_MO
         use dirac_exatensor

         implicit none
         integer, intent(in) :: nocc,nvirt   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in) :: mo_occ(:),mo_virt(:) ! the list of occupied and virtual orbitals

         integer   :: ierr, i
         complex(8):: zero=(0.D0,0.D0),one=(1.D0,0.D0),minus_one=(-1.D0,0.D0)
         complex(8), pointer, contiguous:: ao_tens(:,:,:,:),mo_tens(:,:,:,:)
         type(C_PTR):: body_p
         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: host_arg_max

         integer             :: nshells, number_aoshells   ! number of AO shells
         type(basis_func_info_t), allocatable :: shell(:)  ! array with basis function information
         integer             :: nao, number_aos

         type(talsh_tens_t):: aoint_tensor
         type(talsh_tens_t):: moint_tensor
         type(talsh_tens_t):: vvoo_tensor
         type(talsh_tens_t):: t2_tensor
         type(talsh_tens_t):: one_tensor
         type(talsh_tens_t):: mone_tensor
         integer:: aoint_dims(1:4),moint_dims(1:4),vvoo_dims(1:4),t2_dims(1:4),scalar_dims(1)
         integer:: nmo(4)
         integer, allocatable :: mo_list(:)

!        Make some noise so that we know we are here
         print*, " Arrived in cc_driver routine"

!        Initialize libraries
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         print*, " Initialized talsh library"
         write(*,'(" Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
         call interest_initialize
         print*, " Initialized interest library"

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)
         nshells = number_aoshells()
         allocate(shell(nshells))
         print*, " Allocated ",nshells," shell variables"
         call initialize_shells(shell)
         do i = 1, nshells
            if (compress_unc_shell(shell(i)) == 1) call quit ('This version can not handle contracted functions')
         enddo
         print*, " Retrieved basis set information from DIRAC common blocks"
         
!        Get tensor with AO integrals
         aoint_dims      = nao
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=zero)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens,nshells,shell)

!        Get vvoo tensor (assuming closed shells for the moment)
         nmo(1) = nocc
         nmo(2) = nvirt
         nmo(3) = nocc
         nmo(4) = nvirt
         moint_dims(1) = nocc * 2
         moint_dims(2) = nvirt * 2
         moint_dims(3) = nocc * 2
         moint_dims(4) = nvirt * 2
         mo_list(1:nocc)                        = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvirt)             = mo_virt(1:nvirt)
         mo_list(nocc+nvirt+1:2*nocc+nvirt)     = mo_occ(1:nvirt)
         mo_list(2*nocc+nvirt+1:2*nocc+2*nvirt) = mo_virt(1:nvirt)
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=zero)
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
!        Switch to Dirac order and antisymmetrize
         vvoo_dims(1) = nvirt * 2
         vvoo_dims(2) = nvirt * 2
         vvoo_dims(3) = nocc * 2
         vvoo_dims(4) = nocc * 2
         scalar_dims  = 1
         ierr=talsh_tensor_construct(vvoo_tensor,C8,vvoo_dims,init_val=zero)
         ierr=talsh_tensor_construct(one_tensor,C8,scalar_dims,init_val=one)
         ierr=talsh_tensor_construct(mone_tensor,C8,scalar_dims,init_val=minus_one)
         ierr=talsh_tensor_contract("D(a,b,i,j)+=L(i,a,j,b)*R()",vvoo_tensor,moint_tensor,one_tensor)
         ierr=talsh_tensor_contract("D(a,b,i,j)+=L(j,a,i,b)*R()",vvoo_tensor,moint_tensor,mone_tensor)
         ierr=talsh_tensor_destruct(moint_tensor)

!        Make MP2 wave function by dividing by denominators and calculate energy
         t2_dims = vvoo_dims
         ierr=talsh_tensor_construct(t2_tensor,C8,t2_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(a,b,i,j)+=L(a,b,i,j)*R()",t2_tensor,vvoo_tensor,one_tensor)

!        Destruct remaining tensors and shut down library
         ierr=talsh_tensor_destruct(aoint_tensor)
         ierr=talsh_tensor_destruct(vvoo_tensor)
         ierr=talsh_tensor_destruct(t2_tensor)
         ierr = talsh_shutdown()

!        Make some noise so that we know we are leaving
         print*, " Leaving cc_driver routine"
         return
        end subroutine cc_driver

       end module exacorr_cc
