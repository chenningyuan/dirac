!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.
       module exacorr_ao_to_mo


        implicit none
        private
!PARAMETERS:
        integer, parameter, private:: REALD=8 !default real for this module
!DERIVED TYPES:
!VISIBILITY:
        public compute_mo_integrals

       contains
!-------------------------------------------------------------------------------------------
        subroutine compute_mo_integrals (nao, nmo, mo_list)
!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!
         use dirac_MO
         use tensor_algebra
         use talsh
         use, intrinsic:: ISO_C_BINDING

         implicit none
         integer, intent(in) :: nao 
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: ifrp = 1     ! we do not use symmetry, so ifrp is always 1
         integer             :: nz = 4       ! we do not use symmetry, so nz is always 4
         integer             :: i, imo, ierr
         type(MO)            :: qmo
         integer             :: ipqtoq(4)   ! only needed for the debug print, don't forget to remove later

!        we will test the construction and contraction of two tensors, needs to be moved to a proper subroutine later on
         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         type(talsh_tens_t):: aoint_tensor, mocoef_tensor,a1int_tensor
         integer(C_INT):: aoint_dims(1:4),mocoef_dims(1:2),a1int_dims(1:4)
         type(talsh_task_t):: tsk
         real(C_DOUBLE):: flops,tm,tmc,tmi,tmo,cn1
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: sts,host_arg_max

         real(8)   :: vd,vl,vr
         complex(8):: zero=(0.0,0.0), aoint_value=(1.D-6,0.0), mocoef_value=(1.D-3,0.0)


!        Make some noise so that we know we are here
         print*, " Arrived in compute_mo_integrals routine"
!        Initialize library
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         write(*,'("Status ",i11,": Size (Bytes) = ",i13, &
         ": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
         
!        Create and initialize tensors
         aoint_dims = nao
         mocoef_dims(1) = nao
         mocoef_dims(2) = nmo(1)
         a1int_dims(1)  = nmo(1)
         a1int_dims(2:4) = nao
         ierr=talsh_tensor_construct(aoint_tensor,R8,aoint_dims,init_val=aoint_value)
         ierr=talsh_tensor_construct(a1int_tensor,R8,a1int_dims,init_val=zero)
         ierr=talsh_tensor_construct(mocoef_tensor,R8,mocoef_dims,init_val=mocoef_value)
!        ierr=talsh_tensor_construct(mocoef_tensor,R8,mocoef_dims,init_method=mocoef_init(mo_list(1:nmo(i)))) ! here I need information on how to write such an intialization method....

!        Do dummy 1st quarter-transformation
         ierr=talsh_tensor_contract("D(i,q,r,s)+=L(p,q,r,s)*R(p,i)",a1int_tensor,aoint_tensor,mocoef_tensor,talsh_task=tsk)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=14; return; endif
         ierr=talsh_task_wait(tsk,sts); if(ierr.ne.TALSH_SUCCESS.or.sts.ne.TALSH_TASK_COMPLETED) then; ierr=15; return; endif
         ierr=talsh_task_time(tsk,tm,tmc,tmi,tmo)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=16; return; endif

!        Calculate the speed with which this was done
         vd=1_C_SIZE_T; do i=1,4; vd=vd*a1int_dims(i); enddo !volume of the destination tensor
         vl=1_C_SIZE_T; do i=1,4; vl=vl*aoint_dims(i); enddo !volume of the left tensor
         vr=1_C_SIZE_T; do i=1,2; vr=vr*mocoef_dims(i); enddo !volume of the right tensor
         flops=dsqrt(dble(vd)*dble(vl)*dble(vr))*2d0 !number of floating point operations
         write(*,'(1x,A,D8.2)') 'Contraction speed in flops = ',flops/tm

!        Initialize tensors with MO coefficients

         do i = 1, 4
            ipqtoq(i) = i
         end do

         imo = 1 
         do i = 1, 4
            call create (qmo, nao, nmo(i), nz)
            call get_dirac_mo_coefficients (ifrp,qmo,mo_list(imo:imo+nmo(i)-1),nmo(i))

!           test that we have read the coefficients correctly
            print*,"Mo coefficients for index ",i
            CALL PRQMAT(qmo%coeff,nao,nmo(i),nao,nmo(i),nz,ipqtoq,6)

            call delete (qmo)
            imo = imo + nmo(i)
         end do

!        Initialize AO integral tensors (to be done after we know how to make an initialization method)

!        Shut down library
         ierr = talsh_shutdown()

!        Make some noise so that we know we are leaving
         print*, " Leaving compute_mo_integrals routine"
         return
        end subroutine compute_mo_integrals

       end module exacorr_ao_to_mo
