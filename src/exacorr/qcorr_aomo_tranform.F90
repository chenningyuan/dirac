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
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
!Note that all integrals are printed to standard output, should be replaced by checking the MP2 energy or some other printable quantity
!
         use dirac_MO
         use dirac_exatensor
         use tensor_algebra
         use talsh
         use, intrinsic:: ISO_C_BINDING
         use module_interest_interface
         use module_interest_eri


         implicit none
         integer, intent(in) :: nao 
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: ifrp = 1     ! we do not use symmetry, so ifrp is always 1
         integer             :: nz = 4       ! we do not use symmetry, so nz is always 4
         integer             :: i, j, k, l, imo, ierr
         type(MO)            :: qmo
         integer             :: ipqtoq(4)   ! only needed for the debug print, don't forget to remove later
         integer             :: nshells, number_aoshells   ! number of AO shells
         type(basis_func_info_t), allocatable :: shell(:)  ! array with basis function information

         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         type(talsh_tens_t):: aoint_tensor, mocoef_tensor,a1int_tensor
         type(talsh_tens_t):: a2int_tensor, a3int_tensor, moint_tensor
         integer(C_INT):: aoint_dims(1:4),mocoef_dims(1:2),a1int_dims(1:4)
         integer(C_INT):: a2int_dims(1:4),a3int_dims(1:4),moint_dims(1:4)
         type(talsh_task_t):: tsk
         real(C_DOUBLE):: flops,tm,tmc,tmi,tmo,cn1
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: sts,host_arg_max

         real(8)   :: vd,vl,vr
         complex(8):: zero=(0.0,0.0), aoint_value=(1.D-6,0.0), mocoef_value=(1.D-3,0.0)
         real(8), pointer, contiguous:: camo_tens(:,:),ao_tens(:,:,:,:),mo_tens(:,:,:,:)
         type(C_PTR):: body_p

!        Make some noise so that we know we are here
         print*, " Arrived in compute_mo_integrals routine"

!        Initialize libraries
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         print*, " Initialized talsh library"
         write(*,'(" Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
         call interest_initialize
         print*, " Initialized interest library"

!        Retrieve basis set information from DIRAC common blocks
         nshells = number_aoshells()
         allocate(shell(nshells))
         print*, " Allocated ",nshells," shell variables"
         call initialize_shells(shell)
         do i = 1, nshells
            if (compress_unc_shell(shell(i)) == 1) call quit ('This version can not handle contracted functions')
         enddo
         print*, " Retrieved basis set information from DIRAC common blocks"
         
!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         mocoef_dims(1)  = nao
         mocoef_dims(2)  = nmo(1)
         aoint_dims      = nao
         a1int_dims(1)   = nmo(1)
         a1int_dims(2:4) = nao
         a2int_dims(1:2) = nmo(1:2)
         a2int_dims(3:4) = nao
         a3int_dims(1:3) = nmo(1:3)
         a3int_dims(4)    = nao
         moint_dims(1:4)  = nmo(1:4)

!        Initialize tensors with MO coefficients
         ierr=talsh_tensor_construct(mocoef_tensor,R8,mocoef_dims,init_val=zero)
         ierr=talsh_tensor_get_body_access(mocoef_tensor,body_p,R8,0,DEV_HOST)
         call c_f_pointer(body_p,camo_tens,mocoef_dims(1:2)) ! to use camo_tens as a regular Fortran 2d array

!        This array is only needed to print the MO coefficients using the DIRAC routine PRQMAT
         do i = 1, 4
            ipqtoq(i) = i
         end do
         imo = 1
         do i = 1, 1  ! CAUTION: should be generalized to different coefficient ranges for each index (e.g. for MP2)
            call create (qmo, nao, nmo(i), nz)
            call get_dirac_mo_coefficients (ifrp,qmo,mo_list(imo:imo+nmo(i)-1),nmo(i))
!           uncommenting this can be useful for debugging, this piece of code is to be deleted later
!           test that we have read the coefficients correctly
!           print*,"Mo coefficients for index ",i
!           CALL PRQMAT(qmo%coeff,nao,nmo(i),nao,nmo(i),nz,ipqtoq,6)
!           copy the real coefficients into the tensor (these are the alpha coefficients)
!           in the transformation below we assume real coefficients and furthermore (Calpha = Cbeta) so this
!           will only work for non-relativistic (levy-Leblond) or X2C spinfree calculations
            camo_tens = qmo%coeff(:,:,1)
            call delete (qmo)
            imo = imo + nmo(i)
         end do

!        Get tensor for AO integrals
         ierr=talsh_tensor_construct(aoint_tensor,R8,aoint_dims,init_val=zero)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,R8,0,DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) !now you can use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens,nshells,shell)

!        Do 1st quarter-transformation
         ierr=talsh_tensor_construct(a1int_tensor,R8,a1int_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(i,q,r,s)+=L(p,q,r,s)*R(p,i)",a1int_tensor,aoint_tensor,mocoef_tensor,talsh_task=tsk)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_task_wait(tsk,sts); if(ierr.ne.TALSH_SUCCESS.or.sts.ne.TALSH_TASK_COMPLETED) then; ierr=12; return; endif
         ierr=talsh_task_time(tsk,tm,tmc,tmi,tmo)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=13; return; endif
!        Calculate the speed with which this was done
         vd=1_C_SIZE_T; do i=1,4; vd=vd*a1int_dims(i); enddo !volume of the destination tensor
         vl=1_C_SIZE_T; do i=1,4; vl=vl*aoint_dims(i); enddo !volume of the left tensor
         vr=1_C_SIZE_T; do i=1,2; vr=vr*mocoef_dims(i); enddo !volume of the right tensor
         flops=dsqrt(dble(vd)*dble(vl)*dble(vr))*2d0 !number of floating point operations
         write(*,'(1x,A,D8.2)') 'Contraction speed 1st quarter in flops = ',flops/tm
!        Destruct task and ao integrals:
         ierr=talsh_task_destruct(tsk); if(ierr.ne.TALSH_SUCCESS) then; ierr=14; return; endif
         ierr=talsh_tensor_destruct(aoint_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=15; return; endif

!        Do 2nd quarter-transformation
         ierr=talsh_tensor_construct(a2int_tensor,R8,a2int_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoef_tensor,talsh_task=tsk)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         ierr=talsh_task_wait(tsk,sts); if(ierr.ne.TALSH_SUCCESS.or.sts.ne.TALSH_TASK_COMPLETED) then; ierr=22; return; endif
         ierr=talsh_task_time(tsk,tm,tmc,tmi,tmo)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=23; return; endif
!        Calculate the speed with which this was done
         vd=1_C_SIZE_T; do i=1,4; vd=vd*a2int_dims(i); enddo !volume of the destination tensor
         vl=1_C_SIZE_T; do i=1,4; vl=vl*a1int_dims(i); enddo !volume of the left tensor
         vr=1_C_SIZE_T; do i=1,2; vr=vr*mocoef_dims(i); enddo !volume of the right tensor
         flops=dsqrt(dble(vd)*dble(vl)*dble(vr))*2d0 !number of floating point operations
         write(*,'(1x,A,D8.2)') 'Contraction speed 2nd quarter in flops = ',flops/tm
!        Destruct task and a1 integrals:
         ierr=talsh_task_destruct(tsk); if(ierr.ne.TALSH_SUCCESS) then; ierr=24; return; endif
         ierr=talsh_tensor_destruct(a1int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=25; return; endif

!        Do 3rd quarter-transformation
         ierr=talsh_tensor_construct(a3int_tensor,R8,a3int_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(i,j,k,s)+=L(i,j,r,s)*R(r,k)",a3int_tensor,a2int_tensor,mocoef_tensor,talsh_task=tsk)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=talsh_task_wait(tsk,sts); if(ierr.ne.TALSH_SUCCESS.or.sts.ne.TALSH_TASK_COMPLETED) then; ierr=32; return; endif
         ierr=talsh_task_time(tsk,tm,tmc,tmi,tmo)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=33; return; endif
!        Calculate the speed with which this was done
         vd=1_C_SIZE_T; do i=1,4; vd=vd*a3int_dims(i); enddo !volume of the destination tensor
         vl=1_C_SIZE_T; do i=1,4; vl=vl*a2int_dims(i); enddo !volume of the left tensor
         vr=1_C_SIZE_T; do i=1,2; vr=vr*mocoef_dims(i); enddo !volume of the right tensor
         flops=dsqrt(dble(vd)*dble(vl)*dble(vr))*2d0 !number of floating point operations
         write(*,'(1x,A,D8.2)') 'Contraction speed 3rd quarter in flops = ',flops/tm
!        Destruct task and a2 integrals:
         ierr=talsh_task_destruct(tsk); if(ierr.ne.TALSH_SUCCESS) then; ierr=34; return; endif
         ierr=talsh_tensor_destruct(a2int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=35; return; endif

!        Do final quarter-transformation
         ierr=talsh_tensor_construct(moint_tensor,R8,moint_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoef_tensor,talsh_task=tsk)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=talsh_task_wait(tsk,sts); if(ierr.ne.TALSH_SUCCESS.or.sts.ne.TALSH_TASK_COMPLETED) then; ierr=42; return; endif
         ierr=talsh_task_time(tsk,tm,tmc,tmi,tmo)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=43; return; endif
!        Calculate the speed with which this was done
         vd=1_C_SIZE_T; do i=1,4; vd=vd*moint_dims(i); enddo !volume of the destination tensor
         vl=1_C_SIZE_T; do i=1,4; vl=vl*a3int_dims(i); enddo !volume of the left tensor
         vr=1_C_SIZE_T; do i=1,2; vr=vr*mocoef_dims(i); enddo !volume of the right tensor
         flops=dsqrt(dble(vd)*dble(vl)*dble(vr))*2d0 !number of floating point operations
         write(*,'(1x,A,D8.2)') 'Contraction speed 4rd quarter in flops = ',flops/tm
!        Destruct task and a3 integrals:
         ierr=talsh_task_destruct(tsk); if(ierr.ne.TALSH_SUCCESS) then; ierr=44; return; endif
         ierr=talsh_tensor_destruct(a3int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=45; return; endif

!        Print transformed integrals to check correctness
         ierr=talsh_tensor_get_body_access(moint_tensor,body_p,R8,0,DEV_HOST)
         call c_f_pointer(body_p,mo_tens,moint_dims(1:4))

!        uncommenting this can be useful for debugging, this piece of code is to be deleted later
!         print*,"Transformed integrals larger than 1.D-15"
!         do l = 1, nmo(4)
!            do k = 1, nmo(3)
!               do j = 1, nmo(2)
!                  do i = 1, nmo(1)
!                     if (abs(mo_tens(i,j,k,l)) > 1.D-15) print*,i,j,k,l,mo_tens(i,j,k,l)
!                  end do
!               end do
!            end do
!         end do

!        Write to old format MO-integrals file to allow for easy debugging to the integral transform
         call write_to_mdcint (nmo,mo_tens)
     
!        Shut down library
         ierr = talsh_shutdown()

!        Make some noise so that we know we are leaving
         print*, " Leaving compute_mo_integrals routine"
         return
        end subroutine compute_mo_integrals

        subroutine compute_ao_integrals (nao,ao_tens,nshells,gto)
!This subroutine drives the generation of AO integrals
!
         use dirac_exatensor
         use module_interest_eri


         implicit none
         integer, intent(in) :: nao, nshells ! number of basis functions and basis function shells
         real(8),intent(out) :: ao_tens(nao,nao,nao,nao)
         type(basis_func_info_t), intent(in) :: gto(:)  ! array with basis function information

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         real(8), parameter :: fijkl = 1.0
         integer, parameter :: l_max = 6
         integer            :: ncartesian(l_max)

!        calculate number of Cartesian function for each l-value
         ncartesian(1) = 1
         do i = 2, l_max
            ncartesian(i) = ncartesian(i-1) + i
         end do

         loff = 0
         do lsh = 1, nshells
            ll   =  gto(lsh)%orb_momentum
            el   =  gto(lsh)%exponent(1)
            xl   =  gto(lsh)%coord(1)
            yl   =  gto(lsh)%coord(2)
            zl   =  gto(lsh)%coord(3)
            cl   =  gto(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            koff = 0
            do ksh = 1, nshells
               lk   =  gto(ksh)%orb_momentum
               ek   =  gto(ksh)%exponent(1)
               xk   =  gto(ksh)%coord(1)
               yk   =  gto(ksh)%coord(2)
               zk   =  gto(ksh)%coord(3)
               ck   =  gto(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               joff = 0
               do jsh = 1, nshells
                  lj   =  gto(jsh)%orb_momentum
                  ej   =  gto(jsh)%exponent(1)
                  xj   =  gto(jsh)%coord(1)
                  yj   =  gto(jsh)%coord(2)
                  zj   =  gto(jsh)%coord(3)
                  cj   =  gto(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  ioff = 0
                  do ish = 1, nshells
                     li   =  gto(ish)%orb_momentum
                     ei   =  gto(ish)%exponent(1)
                     xi   =  gto(ish)%coord(1)
                     yi   =  gto(ish)%coord(2)
                     zi   =  gto(ish)%coord(3)
                     ci   =  gto(ish)%coefficient(1)
                     ni   =  ncartesian(li)
!                    output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                     ijkl = 0
                     do l = 1, nl
                        do k = 1, nk
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = ijkl + 1
                                 ao_tens(i+ioff,j+joff,k+koff,l+loff) = gout(ijkl)
                              end do
                           end do
                        end do
                     end do
                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do
         
         end subroutine compute_ao_integrals

       end module exacorr_ao_to_mo
