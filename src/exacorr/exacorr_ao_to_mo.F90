!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.
       module exacorr_ao_to_mo


        implicit none
        private
        public compute_ao_integrals
        public compute_mo_integrals
        public get_integral_tensor

       contains
!-------------------------------------------------------------------------------------------
        subroutine compute_mo_integrals (nmo, mo_list)

!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
!The integrals are written to MDCINT, to allow for testing with unmodified DIRAC correlation modules
 
!        Written by Lucas Visscher, winter 2016/2017
 
         use dirac_MO
         use exacorr_datatypes
         use tensor_algebra
         use talsh
         use, intrinsic:: ISO_C_BINDING
         use module_interest_interface
         use module_interest_eri


         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: i, ierr
         integer             :: nao, number_aos            ! the size of the atomic orbital basis
         integer             :: nshells, number_aoshells   ! number of AO shells
         type(basis_func_info_t), allocatable :: shell(:)  ! array with basis function information

         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         type(talsh_tens_t):: aoint_tensor
         type(talsh_tens_t):: moint_tensor
         integer(C_INT):: aoint_dims(1:4),moint_dims(1:4)
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: host_arg_max

         real(8)   :: vd,vl,vr
         complex(8):: zero=(0.D0,0.D0)
         complex(8), pointer, contiguous:: ao_tens(:,:,:,:),mo_tens(:,:,:,:)
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
         nao     = number_aos(1)
         nshells = number_aoshells()
         allocate(shell(nshells))
         print*, " Allocated ",nshells," shell variables"
         call initialize_shells(shell)
         do i = 1, nshells
            call compress_unc_shell(shell(i),ierr)
            if (ierr == 1) call quit ('This version can not handle contracted functions')
         enddo
         print*, " Retrieved basis set information from DIRAC common blocks"
         
!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         aoint_dims      = nao
         moint_dims(1:4) = 2 * nmo(1:4)

!        Get tensor for AO integrals
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=zero)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens,nshells,shell)

!        Transform to MO basis
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=zero)
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)

!        Write to old format MO-integrals file to allow for easy debugging of the integral transform
         ierr=talsh_tensor_get_body_access(moint_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,mo_tens,moint_dims(1:4))
         call write_to_mdcint (nmo,2,mo_tens)
     
!        Destruct remaining tensors and shut down library
         ierr=talsh_tensor_destruct(aoint_tensor)
         ierr=talsh_tensor_destruct(moint_tensor)
         ierr = talsh_shutdown()

!        Make some noise so that we know we are leaving
         print*, " Leaving compute_mo_integrals routine"
         return
        end subroutine compute_mo_integrals

        subroutine get_mocoef_tensor (nao, nmo, mo_list, mocoefa_tensor, mocoefb_tensor)

!This subroutine reads the desired subset of mo coefficients from DIRAC
 
!        Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 30 C)
 
         use dirac_MO
         use tensor_algebra
         use talsh

         implicit none
         integer, intent(in) :: nao          ! the length of the ao basis
         integer, intent(in) :: nmo          ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         type(talsh_tens_t), intent(inout) :: mocoefa_tensor,mocoefb_tensor

         integer             :: j, k, ierr
         type(MO)            :: qmo
         integer             :: ifrp = 1     ! we do not use symmetry, so ifrp (g/u) is always 1
         integer             :: nz = 4       ! we do not use symmetry, so nz is always 4
         integer             :: tens_rank, mocoef_dims(1:2)

         complex(8), pointer, contiguous:: camo_tens(:,:),cbmo_tens(:,:)
         type(C_PTR):: body_p


         ierr=talsh_tensor_dimensions(mocoefa_tensor,tens_rank,mocoef_dims)
         ierr=talsh_tensor_get_body_access(mocoefa_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,camo_tens,mocoef_dims(1:2)) ! to use camo_tens as a regular Fortran 2d array
         ierr=talsh_tensor_get_body_access(mocoefb_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,cbmo_tens,mocoef_dims(1:2)) ! to use cbmo_tens as a regular Fortran 2d array

         call get_dirac_mo_coefficients (ifrp,qmo,mo_list,nmo)

!        Copy the coefficients into the alpha and beta tensors, using the quaternion relations as defined in 
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         do k = 1, nmo
            do j = 1, nao
                camo_tens(j,k) = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
                cbmo_tens(j,k) = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
                camo_tens(j,k+nmo) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
                cbmo_tens(j,k+nmo) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
            end do
         end do
         call delete (qmo)

         return
        end subroutine get_mocoef_tensor

        subroutine compute_ao_integrals (nao,ao_tens,nshells,gto)

!This subroutine drives the generation of AO integrals by InteRest
!
         use exacorr_datatypes
         use module_interest_eri

         implicit none
         integer, intent(in) :: nao, nshells ! number of basis functions and basis function shells
         complex(8),intent(out) :: ao_tens(nao,nao,nao,nao)
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
                                 ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
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

        subroutine tens_init_ao_integrals(ao_tens,nao,first_ao,ierr)

         use exacorr_datatypes
         use module_interest_eri

!Initialization routine that is used to fill a tensor with AO integrals

         implicit none

         integer, intent(in):: nao(4)
         complex(8),intent(out) :: ao_tens(nao(1),nao(2),nao(3),nao(4))
         integer, intent(in):: first_ao(1:4)
         integer, intent(out), optional:: ierr

         integer :: nshells(4)
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         real(8), parameter :: fijkl = 1.0
         integer, parameter :: l_max = 6
         integer            :: ncartesian(l_max)

!        calculate number of Cartesian functions for each l-value
         ncartesian(1) = 1
         do i = 2, l_max
            ncartesian(i) = ncartesian(i-1) + i
         end do

!        Allocate and get shell information for all shells that we will treat
         call fetch_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call fetch_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call fetch_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call fetch_gtos(first_ao(4),nao(4),gto4,nshells(4))

!        Intranode parallelization should be possible in the loops below because Interest is thread safe

         loff = 0
         do lsh = 1, nshells(4)
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            koff = 0
            do ksh = 1, nshells(3)
               lk   =  gto3(ksh)%orb_momentum
               ek   =  gto3(ksh)%exponent(1)
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               ck   =  gto3(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               joff = 0
               do jsh = 1, nshells(2)
                  lj   =  gto2(jsh)%orb_momentum
                  ej   =  gto2(jsh)%exponent(1)
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  cj   =  gto2(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  ioff = 0
                  do ish = 1, nshells(1)
                     li   =  gto1(ish)%orb_momentum
                     ei   =  gto1(ish)%exponent(1)
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ci   =  gto1(ish)%coefficient(1)
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
                                 ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
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

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)
         
        end subroutine tens_init_ao_integrals

        subroutine aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)

!        Transform integral tensor from ao-basis to mo-basis 
 
         use tensor_algebra
         use talsh

         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         type(talsh_tens_t), intent(inout)  :: aoint_tensor, moint_tensor

         type(talsh_tens_t):: mocoefa_tensor(4),mocoefb_tensor(4)
         type(talsh_tens_t):: a1int_tensor,a2int_tensor, a3int_tensor
         integer:: aoint_dims(1:4),mocoef_dims(1:2),a1int_dims(1:4)
         integer:: a2int_dims(1:4),a3int_dims(1:4),moint_dims(1:4)

         complex(8):: zero=(0.D0,0.D0)
         integer   :: ierr, tens_rank, i, imo, nao

!        Get tensor dimension of input and output tensors
         ierr = talsh_tensor_dimensions(aoint_tensor,tens_rank,aoint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: aoint tensor corrupted'
         ierr = talsh_tensor_dimensions(moint_tensor,tens_rank,moint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: moint tensor corrupted'

!        Initialize tensors with MO coefficients
         nao = aoint_dims(1)
         imo = 0
         do i = 1, 4
            mocoef_dims(1) = aoint_dims(i)
            if (mocoef_dims(1) .ne. nao)      stop 'program error: inconsistent dimensions'
            mocoef_dims(2) = moint_dims(i)
            if (mocoef_dims(2) .ne. 2*nmo(i)) stop 'program error: inconsistent dimensions'
            ierr=talsh_tensor_construct(mocoefa_tensor(i),C8,mocoef_dims,init_val=zero)
            ierr=talsh_tensor_construct(mocoefb_tensor(i),C8,mocoef_dims,init_val=zero)
            call get_mocoef_tensor (nao, nmo(i), mo_list(imo+1:imo+nmo(i)), mocoefa_tensor(i), mocoefb_tensor(i))
            imo = imo + nmo(i)
         end do

!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         a1int_dims(1)   = moint_dims(1)
         a1int_dims(2:4) = aoint_dims(2:4)
         a2int_dims(1:2) = moint_dims(1:2)
         a2int_dims(3:4) = aoint_dims(3:4)
         a3int_dims(1:3) = moint_dims(1:3)
         a3int_dims(4)   = aoint_dims(4)

!        Do 1st and 2nd quarter-transformation
         ierr=talsh_tensor_construct(a1int_tensor,C8,a1int_dims,init_val=zero)
         ierr=talsh_tensor_construct(a2int_tensor,C8,a2int_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)",a1int_tensor,mocoefa_tensor(1),aoint_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefa_tensor(2))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         ierr=talsh_tensor_init(a1int_tensor,zero)
         ierr=talsh_tensor_contract("D(i,q,r,s)=R+(p,i)*L(p,q,r,s)",a1int_tensor,mocoefb_tensor(1),aoint_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefb_tensor(2))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=22; return; endif

!        Destruct a1 and a1 integrals:
         ierr=talsh_tensor_destruct(a1int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=25; return; endif

!        Do 3rd and 4th quarter-transformation
         ierr=talsh_tensor_construct(a3int_tensor,C8,a3int_dims,init_val=zero)
         ierr=talsh_tensor_contract("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefa_tensor(3),a2int_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefa_tensor(4))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=talsh_tensor_init(a3int_tensor,zero)
         ierr=talsh_tensor_contract("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefb_tensor(3),a2int_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=32; return; endif
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefb_tensor(4))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=42; return; endif
 
!        Destruct a2 and a3 integrals:
         ierr=talsh_tensor_destruct(a2int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=35; return; endif
         ierr=talsh_tensor_destruct(a3int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=45; return; endif

         return

        end subroutine aomo_transform

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,nmo,mo_list)
 

         use tensor_algebra
         use talsh
       
         implicit none

         complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         type(talsh_tens_t), intent(inout)  :: aoint_tensor,integral_tensor
         integer, intent(in)                :: nmo(4), mo_list(:)
         logical, intent(in)                :: antisymmetrize

         integer                            :: ierr, moint_dims(4)
         type(talsh_tens_t)                 :: moint_tensor

         moint_dims = 2 * nmo
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
!        Transform integrals
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=talsh_tensor_destruct(moint_tensor)

        end subroutine get_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         use tensor_algebra
         use talsh
       
         implicit none

         complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         type(talsh_tens_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         logical, intent(in)                :: antisymmetrize

         integer:: mulliken_dims(1:4),dirac_dims(1:4)

         integer   :: ierr, tens_rank, mismatch

!        work-around: addition with permutation is not implemented
         type(talsh_tens_t):: one_tensor, minusone_tensor
         integer:: scalar_dims(1)
         scalar_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,scalar_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(minusone_tensor,C8,scalar_dims(1:0),init_val=MINUS_ONE)
!        end work-round

!        get dimensions of the tensors and check consistency
         ierr = talsh_tensor_dimensions(mulliken_tensor,tens_rank,mulliken_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input tensor corrupted in mulliken_to_dirac_sort'
         ierr = talsh_tensor_dimensions(dirac_tensor,tens_rank,dirac_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'output tensor corrupted in mulliken_to_dirac_sort'
         mismatch = 0
         if (mulliken_dims(1) .ne. dirac_dims(1)) mismatch = mismatch + 1
         if (mulliken_dims(2) .ne. dirac_dims(3)) mismatch = mismatch + 1
         if (mulliken_dims(3) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (mulliken_dims(4) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (antisymmetrize) then
            if (dirac_dims(1) .ne. dirac_dims(2) .and. dirac_dims(3) .ne. dirac_dims(4)) mismatch = mismatch + 1
         end if
         if (mismatch.gt.0) stop 'dimension mismatch in mulliken_to_dirac_sort'

         ierr=talsh_tensor_init(dirac_tensor)
!shorter ierr=talsh_tensor_add("D(p,q,r,s)+=M(p,r,q,s)",dirac_tensor,mulliken_tensor)
         ierr=talsh_tensor_contract("D(p,q,r,s)+=M(p,r,q,s)*C()",dirac_tensor,mulliken_tensor,one_tensor)

         if (antisymmetrize) then
            if (dirac_dims(1) .eq. dirac_dims(2)) then
!shorter        ierr=talsh_tensor_add("D(p,q,r,s)+=M(q,r,p,s)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
                ierr=talsh_tensor_contract("D(p,q,r,s)+=M(q,r,p,s)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
            else
!shorter        ierr=talsh_tensor_add("D(p,q,r,s)+=M(p,s,q,r)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
                ierr=talsh_tensor_contract("D(p,q,r,s)+=M(p,s,q,r)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
            end if
         end if

!        work-around: addition with permutation is not implemented
         ierr=talsh_tensor_destruct(one_tensor)
         ierr=talsh_tensor_destruct(minusone_tensor)
!        end work-round

        end subroutine mulliken_to_dirac_sort

       end module exacorr_ao_to_mo
