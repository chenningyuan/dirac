!This module contains the routines needed to calculate ae coupled cluster wave function
       module exacorr_cc

        implicit none

        complex(8), parameter :: zero=(0.D0,0.D0),one=(1.D0,0.D0),minus_one=(-1.D0,0.D0)

        private
        public cc_driver

       contains
!-------------------------------------------------------------------------------------------
        subroutine cc_driver (nocc,nvir,mo_occ,mo_vir)

!This subroutine drives the ccd energy calculation
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
 
!        Written by Lucas Visscher, winter 2016/2017
 
         use tensor_algebra
         use talsh
         use dirac_exatensor
         use exacorr_ao_to_mo

         implicit none
         integer, intent(in) :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in) :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer   :: ierr, i
         complex(8), pointer, contiguous:: ao_tens(:,:,:,:)
         complex(8), pointer :: energy(:)
         type(C_PTR):: body_p
         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: host_arg_max

         integer             :: nshells, number_aoshells   ! number of AO shells
         type(basis_func_info_t), allocatable :: shell(:)  ! array with basis function information
         integer             :: nao, number_aos

!        atomic integral tensor
         type(talsh_tens_t):: aoint_tensor
         integer:: aoint_dims(1:4)
!        sorted integral tensors
         type(talsh_tens_t):: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor
         integer:: oovv_dims(1:4),vvvv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4)
!        intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t):: goint_tensor, gvint_tensor
         type(talsh_tens_t):: aint_tensor, bint_tensor, hint_tensor
         integer:: goint_dims(1:2),gvint_dims(1:2)
         integer:: aint_dims(1:4),bint_dims(1:4),hint_dims(1:4)
!        solution tensors
         type(talsh_tens_t):: t2_tensor, s2_tensor
         integer:: t2_dims(1:4), s2_dims(1:4)
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t):: one_tensor,mone_tensor,energy_tensor
         integer:: scalar_dims(1)

         integer:: nmo(4)
         integer, allocatable :: mo_list(:)
         real(8), allocatable :: eps_occ(:),eps_vir(:)

!        CCD control variables
         integer :: iccd, nccd

!        These declarations are needed for work-arounds
         real(8) :: mp2_energy, ccd_energy
         complex(8), pointer, contiguous:: oovv_tens(:,:,:,:),t2_tens(:,:,:,:),s2_tens(:,:,:,:)
         integer :: j,a,b

!        Make some noise so that we know we are here
         print*, " Arrived in cc_driver routine"

!        Initialize libraries
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         print*, " Initialized talsh library"
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
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

!        Get orbital energies
         allocate (eps_occ(2*nocc))
         allocate (eps_vir(2*nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)

!        Get oovv tensor (assuming closed shells for the moment)
         oovv_dims(1:2) = nocc  * 2
         oovv_dims(3:4) = nvir * 2
         ierr=talsh_tensor_construct(oovv_tensor,C8,oovv_dims,init_val=zero)
         allocate (mo_list(2*nocc+2*nvir))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nocc
         nmo(4) = nvir
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nocc+nvir+1:2*nocc+nvir)     = mo_occ(1:nvir)
         mo_list(2*nocc+nvir+1:2*nocc+2*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,oovv_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Initialize t2 amplitudes with the MP2 values
         t2_dims(1:2) = nvir * 2
         t2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(t2_tensor,C8,t2_dims,init_val=zero)
!        This is a work-around, this transpose and copy should be cast in talsh calls
         ierr=talsh_tensor_get_body_access(oovv_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,oovv_tens,oovv_dims)
         ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,t2_tens,t2_dims)
         do j = 1, 2 * nocc
           do i = 1, 2 * nocc
             do b = 1, 2 * nvir
               do a = 1, 2 * nvir
                  t2_tens(a,b,i,j)   = dconjg(oovv_tens(i,j,a,b))
               end do 
             end do 
           end do 
         end do 
         call scale_with_denominators (t2_tensor,eps_occ,eps_vir)

!        Calculate and print the MP2 energy
!        This is a work-around, should be cast as a complete contraction of tensors
         mp2_energy = zero
         do j = 1, 2 * nocc
           do i = 1, j
             do b = 1, 2 * nvir
               do a = 1, b
                  mp2_energy = mp2_energy + oovv_tens(i,j,a,b) * t2_tens(a,b,i,j)
               end do 
             end do 
           end do 
         end do 
         print*," MP2 energy = ",mp2_energy

!        Prepare for CC by getting the other classes of integrals

!        Get oooo tensor
         oooo_dims = nocc * 2
         ierr=talsh_tensor_construct(oooo_tensor,C8,oooo_dims,init_val=zero)
         allocate (mo_list(4*nocc))
         nmo = nocc 
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,oooo_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Get vovo tensor
         vovo_dims(1) = nvir * 2
         vovo_dims(2) = nocc * 2
         vovo_dims(3) = nvir * 2
         vovo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(vovo_tensor,C8,vovo_dims,init_val=zero)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,vovo_tensor,.false.,nmo,mo_list)
         deallocate (mo_list)

!        Get vvvv tensor
         vvvv_dims = nvir * 2
         ierr=talsh_tensor_construct(vvvv_tensor,C8,vvvv_dims,init_val=zero)
         allocate (mo_list(4*nvir))
         nmo = nvir
         mo_list(1:nvir)          = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
         mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
         mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,vvvv_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Start CCD iterations

         nccd = 30 ! hard wired for the moment
!        Start CCD iterations

         nccd = 30 ! hard wired for the moment

         goint_dims = nocc  * 2
         ierr=talsh_tensor_construct(goint_tensor,C8,goint_dims,init_val=zero)
         gvint_dims = nvir * 2
         ierr=talsh_tensor_construct(gvint_tensor,C8,gvint_dims,init_val=zero)
         aint_dims  = nocc  * 2
         ierr=talsh_tensor_construct(aint_tensor,C8,aint_dims,init_val=zero)
         hint_dims(1) = nvir * 2
         hint_dims(2) = nocc * 2
         hint_dims(3) = nvir * 2
         hint_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(hint_tensor,C8,hint_dims,init_val=zero)
         s2_dims(1:2) = nvir * 2
         s2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(s2_tensor,C8,s2_dims,init_val=zero)

         do iccd = 1, nccd

!        calculate intermediaries
         ierr=talsh_tensor_init(goint_tensor)
         ierr=talsh_tensor_contract("G(k,i)+=V(k,l,c,d)*T(c,d,i,l)",goint_tensor,oovv_tensor,t2_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_tensor_init(gvint_tensor)
         ierr=talsh_tensor_contract("G(a,c)+=V(k,l,c,d)*T(a,d,k,l)",gvint_tensor,oovv_tensor,t2_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=talsh_tensor_init(aint_tensor)
         ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",aint_tensor,oooo_tensor)
         ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",aint_tensor,oovv_tensor,t2_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=13; return; endif
         ierr=talsh_tensor_init(hint_tensor)
         ierr=talsh_tensor_add("H(a,k,i,c)+=V(a,k,i,c)",hint_tensor,vovo_tensor)
!        ierr=talsh_tensor_contract("H(a,k,i,c)+=V(k,l,c,d)*T(a,d,i,l)",hint_tensor,oovv_tensor,t2_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=14; return; endif

!        calculate s2 tensor
!        ierr=talsh_tensor_init(s2_tensor)
!        ierr=talsh_tensor_add("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!        work-around
         ierr=talsh_tensor_get_body_access(s2_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,s2_tens,s2_dims)
         do j = 1, 2 * nocc
           do i = 1, 2 * nocc
             do b = 1, 2 * nvir
               do a = 1, 2 * nvir
                  s2_tens(a,b,i,j)   = dconjg(oovv_tens(i,j,a,b))
               end do 
             end do 
           end do 
         end do 
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,j,k)",s2_tensor,goint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,gvint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,aint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,k,i,c)*T(b,c,j,k)",s2_tensor,hint_tensor,t2_tensor)
!        ierr=talsh_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,vvvv_tensor,t2_tensor)

!        prepare for next iteration
         ierr=talsh_tensor_init(t2_tensor)
         ierr=talsh_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",t2_tensor,s2_tensor)
         call scale_with_denominators (t2_tensor,eps_occ,eps_vir)

         ccd_energy = zero
         do j = 1, 2 * nocc
           do i = 1, j
             do b = 1, 2 * nvir
               do a = 1, b
                  ccd_energy = ccd_energy + oovv_tens(i,j,a,b) * t2_tens(a,b,i,j)
               end do 
             end do 
           end do 
         end do 
         print*," CCD energy = ",ccd_energy

         end do

!        Destruct tensors and shut down library
         ierr=talsh_tensor_destruct(aoint_tensor)
         ierr=talsh_tensor_destruct(oovv_tensor)
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(t2_tensor)
         ierr=talsh_tensor_destruct(goint_tensor)
         ierr=talsh_tensor_destruct(gvint_tensor)
         ierr=talsh_tensor_destruct(aint_tensor)
         ierr=talsh_tensor_destruct(hint_tensor)
         ierr = talsh_shutdown()
         deallocate (eps_occ)
         deallocate (eps_vir)

!        Make some noise so that we know we are leaving
         print*, " Leaving cc_driver routine"
         return
        end subroutine cc_driver

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,nmo,mo_list)

         use tensor_algebra
         use talsh
         use exacorr_ao_to_mo
       
         implicit none

         type(talsh_tens_t), intent(inout)  :: aoint_tensor,integral_tensor
         integer, intent(in)                :: nmo(4), mo_list(:)
         logical, intent(in)                :: antisymmetrize

         integer                            :: ierr, moint_dims(4)
         type(talsh_tens_t)                 :: moint_tensor

         moint_dims = 2 * nmo
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=zero)
!        Transform integrals
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)

        end subroutine get_integral_tensor

        subroutine mulliken_to_dirac_sort (mulliken_tensor,dirac_tensor,antisymmetrize)

         use tensor_algebra
         use talsh
       
         implicit none

         type(talsh_tens_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         logical, intent(in)                :: antisymmetrize

         complex(8), pointer, contiguous:: mulliken_tens(:,:,:,:),dirac_tens(:,:,:,:)
         type(C_PTR):: body_p

         integer:: mulliken_dims(1:4),dirac_dims(1:4)

         integer   :: ierr, tens_rank, p, q, r, s, mismatch

!        Make input tensor accessible
         ierr = talsh_tensor_dimensions(mulliken_tensor,tens_rank,mulliken_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: input tensor corrupted'
         ierr=talsh_tensor_get_body_access(mulliken_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,mulliken_tens,mulliken_dims(1:4)) 

!        Make output tensor accessible
         ierr = talsh_tensor_dimensions(dirac_tensor,tens_rank,dirac_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: input tensor corrupted'
         ierr=talsh_tensor_get_body_access(dirac_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,dirac_tens,dirac_dims(1:4)) 

!        Check consistency
         mismatch = 0
         if (mulliken_dims(1) .ne. dirac_dims(1)) mismatch = mismatch + 1
         if (mulliken_dims(2) .ne. dirac_dims(3)) mismatch = mismatch + 1
         if (mulliken_dims(3) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (mulliken_dims(4) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (antisymmetrize) then
            if (dirac_dims(1) .ne. dirac_dims(2) .and. dirac_dims(3) .ne. dirac_dims(4)) mismatch = mismatch + 1
         end if
         if (mismatch.gt.0) stop 'dimension mismatch in mulliken_to_dirac_sort'

         do s = 1, dirac_dims(4)
           do r = 1, dirac_dims(3)
             do q = 1, dirac_dims(2)
               do p = 1, dirac_dims(1)
                  dirac_tens(p,q,r,s) = mulliken_tens(p,r,q,s)
               end do
             end do
           end do
         end do

         if (antisymmetrize) then
            if (dirac_dims(1) .eq. dirac_dims(2)) then
              do s = 1, dirac_dims(4)
                do r = 1, dirac_dims(3)
                  do q = 1, dirac_dims(2)
                    do p = 1, dirac_dims(1)
                     dirac_tens(p,q,r,s) = dirac_tens(p,q,r,s) - mulliken_tens(q,r,p,s)
                    end do
                  end do
                end do
              end do
            else
              do s = 1, dirac_dims(4)
                do r = 1, dirac_dims(3)
                  do q = 1, dirac_dims(2)
                    do p = 1, dirac_dims(1)
                     dirac_tens(p,q,r,s) = dirac_tens(p,q,r,s) - mulliken_tens(p,s,q,r)
                    end do
                  end do
                end do
              end do
            end if
         end if

        end subroutine mulliken_to_dirac_sort
        
        subroutine scale_with_denominators (t2_tensor,eps_occ,eps_vir)

         use tensor_algebra
         use talsh
       
         implicit none

         type(talsh_tens_t), intent(inout)  :: t2_tensor
         complex(8), pointer, contiguous:: t2_tens(:,:,:,:)
         type(C_PTR):: body_p
         real(8) :: eps_occ(:),eps_vir(:),denominator

         integer:: t2_dims(1:4)

         integer   :: ierr, tens_rank, a, b, i, j

!        Get tensor dimension of input and output tensors
         ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: t2 tensor corrupted'
         ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,t2_tens,t2_dims(1:4)) ! to use <t2_tens> as a regular Fortran 4d array

         do j = 1, t2_dims(4)
           do i = 1, t2_dims(3)
             do b = 1, t2_dims(2)
               do a = 1, t2_dims(1)
                  denominator = eps_occ(i) + eps_occ(j) - eps_vir(a) - eps_vir(b)
                  t2_tens(a,b,i,j) = t2_tens(a,b,i,j) / denominator
               end do
             end do
           end do
         end do
        
       end subroutine scale_with_denominators

       subroutine get_orbital_energies (nmo, mo_list, eps)

!This subroutine reads the desired subset of mo energies from DIRAC

!        Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 25 C)

         use dirac_MO
         use tensor_algebra
         use talsh

         implicit none
         integer, intent(in) :: nmo          ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         real(8) :: eps(:)

         integer             :: j, k, ierr
         type(MO)            :: qmo
         integer             :: ifrp = 1     ! we do not use symmetry, so ifrp (g/u) is always 1
         integer             :: nz = 4       ! we do not use symmetry, so nz is always 4

         call get_dirac_mo_coefficients (ifrp,qmo,mo_list,nmo)

         do k = 1, nmo
            eps(k)     = qmo%energy(k)
            eps(k+nmo) = qmo%energy(k)
         end do
         call delete (qmo)

        end subroutine get_orbital_energies

       end module exacorr_cc
