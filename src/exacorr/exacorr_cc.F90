!This module contains the routines needed to calculate ae coupled cluster wave function
       module exacorr_cc

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),ONE_QUARTER=(0.25D0,0.D0)

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

         implicit none
         integer, intent(in) :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in) :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer   :: ierr, i
         complex(8), pointer :: energy(:),result_tens(:)
         type(C_PTR):: body_p
         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: host_arg_max

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_tens_t):: fo_tensor, fv_tensor, fov_tensor
         type(talsh_tens_t):: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor, vooo_tensor, vovv_tensor
!        intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t):: goint_tensor, gvint_tensor, hoint_tensor, hvint_tensor, hovint_tensor
         type(talsh_tens_t):: aint_tensor, bint_tensor, hint_tensor
         integer:: goint_dims(1:2),gvint_dims(1:2),hoint_dims(1:2),hvint_dims(1:2),hovint_dims(1:2)
         integer:: aint_dims(1:4),bint_dims(1:4),hint_dims(1:4)
!        solution tensors
         type(talsh_tens_t):: t1_tensor, t2_tensor, tau_tensor, s1_tensor, s2_tensor
         integer:: t1_dims(1:2), t2_dims(1:4), tau_dims(1:4), s1_dims(1:2), s2_dims(1:4)
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t):: result_tensor, one_tensor
         integer:: result_dims(1),one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CCD control and result variables
         integer :: iccd, nccd
         real(8) :: mp2_energy, ccd_energy

!        Make some noise so that we know we are here
         print*, " Arrived in cc_driver routine"

!        Initialize libraries
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         print*, " Initialized talsh library"
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
         call interest_initialize
         print*, " Initialized interest library"

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         result_dims(1) = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(2*nocc))
         allocate (eps_vir(2*nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         
!        Get fixed one-body tensors
         call get_CC_fockmatrix (nocc,nvir,eps_occ,eps_vir, &
                                 fo_tensor,fov_tensor,fv_tensor)

!        Get fixed two-body tensors
         call get_CC_integrals (nocc,nvir,mo_occ,mo_vir, &
                                oooo_tensor,vooo_tensor,oovv_tensor, &
                                vovo_tensor,vovv_tensor,vvvv_tensor)

!        Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
         t1_dims(1) = nvir * 2
         t1_dims(2) = nocc * 2
         ierr=talsh_tensor_construct(t1_tensor,C8,t1_dims,init_val=ZERO)

!        Initialize t2 amplitudes with the MP2 values
         t2_dims(1:2) = nvir * 2
         t2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(t2_tensor,C8,t2_dims,init_val=ZERO)
!        ierr=talsh_tensor_add("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor,one_tensor)
!        end of work-around
         call scale_with_denominators (t2_tensor,eps_occ,eps_vir)

!        Calculate and print the MP2 energy
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor)
         mp2_energy=result_tens(1) * ONE_QUARTER
         print*," MP2 energy = ",mp2_energy

!        Start CCSD iterations

         nccd = 30 ! hard wired for the moment

         goint_dims = nocc  * 2
         ierr=talsh_tensor_construct(goint_tensor,C8,goint_dims,init_val=ZERO)
         gvint_dims = nvir * 2
         ierr=talsh_tensor_construct(gvint_tensor,C8,gvint_dims,init_val=ZERO)
         hoint_dims = nocc  * 2
         ierr=talsh_tensor_construct(hoint_tensor,C8,hoint_dims,init_val=ZERO)
         hvint_dims = nvir * 2
         ierr=talsh_tensor_construct(hvint_tensor,C8,hvint_dims,init_val=ZERO)
         hovint_dims(1) = nocc * 2
         hovint_dims(2) = nvir * 2
         ierr=talsh_tensor_construct(hovint_tensor,C8,hovint_dims,init_val=ZERO)
         aint_dims  = nocc  * 2
         ierr=talsh_tensor_construct(aint_tensor,C8,aint_dims,init_val=ZERO)
         hint_dims(1) = nvir * 2
         hint_dims(2) = nocc * 2
         hint_dims(3) = nvir * 2
         hint_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(hint_tensor,C8,hint_dims,init_val=ZERO)
         s1_dims(1) = nvir * 2
         s1_dims(2) = nocc * 2
         ierr=talsh_tensor_construct(s1_tensor,C8,s1_dims,init_val=ZERO)
         s2_dims(1:2) = nvir * 2
         s2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(s2_tensor,C8,s2_dims,init_val=ZERO)
         print*,"allocated all tensors"

         do iccd = 1, nccd

!        calculate intermediaries
         ierr=talsh_tensor_init(hoint_tensor)
         ierr=talsh_tensor_add("H(k,j)+=F(k,j)",hoint_tensor,fo_tensor)
         ierr=talsh_tensor_contract("H(k,j)+=V(k,l,c,d)*T(c,d,l,j)",hoint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         ierr=talsh_tensor_init(goint_tensor)
         ierr=talsh_tensor_add("G(k,j)+=H(k,j)",goint_tensor,hoint_tensor)
         ierr=talsh_tensor_contract("G(k,j)+=F(k,c)*T(c,j)",goint_tensor,fov_tensor,t1_tensor)
!check   ierr=talsh_tensor_contract("G(k,j)+=V+(i,c,k,l)*T(c,l)",goint_tensor,ooov_tensor,t1_tensor, &
!check                               scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_tensor_init(gvint_tensor)
         ierr=talsh_tensor_contract("G(a,c)+=V(k,l,c,d)*T(d,a,k,l)",gvint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=talsh_tensor_init(aint_tensor)
         ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",aint_tensor,oooo_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=13; return; endif
         ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",aint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=13; return; endif
         ierr=talsh_tensor_init(hint_tensor)
         ierr=talsh_tensor_add("H(a,i,c,k)+=V(a,i,c,k)",hint_tensor,vovo_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=14; return; endif
         ierr=talsh_tensor_contract("H(a,i,c,k)+=V(k,l,c,d)*T(a,d,l,i)",hint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=14; return; endif

!        For debugging
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=T(a,b,i,j)*T(a,b,i,j)",result_tensor,t2_tensor,t2_tensor, &
                                     scale=ONE_QUARTER)
         print*, "T2amp = ",result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=G(k,j)*G(k,j)",result_tensor,goint_tensor,goint_tensor)
         print*, "goint = ",result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=G(a,c)*G(a,c)",result_tensor,gvint_tensor,gvint_tensor)
         print*, "gvint = ",result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=A(k,l,i,j)*A(k,l,i,j)",result_tensor,aint_tensor,aint_tensor)
         print*, "a int = ",result_tens(1)/4.D0
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=H(a,i,c,k)*H(a,i,c,k)",result_tensor,hint_tensor,hint_tensor)
         print*, "h int = ",result_tens(1)

!        calculate s2 tensor
         ierr=talsh_tensor_init(s2_tensor)
!        ierr=talsh_tensor_add("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor,one_tensor)
!        end of work-around

!        This terms appears twice as we need to permute i and j
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,goint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,goint_tensor,t2_tensor)

!        This terms appears twice as we need to permute a and b
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,gvint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,gvint_tensor,t2_tensor)

!        The oooo contribution
         ierr=talsh_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,aint_tensor,t2_tensor, &
                                    scale=ONE_HALF)

!        The vvvv contribution
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,vvvv_tensor,t2_tensor, &
                                     scale=ONE_HALF)

!        This terms appears four times as we need to permute both i and j and a and b
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,i,c,k)*T(c,b,j,k)",s2_tensor,hint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,j,c,k)*T(c,b,k,i)",s2_tensor,hint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,i,c,k)*T(a,c,j,k)",s2_tensor,hint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,j,c,k)*T(a,c,k,i)",s2_tensor,hint_tensor,t2_tensor)

!        prepare for next iteration
         ierr=talsh_tensor_init(t2_tensor)
         ierr=talsh_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",t2_tensor,s2_tensor)
         call scale_with_denominators (t2_tensor,eps_occ,eps_vir)

         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor)
         ccd_energy=result_tens(1) * ONE_QUARTER
         print*," CCD energy = ",ccd_energy

         end do

!        Destruct tensors and shut down library
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

        subroutine get_CC_fockmatrix (nocc,nvir,eps_occ,eps_vir,fo_tensor,fov_tensor,fv_tensor)

!        Routine to get tensors with Fock matrix elements (with on the diagonal the orbital energies subtracted).
!        In this implementation we assume canonical orbitals, so we may simply construct zero tensors.
!        We also assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017
 
         use tensor_algebra
         use talsh

         implicit none
         integer, intent(in)               :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         real(8), intent(in)               :: eps_occ(:),eps_vir(:) ! the orbital energies
!        target tensors to be filled and given back to the caller
         type(talsh_tens_t), intent(inout) :: fo_tensor, fv_tensor, fov_tensor
         integer:: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)

         integer   :: ierr, i

!        Get foo tensor
         fo_dims(1:2) = nocc * 2
         ierr=talsh_tensor_construct(fo_tensor,C8,fo_dims,init_val=ZERO)

!        Get fvv tensor
         fv_dims(1:2) = nvir * 2
         ierr=talsh_tensor_construct(fv_tensor,C8,fv_dims,init_val=ZERO)

!        Get fo tensor
         fov_dims(1) = nocc * 2
         fov_dims(2) = nvir * 2
         ierr=talsh_tensor_construct(fov_tensor,C8,fov_dims,init_val=ZERO)

        end subroutine get_CC_Fockmatrix

        subroutine get_CC_integrals (nocc,nvir,mo_occ,mo_vir, &
                                     oooo_tensor,vooo_tensor,oovv_tensor, &
                                     vovo_tensor,vovv_tensor,vvvv_tensor)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017
 
         use tensor_algebra
         use talsh
         use dirac_exatensor
         use exacorr_ao_to_mo

         implicit none
         integer, intent(in)               :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)               :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(talsh_tens_t), intent(inout) :: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor, vooo_tensor, vovv_tensor
         integer:: oovv_dims(1:4),vvvv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4),vooo_dims(1:4),vovv_dims(1:4)

         integer   :: ierr, i
         complex(8), pointer, contiguous:: ao_tens(:,:,:,:)
         type(C_PTR):: body_p

         integer             :: nshells, number_aoshells   ! number of AO shells
         type(basis_func_info_t), allocatable :: shell(:)  ! array with basis function information
         integer             :: nao, number_aos

!        auxilliary integral tensors (only needed inside this routine)
         type(talsh_tens_t):: aoint_tensor, ovvo_tensor
         integer:: aoint_dims(1:4), ovvo_dims(1:4)

!        arrays needed to communicate with the MO integral generator
         integer:: nmo(4)
         integer, allocatable :: mo_list(:)
         real(8), allocatable :: eps_occ(:),eps_vir(:)

!        scalars (need to be defined as tensor types)
         type(talsh_tens_t):: minusone_tensor
         integer:: minusone_dims(1)

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)      ! number of basis functions
         nshells = number_aoshells()  ! basis function are organized in shells, this is their number
         allocate(shell(nshells))
         call initialize_shells(shell)
         do i = 1, nshells
            if (compress_unc_shell(shell(i)) == 1) call quit ('This version can not handle contracted functions')
         enddo
         print*, " Retrieved basis set information from DIRAC common blocks"

!        Get tensor with AO integrals
         aoint_dims      = nao
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens,nshells,shell)

!        Get antisymmetric oovv tensor
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         oovv_dims(1:2) = nocc * 2
         oovv_dims(3:4) = nvir * 2
         ierr=talsh_tensor_construct(oovv_tensor,C8,oovv_dims,init_val=ZERO)
!        whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
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

!        Get antisymmetric oooo tensor
         oooo_dims = nocc * 2
         ierr=talsh_tensor_construct(oooo_tensor,C8,oooo_dims,init_val=ZERO)
         allocate (mo_list(4*nocc))
         nmo = nocc 
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,oooo_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         vovo_dims(1) = nvir * 2
         vovo_dims(2) = nocc * 2
         vovo_dims(3) = nvir * 2
         vovo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(vovo_tensor,C8,vovo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,vovo_tensor,.false.,nmo,mo_list)
         deallocate (mo_list)
!        Get ovvo tensor
         ovvo_dims(1) = nocc * 2
         ovvo_dims(2) = nvir * 2
         ovvo_dims(3) = nvir * 2
         ovvo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nvir
         nmo(4) = nocc
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nvir+nocc+1:2*nvir+nocc)     = mo_vir(1:nvir)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,ovvo_tensor,.false.,nmo,mo_list)
         deallocate (mo_list)
!        Make anti-symmetrized vovo tensor
         minusone_dims(1) = 1
         ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",vovo_tensor,ovvo_tensor,minusone_tensor)
         ierr=talsh_tensor_destruct(ovvo_tensor)

!        Get vvvv tensor
         vvvv_dims = nvir * 2
         ierr=talsh_tensor_construct(vvvv_tensor,C8,vvvv_dims,init_val=ZERO)
         allocate (mo_list(4*nvir))
         nmo = nvir
         mo_list(1:nvir)          = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
         mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
         mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,vvvv_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Get vooo tensor
         vooo_dims(1)   = nvir * 2
         vooo_dims(2:4) = nocc * 2
         ierr=talsh_tensor_construct(vooo_tensor,C8,vooo_dims,init_val=ZERO)
         allocate (mo_list(nvir+3*nocc))
         nmo(1)   = nvir
         nmo(2:4) = nocc
         mo_list(1:nvir)                    = mo_vir(1:nvir)
         mo_list(nvir+1:nvir+nocc)          = mo_occ(1:nocc)
         mo_list(nvir+nocc+1:nvir+2*nocc)   = mo_occ(1:nocc)
         mo_list(nvir+2*nocc+1:nvir+3*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,vooo_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Get vovv tensor
         vovv_dims(1)   = nvir * 2
         vovv_dims(2)   = nocc * 2
         vovv_dims(3:4) = nvir * 2
         ierr=talsh_tensor_construct(vovv_tensor,C8,vovv_dims,init_val=ZERO)
         allocate (mo_list(3*nvir+nocc))
         nmo(1)   = nvir
         nmo(2)   = nvir
         nmo(3)   = nocc
         nmo(4)   = nvir
         mo_list(1:nvir)                    = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)             = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)      = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:3*nvir+nocc) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,vovv_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Clean-up
         ierr=talsh_tensor_destruct(aoint_tensor)

        end subroutine get_CC_integrals


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
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
!        Transform integrals
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)

        end subroutine get_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         use tensor_algebra
         use talsh
       
         implicit none

         type(talsh_tens_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         logical, intent(in)                :: antisymmetrize

         integer:: mulliken_dims(1:4),dirac_dims(1:4)

         integer   :: ierr, mismatch

!        work-around: addition with permutation is not implemented
         type(talsh_tens_t):: one_tensor, minusone_tensor
         integer:: scalar_dims(1)
         scalar_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,scalar_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(minusone_tensor,C8,scalar_dims(1:0),init_val=MINUS_ONE)
!        end work-round

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

!--- Obsolete, remove after CCSD is working ---
        subroutine ccd_driver (nocc,nvir,mo_occ,mo_vir)

!This subroutine drives the ccd energy calculation
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
 
!        Written by Lucas Visscher, winter 2016/2017
 
         use tensor_algebra
         use talsh

         implicit none
         integer, intent(in) :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in) :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer   :: ierr, i
         complex(8), pointer :: energy(:),result_tens(:)
         type(C_PTR):: body_p
         integer(C_SIZE_T), parameter:: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_SIZE_T):: host_buf_size
         integer(C_INT):: host_arg_max

!        two-electron integral tensors
         type(talsh_tens_t):: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor
!        intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t):: goint_tensor, gvint_tensor
         type(talsh_tens_t):: aint_tensor, bint_tensor, hint_tensor
         integer:: goint_dims(1:2),gvint_dims(1:2)
         integer:: aint_dims(1:4),bint_dims(1:4),hint_dims(1:4)
!        solution tensors
         type(talsh_tens_t):: t2_tensor, s2_tensor
         integer:: t2_dims(1:4), s2_dims(1:4)
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t):: result_tensor, one_tensor
         integer:: result_dims(1),one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CCD control and result variables
         integer :: iccd, nccd
         real(8) :: mp2_energy, ccd_energy

!        Make some noise so that we know we are here
         print*, " Arrived in cc_driver routine"

!        Initialize libraries
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         print*, " Initialized talsh library"
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
         call interest_initialize
         print*, " Initialized interest library"

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         result_dims(1) = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(2*nocc))
         allocate (eps_vir(2*nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         
!        Get integral tensors
         call get_CCD_integrals (nocc,nvir,mo_occ,mo_vir,oooo_tensor,oovv_tensor,vvvv_tensor,vovo_tensor)

!        Initialize t2 amplitudes with the MP2 values
         t2_dims(1:2) = nvir * 2
         t2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(t2_tensor,C8,t2_dims,init_val=ZERO)
!        ierr=talsh_tensor_add("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor,one_tensor)
!        end of work-around
         call scale_with_denominators (t2_tensor,eps_occ,eps_vir)

!        Calculate and print the MP2 energy
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor)
         mp2_energy=result_tens(1) * ONE_QUARTER
         print*," MP2 energy = ",mp2_energy

!        Start CCD iterations

         nccd = 30 ! hard wired for the moment

         goint_dims = nocc  * 2
         ierr=talsh_tensor_construct(goint_tensor,C8,goint_dims,init_val=ZERO)
         gvint_dims = nvir * 2
         ierr=talsh_tensor_construct(gvint_tensor,C8,gvint_dims,init_val=ZERO)
         aint_dims  = nocc  * 2
         ierr=talsh_tensor_construct(aint_tensor,C8,aint_dims,init_val=ZERO)
         hint_dims(1) = nvir * 2
         hint_dims(2) = nocc * 2
         hint_dims(3) = nvir * 2
         hint_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(hint_tensor,C8,hint_dims,init_val=ZERO)
         s2_dims(1:2) = nvir * 2
         s2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(s2_tensor,C8,s2_dims,init_val=ZERO)

         do iccd = 1, nccd

!        calculate intermediaries
         ierr=talsh_tensor_init(goint_tensor)
         ierr=talsh_tensor_contract("G(k,j)+=V(k,l,c,d)*T(c,d,l,j)",goint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_tensor_init(gvint_tensor)
         ierr=talsh_tensor_contract("G(a,c)+=V(k,l,c,d)*T(d,a,k,l)",gvint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=talsh_tensor_init(aint_tensor)
         ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",aint_tensor,oooo_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=13; return; endif
         ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",aint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=13; return; endif
         ierr=talsh_tensor_init(hint_tensor)
         ierr=talsh_tensor_add("H(a,i,c,k)+=V(a,i,c,k)",hint_tensor,vovo_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=14; return; endif
         ierr=talsh_tensor_contract("H(a,i,c,k)+=V(k,l,c,d)*T(a,d,l,i)",hint_tensor,oovv_tensor,t2_tensor, &
                                     scale=ONE_HALF)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=14; return; endif

!        For debugging
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=T(a,b,i,j)*T(a,b,i,j)",result_tensor,t2_tensor,t2_tensor, &
                                     scale=ONE_QUARTER)
         print*, "T2amp = ",result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=G(k,j)*G(k,j)",result_tensor,goint_tensor,goint_tensor)
         print*, "goint = ",result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=G(a,c)*G(a,c)",result_tensor,gvint_tensor,gvint_tensor)
         print*, "gvint = ",result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=A(k,l,i,j)*A(k,l,i,j)",result_tensor,aint_tensor,aint_tensor)
         print*, "a int = ",result_tens(1)/4.D0
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("S()+=H(a,i,c,k)*H(a,i,c,k)",result_tensor,hint_tensor,hint_tensor)
         print*, "h int = ",result_tens(1)

!        calculate s2 tensor
         ierr=talsh_tensor_init(s2_tensor)
!        ierr=talsh_tensor_add("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor,one_tensor)
!        end of work-around

!        This terms appears twice as we need to permute i and j
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,goint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,goint_tensor,t2_tensor)

!        This terms appears twice as we need to permute a and b
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,gvint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,gvint_tensor,t2_tensor)

!        The oooo contribution
         ierr=talsh_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,aint_tensor,t2_tensor, &
                                    scale=ONE_HALF)

!        The vvvv contribution
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,vvvv_tensor,t2_tensor, &
                                     scale=ONE_HALF)

!        This terms appears four times as we need to permute both i and j and a and b
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,i,c,k)*T(c,b,j,k)",s2_tensor,hint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,j,c,k)*T(c,b,k,i)",s2_tensor,hint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,i,c,k)*T(a,c,j,k)",s2_tensor,hint_tensor,t2_tensor)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,j,c,k)*T(a,c,k,i)",s2_tensor,hint_tensor,t2_tensor)

!        prepare for next iteration
         ierr=talsh_tensor_init(t2_tensor)
         ierr=talsh_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",t2_tensor,s2_tensor)
         call scale_with_denominators (t2_tensor,eps_occ,eps_vir)

         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor)
         ccd_energy=result_tens(1) * ONE_QUARTER
         print*," CCD energy = ",ccd_energy

         end do

!        Destruct tensors and shut down library
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
        end subroutine ccd_driver

        subroutine get_CCD_integrals (nocc,nvir,mo_occ,mo_vir,oooo_tensor,oovv_tensor,vvvv_tensor,vovo_tensor)

!        Routine to get integrals needed in CCD in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017
 
         use tensor_algebra
         use talsh
         use dirac_exatensor
         use exacorr_ao_to_mo

         implicit none
         integer, intent(in)               :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)               :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(talsh_tens_t), intent(inout) :: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor
         integer:: oovv_dims(1:4),vvvv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4)

         integer   :: ierr, i
         complex(8), pointer, contiguous:: ao_tens(:,:,:,:)
         type(C_PTR):: body_p

         integer             :: nshells, number_aoshells   ! number of AO shells
         type(basis_func_info_t), allocatable :: shell(:)  ! array with basis function information
         integer             :: nao, number_aos

!        auxilliary integral tensors (only needed inside this routine)
         type(talsh_tens_t):: aoint_tensor, ovvo_tensor
         integer:: aoint_dims(1:4), ovvo_dims(1:4)

!        arrays needed to communicate with the MO integral generator
         integer:: nmo(4)
         integer, allocatable :: mo_list(:)
         real(8), allocatable :: eps_occ(:),eps_vir(:)

!        scalars (need to be defined as tensor types)
         type(talsh_tens_t):: minusone_tensor
         integer:: minusone_dims(1)

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)      ! number of basis functions
         nshells = number_aoshells()  ! basis function are organized in shells, this is their number
         allocate(shell(nshells))
         call initialize_shells(shell)
         do i = 1, nshells
            if (compress_unc_shell(shell(i)) == 1) call quit ('This version can not handle contracted functions')
         enddo
         print*, " Retrieved basis set information from DIRAC common blocks"

!        Get tensor with AO integrals
         aoint_dims      = nao
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens,nshells,shell)

!        Get antisymmetric oovv tensor
         oovv_dims(1:2) = nocc * 2
         oovv_dims(3:4) = nvir * 2
         ierr=talsh_tensor_construct(oovv_tensor,C8,oovv_dims,init_val=ZERO)
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

!        Get antisymmetric oooo tensor
         oooo_dims = nocc * 2
         ierr=talsh_tensor_construct(oooo_tensor,C8,oooo_dims,init_val=ZERO)
         allocate (mo_list(4*nocc))
         nmo = nocc 
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,oooo_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         vovo_dims(1) = nvir * 2
         vovo_dims(2) = nocc * 2
         vovo_dims(3) = nvir * 2
         vovo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(vovo_tensor,C8,vovo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,vovo_tensor,.false.,nmo,mo_list)
         deallocate (mo_list)
!        Get ovvo tensor
         ovvo_dims(1) = nocc * 2
         ovvo_dims(2) = nvir * 2
         ovvo_dims(3) = nvir * 2
         ovvo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nvir
         nmo(4) = nocc
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nvir+nocc+1:2*nvir+nocc)     = mo_vir(1:nvir)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,ovvo_tensor,.false.,nmo,mo_list)
         deallocate (mo_list)
!        Make anti-symmetrized vovo tensor
         minusone_dims(1) = 1
         ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",vovo_tensor,ovvo_tensor,minusone_tensor)

!        Get vvvv tensor
         vvvv_dims = nvir * 2
         ierr=talsh_tensor_construct(vvvv_tensor,C8,vvvv_dims,init_val=ZERO)
         allocate (mo_list(4*nvir))
         nmo = nvir
         mo_list(1:nvir)          = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
         mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
         mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,vvvv_tensor,.true.,nmo,mo_list)
         deallocate (mo_list)

!        Clean-up
         ierr=talsh_tensor_destruct(aoint_tensor)
         ierr=talsh_tensor_destruct(ovvo_tensor)

        end subroutine get_CCD_integrals

       end module exacorr_cc
