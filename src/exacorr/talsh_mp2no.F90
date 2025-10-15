module talsh_mp2no
!This module contains the routines needed to calculate the frozen mp2 natural orbitals

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils
        use talsh_common_routines
        use exacorr_mo
        use exacorr_global
        use mp2no_matrixtransforms
        use mp2no_fileinterfaces
        use talsh_cis
        use talsh_ccrsp_common
        use, intrinsic:: ISO_C_BINDING

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
        real(8), parameter    :: ONE_QUARTER=0.25D0

        type(talsh_tens_t) :: one_tensor

        private

        public talsh_mp2no_driver, get_new_orbital, get_transform_Fock

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine talsh_mp2no_driver (exa_input,n_cano)
!       This module is used to construct the MP2 density matrix by talsh.
!       The complex density matrix are saved on the file 'DM_complex'.
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new reduced canonical orbitals are written into 'MP2NOs_AO' by Module3.
!       Written by Xiang Yuan, modified by Johann Pototschnig

        use talsh_mp2lap

        type(exacc_input), intent(in) :: exa_input
        integer, intent(inout)        :: N_cano                  ! number of new canonical orbital


!       Point connected with tensors
        integer                :: ierr
        real(8),pointer        :: NOs_AO(:,:,:)           ! new quaternion orbital in AO basis
        real(kind=8), pointer  :: N_OCCUPATION(:)         ! occupation number from qdia of DM
        real(kind=8), pointer  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
        logical                :: alive
        real(kind=8)           :: threshold
        integer                :: i, nocc
        character(len=10)      :: fn

        call print_date("Entering MP2_NO Module")

        threshold = exa_input%t_mp2no
        !LV: the line below should be improved, will not work if a user does not give the HOMO as the last orbital
        nocc = exa_input%mokr_occ(size(exa_input%mokr_occ))

        inquire(file="NOs_MO", exist=alive)

        if (alive) then
            call print_date('File NOs_MO exists')
        else
            call print_date('File NOs_MO does not exist. Start constructing MP2 density matrix')
            if (exa_input%mp2no<9) then
                  call talsh_generate_MP2_VNO_MObasis(exa_input)
            else
                  call mp2lap_get_DM(exa_input)
            endif
            !call get_transform_Fock(exa_input,'DM_complex')
            call get_transform_Fock(exa_input)
        end if

        call get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION, &
                              New_orbital_energy,'NOs_MO',N_cano)
        call write_nos(nocc,NOs_AO,New_orbital_energy,'scf','mbpt')

        call print_date('Leaving MP2_NO Module')
      end subroutine talsh_mp2no_driver


!------------------------------------------------------------------------------------------------------------------------------------
      subroutine talsh_generate_MP2_VNO_MObasis (exa_input)
!       This subroutine obtains the complex mp2 density matrix in MO basis
!       The non-diagonal complex DM is written into file DM_complex.
!       Written by Xiang Yuan

         use checkpoint

         type(exacc_input), intent(in)          :: exa_input
         integer                                :: nocc, nvir              ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable                  :: mo_vir(:)               ! the list index of virtual spinors
         integer,  allocatable                  :: mo_occ(:)               ! the list index of occupied spinors
         real(8), allocatable                   :: eps_occ(:),eps_vir(:)   ! HF orbital energies

!        Intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t)          :: MP2_vv_density, MP2_vv_density_tmp
         type(talsh_tens_t)          :: CIS_vv_effden
         type(talsh_tens_t)          :: result_tensor
         integer(INTD), dimension(4) :: vvoo_dims, oovv_dims
         integer(INTD), dimension(2) :: vv_dims, ov_dims, oo_dims, vo_dims
         integer(C_INT)              :: one_dims(1)

!        Common tensors
         type(talsh_comm_tens_t) :: comm_t

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens_t)       :: int_t

!        property integrals for PNO
         type(talsh_orb_prop_integ_tens_t)   :: orb_int_tensor1
         type(talsh_rsp_interm_tens_t)       :: interm_t_prop1
         !type(talsh_tens_shape_t)            :: tens_shape
         integer(INTD), dimension(:), pointer:: tens_shape_ptr
         complex(8), pointer   :: hbar_oo_ptr(:,:), hbar_vv_ptr(:,:)
         complex(8), pointer   :: t1_ptr(:,:), t2_ptr(:,:,:,:)
         complex(8), pointer   :: MP2_vv_density_array_ptr(:,:)
         complex(8), pointer   :: MP2_vv_density_array_tmp(:,:)
         !integer(C_SIZE_T)           :: volume
         !integer(C_INT)           :: rankk

!        Point connected with tensors
         integer(C_INT)              :: ierr
         complex(8), pointer         :: result_tens(:)
         integer                     :: initialize_talsh
         real(8)                     :: mp2_energy,scf_energy

         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max

!        Scalars (need to be defined as tensor types)
         real(kind=8)  :: energy_shift(2) = 0.0d0
         integer       :: i, j, k, a, c

         integer :: num_gpus = 0, igpu
         integer :: i_mp2fno_prop 
         character(len=20)                  :: fn, istr

         logical :: fock_flag ! used for FNO++ orbitals

         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir


!        Make some noise so that we know we are here
         call print_date('Entered mp2_driver routine')

         call talsh_get_number_gpus(num_gpus)

!        Initialize libraries
         initialize_talsh=1
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max,gpu_list=(/(igpu,igpu=0,num_gpus-1)/))
         call print_date('Initialized talsh library')
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

         call print_date('Initializing interest library')
         call interest_initialize(.true.)
         call print_date('Initialized interest library')

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         vv_dims(1) = nvir
         vv_dims(2) = nvir

         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(MP2_vv_density,C8,vv_dims,init_val=ZERO)

!        Get HF orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("retrieved orbital energies")

!        Get oovv integrals
         fock_flag = exa_input%mp2no.eq.3
         call get_MP2_integrals(nocc,nvir,mo_occ,mo_vir,eps_occ,eps_vir,int_t,exa_input%print_level, fock_flag)

!        Make some noise so that we know we are proceeding
         call print_date('started forming VV block of MP2 density matrix')
!        Initialize scale with denominators of oovv integrals
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc

         ierr=talsh_tensor_construct(comm_t%t2,C8,vvoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,one_dims)
!        Initialize t2 amplitudes with the MP2 values
         ierr=talsh_tensor_init(comm_t%t2)
         ierr=talsh_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",comm_t%t2,int_t%oovv,one_tensor)
         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
!        Calculate and print the MP2 energy
!        Note: here we don't consider the t1 amplitude that is different with excorr_mp2no case.
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv) ! o(2)v(2)
         mp2_energy=result_tens(1)*ONE_QUARTER
         write(*,*) ""
         write(*,*) "    MP2 energy = ", mp2_energy
         write(*,*) ""

         call checkpoint_write ('/result/wavefunctions/mbpt/method',sdata="MP2")
         call checkpoint_write ('/result/wavefunctions/mbpt/e_corr',rdata=mp2_energy)
         call checkpoint_read ('/result/wavefunctions/scf/energy',rdata=scf_energy)
         call checkpoint_write ('/result/wavefunctions/mbpt/energy',rdata=scf_energy+mp2_energy)
         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=int_t%oovv)

!        Calculate and print the MP2 density matrix of virtual-virtual block Dvv
         ierr=talsh_tensor_init(MP2_vv_density)
         if(exa_input%mp2no.eq.3) then

            write(*,*) ""
            call print_date("Computing MP2FNO++ density matrix..")
            write(*,*) "omega = ", exa_input%omega

            ! FNO++
            ! Notations from:
            ! Reduced-scaling coupled cluster response theory: Challenges and opportunities (2018)
            !
            !---------------
            ! constructing tensors
            !---------------
            ! Hbar in comm_t%hoo and comm_t%hvv (not to confuse with actual hoo and hvv)
            oo_dims = nocc

            ierr=talsh_tensor_construct(comm_t%hoo,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_add("H(i,j)+=F(i,j)",comm_t%hoo,int_t%foo)
            ierr=talsh_tensor_contract("H(i,j)+=T(e,f,i,n)*V(j,n,e,f)",comm_t%hoo,comm_t%t2,int_t%oovv,scale=ONE_HALF)

            ierr=talsh_tensor_construct(comm_t%hvv,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_add("H(a,b)+=F(a,b)",comm_t%hvv,int_t%fvv)
            ierr=talsh_tensor_contract("H(a,b)+=T(f,a,m,n)*V(m,n,f,b)",comm_t%hvv,comm_t%t2,int_t%oovv,scale=MINUS_ONE_HALF)

            ierr=talsh_tensor_construct(comm_t%t1,C8,vo_dims,init_val=ZERO)
            do i_mp2fno_prop=1,exa_input%n_mp2fno_prop
               write(*,*) "Computing density matrix of operator ", i_mp2fno_prop
               ! Abar ai
               call get_prop_integral(exa_input,orb_int_tensor1,exa_input%mp2fno_prop(i_mp2fno_prop))
               interm_t_prop1%property_label = orb_int_tensor1%property_label

               ov_dims(1) = nocc
               ov_dims(2) = nvir
               call get_prop_bar_ov(interm_t_prop1%prop_bar_ov,orb_int_tensor1%prop_ov,comm_t%t2,ov_dims)

               ! Abar acij
               oovv_dims(1:2) = nocc
               oovv_dims(3:4) = nvir
               call get_prop_bar_oovv(interm_t_prop1%prop_bar_oovv,orb_int_tensor1,comm_t%t2,oovv_dims)

               !---------------
               ! approximated amplitudes
               !---------------
               ! t1
               vo_dims(1) = nvir
               vo_dims(2) = nocc

               ierr=talsh_tensor_scale(comm_t%t1,val=ZERO)
               ierr=talsh_tensor_add("T(a,i)+=A(i,a)",comm_t%t1,interm_t_prop1%prop_bar_ov)

               ierr=talsh_tensor_get_body_access(comm_t%t1,body_p,C8,int(0,C_INT),DEV_HOST)
               call c_f_pointer(body_p,t1_ptr,vo_dims)
               ierr=talsh_tensor_get_body_access(comm_t%hoo,body_p,C8,int(0,C_INT),DEV_HOST)
               call c_f_pointer(body_p,hbar_oo_ptr,oo_dims)
               ierr=talsh_tensor_get_body_access(comm_t%hvv,body_p,C8,int(0,C_INT),DEV_HOST)
               call c_f_pointer(body_p,hbar_vv_ptr,vv_dims)

               do i=1,nocc
                  do a=1,nvir
                      t1_ptr(a,i) = t1_ptr(a,i)/(hbar_vv_ptr(a,a) - hbar_oo_ptr(i,i) + exa_input%omega)
                  enddo
               enddo

               ! t2
               ierr=talsh_tensor_scale(comm_t%t2,val=ZERO)
               ierr=talsh_tensor_add("T(a,c,i,j)+=A(i,j,a,c)",comm_t%t2,interm_t_prop1%prop_bar_oovv)
               ierr=talsh_tensor_get_body_access(comm_t%t2,body_p,C8,int(0,C_INT),DEV_HOST)
               call c_f_pointer(body_p,t2_ptr,vvoo_dims)

               do i=1,nocc
                  do j=1,nocc
                     do a=1,nvir
                        do c=1,nvir
                           t2_ptr(a,c,i,j) = t2_ptr(a,c,i,j)/(hbar_vv_ptr(a,a) + hbar_vv_ptr(c,c)&
                              - hbar_oo_ptr(i,i) - hbar_oo_ptr(j,j) + exa_input%omega)
                        enddo
                     enddo
                  enddo
               enddo

               !---------------
               ! MP2 density matrix
               !---------------
               ierr=talsh_tensor_contract("D(a,b)+=T(a,c,i,j)*T(b,c,i,j)",MP2_vv_density,&
                                           comm_t%t2, comm_t%t2,scale=ONE_HALF)
               ierr=talsh_tensor_contract("D(a,b)+=T(a,i)*T(b,i)",MP2_vv_density,&
                                           comm_t%t1, comm_t%t1)
               call write_talsh_matrix(nvir,MP2_vv_density,&
                    trim(get_filename_with_index('DM_complex',i_mp2fno_prop)))

               ierr=talsh_tensor_scale(MP2_vv_density,val=ZERO)
            enddo

            if(exa_input%do_sum_mp2fno_prop_density) then
               ! initialize arrays
               ierr=talsh_tensor_construct(MP2_vv_density_tmp,C8,vv_dims,init_val=ZERO)
               allocate(MP2_vv_density_array_tmp(nvir,nvir))

               do i_mp2fno_prop=1,exa_input%n_mp2fno_prop
                  ! read the matrix
                  call read_matrix_comp(MP2_vv_density_array_tmp,&
                      trim(get_filename_with_index('DM_complex',i_mp2fno_prop)))

                  ! assign the values of an array to a tensor
                  ierr=talsh_tensor_get_body_access(MP2_vv_density_tmp,body_p,C8,int(0,C_INT),DEV_HOST)
                  call c_f_pointer(body_p,MP2_vv_density_array_ptr,vv_dims)
                  do a=1,nvir
                     do c=1,nvir
                        MP2_vv_density_array_ptr(a,c) = MP2_vv_density_array_tmp(a,c)
                     enddo
                  enddo

                  ! contract
                  ierr=talsh_tensor_add("D(a,b)+=D_tmp(a,b)",MP2_vv_density,MP2_vv_density_tmp,&
                      scale=dcmplx(exa_input%mp2fno_prop_weights(i_mp2fno_prop),0.0))
               enddo
               call write_talsh_matrix(nvir,MP2_vv_density,'DM_complex')

               ! destroy arrays
               deallocate(MP2_vv_density_array_tmp)
               ierr=talsh_tensor_destruct(MP2_vv_density_tmp)
            endif

            call print_date("Finished MP2FNO++ density matrix!")
            return
         else
            ierr=talsh_tensor_contract("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_vv_density,&
                                        int_t%oovv,int_t%oovv,scale=ONE_HALF) ! o(2)v(2)
         endif

!        Make some noise so that we know we are proceeding
         call print_date('finished MP2_DM_VV calculation and save it to file DM_complex')

         if (exa_input%do_neo) then
            call generate_effec_densive(exa_input,MP2_vv_density,CIS_vv_effden)
         end if

         call write_talsh_matrix(nvir,MP2_vv_density,'DM_complex')

        end subroutine talsh_generate_MP2_VNO_MObasis

!------------------------------------------------------------------------------------------------------------------------------------
    subroutine get_MP2_integrals (nocc,nvir,mo_occ,mo_vir,eps_occ,eps_vir,int_t,print_level,fock_flag)

!    Routine to get oovv integral in the form of antisymmetric tensors.
!    In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!    full integral list.
!    copy from talsh_cc.F90 routine: get_CC_integrals

      use talsh_ao_to_mo
      use talsh_cc, only : init_fock_talsh

      integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
      integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
      real(8), intent(in)    :: eps_occ(:),eps_vir(:)
!     target integral tensors to be filled and given back to the caller
      type(talsh_intg_tens_t),intent(inout)  :: int_t
      integer, intent(in)                  :: print_level
      logical                              :: fock_flag
      integer(INTD)                        :: oovv_dims(1:4)
      integer(INTD)                        :: ierr
      complex(8), pointer, contiguous      :: ao_tens(:,:,:,:)
      type(C_PTR)                          :: body_p
      integer                              :: nao

!     auxilliary integral tensors (only needed inside this routine)
      type(talsh_tens_t) :: aoint_tensor, ovvo_tensor
      complex(8), pointer, contiguous      :: h_tens(:,:)
      integer(INTD)      :: aoint_dims(1:4), ovvo_dims(1:4), oooo_dims(1:4), vovo_dims(1:4)
      integer(INTD), dimension(2) :: vv_dims, oo_dims, vo_dims

!     arrays needed to communicate with the MO integral generator
      integer:: nmo(4)
      integer, allocatable :: mo_list(:)

!     scalars (need to be defined as tensor types)
      type(talsh_tens_t) :: minusone_tensor
      integer(INTD)      :: minusone_dims(1)

      real(8) :: scf_energy

      call print_date('entering integral transformation, scheme = 2')

!     Retrieve basis set information
      nao = get_nao()      ! number of basis functions

!     Get tensor with AO integrals
      aoint_dims = nao
      ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
      if (ierr /= 0) then
            print*, "Could not construct ao integral tensor, increasing TALSH_BUFF may help"
            call quit ('error constructing aoint_tensor')
      end if
      call print_date('allocated AO integral tensor')
!     Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
      ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
      call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
      call print_date('finished ao integral calculation')

!     Get antisymmetric oovv tensor
!     Note that the dimensions for the tensor are given in physicists notation <12||34>
      oovv_dims(1:2) = nocc
      oovv_dims(3:4) = nvir
      ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
!     whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
      allocate (mo_list(2*nocc+2*nvir))
      nmo(1) = nocc
      nmo(2) = nvir
      nmo(3) = nocc
      nmo(4) = nvir
      mo_list(1:nocc)                      = mo_occ(1:nocc)
      mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
      mo_list(nocc+nvir+1:2*nocc+nvir)     = mo_occ(1:nocc)
      mo_list(2*nocc+nvir+1:2*nocc+2*nvir) = mo_vir(1:nvir)
      call get_integral_tensor (aoint_tensor,int_t%oovv,12,nmo,mo_list)
      deallocate (mo_list)
      call print_date('finished oovv integrals')

      if(fock_flag) then
         ! construct and init tensors
         oo_dims = nocc
         vv_dims = nvir
         oooo_dims = nocc
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ovvo_dims(1) = nocc
         ovvo_dims(4) = nocc
         ovvo_dims(2:3) = nvir
         ierr=talsh_tensor_construct(int_t%foo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(int_t%fvv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(int_t%oooo,C8,oooo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(int_t%vovo,C8,vovo_dims,init_val=ZERO)

         ! compute integrals
         allocate (mo_list(4*nocc))
         nmo = nocc
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,int_t%oooo,12,nmo,mo_list)
         deallocate (mo_list)
         call print_date('finished oooo integrals')

         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,int_t%vovo,0,nmo,mo_list)
         deallocate (mo_list)
!        Get ovvo tensor
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
         call get_integral_tensor (aoint_tensor,ovvo_tensor,0,nmo,mo_list)
         deallocate (mo_list)
!        Make anti-symmetrized vovo tensor
         minusone_dims(1) = 1
         ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,minusone_tensor)
         ierr=talsh_tensor_destruct(ovvo_tensor)
         call print_date('finished vovo integrals')

         ! f integrals
         call init_fock_talsh(2, mo_occ, mo_vir, eps_occ, eps_vir, int_t, scf_energy, print_level)
         call print_date('Finished foo and fvv integrals!')

      endif
   end subroutine get_MP2_integrals

!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_Fock (exa_input)
         ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
         ! The new basis set comes from the digonalization of density matrix.
         ! The all natural orbitals and occupation numbers are saved in file NOs_MO
         ! The transformed Fock matrix are saved in file FM_in_NatOrb
         ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix.
         ! This routine doesn't truncate the new basis set.
         ! Written by Xiang Yuan


                  !character(len=*),intent(in)        :: Matrix_file
                  type(exacc_input), intent(in)      :: exa_input

                  type(cmo)                          :: qorbital
                  type(basis_set_info_t)             :: ao_basis
                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  real(kind=8), pointer              :: N_OCCUPATION_tmp(:)     ! occupation number from qdia of DM
                  complex(kind=8), pointer           :: MP2_density_comp(:,:)   ! MP2 density matrix from DM_complex
                  real(8), pointer                   :: MP2_density_quat(:,:,:) ! quaternion MP2 density matrix
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: NOs_MO_tmp(:,:,:)       ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix
                  real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
                  !integer                            :: sz,flb,sz_occ           ! size of matrix and flag for block
                  integer                            :: sz1, sz2                ! size of matrix and flag for block
                  real(kind=8)                       :: Sum_N_OCCUPATION_total, Sum_N_OCCUPATION_total_tmp
                  integer                            :: i, ierr, nocc, nvirt
                  integer                            :: i_mp2fno_prop, counter
                  integer                            :: N_reduced, N_reduced_tmp
                  character(len=20)                  :: fn, istr

                  if((exa_input%mp2no==3).and.exa_input%do_sum_mp2fno_prop_density) then
                      ! Generating the FNO and save them in files
                      sz1 = 0
                      sz2 = 0
                      do i_mp2fno_prop = 1,exa_input%n_mp2fno_prop
                         call print_date('Diagonalizing perturbed density...')
                         call read_matrix_comp(MP2_density_comp,&
                             trim(get_filename_with_index('DM_complex',i_mp2fno_prop)))
                         call print_date('read density matrix (DM) in MO basis (spinor) from DM_complex')

                         call complex_to_quaternion(MP2_density_comp,MP2_density_quat)
                         call print_date('changed MO DM representation (spinor->quaternion)')

                         call diagonal_matrix_quater(MP2_density_quat,N_OCCUPATION_tmp,NOs_MO_tmp)
                         call print_date('diagonalized MO DM (quaternion)')

                         call quater_write(NOs_MO_tmp,N_OCCUPATION_tmp,&
                             trim(get_filename_with_index('NOs_MO',i_mp2fno_prop)))
                         call print_date('wrote all natural orbital and occupation number in file NOs_MO')

                         sz1 = sz1 + size(NOs_MO_tmp,1)
                         sz2 = sz2 + size(NOs_MO_tmp,2)

                         deallocate(MP2_density_quat) ! allocation done in complex_to_quaternion
                         deallocate(NOs_MO_tmp)       ! allocation done in diagonal_matrix_quater
                         deallocate(N_OCCUPATION_tmp) ! allocation done in diagonal_matrix_quater
                      enddo

                      ! Read the FNO from files and concatenate them in arrays 
                      allocate(NOs_MO(sz1,sz2,4))
                      allocate(N_OCCUPATION(sz2))

                      Sum_N_OCCUPATION_total = 0.0
                      N_reduced = 0
                      counter = 0
                      do i_mp2fno_prop = 1,exa_input%n_mp2fno_prop
                         call quater_read_orbital(NOs_MO_tmp,N_OCCUPATION_tmp,Sum_N_OCCUPATION_total_tmp,exa_input%t_mp2no,N_reduced_tmp,&
                             trim(get_filename_with_index('NOs_MO',i_mp2fno_prop)))
                         
                         forall (i=1:size(N_OCCUPATION_tmp)) N_OCCUPATION(i+counter) = N_OCCUPATION_tmp(i)
                         do i = 1,size(N_OCCUPATION_tmp)
                            NOs_MO(:,i+counter,:) = NOs_MO_tmp(:,i,:)
                         enddo

                         Sum_N_OCCUPATION_total = Sum_N_OCCUPATION_total + Sum_N_OCCUPATION_total_tmp
                         N_reduced = N_reduced + N_reduced_tmp
                         counter = counter + size(N_OCCUPATION_tmp)

                         deallocate(NOs_MO_tmp)       ! allocation done in quater_read_orbital
                         deallocate(N_OCCUPATION_tmp) ! allocation done in quater_read_orbital
                      enddo

                      call print_mp2fvno_occupation_number_details(exa_input%t_mp2no, N_reduced, N_OCCUPATION, Sum_N_OCCUPATION_total)
                  else
                     call read_matrix_comp(MP2_density_comp,'DM_complex')
                     call print_date('read density matrix (DM) in MO basis (spinor) from DM_complex')

                     call complex_to_quaternion(MP2_density_comp,MP2_density_quat)
                     call print_date('changed MO DM representation (spinor->quaternion)')

                     call diagonal_matrix_quater(MP2_density_quat,N_OCCUPATION,NOs_MO)
                     call print_date('diagonalized MO DM (quaternion)')
                  endif

                  call sort_values_vectors (N_OCCUPATION,NOs_MO)
                  call print_date('Sorted natural orbitals by occupation number')

                  if(exa_input%print_occupation_number.eq.1) then
                      print*, "-------------------------------"
                      print*, "Printing all occupation numbers"
                      print*, "-------------------------------"
                      do i=1,size(N_OCCUPATION)
                         write(*,'(I6, A2, E15.4)') i,') ', N_OCCUPATION(i)
                      enddo
                  endif

                  call quater_write(NOs_MO,N_OCCUPATION,'NOs_MO')
                  call print_date('wrote all natural orbital and occupation number in file NOs_MO')

                  ! we initialize the (diagonal) vv block of the Fock matrix with the orbital energies
                  call read_from_dirac(ao_basis,qorbital,'scf',ierr)
                  call convert_mo_to_quaternion(qorbital)
                  nvirt = exa_input%nkr_vir        ! the number of virtual orbitals
                  nocc  = exa_input%mokr_occ(size(exa_input%mokr_occ))
                  allocate(Fock_HF(nvirt,nvirt,4))
                  Fock_HF = 0.0
                  forall (i=1:nvirt) Fock_HF(i,i,1) = qorbital%energy(nocc+i)
                  call print_date('read Fock matrix from checkpoint')

                  call basis_transform_matrix(Fock_HF,NOs_MO,Fock_NO)
                  call print_date('got full VV Fock matrix in natural orbital basis F_VV_NO')

                  call quater_write(Fock_NO,'FM_inNatOrb')
                  call print_date('wrote F_VV_NO in file FM_in_NatOrb')

      end subroutine get_transform_Fock


!------------------------------------------------------------------------------------------------------------------------------
      subroutine get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION_reduced,New_orbital_energy,file_orb,N_reduced)
      ! This routine is used for getting new truncated canonical orbital in AO basis
      ! and the corresponding orbital energy.


            type(exacc_input), intent(in)      :: exa_input
            real(kind=8),intent(in)            :: threshold
            real(kind=8), pointer,intent(out)  :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis
            real(kind=8), pointer,intent(out)  :: N_OCCUPATION_reduced(:) ! occupation number from qdia of DM
            real(kind=8), pointer,intent(out)  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
            character(len=*),intent(in)        :: file_orb
            integer,intent(inout)              :: N_reduced               ! number of reduced natural orbital

            type(cmo)                          :: qorbital
            type(basis_set_info_t)             :: ao_basis
            real(kind=8)                       :: Sum_N_OCCUPATION_total
            real(kind=8)                       :: Sum_N_OCCUPATION_reduced
            real(kind=8)                       :: Perc_OCCUPATION_reduced
            real(kind=8), pointer              :: NOs_MO_reduced(:,:,:)   ! natural orbitals which have occupation number larger than threshold
            real(kind=8), pointer              :: New_cano_MO(:,:,:)      ! new canonical orbital in AO basis
            real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
            real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix
            real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital
            integer                            :: nocc, nvirt             ! number of occupied and virtuals
            integer                            :: ierr,i,j,k

            call quater_read(NOs_MO_reduced,N_OCCUPATION_reduced,Sum_N_OCCUPATION_total,&
                              threshold,N_reduced,file_orb)
            call print_date('got the truncated NO and occupation number from NOs_MO')

            nvirt = exa_input%nkr_vir        ! the number of virtual orbitals
            nocc  = exa_input%mokr_occ(size(exa_input%mokr_occ))
            call read_from_dirac(ao_basis,qorbital,'scf',ierr)
            call convert_mo_to_quaternion(qorbital)
            allocate(HF_AO(qorbital%nao,nvirt,4))
            HF_AO = qorbital%coeff_q(:,nocc+1:nocc+nvirt,:)
            call print_date('read HF orbital coefficient from checkpoint')

            ! we initialize the (diagonal) vv block of the Fock matrix with the orbital energies
            allocate(Fock_HF(nvirt,nvirt,4))
            Fock_HF = 0.0
            forall (i=1:nvirt) Fock_HF(i,i,1) = qorbital%energy(nocc+i)
            call print_date('read Fock matrix from checkpoint')

            call basis_transform_matrix(Fock_HF,NOs_MO_reduced,Fock_NO)
            call print_date('transformed Fock matrix into reduced NO basis')

            call get_Cano_Orbital(NOs_MO_reduced,New_cano_MO,Fock_NO,New_orbital_energy)
            call print_date('got new Canonical orbital in MO basis')

            call basis_transform_vector(New_cano_MO,HF_AO,NOs_AO)
            call print_date('got new orbitals in AO basis')

            call print_rMP2FVNO(New_orbital_energy,N_OCCUPATION_reduced)
      end subroutine

!------------------------------------------------------------------------------------------------------------------------------
      subroutine generate_effec_densive(exa_input,MP2_vv_density,CIS_vv_effden)

            type(exacc_input), intent(in) :: exa_input
            type(talsh_tens_t)          :: MP2_vv_density
            type(talsh_tens_t)          :: CIS_vv_effden
            integer(INTD)               :: ierr

            call print_date('Start computing CIS effective density matrix')
            call talsh_cis_driver(exa_input,CIS_vv_effden,initialize_tensor_library=.false.,initialize_interest=.false.)
            ierr=talsh_tensor_add('D(a,b)+=DC(a,b)',MP2_vv_density,CIS_vv_effden,scale=ONE)
            call print_date('Add CIS effective density matrix into MP2 desnity matrix')

      end subroutine

      subroutine get_prop_bar_ov(prop_bar_ov,prop_ov,t2,ov_dims)
         type(talsh_tens_t), intent(inout)       :: prop_bar_ov, prop_ov, t2
         integer(INTD), dimension(2), intent(in) :: ov_dims
         integer(INTD)                           :: ierr

         ierr=talsh_tensor_construct(prop_bar_ov,C8,ov_dims,init_val=ZERO)
         ierr=talsh_tensor_add("A(i,a)+=B(i,a)",prop_bar_ov,prop_ov)
         ierr=talsh_tensor_contract("A(i,a)+=T(a,b,i,j)*B(j,b)",prop_bar_ov,t2,prop_ov)
      end subroutine

      subroutine get_prop_bar_oovv(prop_bar_oovv,orb_int_tensor1,t2,oovv_dims)
         type(talsh_tens_t), intent(inout)               :: prop_bar_oovv, t2
         type(talsh_orb_prop_integ_tens_t), intent(inout):: orb_int_tensor1
         integer(INTD), dimension(4), intent(in)         :: oovv_dims
         integer(INTD)                                   :: ierr

         ierr=talsh_tensor_construct(prop_bar_oovv,C8,oovv_dims,init_val=ZERO)

         ! First term with P(a,c)
         ierr=talsh_tensor_contract("A(i,j,a,c)+=T(e,c,i,j)*B(a,e)",prop_bar_oovv,t2,orb_int_tensor1%prop_vv)
         ierr=talsh_tensor_contract("A(i,j,a,c)+=T(e,a,i,j)*B(c,e)",prop_bar_oovv,t2,orb_int_tensor1%prop_vv,scale=MINUS_ONE)

         ! Second term P(i,j)
         ierr=talsh_tensor_contract("A(i,j,a,c)+=T(a,c,m,j)*B(m,i)",prop_bar_oovv,t2,orb_int_tensor1%prop_oo,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("A(i,j,a,c)+=T(a,c,m,i)*B(m,j)",prop_bar_oovv,t2,orb_int_tensor1%prop_oo)
      end subroutine

end module talsh_mp2no

