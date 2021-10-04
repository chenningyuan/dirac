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

        public talsh_mp2no_driver, get_new_orbital, get_transform_Fock, write_talsh_matrix

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine talsh_mp2no_driver (exa_input,n_cano)
!       This module is used to construct the MP2 density matrix by talsh.
!       The complex density matrix are saved on the file 'DM_complex'.  
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new reduced canonical orbitals are written into 'MP2NOs_AO' by Module3.
!       Written by Xiang Yuan

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

        call print_date("Entering MP2_NO Module")
 
        threshold = exa_input%t_mp2no
        !LV: the line below should be improved, will not work if a user does not give the HOMO as the last orbital
        nocc = exa_input%mokr_occ(size(exa_input%mokr_occ))

        inquire(file="NOs_MO", exist=alive)
        if (alive) then
            call print_date('File NOs_MO exists')
        else
            call print_date('File NOs_MO does not exist. Start constructing MP2 density matrix')
            call talsh_generate_MP2_VNO_MObasis(exa_input)   
            call get_transform_Fock(exa_input,'DM_complex')
        end if
        call get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION, &
                              New_orbital_energy,'NOs_MO',N_cano)
        call write_to_newfile(nocc,NOs_AO,New_orbital_energy,'DFCOEF','MP2NOs_AO')

        call print_date('Leaving MP2_NO Module')
      end subroutine talsh_mp2no_driver


!------------------------------------------------------------------------------------------------------------------------------------
      subroutine talsh_generate_MP2_VNO_MObasis (exa_input)
!       This subroutine obtains the complex mp2 density matrix in MO basis
!       The non-diagonal complex DM is written into file DM_complex.
!       Written by Xiang Yuan

         type(exacc_input), intent(in)          :: exa_input
         integer                                :: nocc, nvir              ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable                  :: mo_vir(:)               ! the list index of virtual spinors
         integer,  allocatable                  :: mo_occ(:)               ! the list index of occupied spinors
         real(8), allocatable                   :: eps_occ(:),eps_vir(:)   ! HF orbital energies

!        Intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t)          :: MP2_vv_density
         type(talsh_tens_t)          :: result_tensor
         integer(INTD), dimension(4) :: vvoo_dims
         integer(INTD), dimension(2) :: vv_dims
         integer(C_INT)              :: one_dims(1)

!        Common tensors
         type(talsh_comm_tens_t) :: comm_t

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens_t)       :: int_t

!        Point connected with tensors 
         integer(C_INT)              :: ierr
         complex(8), pointer         :: result_tens(:)
         integer                     :: initialize_talsh
         real(8)                     :: mp2_energy

         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max
         
!        Scalars (need to be defined as tensor types)
         real(kind=8)  :: energy_shift(2) = 0.0d0
         integer       :: i, j, k

         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir


!        Make some noise so that we know we are here
         call print_date('Entered mp2_driver routine')

!        Initialize libraries
         initialize_talsh=1
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)
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
         call get_MP2_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%print_level)

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

         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=int_t%oovv)

!        Calculate and print the MP2 density matrix of virtual-virtual block Dvv 
         ierr=talsh_tensor_init(MP2_vv_density)
         ierr=talsh_tensor_contract("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_vv_density,&
                                        int_t%oovv,int_t%oovv,scale=ONE_HALF) ! o(2)v(2)

!        Make some noise so that we know we are proceeding 
         call print_date('finished MP2_DM_VV calculation and save it to file DM_complex')
         call write_talsh_matrix(nvir,MP2_vv_density,'DM_complex')

        end subroutine talsh_generate_MP2_VNO_MObasis

!------------------------------------------------------------------------------------------------------------------------------------
    subroutine get_MP2_integrals (nocc,nvir,mo_occ,mo_vir,int_t,print_level)
                
!    Routine to get oovv integral in the form of antisymmetric tensors.
!    In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!    full integral list.
!    copy from talsh_cc.F90 routine: get_CC_integrals
         
      use talsh_ao_to_mo
         
      integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
      integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!     target integral tensors to be filled and given back to the caller
      type(talsh_intg_tens_t),intent(inout)  :: int_t
      integer, intent(in)                  :: print_level 
      integer(INTD)                        :: oovv_dims(1:4)
      integer(INTD)                        :: ierr
      complex(8), pointer, contiguous      :: ao_tens(:,:,:,:)
      type(C_PTR)                          :: body_p 
      integer                              :: nao
         
!     auxilliary integral tensors (only needed inside this routine)
      type(talsh_tens_t) :: aoint_tensor, ovvo_tensor
      integer(INTD)      :: aoint_dims(1:4), ovvo_dims(1:4)
         
!     arrays needed to communicate with the MO integral generator
      integer:: nmo(4)
      integer, allocatable :: mo_list(:)
         
!     scalars (need to be defined as tensor types)
      type(talsh_tens_t) :: minusone_tensor
      integer(INTD)      :: minusone_dims(1)
         
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
         
   end subroutine get_MP2_integrals
        
!------------------------------------------------------------------------------------------------------------------------------------
      subroutine write_talsh_matrix(nvir,dm_vv,file_out)
! Routine for writing complex density matrix into the file'DM_complex'. 


! input variables
         integer,intent(in)                  :: nvir
         type(talsh_tens_t), intent(inout)   :: dm_vv
         character(len=*),intent(in)         :: file_out

! body pointers
         type(C_PTR)                     :: body_p2
         complex(8), pointer, contiguous :: dm_block(:,:)
         complex(8), allocatable         :: dm_mo_fort(:,:)
! error code
         integer :: ierr
         integer :: mp2dm

         allocate (dm_mo_fort(nvir,nvir))
         dm_mo_fort = ZERO
         
         ierr=talsh_tensor_get_body_access(dm_vv,body_p2,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p2,dm_block,(/nvir,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
         dm_mo_fort(1:nvir,1:nvir) = dm_block
         
         call get_free_fileunit(mp2dm)
         open (mp2dm,file=file_out,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
         ! write dimension of the matrix 
         write(mp2dm) nvir
         ! write matrix element 
         write(mp2dm) dm_mo_fort

         close(mp2dm,status='KEEP')

        end subroutine write_talsh_matrix

!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_Fock (exa_input,Matrix_file)
         ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
         ! The new basis set comes from the digonalization of density matrix.
         ! The all natural orbitals and occupation numbers are saved in file NOs_MO 
         ! The transformed Fock matrix are saved in file FM_in_NatOrb
         ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix. 
         ! This routine doesn't truncate the new basis set.    
         ! Written by Xiang Yuan   


                  character(len=*),intent(in)        :: Matrix_file
                  type(exacc_input), intent(in)      :: exa_input

                  type(cmo)                          :: qorbital
                  type(basis_set_info_t)             :: ao_basis
                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  complex(kind=8), pointer           :: MP2_density_comp(:,:)   ! MP2 density matrix from DM_complex
                  real(8), pointer                   :: MP2_density_quat(:,:,:) ! quaternion MP2 density matrix 
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
                  real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
                  integer                            :: sz,flb,sz_occ           ! size of matrix and flag for block
                  integer                            :: i, ierr, nocc, nvirt

                  call read_matrix_comp(MP2_density_comp,Matrix_file)
                  call print_date('read density matrix (DM) in MO basis (spinor) from DM_complex')

                  call complex_to_quaternion(MP2_density_comp,MP2_density_quat)
                  call print_date('changed MO DM representation (spinor->quaternion)')

                  call diagonal_matrix_quater(MP2_density_quat,N_OCCUPATION,NOs_MO)    
                  call print_date('diagonalized MO DM (quaternion)')

                  call sort_values_vectors (N_OCCUPATION,NOs_MO)
                  call print_date('Sorted natural orbitals by occupation number')

                  call quater_write(NOs_MO,N_OCCUPATION,'NOs_MO')
                  call print_date('wrote all natural orbital and occupation number in file NOs_MO')
                    
                  ! we initialize the (diagonal) vv block of the Fock matrix with the orbital energies
                  call read_from_dirac(ao_basis,qorbital,'DFCOEF',ierr)
                  call convert_mo_to_quaternion(qorbital)
                  nvirt = exa_input%nkr_vir        ! the number of virtual orbitals
                  nocc  = exa_input%mokr_occ(size(exa_input%mokr_occ))
                  allocate(Fock_HF(nvirt,nvirt,4))
                  Fock_HF = 0.0
                  forall (i=1:nvirt) Fock_HF(i,i,1) = qorbital%energy(nocc+i)
                  call print_date('read Fock matrix from DFCOEF')

                  call basis_transform_matrix(Fock_HF,NOs_MO,Fock_NO)
                  call print_date('got full VV Fock matrix in natural orbital basis F_VV_NO')

                  call quater_write(Fock_NO,"FM_in_NatOrb")
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

            Sum_N_OCCUPATION_reduced = 0.0 
            do i=1, size(N_OCCUPATION_reduced)
                  Sum_N_OCCUPATION_reduced = Sum_N_OCCUPATION_reduced + N_OCCUPATION_reduced(i)
            end do

            Perc_OCCUPATION_reduced = Sum_N_OCCUPATION_reduced / Sum_N_OCCUPATION_total
            Perc_OCCUPATION_reduced = Perc_OCCUPATION_reduced * 100

            write(*,*) ""
            write(*,"(A50,2X,D10.4)") 'MP2 occupation number (o.n.) threshold (t) : ', threshold
            write(*,"(A50,2X,I10)")    'Number of MP2 FVNOs retained (o.n. > t)    : ', N_reduced
            write(*,*) ""
            write(*,"(A50,2X,D10.4)") 'Sum of MP2 FVNOs o.n., retained NOs (a)    : ', Sum_N_OCCUPATION_reduced
            write(*,"(A50,2X,D10.4)") 'Sum of MP2 FVNOs o.n., all NOs (b)         : ', Sum_N_OCCUPATION_total
            write(*,"(A50,2X,D10.4)") 'Percentage of occupation retained (100*a/b): ', Perc_OCCUPATION_reduced
            write(*,*) ""

            nvirt = exa_input%nkr_vir        ! the number of virtual orbitals
            nocc  = exa_input%mokr_occ(size(exa_input%mokr_occ))
            call read_from_dirac(ao_basis,qorbital,'DFCOEF',ierr)
            call convert_mo_to_quaternion(qorbital)
            allocate(HF_AO(qorbital%nao,nvirt,4))
            HF_AO = qorbital%coeff_q(:,nocc+1:nocc+nvirt,:)
            call print_date('read HF orbital coefficient from DFCOEF')
            
            ! we initialize the (diagonal) vv block of the Fock matrix with the orbital energies
            allocate(Fock_HF(nvirt,nvirt,4))
            Fock_HF = 0.0
            forall (i=1:nvirt) Fock_HF(i,i,1) = qorbital%energy(nocc+i)
            call print_date('read Fock matrix from DFCOEF')

            call basis_transform_matrix(Fock_HF,NOs_MO_reduced,Fock_NO)
            call print_date('transformed Fock matrix into reduced NO basis')
            
            call get_Cano_Orbital(NOs_MO_reduced,New_cano_MO,Fock_NO,New_orbital_energy)
            call print_date('got new Canonical orbital in MO basis')

            call basis_transform_vector(New_cano_MO,HF_AO,NOs_AO)
            call print_date('got new orbitals in AO basis')

            call print_rMP2FVNO(New_orbital_energy,N_OCCUPATION_reduced)
      end subroutine

end module talsh_mp2no

