module talsh_mp2no
!This module contains the routines needed to calculate the mp2 natural orbitals 

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils
        use intermediates
        use talsh_common_routines
        use exacorr_mo
        use exacorr_global,only: nfunctions
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

        public talsh_mp2no_driver

!        interface scale_with_denominators
!           module procedure scale_with_denominators
!        end interface
!

       contains
!-------------------------------------------------------------------------------------------------------------------------------------

       subroutine talsh_mp2no_driver (exa_input,n_cano)
!       This module is used to construct the MP2 density matrix by talsh.
!       The complex density matrix are saved on the file 'DM_complex'.  
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new canonical orbitals are written into new file by Module3.
!       Written by Xiang Yuan
        use exacorr_global
        type(exacc_input), intent(in) :: exa_input
        integer, intent(inout)        :: N_cano                  ! number of new canonical orbital
!        real(kind=8),intent(in)       :: threshold=0.0d-0


!       Point connected with tensors 
        integer                :: ierr
        real(kind=8), pointer  :: N_OCCUPATION_comp(:)    ! natural orbital occupation number in complex diagonalization. For checking 
        real(8),pointer        :: NOs_AO(:,:,:)           ! new quaternion orbital in AO basis
        real(kind=8), pointer  :: N_OCCUPATION_quat(:)    ! occupation number from qdia of DM
        real(kind=8), pointer  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
        logical                :: alive
        real(kind=8)           :: threshold
        integer                :: i

        
        threshold = exa_input%t_mp2no
        print *, "MP2 occupation number Threshold: ", threshold

        inquire(file="NOs_MO", exist=alive)
        if (alive) then
            call print_date('File NOs_MO existed')
        else
            call print_date('File NOs_MO is not existed. Start constructing MP2 density matrix')
            call talsh_generate_MP2_VNO_MObasis(exa_input,N_OCCUPATION_comp)   
            call get_transform_Fock(exa_input,'DM_complex')
!               print *,'Occupation numbers of all natural orbital in complex diago'
!               do i = 1, size(NO_OCC_1,1)
!                     print "(4x,i4,a,2e14.6)",i,"    ",NO_OCC_1(i)
!               end do
            deallocate(N_OCCUPATION_comp)
        end if
        call get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION_quat, &
                              New_orbital_energy,'NOs_MO',N_cano)
        call write_to_newfile(exa_input,NOs_AO,N_OCCUPATION_quat,New_orbital_energy,'DFCOEF','MP2NOs_AO')

        call print_date('Leave MP2_NO Module')
      end subroutine talsh_mp2no_driver


!-----------------------------------------------------------------------------------------------------------------------------------------------

      subroutine talsh_generate_MP2_VNO_MObasis (exa_input,N_OCCUPATION_comp)
!       This subroutine obtains the complex mp2 natural orbital in MO basis
!       Do complex diag to obtain complex occupation number. However, we don't use these orbital and occupation number. 
!       The non-diagonal complex DM is written into file DM_complex.
!       Written by Xiang Yuan


         use exacorr_global
!         use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
         use talsh_gradient


         type(exacc_input), intent(in)          :: exa_input
         real(kind=8), pointer,intent(inout)    :: N_OCCUPATION_comp(:)    ! natural orbital occupation number in complex diagonalization. For checking
         integer                                :: nocc, nvir              ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable                  :: mo_vir(:)               ! the list index of virtual spinors
         integer,  allocatable                  :: mo_occ(:)               ! the list index of occupied spinors
         real(8), allocatable                   :: eps_occ(:),eps_vir(:)   ! HF orbital energies
         complex(kind=8), pointer               :: NOs_MO(:,:)             ! complex natural orbital

!        Intermediate tensors (same notation as in relccsd)
         type(talsh_dm_tens)         :: MP2_density_tensor
         integer(INTD), dimension(4) :: vvoo_dims
         integer(INTD), dimension(2) :: vv_dims
         integer(C_INT)              :: one_dims(1)

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens)       :: int_t

!        Point connected with tensors 
         integer(C_INT)              :: ierr
         complex(8), pointer         :: MP2_density_tens(:,:)
         integer                     :: initialize_talsh

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
         ierr=talsh_tensor_construct(MP2_density_tensor%prime_vv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(MP2_density_tensor%prime_vv,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,MP2_density_tens,vv_dims)
         call print_exainput(exa_input,get_nao())
         call print_date("Getting orbital energies")


!        Get HF orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("Retrieved orbital energies")

!        Get oovv integrals 
         call get_MP2_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%print_level)

!        Make some noise so that we know we are proceeding 
         call print_date('Forming VV block of MP2 density matrix')
!        Initialize scale with denominators of oovv integrals
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=int_t%oovv)

!        Calculate and print the MP2 density matrix of virtual-virtual block Dvv 
         ierr=talsh_tensor_init(MP2_density_tensor%prime_vv)
         ierr=talsh_tensor_contract("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_density_tensor%prime_vv,&
                                        int_t%oovv,int_t%oovv,scale=ONE_HALF) ! o(2)v(2)

!        Make some noise so that we know we are proceeding 
         call print_date('Save complex non-digonal DM to file DM_complex')
         call write_matrix(exa_input,MP2_density_tensor,'DM_complex')



!        Diagoliazation of Density matrix
!         allocate(NO_OCC(2,size(MP2_density_tens,2)))
         allocate(N_OCCUPATION_comp(size(MP2_density_tens,2)))
!         call diagoanl_matrix_complex(MP2_density_tens,NO_OCC,NOs_MO) 
         N_OCCUPATION_comp = 0
!        N_OCCUPATION_comp(:) = NO_OCC(1,:)


        end subroutine talsh_generate_MP2_VNO_MObasis

!----------------------------------------------------------------------------------------------------------------------------------------------

    subroutine get_MP2_integrals (nocc,nvir,mo_occ,mo_vir,int_t,print_level)
                
!    Routine to get oovv integral in the form of antisymmetric tensors.
!    In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!    full integral list.
!    copy from talsh_cc.F90 routine: get_CC_integrals
         
      use exacorr_global
      use talsh_ao_to_mo
         
      integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
      integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!     target integral tensors to be filled and given back to the caller
      type(talsh_intg_tens),intent(inout)  :: int_t
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
        
!---------------------------------------------------------------------------------------------------------------------------------------------------

      subroutine write_matrix(exa_input,gamma,file_out)
! Routine for writing complex density matrix into the file'DM_complex'. 

         use exacorr_global, only : get_nsolutions
         use exacorr_utils

! input variables
         type(exacc_input),   intent(in   )  :: exa_input
         type(talsh_dm_tens), intent(inout)  :: gamma
         character(len=*),intent(in)         :: file_out

! matrix dimensions
         integer :: nocc, nvir
! body pointers
         type(C_PTR)                     :: body_p2
         complex(8), pointer, contiguous :: dm_block(:,:)
         complex(8), allocatable         :: dm_mo_fort(:,:)
! error code
         integer :: ierr
! loop variables
         integer :: i, j, k
         integer :: mp2dm
! debug
         logical :: debug=.false.

         
      !***************
      !     set dims *
      !***************
            nocc    = 0
            nvir    = exa_input%nvir
         
      !***************************************
      !     copy DM blocks into one large DM *
      !***************************************
         
            allocate (dm_mo_fort(nvir,nvir))
            dm_mo_fort = ZERO
         
      !-----------
      !| OO | OV |
      !-----------
      !| VO | VV |
      !-----------
                  
      !vv
            ierr=talsh_tensor_get_body_access(gamma%prime_vv,body_p2,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p2,dm_block,(/nvir,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
            dm_mo_fort(1:nvir,1:nvir) = dm_block
         
            call get_free_fileunit(mp2dm)
            open (mp2dm,file=file_out,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
      ! write dimension of the matrix 
            write(mp2dm) nocc+nvir
      ! write matrix element 
            write(mp2dm) dm_mo_fort

            if (debug) then
              print *,'<<<<< DM, complete, before formatting to DM_complex'
                 do j=1,(nocc+nvir)
                    do i=1,(nocc+nvir)
                       write(mp2dm,'(2x,2F28.16,2x,2I4)') dm_mo_fort(i,j),i,j
                    end do
                  end do
                  print *,'>>>>> DM, complete, after formatting to DM_complex'
            end if

            close(mp2dm,status='KEEP')

        end subroutine write_matrix


!------------------------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_Fock (exa_input,Matrix_file)
         ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
         ! The new basis set comes from the digonalization of density matrix.
         ! The all natural orbitals and occupation numbers are saved in file NOs_MO 
         ! The transformed Fock matrix are saved in file FM_in_NatOrb
         ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix. 
         ! This routine doesn't truncate the new basis set.    
         ! Written by Xiang Yuan   

            use matrix_complex_to_quaternion 
            use write_orbital_to_file

                  character(len=*),intent(in)        :: Matrix_file
                  type(exacc_input), intent(in)      :: exa_input

                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  complex(kind=8), pointer           :: MP2_density_comp(:,:)   ! MP2 density matrix from DM_complex
                  real(8), pointer                   :: MP2_density_quat(:,:,:) ! quaternion MP2 density matrix 
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
                  real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
                  integer                            :: sz,flb                  ! size of matrix and flag for block
                  character                          :: mat                     ! flag used for matrix basis transform
                  integer                            :: i

                                                        
                  call print_date('Read Density Matrix from DM_complex')
                  call read_matrix_comp(MP2_density_comp,Matrix_file)

                  call print_date('Transform DM from complex into quaternion')
                  call complex_to_quaternion(exa_input,MP2_density_comp,MP2_density_quat)

                  call print_date('Qdiag Density matrix')
                  call diagoanl_matrix_quater(MP2_density_quat,N_OCCUPATION,NOs_MO)           

                  call print_date('Sort quaterbion natural orbital')
                  call sort_values_vectors (N_OCCUPATION,NOs_MO)

                  call print_date('Write the all natural orbital and occupation number in file NOs_MO')
                  call quater_write(NOs_MO,N_OCCUPATION,'NOs_MO')
                    
                  call print_date('Read Fock matrix from DFCOEF')
                  flb = 0                     ! we read vv block of Fock matrix
                  sz = exa_input%nvir/2       ! the size of the FM is the same as DM
                  call quater_read_from_dirac_fock(Fock_HF,sz,flb)

                  mat = 'Matrix'
                  call basis_transform(Fock_HF,NOs_MO,Fock_NO,mat)
                  call print_date('Get full vv Fock matrix in natural orbital basis')
                  call print_date('Write the Full Fock matrix of natural orbital in file FM_in_NatOrb')
                  call quater_write(Fock_NO,"FM_in_NatOrb")
                        
      end subroutine get_transform_Fock

         
!------------------------------------------------------------------------------------------------------------------------------
      subroutine get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION_reduced,New_orbital_energy,file_orb,N_reduced)
      ! This routine is used for getting new truncated canonical orbital in AO basis
      ! and the corresponding orbital energy.

            use matrix_complex_to_quaternion 
            use write_orbital_to_file


            type(exacc_input), intent(in)      :: exa_input
            real(kind=8),intent(in)            :: threshold
            real(kind=8), pointer,intent(out)  :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis 
            real(kind=8), pointer,intent(out)  :: N_OCCUPATION_reduced(:) ! occupation number from qdia of DM
            real(kind=8), pointer,intent(out)  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
            character(len=*),intent(in)        :: file_orb
            integer,intent(inout)              :: N_reduced               ! number of reduced natural orbital 


            real(kind=8), pointer              :: NOs_MO_reduced(:,:,:)   ! natural orbitals which have occupation number larger than threshold
            real(kind=8), pointer              :: New_cano_orb(:,:,:)     ! new canonical orbital 
            real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
            real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
            real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital 
            integer                            :: sz,flb, fl              ! size of matrix and flag for block
            character                          :: mat                     ! flag used for matrix basis transform
            integer                            :: i,j,k
            

            call print_date('Read and truncate the natural orbital and occupation number from NOs_MO')
            call quater_read(NOs_MO_reduced,N_OCCUPATION_reduced,threshold,N_reduced,file_orb)


            call print_date('Read HF orbital coefficient from DFCOEF')
            flb = 0                       ! we only read virtual orbital 
            sz = exa_input%nvir/2         ! the number of virtual orbital 
            call quater_read_from_dirac_orbital(HF_AO,sz,flb)


            call print_date('Read Fock matrix from DFCOEF')
            flb = 0                         ! we read vv block of Fock matrix
            sz = exa_input%nvir/2           ! the size of the FM is the same as DM
            call quater_read_from_dirac_fock(Fock_HF,sz,flb)

            mat = 'Matrix'
            call basis_transform(Fock_HF,NOs_MO_reduced,Fock_NO,mat)
            call print_date('Get reduced vv Fock Matrix in Natural Orbital Basis')

            call print_date('Get new Canonical orbital in MO basis')
            call get_Cano_Orbital(NOs_MO_reduced,New_cano_orb,Fock_NO,New_orbital_energy)


            print *, "Print the first 10 New orbital Energy"
            do i=1, 10
               print *, 'Orbital',i,':Energy: ', New_orbital_energy(i)
            end do

!      Debug part, write HF orbital into the new file to test the code.
!*********************************************** 
!            New_cano_orb = 0.0                
!            do i=1, size(New_cano_orb,2)      
!                  New_cano_orb(i,i,1) = 1.0   
!            end do                            
!***********************************************

            call print_date('Get new orbital in AO basis')
            fl = 0 
            call basis_transform(New_cano_orb,HF_AO,NOs_AO,fl)

      end subroutine


!-----------------------------------------------------------------------------------------------------------------------------
      subroutine write_to_newfile (exa_input,New_cano_orb,N_OCCUPATION,New_orbital_energy,file_old,file_new)
  !      write new orbital information into new file 'MP2NOs_AO', which is similar to DFCOEF 
  !      copy the basis information from DFCOEF and use new orbital energy and coefficient to replace the old orbital.   
  !      revise subroutine write_mos_to_dirac of exacorr_dirac.F

         use exacorr_mo
         use exacorr_datatypes
         use exacorr_utils
         use write_orbital_to_file
  
         implicit none
    
         type(exacc_input), intent(in)         :: exa_input
         character(len=*),intent(in)           :: file_old,file_new
         real(8),pointer,intent(in)            :: New_cano_orb(:,:,:)    ! new canonical orbital in AO basis
         real(kind=8),pointer,intent(in)       :: N_OCCUPATION(:)        ! occupation number
         real(kind=8),pointer,intent(in)       :: New_orbital_energy(:)  ! new orbital energy 


         type(basis_set_info_t)   :: ao_basis
         type(qumo)               :: qorbital        ! save the new canonical orbital 
         type(qumo)               :: qorbital_hf     ! save the original HF canonical orbital 

         integer                  :: nocc            ! occupied orbital 
         integer                  :: nao, nmo, nno   ! number of atomic orbital; molecular orbital; natural orbital
         integer,allocatable      :: index(:)
         integer                  :: i,j
         integer                  :: ierr
         real(kind=8)             :: diff_1, diff_2, diff_3, diff_4


!        read information of from DFCOEF
         call quater_read_from_dirac_all(ao_basis,qorbital_hf,file_old,ierr)

         nocc = (exa_input%nocc)/2
         nao = qorbital_hf%nao
         nno = size(New_cano_orb,2)
         nmo = nocc+nno
         allocate(index(nmo))
         do i = 1, nmo
            index(i) = i ! we set the index compatible with the fact that we store all mo's
         end do

         call alloc_mo(qorbital,nao,nmo,4)
         qorbital%allocated = qorbital_hf%allocated
         qorbital%restricted = qorbital_hf%restricted
         qorbital%complex_algebra = qorbital_hf%complex_algebra
         qorbital%total_energy = qorbital_hf%total_energy
         qorbital%fermion_irrep = qorbital_hf%fermion_irrep
         qorbital%boson_irrep = 0
         qorbital%index = index


         do i = 1, nocc
            qorbital%energy(i)=qorbital_hf%energy(i)
         end do

         do i = nocc+1, qorbital%nmo
            qorbital%energy(i)=New_orbital_energy(i-nocc)
!            qorbital%energy(i)=qorbital_hf%energy(i)
         end do

         do j = 1, nocc
            do  i = 1, qorbital%nao
                  qorbital%coeff(i,j,1)=qorbital_hf%coeff(i,j,1)
                  qorbital%coeff(i,j,2)=qorbital_hf%coeff(i,j,2)
                  qorbital%coeff(i,j,3)=qorbital_hf%coeff(i,j,3)
                  qorbital%coeff(i,j,4)=qorbital_hf%coeff(i,j,4)
            enddo
         enddo


         do  j = nocc+1, qorbital%nmo
            do  i = 1, qorbital%nao
                  qorbital%coeff(i,j,1)=New_cano_orb(i,j-nocc,1)
                  qorbital%coeff(i,j,2)=New_cano_orb(i,j-nocc,2)
                  qorbital%coeff(i,j,3)=New_cano_orb(i,j-nocc,3)
                  qorbital%coeff(i,j,4)=New_cano_orb(i,j-nocc,4)
            enddo
         enddo


         call quater_write_mos_to_dirac(ao_basis,qorbital,file_new)

        end subroutine write_to_newfile


!---------------------------------------------------------------------------------------------------------------------------------------
      subroutine diagoanl_matrix_complex(M,eVal,eVectorsR)
      ! This routine is used to do complex digonalization in debug. 
      ! Set parameter needed in lapack

            integer :: lwork, info,i,j
            real(kind=8), allocatable :: work(:)
            real(kind=8), allocatable :: rwork(:)
            integer :: vectors_L, vectors_R
            complex(kind=8), intent(inout), pointer :: M(:,:)
            real(kind=8), pointer :: eVal(:,:) 
            complex(kind=8), pointer :: eVectorsR(:,:)
            complex(kind=8), pointer :: eVectorsL(:,:)
            character(len=1) :: do_left = 'N', do_right = 'V'      

            vectors_L = size(M,1)
            vectors_R = size(M,2)
            lwork = 1000000*vectors_R	
            allocate(eVectorsL(vectors_L,vectors_R))
            allocate(eVectorsR(vectors_L,vectors_R))
            allocate(eVal(2,vectors_R))
            allocate(work(MAX(1,LWORK)))
            allocate(rwork(2*vectors_R))

            call zgeev (do_left, do_right, vectors_R, M, vectors_L, &
                     eVal, eVectorsL, vectors_L, eVectorsR, vectors_R, &
                     WORK, LWORK, rwork, INFO)

            deallocate(eVectorsL)
            deallocate(work)
            deallocate(rwork)

            call print_date('After diagonal_matrix_complex')
         end subroutine diagoanl_matrix_complex
!------------------------------------------------------------------------------------------------------------

end module talsh_mp2no

