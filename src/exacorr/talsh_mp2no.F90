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
!-------------------------------------------------------------------------------------------
        subroutine talsh_mp2no_driver (exa_input)

!       This subroutine drives the mp2 natural orbital calculation
!       Written by Xiang Yuan, July 2020

         use exacorr_global
         use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
         use talsh_gradient
         use davidson

         type(exacc_input), intent(in) :: exa_input

         integer                :: nocc, nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         real(8), allocatable   :: eps_occ(:),eps_vir(:) ! HF orbital energies

!        Point connected with tensors 
         integer(C_INT)       :: ierr
         complex(8), pointer :: MP2_density_tens(:,:)
         complex(kind=8), pointer :: NOs_MO(:,:)
         complex(kind=8), pointer :: test_tens(:,:)
         real(kind=8), pointer :: eVectors_Real(:,:)  ! natural orbital in real number used for sorted routine 
         real(kind=8), pointer :: NO_OCC(:,:)         ! natural orbital occupation number 
         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max

!        Common tensors
         type(talsh_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens) :: int_t
         
!        Intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t) :: result_tensor
         type(talsh_dm_tens) :: MP2_density_tensor
         type(talsh_tens_t) :: NOs_MO_tensor  ! natural orbital in MO basis
         integer(INTD), dimension(4) :: oooo_dims, oovv_dims, vvoo_dims, vovo_dims
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vv_dims
         integer(C_INT)              :: one_dims(1)
         real(kind=8)                :: energy_shift(2)
         type(basis_set_info_t)      :: ao_basis ! information about the ao basis
         type(cmo)                   :: cspinor_all ! Stores all MO-coefficients and energies obtained from the parent code

!        Scalars (need to be defined as tensor types)
         real(kind=8)       :: threshold
         integer            :: i,j
         integer            :: num_no_red       ! number of reduced natural orbital
         integer            :: sorting_order_local
         logical            :: left_local, right_local, apply_evalues_shift  !variable in davidson module controling sorted routine 

         


!        Make some noise so that we know we are here
         call print_date('Entered mp2_driver routine')

!        Initialize libraries
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)
         call print_date('Initialized talsh library')
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max
         call interest_initialize(.true.)
         call print_date('Initialized interest library')

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         print *,'number of virtual orbital', nvir
         print *,'number of occupied orbital', nocc
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir
         vv_dims(1) = nvir
         vv_dims(2) = nvir
         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(MP2_density_tensor%prime_vv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(MP2_density_tensor%prime_vv,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,MP2_density_tens,vv_dims)
         print *,'construct mp2_tensor',ierr
         print *,'mo_occ', mo_occ
         print *,'mo_vir', mo_vir  
         call print_exainput(exa_input,get_nao())
         call print_date("Getting orbital energies")


!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("Retrieved orbital energies")

!        Get oovv integrals 
         call get_MP2_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%print_level)

!        Get fixed one-body tensors, needs to follow the two body tensors
!        call get_MP2_fockmatrix (nocc,nvir,mo_occ,mo_vir,diag_foo_tensor, diag_fvv_tensor,int_t, exa_input%print_level)


!        Initialize t2 amplitudes with the MP2 values
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc

         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=int_t%oovv)


!        Calculate and print the MP2 density matrix within virtual molecular basis
         ierr=talsh_tensor_init(MP2_density_tensor%prime_vv)
         ierr=talsh_tensor_contract("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_density_tensor%prime_vv,&
                                        int_t%oovv,int_t%oovv,scale=ONE_HALF) ! o(2)v(2)
         print *,'after constructing density matrix',ierr
!         print *,'Row_DM', size(MP2_density_tens,1)
!         print *,'Col_DM', size(MP2_density_tens,2)


!        Make some noise so that we know we are leaving
         write(*,*) ""
         call print_date('Go to natural orbital')


!        Diagoliazation of Density matrix
         allocate(NO_OCC(2,size(MP2_density_tens,2)))
         NO_OCC = 0.0
         call diagoanl_DMP2_matrix(MP2_density_tens,NO_OCC,NOs_MO)
         !do i=1,size(NO_OCC,2),1
         !       print *,'occupation number',NO_OCC(1,i)
         !end do


!        Sorted occupation number and natural orbital 
         allocate(eVectors_Real(size(NOs_MO,1),size(NOs_MO,2))) 
         energy_shift =0.0
         sorting_order_local = 1
         apply_evalues_shift = .false.
         left_local = .false. 
         right_local = .true.
         do i=1,size(NOs_MO,1),1
                eVectors_Real (:,i)= real(NOs_MO(:,i))
         end do
         call sort_values_vectors(sorting_order_local, NO_OCC, eVectors_Real, eVectors_Real, &
                                     right_local, left_local,apply_evalues_shift,energy_shift)

         do i=1,size(NOs_MO,1),1
                do j=1,size(NOs_MO,2),1            
                     NOs_MO(j,i)=complex(eVectors_Real(j,i),0.0)  
                end do 
         end do

         call MP2_read_write_from_dirac (ao_basis,cspinor_all,'DFCOEF',ierr)
!---------------------------------------------------------------------------------------
        ! threshold = 1.0d-5
        if (threshold /= 0) then
            call get_MP2_FNOs (threshold,nocc,nvir,mo_vir,eVectors_Real,NO_OCC)
        end if 

!--------------------------------------------Test part-----------------------------------------------------
   
          do i=1,size(NO_OCC,2),1
             print *,'occupation number_order',NO_OCC(1,i)
          end do
          
          do j=1,size(eVectors_Real,2),1
             print *,'natural orbital in MO basis_order_in complex',NOs_MO(:,j)
          end do
!         print *, '!--------------------------------------------------------------'
!         ierr=talsh_tensor_get_body_access(mocoefb_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
!         call c_f_pointer(body_p,test_tens,mocoef_dims) ! to use as a regular Fortran array
!         do i=1,nmo,1 
!                print *,'MO_AO coeff matrix_b ',test_tens(:,i)
                !print *, 'size row', size(test_tens,1)
                !print *, 'size col', size(test_tens,2)
                !print *,'MO_AO coeff matrix ',test_tens
!         end do
!-------------------------------------------------------------------------------------------------------------
        end subroutine talsh_mp2no_driver



!----------------------------------------------------------------------------------------------------------------------    
        subroutine diagoanl_DMP2_matrix(M,eVal,eVectorsR)

                ! set parameter needed in lapack
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
                    print *,'vector_R', vectors_R
                    print *,'vector_L', vectors_L
                    lwork = 10000000*vectors_R	
                    allocate(eVectorsL(vectors_L,vectors_R))
                    allocate(eVectorsR(vectors_L,vectors_R))
                    allocate(eVal(2,vectors_R))
                    allocate(work(MAX(1,LWORK)))
                    allocate(rwork(2*vectors_R))
                    call zgeev (do_left, do_right, vectors_R, M, vectors_L, &
                    eVal, eVectorsL, vectors_L, eVectorsR, vectors_R, &
                    WORK, LWORK, rwork, INFO)
                    !print *,'INFO',info
                      
        end subroutine diagoanl_DMP2_matrix
            
            
!----------------------------------------------------------------------------------------------------------------

        subroutine get_MP2_mocoef_tensor (nao, nmo, mo_list, mocoefa_tensor, mocoefb_tensor)

                !This subroutine reads the desired subset of mo coefficients from DIRAC
                
                !        Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 30 C)
                
                         use exacorr_mo
                         use exacorr_global
                
                         implicit none
                         integer, intent(in) :: nao          ! the length of the ao basis
                         integer, intent(in) :: nmo          ! the length of the mo basis
                         integer, intent(in) :: mo_list(:)   ! and their indices
                         type(talsh_tens_t), intent(inout) :: mocoefa_tensor,mocoefb_tensor
                
                         integer(INTD)       :: ierr, tens_rank
                         type(cmo)           :: cspinor
                         integer(INTD)       :: mocoef_dims(1:2)
                
                         complex(8), pointer, contiguous:: camo_tens(:,:),cbmo_tens(:,:)
                         type(C_PTR):: body_p
                
                         ierr=talsh_tensor_dimensions(mocoefa_tensor,tens_rank,mocoef_dims)
                         ierr=talsh_tensor_get_body_access(mocoefa_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
                         call c_f_pointer(body_p,camo_tens,mocoef_dims(1:2)) ! to use camo_tens as a regular Fortran 2d array
                         ierr=talsh_tensor_get_body_access(mocoefb_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
                         call c_f_pointer(body_p,cbmo_tens,mocoef_dims(1:2)) ! to use cbmo_tens as a regular Fortran 2d array
                
                         call get_mo_coefficients (cspinor,mo_list,nmo)
                
                !        Copy the coefficients into the alpha and beta tensors, using the quaternion relations as defined in
                !        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
                !        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
                !        using Kramers' symmetry.
                
                         camo_tens = cspinor%coeff(:,:,1)
                         cbmo_tens = cspinor%coeff(:,:,2)
                
        end subroutine get_MP2_mocoef_tensor



!--------------------------------------------------------------------------------------------------------------------

        subroutine get_MP2_integrals (nocc,nvir,mo_occ,mo_vir,int_t,print_level)

  !        Routine to get integrals needed in CC in the form of antisymmetric tensors.
  !        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
  !        full integral list.
  
  !        Written by Lucas Visscher, summer 2017
  !        Revised by Xiang Yuan, July 2020
  
           use exacorr_global
           use talsh_ao_to_mo
  
           integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
           integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
  !        target integral tensors to be filled and given back to the caller
           type(talsh_intg_tens),intent(inout)  :: int_t
           integer, intent(in)                  :: print_level
  
           integer(INTD) :: oovv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4)
  
           integer(INTD)                      :: ierr
           complex(8), pointer, contiguous    :: ao_tens(:,:,:,:)
           type(C_PTR)                        :: body_p
  
           integer :: nao
  
  !        auxilliary integral tensors (only needed inside this routine)
           type(talsh_tens_t) :: aoint_tensor, ovvo_tensor
           integer(INTD)      :: aoint_dims(1:4), ovvo_dims(1:4)
  
  !        arrays needed to communicate with the MO integral generator
           integer:: nmo(4)
           integer, allocatable :: mo_list(:)
  
  !        scalars (need to be defined as tensor types)
           type(talsh_tens_t) :: minusone_tensor
           integer(INTD)      :: minusone_dims(1)
  
           call print_date('entering integral transformation, scheme = 2')
  
  !        Retrieve basis set information
           nao     = get_nao()      ! number of basis functions
  
  !        Get tensor with AO integrals
           aoint_dims      = nao
           ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
           if (ierr /= 0) then
              print*, "Could not construct ao integral tensor, increasing TALSH_BUFF may help"
              call quit ('error constructing aoint_tensor')
           end if
           call print_date('allocated AO integral tensor')
  !        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
           ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
           call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
           call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
           call print_date('finished ao integral calculation')
  
  !        Get antisymmetric oovv tensor
  !        Note that the dimensions for the tensor are given in physicists notation <12||34>
           oovv_dims(1:2) = nocc
           oovv_dims(3:4) = nvir
           ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
  !        whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
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




!------------------------------------------------------------------------------------------------------------

     subroutine get_MP2_FNOs (threshold,nocc,nvir,mo_vir,eVectors_Real,NO_OCC)
!        Routine to et theshold and reduce natural orbital space including transform from MO basis to AO basis 

        integer, intent(in)                :: nocc,nvir      ! the size of the mo basis for occupied and virtual spinors
        integer, intent(in)                :: mo_vir(:)      ! the list of virtual orbitals
        real(kind=8), pointer, intent(in)  :: eVectors_Real(:,:)  ! natural orbital in real number used for sorted routine 
        real(kind=8), intent(in)           :: threshold
        real(kind=8), pointer,intent(in)   :: NO_OCC(:,:)         ! natural orbital occupation number 

        integer(C_INT)                     :: ierr
        type(C_PTR)                        :: body_p
        integer(C_SIZE_T)                  :: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT)                     :: host_arg_max
        type(talsh_tens_t)                 :: mocoefa_tensor, mocoefb_tensor ! HF molecular orbital in AO basis 
        type(talsh_tens_t)                 :: NOs_MO_red_tensor  ! reduced natural orbital in MO basis
        type(talsh_tens_t)                 :: NOs_AO_tensor  ! reduced natural orbital in AO basis
        complex(kind=8), pointer           :: NOs_MO_red(:,:)
        integer, allocatable               :: mp2no_mo_list(:)   ! molecular orbital indices
        integer                            :: num_no_red
        integer(INTD), dimension(2)        :: no_red_dims, no_ao_dims
        integer(INTD)                      :: mocoef_dims(2)
        integer                            :: nao, nmo     ! the size of the ao & mo basis for natural orbital transformation 
        integer                            :: i,j
        complex(kind=8), pointer           :: test_tens(:,:)

        num_no_red=0
        do i=1,size(eVectors_Real,1),1
               if (NO_OCC(1,i) > threshold) then
                       num_no_red = num_no_red +1
               end if 
        end do
        print *,'number of NOs kept in following',num_no_red
        no_red_dims(1)=nvir
        no_red_dims(2)=num_no_red
        ierr=talsh_tensor_construct(NOs_MO_red_tensor,C8,no_red_dims,init_val=ZERO)
        ierr=talsh_tensor_get_body_access(NOs_MO_red_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,NOs_MO_red,no_red_dims) ! to use as a regular Fortran array
        do i=1,size(NOs_MO_red,2),1
               do j=1,size(NOs_MO_red,1),1            
                        NOs_MO_red(j,i)=complex(eVectors_Real(j,i+nvir-num_no_red),0.0)  
               end do 
        end do
      
        do j=1,size(NOs_MO_red,2),1
               print *,'reduced natural orbital in MO basis',NOs_MO_red(:,j)
        end do


!        Extracting virtual MO coeffcient in AO basis set
        nao = (nocc + nvir)/2
        nmo = nvir 
        mocoef_dims(1) = nao
        mocoef_dims(2) = nmo
        allocate(mp2no_mo_list(nvir))
        mp2no_mo_list = mo_vir
        ierr=talsh_tensor_construct(mocoefa_tensor,C8,mocoef_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(mocoefb_tensor,C8,mocoef_dims,init_val=ZERO)
        call get_MP2_mocoef_tensor (nao, nmo, mp2no_mo_list, mocoefa_tensor, mocoefb_tensor) 
        ierr=talsh_tensor_add("moco_a(a,b)+=moco_b(a,b)",mocoefa_tensor,mocoefb_tensor)


!        Transform reduced natural orbital from MO basis into AO basis
        no_ao_dims(1)=nao
        no_ao_dims(2)=num_no_red
        ierr=talsh_tensor_construct(NOs_AO_tensor,C8,no_ao_dims,init_val=ZERO)
        ierr=talsh_tensor_get_body_access(NOs_AO_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        ierr=talsh_tensor_contract("NO_AO(a,b)+=MO_AO(a,i)*NO_MO(i,b)",NOs_AO_tensor,mocoefa_tensor,NOs_MO_red_tensor,scale=ONE)

        
!--------------------------------------------Test part--------------------------------------------------------------------------
        print *, '!--------------------------------------------------------------'
        !ierr=talsh_tensor_get_body_access(NOs_AO_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,test_tens,no_ao_dims) ! to use as a regular Fortran array
        print *, 'size test_tens', size(test_tens,2)
        do i=1,num_no_red,1 
                print *,'reduced natural orbital in AO basis',test_tens(:,i)
        end do

        print *, '!--------------------------------------------------------------'
        ierr=talsh_tensor_get_body_access(mocoefa_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,test_tens,mocoef_dims) ! to use as a regular Fortran array
        print *, 'size test_tens', size(test_tens,2)
        do i=1,nmo,1 
                !print *,'MO_AO coeff matrix_a ',test_tens(:,i)
                !print *, 'size row', size(test_tens,1)
                !print *, 'size col', size(test_tens,2)
                print *,'MO in AO basis coeff matrix ',test_tens(:,i)
        end do
!--------------------------------------------Test part--------------------------------------------------------------------------
      end subroutine get_MP2_FNOs

!----------------------------------------------------------------------------------------------------------------------------
      subroutine MP2_read_write_from_dirac (ao_basis,cspinor,file_in_name,ierr)

        use exacorr_mo
        use exacorr_datatypes
        use exacorr_global,only: nfunctions
        use exacorr_utils
 
        implicit none

 !      read all positive energy mo coefficients present on DFCOEF file
 !      check for imcompatibiklities (use of symmetry) as well and continue
 !      by reading also the basis set information into gtos
 
        character(len=*),intent(in) :: file_in_name
        type(basis_set_info_t), intent(inout) :: ao_basis
        type(cmo), intent(inout):: cspinor
        integer, intent(out)   :: ierr 
        integer, parameter     :: NOERROR  = 0 ! no error
        integer, parameter     :: NOTFOUND = 1 ! file not found
        integer, parameter     :: CORRUPT  = 2 ! file unreadable
        integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
        integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old DFCOEF?)
        integer, parameter     :: CONTRAC  = 5 ! basis set is contracted (Interest can not handle this)
        integer, parameter     :: INCONSIS = 6 ! DFCOEF does not contain consistent information wrt number of basis functions
        character              :: text*74      ! holds text read from DFCOEF, only used for printing this
        integer                :: nsym
        integer                :: idim(3,2)
        real(8)                :: total_energy
        logical                :: tobe
        integer                :: lucoef
        integer                :: inf_mp2no
        integer                :: i, j, k, ishell, iprimstart, iprimend
        integer                :: nao, nao_bas, nshells, iao, ibas
        integer                :: nmo, npo, nto, npriexp, ncnt
        integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
        real(8), allocatable   :: coeff(:,:,:),eig(:), seig(:)
        integer, allocatable   :: ibeig(:),index(:)
        integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
        integer, allocatable   :: jstrt(:)
        real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)
        complex(8), allocatable:: scoeff(:,:,:)
        integer, allocatable   :: ao_map(:,:)
        integer                :: bi
 
        ierr = NOERROR
 
        ! open file and read dimension info
        inquire (file=file_in_name,exist=tobe)
        bi=0
        if (.not.tobe) then
           print*, "error: MO coefficient file"//file_in_name//" not found"
           ierr = NOTFOUND
           return
        else
           !call get_free_fileunit(lucoef)
           !call get_free_fileunit(inf_mp2no)
           lucoef=100
           inf_mp2no = lucoef+10
           open (lucoef,file=file_in_name,status='OLD',FORM='UNFORMATTED',&     
                access='SEQUENTIAL')
           !open (inf_mp2no,file='MP2NOs',FORM='UNFORMATTED',&     
                !access='SEQUENTIAL')
           open (inf_mp2no,file='MP2NOs',FORM='FORMATTED',&     
                access='SEQUENTIAL')
           !open (unit=99,file='test.txt')
           read (lucoef,end=10,err=10) text,nsym,&         
                ((idim(i,j),i=1,3),j=1,nsym),total_energy
        end if
 
        ! check symmetry and get dimensions of the coefficient array
        if (nsym > 1) then
           ierr = SYMMETRY
           close (lucoef,status='KEEP')
           close (inf_mp2no,status='KEEP')
           return
        else
           nao = idim(3,1) ! number of basis functions
           npo = idim(1,1) ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
           nmo = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
           nto = npo+nmo   ! total number of solutions on this file
           if (npo > 0) print*, "warning: reading from 4C DFCOEF file" ! this warning will go once we have the full index transform implemented
        end if
 
        allocate(coeff(nao,nto,nz))
        allocate(eig(nto))
        allocate(ibeig(nto))
        allocate(index(2*nmo))
         
        ! read the mo coefficients and orbital eigenvalues
        read (lucoef,end=10,err=10) coeff
        read (lucoef,end=10,err=10) eig
        read (lucoef,end=10,err=10) ibeig
        ! fill in index values
        do i = 1, 2*nmo
           index(i) = i ! we set the index compatible with the fact that we store all mo's
        end do
         !WRITE(inf_mp2no) TEXT,NSYM,((IDIM(I,J),I = 1,3),J=1,NSYM),total_energy
         !WRITE(inf_mp2no) coeff
         !WRITE(inf_mp2no) EIG
         !WRITE(inf_mp2no) IBEIG
         WRITE(inf_mp2no,'(A74)') TEXT
         WRITE(inf_mp2no,'(7(X,I0))') NSYM,((IDIM(I,J),I = 1,3),J=1,NSYM)
         WRITE(inf_mp2no,'(E24.16)') total_energy
         WRITE(inf_mp2no,'(6F22.16)') coeff
         WRITE(inf_mp2no,'(6E22.12)') EIG
         WRITE(inf_mp2no,'(66(X,I0))') IBEIG
 !      switch to unrestricted spinor picture
        allocate(seig(2*nmo))
        allocate(scoeff(nao,2*nmo,2)) ! nao, 2*nmo, spin
        call convert_quaternion_to_complex(coeff,scoeff,eig,seig,&                                              
             nao,nmo,npo)
 
        ! put in cmo
        call alloc_mo (cspinor,nao,2*nmo)
        cspinor%total_energy = total_energy
        cspinor%index        = index
        cspinor%coeff        = scoeff
        cspinor%energy       = seig
        cspinor%boson_irrep  = 0 ! because nosym. if there is, include ibeig in converter
 
        deallocate(coeff)
        deallocate(eig)
        deallocate(ibeig)
        deallocate(index)
        deallocate(seig)
        deallocate(scoeff)
 
        ! read the basis set information and check for consistency
        read (LUCOEF,err=20,end=20) nao_bas,nshells,npriexp,ncnt 
        if (nao < nao_bas) then
           ierr = INCONSIS
           close (lucoef,status='KEEP')
           return
        elseif (nao > nao_bas) then
           print*,"WARNING: assuming Levy-Leblond Hamiltonian"
           print*,"Reading fewer basis functions than MO-coefficients"
           allocate (ao_map(2,nao_bas))
           forall (i=1:nao_bas) ao_map(:,i)=i
           call compress_mo (cspinor,ao_map)
        end if
 
        ! allocate the arrays to read the basis set information in DIRAC format
        allocate(ncent(nshells))
        allocate(nhkt(nshells))
        allocate(nuco(nshells))
        allocate(numcf(nshells))
        allocate(jstrt(nshells))
        allocate(cent(nshells,3))
        allocate(priexp(npriexp))
        allocate(priccf(npriexp,ncnt))
 
 !      orb_momentum, atom_number, n_primives for all shells
        read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
 !      start of exponent values, number of contracted for all shells
        read (LUCOEF,err=20,end=20) jstrt,numcf
 !      coordinates (x, y, z) of center for all shells
        read (LUCOEF,err=20,end=20) cent
 !      values of the exponents for all shells (packed into a 1-d array)
        read (LUCOEF,err=20,end=20) priexp
 !      values of the contraction coefficients for all shells (packed into a 1-d array)
        read (LUCOEF,err=20,end=20) priccf
 
        !write(inf_mp2no) nao_bas,nshells,npriexp,ncnt 
        !WRITE(inf_mp2no) nhkt,ncent,nuco
        !WRITE(inf_mp2no) jstrt,numcf
        !WRITE(inf_mp2no) cent
        !WRITE(inf_mp2no) priexp
        !WRITE(inf_mp2no) priccf

        !write(inf_mp2no,*) nao_bas,nshells,npriexp,ncnt 
        !WRITE(inf_mp2no,*) nhkt,ncent,nuco
        !WRITE(inf_mp2no,*) jstrt,numcf
        !WRITE(inf_mp2no,*) cent
        !WRITE(inf_mp2no,*) priexp
        !WRITE(inf_mp2no,*) priccf

        ! put in gtos
        ao_basis%nshells       = nshells
        ao_basis%nao           = nao_bas
        ao_basis%basis_angular = 1
        nullify(ao_basis%gtos)
        allocate(ao_basis%gtos(nshells))
        do ishell = 1, nshells
           ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
           ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
           ao_basis%gtos(ishell)%n_primitives = nuco(ishell)
           nullify  (ao_basis%gtos(ishell)%exponent) 
           allocate (ao_basis%gtos(ishell)%exponent(nuco(ishell)))
           nullify  (ao_basis%gtos(ishell)%coefficient)
           allocate (ao_basis%gtos(ishell)%coefficient(nuco(ishell)))
           iPrimStart = jstrt(ishell) + 1
           iPrimEnd   = jstrt(ishell) + nuco(ishell)
           k          = numcf(ishell)
           ao_basis%gtos(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
           ao_basis%gtos(ishell)%coefficient = &                                   
                        priccf(iPrimStart:iPrimEnd,k)
           ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
 !         In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
           call compress_unc_shell (ao_basis%gtos(ishell),ierr)
           if (ierr == 1) ierr = CONTRAC
        end do
 
        deallocate(ncent)
        deallocate(nhkt)
        deallocate(nuco)
        deallocate(numcf)
        deallocate(jstrt)
        deallocate(cent)
        deallocate(priexp)
        deallocate(priccf)
 
        ! make array with pointers from a basis function index to its shell
        nullify(ao_basis%shell_indices)
        allocate(ao_basis%shell_indices(nao))
        iao = 0
        do ishell = 1, nshells
           do ibas = 1,&       
                nfunctions(&       
                ao_basis%gtos(ishell)%orb_momentum,ao_basis%basis_angular)
              iao = iao + 1
              if (iao.gt.nao) then
                 ierr = INCONSIS
                 call quit ('Error in read_from_dirac')
              endif
              ao_basis%shell_indices(iao) = ishell
           end do
        end do
        if (iao.gt.nao) then
           ierr = INCONSIS
           close (lucoef,status='KEEP')
           close (inf_mp2no,status='KEEP')
           call quit ('Error in read_from_dirac')
        endif
 
        close (lucoef,status='KEEP')
        close (inf_mp2no,status='KEEP')
 
        call print_date("read from DFCOEF:"//text)
        return
 
   10   ierr = CORRUPT
        close (lucoef,status='KEEP')
        return
 
   20   ierr = NOBASIS
        close (lucoef,status='KEEP')
        return


        end subroutine MP2_read_write_from_dirac



!------------------------------------------------------------------------------------------------------------------




        subroutine get_MP2_fockmatrix (nocc,nvir,mo_occ,mo_vir,diag_mat_occ,diag_mat_vir,int_t, print_level)

!        Written by Lucas Visscher, summer 2017
!        Revised by Xiang Yuan, July 2020
!        Maybe used one day.....
         use exacorr_global

         implicit none

         integer, intent(in)                  :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)                  :: print_level
         integer, intent(in)                  :: mo_occ(:),mo_vir(:) ! the orbital energies
   


!        target tensors to be filled and given back to the caller
         type(talsh_tens_t), intent(inout)    :: diag_mat_occ, diag_mat_vir
         type(talsh_tens_t)                   :: fo_tensor, fv_tensor
         type(talsh_intg_tens), intent(inout) :: int_t

         type(talsh_tens_t)              :: id_mat_occ
         type(talsh_tens_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2)
         integer(INTD)                   :: ierr
         logical                         :: one_el_exist
         real(8), allocatable            :: eps_occ(:),eps_vir(:) ! orbital energies
         type(one_el_t)                  :: one_int
         real(8)                         :: e_core
         integer                         :: nkp, min_occ
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer                         :: i, j
         type(talsh_tens_t)              :: result_tensor
         integer(INTD)                   :: result_dims(1)
         complex(8), pointer             :: result_tens(:)


         real(8), parameter :: THRESHOLD=1.0D-8

         min_occ=minval(mo_occ)
         result_dims(1)=1

!        Get foo tensor
         fo_dims(1:2) = nocc
         ierr=talsh_tensor_construct(fo_tensor,C8,fo_dims,init_val=ZERO)

!        Get fvv tensor
         fv_dims(1:2) = nvir
         ierr=talsh_tensor_construct(fv_tensor,C8,fv_dims,init_val=ZERO)


!        get one electron integrals
         one_el_exist = exist_one_el(e_core)

!        write one electron integrals to tensor
         if (one_el_exist) then

!           create density (identity) matrix for the occupied space
            ierr=talsh_tensor_construct(id_mat_occ,C8,fo_dims,init_val=ZERO)
            ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, nocc
              h_tens(i,i) = ONE
            enddo

!           now write one electron integrals

            ierr=talsh_tensor_get_body_access(fo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_occ,nocc,min_occ)


            ierr=talsh_tensor_get_body_access(fv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_vir,nvir,mo_vir,nvir,min_occ)

!           Get access to result
            result_dims(1) = 1
            ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
            ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,result_tens,result_dims)

            call print_date('one electron integrals added to Fock Matrix')

            ierr=talsh_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",fo_tensor,int_t%oooo,id_mat_occ, ONE)
            ierr=talsh_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",fv_tensor,int_t%vovo,id_mat_occ, ONE)
            call print_date('two electron integrals added to Fock Matrix')
            
            if (print_level >10) then
              call print_tensor(fo_tensor, 1.0D-6, 'fo_tensor')
              call print_tensor(fv_tensor, 1.0D-6, 'fv_tensor')
            end if

            !now we can get the correct epsilon values
            ierr=talsh_tensor_get_body_access(fo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            do i=1,fo_dims(1)
              eps_occ(i)=h_tens(i,i)
            end do

            ierr=talsh_tensor_get_body_access(fv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            do i=1,fv_dims(1)
              eps_vir(i)=h_tens(i,i)
            end do


            ierr=talsh_tensor_construct(diag_mat_occ,C8,fo_dims,init_val=ZERO)
            ierr=talsh_tensor_get_body_access(diag_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, nocc
                h_tens(i,i) = eps_occ(i)
            end do

            ierr=talsh_tensor_construct(diag_mat_vir,C8,fv_dims,init_val=ZERO)
            ierr=talsh_tensor_get_body_access(diag_mat_vir,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, nvir
                h_tens(i,i) = eps_vir(i)
            end do


            ierr=talsh_tensor_destruct(id_mat_occ)
            ierr=talsh_tensor_destruct(result_tensor)
       
         if (print_level >0) call print_date('Fock matrix initialized')
        end if 
        end subroutine get_MP2_Fockmatrix
end module talsh_mp2no
        !------------------------------------------------------------------------------------------------
