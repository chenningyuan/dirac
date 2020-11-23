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

       subroutine talsh_mp2no_driver (exa_input)
!       In principle, current code will generate two files:
!       MP2NOs_MO: MP2 natural orbital in MO basis & MP2NOs_AO: MP2 natural orbital in AO basis
!       Written by Xiang Yuan, Oct 2020.
        use exacorr_global
        use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
        use talsh_gradient
        use davidson
        type(exacc_input), intent(in) :: exa_input
        integer                :: nocc, nvir   ! the size of the mo basis for occupied and virtual spinors
        integer,  allocatable  :: mo_vir(:)    ! the list index of virtual orbitals

!       Point connected with tensors 
        integer(C_INT)              :: ierr
        real(kind=8), pointer       :: NO_OCC_1(:)         ! natural orbital occupation number saved in new files
        complex(kind=8), pointer    :: NOs_MO(:,:)  ! sorted natural orbital in MO basis
        complex(kind=8), pointer    :: NOs_AO_tens(:,:)  ! natural orbital in AO basis
        type(basis_set_info_t)      :: ao_basis ! information about the ao basis
        type(cmo)                   :: cspinor_all ! Stores all MO-coefficients and energies obtained from the parent code
        character                   :: text*74      ! holds text read from DFCOEF, only used for printing this
        integer, parameter          :: CORRUPT  = 2 ! file unreadable
        integer                     :: nsym
        integer                     :: mp2out
        integer                     :: err
        integer                     :: initialize_talsh
        logical                     :: mp2no_present 

!       Scalars (need to be defined as tensor types)
!       Normally, we don't generate MP2NOs_MO, but only MP2NOs_AO
        real(kind=8)       :: threshold=0.0d0
        integer            :: i,j

        nocc = exa_input%nocc
        nvir = exa_input%nvir
        mo_vir = exa_input%mo_vir
        INQUIRE(FILE="MP2NOs_MO",EXIST=mp2no_present)
        if (mp2no_present) then
           initialize_talsh = 0
           open(mp2out, file = 'MP2NOs_MO', status='OLD',FORM='UNFORMATTED',&     
           access='SEQUENTIAL',iostat=err) 
           allocate(NOs_MO(nvir,nvir))
           allocate(NO_OCC_1(nvir))
           read(mp2out,iostat=err) text,nsym
           read(mp2out,iostat=err) NOs_MO
           read(mp2out,iostat=err) NO_OCC_1
           print *, 'foobaz'

        else
           print *, 'foobar'
           initialize_talsh=1
! generate MP2NOs_MO file containing all natural orbital and occupation number in MO basis
           call talsh_generate_MP2_VNO_MObasis (exa_input,NOs_MO,NO_OCC_1,nocc,nvir,mo_vir)   
        end if 

!--------------------------------------------Test part-------------------------------------------------------------------------------------------
       print *,'MP2 virtual spinor occupation numbers'
       do i = 1, size(NO_OCC_1,1)
          print "(4x,i4,a,2e14.6)",i,"    ",NO_OCC_1(i)
       end do

!      call print_values_R_vectors_C(NO_OCC_1,NOs_MO)
       deallocate(NO_OCC_1)

      end subroutine talsh_mp2no_driver

 !-----------------------------------------------------------------------------------------------------------------------------------------------

      subroutine talsh_generate_MP2_VNO_MObasis (exa_input,NOs_MO,NO_OCC_1,nocc,nvir,mo_vir)

!       This subroutine obtains the mp2 natural orbital in MO basis in complex representation 
!       Written by Xiang Yuan, July 2020

         use exacorr_global
         use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
         use talsh_gradient
         use davidson

         type(exacc_input), intent(in) :: exa_input

         complex(kind=8), pointer,intent(inout)     :: NOs_MO(:,:)
         real(kind=8), pointer,intent(inout)        :: NO_OCC_1(:)         ! natural orbital occupation number and the one saved in new files
         integer,intent(inout)                      :: nocc, nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable,intent(inout)        :: mo_vir(:) ! the list of virtual orbitals

         integer,  allocatable  :: mo_occ(:)               ! the list of occupied orbitals
         real(8), allocatable   :: eps_occ(:),eps_vir(:) ! HF orbital energies

!        Point connected with tensors 
         integer(C_INT)              :: ierr
         complex(8), pointer         :: MP2_density_tens(:,:)
         real(kind=8), pointer       :: NOs_MO_R(:,:) ! natural orbitals, complex values represented as reals 
         complex(kind=8), pointer :: test_tens(:,:)
        
         real(kind=8), pointer :: NO_OCC(:,:)       ! natural orbital occupation number and the one saved in new files
         complex(kind=8), pointer :: NOs_AO_tens(:,:)  ! natural orbital in AO basis
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
         integer(INTD), dimension(4) :: vvoo_dims
         integer(INTD), dimension(2) :: vv_dims
         integer(C_INT)              :: one_dims(1)
         type(basis_set_info_t)      :: ao_basis ! information about the ao basis
         type(cmo)                   :: cspinor_all ! Stores all MO-coefficients and energies obtained from the parent code

!        Scalars (need to be defined as tensor types)
         real(kind=8)  :: threshold
         real(kind=8)  :: energy_shift(2) = 0.0d0
         integer       :: i,j, k, n_mp2nos
         integer       :: num_no_red       ! number of reduced natural orbital
         integer       :: sorting_order = -1
         logical       :: left  = .false.  !variable in davidson module controling sorted routine 
         logical       :: right = .true.   !variable in davidson module controling sorted routine 
         logical       :: apply_evalues_shift = .false.  !variable in davidson module controling sorted routine 


!        Make some noise so that we know we are here
         call print_date('Entered mp2_driver routine')

!        Initialize libraries
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)
         call print_date('Initialized talsh library')
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

         call print_date('Initializing interest library')
         call interest_initialize(.true.)
         call print_date('Initialized interest library')

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         

         print *,'number of virtual orbital', nvir
         print *,'number of occupied orbital', nocc
         allocate(mo_occ(nocc))
         mo_occ = exa_input%mo_occ


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

!  Make some noise so that we know we are proceeding 
         call print_date('Forming MP2 VV density matrix')
!        Initialize scale with denominators of oovv integrals
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=int_t%oovv)

!        Calculate and print the MP2 density matrix of virtual-virtual block Dvv 
         ierr=talsh_tensor_init(MP2_density_tensor%prime_vv)
         ierr=talsh_tensor_contract("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_density_tensor%prime_vv,&
                                        int_t%oovv,int_t%oovv,scale=ONE_HALF) ! o(2)v(2)
!        print *,'after constructing density matrix',ierr

         call write_DenMatrix_complex_into_DMComplefile(MP2_density_tens)

!  Make some noise so that we know we are proceeding 
         call print_date('Diagonalizing MP2 VV density matrix')
!        Diagoliazation of Density matrix
         allocate(NO_OCC(2,size(MP2_density_tens,2)))
         allocate(NO_OCC_1(size(MP2_density_tens,2)))
         NO_OCC = 0.0d0
         call diagoanl_DMP2_matrix(MP2_density_tens,NO_OCC,NOs_MO)

         n_mp2nos = size(NOs_MO,1) ! MP2NOs is a square matrix, so we only leed one dimension

!        print *,'unsorted MP2 NOs and occupation numbers'
!        call print_values_R_vectors_C(NO_OCC,NOs_MO)

!  Make some noise so that we know we are proceeding 
         call print_date('Sort natural orbital (using RELCCSD davidson auxiliary routines)')

! temporaryly using the relccsd davidson sort_values_vectors, which in this case
! means we need to setup the complex mode in the davidson module, and with that
! we can handle complex eigenvalues and eigenvectors
         call davidson_setup(complex_mode=.true.,unit_output=6,solve_right=.false.)

! allocate first dimension as 2*n_mp2nos since representing the complex numbers as pairs of reals
            allocate(NOs_MO_R(2*n_mp2nos,n_mp2nos)) 

! copy over the MP2 NOs from the original array of complex variables to the auxiliary real array
            do i = 1, n_mp2nos
               j = 2*i - 1 
               do k = 1, n_mp2nos
                  NOs_MO_R(j,k)   = dreal(NOs_MO(i,k))
                  NOs_MO_R(j+1,k) = dimag(NOs_MO(i,k))
               enddo
            enddo

            call sort_values_vectors(sorting_order, NO_OCC, NOs_MO_R, NOs_MO_R, &
                                     right,left,apply_evalues_shift,energy_shift)

!           call davidson_print_vectors(NO_OCC,   label="MP2 NOs occupation numbers    ")
!           call davidson_print_vectors(NOs_MO_R, label="MP2 Natural orbitals (RHS)    ", complex_dimension=1)

!  as we had to use a representation of complex numbers as pairs of reals, we reconstruct now the actual array of complex variables
            do i = 1, n_mp2nos
               j = 2*i - 1
               do k = 1, n_mp2nos
                  NOs_MO(i,k) = dcmplx(NOs_MO_R(j,k),NOs_MO_R(j+1,k))
               end do
            enddo
            deallocate(NOs_MO_R)

!           print *,'sorted MP2 NOs and occupation numbers'
!           call print_values_R_vectors_C(NO_OCC,NOs_MO)

            NO_OCC_1(:) = NO_OCC(1,:)
            deallocate(NO_OCC)

        end subroutine talsh_generate_MP2_VNO_MObasis

!--------------------------------------------------------------------------------------------------------------------------------------------------------------
                
        subroutine print_values_R_vectors_C(values,vectors)
           real(kind=8),    intent(in) :: values(:,:)
           complex(kind=8), intent(in) :: vectors(:,:)
           integer :: ndim, i, k
   
           ndim = size(values,2)
 
           print "(4x,a)", "Eigenvalues"
           print "(4x,a,$)", ""
           do i = 1, ndim
              print "(a,2e14.6,$)","    ",values(:,i)
           enddo
           print "(4x,a)", ""
           print "(4x,a)", "Eigenvectors"
           do i = 1, ndim
              print "(4x,a,$)", ""
              do k = 1, ndim
                 print "(a,2e14.6,$)","    ",vectors(i,k)
              enddo
              print *, ""
           enddo
        end subroutine

!--------------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine diagoanl_DMP2_matrix(M,eVal,eVectorsR)

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
  
        end subroutine diagoanl_DMP2_matrix
                                                        
                                              
!------------------------------------------------------------------------------------------------------------------------------------------------------------------

      subroutine write_DenMatrix_complex_into_DMComplefile(DM)
         complex(8), pointer         :: DM(:,:)
         integer                     :: mp2dm,err
         open(mp2dm, file = 'MP2DMs', status='NEW',FORM='UNFORMATTED',&     
         access='SEQUENTIAL',iostat=err) 
         write(mp2dm) DM
         write(*,*) err
      end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine get_MP2_integrals (nocc,nvir,mo_occ,mo_vir,int_t,print_level)
                
                  !        Routine to get oovv integral in the form of antisymmetric tensors.
                  !        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
                  !        full integral list.
                  !        Revised by Xiang Yuan, July 2020
                  
                           use exacorr_global
                           use talsh_ao_to_mo
                  
                           integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
                           integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
                  !        target integral tensors to be filled and given back to the caller
                           type(talsh_intg_tens),intent(inout)  :: int_t
                           integer, intent(in)                  :: print_level 
                           integer(INTD)                        :: oovv_dims(1:4)
                           integer(INTD)                        :: ierr
                           complex(8), pointer, contiguous      :: ao_tens(:,:,:,:)
                           type(C_PTR)                          :: body_p 
                           integer                              :: nao
                  
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
                           nao = get_nao()      ! number of basis functions
                  
                  !        Get tensor with AO integrals
                           aoint_dims = nao
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
                 

end module talsh_mp2no

