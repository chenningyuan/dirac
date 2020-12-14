module matrix_complex_to_quaternion
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

        public talsh_mp2no_driver,talsh_create_AO_to_MNO_trans

!        interface scale_with_denominators
!           module procedure scale_with_denominators
!        end interface
!

       contains
!-------------------------------------------------------------------------------------------------------------------------------------

       subroutine read_matrix (Matrix, Matrix_name)
!       Written by Xiang Yuan, Dec 2020.

        use davidson
        complex(kind=8), pointer, intent(inout)  :: Matrix(:,:)   ! MP2 density matrix from DM_complex
        character(len=*),intent(in)              :: Matrix_name




!--------------------------------------------Test part-------------------------------------------------------------------------------------------
!        print *,'Occupation numbers of all natural orbital'
!       do i = 1, size(NO_OCC_1,1)
!          print "(4x,i4,a,2e14.6)",i,"    ",NO_OCC_1(i)
!       end do
!       deallocate(NO_OCC_1)

!      end subroutine read_matrix

 !-----------------------------------------------------------------------------------------------------------------------------------------------

      subroutine sort_values_vectors (sortingorder, NOs_OCC_quat, NOs_MO_quat, threshold)

!       This subroutine sorts and choose the natural orbital used in following computation 
!       Written by Xiang Yuan, July 2020

         use exacorr_global
         use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
         use talsh_gradient
         use davidson

         type(exacc_input), intent(in) :: exa_input

         real(kind=8), pointer       :: NO_OCC_quat(:)          ! occupation number from qdia of DM
         real(kind=8), pointer       :: NOs_MO_quat(:,:,:)      ! natural orbital from qdia of DM
         real(kind=8)                :: threshold

         complex(kind=8), pointer,intent(inout)     :: NOs_MO(:,:)

!        Point connected with tensors 
         integer(C_INT)              :: ierr
         real(kind=8), allocatable   :: NOs_MO_R(:,:) ! auxilary array of quaternion natural orbital used in sort routine        
         real(kind=8), allocatable   :: NOs_MO_L(:,:) ! auxilary array of quaternion natural orbital used in sort routine   
         real(kind=8), allocatable   :: NO_OCC(:,:)   ! auxilary array of occupation number used in sort routine 

!        Scalars (need to be defined as tensor types)
         real(kind=8)  :: energy_shift(2) = 0.0d0
         integer       :: i, j, k, n_mp2nos
         integer       :: num_no_red       ! number of reduced natural orbital
         integer       :: sorting_order = -1
         logical       :: left  = .true.  !variable in davidson module controling sorted routine 
         logical       :: right = .true.   !variable in davidson module controling sorted routine 
         logical       :: apply_evalues_shift = .false.  !variable in davidson module controling sorted routine 


         n_mp2nos = size(NOs_MO_quat,1) ! MP2NOs is a square matrix, so we only need one dimension


!  Make some noise so that we know we are proceeding 
         call print_date('Sort natural orbital (using RELCCSD davidson auxiliary routines)')

! temporaryly using the relccsd davidson sort_values_vectors, which in this case
! means we need to setup the complex mode in the davidson module, and with that
! we can handle complex eigenvalues and eigenvectors
         call davidson_setup(complex_mode=.true.,unit_output=6,solve_right=.false.)

! allocate first dimension as 2*n_mp2nos since representing the complex numbers as pairs of reals
         print *, "n_mp2no", n_mp2nos
            allocate(NOs_MO_R(n_mp2nos,n_mp2nos)) 
            allocate(NOs_MO_L(n_mp2nos,n_mp2nos))
            allocate(NO_OCC(2,size(NO_OCC_quat,1))) 

! quaternions have 4 component, so we sort two times, sort 1, 2 component 

            NO_OCC = 0.0d0
            NO_OCC(1,:) = NO_OCC_quat(:)
            NO_OCC(2,:) = NO_OCC_quat(:)


            do i = 1, n_mp2nos
               do k = 1, n_mp2nos
                  NOs_MO_R(i,k) = (NOs_MO_quat(i,k,1))
                  NOs_MO_L(i,k) = (NOs_MO_quat(i,k,2))
               enddo
            enddo


            call sort_values_vectors(sorting_order, NO_OCC, NOs_MO_R, NOs_MO_L, &
                                     right,left,apply_evalues_shift,energy_shift)


            do i = 1, n_mp2nos
               do k = 1, n_mp2nos
                     NOs_MO_quat(i,k,1) = (NOs_MO_R(i,k))
                     NOs_MO_quat(i,k,2) = (NOs_MO_L(i,k))
               enddo
            enddo

!  sort 3, 4 component 

            NO_OCC = 0.0d0
            NO_OCC(1,:) = NO_OCC_quat(:)
            NO_OCC(2,:) = NO_OCC_quat(:)

            do i = 1, n_mp2nos
               do k = 1, n_mp2nos
                     NOs_MO_R(i,k) = (NOs_MO_quat(i,k,3))
                     NOs_MO_L(i,k) = (NOs_MO_quat(i,k,4))
               enddo
            enddo


            call sort_values_vectors(sorting_order, NO_OCC, NOs_MO_R, NOs_MO_L, &
                                     right,left,apply_evalues_shift,energy_shift)


            do i = 1, n_mp2nos
               do k = 1, n_mp2nos
                     NOs_MO_quat(i,k,3) = (NOs_MO_R(i,k))
                     NOs_MO_quat(i,k,4) = (NOs_MO_L(i,k))
               enddo
            enddo


            NO_OCC_quat(i) = NO_OCC_(2,i)

!           call davidson_print_vectors(NO_OCC,   label="MP2 NOs occupation numbers    ")
!           call davidson_print_vectors(NOs_MO_R, label="MP2 Natural orbitals (RHS)    ", complex_dimension=1)

            deallocate(NOs_MO_R)
            deallocate(NOs_MO_L)
            deallocate(NOs_OCC)            

!            NO_OCC_1(:) = NO_OCC(1,:)
!            deallocate(NO_OCC_quater)

        end subroutine sort_values_vectors

!-------------------------------------------------------------------------------------------------------------------------------------------

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
!           print "(4x,a)", ""
!           print "(4x,a)", "Eigenvectors"
!           do i = 1, ndim
!              print "(4x,a,$)", ""
!              do k = 1, ndim
!                 print "(a,2e14.6,$)","    ",vectors(i,k)
!              enddo
!              print *, ""
!           enddo
        end subroutine

!----------------------------------------------------------------------------------------------------------------------    

        subroutine basis_transform (exa_input,initialize_talsh,threshold,nocc,nvir,mo_vir,NOs_MO,NO_OCC,NOs_AO_tens)
!        Routine to set theshold and reduce natural orbital space including transform from MO basis to AO basis 
!        Now theshold is 0
                
           type(exacc_input), intent(in)          :: exa_input
           integer, intent(in)                    :: nocc,nvir      ! the size of the mo basis for occupied and virtual spinors
           integer, intent(in)                    :: mo_vir(:)      ! the list of virtual orbitals
           integer, intent(in)                    :: initialize_talsh         ! key to initial the lib
           complex(kind=8), pointer, intent(in)   :: NOs_MO(:,:,:)    ! sorted natural orbital in MO basis
           real(kind=8), intent(in)               :: threshold
           real(kind=8), pointer,intent(in)       :: NO_OCC(:)         ! natural orbital occupation number 
           complex(kind=8), pointer,intent(out)   :: NOs_AO_tens(:,:)
           integer(C_INT)                     :: ierr
           type(C_PTR)                        :: body_p
           integer(C_SIZE_T)                  :: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
           integer(C_INT)                     :: host_arg_max
           type(talsh_tens_t)                 :: mocoefa_tensor, mocoefb_tensor ! HF molecular orbital in AO basis 
           type(talsh_tens_t)                 :: NOs_MOa_red_tensor  ! reduced natural orbital spin_1 in MO basis
           type(talsh_tens_t)                 :: NOs_MOb_red_tensor  ! reduced natural orbital spin_2 in MO basis
           type(talsh_tens_t)                 :: NOs_AO_tensor  ! reduced natural orbital in AO basis
           type(talsh_tens_t)                 :: NOs_AOb_tensor  ! reduced natural orbital in AO basis
           complex(kind=8), pointer           :: NOs_MO_red(:,:)
           complex(kind=8), pointer           :: show_result(:,:)
           integer, allocatable               :: mp2no_mo_list(:)   ! molecular orbital indices
           integer                            :: num_no_red
           integer(INTD), dimension(2)        :: no_red_dims, no_ao_dims
           integer(INTD)                      :: mocoef_dims(2)
           integer                            :: nao, nmo     ! the size of the ao & mo basis for natural orbital transformation 
           integer                            :: i,j
                      
                
!       Initialize libraries
                
           if (initialize_talsh == 0) then
              buf_size=exa_input%talsh_buff*buf_size
              ierr=talsh_init(buf_size,host_arg_max)
              call print_date('Initialized talsh library')
              write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max
              call interest_initialize(.true.)
              call print_date('Initialized interest library')
           end if
                
           num_no_red = 0
           do i=1,size(NO_OCC,1),1
              if (NO_OCC(i) > threshold) then
                 num_no_red = num_no_red +1
              end if 
           end do

!           print *, size(NOs_MO,1)
!           print *, size(NO_OCC,1)
           no_red_dims(1)=nvir
           no_red_dims(2)=num_no_red
!          print *,'size no_red_dims(1)',no_red_dims(1)
!          print *,'size no_red_dims(2)',no_red_dims(2)
           ierr=talsh_tensor_construct(NOs_MOa_red_tensor,C8,no_red_dims,init_val=ZERO)
           ierr=talsh_tensor_get_body_access(NOs_MOa_red_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
           call c_f_pointer(body_p,NOs_MO_red,no_red_dims) ! to use as a regular Fortran array
           do i=1,size(NOs_MO_red,2),1
                do j=1,size(NOs_MO_red,1),1            
                        NOs_MO_red(j,i)=NOs_MO(j,i+nvir-num_no_red,1)
                end do 
           end do

           ierr=talsh_tensor_construct(NOs_MOb_red_tensor,C8,no_red_dims,init_val=ZERO)
           ierr=talsh_tensor_get_body_access(NOs_MOb_red_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
           call c_f_pointer(body_p,NOs_MO_red,no_red_dims) ! to use as a regular Fortran array
           do i=1,size(NOs_MO_red,2),1
                do j=1,size(NOs_MO_red,1),1            
                        NOs_MO_red(j,i)=NOs_MO(j,i+nvir-num_no_red,2)
                end do 
           end do

!       print *,'NOs_MO_red',ierr
                  
        
           call print_date('Extracting virtual MO coeffcient in AO basis set')
!       Extracting virtual MO coeffcient in AO basis set
           nao = (nocc + nvir)/2
           nmo = nvir 
           mocoef_dims(1) = nao
           mocoef_dims(2) = nmo
           allocate(mp2no_mo_list(nvir))
           mp2no_mo_list = mo_vir
           ierr=talsh_tensor_construct(mocoefa_tensor,C8,mocoef_dims,init_val=ZERO)
           ierr=talsh_tensor_construct(mocoefb_tensor,C8,mocoef_dims,init_val=ZERO)
           call get_MP2_mocoef_tensor (nao, nmo, mp2no_mo_list, mocoefa_tensor, mocoefb_tensor) 
!           ierr=talsh_tensor_add("moco_a(a,b)+=moco_b(a,b)",mocoefa_tensor,mocoefb_tensor)
          Print *,'Extracting virtual MO coeffcient in AO basis set',ierr
           ierr=talsh_tensor_get_body_access(mocoefa_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
           call c_f_pointer(body_p,show_result,mocoef_dims)
           print *,'coeff of HF Molecular Orbital'
           do i = 1, nao
              print "(4x,a,$)", ""
              do j = 1, nmo
                 print "(a,2e14.6,$)","    ",show_result(i,j)
              enddo
              print *, ""
           enddo
                
                
           call print_date('Transform reduced natural orbital from MO basis into AO basis')
!          Transform reduced natural orbital from MO basis into AO basis
           no_ao_dims(1)=nao
           no_ao_dims(2)=num_no_red
           ierr=talsh_tensor_construct(NOs_AO_tensor,C8,no_ao_dims,init_val=ZERO)
           ierr=talsh_tensor_construct(NOs_AOb_tensor,C8,no_ao_dims,init_val=ZERO)
           ierr=talsh_tensor_get_body_access(NOs_AO_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
           ierr=talsh_tensor_contract("NO_AO(a,b)+=MO_AO(a,i)*NO_MO(i,b)",NOs_AO_tensor,mocoefa_tensor,NOs_MOa_red_tensor,scale=ONE)
          ierr=talsh_tensor_contract("NO_AO(a,b)+=MO_AO(a,i)*NO_MO(i,b)",NOs_AOb_tensor,mocoefb_tensor,NOs_MOb_red_tensor,scale=ONE)
           ierr=talsh_tensor_add("NO_AO(a,b)+=NO_AO2(a,b)",NOs_AO_tensor,NOs_AOb_tensor)
           call c_f_pointer(body_p,NOs_AO_tens,no_ao_dims) ! to use as a regular Fortran array
           print *,'coeff of natural Orbital in AO basis'
           do i = 1, nao
              print "(4x,a,$)", ""
              do j = 1, num_no_red
                 print "(a,2e14.6,$)","    ",NOs_AO_tens(i,j)
              enddo
              print *, ""
           enddo
                          
        end subroutine basis_transform
!-----------------------------------------------------------------------------------------------------------------------
                                                  
        subroutine diagoanl_matrix_quater(A,EIG,VEC)

         ! Set parameter needed in qdiag
                 
                 real(kind=8), intent(in), pointer :: A(:,:,:)
                 real(kind=8), pointer :: EIG(:)
                 real(kind=8), pointer :: VEC(:,:,:)                 
                 integer :: LRA, LCA, LRV, LCV
                 integer :: NZ, N, MATZ, LWORK, IERR
                 integer :: i, j
                 real(kind=8), allocatable :: WORK(:)
                 
               
                 LRA = size(A,1)
                 print *, "LRA?", LRA
                 LCA = size(A,2)
                 print *, "LCA?", LCA
                 LRV = LRA
                 print *, "LRV?", LRV
                 LCV = LCA
                 print *, "LCV?", LCV
                 lwork = 10000000*LRA
                 NZ = 4
                 N = LRA
                 MATZ = 1
                 allocate(EIG(N))
                 allocate(VEC(LRV,LCV,NZ))
                 allocate(work(MAX(1,LWORK)))


                 call qdiag (NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,&
                             WORK,LWORK,IERR)

                  print *, "occupation number in quaternion"
                  do i = 1, LRA
                     print *, EIG(i)
                  end do
                  deallocate(work)
   
                print *, "it works in quaternion_to_complex?"
      end subroutine diagoanl_matrix_quater
                
                
!---------------------------------------------------------------------------------------------------------------------------------------------------
            
      subroutine complex_to_quaternion(exa_input,dm_mo_fort,dm_diracformat)
         
               use exacorr_global, only : get_nsolutions

               !input variables
               type(exacc_input),   intent(in)  :: exa_input
               !matrix dimensions
               integer :: nocc, nvir, nkr_occ, nkr_vir, nesh
               complex(8), pointer         :: dm_mo_fort(:,:)
               real(8),pointer             :: dm_diracformat(:,:,:)
               !complete mo list
               integer, allocatable :: mo_list(:)
               !error code
               integer :: ierr
               !loop variables
               integer :: i, j, norbt
               integer :: i_unbar, j_unbar, j_bar
               !debug
               logical :: debug=.false.
               integer :: k


!***************
!     set dims *
!***************
 
      nocc    = exa_input%nocc
      nvir    = exa_input%nvir
      nkr_occ = exa_input%nkr_occ
      nkr_vir = exa_input%nkr_vir


!*****************************************************************************************************************************
!     Determine where this should go in the full list of MOs (NB: DIRAC assumes restricted densities and uses Kramers pairs) *
!*****************************************************************************************************************************
      if (debug) then
        print *,'<<<<< DM, complete, before formatting to dirac'
        k = 1
        do i=1,(nocc+nvir)
           do j=1,(nocc+nvir)
              write(*,'(3X,I4,2x,2F28.10,2x,4I4)') k, dm_mo_fort(i,j),i,j
              k = k + 1
           end do
         end do
         print *,'>>>>> DM, complete, before formatting to dirac'
      end if

!      norbt = nkr_occ + nkr_vir
      norbt = nkr_vir
      print *,"norbt", norbt 
      allocate(mo_list(norbt))
      mo_list(1:nkr_occ)                 = exa_input%mokr_occ
      mo_list(1+nkr_occ:nkr_occ+nkr_vir) = exa_input%mokr_vir

!      nesh = get_nsolutions() / 2
      nesh = norbt
      print *,"nesh", nesh 
      allocate (dm_diracformat(nesh,nesh,4))
      dm_diracformat = real(ZERO)
      do j = 1, norbt
        j_unbar = 2*j-1
        j_bar   = 2*j
        do i = 1, norbt
          ! We copy the unbar-unbar part of the density matrix, as well as the unbar-bar part
          ! This follows from eq. 27 of Shee, Visscher, Saue JCP 145 184107, 2016
          ! it would be good to test whether the Kramers-restricted assumption holds
          ! Or even to average over the two density matrices and enforce a Kramers-restricted density
          i_unbar = 2*i-1
          dm_diracformat(mo_list(i),mo_list(j),1) = real(dm_mo_fort(i_unbar,j_unbar))
          dm_diracformat(mo_list(i),mo_list(j),2) = imag(dm_mo_fort(i_unbar,j_unbar))
          dm_diracformat(mo_list(i),mo_list(j),3) = real(dm_mo_fort(i_unbar,j_bar))
          dm_diracformat(mo_list(i),mo_list(j),4) = imag(dm_mo_fort(i_unbar,j_bar))
        end do
      end do

!      deallocate (dm_diracformat)
      deallocate (mo_list)
!      deallocate (dm_mo_fort)
      end subroutine complex_to_quaternion
 !-------------------------------------------------------------------------------------------------
              
end module matrix_complex_to_quaternion

