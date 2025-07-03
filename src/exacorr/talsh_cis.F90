module talsh_cis
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use talsh_cis_property
    use talsh_cis_density
    

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor
    private

    public talsh_cis_driver

    contains
!---------------------------------------------------------------------------------------------------------------

    subroutine talsh_cis_driver (exa_input,CIS_vv_effden,initialize_tensor_library,initialize_interest)
    ! This routine drives CIS calculation including construction of CIS Hamiltonian and diagnolization
    ! Written by Xiang Yuan, June 2020
        use exacorr_global
        use talsh_ao_to_mo
        use talsh_rt_solver


        type(exacc_input), intent(in) :: exa_input
        type(talsh_tens_t), optional   :: CIS_vv_effden
        logical, optional :: initialize_tensor_library
        logical, optional :: initialize_interest

         integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         integer(C_INT)       :: ierr
         complex(8), pointer :: result_tens(:,:,:,:),sigma1(:,:),check_tens(:,:)
         integer(INTD)       :: fo_dims(1:2),fv_dims(1:2)
         complex(kind=8), pointer :: eValues(:) 
         complex(kind=8), pointer :: eVectorsR(:,:)
         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer :: i,j,a,b,n
         integer :: bj,ai
         logical :: cisrt


!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens_t) :: int_t
!        single-body and body_two tensors
         type(talsh_tens_t) :: body_one_tensor, diag_foo_tensor, diag_fvv_tensor
         type(talsh_tens_t) :: id_mat_occ, id_mat_vir
!        intermediate tensors (same notation as in relccsd)
         integer(INTD), dimension(2) :: oo_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims, vovo_dims, voov_dims
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor,sigma_tensor
         integer(C_INT)     :: one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CIS control and result variables
         real(8) :: tem
         integer :: ncycles
         integer :: tar_state
         logical :: initialize_tensor_library_local
         logical :: initialize_interest_local

         initialize_tensor_library_local = .true.
         initialize_interest_local  = .true.

         if (present(initialize_tensor_library)) initialize_tensor_library_local = initialize_tensor_library 
         if (present(initialize_interest)) initialize_interest_local = initialize_interest
         
         if (initialize_tensor_library_local) then 
!        Initialize libraries
            buf_size=exa_input%talsh_buff*buf_size
            ierr=talsh_init(buf_size,host_arg_max)
            call print_date('Initialized talsh library in CIS driver')
            write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max
         else
            call print_date('NOTE: Tensor libraries NOT initialized in CIS driver')
         end if

         call interest_initialize(initialize_interest_local)
         if (initialize_interest_local) then 
            call print_date('Initialized Interest library in CIS driver')
         else
            call print_date('NOTE: Interest library NOT initialized in CIS driver')
         end if
        
!        DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         allocate(sigma1(nocc*nvir,nocc*nvir))
         allocate(eValues(nocc*nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir


!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         oo_dims(1:2) = nocc
         vv_dims(1:2) = nvir
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         voov_dims(1) = nvir
         voov_dims(2:3) = nocc
         voov_dims(4) = nvir


!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))

!        Construct necessary tensor 
         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(result_tensor,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,vovo_dims)

 !       Create identity matrix for the Occupied space
         ierr=talsh_tensor_construct(id_mat_occ,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,oo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
           h_tens(i,i) = ONE
         end do

!        Create identity matrix for the Virtual space
         ierr=talsh_tensor_construct(id_mat_vir,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(id_mat_vir,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,vv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
           h_tens(i,i) = ONE
         end do

         call get_CIS_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%print_level,&
                        perturbation=exa_input%do_cis_d)  !        Get molecular integral used in CIS 


!        Constrcut one body integral tensor
         ierr=talsh_tensor_construct(body_one_tensor,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(diag_foo_tensor,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(diag_fvv_tensor,C8,vv_dims,init_val=ZERO)

!        Get orbital energy
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)

         fo_dims(1:2) = nocc
         fv_dims(1:2) = nvir

         ierr=talsh_tensor_construct(diag_foo_tensor,C8,fo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(diag_foo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
             h_tens(i,i) = eps_occ(i)
         end do
     
         ierr=talsh_tensor_construct(diag_fvv_tensor,C8,fv_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(diag_fvv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
             h_tens(i,i) = eps_vir(i)
         end do

!        Set one body integral tensor element 
         ierr=talsh_tensor_contract("OB(b,j,a,i)+=F(i,j)*I(b,a)",body_one_tensor,diag_foo_tensor,id_mat_vir,scale=MINUS_ONE)
         ierr=talsh_tensor_contract("OB(b,j,a,i)+=F(b,a)*I(i,j)",body_one_tensor,diag_fvv_tensor,id_mat_occ)



!        Constrcut two body integral tensor and Create CIS Hamiltonian tensor 
         ierr=talsh_tensor_contract("HCIS(b,j,a,i)+=TB(j,b,i,a)*M()",result_tensor,int_t%ovov,one_tensor,scale=MINUS_ONE)
         ierr=talsh_tensor_add("HCIS(b,j,a,i)+=OB(b,j,a,i)",result_tensor,body_one_tensor)


        do j=1,nocc,1
            do b=1,nvir,1
                do i=1,nocc,1
                    do a=1,nvir,1
                        bj=(j-1)*nvir+b
                        ai=(i-1)*nvir+a
                        sigma1(ai,bj)=result_tens(b,j,a,i)
                    end do
                end do
            end do 
         end do 


        call diagonalize_cis_matrix(sigma1,eValues,eVectorsR) ! Get CIS eigenvector and excitation energy
!   Sort the eigenvector according to the excitation energy 
        call sort_cis_eig(eValues,eVectorsR)
!       call talsh_cis_pro_driver(exa_input,eValues,eVectorsR)

        n = size(eValues)
        call print_energy_CIS(eValues,size(eValues),n)

        if (present(CIS_vv_effden)) then
            call generate_cis_effden(exa_input,CIS_vv_effden,eVectorsR,eValues)
        end if

         ierr=talsh_tensor_destruct(id_mat_occ)
         ierr=talsh_tensor_destruct(id_mat_vir)
         ierr=talsh_tensor_destruct(result_tensor)

        if(exa_input%do_cis_d .eqv. .true.) then
            call print_date('Starting CIS(D) computation')
            call talsh_cis_d_driver(exa_input,int_t,eVectorsR,eValues,nocc,nvir,eps_occ,eps_vir)
        end if 

        cisrt = exa_input%do_cisrt
        if(cisrt) then
            call print_date('Starting Real-time CIS computation')
            call talsh_cis_rt_solver_driver(exa_input,eVectorsR,eValues)
        end if 

    end subroutine talsh_cis_driver 


    
!----------------------------------------------------------------------------------------------------------------------    
    subroutine diagonalize_cis_matrix(M,eVal,eVectorsR)

!       Set parameter needed in zgeev
        integer :: lwork, info,i,j
        real(kind=8), allocatable :: work(:)
        real(kind=8), allocatable :: rwork(:)
        integer :: vectors_L, vectors_R
        complex(kind=8), intent(inout), pointer :: M(:,:)
        complex(kind=8), pointer :: eVal(:) 
        complex(kind=8), pointer :: eVectorsR(:,:)
        complex(kind=8), pointer :: eVectorsL(:,:)
        character(len=1) :: do_left = 'N', do_right = 'V'


        vectors_R = size(M,1)
        vectors_L = size(M,2)
        lwork = 100000*vectors_R
        allocate(eVectorsL(vectors_L,vectors_R))
        allocate(eVectorsR(vectors_L,vectors_R))
        allocate(eVal(vectors_R))
        allocate(work(MAX(1,LWORK)))
        allocate(rwork(2*vectors_R))
        call zgeev (do_left, do_right, vectors_R, M, vectors_L, &
        eVal, eVectorsL, vectors_L, eVectorsR, vectors_R, &
        WORK, LWORK, rwork, INFO)


        
    end subroutine diagonalize_cis_matrix

!-----------------------------------------------------------------------------------------------------------------------
subroutine sort_cis_eig(eValues,eVectors)
    complex(kind=8), pointer :: eValues(:) 
    complex(kind=8), pointer :: eVectors(:,:)
    complex(kind=8), ALLOCATABLE :: eValues_tem(:) 
    complex(kind=8), ALLOCATABLE :: eValues_tem_2(:) 
    complex(kind=8), ALLOCATABLE :: eVectors_tem(:,:)
    integer, pointer :: index_cis_eig(:)
    real(8) :: tem, tem_internal
    integer :: i,j,index_tem,dege


    allocate(index_cis_eig(size(eValues)))
    allocate(eValues_tem(size(eValues)))
    allocate(eValues_tem_2(size(eValues)))
    allocate(eVectors_tem(size(eVectors,1),size(eVectors,2)))
    eValues_tem = eValues
    eVectors_tem = eVectors
    do i = 1, size(index_cis_eig),1
        index_cis_eig(i) = i
    end do

!        print *, "original index", index_cis_eig
    dege = 1
    do i=1,size(eValues)-1,1
        do j=i+1,size(eValues),1
        tem = abs(eValues(i))
        if (abs(eValues(i)) > abs(eValues(j))) then
            tem=eValues(i)
            index_tem = index_cis_eig(i)
            eValues(i)=eValues(j)
            index_cis_eig(i) = index_cis_eig(j)
            eValues(j)=tem
            index_cis_eig(j) =index_tem
        end if 
        end do 
        end do

    do i=1,size(index_cis_eig),1
        eVectors(:,i) = eVectors_tem(:,index_cis_eig(i))
    end do

    deallocate(eVectors_tem)
    deallocate(eValues_tem)
    deallocate(index_cis_eig)
end subroutine   

!-------------------------------------------------------------------------------------------------------------------------------

    subroutine get_CIS_integrals (nocc,nvir,mo_occ,mo_vir,int_t,print_level,perturbation)

!            Routine to get integrals needed in CIS in the form of antisymmetric tensors.
!            In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!            full integral list.
!            Written by Lucas Visscher, summer 2017
!            Revised by Xiang Yuan, June 2020
!            This routine now only contains the integrals int_t%oooo, int_t%vovo needed in CIS.
    
             use exacorr_global
             use talsh_ao_to_mo
    
             integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
             integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
             type(talsh_intg_tens_t),intent(inout)  :: int_t ! target integral tensors to be filled and given back to the caller
             integer, intent(in)                  :: print_level
    
             integer(INTD) :: oooo_dims(1:4),vovo_dims(1:4),ovov_dims(1:4),oovv_dims(1:4)
             integer(INTD) :: oovo_dims(1:4), ovvv_dims(1:4), voov_dims(1:4)

             integer(INTD)                      :: ierr
             complex(8), pointer, contiguous    :: ao_tens(:,:,:,:)
             type(C_PTR)                        :: body_p
    
             integer :: nao

             logical, optional :: perturbation
    
!            Auxilliary integral tensors (only needed inside this routine)
             type(talsh_tens_t) :: aoint_tensor, ovvo_tensor
             integer(INTD)      :: aoint_dims(1:4), ovvo_dims(1:4)
    
!            Arrays needed to communicate with the MO integral generator
             integer:: nmo(4)
             integer, allocatable :: mo_list(:)
    
!            Scalars (need to be defined as tensor types)
             type(talsh_tens_t) :: minusone_tensor
             integer(INTD)      :: minusone_dims(1)
    
             call print_date('entering integral transformation, scheme = 2')
    
!            Retrieve basis set information
             nao     = get_nao()      ! number of basis functions
    
!            Get tensor with AO integrals
             aoint_dims      = nao
             ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
             if (ierr /= 0) then
                print*, "Could not construct ao integral tensor, increasing TALSH_BUFF may help"
                call quit ('error constructing aoint_tensor')
             end if
             call print_date('allocated AO integral tensor')
!            Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
             ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
             call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
             call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
             call print_date('finished ao integral calculation')
        


!            Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.

!            Get ovov tensor
             ovov_dims(1) = nocc
             ovov_dims(2) = nvir
             ovov_dims(3) = nocc
             ovov_dims(4) = nvir
             ierr=talsh_tensor_construct(int_t%ovov,C8,ovov_dims,init_val=ZERO)

             allocate (mo_list(2*nvir+2*nocc))

             nmo(1:2) = nocc
             nmo(3:4) = nvir
             mo_list(1:nocc)                      = mo_occ(1:nocc)
             mo_list(nocc+1:2*nocc)               = mo_occ(1:nocc)
             mo_list(2*nocc+1:2*nocc+nvir)        = mo_vir(1:nvir)
             mo_list(2*nocc+nvir+1:2*nvir+2*nocc) = mo_vir(1:nvir)
             call get_integral_tensor (aoint_tensor,int_t%ovov,0,nmo,mo_list)
             deallocate (mo_list)

            !Get voov tensor
             voov_dims(1) = nvir
             voov_dims(2) = nocc
             voov_dims(3) = nocc
             voov_dims(4) = nvir
             ierr=talsh_tensor_construct(int_t%voov,C8,voov_dims,init_val=ZERO)
             allocate (mo_list(2*nvir+2*nocc))
             nmo(1) = nvir
             nmo(2) = nocc
             nmo(3) = nocc
             nmo(4) = nvir
             mo_list(1:nvir)                        = mo_vir(1:nvir) 
             mo_list(nvir+1:nocc+nvir)              = mo_occ(1:nocc)
             mo_list(nocc+nvir+1:2*nocc+nvir)       = mo_occ(1:nocc)
             mo_list(2*nocc+nvir+1:2*nvir+2*nocc)   = mo_vir(1:nvir)
             call get_integral_tensor (aoint_tensor,int_t%voov,0,nmo,mo_list)


             deallocate (mo_list)
!            Make anti-symmetrized vovo tensor
             minusone_dims(1) = 1
             ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)

             ierr=talsh_tensor_contract("V(i,a,j,b)+=V(a,j,i,b)*C()",int_t%ovov,int_t%voov,minusone_tensor)
             ierr=talsh_tensor_destruct(int_t%voov)
             call print_date('finished ovov integrals')
            

             if(present(perturbation)) then
                if(perturbation .eqv. .true.) then
                    !get oovo tensor
                    oovo_dims(1:2)   = nocc
                    oovo_dims(3) = nvir
                    oovo_dims(4) = nocc
                    ierr = talsh_tensor_construct(int_t%oovo, C8, oovo_dims, init_val=ZERO)
                    ierr = talsh_tensor_init(int_t%oovo, ZERO)
                    allocate (mo_list(3*nvir+nocc))
                    nmo(1)   = nocc
                    nmo(2)   = nvir
                    nmo(3)   = nocc
                    nmo(4)   = nocc
                    mo_list(1:nocc)                    = mo_occ(1:nocc)
                    mo_list(nocc+1:nvir+nocc)             = mo_vir(1:nvir)
                    mo_list(nvir+nocc+1:2*nocc+nvir)      = mo_occ(1:nocc)
                    mo_list(2*nocc+nvir+1:3*nocc+nvir) = mo_occ(1:nocc)
                    call get_integral_tensor (aoint_tensor,int_t%oovo,12,nmo,mo_list)
                    deallocate (mo_list)
                    call print_date('finished oovo integrals')


                    !Compute ovvv integrals
                    ovvv_dims(1)   = nocc
                    ovvv_dims(2:4) = nvir
                    ierr=talsh_tensor_construct(int_t%ovvv,C8,ovvv_dims,init_val=ZERO)
                    ierr = talsh_tensor_init(int_t%ovvv, ZERO)

                    allocate (mo_list(3*nvir+nocc))
                    nmo(1)   = nocc
                    nmo(2:4)   = nvir
                    mo_list(1:nocc)                    = mo_occ(1:nocc)
                    mo_list(nocc+1:nvir+nocc)          = mo_vir(1:nvir)
                    mo_list(nvir+nocc+1:2*nvir+nocc)   = mo_vir(1:nvir)
                    mo_list(2*nvir+nocc+1:3*nvir+nocc) = mo_vir(1:nvir)
                    call get_integral_tensor (aoint_tensor,int_t%ovvv,34,nmo,mo_list)
                    deallocate (mo_list)
                    call print_date('finished ovvv integrals')



                    !get MP2 integrals, code from talsh_mp2no.F90, routine: get_MP2_integrals
                    !get oovv integrals
                    oovv_dims(1:2) = nocc
                    oovv_dims(3:4) = nvir
                    ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
                    ierr = talsh_tensor_init(int_t%oovv, ZERO)
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
                endif
            endif

    
!            Clean-up
             ierr=talsh_tensor_destruct(aoint_tensor)
    
           end subroutine get_CIS_integrals

!-----------------------------------------------------------------------------------------------------------

    subroutine generate_cis_effden(exa_input,CIS_vv_effden,eVectorsR,eValues)
        type(exacc_input), intent(in) :: exa_input
        type(talsh_tens_t)   :: CIS_vv_effden
        complex(kind=8), pointer :: eVectorsR(:,:)
        complex(kind=8), pointer :: eValues(:)
        
        type(talsh_tens_t), pointer :: CIS_eigen_tensors_right(:) ! the CIS matrix is Hermitian so no need for left eigenvectors 
        integer            :: N_states, i_state
        integer            :: nvir, nocc
        integer            :: ierr
        integer(INTD), dimension(2) :: vv_dims 
        integer(INTD), dimension(2) :: cis_eig_dims
        complex(kind=8)    :: C_cocc

        nvir = exa_input%nvir
        nocc = exa_input%nocc
        cis_eig_dims(1) = nvir
        cis_eig_dims(2) = nocc
        vv_dims = nvir


        N_states = exa_input%neo_cis_nstates
        allocate(CIS_eigen_tensors_right(N_states))

        ierr=talsh_tensor_construct(CIS_vv_effden,C8,vv_dims,init_val=ZERO)

        do i_state=1, N_states
            ierr=talsh_tensor_construct(CIS_eigen_tensors_right(i_state),C8,cis_eig_dims,init_val=ZERO)
            call reshape_cis(eVectorsR,CIS_eigen_tensors_right(i_state),cis_eig_dims,i_state)
        end do

        do i_state=1, N_states
            C_cocc = ONE
            ierr=talsh_tensor_contract("D(a,b)+=CISL+(a,i)*CISR(b,i)",CIS_vv_effden,&
                CIS_eigen_tensors_right(i_state),CIS_eigen_tensors_right(i_state),scale=C_cocc)
        end do

        do i_state=1, N_states
            ierr=talsh_tensor_destruct(CIS_eigen_tensors_right(i_state))
        end do
        deallocate(CIS_eigen_tensors_right)

    end subroutine 


!-----------------------------------------------------------------------------------------------------------
    subroutine print_energy_CIS(eValues,size,n,pmin,perturbation)
        complex(kind=8) :: eValues(:)
        integer :: size, n, i  
        logical, optional :: perturbation
        integer, optional :: pmin
        if(present(perturbation)) then
            if (perturbation .eqv. .true.) then
                print *,'         -     CIS(D) excitation energy     -'
            else 
                print *,'         -     CIS excitation energy     -'
            endif
        else 
            print *,'         -     CIS excitation energy     -'
        endif
        
        print *,'-----------------------------------------------------------------------------------------'
        print *,'state | rel. eigenvalue (a.u.) | rel. eigenvalue (eV) | rel. eigenvalue (cm-1) |'
        if(present(pmin)) then
            do i=1,size,1
                print "(I5,(F18.5),(F23.5),(F27.5))",i+pmin-1,abs(eValues(i)),abs(eValues(i))*27.211386245988,abs(eValues(i))*219474.63137
            end do
        else
            if (n > size) then
                print *,'wrong number of state, please revise n'
            end if 
            do i=1,size,1
                print "(I5,(F18.5),(F23.5),(F27.5))",i,abs(eValues(i)),abs(eValues(i))*27.211386245988,abs(eValues(i))*219474.63137
            end do
        end if
        print *,'-----------------------------------------------------------------------------------------'

    end subroutine print_energy_CIS

!---------------------------------------------------------------------------------------------------------
    subroutine talsh_cis_d_driver(exa_input,int_t,eVectorsR,eValues,nocc,nvir,eps_occ,eps_vir)

        type(exacc_input), intent(in) :: exa_input
        !        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
        type(talsh_intg_tens_t), intent(inout) :: int_t
        complex(kind=8), pointer, intent(inout) :: eVectorsR(:,:)
        complex(kind=8), pointer, intent(inout) :: eValues(:)
        integer,intent(inout)                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
        !        orbital energies
        real(8), allocatable,intent(inout)     :: eps_occ(:),eps_vir(:) 

        
        integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
        integer(C_INT)       :: ierr
        complex(8), pointer :: result_tens(:,:,:,:),sigma1(:,:),check_tens(:,:)
        integer(INTD)       :: fo_dims(1:2),fv_dims(1:2)
        
       
        type(C_PTR):: body_p
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max
        complex(8), pointer, contiguous :: h_tens(:,:)
        integer :: i,j,a,b,n,k,p
        integer :: bj,ai



!        single-body and body_two tensors
        type(talsh_tens_t) :: body_one_tensor, diag_foo_tensor, diag_fvv_tensor
        type(talsh_tens_t) :: id_mat_occ, id_mat_vir
!        intermediate tensors (same notation as in relccsd)
        integer(INTD), dimension(2) :: oo_dims, vv_dims
        integer(INTD), dimension(4) :: oooo_dims, vovo_dims, voov_dims
!        scalars (need to be defined as tensor types)
        type(talsh_tens_t) :: result_tensor,sigma_tensor
        integer(C_INT)     :: one_dims(1)

!        CIS control and result variables
        real(8) :: tem
        integer :: ncycles
        integer :: tar_state

!        CIS(D) variables
       complex(8), pointer, contiguous :: s_tens(:)
       complex(8), pointer, contiguous :: t_tens(:,:)
       complex(8), pointer, contiguous :: q_tens(:,:,:,:)

       integer(INTD)                   :: oovv_dims(1:4)
       integer(INTD)                   :: vvoo_dims(1:4)
       integer(INTD)                   :: ov_dims(1:2)
       integer                         :: pmin, pmax

       type(talsh_tens_t)              :: coefs_tensor
       type(talsh_tens_t)              :: t2_tensor
       type(talsh_tens_t)              :: nu_tensor
       type(talsh_tens_t)              :: temp_tensor1
       type(talsh_tens_t)              :: temp_tensor2
       type(talsh_tens_t)              :: temp_tensor3
       type(talsh_tens_t)              :: u_tensor
       type(talsh_tens_t)              :: e_tensor
       type(talsh_tens_t)              :: b_tensor
       type(talsh_tens_t)              :: eps_diag_tensor
       complex(8), dimension(:),  allocatable :: e_vector

        pmin = exa_input%cis_d_pmin
        pmax = exa_input%cis_d_pmax

        !Defining dimensions arrays
        oovv_dims(1:2) = nocc
        oovv_dims(3:4) = nvir
        ov_dims(:)     = [nocc, nvir]
        vvoo_dims(1:2) = nvir 
        vvoo_dims(3:4) = nocc
        vv_dims(:)     = [nvir, nvir]
        oo_dims(:)     = [nocc, nocc]

        allocate(e_vector(pmax-pmin+1))
!           Initializing tensors before the loop
        ierr = talsh_tensor_construct(t2_tensor, C8, vvoo_dims) 
        ierr = talsh_tensor_construct(coefs_tensor, C8, ov_dims)
        ierr = talsh_tensor_construct(nu_tensor, C8, ov_dims)
        ierr = talsh_tensor_construct(e_tensor, C8, one_dims(1:0))
        ierr = talsh_tensor_construct(b_tensor, C8, oovv_dims)
        ierr = talsh_tensor_construct(u_tensor, C8, oovv_dims)
        ierr = talsh_tensor_construct(temp_tensor1, C8, oo_dims)
        ierr = talsh_tensor_construct(temp_tensor2, C8, vv_dims)
        ierr = talsh_tensor_construct(temp_tensor3, C8, ov_dims)


        ierr = talsh_tensor_init(t2_tensor, ZERO)
        ierr = talsh_tensor_init(temp_tensor1, ZERO)
        
        e_vector(:) = ZERO
        

!            MP2 amplitudes
        ierr = talsh_tensor_contract("T(a,b,i,j)+=I+(i,j,a,b)*C()", t2_tensor, int_t%oovv, one_tensor, SCALE=MINUS_ONE)
        call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=t2_tensor)

!           Initializing temporary tensor that will be use in tensor contraction
        
        ierr = talsh_tensor_contract("TMP(j,i)+=I(j,k,b,c)*T2(b,c,i,k)", temp_tensor1, int_t%oovv, t2_tensor)
        
        
!           Loop on needed excitations level
        do p = pmin, pmax
            ierr = talsh_tensor_init(coefs_tensor, ZERO)
            ierr = talsh_tensor_init(nu_tensor, ZERO)
            ierr = talsh_tensor_init(e_tensor, ZERO)
            ierr = talsh_tensor_init(u_tensor, ZERO)
            ierr = talsh_tensor_init(b_tensor, ZERO)
            ierr = talsh_tensor_init(temp_tensor2, ZERO)
            ierr = talsh_tensor_init(temp_tensor3, ZERO)

            ierr = talsh_tensor_get_body_access(coefs_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,t_tens,ov_dims) ! to use as a regular Fortran array
            t_tens = ZERO
            j = 1
            do i=1, nocc, 1
                t_tens(i,:) = eVectorsR(j:j+nvir-1,p)
                j = j+nvir
            end do     

            ierr = talsh_tensor_contract("TMP(a,b)+=I(j,k,b,c)*T2(a,c,j,k)", temp_tensor2, int_t%oovv, t2_tensor)
            ierr = talsh_tensor_contract("NU(i,a)+=TMP+(a,b)*C(i,b)", nu_tensor, temp_tensor2, coefs_tensor, scale=ONE_HALF)

            ierr = talsh_tensor_contract("NU(i,a)+=TMP+(j,i)*C(j,a)", nu_tensor, temp_tensor1, coefs_tensor, scale=ONE_HALF)

            ierr = talsh_tensor_contract("TMP(k,c)+=I(j,k,b,c)*C(j,b)", temp_tensor3, int_t%oovv, coefs_tensor)
            ierr = talsh_tensor_contract("NU(i,a)+=TMP(k,c)*T2+(c,a,i,k)", nu_tensor, temp_tensor3, t2_tensor)
            
            ierr = talsh_tensor_contract("E()+=C+(i,a)*NU(i,a)", e_tensor, coefs_tensor, nu_tensor)

            ierr = talsh_tensor_contract("U(i,j,a,b)+=I(i,c,a,b)*B(j,c)", u_tensor, int_t%ovvv, coefs_tensor)
            ierr = talsh_tensor_contract("U(i,j,a,b)+=I(j,c,a,b)*B(i,c)", u_tensor, int_t%ovvv, coefs_tensor, scale=MINUS_ONE)
            ierr = talsh_tensor_contract("U(i,j,a,b)+=I(i,j,a,k)*B(k,b)", u_tensor, int_t%oovo, coefs_tensor, scale=MINUS_ONE)
            ierr = talsh_tensor_contract("U(i,j,a,b)+=I(i,j,b,k)*B(k,a)", u_tensor, int_t%oovo, coefs_tensor, scale=ONE)

            ierr = talsh_tensor_copy("B(i,j,a,b)=U(i,j,a,b)", b_tensor, u_tensor)

            !Scale with denominator:
            ierr = talsh_tensor_get_body_access(b_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,q_tens,oovv_dims) ! to use as a regular Fortran array
            do i = 1, nocc
                do j = 1, nocc
                    do a = 1, nvir
                        do b = 1, nvir
                            q_tens(i,j,a,b) = -q_tens(i,j,a,b)/(eps_vir(a) + eps_vir(b) - eps_occ(i) - eps_occ(j) - eValues(p))
                        end do
                    end do
                end do
            end do

            ierr = talsh_tensor_contract("E()+=B+(i,j,a,b)*U(i,j,a,b)", e_tensor, b_tensor, u_tensor, scale=ONE_QUARTER_C)
            ierr = talsh_tensor_get_body_access(e_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,s_tens,one_dims) ! to use as a regular Fortran array
            e_vector(p-pmin+1) = s_tens(1) + eValues(p)
        end do

!       clean up tensors
        ierr = talsh_tensor_destruct(coefs_tensor)
        ierr = talsh_tensor_destruct(t2_tensor)
        ierr = talsh_tensor_destruct(nu_tensor)
        ierr = talsh_tensor_destruct(u_tensor)
        ierr = talsh_tensor_destruct(e_tensor)
        ierr = talsh_tensor_destruct(b_tensor)
        call print_energy_CIS(e_vector,size(e_vector),n, pmin, .true.)
    end subroutine

end module
