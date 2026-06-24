module eti_cis
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy

    use tensor_algebra
    use exacorr_tensor_interface
    use exacorr_datatypes
    use exacorr_cfg, only: exacc_input
    use exacorr_utils
    use exacorr_scaling_parameters
    use eti_common_routines
    use eti_cis_property
    use eti_cis_density
    use eti_rt_solver, only: eti_cis_rt_solver_driver
    use eti_ao_to_mo

    implicit none

    type(tensor_t) :: one_tensor
    private

    public eti_cis_driver

    contains
!---------------------------------------------------------------------------------------------------------------

    subroutine eti_cis_driver (exa_input,CIS_vv_effden,initialize_tensor_library,initialize_interest)
    ! This routine drives CIS calculation including construction of CIS Hamiltonian and diagnolization
    ! Written by Xiang Yuan, June 2020, extended since by others

        use eti_integrals

        type(exacc_input), intent(in) :: exa_input
        type(tensor_t), optional   :: CIS_vv_effden
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
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer :: i,j,a,b,n
         integer :: bj,ai
         logical :: cisrt


!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_intg_tens_t) :: int_t
!        single-body and body_two tensors
         type(tensor_t) :: body_one_tensor, diag_foo_tensor, diag_fvv_tensor
         type(tensor_t) :: id_mat_occ, id_mat_vir
!        intermediate tensors (same notation as in relccsd)
         integer(INTD), dimension(2) :: oo_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims, vovo_dims, voov_dims
!        scalars (need to be defined as tensor types)
         type(tensor_t) :: result_tensor,sigma_tensor
         integer(C_INT)     :: one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CIS control and result variables
         real(8) :: tem
         integer :: ncycles
         integer :: tar_state
         logical :: initialize_tensor_library_local
         logical :: initialize_interest_local

         call print_date('Entered cis_driver routine')
!        The user may have modified the default executor via input
         if (exa_input%executor /= 0) then ! 0 means no input value is set
            ierr = set_default_executor(exa_input%executor)
            if (ierr /= 0) write (*,*) ' WARNING: incompatible executor specified in input, keeping default'
            select case(get_default_executor())
            case(1)
                write(*,*) ' User set TAPP (executor 1) via input as the default for tensor operations'
            case(2)
                write(*,*) ' User set TALSH (executor 2) via input as the default for tensor operations'
            case(3)
                write(*,*) ' User set ExaTensor (executor 3) via input as the default for tensor operations'
            end select
         end if

!        Initialize libraries
         call eti_executor_initialize(exa_input%talsh_buff,get_default_executor())
         call interest_initialize(.true.)
         call print_date('Initialized interest library')
        
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
         ierr=construct_tensor_t(one_tensor,exa_input%eti_precision%storage_t,one_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(result_tensor,exa_input%eti_precision%storage_t,voov_dims,init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,voov_dims)

 !       Create identity matrix for the Occupied space
         ierr=construct_tensor_t(id_mat_occ,exa_input%eti_precision%storage_t,oo_dims,init_val=ZERO)
         ierr=direct_data_pointer(id_mat_occ,body_p)
         call c_f_pointer(body_p,h_tens,oo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
           h_tens(i,i) = ONE
         end do

!        Create identity matrix for the Virtual space
         ierr=construct_tensor_t(id_mat_vir,exa_input%eti_precision%storage_t,vv_dims,init_val=ZERO)
         ierr=direct_data_pointer(id_mat_vir,body_p)
         call c_f_pointer(body_p,h_tens,vv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
           h_tens(i,i) = ONE
         end do

         ! Get molecular integrals used in CIS
         int_t%vovo_active = .true.
         if (exa_input%do_cis_d) then 
            int_t%vovv_active = .true.
            int_t%ooov_active = .true.
            int_t%oovv_active = .true.
         else
            int_t%vovv_active = .false.
            int_t%ooov_active = .false.
            int_t%oovv_active = .false.
         end if
         int_t%oooo_active = .false.
         int_t%vvvv_active = .false.
         
         call get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t,1,exa_input%print_level,exa_input%eti_precision)

!        Get orbital energy
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)

         fo_dims(1:2) = nocc
         fv_dims(1:2) = nvir

         ierr=construct_tensor_t(diag_foo_tensor,exa_input%eti_precision%storage_t,fo_dims,init_val=ZERO)
         ierr=direct_data_pointer(diag_foo_tensor,body_p)
         call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
             h_tens(i,i) = eps_occ(i)
         end do
     
         ierr=construct_tensor_t(diag_fvv_tensor,exa_input%eti_precision%storage_t,fv_dims,init_val=ZERO)
         ierr=direct_data_pointer(diag_fvv_tensor,body_p)
         call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
             h_tens(i,i) = eps_vir(i)
         end do

!        Construct one body integral tensor
         ierr=construct_tensor_t(body_one_tensor,exa_input%eti_precision%storage_t,voov_dims,init_val=ZERO)

!        Set one body integral tensor element 
         ierr=t_einsum("OB(b,i,j,a)+=F(i,j)*I(b,a)",body_one_tensor,diag_foo_tensor,id_mat_vir,scale=MINUS_ONE)
         ierr=t_einsum("OB(b,i,j,a)+=F(b,a)*I(i,j)",body_one_tensor,diag_fvv_tensor,id_mat_occ)
         ierr=destroy_tensor_t(id_mat_occ)
         ierr=destroy_tensor_t(id_mat_vir)
         ierr=destroy_tensor_t(diag_foo_tensor)
         ierr=destroy_tensor_t(diag_fvv_tensor)

!        Construct two body integral tensor and Create CIS Hamiltonian tensor 
         ierr=t_einsum("HCIS(b,i,j,a)+=TB(b,i,a,j)*M()",result_tensor,int_t%vovo,one_tensor,scale=MINUS_ONE)
         ierr=t_einsum("HCIS(b,i,j,a)+=OB(b,i,j,a)",result_tensor,body_one_tensor)

!       Copy data into a regular fortran array so that we can diagonalize this.
        do j=1,nocc,1
            do b=1,nvir,1
                do i=1,nocc,1
                    do a=1,nvir,1
                        bj=(j-1)*nvir+b
                        ai=(i-1)*nvir+a
                        sigma1(ai,bj)=result_tens(b,i,j,a)
                    end do
                end do
            end do 
         end do 

        ierr=destroy_tensor_t(result_tensor)
        ierr=destroy_tensor_t(body_one_tensor)

        call diagonalize_cis_matrix(sigma1,eValues,eVectorsR) ! Get CIS eigenvector and excitation energy
!   Sort the eigenvector according to the excitation energy 
        call sort_cis_eig(eValues,eVectorsR)
!       call eti_cis_pro_driver(exa_input,eValues,eVectorsR)

        n = size(eValues)
        call print_energy_CIS(eValues,size(eValues),n)

        if (present(CIS_vv_effden)) then
            call generate_cis_effden(exa_input,CIS_vv_effden,eVectorsR,eValues)
        end if

        if(exa_input%do_cis_d .eqv. .true.) then
            call print_date('Starting CIS(D) computation')
            call eti_cis_d_driver(exa_input,int_t,eVectorsR,eValues,nocc,nvir,eps_occ,eps_vir)
        end if 

        cisrt = exa_input%do_cisrt
        if(cisrt) then
            call print_date('Starting Real-time CIS computation')
            call eti_cis_rt_solver_driver(exa_input,eVectorsR,eValues)
        end if 

    end subroutine eti_cis_driver 
    
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

    subroutine generate_cis_effden(exa_input,CIS_vv_effden,eVectorsR,eValues)
        type(exacc_input), intent(in) :: exa_input
        type(tensor_t)   :: CIS_vv_effden
        complex(kind=8), pointer :: eVectorsR(:,:)
        complex(kind=8), pointer :: eValues(:)
        
        type(tensor_t), pointer :: CIS_eigen_tensors_right(:) ! the CIS matrix is Hermitian so no need for left eigenvectors 
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

        ierr=construct_tensor_t(CIS_vv_effden,exa_input%eti_precision%storage_t,vv_dims,init_val=ZERO)

        do i_state=1, N_states
            ierr=construct_tensor_t(CIS_eigen_tensors_right(i_state),exa_input%eti_precision%storage_t,&
                                  cis_eig_dims,init_val=ZERO)
            call reshape_cis(eVectorsR,CIS_eigen_tensors_right(i_state),cis_eig_dims,i_state)
        end do

        do i_state=1, N_states
            C_cocc = ONE
            ierr=t_einsum("D(a,b)+=CISL+(a,i)*CISR(b,i)",CIS_vv_effden,&
                CIS_eigen_tensors_right(i_state),CIS_eigen_tensors_right(i_state),scale=C_cocc)
        end do

        do i_state=1, N_states
            ierr=destroy_tensor_t(CIS_eigen_tensors_right(i_state))
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
        print *,'state | rel. eigenvalue (a.u.) | rel. eigenvalue (eV) | rel. eigenvalue (cm-1)'
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
    subroutine eti_cis_d_driver(exa_input,int_t,eVectorsR,eValues,nocc,nvir,eps_occ,eps_vir)

        type(exacc_input), intent(in) :: exa_input
        !        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
        type(eti_intg_tens_t), intent(inout) :: int_t
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
        complex(8), pointer, contiguous :: h_tens(:,:)
        integer :: i,j,a,b,n,k,p
        integer :: bj,ai



!        single-body and body_two tensors
        type(tensor_t) :: body_one_tensor, diag_foo_tensor, diag_fvv_tensor
        type(tensor_t) :: id_mat_occ, id_mat_vir
!        intermediate tensors (same notation as in relccsd)
        integer(INTD), dimension(2) :: oo_dims, vv_dims
        integer(INTD), dimension(4) :: oooo_dims, vovo_dims, voov_dims
!        scalars (need to be defined as tensor types)
        type(tensor_t) :: result_tensor,sigma_tensor
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

       type(tensor_t)              :: coefs_tensor
       type(tensor_t)              :: t2_tensor
       type(tensor_t)              :: nu_tensor
       type(tensor_t)              :: temp_tensor1
       type(tensor_t)              :: temp_tensor2
       type(tensor_t)              :: temp_tensor3
       type(tensor_t)              :: u_tensor
       type(tensor_t)              :: e_tensor
       type(tensor_t)              :: b_tensor
       type(tensor_t)              :: eps_diag_tensor
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
        ierr = construct_tensor_t(t2_tensor, C8, vvoo_dims) 
        ierr = construct_tensor_t(coefs_tensor, C8, ov_dims)
        ierr = construct_tensor_t(nu_tensor, C8, ov_dims)
        ierr = construct_tensor_t(e_tensor, C8, one_dims(1:0))
        ierr = construct_tensor_t(b_tensor, C8, oovv_dims)
        ierr = construct_tensor_t(u_tensor, C8, oovv_dims)
        ierr = construct_tensor_t(temp_tensor1, C8, oo_dims)
        ierr = construct_tensor_t(temp_tensor2, C8, vv_dims)
        ierr = construct_tensor_t(temp_tensor3, C8, ov_dims)


        ierr = initialize_tensor_t(t2_tensor, ZERO)
        ierr = initialize_tensor_t(temp_tensor1, ZERO)
        
        e_vector(:) = ZERO
        

!            MP2 amplitudes
        ierr = t_einsum("T(a,b,i,j)+=I+(i,j,a,b)*C()", t2_tensor, int_t%oovv, one_tensor, SCALE=MINUS_ONE)
        call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=t2_tensor)

!           Initializing temporary tensor that will be use in tensor contraction
        
        ierr = t_einsum("TMP(j,i)+=I(j,k,b,c)*T2(b,c,i,k)", temp_tensor1, int_t%oovv, t2_tensor)
        
        
!           Loop on needed excitations level
        do p = pmin, pmax
            ierr = initialize_tensor_t(coefs_tensor, ZERO)
            ierr = initialize_tensor_t(nu_tensor, ZERO)
            ierr = initialize_tensor_t(e_tensor, ZERO)
            ierr = initialize_tensor_t(u_tensor, ZERO)
            ierr = initialize_tensor_t(b_tensor, ZERO)
            ierr = initialize_tensor_t(temp_tensor2, ZERO)
            ierr = initialize_tensor_t(temp_tensor3, ZERO)

            ierr = direct_data_pointer(coefs_tensor,body_p)
            call c_f_pointer(body_p,t_tens,ov_dims) ! to use as a regular Fortran array
            t_tens = ZERO
            j = 1
            do i=1, nocc, 1
                t_tens(i,:) = eVectorsR(j:j+nvir-1,p)
                j = j+nvir
            end do     

            ierr = t_einsum("TMP(a,b)+=I(j,k,b,c)*T2(a,c,j,k)", temp_tensor2, int_t%oovv, t2_tensor)
            ierr = t_einsum("NU(i,a)+=TMP+(a,b)*C(i,b)", nu_tensor, temp_tensor2, coefs_tensor, scale=ONE_HALF)

            ierr = t_einsum("NU(i,a)+=TMP+(j,i)*C(j,a)", nu_tensor, temp_tensor1, coefs_tensor, scale=ONE_HALF)

            ierr = t_einsum("TMP(k,c)+=I(j,k,b,c)*C(j,b)", temp_tensor3, int_t%oovv, coefs_tensor)
            ierr = t_einsum("NU(i,a)+=TMP(k,c)*T2+(c,a,i,k)", nu_tensor, temp_tensor3, t2_tensor)
            
            ierr = t_einsum("E()+=C+(i,a)*NU(i,a)", e_tensor, coefs_tensor, nu_tensor)

            ierr = t_einsum("U(i,j,a,b)+=I(c,j,a,b)*B(i,c)", u_tensor, int_t%vovv, coefs_tensor)
            ierr = t_einsum("U(i,j,a,b)+=I(c,i,a,b)*B(j,c)", u_tensor, int_t%vovv, coefs_tensor, scale=MINUS_ONE)
            ierr = t_einsum("U(i,j,a,b)+=I(i,j,k,a)*B(k,b)", u_tensor, int_t%ooov, coefs_tensor, scale=ONE)
            ierr = t_einsum("U(i,j,a,b)+=I(i,j,k,b)*B(k,a)", u_tensor, int_t%ooov, coefs_tensor, scale=MINUS_ONE)

            ierr = t_einsum("B(i,j,a,b)=U(i,j,a,b)", b_tensor, u_tensor)

            !Scale with denominator:
            ierr = direct_data_pointer(b_tensor,body_p)
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

            ierr = t_einsum("E()+=B+(i,j,a,b)*U(i,j,a,b)", e_tensor, b_tensor, u_tensor, scale=ONE_QUARTER)
            ierr = direct_data_pointer(e_tensor,body_p)
            call c_f_pointer(body_p,s_tens,one_dims) ! to use as a regular Fortran array
            e_vector(p-pmin+1) = s_tens(1) + eValues(p)
        end do

!       clean up tensors
        ierr = destroy_tensor_t(coefs_tensor)
        ierr = destroy_tensor_t(t2_tensor)
        ierr = destroy_tensor_t(nu_tensor)
        ierr = destroy_tensor_t(u_tensor)
        ierr = destroy_tensor_t(e_tensor)
        ierr = destroy_tensor_t(b_tensor)
        call print_energy_CIS(e_vector,size(e_vector),n, pmin, .true.)
    end subroutine

end module
