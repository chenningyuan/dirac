module talsh_cis
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use intermediates
    use talsh_common_routines
!    use davidson
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

    subroutine talsh_cis_driver (exa_input)
    ! This routine drives CIS calculation including construction of CIS Hamiltonian and diagnolization
    !  Written by Xiang Yuan, June 2020
        use exacorr_global
        use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
        use talsh_gradient


        type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer(C_INT)       :: ierr
         complex(8), pointer :: result_tens(:,:,:,:),sigma1(:,:),check_tens(:,:)
         complex(kind=8), pointer :: eValues(:) 
         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer :: i,j,a,b
         integer :: bj,ai


!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens) :: int_t
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
         real(8) :: t_target_precision
         integer :: ncycles

!        Initialize libraries
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)
         call print_date('Initialized talsh library')
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max
         call interest_initialize(.true.)
         call print_date('Initialized interest library')

        !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         allocate(sigma1(nocc*nvir,nocc*nvir))
         allocate(eValues(nocc*nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir


         !Initialize scalars that are to be used as tensors in contractions
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


         call print_date("Getting orbital energies")
         ! Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("Retrieved orbital energies")


 !      create identity matrix for the occupied space
         ierr=talsh_tensor_construct(id_mat_occ,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,oo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
           h_tens(i,i) = ONE
         end do
!        create identity matrix for the virtual space
         ierr=talsh_tensor_construct(id_mat_vir,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(id_mat_vir,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,vv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
           h_tens(i,i) = ONE
         end do

         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(result_tensor,C8,voov_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,voov_dims)

         print *,'before get_CIS_integrals'
         call get_CIS_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%print_level)
         print *,'after get_CIS_integrals'

         !Constrcut body_one tensor
         ierr=talsh_tensor_construct(body_one_tensor,C8,voov_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(diag_foo_tensor,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(diag_fvv_tensor,C8,vv_dims,init_val=ZERO)
         !compute fock matrix tensor 
         print *,'before get_CIS_fockmatrix'
         call get_CIS_fockmatrix (nocc,nvir,mo_occ,mo_vir,diag_foo_tensor,diag_fvv_tensor, &
                                  int_t, exa_input%print_level)


         print *,'after get_CIS_fockmatrix'
!        set body_one tensor element 
         ierr=talsh_tensor_contract("OB(b,i,j,a)+=F(i,j)*I(b,a)",body_one_tensor,diag_foo_tensor,id_mat_vir,scale=MINUS_ONE)
         print *,'OB(b,i,j,a)+=F(i,j)*I(b,a) ierr',ierr
         ierr=talsh_tensor_contract("OB(b,i,j,a)+=F(b,a)*I(i,j)",body_one_tensor,diag_fvv_tensor,id_mat_occ)
         print *,'OB(b,i,j,a)+=F(b,a)*I(i,j) ierr',ierr

         !ierr=talsh_tensor_get_body_access(body_one_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         !call c_f_pointer(body_p,result_tens,voov_dims)
         !print *,'int_t%vovo(2,1,2,1)',result_tens(2,1,2,1)
         !print *,'OB(b,i,j,a)',result_tens(1,1,1,1)



          !create CIS hamiltonian tensor 
         ierr=talsh_tensor_contract("HCIS(b,i,j,a)+=TB(b,i,a,j)*M()",result_tensor,int_t%vovo,one_tensor,scale=MINUS_ONE)
         print *,'HCIS(b,i,j,a)+=TB(b,i,a,j)*M() ierr',ierr
         !print *,'result(1,1,1,1)',result_tens(2,1,1,1)
         !print *,'result(2,1,1,2)',result_tens(2,1,1,2)
         !do i=1,nvir,1
            !print *,'two-electron',result_tens(i,1,1,i)
         !end do
         ierr=talsh_tensor_add("HCIS(b,i,j,a)+=OB(b,i,j,a)",result_tensor,body_one_tensor)
         print *,'HCIS(b,i,j,a)+=OB(b,i,j,a) ierr',ierr
         print *,'nocc',nocc
         print *,'nvir',nvir
         print *,'result_tens',size(result_tens)

         !print *,'result(1,1,1,1)',result_tens(1,1,1,1)
         !print *,'result(2,1,1,2)',result_tens(2,1,1,2)


         ! transform tensor to matrix sigma1 and digonal it

        do j=1,nocc,1
            do b=1,nvir,1
                do i=1,nocc,1
                    do a=1,nvir,1
                        bj=(j-1)*nvir+b
                        ai=(i-1)*nvir+a
                        !if (ABS(result_tens(b,i,j,a))>0.001) then
                        sigma1(bj,ai)=result_tens(b,i,j,a)
                            !sigma1(bj,ai)=(1.0000000000000000,0.0000000000000000)
                        !else 
                            !sigma1(bj,ai)=(0.0000000000000000,0.0000000000000000)
                        !end if 
                    end do
                end do
            end do 
         end do 

         !do i=1,nocc*nvir,1
            !print *,'diagonal element_sigma1',sigma1(1,i)
         !end do



         call diagoanl_cis_matrix(sigma1,eValues)
         print *,'after diagoanl_cis_matrix',eValues


         ierr=talsh_tensor_destruct(id_mat_occ)
         ierr=talsh_tensor_destruct(id_mat_vir)
         ierr=talsh_tensor_destruct(result_tensor)
    end subroutine talsh_cis_driver 


    
!----------------------------------------------------------------------------------------------------------------------    
    subroutine diagoanl_cis_matrix(M,eVal)

    ! set parameter needed in lapack
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
        print *,'vectors_R',vectors_R
        print *,'vectors_L',vectors_L
        !print *,'sigma(3,4)_dia',M(3,4)
        !print *,'sigma(4,3)_dia',M(4,3)
        lwork = 10*vectors_R	
        allocate(eVectorsL(vectors_L,vectors_R))
        allocate(eVectorsR(vectors_L,vectors_R))
        allocate(eVal(vectors_R))
        !eVectorsR = 0.0d0
        allocate(work(MAX(1,LWORK)))
        allocate(rwork(2*vectors_R))
        !print *,'Matrix',M
        call zgeev (do_left, do_right, vectors_R, M, vectors_L, &
        eVal, eVectorsL, vectors_L, eVectorsR, vectors_R, &
        WORK, LWORK, rwork, INFO)
        print *,'INFO',info


    end subroutine diagoanl_cis_matrix

!-----------------------------------------------------------------------------------------------------------------------

    subroutine get_CIS_fockmatrix (nocc,nvir,mo_occ,mo_vir,diag_mat_occ, diag_mat_vir, &
        int_t, print_level)

!        Written by Lucas Visscher, summer 2017
!        Revised by Xiang Yuan, June 2020
    use exacorr_global

    implicit none

    integer, intent(in)                  :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
    integer, intent(in)                  :: print_level
    integer, intent(in)                  :: mo_occ(:),mo_vir(:) ! the orbital energies


!        target tensors to be filled and given back to the caller
    type(talsh_tens_t), intent(inout)    :: diag_mat_occ, diag_mat_vir
    type(talsh_intg_tens), intent(inout) :: int_t

    type(talsh_tens_t)              :: fo_tensor, fv_tensor
    type(talsh_tens_t)              :: id_mat_occ 
    integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2)
    integer(INTD)                   :: ierr
    logical                         :: one_el_exist
    real(8), allocatable            :: eps_occ(:),eps_vir(:)
    real(8)                         :: e_core
    integer                         :: min_occ
    type(C_PTR)                     :: body_p
    complex(8), pointer, contiguous :: h_tens(:,:)
    integer                         :: i, j, iv, jv, imo, jmo


    min_occ=minval(mo_occ)
    allocate (eps_occ(nocc))
    allocate (eps_vir(nvir))

!        Get foo tensor
    fo_dims(1:2) = nocc
    ierr=talsh_tensor_construct(fo_tensor,C8,fo_dims,init_val=ZERO)

!        Get fvv tensor
    fv_dims(1:2) = nvir
    ierr=talsh_tensor_construct(fv_tensor,C8,fv_dims,init_val=ZERO)

    !   get one electron integrals
    one_el_exist = exist_one_el(e_core)

!           create density (identity) matrix for the occupied space
    ierr=talsh_tensor_construct(id_mat_occ,C8,fo_dims,init_val=ZERO)
    ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
    h_tens = ZERO
    do i = 1, nocc
    h_tens(i,i) = ONE
    end do

!           create fock matrix for the occupied and virtual space
    ierr=talsh_tensor_get_body_access(fo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
    call get_one_el(h_tens,mo_occ,nocc,mo_occ,nocc,min_occ)

    ierr=talsh_tensor_get_body_access(fv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
    call get_one_el(h_tens,mo_vir,nvir,mo_vir,nvir,min_occ)
    call print_date('one electron integrals added to Fock Matrix')

    ierr=talsh_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",fo_tensor,int_t%oooo,id_mat_occ, ONE)
    ierr=talsh_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",fv_tensor,int_t%vovo,id_mat_occ, ONE)
    call print_date('two electron integrals added to Fock Matrix')

    if (print_level >10) then
        call print_tensor(fo_tensor, 1.0D-6, 'fo_tensor')
        call print_tensor(fv_tensor, 1.0D-6, 'fv_tensor')
    end if


!set diag element of fock matrix 
    ierr=talsh_tensor_get_body_access(fo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
    do i=1,fo_dims(1)
        eps_occ(i)=h_tens(i,i)
        print *,'eps_occ(i)',eps_occ(i)
    end do

    ierr=talsh_tensor_get_body_access(fv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
    do i=1,fv_dims(1)
        eps_vir(i)=h_tens(i,i)
        print *,'eps_vir(i)',eps_vir(i)
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
    if (print_level >0) call print_date('Fock matrix initialized')  
    end subroutine get_CIS_fockmatrix

!-----------------------------------------------------------------------------------------------------------------


    subroutine get_CIS_integrals (nocc,nvir,mo_occ,mo_vir,int_t,print_level)

    !        Routine to get integrals needed in CC in the form of antisymmetric tensors.
    !        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
    !        full integral list.
    
    !        Written by Lucas Visscher, summer 2017
    !        Revised by Xiang Yuan, June 2020
    !        This routine now only contains the integrals int_t%oooo, int_t%vovo needed in CIS.
    
             use exacorr_global
             use talsh_ao_to_mo
    
             integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
             integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
    !        target integral tensors to be filled and given back to the caller
             type(talsh_intg_tens),intent(inout)  :: int_t
             integer, intent(in)                  :: print_level
    
             integer(INTD) :: oooo_dims(1:4),vovo_dims(1:4)
    
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
    
    !        Get antisymmetric oooo tensor
             oooo_dims = nocc
             ierr=talsh_tensor_construct(int_t%oooo,C8,oooo_dims,init_val=ZERO)
             allocate (mo_list(4*nocc))
             nmo = nocc
             mo_list(1:nocc)          = mo_occ(1:nocc)
             mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
             mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
             mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
             call get_integral_tensor (aoint_tensor,int_t%oooo,12,nmo,mo_list)
             deallocate (mo_list)
             call print_date('finished oooo integrals')
    
    !        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
    !        Get vovo tensor
             vovo_dims(1) = nvir
             vovo_dims(2) = nocc
             vovo_dims(3) = nvir
             vovo_dims(4) = nocc
             ierr=talsh_tensor_construct(int_t%vovo,C8,vovo_dims,init_val=ZERO)
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
             ovvo_dims(1) = nocc
             ovvo_dims(2) = nvir
             ovvo_dims(3) = nvir
             ovvo_dims(4) = nocc
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
    
    !        Clean-up
             ierr=talsh_tensor_destruct(aoint_tensor)
    
           end subroutine get_CIS_integrals

end module
