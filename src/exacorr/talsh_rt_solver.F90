module talsh_rt_solver
! This module contains the routines needed to solve real-time equations
! Written by Xiang Yuan, Shanghai, June 2024
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    
    
    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)

    public talsh_cis_rt_solver_driver

    contains

!---------------------------------------------------------------------------------------------------------------

subroutine talsh_cis_rt_solver_driver(exa_input,eVectorsR,eValues)

    ! use talsh_ao_to_mo
    use exacorr_global
    use talsh_cis_density,   only       : reshape_cis
    use talsh_ccrsp_common,  only       : get_prop_integral


    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), intent(inout), pointer :: eVectorsR(:,:)
    complex(kind=8), intent(inout), pointer :: eValues(:)

!   CIS tensors
    type(talsh_tens_t), pointer :: CIS_eigen_tensors(:)
    integer   :: N_states, N_props
    integer   :: i_state, j_state
    integer(INTD), dimension(2) :: cis_eig_dims
    integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
    integer :: i,j,a,b,N_steps    

!   Property integral and external field
    type(talsh_orb_prop_integ_tens_t), pointer  :: orb_int_tensors(:)   ! property integral in MO 
    real(kind=8), pointer       :: Fs(:) ! external field    
    integer            :: i_prop
    logical            :: exchange     
    integer            :: field_type   
    type(talsh_tens_t) :: M_prop_inter, M_prop_inter_notime    
    ! type(talsh_tens_t), intent(inout) :: DiaM_prop_inter
    type(talsh_tens_t) :: ECIS_tensor    
    integer(INTD), dimension(1) :: R_dim

!   CIS one-electron density matrix
    type(talsh_dm_tens_t), pointer   :: Tra_Den_exe_exe_tensors(:,:)
    type(talsh_dm_tens_t)  :: Tra_Den_exe_exe_tensor
    type(talsh_dm_tens_t), pointer   :: Tra_Den_exe_grd_tensors(:)
    type(talsh_dm_tens_t)  :: Tra_Den_exe_grd_tensor
    integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims
    integer(INTD), dimension(2) :: M_dims

!   Propagate
    real(kind=8) t, t0, dt
    type(talsh_tens_t) :: R_tensor

!   Talsh 
    integer(C_INT)       :: ierr
    type(C_PTR):: body_p
    integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
    integer(C_INT):: host_arg_max

    nvir = exa_input%nvir
    nocc = exa_input%nocc
    cis_eig_dims(1) = nvir
    cis_eig_dims(2) = nocc
    vv_dims = nvir
    oo_dims = nocc
    vo_dims(1) = nvir
    vo_dims(2) = nocc


    N_states = size(eVectorsR,2)
    N_props = 1

!   Constrcut talsh cis_eigenvectors
    allocate(CIS_eigen_tensors(N_states))
    do i_state=1, N_states
        ierr=talsh_tensor_construct(CIS_eigen_tensors(i_state),C8,cis_eig_dims,init_val=ZERO)
        ! print *, "ierr_CIS_eigen: ", ierr
        call reshape_cis(eVectorsR,CIS_eigen_tensors(i_state),cis_eig_dims,i_state)
        ! print *, "issue in reshape?"
    end do

!   Get property integrals
    allocate(orb_int_tensors(N_props))
    allocate(Fs(N_props))
    i_prop = 3
    exchange = .false.
    do i=1, N_props
        call get_prop_integral(exa_input,orb_int_tensors(i),i_prop,exchange)
        ! print *, "issues?"
    end do 

!   Constrcut talsh cis_density_matrix
    allocate(Tra_Den_exe_grd_tensors(N_states))
    allocate(Tra_Den_exe_exe_tensors(N_states,N_states))
    ierr=talsh_tensor_construct(Tra_Den_exe_grd_tensor%prime_vo,C8,vo_dims,init_val=ZERO)
    ! print *, "ierrtra_exe_grd: ", ierr
    ierr=talsh_tensor_construct(Tra_Den_exe_exe_tensor%prime_vv,C8,vv_dims,init_val=ZERO)
    ! print *, "ierrtra_exe_exe1: ", ierr
    ierr=talsh_tensor_construct(Tra_Den_exe_exe_tensor%prime_oo,C8,oo_dims,init_val=ZERO)
    ! print *, "ierrtra_exe_exe2: ", ierr

    do i_state = 1, N_states
        call get_cis_den_exe_grd(exa_input,CIS_eigen_tensors(i_state),&
                                    Tra_Den_exe_grd_tensor, vo_dims)
        Tra_Den_exe_grd_tensors(i_state) = Tra_Den_exe_grd_tensor
    end do 

    do j_state = 1, N_states
        do i_state = j_state, N_states
            call get_cis_den_exe_exe(exa_input,CIS_eigen_tensors(i_state),&
                                        CIS_eigen_tensors(j_state),Tra_Den_exe_exe_tensor,&
                                        oo_dims, vv_dims)
            Tra_Den_exe_exe_tensors(i_state,j_state) = Tra_Den_exe_exe_tensor

            if (i_state > j_state) then
                Tra_Den_exe_exe_tensors(j_state,i_state) = Tra_Den_exe_exe_tensors(i_state,j_state)
            end if 
        end do
    end do

!   Constrcut talsh time-dependent coeffcients
    R_dim = N_states + 1
    M_dims = N_states + 1
    call form_R_tensor(R_dim,R_tensor)
    call form_ECIS_tensor(M_dims,eValues,ECIS_tensor)
    ierr=talsh_tensor_construct(M_prop_inter,C8,M_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(M_prop_inter_notime,C8,M_dims,init_val=ZERO)


    t0 = 0.0
    t = 0.0
    dt = exa_input%rt_dt  
    N_steps = exa_input%rt_Nsteps
    field_type = exa_input%rt_fieldtype

    i_prop = 1    !to calculaye the time-dependent ith prop
    
    call get_submatrix_notime(i_prop,orb_int_tensors,&
            Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
            oo_dims,M_dims,M_prop_inter_notime)

!   Propagate R
    do i=1, N_steps
        call get_submatrix(t,Fs,orb_int_tensors,&
                        Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
                        oo_dims,M_dims,M_prop_inter,field_type)
        print *,'rounds', i
        call propagate_R(dt,ECIS_tensor,M_prop_inter,R_dim,R_tensor)
        t = t0+dt
        t0 = t

        call gen_time_prop(M_prop_inter_notime,R_tensor,R_dim)
    end do

    call print_date('Finishing Real-time CIS computation')

end subroutine
!---------------------------------------------------------------------------------------------------------------

    subroutine form_R_tensor(R_dim,R_tensor)
        integer(INTD), dimension(1) :: R_dim
        type(talsh_tens_t), intent(inout) :: R_tensor

        
        complex(8), pointer, contiguous :: h_tens(:)
        integer(C_INT)       :: ierr
        type(C_PTR):: body_p
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max
        

        ierr=talsh_tensor_construct(R_tensor,C8,R_dim,init_val=ZERO)
        ! print *, "ierr R1:", ierr

!       Set HF as the initial state at t=0
        ierr=talsh_tensor_get_body_access(R_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        ! print *, "ierr R2:", ierr
        call c_f_pointer(body_p,h_tens,R_dim) ! to use as a regular Fortran array
        h_tens(1) = 1.0

    end subroutine

!---------------------------------------------------------------------------------------------------------------

    subroutine form_ECIS_tensor(M_dims,eValues,ECIS_tensor)
        integer(INTD), dimension(2) :: M_dims
        complex(kind=8), intent(inout), pointer :: eValues(:)
        type(talsh_tens_t), intent(inout) :: ECIS_tensor


        integer :: i
        complex(8), pointer, contiguous :: h_tens(:,:)
        integer(C_INT)       :: ierr
        type(C_PTR):: body_p
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max
        complex(kind=8)  :: e_hf

    
        ierr=talsh_tensor_construct(ECIS_tensor,C8,M_dims,init_val=ZERO)
        ! print *, "ierr ECIS1:", ierr

        ierr=talsh_tensor_get_body_access(ECIS_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        ! print *, "ierr ECIS2:", ierr
        call c_f_pointer(body_p,h_tens,M_dims) ! to use as a regular Fortran array

        e_hf = -7.985
        h_tens(1,1) = 0.0 + e_hf
        do i=2, M_dims(1)
            h_tens(i,i) =eValues(i-1) + e_hf
        end do

    end subroutine

!---------------------------------------------------------------------------------------------------------------
    function fieldfunc(t,field_type)
        implicit none
        real(kind=8) :: fieldfunc
        real(kind=8), intent(in) :: t
        integer, intent(in)     :: field_type   

        if (field_type == 1) then
            fieldfunc = 0.0005
        end if
    end function fieldfunc
!---------------------------------------------------------------------------------------------------------------
    subroutine get_cis_den_exe_grd(exa_input,CIS_eigen_tensor,&
                                    Tra_Den_exe_grd_tensor,vo_dims)

    type(exacc_input), intent(in) :: exa_input
    type(talsh_tens_t),intent(inout) :: CIS_eigen_tensor
    type(talsh_dm_tens_t),intent(inout)  :: Tra_Den_exe_grd_tensor
    integer(INTD), dimension(2) :: vo_dims
    

    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    complex(8), pointer, contiguous :: Tra_Den_exe_grd_tens(:,:)
    
    type(talsh_tens_t) :: one_tensor
    integer(C_INT)     :: one_dims(1)

    integer :: i, a

    ierr=talsh_tensor_init(Tra_Den_exe_grd_tensor%prime_vo)
    ! print *, "ierr_get_tra_exe_grd1: ", ierr
    one_dims(1) = 1
    ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
    ! print *, "ierr_get_tra_exe_grd2: ", ierr
    
    ierr=talsh_tensor_contract("TDM(a,i)+=CISR(a,i)*I()",Tra_Den_exe_grd_tensor%prime_vo,&
                                CIS_eigen_tensor,one_tensor,scale=ONE)

    ierr=talsh_tensor_destruct(one_tensor)

    end subroutine 
!---------------------------------------------------------------------------------------------------------------

    subroutine get_cis_den_exe_exe(exa_input,CIS_eigen_left_tensor,& 
                            CIS_eigen_right_tensor,Tra_Den_exe_exe_tensor,&
                            oo_dims, vv_dims)

        type(exacc_input), intent(in) :: exa_input
        type(talsh_tens_t),intent(inout) :: CIS_eigen_left_tensor
        type(talsh_tens_t),intent(inout) :: CIS_eigen_right_tensor
        type(talsh_dm_tens_t),intent(inout)  :: Tra_Den_exe_exe_tensor
        integer(INTD), dimension(2) :: oo_dims, vv_dims
    

        type(C_PTR):: body_p
        integer(C_INT)       :: ierr
        complex(8), pointer, contiguous :: Tra_Den_exe_exe_oo_tens(:,:)
        complex(8), pointer, contiguous :: Tra_Den_exe_exe_vv_tens(:,:)

        integer :: i, j, a, b

        ierr=talsh_tensor_init(Tra_Den_exe_exe_tensor%prime_oo)
        ! print *, "ierr_get_tra_exe_exe1: ", ierr
        ierr=talsh_tensor_init(Tra_Den_exe_exe_tensor%prime_vv)
        ! print *, "ierr_get_tra_exe_exe2: ", ierr
    
        ierr=talsh_tensor_contract("TDM(a,b)+=CISL+(a,i)*CISR(b,i)",Tra_Den_exe_exe_tensor%prime_vv,&
                                CIS_eigen_left_tensor,CIS_eigen_right_tensor,scale=ONE)

        ! print *, "ierr_get_tra_exe_exe3: ", ierr
        ierr=talsh_tensor_contract("TDM(i,j)+=CISL+(a,i)*CISR(a,j)",Tra_Den_exe_exe_tensor%prime_oo,&
                                CIS_eigen_left_tensor,CIS_eigen_right_tensor,scale=ONE)

    end subroutine 
!---------------------------------------------------------------------------------------------------------------

    subroutine get_submatrix(t,Fs,orb_int_tensors,&
                        Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
                        oo_dims,M_dims,M_prop_inter,field_type)

        real(kind=8), intent(in) :: t
        real(kind=8)       :: Fs(:) ! external field   
        type(talsh_orb_prop_integ_tens_t), pointer  :: orb_int_tensors(:)   ! property integral in MO 
        type(talsh_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_exe_tensors(:,:)
        type(talsh_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_grd_tensors(:)
        integer(INTD), dimension(2) :: oo_dims
        integer(INTD), dimension(2) :: M_dims
        type(talsh_tens_t), intent(inout) :: M_prop_inter
        integer, intent(in)            :: field_type   



        integer     :: N_props, N_states
        integer     :: i_state, j_state
        integer     :: i, j

        integer     :: ierr
        type(C_PTR) :: body_p, body_p2, body_p3
       
        type(talsh_tens_t)          :: result_tensor
        complex(8), pointer         :: result_tens(:)
        integer(C_INT)              :: one_dims(1)

        complex(8)                  :: com_f
        complex(8), pointer         :: orb_int_tens(:,:)
        complex(8), pointer, contiguous :: M_prop_inter_tens(:,:)

        ierr=talsh_tensor_init(M_prop_inter)
        ierr=talsh_tensor_get_body_access(M_prop_inter,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,M_prop_inter_tens,M_dims)

        N_states = size(M_prop_inter_tens,1)
        N_props = size(Fs, 1)


        one_dims(1) = 1
        ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
        ierr=talsh_tensor_get_body_access(result_tensor,body_p2,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p2,result_tens,one_dims)

!       CIS-CIS block
        do j_state = 1, N_states-1
            do i_state = j_state, N_states-1        
                ierr=talsh_tensor_init(result_tensor)

                do i=1, N_props
                    Fs(i) = fieldfunc(t,field_type)
                    com_f = complex(Fs(i), 0.0)
                    ierr=talsh_tensor_contract("M()+=Prop(a,b)*DM(a,b)",result_tensor,&
                             orb_int_tensors(i)%prop_vv,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_vv,&
                             scale=com_f)
                
                    ierr=talsh_tensor_contract("M()+=Prop(i,j)*DM(i,j)",result_tensor,&
                             orb_int_tensors(i)%prop_oo,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_oo,&
                             scale=-com_f)

                end do

                M_prop_inter_tens(i_state+1, j_state+1) = result_tens(1)
                M_prop_inter_tens(j_state+1, i_state+1) = conjg(result_tens(1))
                
            end do
        end do

!       CIS-HF block
        do i_state = 1, N_states-1
            ierr=talsh_tensor_init(result_tensor)
            do i=1, N_props
                Fs(i) = fieldfunc(t,field_type)
                com_f = complex(Fs(i), 0.0)
                ierr=talsh_tensor_contract("M()+=Prop(a,i)*DM(a,i)",result_tensor,&
                    orb_int_tensors(i)%prop_vo,&
                    Tra_Den_exe_grd_tensors(i_state)%prime_vo,&
                    scale=com_f)
            end do
            ! print *, "M: ", result_tens(1)
            M_prop_inter_tens(i_state+1, 1) = result_tens(1)
            M_prop_inter_tens(1, i_state+1) = conjg(result_tens(1))
                
        end do

!       HF-HF block
        do i=1, N_props
            Fs(i) = fieldfunc(t,field_type)
            com_f = complex(Fs(i), 0.0)
            ierr=talsh_tensor_get_body_access(orb_int_tensors(i)%prop_oo,body_p3,&
                                            C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p3,orb_int_tens,oo_dims)
          do j = 1, oo_dims(1)
            M_prop_inter_tens(1, 1) = M_prop_inter_tens(1, 1) + orb_int_tens(j,j) * com_f
          end do 
        end do

        ! call print_tensor(M_prop_inter,1.0d-6,"M_prop_inter")

    end subroutine 
!---------------------------------------------------------------------------------------------------------------

    subroutine get_submatrix_notime(i_prop,orb_int_tensors,&
                        Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
                        oo_dims,M_dims,M_prop_inter_notime)

        integer  :: i_prop
        type(talsh_orb_prop_integ_tens_t), pointer  :: orb_int_tensors(:)   ! property integral in MO 
        type(talsh_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_exe_tensors(:,:)
        type(talsh_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_grd_tensors(:)
        integer(INTD), dimension(2) :: oo_dims
        integer(INTD), dimension(2) :: M_dims
        type(talsh_tens_t), intent(inout) :: M_prop_inter_notime



        integer     :: N_states
        integer     :: i_state, j_state
        integer     :: j

        integer     :: ierr
        type(C_PTR) :: body_p, body_p2, body_p3
       
        type(talsh_tens_t)          :: result_tensor
        complex(8), pointer         :: result_tens(:)
        integer(C_INT)              :: one_dims(1)

        complex(8)                  :: com_f
        complex(8), pointer         :: orb_int_tens(:,:)
        complex(8), pointer, contiguous :: M_prop_inter_tens(:,:)

        ierr=talsh_tensor_init(M_prop_inter_notime)
        ierr=talsh_tensor_get_body_access(M_prop_inter_notime,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,M_prop_inter_tens,M_dims)

        N_states = size(M_prop_inter_tens,1)

        one_dims(1) = 1
        ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
        ierr=talsh_tensor_get_body_access(result_tensor,body_p2,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p2,result_tens,one_dims)
        

!       CIS-CIS block
        do j_state = 1, N_states-1
            do i_state = j_state, N_states-1
                ierr=talsh_tensor_init(result_tensor)

                    ierr=talsh_tensor_contract("M()+=Prop(a,b)*DM(a,b)",result_tensor,&
                             orb_int_tensors(i_prop)%prop_vv,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_vv,&
                             scale=ONE)
                    
                    ierr=talsh_tensor_contract("M()+=Prop(i,j)*DM(i,j)",result_tensor,&
                             orb_int_tensors(i_prop)%prop_oo,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_oo,&
                             scale=MINUS_ONE)


                M_prop_inter_tens(i_state+1, j_state+1) = result_tens(1)
                M_prop_inter_tens(j_state+1, i_state+1) = conjg(result_tens(1))
                
            end do
        end do 


!       CIS-HF block
        do i_state = 1, N_states-1
            ierr=talsh_tensor_init(result_tensor)
                ierr=talsh_tensor_contract("M()+=Prop(a,i)*DM(a,i)",result_tensor,&
                    orb_int_tensors(i_prop)%prop_vo,&
                    Tra_Den_exe_grd_tensors(i_state)%prime_vo,&
                    scale=ONE)
            M_prop_inter_tens(i_state+1, 1) = result_tens(1)
            M_prop_inter_tens(1, i_state+1) = conjg(result_tens(1))
                
        end do

!       HF-HF block
        ierr=talsh_tensor_get_body_access(orb_int_tensors(i_prop)%prop_oo,body_p3,&
                                            C8,int(0,C_INT),DEV_HOST)

        call c_f_pointer(body_p3,orb_int_tens,oo_dims)
        do j = 1, oo_dims(1)
            M_prop_inter_tens(1, 1) = M_prop_inter_tens(1, 1) + orb_int_tens(j,j) 
        end do 

    end subroutine

!------------------------------------------------------------------------------------------------------------
    subroutine diagonal_submatrix(r_dims,M_prop_inter_tens,DiaM_prop_inter_tens)

        integer(INTD), dimension(2), intent(in) :: r_dims
        complex(8), pointer, contiguous, intent(inout) :: M_prop_inter_tens(:,:)
        complex(8), pointer, contiguous :: DiaM_prop_inter_tens(:,:)

        integer     :: ierr
        type(C_PTR) :: body_p
        
        integer :: i,j


        integer :: lwork, info
        real(kind=8), allocatable :: work(:)
        real(kind=8), allocatable :: rwork(:)
        integer :: vectors_L, vectors_R
        complex(kind=8), pointer :: eVal(:) 
        complex(kind=8), pointer :: eVectorsR(:,:)
        complex(kind=8), pointer :: eVectorsL(:,:)
        character(len=1) :: do_left = 'N', do_right = 'V'


        DiaM_prop_inter_tens = complex(0.0,0.0)


        vectors_R = size(M_prop_inter_tens,1)
        vectors_L = size(M_prop_inter_tens,2)
        lwork = 100000*vectors_R
        allocate(eVectorsL(vectors_L,vectors_R))
        allocate(eVectorsR(vectors_L,vectors_R))
        allocate(eVal(vectors_R))
        allocate(work(MAX(1,LWORK)))
        allocate(rwork(2*vectors_R))
        call zgeev (do_left, do_right, vectors_R, M_prop_inter_tens, vectors_L, &
        eVal, eVectorsL, vectors_L, eVectorsR, vectors_R, &
        WORK, LWORK, rwork, INFO)

        do i=1, r_dims(1)
            DiaM_prop_inter_tens(i,i) = eVal(i)
        end do 


    end subroutine 
!---------------------------------------------------------------------------------------------------------------

    subroutine propagate_R(dt,ECIS_tensor,M_prop_inter,R_dim,R_tensor)

        type(talsh_tens_t), intent(inout) :: M_prop_inter
        type(talsh_tens_t), intent(inout) :: ECIS_tensor
        type(talsh_tens_t), intent(inout) :: R_tensor
        real(kind=8) dt
        integer(INTD), dimension(1) :: R_dim

        integer :: i,j

        integer     :: ierr
        type(C_PTR) :: body_p
        complex(8), pointer         :: R_tens(:)
        

        type(talsh_tens_t)   :: R_temp_tensor
        complex(8)           :: com_f

        ierr=talsh_tensor_get_body_access(R_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,R_tens,R_dim)
        ! print *, 'ierr1: ', ierr

        com_f = complex(dt, 0.0)
        ierr=talsh_tensor_construct(R_temp_tensor,C8,R_dim,init_val=ZERO)
        ! print *, 'ierr2: ', ierr

        ierr=talsh_tensor_contract("RTemp(i)+=ECIS(i,j)*R(j)",R_temp_tensor,&
                             ECIS_tensor,R_tensor,scale=com_f)


        ierr=talsh_tensor_contract("R(i)+=M(i,j)*RTemp(j)",R_tensor,&
                             M_prop_inter,R_temp_tensor,scale=com_f)

       
        ierr=talsh_tensor_destruct(R_temp_tensor)

    end subroutine 
!---------------------------------------------------------------------------------------------------------------

    subroutine gen_time_prop(M_prop_inter_notime,R_tensor,R_dim)

        type(talsh_tens_t) :: M_prop_inter_notime
        type(talsh_tens_t) :: R_tensor 
        integer(INTD), dimension(1) :: R_dim   


        integer     :: ierr
        type(C_PTR) :: body_p
        complex(8), pointer         :: result_tens(:)
        integer(C_INT)     :: one_dims(1)
        type(talsh_tens_t) :: temp_tensor, result_tensor


        one_dims(1) = 1
        ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ONE)
        ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,result_tens,one_dims)

        ierr=talsh_tensor_construct(temp_tensor,C8,R_dim,init_val=ZERO)

        ierr=talsh_tensor_contract("Temp(i)+=M(i,j)*R(j)",temp_tensor,&
                            M_prop_inter_notime,R_tensor,scale=ONE)

        ierr=talsh_tensor_contract("Result()+=R(i)*Temp(i)",result_tensor,&
                            R_tensor,temp_tensor,scale=ONE)

        print *, 'Time-dependent prop: ', result_tens(1)

        ierr=talsh_tensor_destruct(temp_tensor)
        ierr=talsh_tensor_destruct(result_tensor)

    end subroutine


end module
