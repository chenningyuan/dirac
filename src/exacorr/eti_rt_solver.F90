module eti_rt_solver
! This module contains the routines needed to solve real-time equations
! basic features and CIS code Written by Xiang Yuan, Shanghai, June 2024

    use tensor_algebra
    use exacorr_tensor_interface
    use exacorr_datatypes
    use exacorr_cfg
    use exacorr_utils
    use exacorr_scaling_parameters
    use eti_common_routines
    
    implicit none

    public eti_cis_rt_solver_driver

    contains
!---------------------------------------------------------------------------------------------------------------
    subroutine eti_cis_rt_solver_driver(exa_input,eVectorsR,eValues)

        use exacorr_global
        use eti_cis_density,   only       : reshape_cis
        use eti_ccrsp_common,  only       : get_prop_integral


        type(exacc_input), intent(in) :: exa_input
        complex(kind=8), intent(inout), pointer :: eVectorsR(:,:)
        complex(kind=8), intent(inout), pointer :: eValues(:)

!   CIS tensors
        type(tensor_t), pointer :: CIS_eigen_tensors(:)
        integer   :: N_states, N_props
        integer   :: i_state, j_state
        integer(INTD), dimension(2) :: cis_eig_dims
        integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
        integer :: i,j,a,b,N_steps    

!   Property integral and external field
        type(eti_orb_prop_integ_tens_t), pointer  :: orb_int_tensors(:)   ! property integral in MO 
        real(kind=8), pointer       :: Fs(:) ! external field    
        integer            :: i_prop
        logical            :: exchange     
        integer            :: field_type   
        type(tensor_t) :: M_prop_inter
        type(tensor_t), pointer :: M_prop_inter_notime(:)  
        type(tensor_t) :: ECIS_tensor    
        integer(INTD), dimension(1) :: R_dim

!   CIS one-electron density matrix
        type(eti_dm_tens_t), pointer   :: Tra_Den_exe_exe_tensors(:,:)
        type(eti_dm_tens_t)  :: Tra_Den_exe_exe_tensor
        type(eti_dm_tens_t), pointer   :: Tra_Den_exe_grd_tensors(:)
        type(eti_dm_tens_t)  :: Tra_Den_exe_grd_tensor
        integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims
        integer(INTD), dimension(2) :: M_dims

!   Propagate
        real(kind=8) t, t0, dt
        type(tensor_t) :: R_tensor
        complex(8) :: td_prop   

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
        N_props = 3    ! Number of properties to be calculated

!   Constrcut cis_eigenvectors
        allocate(CIS_eigen_tensors(N_states))
        do i_state=1, N_states
          ierr=construct_tensor_t(CIS_eigen_tensors(i_state),exa_input%eti_precision%storage_t,cis_eig_dims,init_val=ZERO)
          call reshape_cis(eVectorsR,CIS_eigen_tensors(i_state),cis_eig_dims,i_state)
        end do

!   Get property integrals
        allocate(orb_int_tensors(N_props))
        allocate(Fs(N_props))
        exchange = .false.
        do i=1, N_props
          call get_prop_integral(exa_input,orb_int_tensors(i),i,exchange)
        end do 

!   Constrcut cis_density_matrix
        allocate(Tra_Den_exe_grd_tensors(N_states))
        allocate(Tra_Den_exe_exe_tensors(N_states,N_states))

        do i_state=1, N_states
          call get_cis_den_exe_grd(exa_input,CIS_eigen_tensors(i_state),&
                Tra_Den_exe_grd_tensors(i_state), vo_dims,exa_input%eti_precision)
        end do 

        do j_state=1, N_states
            do i_state = 1, N_states
              call get_cis_den_exe_exe(exa_input,CIS_eigen_tensors(i_state),&
                    CIS_eigen_tensors(j_state),Tra_Den_exe_exe_tensors(i_state,j_state),&
                    oo_dims, vv_dims,exa_input%eti_precision)
            end do
        end do

!   Constrcut time-dependent coeffcients
        R_dim = N_states + 1
        M_dims = N_states + 1
        call form_R_tensor(R_dim,R_tensor,exa_input%eti_precision)
        call form_ECIS_tensor(M_dims,eValues,ECIS_tensor,exa_input%eti_precision)
        ierr=construct_tensor_t(M_prop_inter,exa_input%eti_precision%storage_t,M_dims,init_val=ZERO)
    
        t0 = 0.0
        t = 0.0
        dt = exa_input%rt_dt  
        N_steps = exa_input%rt_Nsteps
        field_type = exa_input%rt_fieldtype
        allocate(M_prop_inter_notime(N_props))
    
        do i_prop=1, N_props
          ierr=construct_tensor_t(M_prop_inter_notime(i_prop),exa_input%eti_precision%storage_t,M_dims,init_val=ZERO)
          call get_submatrix_notime(i_prop,orb_int_tensors,&
            Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
            oo_dims,M_dims,M_prop_inter_notime(i_prop),exa_input%eti_precision)
        end do

!   Propagate R(t)
        open(unit=10, file='TDProp.out', status='new', action='write', iostat=i)

        do i=1, N_steps
          call get_submatrix(t,Fs,orb_int_tensors,&
                        Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
                        oo_dims,M_dims,M_prop_inter,field_type,exa_input%eti_precision)
        
          call propagate_R(dt,ECIS_tensor,M_prop_inter,R_dim,R_tensor,exa_input%eti_precision)
          t = t0+dt
          t0 = t

        !   do i_prop=1, N_props
            i_prop = 3 !calculate one Z-component of the property only
            call gen_time_prop(M_prop_inter_notime(i_prop),R_tensor,R_dim, td_prop)
            write(10, '(F14.12, A, F14.12)')  real(td_prop), ' + ', aimag(td_prop)
        !   end do


          if (i<=10) then
            print *,'rounds', i
            ! do i_prop=1, N_props
                i_prop = 3 !calculate one Z-component of the property only
                print *, 'Time-dependent prop: ', td_prop
            ! end do 
          end if 

        end do

        close(10)

        call print_date('Finishing Real-time CIS computation')

    end subroutine

!---------------------------------------------------------------------------------------------------------------
    subroutine form_R_tensor(R_dim,R_tensor,eti_precision)
        integer(INTD), dimension(1) :: R_dim
        type(tensor_t), intent(inout) :: R_tensor
        type(eti_precision_t), intent(in) :: eti_precision

        complex(8), pointer, contiguous :: h_tens(:)
        integer(C_INT)       :: ierr
        type(C_PTR):: body_p
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max
        
        ierr=construct_tensor_t(R_tensor,eti_precision%storage_t,R_dim,init_val=ZERO)

!       Set HF as the initial state at t=0
        ierr=direct_data_pointer(R_tensor,body_p)
        call c_f_pointer(body_p,h_tens,R_dim) ! to use as a regular Fortran array
        h_tens(1) = 1.0

    end subroutine

!---------------------------------------------------------------------------------------------------------------
    subroutine form_ECIS_tensor(M_dims,eValues,ECIS_tensor,eti_precision)
        integer(INTD), dimension(2) :: M_dims
        complex(kind=8), intent(inout), pointer :: eValues(:)
        type(tensor_t), intent(inout) :: ECIS_tensor
        type(eti_precision_t), intent(in) :: eti_precision

        integer :: i
        complex(8), pointer, contiguous :: h_tens(:,:)
        integer(C_INT)       :: ierr
        type(C_PTR):: body_p
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max


        ierr=construct_tensor_t(ECIS_tensor,eti_precision%storage_t,M_dims,init_val=ZERO)
        ierr=direct_data_pointer(ECIS_tensor,body_p)
        call c_f_pointer(body_p,h_tens,M_dims) ! to use as a regular Fortran array

        h_tens(1,1) = 0.0 
        do i=2, M_dims(1)
            h_tens(i,i) =eValues(i-1) 
        end do

    end subroutine

!---------------------------------------------------------------------------------------------------------------
    function fieldfunc(t,field_type)
        implicit none
        real(kind=8) :: fieldfunc
        real(kind=8), intent(in) :: t
        integer, intent(in)     :: field_type   

        if (field_type == 1) then
            if(t<= 0.00001) then
                fieldfunc = 0.0005
            else 
                fieldfunc = 0.0
            end if
        end if
    end function fieldfunc

!---------------------------------------------------------------------------------------------------------------
    subroutine get_cis_den_exe_grd(exa_input,CIS_eigen_tensor,&
                                    Tra_Den_exe_grd_tensor,vo_dims,eti_precision)

    type(exacc_input), intent(in) :: exa_input
    type(tensor_t),intent(inout) :: CIS_eigen_tensor
    type(eti_dm_tens_t),intent(inout)  :: Tra_Den_exe_grd_tensor
    integer(INTD), dimension(2) :: vo_dims
    type(eti_precision_t), intent(in) :: eti_precision
    

    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    complex(8), pointer, contiguous :: Tra_Den_exe_grd_tens(:,:)
    
    type(tensor_t) :: one_tensor
    integer(C_INT)     :: one_dims(1)

    integer :: i, a

    ierr=construct_tensor_t(Tra_Den_exe_grd_tensor%prime_vo,eti_precision%storage_t,vo_dims,init_val=ZERO)
    ierr=initialize_tensor_t(Tra_Den_exe_grd_tensor%prime_vo)
    one_dims(1) = 1
    ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ONE)
    
    ierr=t_einsum("TDM(a,i)+=CISR(a,i)*I()",Tra_Den_exe_grd_tensor%prime_vo,&
                                CIS_eigen_tensor,one_tensor,scale=ONE)

    ierr=destroy_tensor_t(one_tensor)

    end subroutine 

!---------------------------------------------------------------------------------------------------------------
    subroutine get_cis_den_exe_exe(exa_input,CIS_eigen_left_tensor,& 
                            CIS_eigen_right_tensor,Tra_Den_exe_exe_tensor,&
                            oo_dims, vv_dims,eti_precision)

        type(exacc_input), intent(in) :: exa_input
        type(tensor_t),intent(inout) :: CIS_eigen_left_tensor
        type(tensor_t),intent(inout) :: CIS_eigen_right_tensor
        type(eti_dm_tens_t),intent(inout)  :: Tra_Den_exe_exe_tensor
        type(eti_precision_t), intent(in) :: eti_precision
        integer(INTD), dimension(2) :: oo_dims, vv_dims
    
        type(C_PTR):: body_p
        integer(C_INT)       :: ierr
        complex(8), pointer, contiguous :: Tra_Den_exe_exe_oo_tens(:,:)
        complex(8), pointer, contiguous :: Tra_Den_exe_exe_vv_tens(:,:)

        integer :: i, j, a, b

        ierr=construct_tensor_t(Tra_Den_exe_exe_tensor%prime_vv,eti_precision%storage_t,vv_dims,init_val=ZERO)
        ierr=construct_tensor_t(Tra_Den_exe_exe_tensor%prime_oo,eti_precision%storage_t,oo_dims,init_val=ZERO)

        ierr=initialize_tensor_t(Tra_Den_exe_exe_tensor%prime_oo)
        ierr=initialize_tensor_t(Tra_Den_exe_exe_tensor%prime_vv)
    
        ierr=t_einsum("TDM(a,b)+=CISL+(a,i)*CISR(b,i)",Tra_Den_exe_exe_tensor%prime_vv,&
                                CIS_eigen_left_tensor,CIS_eigen_right_tensor,scale=ONE)

        ierr=t_einsum("TDM(i,j)+=CISL+(a,i)*CISR(a,j)",Tra_Den_exe_exe_tensor%prime_oo,&
                                CIS_eigen_left_tensor,CIS_eigen_right_tensor,scale=ONE)

    end subroutine 

!---------------------------------------------------------------------------------------------------------------
    subroutine get_submatrix(t,Fs,orb_int_tensors,&
                        Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
                        oo_dims,M_dims,M_prop_inter,field_type,eti_precision)

        real(kind=8), intent(in) :: t
        real(kind=8)       :: Fs(:) ! external field   
        type(eti_orb_prop_integ_tens_t), pointer  :: orb_int_tensors(:)   ! property integral in MO 
        type(eti_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_exe_tensors(:,:)
        type(eti_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_grd_tensors(:)
        integer(INTD), dimension(2) :: oo_dims
        integer(INTD), dimension(2) :: M_dims
        type(tensor_t), intent(inout) :: M_prop_inter
        integer, intent(in)            :: field_type   
        type(eti_precision_t), intent(in) :: eti_precision

        integer     :: ierr
        type(C_PTR) :: body_p, body_p2, body_p3
        type(tensor_t)          :: result_tensor
        complex(8), pointer         :: result_tens(:)
        integer(C_INT)              :: one_dims(1)
        complex(8)                  :: com_f
        complex(8), pointer         :: orb_int_tens(:,:)
        complex(8), pointer, contiguous :: M_prop_inter_tens(:,:)

        integer     :: N_props, N_states
        integer     :: i_state, j_state
        integer     :: i, j

        ierr=initialize_tensor_t(M_prop_inter)
        ierr=direct_data_pointer(M_prop_inter,body_p)
        call c_f_pointer(body_p,M_prop_inter_tens,M_dims)

        N_states = size(M_prop_inter_tens,1)
        N_props = size(Fs, 1)

        one_dims(1) = 1
        ierr=construct_tensor_t(result_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ZERO)
        ierr=direct_data_pointer(result_tensor,body_p2)
        call c_f_pointer(body_p2,result_tens,one_dims)


!       CIS-CIS block
        do j_state = 1, N_states-1
            do i_state = j_state, N_states-1        
                ierr=initialize_tensor_t(result_tensor)

                do i=1, N_props
                    Fs(i) = fieldfunc(t,field_type)
                    com_f = cmplx(Fs(i), 0.0)
                    ierr=t_einsum("M()+=Prop(a,b)*DM(a,b)",result_tensor,&
                             orb_int_tensors(i)%prop_vv,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_vv,&
                             scale=com_f)
                
                    ierr=t_einsum("M()+=Prop(i,j)*DM(i,j)",result_tensor,&
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
            ierr=initialize_tensor_t(result_tensor)
            do i=1, N_props
                Fs(i) = fieldfunc(t,field_type)
                com_f = cmplx(Fs(i), 0.0)
                ierr=t_einsum("M()+=Prop(a,i)*DM(a,i)",result_tensor,&
                    orb_int_tensors(i)%prop_vo,&
                    Tra_Den_exe_grd_tensors(i_state)%prime_vo,&
                    scale=com_f)
            end do
            M_prop_inter_tens(i_state+1, 1) = result_tens(1)
            M_prop_inter_tens(1, i_state+1) = conjg(result_tens(1))
                
        end do

!       HF-HF block
        do i=1, N_props
            Fs(i) = fieldfunc(t,field_type)
            com_f = cmplx(Fs(i), 0.0)
            ierr=direct_data_pointer(orb_int_tensors(i)%prop_oo,body_p3)
            call c_f_pointer(body_p3,orb_int_tens,oo_dims)
          do j = 1, oo_dims(1)
            M_prop_inter_tens(1, 1) = M_prop_inter_tens(1, 1) + orb_int_tens(j,j) * com_f
          end do 
        end do

        ierr=destroy_tensor_t(result_tensor)

    end subroutine 

!---------------------------------------------------------------------------------------------------------------

    subroutine get_submatrix_notime(i_prop,orb_int_tensors,&
                        Tra_Den_exe_exe_tensors,Tra_Den_exe_grd_tensors,&
                        oo_dims,M_dims,M_prop_inter_notime,eti_precision)

        integer  :: i_prop
        type(eti_orb_prop_integ_tens_t), pointer  :: orb_int_tensors(:)   ! property integral in MO 
        type(eti_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_exe_tensors(:,:)
        type(eti_dm_tens_t), pointer, intent(inout)  :: Tra_Den_exe_grd_tensors(:)
        integer(INTD), dimension(2) :: oo_dims
        integer(INTD), dimension(2) :: M_dims
        type(tensor_t), intent(inout) :: M_prop_inter_notime
        type(eti_precision_t) :: eti_precision

        integer     :: ierr
        type(C_PTR) :: body_p, body_p2, body_p3       
        type(tensor_t)          :: result_tensor
        complex(8), pointer         :: result_tens(:)
        integer(C_INT)              :: one_dims(1)
        complex(8)                  :: com_f
        complex(8), pointer         :: orb_int_tens(:,:)
        complex(8), pointer, contiguous :: M_prop_inter_tens(:,:)


        integer     :: N_states
        integer     :: i_state, j_state
        integer     :: j


        ierr=initialize_tensor_t(M_prop_inter_notime)
        ierr=direct_data_pointer(M_prop_inter_notime,body_p)
        call c_f_pointer(body_p,M_prop_inter_tens,M_dims)

        N_states = size(M_prop_inter_tens,1)

        one_dims(1) = 1
        ierr=construct_tensor_t(result_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ZERO)
        ierr=direct_data_pointer(result_tensor,body_p2)
        call c_f_pointer(body_p2,result_tens,one_dims)
        

!       CIS-CIS block
        do j_state = 1, N_states-1
            do i_state = j_state, N_states-1
                ierr=initialize_tensor_t(result_tensor)

                    ierr=t_einsum("M()+=Prop(a,b)*DM(a,b)",result_tensor,&
                             orb_int_tensors(i_prop)%prop_vv,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_vv,&
                             scale=ONE)
                    
                    ierr=t_einsum("M()+=Prop(i,j)*DM(i,j)",result_tensor,&
                             orb_int_tensors(i_prop)%prop_oo,&
                             Tra_Den_exe_exe_tensors(i_state, j_state)%prime_oo,&
                             scale=MINUS_ONE)

                M_prop_inter_tens(i_state+1, j_state+1) = result_tens(1)
                M_prop_inter_tens(j_state+1, i_state+1) = conjg(result_tens(1))
                
            end do
        end do 


!       CIS-HF block
        do i_state = 1, N_states-1
            ierr=initialize_tensor_t(result_tensor)
                ierr=t_einsum("M()+=Prop(a,i)*DM(a,i)",result_tensor,&
                    orb_int_tensors(i_prop)%prop_vo,&
                    Tra_Den_exe_grd_tensors(i_state)%prime_vo,&
                    scale=ONE)
            M_prop_inter_tens(i_state+1, 1) = result_tens(1)
            M_prop_inter_tens(1, i_state+1) = conjg(result_tens(1))
        end do


!       HF-HF block
        ierr=direct_data_pointer(orb_int_tensors(i_prop)%prop_oo,body_p3)

        call c_f_pointer(body_p3,orb_int_tens,oo_dims)
        do j = 1, oo_dims(1)
            M_prop_inter_tens(1, 1) = M_prop_inter_tens(1, 1) + orb_int_tens(j,j) 
        end do 


        print *, "Time-independent prop(HF):", M_prop_inter_tens(1, 1)
        print "(A60,I0,A11)", "Print Transition matrix between CIS and HF for ", i_prop, "th operator"
        do i_state = 1, N_states-1
            print "(I3,A11,F10.6,A4,F10.6)", i_state, "---> HF = ", real(M_prop_inter_tens(i_state+1, 1)), &
                                " + i ", aimag(M_prop_inter_tens(i_state+1, 1))
        end do

        ierr=destroy_tensor_t(result_tensor)

    end subroutine

!------------------------------------------------------------------------------------------------------------
    subroutine diagonal_submatrix(r_dims,M_prop_inter_tens,DiaM_prop_inter_tens)

        integer(INTD), dimension(2), intent(in) :: r_dims
        complex(8), pointer, contiguous, intent(inout) :: M_prop_inter_tens(:,:)
        complex(8), pointer, contiguous :: DiaM_prop_inter_tens(:,:)
        
        integer :: lwork, info
        real(kind=8), allocatable :: work(:)
        real(kind=8), allocatable :: rwork(:)
        integer :: vectors_L, vectors_R
        complex(kind=8), pointer :: eVal(:) 
        complex(kind=8), pointer :: eVectorsR(:,:)
        complex(kind=8), pointer :: eVectorsL(:,:)
        character(len=1) :: do_left = 'N', do_right = 'V'

        integer :: i

        DiaM_prop_inter_tens = cmplx(0.0,0.0)

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
    subroutine propagate_R(dt,ECIS_tensor,M_prop_inter,R_dim,R_tensor,eti_precision)

        type(tensor_t), intent(inout) :: M_prop_inter
        type(tensor_t), intent(inout) :: ECIS_tensor
        type(tensor_t), intent(inout) :: R_tensor
        type(eti_precision_t), intent(in) :: eti_precision
        real(kind=8) dt
        integer(INTD), dimension(1) :: R_dim

        integer     :: ierr
        type(C_PTR) :: body_p
        complex(8), pointer  :: R_tens(:)
        type(tensor_t)   :: R_temp_tensor, R_temp_tensor2
        complex(8)           :: com_f_real, com_f_imag


        ierr=direct_data_pointer(R_tensor,body_p)
        call c_f_pointer(body_p,R_tens,R_dim)

        com_f_real = cmplx(dt, 0.0)
        com_f_imag = cmplx(0.0, dt)

        ierr=construct_tensor_t(R_temp_tensor,eti_precision%storage_t,R_dim,init_val=ZERO)
        ierr=construct_tensor_t(R_temp_tensor2,eti_precision%storage_t,R_dim,init_val=ZERO)

        ierr=t_einsum("RTemp(i)+=ECIS(i,j)*R(j)",R_temp_tensor,&
                             ECIS_tensor,R_tensor,scale=-com_f_imag)

        ierr=t_einsum("R(i)+=RTemp(i)",R_tensor,R_temp_tensor)


        ierr=initialize_tensor_t(R_temp_tensor)
        ierr=t_einsum("RTemp(i)+=M(i,j)*R(j)",R_temp_tensor,&
                             M_prop_inter,R_tensor,scale=com_f_imag)

        ierr=t_einsum("R(i)+=RTemp(i)",R_tensor,R_temp_tensor)


        ierr=initialize_tensor_t(R_temp_tensor)
        ierr=t_einsum("RTemp(i)+=ECIS(i,j)*R(j)",R_temp_tensor,&
                             ECIS_tensor,R_tensor,scale=com_f_real)


        ierr=t_einsum("R(i)+=M(i,j)*RTemp(j)",R_tensor,&
                             M_prop_inter,R_temp_tensor,scale=com_f_real)


        ierr=initialize_tensor_t(R_temp_tensor)
        ierr=t_einsum("RTemp(i)+=ECIS(i,j)*R(j)",R_temp_tensor,&
                             ECIS_tensor,R_tensor,scale=com_f_real)


        ierr=t_einsum("R(i)+=ECIS(i,j)*RTemp(j)",R_tensor,&
                             ECIS_tensor,R_temp_tensor,scale=com_f_real*0.5)
       
        
        ierr=initialize_tensor_t(R_temp_tensor)
        ierr=t_einsum("RTemp(i)+=M(i,j)*R(j)",R_temp_tensor,&
                             M_prop_inter,R_tensor,scale=com_f_real)


        ierr=t_einsum("R(i)+=M(i,j)*RTemp(j)",R_tensor,&
                             M_prop_inter,R_temp_tensor,scale=com_f_real*0.5)


        ierr=initialize_tensor_t(R_temp_tensor)
        ierr=t_einsum("RTemp(i)+=ECIS(i,j)*R(j)",R_temp_tensor,&
                             ECIS_tensor,R_tensor,scale=com_f_real)
        ierr=t_einsum("RTemp2(i)+=ECIS(i,j)*Rtemp1(j)",R_temp_tensor2,&
                             ECIS_tensor,R_temp_tensor,scale=com_f_real)

        ierr=t_einsum("R(i)+=M(i,j)*RTemp2(j)",R_tensor,&
                             M_prop_inter,R_temp_tensor2,scale=-(com_f_imag *0.5))


        ierr=initialize_tensor_t(R_temp_tensor)
        ierr=initialize_tensor_t(R_temp_tensor2)
        ierr=t_einsum("RTemp(i)+=ECIS(i,j)*R(j)",R_temp_tensor,&
                             ECIS_tensor,R_tensor,scale=com_f_real)
        ierr=t_einsum("RTemp2(i)+=M(i,j)*Rtemp1(j)",R_temp_tensor2,&
                             M_prop_inter,R_temp_tensor,scale=com_f_real)

        ierr=t_einsum("R(i)+=M(i,j)*RTemp2(j)",R_tensor,&
                             M_prop_inter,R_temp_tensor2,scale=(com_f_imag *0.5))

        ierr=destroy_tensor_t(R_temp_tensor)
        ierr=destroy_tensor_t(R_temp_tensor2)

    end subroutine 

!---------------------------------------------------------------------------------------------------------------
    subroutine gen_time_prop(M_prop_inter_notime,R_tensor,R_dim, td_prop)

        type(tensor_t), intent(inout) :: M_prop_inter_notime
        type(tensor_t), intent(inout) :: R_tensor 
        integer(INTD), dimension(1) :: R_dim
        complex(8) :: td_prop   

        integer     :: ierr
        type(C_PTR) :: body_p
        complex(8), pointer         :: result_tens(:)
        integer(C_INT)     :: one_dims(1)
        type(tensor_t) :: temp_tensor, result_tensor
        complex(8) :: norm
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = R_tensor%datatype

        one_dims(1) = 1
        ierr=construct_tensor_t(result_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ZERO)
        ierr=direct_data_pointer(result_tensor,body_p)
        call c_f_pointer(body_p,result_tens,one_dims)

        ierr=t_einsum("Norm()+=R+(i)*R(i)",result_tensor,&
                            R_tensor,R_tensor,scale=ONE)

        norm = result_tens(1)

        ierr=initialize_tensor_t(result_tensor)
        ierr=construct_tensor_t(temp_tensor,eti_precision%storage_t,R_dim,init_val=ZERO)
        
        ierr=t_einsum("Temp(i)+=M(i,j)*R(j)",temp_tensor,&
                            M_prop_inter_notime,R_tensor,scale=ONE)

        ierr=t_einsum("Result()+=R+(i)*Temp(i)",result_tensor,&
                            R_tensor,temp_tensor,scale=ONE)

        td_prop = result_tens(1)/norm

        ierr=destroy_tensor_t(temp_tensor)
        ierr=destroy_tensor_t(result_tensor)

    end subroutine


end module
