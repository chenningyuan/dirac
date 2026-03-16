module eti_davidson

use exacorr_tensor_interface
use tensor_algebra
use eti_common_routines
use exacorr_datatypes
use exacorr_scaling_parameters
use mfsolver_cfg

implicit none

public eti_davidson_get_eigenvalue_shift
public eti_davidson_setup
public eti_davidson_destroy
public eti_davidson_driver
public eti_solve_values_vectors
public eti_sort_values_vectors
public eti_davidson_workspace_create
public eti_davidson_workspace_destroy

interface eti_davidson_driver
    module procedure eti_get_values_vectors_davidson
end interface

interface davidson_print_vectors
    module procedure davidson_print_vectors_real
    module procedure davidson_print_vectors_complex
end interface

interface print_values_vectors
     module procedure print_values_vectors_R
     module procedure print_values_vectors_C
end interface

#ifdef NOT_USED
interface davidson_print_vectors_complex
    module procedure davidson_print_vectors_complex_1d
end interface
#endif

public eti_davidson_results_tens_t
public eti_davidson_workspace_tens_t


private

type(mfsolver_config_t), save :: config
integer, parameter :: max_label_length = 30

! derived types

type eti_davidson_results_tens_t
    complex(kind=8), pointer :: values_rhs(:,:) 
    complex(kind=8), pointer :: values_lhs(:,:) ! keep left and right energies to verify whether they are the same 
    type(eti_wavefun_tens_t), pointer :: wf(:)
end type eti_davidson_results_tens_t

! this is the workspace storage for a given instance of the solver
type :: eti_davidson_workspace_tens_t
   integer :: nroots              ! number of solutions we're looking for
   integer, pointer :: N(:)       ! number of rows per excitation level
   integer :: Nmax                ! total number of rows =  dim(T1)+dim(T2)+...
   integer :: Lmax                ! total number of columns = dim(maximum number of trial vectors)
   integer :: L, Lprev            ! counters for keeping tabs on how much of the workspace we've used
                                  ! in most usages N = Nmax and L will grow up to Lmax
                                  ! with Lprev keeping the number of the prior iteration (staring at 0)
! aspg, 24/03/2020
! original idea was to that  to M to be a tensor, but using arrays of tensors as below, better to keep it as a pointer
!          type(tensor_t) :: M        ! subspace matrix, dimension : (Lmax, Lmax); it is in tensor 
!                                         ! form to conform with the operations, but it will be a pointer to it 
!                                         ! that will be passed to the lapack routines 
   complex(8), pointer, contiguous:: M(:,:) !  subspace matrix, dimension : (Lmax, Lmax)
! aspg, added from (q)rsp
   complex(8), pointer, contiguous:: C(:)   !  subspace vector, dimension : (Lmax)
! aspg, added from (q)rsp
! see to change those below to tensors segmented into singly/doubly/... manifolds
   type(eti_sa_tens_t), pointer :: trial(:)             ! workspace for trial vectors
   type(eti_sa_tens_t), pointer :: sigma(:)             ! workspace for sigma vectors
   type(eti_sa_tens_t), pointer :: delta(:)             ! workspace for correction vectors
   type(eti_sa_tens_t), pointer :: precond(:)           ! workspace for preconditioner
   type(eti_sa_tens_t), pointer :: metric(:)            ! workspace for metric, for generalized eigenvalue problem, or linear equation problem 
   type(eti_sa_tens_t)          :: Hdiag                ! workspace for diagonal elements of the similarity-transformed Hamiltonian
   type(eti_sa_tens_t)          :: Norm_prop_grad       ! workspace for norm of each element of property gradient vectors C in LR
   integer                      :: tens1_rank           ! tensors for a given model IP/EA/EE/... have the same rank, keep this info around
   integer                      :: tens2_rank
   integer                      :: tens3_rank
   integer                      :: tens4_rank
   integer                      :: tens5_rank
   integer, pointer             :: tens1_dims(:)        ! tensors for a given model IP/EA/EE/... have the same dimensions, keep this info around
   integer, pointer             :: tens2_dims(:)
   integer, pointer             :: tens3_dims(:)
   integer, pointer             :: tens4_dims(:)
   integer, pointer             :: tens5_dims(:)
end type eti_davidson_workspace_tens_t


contains 

! this should never be used to update an instance, only to create new ones. it's better to have results
! for left and right separately if we determine these in two steps. but we must check whether an 
! element of the derived type has been allocated already... perhaps this should then be changed to a
! function returning a results type, and have the davidson driver return that to its caller
    subroutine davidson_create_results(results, N, nroots, values, left, right)
        type(eti_davidson_results_tens_t), intent(inout) :: results
        integer, intent(in)    :: N, nroots 
        integer                :: N_local
        logical, intent(in)    :: values, left, right

! we don't use complex variables but rather use twice the amount of real variables
! for symmetric matrices we only need one of the vectors, and we choose it to be the right one
! left vectors, if asked, are ignored
        N_local = N

!       results%eValues => NULL()
!       results%eVectorsR => NULL()
!       results%eVectorsL => NULL()


!       if (config%symmetric) then
!           if (values) then 
!               allocate(results%eValues(1,nroots))
!               results%eValues   = 0.0d0
!           end if
!           if (right) then
!               allocate(results%eVectorsR(N_local,nroots))
!               results%eVectorsR = 0.0d0
!           end if
!       else
!           if (values) then
!               allocate(results%eValues(2,nroots))
!               results%eValues   = 0.0d0
!           end if
!           if (right) then
!               allocate(results%eVectorsR(N_local,nroots))
!               results%eVectorsR = 0.0d0
!           else
!               allocate(results%eVectorsR(1,nroots))
!               results%eVectorsR = 0.0d0
!           end if
!           if (left) then
!               allocate(results%eVectorsL(N_local,nroots))
!               results%eVectorsL = 0.0d0
!           else
!               allocate(results%eVectorsL(1,nroots))
!               results%eVectorsL = 0.0d0
!           end if
!       end if
    end subroutine


    subroutine free_davidson_results(results)
        type(eti_davidson_results_tens_t), intent(inout) :: results

!       if (associated(results%eValues)) deallocate(results%eValues)
!       if (associated(results%eVectorsR)) deallocate(results%eVectorsR)
!       if (associated(results%eVectorsL)) deallocate(results%eVectorsL)
    end subroutine


    subroutine eti_davidson_setup(convergence_threshold, &
                              max_subspace_size, &
                              max_iterations,    & 
                              unit_output,       & 
                              unit_sigmas,       & 
                              unit_subspaces,    &
                              unit_results,      &
                              solve_values,      &
                              solve_linear_sys,  &
                              solve_right,       &
                              solve_left,        &
                              symmetric,         &
                              complex_mode,      &
                              use_preconditioner,&
                              refresh_trial_rate,&
                              verbose,           & 
                              overlap_sorting,   &
                              overlap_sorting_scheme,  &
                              save_results,      &
                              save_sigmas,       &
                              save_subspaces,    &
                              energy_shift)

        real(kind=8), intent(in), optional :: convergence_threshold
        real(kind=8), intent(in), optional :: energy_shift(2)
        integer, intent(in), optional :: max_subspace_size 
        integer, intent(in), optional :: max_iterations
        integer, intent(in), optional :: unit_output
        integer, intent(in), optional :: unit_subspaces
        integer, intent(in), optional :: unit_sigmas
        integer, intent(in), optional :: unit_results
        integer, intent(in), optional :: overlap_sorting_scheme
        integer, intent(in), optional :: refresh_trial_rate
        logical, optional :: solve_values
        logical, optional :: solve_linear_sys
        logical, optional :: solve_right
        logical, optional :: solve_left
        logical, optional :: verbose
        logical, optional :: symmetric
        logical, optional :: overlap_sorting
        logical, optional :: save_sigmas 
        logical, optional :: save_subspaces
        logical, optional :: save_results
        logical, optional :: complex_mode
        logical, optional :: use_preconditioner

        call mfsolver_setup(config, & 
                              convergence_threshold, &
                              max_subspace_size, &
                              max_iterations,    &
                              unit_output,       &
                              unit_sigmas,       &
                              unit_subspaces,    &
                              unit_results,      &
                              solve_values,      &
                              solve_linear_sys,  &
                              solve_right,       &
                              solve_left,        &
                              symmetric,         &
                              complex_mode,      &
                              use_preconditioner,&
                              refresh_trial_rate,&
                              verbose,           &
                              overlap_sorting,   &
                              overlap_sorting_scheme,  &
                              save_results,      &
                              save_sigmas,       &
                              save_subspaces,    &
                              energy_shift)
    end subroutine


    subroutine eti_davidson_destroy(results)
        type(eti_davidson_results_tens_t), intent(inout) :: results
        call free_davidson_results(results)
    end subroutine


    subroutine update_trial_vectors(delta,workspace_t,flavor,Index_small_v,isv)
        type(eti_davidson_workspace_tens_t) :: workspace_t
        logical, pointer, intent(inout) :: delta(:)
        character(len=2) :: flavor
        integer,ALLOCATABLE :: Index_small_v(:)
        integer :: isv

        integer :: ierr, i, j
        integer(INTD) :: rank, rank1, rank2
        logical :: debug = .false.

        rank1=workspace_t%tens1_dims(1)
        rank2=workspace_t%tens1_dims(2)

        if (rank1.gt.rank2) then 
            rank = rank2
        else
            rank = rank1
        end if 

        workspace_t%Lprev = workspace_t%L

        do i = 1, size(delta,1)
           if (delta(i)) then
              workspace_t%L = workspace_t%L + 1
              j = workspace_t%L

              ierr=initialize_tensor_t(workspace_t%trial(j)%s1, ZERO)
              ierr=initialize_tensor_t(workspace_t%trial(j)%s2, ZERO)

              ierr=t_einsum("t(a,b)+=d(a,b)",        workspace_t%trial(j)%s1,workspace_t%delta(i)%s1)
              ierr=t_einsum("t(a,b,c,d)+=d(a,b,c,d)",workspace_t%trial(j)%s2,workspace_t%delta(i)%s2)
              
              if (debug) then
                 write (6,*) 'debug >>> correction vector ',i,' added to trial vector ',j
                 call print_tensor(workspace_t%trial(j)%s1,1.0d-8, "new b1")
                 call print_tensor(workspace_t%trial(j)%s2,1.0d-8, "new b2")
                 write (6,*) ''
              end if
           end if 
        end do

        if (debug) then
           write (6,*) 'debug >>> in update_trial_vectors, after  processing correction vectors'
           write (6,*) 'debug >>> Lprev ',workspace_t%Lprev
           write (6,*) 'debug >>> L     ',workspace_t%L
           write (6,*) ''

           write (6,*) 'debug >>> in update_trial_vectors, trial vectors before GS orthonormalization'
           do i = 1, workspace_t%L
              write (6,*) 'debug >>> trial vector ',i
              call print_tensor(workspace_t%trial(i)%s1,1.0d-8, "b1")
              call print_tensor(workspace_t%trial(i)%s2,1.0d-8, "b2")
              write (6,*) ''
           end do
        end if

        call gram_schmidt_orthogonalization_tensors(workspace_t%Lprev, workspace_t%L, &
                                                workspace_t%trial%s1, workspace_t%trial%s2,rank,flavor,Index_small_v,isv)

        if (debug) then
           write (6,*) 'debug >>> in update_trial_vectors, trial vectors after  GS orthonormalization '
           do i = 1, workspace_t%L
              write (6,*) 'debug >>> trial vector ',i
              call print_tensor(workspace_t%trial(i)%s1,1.0d-8, "b1")
              call print_tensor(workspace_t%trial(i)%s2,1.0d-8, "b2")
              if (rank.eq.1) then
                print *, "norm of new trial", sqrt(tensor_norm2(workspace_t%trial(i)%s1) + tensor_norm2(workspace_t%trial(i)%s2)/2)
              else
                print *, "norm of new trial", sqrt(tensor_norm2(workspace_t%trial(i)%s1) + tensor_norm2(workspace_t%trial(i)%s2)/4)
              end if
              write (6,*) ''
           end do
        end if
    end subroutine


    subroutine double_check_trial_vector(L,Lnew,trial1,trial2,rank)
        integer, intent(inout) :: L, Lnew
        type(tensor_t), intent(inout) :: trial1(:), trial2(:)
        integer :: rank

        integer :: Lnew_tem, L_tem
        integer :: i, j, vi
        real(8) :: norm_r
        integer, allocatable :: valid_index(:)
        type(C_PTR) :: body_p, body_p1
        integer :: ierr
        logical :: debug_double_check = .false.

        integer :: t2_dims(4), t2_rank, t1_dims(2), t1_rank
        complex(8), pointer, contiguous:: trial2_tens(:,:,:,:), trial2_tens_tem(:,:,:,:)
        complex(8), pointer, contiguous:: trial1_tens(:,:), trial1_tens_tem(:,:)


        Lnew_tem = Lnew
        L_tem = L
        allocate(valid_index(Lnew_tem-L+1))
        valid_index = 0

        if (debug_double_check) then
            print *, "valid index of vectors before checking", valid_index
        end if 

        do i = L+1, Lnew_tem
            if (rank==1) then
                norm_r = sqrt(tensor_norm2(trial1(i)) + tensor_norm2(trial2(i))/2)
            else
                norm_r = sqrt(tensor_norm2(trial1(i)) + tensor_norm2(trial2(i))/4)
            end if
            if (debug_double_check) then
                print *, "norm in double check", norm_r
            end if 
            if (norm_r.gt.config%antisymme_threshold) then 
                L_tem = L_tem + 1
                valid_index(L_tem-L) = i
            end if 
        end do
        Lnew = L_tem
       
        if (debug_double_check) then
            print *, "valid index of vectors after checking", valid_index
        end if 

        j = 1
        ierr = get_nmodes(trial1(1),t1_rank)
        ierr = get_extents_static(trial1(1),t1_dims)

        ierr = get_nmodes(trial2(1),t2_rank)
        ierr = get_extents_static(trial2(1),t2_dims)

        do i = L+1, Lnew
            vi = valid_index(j)

            if (debug_double_check) then
                print *, "valid index vi:", vi
                print *, "original index i:", i
                call print_tensor(trial1(i),1.0d-18, "trail1(i) in double check before")
                call print_tensor(trial2(i),1.0d-18, "trail2(i) in double check before")

                call print_tensor(trial1(vi),1.0d-18, "trail1(vi) in double check before")
                call print_tensor(trial2(vi),1.0d-18, "trail2(vi) in double check before")
            end if 

            ierr=direct_data_pointer(trial2(vi),body_p)
            ierr=direct_data_pointer(trial2(i),body_p1)
            ! ierr=talsh_init(tensor_position_you_waht_to_overwrite,ZERO)
            ! ierr=talsh_add(replacement+=vector_you_keep)
            call c_f_pointer(body_p,trial2_tens,t2_dims)
            call c_f_pointer(body_p1,trial2_tens_tem,t2_dims)
            trial2_tens_tem = trial2_tens
            
            ierr=direct_data_pointer(trial1(vi),body_p)
            ierr=direct_data_pointer(trial1(i),body_p1)
            call c_f_pointer(body_p,trial1_tens,t1_dims)
            call c_f_pointer(body_p1,trial1_tens_tem,t1_dims)
            trial1_tens_tem = trial1_tens
            
            j = j+1

            if (debug_double_check) then
                call print_tensor(trial1(i),1.0d-18, "trail1(i) in double check after")
                call print_tensor(trial2(i),1.0d-18, "trail2(i) in double check after")
            end if 

        end do
        
        if (Lnew > L) then
            L = L+1
        end if
        print *, "Lnew:", Lnew
        print *, "L:", L
        deallocate(valid_index)

    end subroutine

    subroutine renormalize_vectors(v1,v2,rank)

       type(tensor_t), intent(inout) :: v1, v2
       integer :: rank
       real :: norm_r
       complex(8) :: norm
       integer :: ierr
  

       if (rank.eq.1) then
            norm_r = sqrt(tensor_norm2(v1) + tensor_norm2(v2)/2)
       else 
            norm_r = sqrt(tensor_norm2(v1) + tensor_norm2(v2)/4)
       end if 

    !    norm_r = sqrt(tensor_norm2(v1) + tensor_norm2(v2)/4) !XY 4 is for EE(ab.ij) but we need to take care for IP
       norm = dcmplx(norm_r, 0.0d0)

!       if (norm_r.lt.config%antisymme_threshold) then  ! to avoid numerical noise
!            ierr=initialize_tensor_t(v1)
!            ierr=initialize_tensor_t(v2)
!       else
!            call divide_tensor_elements_by_value(v1, norm)
!            call divide_tensor_elements_by_value(v2, norm)
!       end if

        if (norm_r .ge. config%safe_minimum ) then
            call divide_tensor_elements_by_value(v1, norm)
            call divide_tensor_elements_by_value(v2, norm)
        else
            call quit ('norm of tensors too small for division, quitting')
        end if

    end subroutine


    subroutine gram_schmidt_projection_tensors(factor, u1, u2, v1, v2, rank)
       complex(8), intent(inout) :: factor
       type(tensor_t), intent(inout) :: u1, u2
       type(tensor_t), intent(inout) :: v1, v2
       integer, intent(in) :: rank


       type(tensor_t) :: result_tensor_num, result_tensor_den
       complex(8), pointer :: result_tens_num(:), result_tens_den(:)
       type(C_PTR):: body_p_num, body_p_den
       type(tensor_t) :: one_tensor

       integer(C_INT)     :: one_dims(1)
       integer  :: ierr
       logical :: debug_GS_proj = .false.
       complex(8)    :: contrac_coeff

       one_dims(1) = 1

       ierr=construct_tensor_t(result_tensor_den,C8,one_dims(1:0),init_val=ZERO)
       ierr=direct_data_pointer(result_tensor_den,body_p_den)
       call c_f_pointer(body_p_den,result_tens_den,one_dims)
       ierr=construct_tensor_t(result_tensor_num,C8,one_dims(1:0),init_val=ZERO)
       ierr=direct_data_pointer(result_tensor_num,body_p_num)
       call c_f_pointer(body_p_num,result_tens_num,one_dims)

       ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)

    ! aspg: verify this and generalize to other ranks !       
       if (rank == 1) then
        contrac_coeff = ONE_HALF
       else ! generalize for other ranks !
        contrac_coeff = ONE_QUARTER
       end if 

       if (debug_GS_proj) then
           call print_tensor(u1, 1.0d-8, "u1 tensor")
           call print_tensor(v1, 1.0d-8, "v1 tensor")
       end if 
       ierr=t_einsum("num()+=u+(a,b)*v(a,b)",result_tensor_num,u1,v1)
       ierr=t_einsum("den()+=u+(a,b)*u(a,b)",result_tensor_den,u1,u1)
       ierr=t_einsum("num()+=u+(a,b,c,d)*v(a,b,c,d)",result_tensor_num,u2,v2,scale=contrac_coeff)
       ierr=t_einsum("den()+=u+(a,b,c,d)*u(a,b,c,d)",result_tensor_den,u2,u2,scale=contrac_coeff)

       factor = result_tens_num(1) / result_tens_den(1)

       if (debug_GS_proj) then
          print *, "result_tens_den():", result_tens_den(1) 
          print *, "result_tens_num():", result_tens_num(1)
          print *, "factor: ", factor
       end if 

       ierr=destroy_tensor_t(result_tensor_num)
       ierr=destroy_tensor_t(result_tensor_den)
       ierr=destroy_tensor_t(one_tensor)

    end subroutine


    subroutine gram_schmidt_orthogonalization_tensors(L,Lnew,trial1,trial2,rank,flavor,Index_small_v,isv)
        integer, intent(inout) :: L, Lnew
        type(tensor_t), intent(inout) :: trial1(:), trial2(:)
        integer(INTD) :: rank
        character(len=2) :: flavor
        integer,ALLOCATABLE :: Index_small_v(:)
        integer :: isv

        real :: norm_r
        integer :: ierr, i, j, n, m
        complex(8) :: factor
        type(tensor_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        logical :: debug_GS_orthg = .false.
        logical :: found

        n = size(Index_small_v)

        do j = L+1, Lnew
            if (debug_GS_orthg) then
                call print_tensor(trial1(j), 1.0d-18, "trial1(j) before ")
                call print_tensor(trial2(j), 1.0d-18, "trial2(j) before ")
            end if
            

            call antisymm_t2(trial2(j),flavor)
            
            do i = 1, j-1
                if (debug_GS_orthg) then
                    print *, "I: ", i, "J: ", j
                    call print_tensor(trial1(i), 1.0d-10, "trial1(i)")
                    call print_tensor(trial2(i), 1.0d-10, "trial2(i)")
                end if 
                
                found = is_number_in_array(i, Index_small_v, n)
                if (.not.found) then
                    call gram_schmidt_projection_tensors(factor,trial1(i), trial2(i), trial1(j), trial2(j),rank)
                end if

                ierr=t_einsum("t(a,b)+=p(a,b)",trial1(j),trial1(i),scale=-factor)
                ierr=t_einsum("t(a,b,c,d)+=p(a,b,c,d)",trial2(j),trial2(i),scale=-factor)

            end do


            do i = 1, j-1            
                found = is_number_in_array(i, Index_small_v, n) 
                if (.not.found) then     
                    call gram_schmidt_projection_tensors(factor,trial1(i), trial2(i), trial1(j), trial2(j),rank)
                end if
                ierr=t_einsum("t(a,b)+=p(a,b)",trial1(j),trial1(i),scale=-factor)
                ierr=t_einsum("t(a,b,c,d)+=p(a,b,c,d)",trial2(j),trial2(i),scale=-factor)    
            end do

            if (rank.eq.1) then
                norm_r = sqrt(tensor_norm2(trial1(j))+tensor_norm2(trial2(j))/2)  
            else 
                norm_r = sqrt(tensor_norm2(trial1(j))+tensor_norm2(trial2(j))/4)
            end if

            if (norm_r.lt.config%antisymme_threshold) then 
                ierr=destroy_tensor_t(trial1(j))
                ierr=destroy_tensor_t(trial2(j))
                isv = isv+1
                Index_small_v(isv) = j

            else
                call renormalize_vectors(trial1(j),trial2(j),rank)
            end if 

        end do
        L = L+1


    end subroutine


    function is_number_in_array(a, V, n) result(found)
        integer, intent(in) :: a, n
        integer, intent(in) :: V(n)
        logical :: found
        integer :: i

        found = .false.

        do i = 1, n
            if (V(i) == a) then
                found = .true.
                exit
            end if
        end do
    end function is_number_in_array


    function which_number_in_array(a, V, n) result(i)
        integer, intent(in) :: a, n
        integer, intent(in) :: V(n)
        integer :: i

        do i = 1, n
            if (V(i) == a) then
                exit
            end if
        end do
    end function which_number_in_array


    subroutine antisymm_t2(trial2,flavor)

        type(tensor_t), intent(inout) :: trial2
        character(len=2) :: flavor
    
        integer :: i, j, a, b
        integer :: N_1, N_2, N_3, N_4
        real(8) :: norm_r
        integer, allocatable :: valid_index(:)
        type(C_PTR) :: body_p, body_p1
        integer :: ierr
        logical :: debug_double_check = .false.

        integer :: t2_dims(4), t2_rank
        complex(8), pointer, contiguous:: trial2_tens(:,:,:,:), trial2_tens_tem(:,:,:,:)


        ierr = get_nmodes(trial2,t2_rank)
        ierr = get_extents_static(trial2,t2_dims)

        if (debug_double_check) then
            call print_tensor(trial2,1.0d-12, "trial2 in double check before")
        end if 

        ierr=direct_data_pointer(trial2,body_p)
        call c_f_pointer(body_p,trial2_tens,t2_dims)
        N_1 = size(trial2_tens,1)
        N_2 = size(trial2_tens,2)
        N_3 = size(trial2_tens,3)
        N_4 = size(trial2_tens,4)
        allocate(trial2_tens_tem(N_1,N_2,N_3,N_4))
        trial2_tens_tem = trial2_tens
        

        if (flavor=='IP') then

            do j=1, size(trial2_tens_tem,4)
                do i=j, size(trial2_tens_tem,3)
                    do a=1, size(trial2_tens_tem,1)
                            if (i==j) then 
                                trial2_tens(a,1,i,j) = dcmplx(0.0,0.0)
                            else
                                trial2_tens(a,1,i,j) = trial2_tens_tem(a,1,i,j)
                            end if

                            trial2_tens(a,1,j,i) = -trial2_tens(a,1,i,j)
                    end do
                end do
            end do

        else if (flavor=='EA') then
            do i=1, size(trial2_tens_tem,3)
                do b=1, size(trial2_tens_tem,2)
                    do a=b, size(trial2_tens_tem,1)
                            if (a==b) then 
                                trial2_tens(a,b,i,1) = dcmplx(0.0,0.0)
                            else
                                trial2_tens(a,b,i,1) = trial2_tens_tem(a,b,i,1)
                            end if

                            trial2_tens(b,a,i,1) = -trial2_tens(a,b,i,1)
                    end do
                end do
            end do

        else

            do j=1, size(trial2_tens_tem,4)
                do i=j, size(trial2_tens_tem,3)
                    do b=1, size(trial2_tens_tem,2)
                        do a=b, size(trial2_tens_tem,1)
                            if (i==j.or.a==b) then 
                                trial2_tens(a,b,i,j) = dcmplx(0.0,0.0)
                            else
                                trial2_tens(a,b,i,j) = trial2_tens_tem(a,b,i,j)
                            end if

                            trial2_tens(b,a,i,j) = -trial2_tens(a,b,i,j)
                            trial2_tens(a,b,j,i) = -trial2_tens(a,b,i,j)
                            trial2_tens(b,a,j,i) = trial2_tens(a,b,i,j)
                    end do
                end do
            end do
        end do

        end if 
        

        if (debug_double_check) then
            call print_tensor(trial2,1.0d-12, "trial2 in double check after")
        end if 

        
        deallocate(trial2_tens_tem)
    end subroutine


    subroutine davidson_print_vectors_complex(vectors, label, nvectors)
        complex(kind=8), intent(in) :: vectors(:,:)
        integer, intent(in), optional :: nvectors
        character(len=max_label_length), intent(in), optional :: label

        integer :: L, N, i, j
        real(kind=8) :: vr, vi
        real(kind=8) :: print_vector_threshold = 1.0d-10
        character(len=max_label_length) :: label_local

        N = size(vectors,1)
        L = size(vectors,2)
        if (present(nvectors)) L = nvectors

        if (present(label)) then
            label_local = label
            print *, ""
            print "(a)", label_local(1:max_label_length)
            print *, ""
            print *, ""
        end if

        do i = 1, N
            print "(2x,i6,$)", i
            do j = 1, L
                vr = real(vectors(i,j))
                vi = aimag(vectors(i,j))
                if (abs(vi).lt.print_vector_threshold) vi = 0.0d0
                if (abs(vr).lt.print_vector_threshold) vr = 0.0d0
                print '(4x,1e14.6,2x,SP,1e14.6,"i",$)', dcmplx(vr,vi)
            end do
            print *, ""
        end do
    end subroutine 


    subroutine davidson_print_vectors_real(vectors, label, nvectors, complex_dimension)
        real(kind=8), intent(in) :: vectors(:,:)
        integer, intent(in), optional :: nvectors
        character(len=max_label_length), intent(in), optional :: label 
        integer, intent(in), optional :: complex_dimension

        integer :: L, N, i, j, N_complex, L_complex, N_step, L_step, kr, ki
        real(kind=8) :: vr, vi
! aspg, for debugging, set this threshold to eliminate small values in printout
!       for both real (and imaginary) parts
        real(kind=8) :: print_vector_threshold = 1.0d-10 
        character(len=max_label_length) :: label_local

        N = size(vectors,1)
        L = size(vectors,2)
        if (present(nvectors)) L = nvectors

        if (present(label)) then
            label_local = label
            print *, ""
            print "(a)", label_local(1:max_label_length)
            print *, ""
            print *, ""
        end if
        if (present(complex_dimension)) then
            if (complex_dimension.eq.1) then
                N_step = 2
                L_step = 1
            else if (complex_dimension.eq.2) then
                N_step = 1
                L_step = 2
            end if
        else
           N_step = 2
           L_step = 1
        end if

        if (config%complex_mode_with_reals) then
            N_complex = N / N_step
            L_complex = L / L_step

            do i = 1, N_complex
                print "(2x,i6,$)", i
                do j = 1, L_complex
                    if (N_step.eq.2) then
                        kr = N_step*(i-1) + 1
                        ki = kr + 1
                        vr = vectors(kr,j)
                        vi = vectors(ki,j)
                    else 
                        kr = L_step*(j-1) + 1
                        ki = kr + 1
                        vr = vectors(i, kr)
                        vi = vectors(i, ki)
                    end if
                    if (abs(vi).lt.print_vector_threshold) vi = 0.0d0
                    if (abs(vr).lt.print_vector_threshold) vr = 0.0d0
                    if (vi < 0.0d0) then
                        print "(4x,1e14.6,a,1e13.6,a,$)", vr," -",abs(vi)," i"
                    else
                        print "(4x,1e14.6,a,1e13.6,a,$)", vr," +",abs(vi)," i"
                    end if
                end do
                print *, ""
            end do
        end if 
    end subroutine


#ifdef NOT_USED
    subroutine davidson_print_vectors_complex_1d(vectors, N, L, label, nvectors)
        complex(kind=8), intent(in) :: vectors(:)
        integer, intent(in) :: N, L
        integer, intent(in), optional :: nvectors
        character(len=max_label_length), intent(in), optional :: label

        integer :: i, j, k, kr, ki
        real(kind=8) :: vr, vi
! aspg, for debugging, set this threshold to eliminate small values in printout
!       for both real (and imaginary) parts
        real(kind=8) :: print_vector_threshold = 1.0d-10
        character(len=max_label_length) :: label_local

        if (present(label)) then
            label_local = label
            print *, ""
            print "(a)", label_local(1:max_label_length)
            print *, ""
            print *, ""
        end if

        do i = 1, N*L
           print "(2x,i6,$)", i
!          do j = 1, L
!              
!              k  = N*(i-1) + L
!              print "(8x,1e14.6,1e14.6)", vectors(l)
               k = i
               vi = real(vectors(k))
               vr = aimag(vectors(k)) 
               if (abs(vi).lt.print_vector_threshold) vi = 0.0d0
               if (abs(vr).lt.print_vector_threshold) vr = 0.0d0
               if (vi < 0.0d0) then
                   print "(4x,1e14.6,a,1e13.6,a,$)", vr," -",abs(vi)," i"
               else
                   print "(4x,1e14.6,a,1e13.6,a,$)", vr," +",vi," i"
               end if
!          end do
           print *, ""
        end do
        
    end subroutine
#endif



    subroutine print_values_vectors_C(output_file_unit, eValues, eVectorsR, eVectorsL, left, right, nroots)
        complex(kind=8), intent(in) :: eVectorsL(:,:)
        complex(kind=8), intent(in) :: eVectorsR(:,:)
        complex(kind=8), intent(in) :: eValues(:)
        integer, optional, intent(in) :: nroots
        integer, intent(in) :: output_file_unit

        logical :: left, right
        integer :: vectors_N
        integer :: nroots_local
        integer :: i
        character(len=max_label_length) :: label
        real(kind=8) :: vr, vi
! aspg, for debugging, set this threshold to eliminate small values in printout
!       for both real (and imaginary) parts
        real(kind=8) :: print_values_threshold = 1.0d-10


        vectors_N = size(eValues,1)
        nroots_local = vectors_N
        if (present(nroots)) nroots_local = nroots

        label = "Eigenvalues"

        do i = 1, nroots_local 
            vr = real(eValues(i))
            vi = aimag(eValues(i))
            if (abs(vi).lt.print_values_threshold) vi = 0.0d0
            if (abs(vr).lt.print_values_threshold) vr = 0.0d0

            if (vi .lt. 0.0d0) then
                write (config%output_file_unit,'(4x,a,i6,a,1e20.12,a,1e20.12,a)')  &
                      "Eigenvalue ",i," : ", vr," -",abs(vi)," i"
            else if (vi .gt. 0.0d0) then
                write (config%output_file_unit,'(4x,a,i6,a,1e20.12,a,1e20.12,a)')  &
                      "Eigenvalue ",i," : ", vr," +",abs(vi)," i"
            else
                write (config%output_file_unit,'(4x,a,i6,a,1e20.12)')  &
                      "Eigenvalue ",i," : ", vr
            end if
        end do
        if (left) then
            label = "Left Eigenvectors            "
            call davidson_print_vectors(eVectorsL, label=label) 
        end if
        if (right) then
            label = "Right Eigenvectors           "
            call davidson_print_vectors(eVectorsR, label=label) 
        end if

    end subroutine 

    subroutine print_values_vectors_R(output_file_unit, eValues, eVectorsR, eVectorsL, left, right, nroots)
        real(kind=8), intent(in) :: eVectorsL(:,:)
        real(kind=8), intent(in) :: eVectorsR(:,:)
        real(kind=8), intent(in) :: eValues(:,:)
        integer, optional, intent(in) :: nroots
        integer, intent(in) :: output_file_unit

        logical :: left, right
        integer :: vectors_N
        integer :: nroots_local
        integer :: i
        character(len=max_label_length) :: label
        real(kind=8) :: vr, vi
! aspg, for debugging, set this threshold to eliminate small values in printout
!       for both real (and imaginary) parts
        real(kind=8) :: print_values_threshold = 1.0d-10


        vectors_N = size(eValues,1)
        nroots_local = vectors_N
        if (present(nroots)) nroots_local = nroots

        label = "Eigenvalues"
        if (config%debug) print *, "debug, evalues inside print_values_vectors_R"
        if (config%debug) print *, eValues(1,:)
        if (config%debug) print *, eValues(2,:)

        do i = 1, nroots_local
            vr = eValues(1,i)
            vi = eValues(2,i)
            if (abs(vi).lt.print_values_threshold) vi = 0.0d0
            if (abs(vr).lt.print_values_threshold) vr = 0.0d0

            if (vi .lt. 0.0d0) then
                write (config%output_file_unit,'(4x,a,i6,a,1e14.6,a,1e14.6,a)')  &
                      "Eigenvalue ",i," : ", vr," -",abs(vi)," i"
            else if (vi .gt. 0.0d0) then
                write (config%output_file_unit,'(4x,a,i6,a,1e14.6,a,1e14.6,a)')  &
                      "Eigenvalue ",i," : ", vr," +",abs(vi)," i"
            else
                write (config%output_file_unit,'(4x,a,i6,a,1e14.6)')  &
                      "Eigenvalue ",i," : ", vr
            end if
        end do
        if (left) then
            label = "Left Eigenvectors            "
            call davidson_print_vectors(eVectorsL, label=label)
        end if
        if (right) then
            label = "Right Eigenvectors           "
            call davidson_print_vectors(eVectorsR, label=label)
        end if

    end subroutine

    function vector_norm(v)
        real(kind=8) :: vector_norm
        real(kind=8), intent(in) :: v(:)
        integer :: vector_N, i, vector_N_complex
        double precision dnrm2, dznrm2

        vector_N = size(v,1)
        if (config%complex_mode_with_reals) then
            vector_N_complex = vector_N / 2
            vector_norm = dznrm2(vector_N_complex, v, 1)
        end if
    end function

    subroutine form_r_tensors(k,eigenvectors_subspace,eigenvalues_subspace,workspace_t,Index_small_v,prop_grad)
! we use delta{1,2,...} to store the values of r^k
        integer, intent(in) :: k
        complex(kind=8), intent(in)  :: eigenvectors_subspace(:,:)
        complex(kind=8), intent(in)  :: eigenvalues_subspace(:)
        type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
        integer,ALLOCATABLE :: Index_small_v(:)
! aspg, added from (q)rsp
        ! complex(kind=8), pointer, contiguous, optional, intent(in) :: C_subspace(:)
!       type(eti_sa_tens_t),intent(inout), optional :: prop_grad
        type(eti_comm_tens_t),intent(inout), optional :: prop_grad
! aspg, added from (q)rsp
        complex(kind=8) :: factor1, factor2, factor3
        integer(INTD) :: r_rank(5)
        integer :: ierr, i, i_pri
        logical :: debug_form_r = .false.
        logical :: found 
        type(tensor_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        one_dims(1) = 1
        ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
! we only need to know the rank of delta1 because the rank of r2 will always be r1_rank+2, that is: 
! r1 = 1 (1h for IP, 1p for EA) then r2 = 3 (2h1p for IP, 1h2p for EA)
! r1 = 2 (1h1p for EE, 2h for DIP, 2p for DEA) then r2 = 4 (2h2p for EE, 3h1p for DIP, 1h3p for DEA)
! etc
        ierr = get_nmodes(workspace_t%delta(k)%s1,r_rank(1))
        ierr = get_nmodes(workspace_t%delta(k)%s2,r_rank(2))

        if (r_rank(1).gt.2) call quit('rank not implemented in form_r_tensors')

!       ierr=initialize_tensor_t(workspace_t%delta(k)%s1, ZERO)
!       ierr=initialize_tensor_t(workspace_t%delta(k)%s2, ZERO)

!       if (debug_form_r) then
!          call davidson_print_vectors(eigenvectors_subspace,   label="eigenvectors_subspace in form_r")
!       end if
        
        i_pri = 0
        do i = 1, workspace_t%L
          found = is_number_in_array(i,Index_small_v,size(Index_small_v))
          if (.not.found) then

           i_pri = i_pri+1
        !    factor1 =  eigenvectors_subspace(i,k)
           factor1 = eigenvectors_subspace(i_pri,k)
        !    factor1 = conjg(factor1)
           factor2 = -eigenvalues_subspace(k)*factor1
! aspg, added from (q)rsp
        !    if (present(C_subspace)) factor2 = -C_subspace(i)
        !    if (present(C_subspace)) factor3 = -C_subspace(i)
! aspg, added from (q)rsp

           if (debug_form_r) write(6,*) 'c(',i,',',k,')',factor1,'lambda*c',factor2

           if (r_rank(1).eq.1) then
              
            !   call rotate_eigenvector(workspace_t%sigma(i)%s1,workspace_t%sigma(i)%s2,workspace_t)
              ierr=t_einsum("r(a)+=s(a)*K()",workspace_t%delta(k)%s1,workspace_t%sigma(i)%s1,one_tensor,scale=factor1)
              ierr=t_einsum("r(a)+=t(a)*K()",workspace_t%delta(k)%s1,workspace_t%trial(i)%s1,one_tensor,scale=factor2)
              ierr=t_einsum("r(a,b,c)+=s(a,b,c)*K()",workspace_t%delta(k)%s2,workspace_t%sigma(i)%s2,&
                                & one_tensor,scale=factor1)
              ierr=t_einsum("r(a,b,c)+=t(a,b,c)*K()",workspace_t%delta(k)%s2,workspace_t%trial(i)%s2,&
                                & one_tensor,scale=factor2)

            else if (r_rank(1).eq.2) then
              ierr=t_einsum("r(a,b)+=s(a,b)",workspace_t%delta(k)%s1,workspace_t%sigma(i)%s1,scale=factor1)
              ierr=t_einsum("r(a,b)+=t(a,b)",workspace_t%delta(k)%s1,workspace_t%trial(i)%s1,scale=factor2)
              ierr=t_einsum("r(a,b,c,d)+=s(a,b,c,d)",workspace_t%delta(k)%s2,workspace_t%sigma(i)%s2,scale=factor1)
              ierr=t_einsum("r(a,b,c,d)+=t(a,b,c,d)",workspace_t%delta(k)%s2,workspace_t%trial(i)%s2,scale=factor2)

            !   call print_tensor(workspace_t%delta(k)%s1, 1.0d-10, "delta2")

            !   if (present(C_subspace)) then
                ! ierr=t_einsum("r(a,b)+=t(a,b)",workspace_t%delta(k)%s1,workspace_t%trial(i)%s1,scale=factor3)
                ! ierr=t_einsum("r(a,b,c,d)+=t(a,b,c,d)",workspace_t%delta(k)%s2,workspace_t%trial(i)%s2,scale=factor3)
            !   end if 

           end if

         if (debug_form_r) then
            write (6,*) 'debug >>> in form_r_tensors, accumulated r1 root',k,' for L',i,' >>>'
            call print_tensor_elements(workspace_t%delta(k)%s1)
            write (6,*) 'debug <<< in form_r_tensors, accumulated r1 for L',i,' <<<'
            write (6,*) 'debug >>> in form_r_tensors, accumulated r2 root',k,' for L',i,' >>>'
            call print_tensor_elements(workspace_t%delta(k)%s2)
            write (6,*) 'debug <<< in form_r_tensors, accumulated r2 for L',i,' <<<'
         end if
        end if
           
       end do

       if (present(prop_grad)) then
        ierr=t_einsum("r(a,b)+=p(a,b)*K()",workspace_t%delta(k)%s1,prop_grad%t1,&
        & one_tensor,scale=ONE)
        ierr=t_einsum("r(a,b,c,d)+=p(a,b,c,d)*K()",workspace_t%delta(k)%s2,prop_grad%t2,&
        & one_tensor,scale=ONE)
      end if 

       if (debug_form_r) then
          write (6,*) 'debug >>> in form_r_tensors, r1 root',k,' >>>'
          call print_tensor_elements(workspace_t%delta(k)%s1)
          write (6,*) 'debug <<< in form_r_tensors, r1 <<<'
          write (6,*) 'debug >>> in form_r_tensors, r2 root',k,' >>>'
          call print_tensor_elements(workspace_t%delta(k)%s2)
          write (6,*) 'debug <<< in form_r_tensors, r1 <<<'
       end if

       ierr=destroy_tensor_t(one_tensor)

    end subroutine    



    function correction_to_vectors(output_file_unit, e_small_M, c_small_M, workspace_t, Index_small_v,&
                                    prop_grad, correct_all)
        integer, intent(in)       :: output_file_unit
        type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
        complex(kind=8), intent(in)  :: e_small_M(:)
        complex(kind=8), intent(in)  :: c_small_M(:,:)
        integer,ALLOCATABLE :: Index_small_v(:)
! aspg, added from (q)rsp
        ! complex(kind=8), pointer, contiguous, optional, intent(in) ::  C_small(:)
!       type(eti_sa_tens_t),intent(inout), optional :: prop_grad
        type(eti_comm_tens_t),intent(inout), optional :: prop_grad
! aspg, added from (q)rsp
        logical, optional, intent(in) :: correct_all

        logical, pointer :: correction_to_vectors(:)

        integer :: ierr, k, unconverged_roots, nroots
        integer :: rank, rank1, rank2
        real(kind=8) :: norm_r = 0.0d0
        logical :: correct_all_local
        logical :: debug_correction_vectors = .false.

        if (present(correct_all)) then
            correct_all_local = correct_all
        else
            correct_all_local = .false. 
        end if

        ! if (present(C_small)) then
        if (present(prop_grad)) then
            nroots = 1
        else
            nroots = workspace_t%nroots
        end if

        allocate(correction_to_vectors(nroots))
        correction_to_vectors = .true.

        unconverged_roots = 0
        ! err = get_nmodes(workspace_t%delta%s1(1),rank)
        rank1 = workspace_t%tens1_dims(1) 
        rank2 = workspace_t%tens1_dims(2) 
        if (rank1.gt.rank2) then 
            rank = rank2
        else 
            rank = rank1
        end if 

        do k = 1, nroots
           ierr=initialize_tensor_t(workspace_t%delta(k)%s1, ZERO)
           if (debug_correction_vectors) print *,'debug, init workspace_t%delta(k)%s1 to zero, ierr',ierr
           ierr=initialize_tensor_t(workspace_t%delta(k)%s2, ZERO)
           if (debug_correction_vectors) print *,'debug, init workspace_t%delta(k)%s2 to zero, ierr',ierr
        end do
        do k = 1, nroots
! instead of having a new variable to store r(:,k), we will do collect the results on workspace_t%delta%s{1,2,...} 
        !    if (present(C_small)) then
           if (present(prop_grad)) then
            !   call form_r_tensors(k,c_small_M,e_small_M,workspace_t,C_small)
              call form_r_tensors(k,c_small_M,e_small_M,workspace_t,Index_small_v,prop_grad)
           else
              call form_r_tensors(k,c_small_M,e_small_M,workspace_t,Index_small_v)
           end if
            
           if (rank .eq.1) then
              norm_r = sqrt(tensor_norm2(workspace_t%delta(k)%s1) +  tensor_norm2(workspace_t%delta(k)%s2)/2)
           else 
              norm_r = sqrt(tensor_norm2(workspace_t%delta(k)%s1) +  tensor_norm2(workspace_t%delta(k)%s2)/4)
           end if

! get the norm of the whole intermediate r(:,k) tensor, which means to calculate it over the different excitation ranks 
           if (debug_correction_vectors) then
              print *,'debug, norm of delta1(',k,'):',tensor_norm2(workspace_t%delta(k)%s1)
              print *,'debug, norm of delta2(',k,'):',tensor_norm2(workspace_t%delta(k)%s2)
           end if

           if (norm_r .lt. config%convergence_threshold) then 
              write(output_file_unit, "(4x,a,i6,a,1e14.6,4x,a)") "Root ",k,"  ||q_k,M|| ", norm_r, " <<< converged ! "
              correction_to_vectors(k) = .false.
           else
              write(output_file_unit, "(4x,a,i6,a,1e14.6)")  "Root ",k,"  ||q_k,M|| ", norm_r
              correction_to_vectors(k) = .true.
              unconverged_roots = unconverged_roots + 1
           end if
        end do
        write (config%output_file_unit,*) " "

        if (any(correction_to_vectors)) then
           if (correct_all_local) then
               unconverged_roots = nroots
               correction_to_vectors = .true. 
           end if
! forms delta (the correction vectors which will serve to create new trial vectors
! note here that, as seen above, workspace_t%delta%s{1,2,..} have been used to accumulate r(i,k) 
           do k = 1, nroots
              if (correction_to_vectors(k)) then
                 if (debug_correction_vectors) then
                    write (6,*) 'debug >>> in correction_to_vectors, precond1 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%precond(k)%s1)
                    call print_tensor(workspace_t%precond(k)%s1, 1.0d-10, "workspace_t%precond(k)%s1")
                    write (6,*) 'debug <<< in correction_to_vectors, precond1 <<<'
                    write (6,*) 'debug >>> in correction_to_vectors, delta1 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta(k)%s1)
                    call print_tensor(workspace_t%delta(k)%s1, 1.0d-10, "workspace_t%delta(k)%s1")
                    write (6,*) 'debug <<< in correction_to_vectors, delta1 <<<'
                 end if

                 call multiply_tensor_elements(workspace_t%delta(k)%s1, workspace_t%precond(k)%s1)

                 if (debug_correction_vectors) then
                    write (6,*) 'debug >>> in correction_to_vectors, result delta1 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta(k)%s1)
                    call print_tensor(workspace_t%delta(k)%s1, 1.0d-10, "workspace_t%delta(k)%s1")
                    write (6,*) 'debug <<< in correction_to_vectors, result delta1 <<<'

                    write (6,*) 'debug >>> in correction_to_vectors, precond2 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%precond(k)%s2) 
                    call print_tensor(workspace_t%precond(k)%s2, 1.0d-10, "workspace_t%precond%s2(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, precond2 <<<'
                    write (6,*) 'debug >>> in correction_to_vectors, delta2 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta(k)%s2)
                    call print_tensor(workspace_t%delta(k)%s2, 1.0d-10, "workspace_t%delta(k)%s2")
                    write (6,*) 'debug <<< in correction_to_vectors, delta2 <<<'
                 end if

                !  call print_tensor(workspace_t%delta%s2(1), 1.0d-10, "delta2 before precond")
                !  call print_tensor(workspace_t%precond%s2(1), 1.0d-10, "precond2 before precond")
                 call multiply_tensor_elements(workspace_t%delta(k)%s2, workspace_t%precond(k)%s2)
                !  call print_tensor(workspace_t%delta%s2(1), 1.0d-10, "delta2 after precond")
                !  call print_tensor(workspace_t%precond%s2(1), 1.0d-10, "precond2 after precond")

                 if (debug_correction_vectors) then
                    write (6,*) 'debug >>> in correction_to_vectors, result delta2 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta(k)%s2)
                    call print_tensor(workspace_t%delta(k)%s2, 1.0d-10, "workspace_t%delta(k)%s2")
                    write (6,*) 'debug <<< in correction_to_vectors, result delta2 <<<'
                 end if

                !  norm_r = tensor_norm2(workspace_t%delta(k)%s1) +  tensor_norm2(workspace_t%delta(k)%s2)
                !  print *,'debug, norm of delta+delta2():',norm_r

              else
                 ierr=initialize_tensor_t(workspace_t%delta(k)%s1, ZERO)
                 ierr=initialize_tensor_t(workspace_t%delta(k)%s2, ZERO)
              end if

              if (debug_correction_vectors) then
                 write (6,*) 'debug >>> in correction_to_vectors, delta1 root',k,' >>>'
                !  call print_tensor_elements(workspace_t%delta(k)%s1)
                 call print_tensor(workspace_t%delta(k)%s1, 1.0d-10, "workspace_t%delta(k)%s1")
                 write (6,*) 'debug <<< in correction_to_vectors, delta1 <<<'
                 write (6,*) 'debug >>> in correction_to_vectors, delta2 root',k,' >>>'
                !  call print_tensor_elements(workspace_t%delta(k)%s2)
                 call print_tensor(workspace_t%delta(k)%s2, 1.0d-10, "workspace_t%delta(k)%s2")
                 write (6,*) 'debug <<< in correction_to_vectors, delta2 <<<'
              end if

           end do
        else
            deallocate(correction_to_vectors)
            correction_to_vectors => NULL()
        end if

    end function

    subroutine form_preconditioner(workspace_t, eigenvalues)
        type(eti_davidson_workspace_tens_t) :: workspace_t
        complex(kind=8), intent(in)  :: eigenvalues(:)

        complex(8), pointer, contiguous:: preco_1d_p(:)
        complex(8), pointer, contiguous:: preco_2d_p(:,:) 
        complex(8), pointer, contiguous:: preco_3d_p(:,:,:) 
        complex(8), pointer, contiguous:: preco_4d_p(:,:,:,:) 
        complex(8), pointer, contiguous:: hdiag_1d_p(:)
        complex(8), pointer, contiguous:: hdiag_2d_p(:,:) 
        complex(8), pointer, contiguous:: hdiag_3d_p(:,:,:) 
        complex(8), pointer, contiguous:: hdiag_4d_p(:,:,:,:) 
        type(C_PTR):: preco1_body_p, preco2_body_p
        type(C_PTR):: hdiag1_body_p, hdiag2_body_p

        integer :: iroot, ja, jb, jc, jd, ierr, b
        integer :: nroots
        real(kind=8) :: x2py2, proj_r, proj_i, eigenvalue_r, eigenvalue_i 
        logical :: preconditioner_debug = .false.
        real(kind=8) :: small_number = 1.0d-8 !LV changed to 1.0d-8

! careful with dimension 1 of eigenvalues, 
        nroots    = workspace_t%nroots

! whether or not we want to precondition, we start by setting up the tensor to one
        do iroot = 1, nroots
           ierr=initialize_tensor_t(workspace_t%precond(iroot)%s1, ONE)   
           ierr=initialize_tensor_t(workspace_t%precond(iroot)%s2, ONE)
        end do

        if (config%use_preconditioner) then
! we want to precondition
! first prepare to access the rank-1 and rank-2 diagonal of Hbar 
           ierr = direct_data_pointer(workspace_t%Hdiag%s1,hdiag1_body_p)
           ierr = direct_data_pointer(workspace_t%Hdiag%s2,hdiag2_body_p)

           if (workspace_t%tens1_rank .eq. 1) then
              call c_f_pointer(hdiag1_body_p,hdiag_1d_p,workspace_t%tens1_dims(1:workspace_t%tens1_rank))
              call c_f_pointer(hdiag2_body_p,hdiag_3d_p,workspace_t%tens2_dims(1:workspace_t%tens2_rank))
           else if (workspace_t%tens1_rank .eq. 2) then
              call c_f_pointer(hdiag1_body_p,hdiag_2d_p,workspace_t%tens1_dims(1:workspace_t%tens1_rank))
              call c_f_pointer(hdiag2_body_p,hdiag_4d_p,workspace_t%tens2_dims(1:workspace_t%tens2_rank))
           end if

           if (preconditioner_debug) write (6,*) 'debug eigenvalues:',eigenvalues
! now we get go root by root, and construct the preconditioner for each
           do iroot = 1, nroots

              eigenvalue_r = real(eigenvalues(iroot))
              eigenvalue_i = aimag(eigenvalues(iroot))

              ierr = direct_data_pointer(workspace_t%precond(iroot)%s1,preco1_body_p)
              ierr = direct_data_pointer(workspace_t%precond(iroot)%s2,preco2_body_p)

              if (workspace_t%tens1_rank .eq. 1) then
                 call c_f_pointer(preco1_body_p,preco_1d_p,workspace_t%tens1_dims(1:workspace_t%tens1_rank))
                 call c_f_pointer(preco2_body_p,preco_3d_p,workspace_t%tens2_dims(1:workspace_t%tens2_rank))

                 do ja = 1, workspace_t%tens1_dims(1)
                    proj_r = eigenvalue_r - real(hdiag_1d_p(ja))
                    proj_i = eigenvalue_i - aimag(hdiag_1d_p(ja))
! aspg: added this check to filter out real or imaginary parts that are too small/noisy and would explode as denominators 
!       worst case this will set the preconditioner of this element to one (due to the initialization, see above)
                    if (dabs(proj_r).le.small_number) proj_r = 0.0d0
                    if (dabs(proj_i).le.small_number) proj_i = 0.0d0

                    x2py2 = proj_r*proj_r + proj_i*proj_i
                    if (dabs(x2py2).gt.config%safe_minimum) then
                       proj_r =  proj_r / x2py2
                       proj_i = -proj_i / x2py2
                       preco_1d_p(ja) = dcmplx(proj_r,proj_i)
                    end if
                 end do

                 do jc = 1, workspace_t%tens2_dims(3)
                    do jb = 1, workspace_t%tens2_dims(2)
                       do ja = 1, workspace_t%tens2_dims(1)
                          proj_r = eigenvalue_r - real(hdiag_3d_p(ja,jb,jc))
                          proj_i = eigenvalue_i - aimag(hdiag_3d_p(ja,jb,jc))
! aspg: added this check to filter out real or imaginary parts that are too small/noisy and would explode as denominators 
!       worst case this will set the preconditioner of this element to one (due to the initialization, see above)
                          if (dabs(proj_r).le.small_number) proj_r = 0.0d0
                          if (dabs(proj_i).le.small_number) proj_i = 0.0d0

                          x2py2 = proj_r*proj_r + proj_i*proj_i
                          if (dabs(x2py2).gt.config%safe_minimum) then
                             proj_r =  proj_r / x2py2
                             proj_i = -proj_i / x2py2
                             preco_3d_p(ja,jb,jc) = dcmplx(proj_r,proj_i)
                          end if
                       end do
                    end do
                 end do

              else if (workspace_t%tens1_rank .eq. 2) then
                 call c_f_pointer(preco1_body_p,preco_2d_p,workspace_t%tens1_dims(1:workspace_t%tens1_rank))
                 call c_f_pointer(preco2_body_p,preco_4d_p,workspace_t%tens2_dims(1:workspace_t%tens2_rank))

                 do jb = 1, workspace_t%tens1_dims(2)
                    do ja = 1, workspace_t%tens1_dims(1)
                       proj_r = eigenvalue_r - real(hdiag_2d_p(ja,jb))
                       proj_i = eigenvalue_i - aimag(hdiag_2d_p(ja,jb))
                       if (preconditioner_debug) then
                          write (6,*) 'debug safe minimum',config%safe_minimum
                          write (6,*) 'debug small number',small_number
                          write (6,*) 'debug for root ',iroot
                          write (6,*) 'debug E(',eigenvalue_r,',',eigenvalue_i,')'
                          write (6,*) 'debug Hii',real(hdiag_2d_p(ja,jb)),aimag(hdiag_2d_p(ja,jb))
                          write (6,*) 'debug projr',proj_r,' proji',proj_i
                       end if
! aspg: added this check to filter out real or imaginary parts that are too small/noisy and would explode as denominators 
!       worst case this will set the preconditioner of this element to one (due to the initialization, see above)
                       x2py2 = proj_r*proj_r + proj_i*proj_i
                       if (preconditioner_debug) then
                          write (6,*) 'debug x2py2',x2py2
                       end if
                       if (x2py2.le.small_number) x2py2 = small_number
                       proj_r =  proj_r / x2py2
                       proj_i = -proj_i / x2py2
                       preco_2d_p(ja,jb) = dcmplx(proj_r,proj_i)
                    end do
                 end do

                 do jd = 1, workspace_t%tens2_dims(4)
                    do jc = 1, workspace_t%tens2_dims(3)
                       do jb = 1, workspace_t%tens2_dims(2)
                          do ja = 1, workspace_t%tens2_dims(1)
                             proj_r = eigenvalue_r - real(hdiag_4d_p(ja,jb,jc,jd))
                             proj_i = eigenvalue_i - aimag(hdiag_4d_p(ja,jb,jc,jd))
                       if (preconditioner_debug) then
                          write (6,*) 'debug safe minimum',config%safe_minimum
                          write (6,*) 'debug small number',small_number
                          write (6,*) 'debug for root ',iroot
                          write (6,*) 'debug E(',eigenvalue_r,',',eigenvalue_i,')'
                          write (6,*) 'debug Hii',real(hdiag_4d_p(ja,jb,jc,jd)),aimag(hdiag_4d_p(ja,jb,jc,jd))
                          write (6,*) 'debug projr',proj_r,' proji',proj_i
                       end if

! aspg: added this check to filter out real or imaginary parts that are too small/noisy and would explode as denominators 
!       worst case this will set the preconditioner of this element to one (due to the initialization, see above)
                             if (dabs(proj_r).le.small_number) proj_r = 0.0d0
                             if (dabs(proj_i).le.small_number) proj_i = 0.0d0

                             x2py2 = proj_r*proj_r + proj_i*proj_i
                       if (preconditioner_debug) then
                          write (6,*) 'debug x2py2',x2py2
                       end if
                             if (dabs(x2py2).gt.config%safe_minimum) then
                                proj_r =  proj_r / x2py2
                                proj_i = -proj_i / x2py2
                                preco_4d_p(ja,jb,jc,jd) = dcmplx(proj_r,proj_i)
                             end if
                          end do
                       end do
                    end do
                 end do

              end if

              if (preconditioner_debug) then
                 write (6,*) 'debug >>> precon1_tens  root',iroot,' >>>'
                 call print_tensor(workspace_t%precond(iroot)%s1, 1.0d-18, "precond1")
                 write (6,*) 'debug <<< precon1_tens <<<'
                 write (6,*) 'debug >>> precon2_tens  root',iroot,' >>>'
                 call print_tensor(workspace_t%precond(iroot)%s2, 1.0d-18, "precond2")
                 write (6,*) 'debug <<< precon2_tens <<<'
              end if
           end do ! iroot
       end if ! create preconditioner
       if (preconditioner_debug) write (6,*) 'debug leaving preconditioner'

    end subroutine 


    subroutine davidson_transpose_conjugate_vector(v, vt, complex_dimension_v)
        real(kind=8), intent(in)    :: v(:,:)
        real(kind=8), intent(inout) :: vt(:,:)
        integer, intent(in) :: complex_dimension_v
        integer :: vectors_L, vectors_N
        real(kind=8) :: vr, vi, vtr, vti

        integer :: i, j, kr, ki, lr, li, N_step, L_step, vectors_N_complex, vectors_L_complex

        vectors_N = size(v,1)
        vectors_L = size(v,2)

        if (complex_dimension_v.eq.1) then
            N_step = 2 
            L_step = 1
        else if (complex_dimension_v.eq.2) then
            N_step = 1 
            L_step = 2
        end if
        vectors_N_complex = vectors_N / N_step 
        vectors_L_complex = vectors_L / L_step

        do i = 1, vectors_N_complex 
            do j = 1, vectors_L_complex
                if (complex_dimension_v.eq.1) then
                    kr = N_step*(i-1) + 1
                    ki = kr + 1

                    vr = v(kr,j)
                    vi = v(ki,j)
                    write (config%output_file_unit,'(a,1i4,a,1i4,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                                 'i:',i,' j:',j,' kr',kr,' ki:',ki,'  vr:',vr,'  vi:',vi 

                    vt(j,kr) =  vr
                    vt(j,ki) = -vi

                     write (config%output_file_unit,'(a,1i4,a,1i4,a,1e14.6,a,1e14.6,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           ' v(',i,',',j,')=(',vr,',',vi,') ->  vt(',j,',',i,')=(',vt(j,kr),',',vt(j,ki),')'
                else if (complex_dimension_v.eq.2) then

                    kr = L_step*(j-1) + 1
                    ki = kr + 1

                    vr = v(i, kr)
                    vi = v(i, ki)
                     write (config%output_file_unit,'(a,1i4,a,1i4,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           'i:',i,' j:',j,' kr',kr,' ki:',ki,'  vr:',vr,'  vi:',vi

                    vt(kr, i) =  vr
                    vt(ki, i) = -vi

                     write (config%output_file_unit,'(a,1i4,a,1i4,a,1e14.6,a,1e14.6,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           ' v(',i,',',j,')=(',vr,',',vi,') ->  vt(',j,',',i,')=(',vt(kr,i),',',vt(ki,i),')'

                end if
            end do
        end do

    end subroutine

    subroutine eti_solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left)
        complex(kind=8), intent(inout), pointer :: M(:,:)
        complex(kind=8), intent(inout), pointer :: eVectorsL(:,:)
        complex(kind=8), intent(inout), pointer :: eVectorsR(:,:)
        complex(kind=8), intent(inout), pointer :: eValues(:)
        logical, intent(in) :: left 

        integer :: vectors_N, vectors_L, i

        integer :: lwork, info
        character(len=1) :: do_left = 'N', do_right = 'V'
        real(kind=8), allocatable :: work(:)
        real(kind=8), allocatable :: rwork(:)

        vectors_N = size(M,1)
        vectors_L = size(M,2)

        if (associated(eVectorsL)) deallocate(eVectorsL)
        if (associated(eVectorsR)) deallocate(eVectorsR)
        if (associated(eValues))   deallocate(eValues)

        if (left) then
            do_left = 'V'
            do_right = 'N'
            allocate(eVectorsL(vectors_N,vectors_L))
            eVectorsL = 0.0d0
        else
            allocate(eVectorsL(1,vectors_L))
            eVectorsL = 0.0d0
        end if
        allocate(eVectorsR(vectors_N,vectors_L))
        allocate(eValues(vectors_L))
        eVectorsR = 0.0d0
        eValues   = 0.0d0

! first query optimal lwork value
        lwork = -1
        allocate(work(4*vectors_L))
!       print *, "debug, lwork in diagonalization before query:",lwork,4*vectors_L

! aspg, 17/04/2019
! check whether the factor of two is needed here after all
        allocate(rwork(2*vectors_N))
        call zgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                    eValues, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, rwork, INFO)
        deallocate(rwork)

        lwork = work(1)
!       print *, "debug, lwork in diagonalization after query:",lwork,work(1:4)
        deallocate(work)
        allocate(work(2*lwork))
        allocate(rwork(2*vectors_N))
        call zgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                    eValues, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, rwork, INFO)
        deallocate(rwork)
!       print *, 'in small_solve_eigenvalue: info',info
        if (info .gt. 0) then
            print *, 'in small_solve_eigenvalue: not all eigenvalues converged' 
            print *, 'converged'
            do i = info+1, vectors_L 
                print *, " Re, Im : ",eValues(i)
            end do
        end if
        deallocate(work)

    end subroutine
 
! aspg, added from (q)rsp
    subroutine eti_solve_Linear_systems(left, A, B_input, Solution)
        logical, optional, intent(in) :: left
        complex(kind=8), intent(inout), pointer :: A(:,:)
        complex(kind=8), intent(inout), pointer :: B_input(:)
        complex(kind=8), intent(inout), pointer :: Solution(:,:)

        character ::    FACT
        character ::    TRANS
        integer   ::    N
        integer   ::    NRHS
        integer   ::    LDA
        complex(kind=8),allocatable ::  AF(:,:)
        integer   ::    LDAF
        integer, allocatable   ::       IPIV(:)
        character ::    EQUED
        real(kind=8), allocatable ::    R(:)
        real(kind=8), allocatable ::    C(:)
        complex(kind=8),allocatable  :: B(:,:)
        integer ::      LDB
        complex(kind=8),allocatable ::  X(:,:)
        integer ::      LDX
        real(kind=8) :: RCOND
        real(kind=8), allocatable ::    FERR(:)
        real(kind=8), allocatable :: BERR(:)
        complex(kind=8), allocatable :: WORK(:)
        real(kind=8), allocatable ::    RWORK(:)
        integer ::      INFO
        integer :: dimen_x

        FACT = 'N'
        ! TRANS = 'N'
        N = size(A,1)
        NRHS = 1
        LDA = N
        LDAF = N
        LDX = N
        LDB = N
        EQUED = 'N'
        allocate(AF(LDAF,N))
        allocate(IPIV(N))

        allocate(R(N))
        allocate(C(N))

        allocate(B(LDB,NRHS))
        allocate(X(LDX,NRHS))
        allocate(FERR(NRHS))
        allocate(BERR(NRHS))
        allocate(WORK(2*N))
        allocate(RWORK(2*N))

        if (left) then 
            TRANS = 'T'
        else
            TRANS = 'N'
        end if 
        ! print *, "TRANS:", TRANS

        B(:,NRHS) = B_input(:)
        call zgesvx (FACT, TRANS, N, NRHS, A, LDA, AF, &
                    LDAF, IPIV, EQUED, R, C, B, LDB, &
                    X, LDX, RCOND, FERR, BERR,&
                    WORK, RWORK, INFO)

        ! call zgetrs (TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

        ! call zgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)

        dimen_x = size(X,1)
        allocate(Solution(dimen_x,NRHS))
        if (left) then 
            Solution(:,NRHS) = conjg(X(:,NRHS))
        else
            Solution(:,NRHS) = X(:,NRHS) 
        end if 

        ! dimen_x = size(B,1)
        ! allocate(Solution(dimen_x,NRHS))
        ! Solution(:,NRHS) = B(:,NRHS)

        deallocate(AF)
        deallocate(IPIV)
        deallocate(R)
        deallocate(C)
        deallocate(B)
        deallocate(X)
        deallocate(FERR)
        deallocate(BERR)
        deallocate(rwork)
        deallocate(work)

    end subroutine eti_solve_Linear_systems
! aspg, added from (q)rsp

    subroutine get_tensor_linear_systems(results_t,root_vect,workspace_t,tprop_tensor) ! solve_lhs)
      
        type(eti_davidson_results_tens_t), intent(inout) :: results_t
        complex(kind=8), pointer, contiguous :: root_vect(:,:)
        type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
!       logical, intent(in) :: solve_lhs
        type(eti_comm_tens_t), intent(inout) :: tprop_tensor
        

        type(tensor_t)          :: one_tensor
        integer(C_INT)              :: ierr
        integer                     :: iroot, N
        integer(C_INT)              :: one_dims(1)
        complex(kind=8)             :: factor

        one_dims(1) = 1
        ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)

        N = size(root_vect,1)
        ! print *, "N in get_tensor_linear_systems", N


        ! print *, "dim1", size(root_vect,1)
        ! print *, "dim2", size(root_vect,2)
        ! print *, "roor_vect", root_vect(:,1)
        ! call print_tensor(t_tensor%t1, 1.0d-10, "t_tensor%t1")
        

         do iroot=1, N
            factor = root_vect(iroot,1)
            ierr=t_einsum("t(a,i)+=B1(a,i)*X()",tprop_tensor%t1,workspace_t%trial(iroot)%s1,one_tensor,scale=factor)
            ierr=t_einsum("t(a,b,i,j)+=B2(a,b,i,j)*X()",tprop_tensor%t2,workspace_t%trial(iroot)%s2,one_tensor,scale=factor)
! aspg, todo for moving away from comm_t
!           if (solve_lhs) then 
!              ierr=t_einsum("t(a,i)+=B1(a,i)*X()",results_t%wf(iroot)%l1,workspace_t%trial(iroot)%s1,one_tensor,scale=factor)
!              ierr=t_einsum("t(a,b,i,j)+=B2(a,b,i,j)*X()",results_t%wf(iroot)%l2,workspace_t%trial(iroot)%s2,one_tensor,scale=factor)
!           else
!              ierr=t_einsum("t(a,i)+=B1(a,i)*X()",results_t%wf(iroot)%r1,workspace_t%trial(iroot)%s1,one_tensor,scale=factor)
!              ierr=t_einsum("t(a,b,i,j)+=B2(a,b,i,j)*X()",results_t%wf(iroot)%r2,workspace_t%trial(iroot)%s2,one_tensor,scale=factor)
!           endif
         end do
        !  call print_tensor(t_tensor%t1, 1.0d-10, "t_tensor%t1")
        !  call print_tensor(t_tensor%t2, 1.0d-10, "t_tensor%t2")


    end subroutine get_tensor_linear_systems


    subroutine get_tensor_eigenvector(result_t,root_vect,e_small,workspace_t,left)
      
        type(eti_davidson_results_tens_t), intent(inout) :: result_t
        complex(kind=8), pointer, contiguous :: root_vect(:,:)
        complex(kind=8), pointer, contiguous :: e_small(:)
        type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
        logical                     :: left
        

        type(tensor_t)          :: one_tensor
        integer(C_INT)              :: ierr
        integer                     :: iroot, N, N_state, istate
        integer(C_INT)              :: one_dims(1)
        integer(INTD), dimension(2)   :: ov_dims
        integer(INTD), dimension(4)   :: oovv_dims
        integer(INTD)               ::  ov_rank, oovv_rank
        complex(kind=8)             :: factor

        one_dims(1) = 1
        ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)



        N = size(root_vect,1)
!       N_state = size(root_vect,2)
! aspg: we only transform to the original basis the number of roots we solved for
        N_state = workspace_t%nroots
        ! print *, "N in get_tensor_linear_systems", N


        ! print *, "dim1", size(root_vect,1)
        ! print *, "dim2", size(root_vect,2)
        ! print *, "roor_vect", root_vect(:,1)
        ! call print_tensor(t_tensor%t1, 1.0d-10, "t_tensor%t1")

! aspg TODO
! replace these by fetching the information from the workspace
        ierr = get_nmodes(workspace_t%trial(1)%s1,ov_rank)
        ierr = get_extents_static(workspace_t%trial(1)%s1,ov_dims)

        ierr = get_nmodes(workspace_t%trial(1)%s2,oovv_rank)
        ierr = get_extents_static(workspace_t%trial(1)%s2,oovv_dims)

        allocate(result_t%wf(N_state))
        allocate(result_t%wf(N_state))
        if (left) then 
            allocate(result_t%values_lhs(N_state,1))
        else
            allocate(result_t%values_rhs(N_state,1))
        end if 

        print *, ""
        print *, "Number of eigenvectors          : ",N_state
        print *, "Subspace dimension (L)          : ",N
        print *, "Eigenvector dimension (# Det.)  : ",workspace_t%N
        print *, ""

        do istate=1, N_state
            if (left) then
                result_t%values_lhs(istate,1) = e_small(istate)
                ierr = construct_tensor_t(result_t%wf(istate)%l1,C8,ov_dims,init_val=ZERO)
                ierr = construct_tensor_t(result_t%wf(istate)%l2,C8,oovv_dims,init_val=ZERO)
            else
                result_t%values_rhs(istate,1) = e_small(istate)
                ierr = construct_tensor_t(result_t%wf(istate)%r1,C8,ov_dims,init_val=ZERO)   !for right vector, actually we use vo and vvoo
                ierr = construct_tensor_t(result_t%wf(istate)%r2,C8,oovv_dims,init_val=ZERO)
            end if


          do iroot=1, N
            factor = root_vect(iroot,istate)
            if (left) then 
                ierr=t_einsum("t(i,a)+=B1(i,a)*X()",result_t%wf(istate)%l1,&
                    workspace_t%trial(iroot)%s1,one_tensor,scale=factor)
                ierr=t_einsum("t(i,j,a,b)+=B2(i,j,a,b)*X()",result_t%wf(istate)%l2,&
                    workspace_t%trial(iroot)%s2,one_tensor,scale=factor) 
            else
                ierr=t_einsum("t(a,i)+=B1(a,i)*X()",result_t%wf(istate)%r1,&
                    workspace_t%trial(iroot)%s1,one_tensor,scale=factor)
                ierr=t_einsum("t(a,b,i,j)+=B2(a,b,i,j)*X()",result_t%wf(istate)%r2,&
                    workspace_t%trial(iroot)%s2,one_tensor,scale=factor) 
            end if 
          end do
          print *, " Done with state ",istate
        end do 

        !  call print_tensor(t_tensor%t1, 1.0d-10, "t_tensor%t1")
        !  call print_tensor(t_tensor%t2, 1.0d-10, "t_tensor%t2")


    end subroutine get_tensor_eigenvector



    subroutine eti_davidson_workspace_destroy(workspace_t)
        type(eti_davidson_workspace_tens_t),  intent(inout) :: workspace_t
        integer :: ierr, i

        do i = 1, workspace_t%Lmax
           ierr = destroy_tensor_t(workspace_t%trial(i)%s1)
           ierr = destroy_tensor_t(workspace_t%trial(i)%s2)
           ierr = destroy_tensor_t(workspace_t%sigma(i)%s1)
           ierr = destroy_tensor_t(workspace_t%sigma(i)%s2)
        end do
        deallocate(workspace_t%trial)
        deallocate(workspace_t%trial)
        deallocate(workspace_t%sigma)
        deallocate(workspace_t%sigma)

        do i = 1, workspace_t%nroots
           ierr = destroy_tensor_t(workspace_t%delta(i)%s1)
           ierr = destroy_tensor_t(workspace_t%delta(i)%s2)
           ierr = destroy_tensor_t(workspace_t%precond(i)%s1)
           ierr = destroy_tensor_t(workspace_t%precond(i)%s2)
        end do
        deallocate(workspace_t%delta)
        deallocate(workspace_t%delta)
        deallocate(workspace_t%precond)
        deallocate(workspace_t%precond)

        deallocate(workspace_t%tens1_dims)
        deallocate(workspace_t%tens2_dims)

        ierr = destroy_tensor_t(workspace_t%Hdiag%s1)
        ierr = destroy_tensor_t(workspace_t%Hdiag%s2)

        if (associated(workspace_t%N)) deallocate(workspace_t%N)
        if (associated(workspace_t%M)) deallocate(workspace_t%M)

        write(config%output_file_unit, "(4x,a)") 'leaving eti_davidson_workspace_destroy'
    end subroutine 

    subroutine eti_davidson_workspace_create(workspace_t, nroots, nocc, nvirt, flavor, excitation_level, left_right)
        type(eti_davidson_workspace_tens_t),  intent(inout) :: workspace_t
        character*2 :: flavor
        integer, intent(in) :: nocc, nvirt, nroots, excitation_level
        character*1, OPTIONAL                              :: left_right
        logical :: workspace_create_debug = .false.

        integer :: M_dims(2) 
        integer :: ierr, i
        integer :: n_continuum_occ, n_continuum_virt
        integer :: excitation_manifold_size
        integer :: remainder

        workspace_t%nroots = nroots
! start out by dealing with excitation level
        allocate(workspace_t%N(excitation_level))

! tensor ranks: for simplicity we will always work with even-ranked tensors.
!               in cases such as IP and EA, which could be described by
!               odd-rakned tensors (for single detachment/attachment, or
!               2h1p/1h2p), we define an extra class called "continuum" that
!               will denote the particle that is removed/added. this continnum
!               can have just dimension 1 for IP/EA, 2 for DIP/DEA etc
! singles
        workspace_t%tens1_rank = 2
        allocate(workspace_t%tens1_dims(workspace_t%tens1_rank))
! doubles
        workspace_t%tens2_rank = 4
        allocate(workspace_t%tens2_dims(workspace_t%tens2_rank))


        select case(flavor)
           case("IP")
              n_continuum_virt = 1
! singles
              workspace_t%tens1_dims = (/ n_continuum_virt , nocc /)
              workspace_t%N(1) = nocc*n_continuum_virt
! doubles
              workspace_t%tens2_dims = (/ nvirt, n_continuum_virt, nocc, nocc /)
              workspace_t%N(2) = n_continuum_virt*(nvirt**1)*(nocc**2)
            !   workspace_t%N(2) = nvirt*n_continuum_virt*(nocc*(nocc-1))/2

           case("EA")
              n_continuum_occ = 1
! singles
              workspace_t%tens1_dims = (/ nvirt, n_continuum_occ /)
              workspace_t%N(1) = nvirt*n_continuum_occ 
! doubles
              workspace_t%tens2_dims = (/ nvirt, nvirt, nocc, n_continuum_occ /)
              workspace_t%N(2) = nvirt*nvirt*nocc*n_continuum_occ

           case("EE")
! singles
              if (present(left_right)) then
                workspace_t%tens1_dims = (/ nocc, nvirt /)
              else
                workspace_t%tens1_dims = (/ nvirt, nocc /)
              end if 
              workspace_t%N(1) = nvirt*nocc
! doubles
              if (present(left_right)) then
                workspace_t%tens2_dims = (/ nocc, nocc, nvirt, nvirt /)
              else 
                workspace_t%tens2_dims = (/ nvirt, nvirt, nocc, nocc /)
              end if 
              workspace_t%N(2) = nvirt*nvirt*nocc*nocc

        end select

        excitation_manifold_size = workspace_t%N(1) + workspace_t%N(2)
 !       excitation_manifold_size = workspace_t%N(1) + (workspace_t%N(2)*(workspace_t%N(2)+1)/2)

! now that we know the dimensionality of the excitation manifold (=number of rows in Hbar), start allocating tensors
! the number of solutions (and solution vectors) depends on the user input.
! as was the case with relccsd, if we give a negative number for the number of roots, the code will assume we want
! as many roots as the there are columns in Hbar (=size of the excitation manifold)..

        workspace_t%Lmax = config%max_subspace_size
        remainder = mod(workspace_t%Lmax,workspace_t%nroots)
        if (nroots.lt.0) then
           print *, 'INFO: number of roots reset to the dimension of Hbar'
           workspace_t%nroots = excitation_manifold_size 
           workspace_t%Lmax   = excitation_manifold_size
        else if (nroots.gt.excitation_manifold_size) then
           print *, 'INFO: number of roots > dimension of Hbar, resetting it'
           workspace_t%nroots = min(nroots,excitation_manifold_size)
           print *, "excitation_manifold_size:", excitation_manifold_size
           print *, "max_subspace_size:", config%max_subspace_size
           workspace_t%Lmax   = min(config%max_subspace_size,excitation_manifold_size)
        else if (remainder.ne.0) then
           print *, 'INFO: input max. subspace size  not a multiple of number of roots, adjusting it'
           workspace_t%Lmax   = config%max_subspace_size + workspace_t%nroots - remainder 
        end if
        workspace_t%L      = workspace_t%nroots
        workspace_t%Lprev  = workspace_t%nroots

        allocate(workspace_t%trial(workspace_t%Lmax))
        allocate(workspace_t%sigma(workspace_t%Lmax))
        allocate(workspace_t%delta(workspace_t%nroots))
        allocate(workspace_t%precond(workspace_t%nroots))

!       M_dims     = (/ workspace_t%Lmax, workspace_t%Lmax /)
!       ierr=construct_tensor_t(workspace_t%M,C8,M_dims,init_val=ZERO)
        allocate(workspace_t%M(workspace_t%Lmax,workspace_t%Lmax))
        workspace_t%M = ZERO

        allocate(workspace_t%C(workspace_t%Lmax))
        workspace_t%C = ZERO
! singles
        do i = 1, workspace_t%Lmax
           ierr = construct_tensor_t(workspace_t%trial(i)%s1,C8,workspace_t%tens1_dims,init_val=ZERO)
           ierr = construct_tensor_t(workspace_t%sigma(i)%s1,C8,workspace_t%tens1_dims,init_val=ZERO)
        end do
        do i = 1, workspace_t%nroots
           ierr = construct_tensor_t(workspace_t%delta(i)%s1,C8,workspace_t%tens1_dims,init_val=ZERO)
           ierr = construct_tensor_t(workspace_t%precond(i)%s1,C8,workspace_t%tens1_dims,init_val=ZERO)
        end do
        ierr=construct_tensor_t(workspace_t%Hdiag%s1,C8,workspace_t%tens1_dims,init_val=ZERO)
        ierr=construct_tensor_t(workspace_t%Norm_prop_grad%s1,C8,workspace_t%tens1_dims,init_val=ZERO)

! doubles
        do i = 1, workspace_t%Lmax
           ierr=construct_tensor_t(workspace_t%trial(i)%s2,C8,workspace_t%tens2_dims,init_val=ZERO)
           ierr=construct_tensor_t(workspace_t%sigma(i)%s2,C8,workspace_t%tens2_dims,init_val=ZERO)
        end do
        do i = 1, workspace_t%nroots
           ierr = construct_tensor_t(workspace_t%delta(i)%s2,C8,workspace_t%tens2_dims,init_val=ZERO)
           ierr = construct_tensor_t(workspace_t%precond(i)%s2,C8,workspace_t%tens2_dims,init_val=ZERO)
        end do
        ierr=construct_tensor_t(workspace_t%Hdiag%s2,C8,workspace_t%tens2_dims,init_val=ZERO)
        ierr=construct_tensor_t(workspace_t%Norm_prop_grad%s2,C8,workspace_t%tens2_dims,init_val=ZERO)

        workspace_t%Nmax = workspace_t%N(1) + workspace_t%N(2)


! debug statements
        if (workspace_create_debug) then
           write(config%output_file_unit, "(4x,a)")  'in eti_davidson_workspace_create'
           write(config%output_file_unit, "(4x,a,i6,2x,a,2x,i6)") &
     &  'created tensor for trial vectors, dimensions :',workspace_t%Nmax,"x",workspace_t%Lmax

           do i = 1, workspace_t%Lmax
              write(config%output_file_unit, "(4x,a,i4)") 'info on workspace_t%trial%s1, element',i
              call eti_tensor_print_info(workspace_t%trial(i)%s1)
              write(config%output_file_unit, "(4x,a,i4)") 'info on workspace_t%trial%s2, element',i
              call eti_tensor_print_info(workspace_t%trial(i)%s2)
           end do

           write(config%output_file_unit, "(4x,a)") 'info on workspace_t%Hdiag%s1' 
           call eti_tensor_print_info(workspace_t%Hdiag%s1)
           write(config%output_file_unit, "(4x,a)") 'info on workspace_t%Hdiag%s2' 
           call eti_tensor_print_info(workspace_t%Hdiag%s2)

           write(config%output_file_unit, "(4x,a,i18,2x,a,2x,i6)") &
     &  'created tensor for sigma vectors, dimensions :',workspace_t%Nmax,"x",workspace_t%Lmax

           write(config%output_file_unit, "(4x,a,i18,2x,a,2x,i6)") &
     &  'created tensor for correction vectors, dimensions :',workspace_t%Nmax,"x",workspace_t%nroots
      end if

        
   end subroutine 


   subroutine create_sigma(workspace_t,int_t,interm_t,right,left,Index_small_v)
      use eti_sigma_eomccsd
      type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(eti_intg_tens_t), intent(inout) :: int_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      logical, intent(in) :: left,right
      integer,ALLOCATABLE :: Index_small_v(:)
      logical :: found
      integer :: ierr, iroot
      
      logical :: debug_sigma = .false., debug_print_before_calculating=.false.

    !   do iroot = (workspace_t%L - workspace_t%nroots + 1), workspace_t%L
      do iroot = workspace_t%Lprev, workspace_t%L
         found = is_number_in_array(iroot,Index_small_v,size(Index_small_v))
         if (.not.found) then

            if (debug_sigma) then
                write (6,*) 'root ',iroot
                if (debug_print_before_calculating) then
                    write (6,*) 'debug >>> t_tens before sigma, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    print "(2x,a6,2x,a6,2x,a6,2x,a6,2x,a6)", "idx", "p  ", "q  ","r  ","s  "
                    call print_tensor(workspace_t%trial(iroot)%s1, 0.0d-18, "trial1")
                    call print_tensor(workspace_t%trial(iroot)%s2, 0.0d-18, "trial2")
                    write (6,*) 'debug <<< t_tens before sigma <<<'

                    write (6,*) 'debug >>> s_tens before sigma, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    print "(2x,a6,2x,a6,2x,a6,2x,a6,2x,a6)", "idx", "p  ", "q  ","r  ","s  "
                    call print_tensor(workspace_t%sigma(iroot)%s1, 0.0d-18, "sigma1")
                    call print_tensor(workspace_t%sigma(iroot)%s2, 0.0d-18, "sigma2")
                    write (6,*) 'debug <<< s_tens before sigma <<<'
                else
                    write (6,*) 'debug <<< silenced printout prior to calculation >>>'  
                end if ! silencing printout prior to calculated vectors
            end if
! debug info     
! aspg
            if (right.and..not.left) then
                if (associated(interm_t%sigma_rhs_p)) then
                    call interm_t%sigma_rhs_p(int_t,interm_t,workspace_t%trial(iroot),workspace_t%sigma(iroot))
                else
                    stop 'RHS sigma vector routine has not been defined prior to calling the solver'
                end if
            else if (left.and..not.right) then
                if (associated(interm_t%sigma_lhs_p)) then
                    call interm_t%sigma_lhs_p(int_t,interm_t,workspace_t%trial(iroot),workspace_t%sigma(iroot))
                else
                    stop 'LHS sigma vector routine has not been defined prior to calling the solver'
                end if
            else
               stop  'Calculation of both left and right sigma not implemented'
            end if
        end if
! aspg
! debug info     
         if (debug_sigma) then
            write (6,*) 'debug >>> t_tens, root',iroot,' >>>'
            print "(2x,a6,2x,a6,2x,a6,2x,a6,2x,a6)", "idx", "p  ", "q  ","r  ","s  "
            call print_tensor(workspace_t%trial(iroot)%s1, 1.0d-18, "trial1")
            call print_tensor(workspace_t%trial(iroot)%s2, 1.0d-18, "trial2")
            write (6,*) 'debug <<< t_tens <<<'

            write (6,*) 'debug >>> s_tens, root',iroot,' >>>'
            print "(2x,a6,2x,a6,2x,a6,2x,a6,2x,a6)", "idx", "p  ", "q  ","r  ","s  "
            call print_tensor(workspace_t%sigma(iroot)%s1, 1.0d-18, "sigma1")
            call print_tensor(workspace_t%sigma(iroot)%s2, 1.0d-18, "sigma2")
            write (6,*) 'debug <<< s_tens <<<'
         end if
! debug info     
! aspg

      end do

!     call quit('bb')

   end subroutine

! aspg, added from (q)rsp
   subroutine create_subspace_property_gradient(C_small,workspace_t,prop_grad,right,left,Index_small_v,isv)
      type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t ! here we will keep the M(Lmax,Lmax)
      complex(kind=8), pointer, contiguous, intent(inout) :: C_small(:)   ! this will be a copy of B, B_small(L)
      logical, intent(in) :: right, left
      integer,ALLOCATABLE :: Index_small_v(:)
!     type(eti_sa_tens_t),intent(inout),optional :: prop_grad
      type(eti_comm_tens_t),intent(inout),optional :: prop_grad
      integer             :: isv
       
      integer :: mij_dims(2), mij_rank
      integer(INTD)      :: one_dims(1)
      type(tensor_t) :: result_tensor
      complex(kind=8), pointer ::  ci_result_tens(:)
      type(C_PTR):: ci_body_p

      integer :: m_dims(2), m_rank
      complex(kind=8), pointer, contiguous:: m_tens(:,:)
      type(C_PTR):: m_body_p
      integer :: ierr, i, j,  i_pri
      logical :: found_i
      logical :: debug_subspace_matrix = .false.

      if (associated(C_small)) deallocate(C_small)
      allocate(C_small(workspace_t%L))
      allocate(C_small((workspace_t%L)-isv))
      C_small = ZERO

!     Initialize scalars that are to be used as tensors in contractions
      one_dims(1) = 1
!     ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
      ierr=construct_tensor_t(result_tensor,C8,one_dims(1:0),init_val=ZERO)
      ierr=direct_data_pointer(result_tensor,ci_body_p)
      call c_f_pointer(ci_body_p,ci_result_tens,one_dims)
      

! capital letters here are meant to be the roots
      if (debug_subspace_matrix) then
         write (6,*) 'foobar, before setup M'
         write (6,*) 'foobar, workspace_t%L',workspace_t%L
         write (6,*) 'foobar, workspace_t%Lprev',workspace_t%Lprev
      end if

      i_pri = 0
      do i = workspace_t%Lprev, workspace_t%L
        found_i = is_number_in_array(i,Index_small_v,size(Index_small_v))

        if (debug_subspace_matrix) then
            call print_tensor(workspace_t%trial(i)%s1, 1.0d-18, "trial1")
            call print_tensor(prop_grad%t1, 1.0d-18, "prop_grad%s1")
            call print_tensor(workspace_t%trial(i)%s2, 1.0d-18, "trial2")
            call print_tensor(prop_grad%t2, 1.0d-18, "prop_grad%s2")
        end if 

        if (.not.found_i) then
! this will be a reduction to a scalar at each contractions, and we then will build the 
          if (workspace_t%N(1).ne.0) then
            if (right.and..not.left) then
               ierr=t_einsum("CI()+=Bdag+(p,q)*C(p,q)",             &
   &                                                       result_tensor,        &
   &                                                       workspace_t%trial(i)%s1,&
   &                                                       prop_grad%t1,scale=ONE)
            else if (left.and..not.right) then
               ierr=t_einsum("CI()+=C+(p,q)*B(p,q)",                &
   &                                                       result_tensor,        &
&                                                          prop_grad%t1,           &
   &                                                       workspace_t%trial(i)%s1,scale=ONE)
            end if
          end if ! workspace_t%N(1).ne.0

          if (workspace_t%N(2).ne.0) then
            if (right.and..not.left) then
               ierr=t_einsum("CI()+=Bdag+(p,q,r,s)*C(p,q,r,s)",     &
   &                                                       result_tensor,        &
   &                                                       workspace_t%trial(i)%s2,&
   &                                                       prop_grad%t2,           &
                                                           scale=ONE_QUARTER)
            else if (left.and..not.right) then
               ierr=t_einsum("CI()+=C+(p,q,r,s)*B(p,q,r,s)",       &
   &                                                       result_tensor,       &
   &                                                       prop_grad%t2,          &
                                                           workspace_t%trial(i)%s2,&
                                                           scale=ONE_QUARTER)
            end if
          end if

          if (debug_subspace_matrix) write (6,*) 'i, c_i:', i, ci_result_tens(1)
          workspace_t%C(i) = ci_result_tens(1)*(-1)
        !  workspace_t%C(i) = ci_result_tens(1)
          ci_result_tens(1) = ZERO

        end if
      end do
      C_small = workspace_t%C(1:workspace_t%L)

!     if (debug_subspace_matrix) then
!       if (.true.) then
!          write (6,*) 'debug, M small'
!          do i = 1, workspace_t%L
!            print "(2x,i6,$)", i
!            do j = 1, workspace_t%L
!                call print_complex_vector_element(M_small(i,j))
!                if (mod(j,4) == 0) then 
!                   print *, "  "
! !                 print "(8x,$)",""
!                end if
!            end do
!            print *, ""
!         end do
!      end if
!      if (.true.) call dump_to_file(2222, "C_small_EXACC.bin", C_small)
!      call quit('bb')

   end subroutine
! aspg, added from (q)rsp

   subroutine create_subspace_matrix(M_small,workspace_t,right,left,frequency,Index_small_v,isv)
      type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t ! here we will keep the M(Lmax,Lmax)
      complex(kind=8), pointer, contiguous, intent(inout) :: M_small(:,:)    ! this will be a copy of M, M_small(L,L)
      logical, intent(in) :: right, left
      complex(kind=8), intent(in), optional :: frequency  
      integer,ALLOCATABLE :: Index_small_v(:)
      integer             :: isv

      integer :: mij_dims(2), mij_rank
      integer(INTD)      :: one_dims(1)
      type(tensor_t) :: result_tensor
      complex(8), pointer :: mij_result_tens(:)
      type(C_PTR):: mij_body_p

      integer :: m_dims(2), m_rank
      complex(8), pointer, contiguous:: m_tens(:,:)
      type(C_PTR):: m_body_p
      integer :: ierr, i, j, i_pri, j_pri
      complex(8) :: contrac_coeff
      logical :: debug_subspace_matrix = .false.
      logical :: found_i, found_j

!! we need to fish out the small matrix from workspace_t%M, which is the slice LxL of the full tensor
!! we also need to conjugate transpose the trial vector 
      if (associated(M_small)) deallocate(M_small)
    !   allocate(M_small(workspace_t%L,workspace_t%L))
      allocate(M_small((workspace_t%L)-isv,(workspace_t%L)-isv))
      M_small = ZERO 

!     Initialize scalars that are to be used as tensors in contractions
      one_dims(1) = 1
      ierr=construct_tensor_t(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    !   print *, "ierr construct small M:", ierr 
      ierr=direct_data_pointer(result_tensor,mij_body_p)
      call c_f_pointer(mij_body_p,mij_result_tens,one_dims)

! capital letters here are meant to be the roots
      if (debug_subspace_matrix) then
         write (6,*) ''
         write (6,*) 'debug >>> before setting up subspace Hamiltonian matrix M'
         write (6,*) 'debug >>> Lprev ',workspace_t%Lprev
         write (6,*) 'debug >>> L     ',workspace_t%L
         write (6,*) ''
      end if
    !   do i = workspace_t%Lprev, workspace_t%L
!     print *, "workspce_tensor:", workspace_t%tens1_dims
      if (workspace_t%tens1_dims(1)==1.or.workspace_t%tens1_dims(2)==1) then  !To distinguish IP,EA and EE
        contrac_coeff = ONE_HALF   !It's IP, EA
      else 
        contrac_coeff = ONE_QUARTER  !It's EE
      end if

      i_pri = 0
      
      do i=1, workspace_t%L
        found_i = is_number_in_array(i,Index_small_v,size(Index_small_v))
        if (.not.found_i) i_pri=i_pri+1
        j_pri = 0
        do j =1, workspace_t%L
! this will be a reduction to a scalar at each contractions, and we then will build the
        found_j = is_number_in_array(j,Index_small_v,size(Index_small_v))
        if(.not.found_j) j_pri=j_pri+1
          
          if ((.not.found_i).and.(.not.found_j)) then
            if (workspace_t%N(1).ne.0) then
               if (right.and..not.left) then 
                  ierr=t_einsum("MIJ()+=Bdag+(p,q)*S(p,q)",          &
   &                                                        result_tensor,        &
   &                                                        workspace_t%trial(i)%s1,&
   &                                                        workspace_t%sigma(j)%s1)
               else if (left.and..not.right) then 
                  ierr=t_einsum("MIJ()+=S+(p,q)*B(p,q)",              &
   &                                                        result_tensor,        &
   &                                                        workspace_t%sigma(i)%s1,&
   &                                                        workspace_t%trial(j)%s1)
               end if
            end if

            if (workspace_t%N(2).ne.0) then
               if (right.and..not.left) then
                  ierr=t_einsum("MIJ()+=Bdag+(p,q,r,s)*S(p,q,r,s)",  &
   &                                                        result_tensor,        &
   &                                                        workspace_t%trial(i)%s2,&
   &                                                        workspace_t%sigma(j)%s2,&
   &                                                        scale=contrac_coeff)
               else if (left.and..not.right) then
                  ierr=t_einsum("MIJ()+=S+(p,q,r,s)*B(p,q,r,s)",      &
   &                                                        result_tensor,        &
   &                                                        workspace_t%sigma(i)%s2,&
   &                                                        workspace_t%trial(j)%s2,&
   &                                                        scale=contrac_coeff)
               end if
            end if

            if (debug_subspace_matrix) write (6,*) 'i,j, m_ij:', i,j, mij_result_tens(1)

            if (present(frequency)) then 
                if (i==j) then 
                    mij_result_tens(1) =mij_result_tens(1) - frequency
                end if 
            end if

            ! workspace_t%M(i,j) = mij_result_tens(1) 
            workspace_t%M(i_pri,j_pri) = mij_result_tens(1) 
            mij_result_tens(1) = ZERO
          end if
        end do
      end do

      M_small = workspace_t%M(1:(workspace_t%L-isv),1:(workspace_t%L-isv))

      if (debug_subspace_matrix) then
         write (6,*) ' >>> debug, subspace Hamiltonian matrix M'
         do i = 1, workspace_t%L-isv
           print "(2x,i6,$)", i
           do j = 1, workspace_t%L-isv
               call print_complex_vector_element(M_small(i,j))
               if (mod(j,4) == 0) then 
                  print *, "  "
!                 print "(8x,$)",""
               end if
           end do
           print *, ""
         end do
      end if
      if (.true.) call dump_to_file(2222, "M_small_EXACC.bin", M_small)

    end
    
    subroutine print_complex_vector_element(ve)
       complex(8), intent(in) :: ve

       real(8) :: threshold = 1.0d-10
       real(8) :: vr, vi

       vr = real(ve)
       vi = aimag(ve)

       if (abs(vi).lt.threshold) vi = 0.0d0
       if (abs(vr).lt.threshold) vr = 0.0d0
       if (vi < 0.0d0) then
           print "(4x,1e14.6,a,1e13.6,a,$)", vr," -",abs(vi)," i"
       else
           print "(4x,1e14.6,a,1e13.6,a,$)", vr," +",abs(vi)," i"
       end if

    end subroutine

    subroutine dump_to_file(fileno, filename, matrix)
        integer, intent(in) :: fileno
        character(len=*), intent(in) :: filename
        complex(kind=8) :: matrix(:,:)
        integer :: irec

        inquire (iolength = irec) matrix
        open (unit=fileno, file=filename, form="unformatted", access="direct", recl=irec)
        write(unit=fileno, rec=1) matrix
        close(unit=fileno)
    end subroutine

   subroutine eti_solve_linear_system_davidson(results_t,int_t,interm_t, &
 &            workspace_t, projector_t, frequency, prop_grad)
      type(eti_davidson_results_tens_t) :: results_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(eti_projectors_tens_t), intent(inout) :: projector_t
      type(eti_intg_tens_t), intent(inout) :: int_t
      character(len=2) :: flavor
!     type(eti_sa_tens_t), intent(inout) :: prop_grad
      type(eti_comm_tens_t), intent(inout) :: prop_grad
      complex(kind=8), intent(in) :: frequency

!     do something to populate the property gradient, within workspace_t, before calling the routine below
      call eti_get_values_vectors_davidson(results_t,int_t,interm_t, &
 &          workspace_t, projector_t, flavor,frequency, prop_grad)

   end subroutine



! aspg, 16/04/2019
! for the moment we have a long argument list to get it going, later we try to shorten in
!   subroutine eti_get_values_vectors_davidson(A, workspace, projector)
! removed property gradient since that will be avaliable from within workspace_t
   subroutine eti_get_values_vectors_davidson(results_t,int_t,interm_t, &
 &               workspace_t, projector_t,flavor, frequency, prop_grad, tprop_tensor)
      type(eti_davidson_results_tens_t), intent(inout) :: results_t
      type(eti_intermediates_tens_t),  intent(inout) :: interm_t
      type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(eti_projectors_tens_t), intent(inout) :: projector_t
      type(eti_intg_tens_t), intent(inout) :: int_t
!     type(eti_sa_tens_t), intent(inout), optional :: prop_grad
      type(eti_comm_tens_t), intent(inout), optional :: prop_grad
      character(len=2) :: flavor

      complex(kind=8), intent(in), optional :: frequency  
      type(eti_comm_tens_t), intent(inout), optional    :: tprop_tensor
      complex(kind=8), pointer, contiguous :: x_small(:,:)        => NULL()
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal

!       type(eti_lambda_tens_t), intent(in)   :: A
!       type(eti_davidson_workspace_tens_t),  intent(inout) :: workspace
!       type(eti_projectors), intent(inout) :: projector
!
! nullifying all pointers at the start to give them a state, see discussion on http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html
! "[...] When a pointer is declared its status is undefined, and cannot be safely queried with the associated intrinsic. [...]"

        complex(kind=8), pointer, contiguous :: M_small(:,:)        => NULL()
! aspg, added from (q)rsp
        complex(kind=8), pointer, contiguous :: C_small(:)          => NULL() 
       
! aspg, added from (q)rsp
        complex(kind=8), pointer, contiguous :: c_small_r(:,:)      => NULL()
        complex(kind=8), pointer, contiguous :: c_small_l(:,:)      => NULL()
        complex(kind=8), pointer, contiguous :: e_small(:)          => NULL()

        complex(kind=8), pointer, contiguous :: c_small_r_old(:,:)  => NULL()
        complex(kind=8), pointer, contiguous :: c_small_l_old(:,:)  => NULL()
        logical, pointer :: delta(:) => NULL()
        logical, pointer :: converged(:)

        complex(kind=8) :: frequency_local
        logical :: solve_linear = .false.
        logical :: solve_eigen  = .true.
        integer,ALLOCATABLE :: Index_small_v(:)
        integer :: isv

        integer :: niter
        integer :: i, j, k, l
        integer :: nroots, ierr
        logical :: left_eigen_vec

        ! if (config%solve_linear_system) then
        !     nroots = 1
        ! else
            nroots = workspace_t%nroots
        ! end if
        allocate (converged(nroots))
        converged = .false.

! aspg need to think of a consistency check for whether the property gradient has been set
!      to replace this
!       if (config%solve_linear_system.and..not.present(prop_gradient_t1)) then
!           call quit('bla')
!       end if

        if (present(frequency)) then
           frequency_local = frequency
        else
           frequency_local = (0.0d0,0.0d0) 
        end if

        allocate(Index_small_v(nroots*1000))
        Index_small_v = 0
        isv = 0
!
! davidson_create_results() will allocate storage for results from davidson diagonalization when called a first time. 
! 
!!      call davidson_create_results(results,vectors_N,nroots,&
!!                                   davidson_solve_values,   &
!!                                   davidson_solve_lhs,      &
!!                                   davidson_solve_rhs)

        workspace_t%L = nroots
        workspace_t%Lprev = 1

        write(config%output_file_unit,"(4x,a,i6,a)") "Requested ",nroots," eigenvales"
        if ( nroots .gt. workspace_t%Nmax ) then 
           write(config%output_file_unit, "(4x,a)") "Number of roots requested larger than matrix order !"
           return 
        else if ( nroots .eq. workspace_t%Nmax ) then
           write(config%output_file_unit, "(4x,a)") "Number of roots requested equal to matrix order !"
        endif

       if (config%overlap_sorting) &
           call setup_overlap_sorting(c_small_r_old, c_small_l_old, nroots, config%solve_rhs, config%solve_lhs)


        write(config%output_file_unit,*) ""
        do niter = 1, config%max_iterations
            write(config%output_file_unit,*) "Iteration ", niter
            write(config%output_file_unit,*) ""

!           call apply_projector(workspace_t%L, projector_t, workspace_t%trial%s1,  workspace_t%trial%s2)
            call create_sigma(workspace_t,int_t,interm_t,config%solve_rhs,config%solve_lhs,Index_small_v)
!           call apply_projector(workspace_t%L, projector_t, workspace_t%sigma%s1,  workspace_t%sigma%s2)
            call create_subspace_matrix(M_small,workspace_t,config%solve_rhs,config%solve_lhs,frequency_local,Index_small_v,isv)

! aspg, added from (q)rsp
            if (config%solve_linear_system) then
               call create_subspace_property_gradient(C_small,workspace_t, &
 &                  prop_grad,config%solve_rhs,config%solve_lhs,Index_small_v,isv)
               call calculate_Linear_systems(config%solve_lhs, M_small,C_small,x_small)
! aspg, added from (q)rsp
            else
               call calculate_eigenvalues_and_eigenvectors(config%output_file_unit,    &
                                                           M_small,                      &
                                                           e_small,                      &
                                                           c_small_r, c_small_l,         &
                                                           c_small_r_old, c_small_l_old, &
                                                           config%eigenvalues_shift,   &
                                                           follow_overlap=(config%overlap_sorting.and.(niter.gt.1)),  &
                                                           follow_overlap_scheme=config%overlap_sorting_scheme, &
                                                           apply_evalues_shift=config%apply_eigenvalues_shift, &
                                                           left=config%solve_lhs,      &
                                                           right=config%solve_rhs,     &
                                                           nroots=nroots,                &
                                                           verbose=config%verbose)
!                                                          verbose=.true.)
            end if

            if (config%solve_rhs.and.config%overlap_sorting) then
                c_small_r_old = c_small_r
            end if
            if (config%solve_lhs.and.config%overlap_sorting) then
                c_small_l_old = c_small_l
            end if
 
! aspg, added from (q)rsp
            if (config%solve_linear_system) then
                allocate(e_small(1))
                e_small(1) = frequency_local
            end if
! aspg, added from (q)rsp
            call form_preconditioner(workspace_t, e_small)

! aspg: preparing to create new correction vectors (if any)
!
! first we get deallocate any vectors from the prior iteration
            if (associated(delta)) deallocate(delta)

! now we call correction_to_vectors() to calculate the vectors and check
! whether or not the correction is within the thresholds
            if (config%solve_linear_system) then
               delta => correction_to_vectors(config%output_file_unit, &
   &                                          e_small, x_small, workspace_t, &
                                              Index_small_v,prop_grad)
            else
               if (config%solve_rhs) then
                  delta => correction_to_vectors(config%output_file_unit, &
   &                                             e_small, c_small_r, workspace_t,Index_small_v)
               else 
                  delta => correction_to_vectors(config%output_file_unit, &
   &                                             e_small, c_small_l, workspace_t,Index_small_v)
               end if
            end if

! if there are no new vectors (=delta is not associated)  this means we've converged
! and we exit the iterative loop
            if (.not.associated(delta)) then
                write(config%output_file_unit,*) "All solutions converged, stopping."
                converged = .true.
                exit
            end if

! otherwise...
! not converged yet, if asked to refresh the trial vectors do it now, otherwise create
! new trial vectors
            if (workspace_t%L .ge. workspace_t%Lmax) then
               write(config%output_file_unit,*) "Reached maximum subspace size, stopping."
               converged = .true. 
               do l = 1, nroots
                  if (delta(l)) then
                     converged(l) = .false. 
                     write(config%output_file_unit,*) "WARNING ! solution ",l," is not converged."
                  end if
               end do 
               if (associated(delta)) deallocate(delta)
               exit
            end if

            if (config%refresh_trial_vectors) then
!              if (mod(niter,config%trial_vector_refresh_rate).eq.0) then
!                  print *,'debug : refreshing trial vectors not implemented in this version'
!              end if
               print *,'debug : refreshing trial vectors not implemented in this version'
            else
               call update_trial_vectors(delta,workspace_t,flavor,Index_small_v,isv)
            end if
 
        end do ! davidson iterations

        if (config%solve_rhs.and.config%overlap_sorting) deallocate(c_small_r_old)
        if (config%solve_lhs.and.config%overlap_sorting) deallocate(c_small_l_old)

        write(config%output_file_unit, *) "Final size of reduced subspace : ",workspace_t%L

        if (config%solve_linear_system) then
            write(config%output_file_unit, *) "**************** Solution of linear equations *****************"
            write(config%output_file_unit, *) "        Only solution vectors are printed (if requested)" 
        else
            write(config%output_file_unit, *) "************** Solution of eigenvalue equations ***************"
            call print_values_vectors(config%output_file_unit,e_small, c_small_r, c_small_l, &
                .false.,.false., nroots)
        end if
        write(config%output_file_unit, *) "***************************************************************"

        
!!      call form_eigenvalues_full_problem(results%eValues, nroots, e_small) 

!!      if (config%solve_lhs) then 
!!          call form_eigenvector_full_problem(results%eVectorsL, nroots,b,c_small_l,config%solve_rhs,config%solve_lhs)
!!      end if
!!      if (config%solve_rhs) then 
!!          call form_eigenvector_full_problem(results%eVectorsR, nroots,b,c_small_r,config%solve_rhs,config%solve_lhs)
!!      end if

        write(config%output_file_unit,*) "Generating solution vectors"

        if (config%solve_linear_system) then
             call get_tensor_linear_systems(results_t,x_small,workspace_t,tprop_tensor)
!           if (config%solve_rhs) then
!               call get_tensor_linear_systems(results_t,x_small,workspace_t,.false.)
!           else
!               call get_tensor_linear_systems(results_t,x_small,workspace_t,.true.)
!           endif
        else
            if (config%solve_rhs) then
                left_eigen_vec = .FALSE.
                call get_tensor_eigenvector(results_t,c_small_r,e_small,workspace_t,left_eigen_vec)
            else
                left_eigen_vec = .True.
                call get_tensor_eigenvector(results_t,c_small_l,e_small,workspace_t,left_eigen_vec)
            end if
        end if

!!      if (associated(b)) deallocate(b)
!!      deallocate(M_small)
!!      if (associated(sigma_vector)) deallocate(sigma_vector)
        deallocate(e_small)
        if (associated(c_small_r)) deallocate(c_small_r)
        if (associated(c_small_l)) deallocate(c_small_l)
        if (associated(converged)) deallocate(converged)
        DEALLOCATE(Index_small_v)

    end subroutine




    ! subroutine rotate_eigenvector(sigma1, sigma2, workspace_t)
    !     type(eti_davidson_workspace_tens_t), intent(inout) :: workspace_t
    !     type(tensor_t) sigma1, sigma2


    !     type(C_PTR):: body_p
    !     integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
    !     integer(C_INT):: host_arg_max

    !     integer(C_INT)              :: ierr
    !     integer(INTD), dimension(1) :: one_dims

    !     integer :: i,j,a 
    !     integer :: N_trial1, N_trail2

    !     type(tensor_t)   :: result_tensor
    !     complex(8), pointer  :: result_tens(:), sigma1_tens(:), sigma2_tens(:)


    !     one_dims(1) = 1
    !     N_trail1 = size(workspace_t%trial%s1,1)
    !     N_trail2 = size(workspace_t%trial%s2,1)


    !     ierr=construct_tensor_t(result_tensor,C8,one_dims,init_val=ZERO)

    !     ierr=direct_data_pointer(result_tensor,body_p)
    !     call c_f_pointer(body_p,result_tens,one_dims)


    !     do i=1, N_trail1
    !         ierr=t_einsum("V1()+=B(p,q)*S(p,q)",result_tensor,workspace_t%trial(i)%s1,sigma1) ! 
    !         sigma1_tens(i) = result_tens(1)
    !     end do 


    !     do i=1, N_trail2
    !         ierr=t_einsum("V2()+=B(a,b,c)*S(a,b,c)",result_tensor,workspace_t%trial(i)%s2,sigma2) ! 
    !     end do 



    !     ierr=t_einsum("V1(a,b)+=B(a,m)*V(m,b)",eigenvec_tensor_new1,workspace_t%trial%s1,eigenvec_tensor_old) ! 
    !     ierr=t_einsum("V2(a,b)+=B(a,m)*V(m,b)",eigenvec_tensor_new2,workspace_t%trial%s2,eigenvec_tensor_old) ! 


    !     ierr=direct_data_pointer(eigenvec_tensor_new1,body_p)
    !     call c_f_pointer(body_p,eigenvec_tens1,eigenvec_dims_s)

    !     ierr=direct_data_pointer(eigenvec_tensor_new2,body_p)
    !     call c_f_pointer(body_p,eigenvec_tens2,eigenvec_dims_d)




    !     ierr=destroy_tensor_t(eigenvec_tensor_new1)
    !     ierr=destroy_tensor_t(eigenvec_tensor_new2)

    ! end subroutine


    subroutine setup_overlap_sorting(c_r_old, c_l_old, nroots, solve_rhs, solve_lhs)
        complex(kind=8), intent(inout), pointer :: c_r_old(:,:) 
        complex(kind=8), intent(inout), pointer :: c_l_old(:,:)  
        integer, intent(in) :: nroots
        logical, intent(in) :: solve_rhs, solve_lhs

        integer :: i, factor 

        if (solve_rhs) then
            allocate(c_r_old(nroots,nroots))
            c_r_old = 0.0d0
            do i = 1, nroots
                c_r_old(i,i) = (1.0d0, 0.0d0)
            enddo
        endif
        if (solve_lhs) then
            allocate(c_l_old(nroots,nroots))
            c_l_old = 0.0d0
            do i = 1, nroots
                c_r_old(i,i) = (1.0d0, 0.0d0)
            enddo
        endif
    end subroutine

    subroutine calculate_eigenvalues_and_eigenvectors(output_file_unit, M, eValues, eVectorsR, eVectorsL, &
                                eVectorsR_old, eVectorsL_old, eValuesShift, left, right,    &
                                sorting_order, nroots, verbose, follow_overlap, follow_overlap_scheme, apply_evalues_shift)
        complex(kind=8), intent(inout), pointer   :: M(:,:)
        complex(kind=8), intent(inout), pointer   :: eValues(:)
        complex(kind=8), intent(inout), pointer   :: eVectorsR(:,:), eVectorsL(:,:)
        complex(kind=8), intent(in)   :: eVectorsR_old(:,:), eVectorsL_old(:,:) !old vectors
        real(kind=8), intent(in)   :: eValuesShift(2)
        logical, optional, intent(in) :: left, right, verbose, follow_overlap, apply_evalues_shift
        integer, optional, intent(in) :: sorting_order, nroots, follow_overlap_scheme
        integer, intent(in) :: output_file_unit 

        logical :: left_local, right_local, verbose_local, follow_overlap_local, apply_evalues_shift_local
        integer :: sorting_order_local, M_nroots, nroot, nroots_to_print, follow_overlap_scheme_local
        complex(kind=8), allocatable :: overlap(:,:)

        external zgemm

        if (present(apply_evalues_shift)) then
            apply_evalues_shift_local = apply_evalues_shift
        else
            apply_evalues_shift_local = .false.
        end if

        if (present(follow_overlap)) then
            follow_overlap_local = follow_overlap
        else
            follow_overlap_local = .false.
        end if

        if (present(follow_overlap_scheme)) then
            follow_overlap_scheme_local = follow_overlap_scheme
        else
            follow_overlap_scheme_local = 2 ! 1 = scheme calculating overlap matrix ; 2 = scheme not explicitly calculating overlap matrix
        end if

        if (present(sorting_order)) then
            sorting_order_local = sorting_order
        else
            sorting_order_local = 1 
        end if

        if (present(right)) then
            right_local = right
        else
            right_local = .false.
        end if

        if (present(left)) then
            left_local = left
        else
            left_local = .false. 
        end if

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false. 
        end if

        M_nroots = size(M,1) 

        if (present(nroots)) then
            nroots_to_print = nroots
        else
            nroots_to_print = M_nroots
        end if

        ! call davidson_print_vectors(M,   label="subspace matrix before diagonal")
        call eti_solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left_local)

        if (follow_overlap_local) then
            
            if (follow_overlap_scheme_local .eq. 1) then
! calculate overlap between old and new vectors. Overlap(nroot,mroot) = Old^T(nroot,nroot)*New(1:nroot,mroot)   
! Max(Overlap(iroot,mroot)) for all iroot = 1,nroot 
                nroot = size(eVectorsR_old,1)
                allocate(overlap(nroot,M_nroots))
                overlap = 0.0d0 
                if(right_local) then    
! aspg, check consistency of complex and real modes after real mode is validated
                    call zgemm('C', 'N', nroot, M_nroots, nroot, (1.0d0, 0.0d0), eVectorsR_old, &
                               nroot, eVectorsR, nroot, (1.0d0, 0.0d0), overlap, nroot)
                else if (left_local) then
! aspg, check consistency of complex and real modes after real mode is validated
                     call zgemm('C', 'N', nroot, M_nroots, nroot, (1.0d0, 0.0d0), eVectorsL_old, &
                                nroot, eVectorsL, nroot, (1.0d0, 0.0d0), overlap, nroot)
                endif

                call sort_by_overlap(-1, eValues, overlap, eVectorsR, eVectorsL, &
                                     right_local, left_local) 
                deallocate(overlap)

            else if (follow_overlap_scheme_local .eq. 2) then
! aspg, debug, testing new overlap sorting
                if(right_local) then
                    call sort_by_overlap_of_subspace(eValues, eVectorsR, nroot)
                else if (left_local) then
                    call sort_by_overlap_of_subspace(eValues, eVectorsL, nroot)
                end if
            end if ! selection of overlap sorting scheme 

        else
            call eti_sort_values_vectors(sorting_order_local, eValues, eVectorsR, eVectorsL, &
                                     right_local, left_local, apply_evalues_shift, eValuesShift)
        end if

! debug, print eigenvecturs and eigenvalues from subspace
        ! print *, "eVectorsR(1) after zgemm", eVectorsR(:,1)
        ! print *, "eVectorsR(2) after zgemm", eVectorsR(:,2)
        ! print *, "eVectorsR(3) after zgemm", eVectorsR(:,3)
        ! print *, "eVectorsR(4) after zgemm", eVectorsR(:,4)
        ! print *, "eVectorsR(5) after zgemm", eVectorsR(:,5)
        ! print *, "eVectorsR(6) after zgemm", eVectorsR(:,6)
        ! print *, "eVectorsR(7) after zgemm", eVectorsR(:,7)
        ! print *, "eVectorsR(8) after zgemm", eVectorsR(:,8)
        call print_values_vectors(config%output_file_unit,eValues, eVectorsR, eVectorsL, &
                                    .false.,.false., nroots=nroots_to_print)
!                                   .false.,.true., nroots=nroots_to_print)
!                                   left_local, right_local, nroots=nroots_to_print)

    end subroutine


! aspg, added from (q)rsp
    subroutine calculate_Linear_systems(left, M, B, Solution)

        logical, optional, intent(in) :: left
        complex(kind=8), intent(inout), pointer :: M(:,:)
        complex(kind=8), intent(inout), pointer :: B(:)
        complex(kind=8), intent(inout), pointer :: Solution(:,:)

        call eti_solve_Linear_systems(left, M, B, Solution)
    end subroutine calculate_Linear_systems
! aspg, added from (q)rsp


    subroutine form_eigenvalues_full_problem(eValues, nroots, e_small)
        complex(kind=8), intent(inout) :: eValues(:)
        complex(kind=8), intent(in)    :: e_small(:)
        integer, intent(in) :: nroots

        integer :: i

        do i = 1, nroots
            eValues(i) = e_small(i)
        end do
    end subroutine

    subroutine form_eigenvector_full_problem(eVectorsFull,nroots,b,eVectorsSmall,right_vector,left_vector) 
        complex(kind=8), intent(inout) :: eVectorsFull(:,:)
        complex(kind=8), intent(in)    :: eVectorsSmall(:,:)
        complex(kind=8), intent(in)    :: b(:,:)
        logical, intent(in) :: left_vector, right_vector
        integer, intent(in) :: nroots

        integer :: i
        complex(kind=8), allocatable :: bc(:,:)
        integer :: vectors_L, vectors_N

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (left_vector) then
            allocate(bc(vectors_L,vectors_N))
            bc = 0.0d0
! aspg: verify here whether transpose or conjugate transpose
            print *, 'warning, calling zgemm with T instead of C in form_eigenvector_full_problem'
            call zgemm('T', 'T', vectors_L, vectors_N, vectors_L, (1.0d0, 0.0d0), eVectorsSmall, &
                       vectors_L, b, vectors_N, (1.0d0, 0.0d0), bc, vectors_L)
            do i = 1, nroots
                eVectorsFull(:,i) = bc(i,:) 
            end do
        end if

        if (right_vector) then 
            allocate(bc(vectors_N,vectors_L))
            bc = 0.0d0
            call zgemm('N', 'N', vectors_N, vectors_L, vectors_L, (1.0d0, 0.0d0), b, &
                       vectors_N, eVectorsSmall, vectors_L, (1.0d0, 0.0d0), bc, vectors_N)
            do i = 1, nroots
                eVectorsFull(:,i) = bc(:,i) 
            end do
        end if

        deallocate (bc)
    end subroutine


    function update_sigma_vector(M, b, old_sigma_vector)
        real(kind=8), pointer :: update_sigma_vector(:,:)
        real(kind=8), intent(in)  :: M(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        real(kind=8), pointer :: old_sigma_vector(:,:)

        integer :: vectors_L, vectors_N
        integer :: vectors_L_old
        integer :: new_L, i

        external zgemm

        vectors_N = size(b,1)
        vectors_L = size(b,2)
        vectors_L_old = size(old_sigma_vector,2)

        allocate(update_sigma_vector(vectors_N,vectors_L))
        do i = 1, vectors_L_old 
           update_sigma_vector(:, i) = old_sigma_vector(:, i)  
        end do
        deallocate(old_sigma_vector)

        new_L = vectors_L - vectors_L_old
        if (config%complex_mode_with_reals) then
            vectors_N = vectors_N / 2
            call zgemm('N', 'N', vectors_N, new_L, vectors_N, (1.0d0, 0.0d0), M, vectors_N, &
       &               b(:, vectors_L_old:vectors_L), vectors_N, (0.0d0, 0.0d0),            &
       &               update_sigma_vector(:,vectors_L_old:vectors_L), vectors_N)
         end if
    end function


    function create_sigma_vector_full_matrix(M, b)
        real(kind=8), pointer :: create_sigma_vector_full_matrix(:,:)
        real(kind=8), intent(in)  :: M(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        integer :: vectors_L, vectors_N

        external zgemm

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_sigma_vector_full_matrix(vectors_N,vectors_L))

        if (config%complex_mode_with_reals) then
            vectors_N = vectors_N / 2
            call zgemm('N', 'N', vectors_N, vectors_L, vectors_N, (1.0d0, 0.0d0), M, vectors_N, &
                       b, vectors_N, (0.0d0, 0.0d0), create_sigma_vector_full_matrix, vectors_N)
        end if
    end function

    subroutine eti_sort_values_vectors(order, e_small, c_small_r, c_small_l, right, left, apply_shift, energy_shift) 
        integer, intent(in) :: order ! 1: ascending, -1: descending
        logical, intent(in) :: right ! if true, sort also the right eigenvectors
        logical, intent(in) :: left  ! if true, sort also the left eigenvectors
        logical, intent(in) :: apply_shift ! if true, sort values of |e_small - energy_shift| instead of e_small
        complex(kind=8), intent(inout) :: e_small(:)
        complex(kind=8), intent(inout) :: c_small_r(:,:)
        complex(kind=8), intent(inout) :: c_small_l(:,:)
        real(kind=8), intent(in) :: energy_shift(2)

        complex(kind=8), allocatable :: temp_vector(:,:)
        real(kind=8) :: temp
        integer, allocatable :: indices_sorted_array(:)
        real(kind=8), allocatable :: e_small_real(:,:)
        integer :: vectors_N, vectors_L, i_sorted, i

        vectors_L = size(e_small)
        allocate(e_small_real(2,vectors_L))
        allocate(indices_sorted_array(vectors_L))
        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
            e_small_real(1,i) = real(e_small(i))
            e_small_real(2,i) = aimag(e_small(i))
        end do

        if (apply_shift) then
           call quicksort_vector_complex_as_real_shifted(energy_shift, e_small_real, indices_sorted_array, 1, vectors_L)
        else
           call quicksort_vector_complex_as_real(e_small_real, indices_sorted_array, 1, vectors_L)
        end if

        do i = 1, vectors_L
            e_small(i) = dcmplx(e_small_real(1,i), e_small_real(2,i))
            !e_small(i) = complex(e_small_real(1,i), e_small_real(2,i))
        end do

        allocate(temp_vector(vectors_L,vectors_L))
        temp_vector = dcmplx(0.0d0, 0.0d0)

        if (right) then
            do i = 1, vectors_L
             i_sorted = indices_sorted_array(i)
             temp_vector(:,i) = c_small_r(:,i_sorted)
            end do
            c_small_r = temp_vector
        end if

        if (left) then
            do i = 1, vectors_L
             i_sorted = indices_sorted_array(i)
             temp_vector(:,i) = c_small_l(:,i_sorted)
            end do
            c_small_l = temp_vector
        end if
        deallocate(temp_vector)
        deallocate(indices_sorted_array)
        deallocate(e_small_real)

    end subroutine


    subroutine sort_by_overlap(order, e_small, overlap, c_small_r, c_small_l, right, left) 
        integer, intent(in) :: order ! 1: ascending, -1: descending
        logical, intent(in) :: right ! if true, sort also the right eigenvectors
        logical, intent(in) :: left  ! if true, sort also the left eigenvectors
        complex(kind=8), intent(inout) :: overlap(:,:)
        complex(kind=8), intent(inout) :: c_small_r(:,:)
        complex(kind=8), intent(inout) :: c_small_l(:,:)
        complex(kind=8), intent(inout) :: e_small(:)

        complex(kind=8), allocatable :: temp_vector(:,:)
        real(kind=8), allocatable :: overlap_copy(:,:)
        complex(kind=8),allocatable :: temp(:)
        integer, allocatable :: indices_sorted_array(:), indices_sorted_array_temp(:)
        integer :: vectors_L, i_sorted, i, max_overlap(1), iroot, nroot, jroot, kroot
        real(kind=8) :: vr, vi

        vectors_L = size(overlap,2)
        nroot     = size(overlap,1)

        allocate(indices_sorted_array(vectors_L))
        allocate(indices_sorted_array_temp(vectors_L))
        allocate(temp_vector(vectors_L,vectors_L))
        allocate(overlap_copy(2,vectors_L))

        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
        end do

        do iroot = 1, nroot
         
          do i = 1, vectors_L 
              indices_sorted_array_temp(i) = i
              overlap_copy(1,i) = real(overlap(iroot,i))
              overlap_copy(2,i) = aimag(overlap(iroot,i))
          end do
          call davidson_print_vectors(overlap_copy,   label="copy of overlap matrix before  ")
          call quicksort_vector_complex_as_real(overlap_copy, indices_sorted_array_temp, 1, vectors_L)
          call davidson_print_vectors(overlap_copy,   label="copy of overlap matrix after   ")
          max_overlap(1) =  indices_sorted_array_temp(vectors_L)
          
          i = max_overlap(1)
          indices_sorted_array(iroot) = i
        enddo

        do i = 1, nroot
          if (indices_sorted_array(i) > nroot) then
           indices_sorted_array(indices_sorted_array(i)) = i
          endif
        enddo

        if (right) then
            do i = 1, vectors_L
             i_sorted = indices_sorted_array(i)
             temp_vector(:,i) = c_small_r(:,i_sorted)
            end do
            c_small_r = temp_vector
        end if

        if (left) then
            do i = 1, vectors_L
             i_sorted = indices_sorted_array(i)
             temp_vector(:,i) = c_small_l(:,i_sorted)
            end do
            c_small_l = temp_vector
        end if

        deallocate(temp_vector)

        allocate(temp(vectors_L))
        do i = 1, vectors_L
          i_sorted = indices_sorted_array(i)
          temp(i) = e_small(i_sorted)
        end do
        e_small = temp

        deallocate(temp_vector)
        deallocate(indices_sorted_array)

    end subroutine


    real(kind=8) function shifted_value_complex_as_real(unshifted_value, shift)
       real(kind=8), intent(in)  ::  shift(2), unshifted_value(2)
       real(kind=8)  :: temp(2)

       temp(1) = unshifted_value(1) - shift(1)
       temp(2) = unshifted_value(2) - shift(2)

       shifted_value_complex_as_real = sqrt(temp(1)*temp(1) + temp(2)*temp(2))
    end function shifted_value_complex_as_real


    recursive subroutine quicksort_vector_complex_as_real(v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:,:)
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, t, value
        integer :: i, j, ti, start_index

        start_index = (first + last) / 2 

        x = sqrt(v(1,start_index)*v(1,start_index) + v(2,start_index)*v(2,start_index))

        i = first
        j = last

        do
           do while ((sqrt(v(1,i)*v(1,i) + v(2,i)*v(2,i)) < x))
              i = i + 1
           end do
           do while ((x < sqrt(v(1,j)*v(1,j) + v(2,j)*v(2,j))))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(1,i)
           v(1,i) = v(1,j)
           v(1,j) = t

           t = v(2,i)
           v(2,i) = v(2,j)
           v(2,j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) then 
            call quicksort_vector_complex_as_real(v, indices_sorted, first, i-1)
        end if
        if (j+1 < last)  then
            call quicksort_vector_complex_as_real(v, indices_sorted, j+1, last)
        end if
    end subroutine

    recursive subroutine quicksort_vector_complex_as_real_shifted(shift, v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:,:)
        real(kind=8), intent(in) :: shift(2) 
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, t, value
        integer :: i, j, ti, start_index

        start_index = (first + last) / 2

        x = shifted_value_complex_as_real(v(:,start_index), shift)

        i = first
        j = last

        do
           do while ( shifted_value_complex_as_real(v(:,i), shift) < x)
              i = i + 1
           end do
           do while ( x < shifted_value_complex_as_real(v(:,j), shift))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(1,i)
           v(1,i) = v(1,j)
           v(1,j) = t

           t = v(2,i)
           v(2,i) = v(2,j)
           v(2,j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) then
            call quicksort_vector_complex_as_real_shifted(shift, v, indices_sorted, first, i-1)
        end if
        if (j+1 < last)  then
            call quicksort_vector_complex_as_real_shifted(shift, v, indices_sorted, j+1, last)
        end if
    end subroutine

! taken from http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    character(len=20) function str(k)
!       "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

    subroutine sort_by_overlap_of_subspace(e_small, c_small, nroots)
!
! this routine sorts the eigenvectors with respect to the norm of their projection on
! to the set of trial vectors. Sorting is from high overlap to low overlap,
! so that the first nroots vectors can be taken to extend the space.
!
! Lucas Visscher (2019), based on a routine by ASP Gomes for the initial implementation 
! of EOM in relccsd.  
! ASP Gomes (2022), readapted from the L Visscher code refactor to only use complex 
! variables (and thus not use rcw toggle)
!
        integer, intent(in)         :: nroots
        complex(kind=8), intent(inout) :: c_small(:,:) ! left or right-side reduced subspace eigenvectors (can be real or complex)
        complex(kind=8), intent(inout) :: e_small(:)    ! reduced subspace eigenvalues (complex, as we deal with a non-hermitian eigenvalue problem)

        real(kind=8), allocatable :: overlap_with_subspace(:)
        complex(kind=8), allocatable :: sorted_eigenvalues(:)
        complex(kind=8), allocatable :: sorted_eigenvectors(:,:)

        real(kind=8) :: overlap
        integer :: subspace_dimension, iroot, itrial, i

        subspace_dimension = size(c_small,2)

        allocate(overlap_with_subspace(subspace_dimension))
        overlap_with_subspace = 0.d0
        do iroot = 1, subspace_dimension
           do itrial = 1, nroots
              overlap = abs(c_small(itrial,iroot))
              overlap_with_subspace(iroot) = overlap_with_subspace(iroot) + overlap**2
           end do
        end do

        allocate (sorted_eigenvalues(subspace_dimension))
        allocate (sorted_eigenvectors(subspace_dimension,subspace_dimension))
        do iroot = 1, subspace_dimension
           i = maxloc(overlap_with_subspace,1)
           sorted_eigenvalues(iroot)  = e_small(i)
           sorted_eigenvectors(:,iroot) = c_small(:,i)
           overlap_with_subspace(i) = -max(1.d-12,overlap_with_subspace(i))
           ! with this choice we select the next highest value for the next root without loosing the overlap information, the
           ! max construct is there in case we have multiple vectors with zero overlap.
        end do

        ! copy the sorted vectors back to the original arrays
        e_small = sorted_eigenvalues
        c_small = sorted_eigenvectors

        deallocate (sorted_eigenvalues)
        deallocate (sorted_eigenvectors)
        deallocate (overlap_with_subspace)

    end subroutine sort_by_overlap_of_subspace




    subroutine apply_projector(nvectors, projector_t, target1, target2) 
       integer, intent(in) :: nvectors
       type(eti_projectors_tens_t), intent(inout) :: projector_t
       type(tensor_t), intent(inout) :: target1(:), target2(:)

       integer :: p, q, r, s, ierr, i
       complex(8), pointer, contiguous:: target1_tens(:,:), target2_tens(:,:,:,:)
       type(C_PTR):: target1_body_p, target2_body_p
       complex(8), pointer, contiguous:: proj1_tens(:,:), proj2_tens(:,:,:,:)
       type(C_PTR):: proj1_body_p, proj2_body_p

       logical :: debug = .true.

       ierr = direct_data_pointer(projector_t%p1,proj1_body_p)
       ierr = direct_data_pointer(projector_t%p2,proj2_body_p)
       call c_f_pointer(proj1_body_p,proj1_tens,projector_t%p1_dims(1:projector_t%p1_rank))
       call c_f_pointer(proj2_body_p,proj2_tens,projector_t%p2_dims(1:projector_t%p2_rank))

       do i = 1, nvectors

          ierr = direct_data_pointer(target1(i),target1_body_p)
          ierr = direct_data_pointer(target2(i),target2_body_p)

          call c_f_pointer(target1_body_p,target1_tens,projector_t%p1_dims(1:projector_t%p1_rank))
          call c_f_pointer(target2_body_p,target2_tens,projector_t%p2_dims(1:projector_t%p2_rank))

          if (debug) then
             write (6,*) 'debug >>> in apply_projector, before projection, p1 >>>'
             call print_tensor(projector_t%p1, 1.0d-18, "p1")
             write (6,*) 'debug <<< in apply_projector, before projection, p1 <<<'
             write (6,*) 'debug <<< in apply_projector, before projection, vector1 ',i,' >>>'
             call print_tensor(target1(i), 1.0d-18, "target1")
             write (6,*) 'debug <<< in apply_projector, before projection, vector1 <<<'
          end if
!         ierr=t_einsum("T(p,q)+=T(p,q)*P(p,q)",target1(i),target1(i),p1)
          do q = 1, projector_t%p1_dims(2)
             do p = 1, projector_t%p1_dims(1)
                 target1_tens(p,q) = target1_tens(p,q)*proj1_tens(p,q)
             end do 
          end do 
          if (debug) then
             write (6,*) 'debug <<< in apply_projector, after  projection, vector1 ',i,' >>>'
             call print_tensor(target1(i), 1.0d-18, "target1")
             write (6,*) 'debug <<< in apply_projector, after  projection, vector1 <<<'
          end if

          if (debug) then
             write (6,*) 'debug >>> in apply_projector, before projection, p2 >>>'
             call print_tensor(projector_t%p2, 1.0d-18, "p2")
             write (6,*) 'debug <<< in apply_projector, before projection, p2 <<<'
             write (6,*) 'debug <<< in apply_projector, before projection, vector2 ',i,' >>>'
             call print_tensor(target2(i), 1.0d-18, "target1")
             write (6,*) 'debug <<< in apply_projector, before projection, vector2 <<<'
          end if
!         ierr=t_einsum("T(p,q,r,s)+=T(p,q,r,s)*P(p,q,r,s)",target2(i),target2(i),p2)
          do s = 1, projector_t%p2_dims(4)
             do r = 1, projector_t%p2_dims(3)
                do q = 1, projector_t%p2_dims(2)
                   do p = 1, projector_t%p2_dims(1)
                      target2_tens(p,q,r,s) = target2_tens(p,q,r,s)*proj2_tens(p,q,r,s)
                   end do 
                 end do 
             end do 
          end do 
          if (debug) then
             write (6,*) 'debug <<< in apply_projector, after  projection, vector2 ',i,' >>>'
             call print_tensor(target2(i), 1.0d-18, "target2")
             write (6,*) 'debug <<< in apply_projector, after  projection, vector2 <<<'
          end if

       end do
    end subroutine

    subroutine eti_davidson_get_eigenvalue_shift(apply_shift,shift)
       logical, intent(inout)      :: apply_shift
       real(kind=8), intent(inout) :: shift(2)

       shift = config%eigenvalues_shift
       apply_shift = config%apply_eigenvalues_shift 

    end subroutine

    subroutine divide_tensor_elements_by_value(a_tensor, denominator)
       type(tensor_t), intent(inout)    :: a_tensor
       complex(8), intent(in) :: denominator

       integer(INTD)  ::  a_rank
       integer(INTD)  ::  a_dims1(1)
       integer(INTD)  ::  a_dims2(1:2)
       integer(INTD)  ::  a_dims3(1:3)
       integer(INTD)  ::  a_dims4(1:4)
       integer(INTD)  ::  a_dims5(1:5)
       integer(INTD)  ::  a_dims6(1:6)
       integer(INTD)  ::  a_dims7(1:7)
       integer(INTD)  ::  a_dims8(1:8)
       integer(INTD)  ::  a_dims9(1:9)
       integer(INTD)  ::  a_dims10(1:10)
       integer(INTD)  ::  ierr
       integer  ::  i, j, k, l, m, n, o, p, q, r
       type(C_PTR)         :: a_body
       complex(8), pointer :: a_tens1(:)
       complex(8), pointer :: a_tens2(:,:)
       complex(8), pointer :: a_tens3(:,:,:)
       complex(8), pointer :: a_tens4(:,:,:,:)
       complex(8), pointer :: a_tens5(:,:,:,:,:)
       complex(8), pointer :: a_tens6(:,:,:,:,:,:)
       complex(8), pointer :: a_tens7(:,:,:,:,:,:,:)
#ifndef VAR_NVFORTRAN
       complex(8), pointer :: a_tens8(:,:,:,:,:,:,:,:)
       complex(8), pointer :: a_tens9(:,:,:,:,:,:,:,:,:)
       complex(8), pointer :: a_tens10(:,:,:,:,:,:,:,:,:,:)
#endif

       ierr = get_nmodes(a_tensor,a_rank)

       if (a_rank.eq.1) then 
         ierr = get_extents_static(a_tensor,a_dims1)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens1,a_dims1)

         do i = 1, a_dims2(1)
            a_tens1(i)=a_tens1(i) / denominator 
         end do

       else if (a_rank.eq.2) then
         ierr = get_extents_static(a_tensor,a_dims2)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens2,a_dims2)

         do i = 1, a_dims2(1)
         do j = 1, a_dims2(2)
            a_tens2(i,j)=a_tens2(i,j) / denominator
         end do
         end do

       else if (a_rank.eq.3) then 
         ierr = get_extents_static(a_tensor,a_dims3)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens3,a_dims3)

         do i = 1, a_dims3(1)
         do j = 1, a_dims3(2)
         do k = 1, a_dims3(3)
            a_tens3(i,j,k)=a_tens3(i,j,k) / denominator
         end do
         end do
         end do

       else if (a_rank.eq.4) then 
         ierr = get_extents_static(a_tensor,a_dims4)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens4,a_dims4)

         do i = 1, a_dims4(1)
         do j = 1, a_dims4(2)
         do k = 1, a_dims4(3)
         do l = 1, a_dims4(4)
             a_tens4(i,j,k,l)=a_tens4(i,j,k,l) / denominator 
         end do
         end do
         end do
         end do

      else if (a_rank.eq.5) then
         ierr = get_extents_static(a_tensor,a_dims5)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens5,a_dims5)

         do i = 1, a_dims5(1)
         do j = 1, a_dims5(2)
         do k = 1, a_dims5(3)
         do l = 1, a_dims5(4)
         do m = 1, a_dims5(5)
             a_tens5(i,j,k,l,m)=a_tens5(i,j,k,l,m) / denominator
         end do ! m
         end do ! l
         end do ! k
         end do ! j
         end do ! i

      else if (a_rank.eq.6) then
         ierr = get_extents_static(a_tensor,a_dims6)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens6,a_dims6)

         do i = 1, a_dims6(1)
         do j = 1, a_dims6(2)
         do k = 1, a_dims6(3)
         do l = 1, a_dims6(4)
         do m = 1, a_dims6(5)
         do n = 1, a_dims6(6)
             a_tens6(i,j,k,l,m,n)=a_tens6(i,j,k,l,m,n) / denominator
         end do ! n
         end do ! m
         end do ! l
         end do ! k
         end do ! j
         end do ! i

      else if (a_rank.eq.7) then
         ierr = get_extents_static(a_tensor,a_dims7)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens7,a_dims7)

         do i = 1, a_dims7(1)
         do j = 1, a_dims7(2)
         do k = 1, a_dims7(3)
         do l = 1, a_dims7(4)
         do m = 1, a_dims7(5)
         do n = 1, a_dims7(6)
         do o = 1, a_dims7(7)
             a_tens7(i,j,k,l,m,n,o)=a_tens7(i,j,k,l,m,n,o) / denominator
         end do ! o
         end do ! n
         end do ! m
         end do ! l
         end do ! k
         end do ! j
         end do ! i

#ifndef VAR_NVFORTRAN
      else if (a_rank.eq.8) then
         ierr = get_extents_static(a_tensor,a_dims8)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens8,a_dims8)

         do i = 1, a_dims8(1)
         do j = 1, a_dims8(2)
         do k = 1, a_dims8(3)
         do l = 1, a_dims8(4)
         do m = 1, a_dims8(5)
         do n = 1, a_dims8(6)
         do o = 1, a_dims8(7)
         do p = 1, a_dims8(8)
             a_tens8(i,j,k,l,m,n,o,p)=a_tens8(i,j,k,l,m,n,o,p) / denominator
         end do ! p
         end do ! o
         end do ! n
         end do ! m
         end do ! l
         end do ! k
         end do ! j
         end do ! i
      else if (a_rank.eq.9) then
         ierr = get_extents_static(a_tensor,a_dims9)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens9,a_dims9)

         do i = 1, a_dims9(1)
         do j = 1, a_dims9(2)
         do k = 1, a_dims9(3)
         do l = 1, a_dims9(4)
         do m = 1, a_dims9(5)
         do n = 1, a_dims9(6)
         do o = 1, a_dims9(7)
         do p = 1, a_dims9(8)
         do q = 1, a_dims9(9)
             a_tens9(i,j,k,l,m,n,o,p,q)=a_tens9(i,j,k,l,m,n,o,p,q) / denominator
         end do ! q
         end do ! p
         end do ! o
         end do ! n
         end do ! m
         end do ! l
         end do ! k
         end do ! j
         end do ! i

      else if (a_rank.eq.10) then
         ierr = get_extents_static(a_tensor,a_dims10)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens10,a_dims10)

         do i = 1, a_dims10(1)
         do j = 1, a_dims10(2)
         do k = 1, a_dims10(3)
         do l = 1, a_dims10(4)
         do m = 1, a_dims10(5)
         do n = 1, a_dims10(6)
         do o = 1, a_dims10(7)
         do p = 1, a_dims10(8)
         do q = 1, a_dims10(9)
         do r = 1, a_dims10(10)
             a_tens10(i,j,k,l,m,n,o,p,q,r)=a_tens10(i,j,k,l,m,n,o,p,q,r) / denominator
         end do ! r
         end do ! q
         end do ! p
         end do ! o
         end do ! n
         end do ! m
         end do ! l
         end do ! k
         end do ! j
         end do ! i
 
#endif
       else
          stop 'error in divide_tensor_elements_by_value: only ranks 0 to 10 implemented (0 to 7 for NVFORTRAN)'
       end if



    end subroutine

    subroutine multiply_tensor_elements(a_tensor, b_tensor)
! this routine performs element-wise multiplication of two tensors
! and stores the outcome in the first argument
!
! a(i,j,k,...) = a(i,j,k,..)*b(i,j,k,...)
!
       implicit none
    
       type(tensor_t), intent(inout)    :: a_tensor,b_tensor
    
       integer(INTD)  ::  a_rank, b_rank
       integer(INTD), allocatable  ::  a_dims(:),b_dims(:)
       integer(INTD)  ::  ierr

       type(C_PTR)         :: a_body,b_body
       complex(8), pointer :: a_tens1(:)
       complex(8), pointer :: b_tens1(:)
       complex(8), pointer :: a_tens2(:,:)
       complex(8), pointer :: b_tens2(:,:)
       complex(8), pointer :: a_tens3(:,:,:)
       complex(8), pointer :: b_tens3(:,:,:)
       complex(8), pointer :: a_tens4(:,:,:,:)
       complex(8), pointer :: b_tens4(:,:,:,:)
       complex(8), pointer :: a_tens5(:,:,:,:,:)
       complex(8), pointer :: b_tens5(:,:,:,:,:)
       complex(8), pointer :: a_tens6(:,:,:,:,:,:)
       complex(8), pointer :: b_tens6(:,:,:,:,:,:)
       complex(8), pointer :: a_tens7(:,:,:,:,:,:,:)
       complex(8), pointer :: b_tens7(:,:,:,:,:,:,:)
#ifndef VAR_NVFORTRAN
       complex(8), pointer :: a_tens8(:,:,:,:,:,:,:,:)
       complex(8), pointer :: b_tens8(:,:,:,:,:,:,:,:)
       complex(8), pointer :: a_tens9(:,:,:,:,:,:,:,:,:)
       complex(8), pointer :: b_tens9(:,:,:,:,:,:,:,:,:)
       complex(8), pointer :: a_tens10(:,:,:,:,:,:,:,:,:,:)
       complex(8), pointer :: b_tens10(:,:,:,:,:,:,:,:,:,:)
#endif
       complex(8) :: mult_val 
    
       ierr = get_nmodes(a_tensor, a_rank)
       ierr = get_nmodes(b_tensor, b_rank)
       allocate(a_dims(a_rank))
       allocate(b_dims(b_rank))
       if (a_rank.ne.b_rank) then
          call quit('incompatible ranks for tensors in in multiply_tensor_elements')
       end if

       ierr = get_extents(a_tensor,a_dims)
       if (ierr.ne.0) call quit('error in getting tensor extents in multiply_tensor_elements')
       ierr = get_extents(b_tensor,b_dims)
       if (ierr.ne.0) call quit('error in getting tensor extents in multiply_tensor_elements')
    
       ! only checking a_rank in what follos since we've established a and b have same rank 
       if (a_rank.eq.1) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens1,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens1,b_dims)
    
         a_tens1 = a_tens1 * b_tens1
    
       else if (a_rank.eq.2) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens2,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens2,b_dims)
    
         a_tens2 = a_tens2 * b_tens2
    
       else if (a_rank.eq.3) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens3,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens3,b_dims)
    
         a_tens3 = a_tens3 * b_tens3
    
       else if (a_rank.eq.4) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens4,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens4,b_dims)
    
         a_tens4 =  a_tens4 * b_tens4

       else if (a_rank.eq.5) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens5,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens5,b_dims)

         a_tens5 =  a_tens5 * b_tens5

       else if (a_rank.eq.6) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens6,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens6,b_dims)

         a_tens6 =  a_tens6 * b_tens6

       else if ((a_rank.eq.7).and.(b_rank.eq.7)) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens7,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens7,b_dims)

         a_tens7 =  a_tens7 * b_tens7

#ifndef VAR_NVFORTRAN
       else if (a_rank.eq.8) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens8,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens8,b_dims)

         a_tens8 =  a_tens8 * b_tens8

       else if (a_rank.eq.9) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens9,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens9,b_dims)

         a_tens9 =  a_tens9 * b_tens9

       else if (a_rank.eq.10) then
         ierr=direct_data_pointer(a_tensor,a_body)
         call c_f_pointer(a_body,a_tens10,a_dims)
         ierr=direct_data_pointer(b_tensor,b_body)
         call c_f_pointer(b_body,b_tens10,b_dims)

         a_tens10 =  a_tens10 * b_tens10

#endif
       else
          stop 'error in multiply_tensor_elements: only ranks 0 to 10 implemented (0 to 7 for NVFORTRAN)'
       end if
       deallocate(a_dims)
       deallocate(b_dims)
    
    end subroutine 

    subroutine print_tensor_elements2(tensor)
       implicit none
    
       type(tensor_t), intent(inout)    :: tensor
    
       integer(INTD)  ::  rank
       integer(INTD), allocatable  ::  dims(:)
       integer(INTD)  ::  ierr
       integer  ::  i, j, k, l, b
       type(C_PTR)         :: body
       complex(8), pointer :: tens1(:)
       complex(8), pointer :: tens2(:,:)
       complex(8), pointer :: tens3(:,:,:)
       complex(8), pointer :: tens4(:,:,:,:)
       complex(8), pointer :: a_tens5(:,:,:,:,:)
       complex(8), pointer :: a_tens6(:,:,:,:,:,:)
       complex(8), pointer :: a_tens7(:,:,:,:,:,:,:)
#ifndef VAR_NVFORTRAN
       complex(8), pointer :: a_tens8(:,:,:,:,:,:,:,:)
       complex(8), pointer :: a_tens9(:,:,:,:,:,:,:,:,:)
       complex(8), pointer :: a_tens10(:,:,:,:,:,:,:,:,:,:)
#endif
       real(8)    :: threshold = 1.0e-10
       complex(8) :: value_to_print
    
       ierr = get_nmodes(tensor,rank)
       allocate(dims(rank))
       ierr = get_extents(tensor,dims)
       ierr = direct_data_pointer(tensor,body)
    
       if (rank.eq.1) then
         call c_f_pointer(body,tens1,dims)
    
         b = 0
         do i = 1, dims(1)
            b = b + 1
            value_to_print = tens1(i)
            if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = ZERO
            else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = dcmplx(real(value_to_print),0.0d0) 
            else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = dcmplx(0.0d0,aimag(value_to_print)) 
            end if
            write(6,'(4x,i6,1e14.6,2x,SP,1e14.6,"i")') b, value_to_print
         end do
    
       else if (rank.eq.2) then
         call c_f_pointer(body,tens2,dims)
    
         b = 0
         do j = 1, dims(2)
           do i = 1, dims(1)
            b = b + 1
            value_to_print = tens2(i,j)
            if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = ZERO
            else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = dcmplx(real(value_to_print),0.0d0)
            else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = dcmplx(0.0d0,aimag(value_to_print))
            end if
            write(6,'(4x,i6,2x,i6,2x,i6,2x,1e14.6,2x,SP,1e14.6,"i")') b, i, j, value_to_print

           end do
         end do
    
       else if (rank.eq.3) then
         call c_f_pointer(body,tens3,dims)

         do k = 1, dims(3)
            do j = 1, dims(2)
              do i = 1, dims(1)
                 b = b + 1
                 value_to_print = tens3(i,j,k)
                 if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                    value_to_print = ZERO
                 else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
                    value_to_print = dcmplx(real(value_to_print),0.0d0)
                 else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                    value_to_print = dcmplx(0.0d0,aimag(value_to_print))
                 end if
                 write(6,'(4x,i6,2x,i6,2x,i6,2x,i6,2x,1e14.6,2x,SP,1e14.6,"i")') b, i, j, k, value_to_print

             end do
           end do
         end do
    
       else if (rank.eq.4) then
         call c_f_pointer(body,tens4,dims)

         b = 0
         do l = 1, dims(4)
            do k = 1, dims(3)
              do j = 1, dims(2)
                do i = 1, dims(1)
                   b = b + 1
                   value_to_print = tens4(i,j,k,l)
                   if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                      value_to_print = ZERO
                   else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
                      value_to_print = dcmplx(real(value_to_print),0.0d0)
                   else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                      value_to_print = dcmplx(0.0d0,aimag(value_to_print))
                   end if
                   write(6,'(4x,i6,2x,i6,2x,i6,2x,i6,2x,i6,2x,1e14.6,2x,SP,1e14.6,"i")') b, i, j, k, l, value_to_print

                end do
             end do
           end do
         end do
    
       else
          call quit('unsupported tensor rank in print_tensor_elements')
       end if
       deallocate(dims)
    
    end subroutine 

end module

