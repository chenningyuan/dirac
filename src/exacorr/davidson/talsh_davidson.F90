module talsh_davidson

use talsh
use tensor_algebra
use talsh_common_routines
use exacorr_datatypes
use talsh_common_routines, only : tensor_norm2


implicit none
complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                         ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                         MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                         MINUS_TWO=(-2.D0,0.D0)

public talsh_davidson_get_eigenvalue_shift
public talsh_davidson_setup
public talsh_davidson_destroy
public talsh_davidson_driver
public davidson_print_vectors_complex
public talsh_solve_values_vectors
public talsh_sort_values_vectors
public talsh_davidson_workspace_create
public talsh_davidson_workspace_destroy

interface talsh_davidson_driver
    module procedure talsh_get_values_vectors_davidson
end interface

interface davidson_print_vectors
    module procedure davidson_print_vectors_real
    module procedure davidson_print_vectors_complex
end interface

interface print_values_vectors
     module procedure print_values_vectors_R
     module procedure print_values_vectors_C
end interface

interface davidson_print_vectors_complex
    module procedure davidson_print_vectors_complex_1d
end interface

public talsh_davidson_results_tens_t
public talsh_davidson_workspace_tens_t


private

! derived types

type talsh_davidson_results_tens_t
    complex(kind=8), pointer :: values_rhs(:,:) 
    complex(kind=8), pointer :: values_lhs(:,:) ! keep left and right energies to verify whether they are the same 
    type(talsh_tens_t), pointer :: rhs1(:), rhs2(:) 
    type(talsh_tens_t), pointer :: lhs1(:), lhs2(:) 
end type talsh_davidson_results_tens_t

! this is the workspace storage for a given instance of the solver
type :: talsh_davidson_workspace_tens_t
   integer :: nroots              ! number of solutions we're looking for
   integer, pointer :: N(:)       ! number of rows per excitation level
   integer :: Nmax                ! total number of rows =  dim(T1)+dim(T2)+...
   integer :: Lmax                ! total number of columns = dim(maximum number of trial vectors)
   integer :: L, Lprev            ! counters for keeping tabs on how much of the workspace we've used
                                  ! in most usages N = Nmax and L will grow up to Lmax
                                  ! with Lprev keeping the number of the prior iteration (staring at 0)
! aspg, 24/03/2020
! original idea was to that  to M to be a tensor, but using arrays of tensors as below, better to keep it as a pointer
!          type(talsh_tens_t) :: M        ! subspace matrix, dimension : (Lmax, Lmax); it is in tensor 
!                                         ! form to conform with the operations, but it will be a pointer to it 
!                                         ! that will be passed to the lapack routines 
   complex(8), pointer, contiguous:: M(:,:) !  subspace matrix, dimension : (Lmax, Lmax)
! aspg, added from (q)rsp
   complex(8), pointer, contiguous:: C(:)   !  subspace vector, dimension : (Lmax)
! aspg, added from (q)rsp
! see to change those below to tensors segmented into singly/doubly/... manifolds
   type(talsh_tens_t), pointer :: trial1(:),   trial2(:)   ! workspace for trial vectors
   type(talsh_tens_t), pointer :: sigma1(:), sigma2(:)     ! workspace for sigma vectors
   type(talsh_tens_t), pointer :: delta1(:), delta2(:)     ! workspace for correction vectors
   type(talsh_tens_t), pointer :: precond1(:), precond2(:) ! workspace for preconditioner
   type(talsh_tens_t), pointer :: metric1(:), metric2(:)   ! workspace for metric, for generalized eigenvalue problem, or linear equation problem 
   type(talsh_tens_t) :: Hdiag1, Hdiag2                    ! workspace for diagonal elements of the similarity-transformed Hamiltonian
   type(talsh_tens_t) :: Norm_prop_grad                    ! workspace for norm of each element of property gradient vectors C in LR
! this so far is not really used outside of davidson, if it is, we should create a type
   integer, pointer :: tens1_dims(:), tens2_dims(:)        ! since all tensors for a given model IP/EA/EE/... have the same dimensions, we keep this information
   integer :: tens1_rank, tens2_rank                       ! since all tensors for a given model IP/EA/EE/... have the same rank, we keep this information
end type talsh_davidson_workspace_tens_t



! module parameters

integer, parameter :: max_label_length = 30
integer, parameter :: nr_max_irreps = 32
real(kind=8), parameter :: davidson_eigenvalues_shift_threshold = 1.0e-3

! module variables, should be set with a call to davidson_setup

real(kind=8), save :: safe_minimum = 1.0d-36

real(kind=8), save :: davidson_convergence_threshold 
real(kind=8), save :: davidson_antisymme_threshold
real(kind=8), save :: davidson_eigenvalues_shift(2)

integer, save      :: davidson_max_subspace_size 
integer, save      :: davidson_max_iterations
integer, save      :: davidson_output_file_unit 
integer, save      :: davidson_subspaces_file_unit 
integer, save      :: davidson_sigmas_file_unit 
integer, save      :: davidson_results_file_unit 
integer, save      :: davidson_overlap_sorting_scheme
integer, save      :: davidson_trial_vector_refresh_rate
integer, save      :: davidson_maximum_trial_vector_refresh_rate = 10

logical, save      :: davidson_apply_eigenvalues_shift
logical, save      :: davidson_solve_rhs 
logical, save      :: davidson_solve_lhs 
logical, save      :: davidson_solve_values
logical, save      :: davidson_solve_linear_system
logical, save      :: davidson_verbose 
logical, save      :: davidson_symmetric 
logical, save      :: davidson_complex_mode_with_reals 
logical, save      :: davidson_save_subspaces 
logical, save      :: davidson_save_sigmas 
logical, save      :: davidson_save_results
logical, save      :: davidson_overlap_sorting
logical, save      :: davidson_debug = .false.
logical, save      :: davidson_use_preconditioner = .true.
logical, save      :: davidson_refresh_trial_vectors = .false.

! 
! variables for controlling i/o
!
character(len=16), save :: left_sigma_file_name     = 'LSIGMAV.sym'
character(len=16), save :: left_vectors_file_name   = 'LEVECTO.sym'
character(len=16), save :: left_subspace_file_name  = 'LSUBSPA.sym'

character(len=16), save :: right_sigma_file_name    = 'RSIGMAV.sym'
character(len=16), save :: right_vectors_file_name  = 'REVECTO.sym'
character(len=16), save :: right_subspace_file_name = 'RSUBSPA.sym'

integer, save :: right_sigma_file_unit    = 140
integer, save :: right_vectors_file_unit  = 141
integer, save :: right_subspace_file_unit = 142

integer, save :: left_sigma_file_unit     = 240
integer, save :: left_vectors_file_unit   = 241
integer, save :: left_subspace_file_unit  = 242


contains 

! this should never be used to update an instance, only to create new ones. it's better to have results
! for left and right separately if we determine these in two steps. but we must check whether an 
! element of the derived type has been allocated already... perhaps this should then be changed to a
! function returning a results type, and have the davidson driver return that to its caller
    subroutine davidson_create_results(results, N, nroots, values, left, right)
        type(talsh_davidson_results_tens_t), intent(inout) :: results
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


!       if (davidson_symmetric) then
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
        type(talsh_davidson_results_tens_t), intent(inout) :: results

!       if (associated(results%eValues)) deallocate(results%eValues)
!       if (associated(results%eVectorsR)) deallocate(results%eVectorsR)
!       if (associated(results%eVectorsL)) deallocate(results%eVectorsL)
    end subroutine


    subroutine talsh_davidson_setup(convergence_threshold, &
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

        call set_safe_minimum() 

        if (present(convergence_threshold)) then
            davidson_convergence_threshold = convergence_threshold
        else
            davidson_convergence_threshold = 1.0d-8 
        end if

        davidson_antisymme_threshold = davidson_convergence_threshold * 0.1


        if (present(energy_shift)) then
            davidson_eigenvalues_shift = energy_shift
            davidson_apply_eigenvalues_shift = .true.

            if (abs(davidson_eigenvalues_shift(1)).lt.davidson_eigenvalues_shift_threshold) &
               davidson_eigenvalues_shift(1) = 0.0d0 
 
            if (abs(davidson_eigenvalues_shift(2)).lt.davidson_eigenvalues_shift_threshold) &
               davidson_eigenvalues_shift(2) = 0.0d0 

            if ((abs(davidson_eigenvalues_shift(1)).lt.davidson_eigenvalues_shift_threshold).and. &
                (abs(davidson_eigenvalues_shift(2)).lt.davidson_eigenvalues_shift_threshold)) then
               davidson_apply_eigenvalues_shift = .false.
            end if
        else
            davidson_eigenvalues_shift = (/0.0d0, 0.0d0/) 
            davidson_apply_eigenvalues_shift = .false.
        end if

        if (present(max_subspace_size)) then
            davidson_max_subspace_size = max_subspace_size  
        else
            davidson_max_subspace_size = 2048
        end if

        if (present(max_iterations)) then
            davidson_max_iterations = max_iterations
        else
            davidson_max_iterations = 100
        end if

        if (present(solve_values)) then
            davidson_solve_values = solve_values
        else
            davidson_solve_values = .true. 
        end if

        if (present(solve_linear_sys)) then
            davidson_solve_linear_system = solve_linear_sys
            if (davidson_solve_linear_system) davidson_solve_values = .false.
        else
            davidson_solve_linear_system = .false. 
        end if
        if (davidson_solve_linear_system.and.davidson_solve_values) then
           call quit('inconsistent setup. cannot calculate eigenvalues and solve linear system')
        endif
 

        if (present(solve_left)) then
            davidson_solve_lhs = solve_left 
        else
            davidson_solve_lhs = .false. 
        end if

        if (present(solve_right)) then
            davidson_solve_rhs = solve_right
        else
            davidson_solve_rhs = .true. 
        end if

        if (present(overlap_sorting_scheme)) then
            davidson_overlap_sorting_scheme = overlap_sorting_scheme
        else
            davidson_overlap_sorting_scheme = 1 ! 1 = scheme calculating overlap matrix ; 2 = scheme not explicitly calculating overlap matrix
        end if

        if (present(overlap_sorting)) then
            davidson_overlap_sorting = overlap_sorting
        else
            davidson_overlap_sorting = .false.
        end if

        if (present(use_preconditioner)) then
            davidson_use_preconditioner = use_preconditioner
        else
            davidson_use_preconditioner = .true.
        end if

        if (present(verbose)) then
            davidson_verbose = verbose
            davidson_debug   = verbose
        else
            davidson_verbose = .false.
        end if

        if (present(symmetric)) then
            davidson_symmetric = symmetric
        else
            davidson_symmetric = .false. 
        end if

        if (present(complex_mode)) then
            davidson_complex_mode_with_reals = complex_mode 
        else
            davidson_complex_mode_with_reals = .false. 
        end if

        if (present(refresh_trial_rate)) then
            if (refresh_trial_rate.le.0) then
                davidson_refresh_trial_vectors = .false.
            else 
                davidson_refresh_trial_vectors = .true.
                if (refresh_trial_rate.le.davidson_maximum_trial_vector_refresh_rate) then
                    davidson_trial_vector_refresh_rate = davidson_maximum_trial_vector_refresh_rate 
                else
                    davidson_trial_vector_refresh_rate = refresh_trial_rate
                end if
            end if
        else
            davidson_refresh_trial_vectors = .false.
        end if

        if (present(unit_output)) then
            davidson_output_file_unit = unit_output
        end if

        if (present(save_subspaces)) then
           davidson_save_subspaces = save_subspaces
        else
           davidson_save_subspaces = .false. 
        end if

        if (present(save_sigmas)) then
           davidson_save_sigmas = save_sigmas
        else
           davidson_save_sigmas = .false. 
        end if

        if (present(save_results)) then
           davidson_save_results = save_results
        else
           davidson_save_results = .false. 
        end if

        if (present(unit_subspaces)) then
            davidson_subspaces_file_unit = unit_subspaces
        else
            davidson_subspaces_file_unit = 137 
        end if

        if (present(unit_sigmas)) then
            davidson_sigmas_file_unit = unit_sigmas
        else
            davidson_sigmas_file_unit = 138
        end if

        if (present(unit_results)) then
            davidson_results_file_unit = unit_results
        else
            davidson_results_file_unit = 138
        end if

        if (davidson_complex_mode_with_reals.and.davidson_overlap_sorting) then
            davidson_overlap_sorting_scheme = 2
        end if

        call davidson_show_configurations()
    end subroutine


    subroutine davidson_show_configurations

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "Configuration variables for matrix-free diagonalizer"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  output written to unit        : ", davidson_output_file_unit
        write(davidson_output_file_unit, *) "  verbose output                : ", davidson_verbose

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) " convergence control"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  convergence threshold         : ", davidson_convergence_threshold
        write(davidson_output_file_unit, *) "  maximum subspace size         : ", davidson_max_subspace_size
        write(davidson_output_file_unit, *) "  maximum number of iterations  : ", davidson_max_iterations
        write(davidson_output_file_unit, *) "  refresh trial vectors         : ", davidson_refresh_trial_vectors 
        if (davidson_refresh_trial_vectors) then
            write(davidson_output_file_unit, '(A,I4,A)') &
                         "     set to refresh vectors every ",davidson_trial_vector_refresh_rate," interations"
        end if

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) " restart/data storage"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  save subspaces                :", davidson_save_subspaces
        if (davidson_save_subspaces) &
        write(davidson_output_file_unit, *) "  subspaces written to unit     : ", davidson_subspaces_file_unit
        write(davidson_output_file_unit, *) "  save sigma vectors            : ", davidson_save_sigmas
        if (davidson_save_sigmas) &
        write(davidson_output_file_unit, *) "  sigma vectors written to unit : ", davidson_sigmas_file_unit
        write(davidson_output_file_unit, *) "  save results                  :", davidson_save_results
        if (davidson_save_results) &
        write(davidson_output_file_unit, *) "  results written to unit       : ", davidson_results_file_unit

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) " diagonalization characteristics"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  solve for right eigenvectors  : ", davidson_solve_rhs
        write(davidson_output_file_unit, *) "  solve for left eigenvectors   : ", davidson_solve_lhs
        write(davidson_output_file_unit, *) "  symmetric eigenproblem        : ", davidson_symmetric
        write(davidson_output_file_unit, *) "  complex mode (real variables) : ", davidson_complex_mode_with_reals
        write(davidson_output_file_unit, *) "  root following via overlap    : ", davidson_overlap_sorting
        if (davidson_overlap_sorting) then
            write(davidson_output_file_unit, *) "     overlap following scheme   : ", davidson_overlap_sorting_scheme
        end if
        write(davidson_output_file_unit, *) "  energy shift in eigval sorting: ",davidson_apply_eigenvalues_shift
!       if (davidson_apply_eigenvalues_shift) then
           write(davidson_output_file_unit, *) "     shift value, real  part    : ",davidson_eigenvalues_shift(1)
           write(davidson_output_file_unit, *) "     shift value, imag. part    : ",davidson_eigenvalues_shift(2)
!       end if

        write(davidson_output_file_unit, *) ""

    end subroutine


    subroutine set_safe_minimum()
        double precision dlamch
        safe_minimum = dlamch('S')
    end subroutine 


    subroutine talsh_davidson_destroy(results)
        type(talsh_davidson_results_tens_t), intent(inout) :: results
        call free_davidson_results(results)
    end subroutine


    subroutine update_trial_vectors(delta,workspace_t)
        type(talsh_davidson_workspace_tens_t) :: workspace_t
        logical, pointer, intent(inout) :: delta(:)

        integer :: ierr, i, j
        integer(INTD) :: rank
        logical :: debug = .false.

        ! rank=talsh_tensor_rank(workspace_t%trial1(1))
        rank=workspace_t%tens1_dims(1)
        ! if (rank.gt.2) call quit('rank not implemented in update_trial_vectors')

        if (debug) then
           write (6,*) ''
           write (6,*) 'debug >>> in update_trial_vectors '
           write (6,*) 'debug >>> nroots',workspace_t%nroots
           write (6,*) 'debug >>> Lprev ',workspace_t%Lprev
           write (6,*) 'debug >>> L     ',workspace_t%L
           write (6,*) ''
 
           do i = 1, workspace_t%L
              write (6,*) 'debug >>> in update_trial_vectors, before adding corrections, trial vector ',i
              call print_tensor(workspace_t%trial1(i),1.0d-18, "b1")
              call print_tensor(workspace_t%trial2(i),1.0d-18, "b2")
              write (6,*) ''
           end do
        end if

        workspace_t%Lprev = workspace_t%L

        if (debug) then
           write (6,*) 'debug >>> in update_trial_vectors, before processing correction vectors'
           write (6,*) 'debug >>> Lprev ',workspace_t%Lprev
           write (6,*) 'debug >>> L     ',workspace_t%L
           write (6,*) ''
        end if

        do i = 1, size(delta,1)
           if (delta(i)) then
              workspace_t%L = workspace_t%L + 1
              j = workspace_t%L
              if (debug) then
                 write (6,*) 'debug >>> in update_trial_vectors, processing correction vector ',i
                 write (6,*) 'debug >>> j             ',j
                 write (6,*) ''
              end if

              ierr=talsh_tensor_init(workspace_t%trial1(j), ZERO)
              ierr=talsh_tensor_init(workspace_t%trial2(j), ZERO)

              ierr=talsh_tensor_add("t(a,b)+=d(a,b)",        workspace_t%trial1(j),workspace_t%delta1(i))
              ierr=talsh_tensor_add("t(a,b,c,d)+=d(a,b,c,d)",workspace_t%trial2(j),workspace_t%delta2(i))
              
              if (debug) then
                 write (6,*) 'debug >>> correction vector ',i,' added to trial vector ',j
                 call print_tensor(workspace_t%trial1(j),1.0d-18, "new b1")
                 call print_tensor(workspace_t%trial2(j),1.0d-18, "new b2")
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
              call print_tensor(workspace_t%trial1(i),1.0d-18, "b1")
              call print_tensor(workspace_t%trial2(i),1.0d-18, "b2")
              write (6,*) ''
           end do
        end if

        call gram_schmidt_orthogonalization_tensors(workspace_t%Lprev, workspace_t%L, workspace_t%trial1, workspace_t%trial2,rank)

        ! if (debug) then
        !     write (6,*) 'debug >>> in update_trial_vectors, trial vectors after  GS orthonormalization'
        !     do i = 1, workspace_t%L
        !        write (6,*) 'debug >>> trial vector ',i
        !        call print_tensor(workspace_t%trial1(i),1.0d-18, "b1")
        !        call print_tensor(workspace_t%trial2(i),1.0d-18, "b2")
        !        print *, "norm of new trial", sqrt(tensor_norm2(workspace_t%trial1(i)) + tensor_norm2(workspace_t%trial2(i))/4)
        !        write (6,*) ''
        !     end do
        !  end if

        call double_check_trial_vector(workspace_t%Lprev, workspace_t%L, workspace_t%trial1, workspace_t%trial2,rank)

        if (debug) then
           write (6,*) 'debug >>> in update_trial_vectors, trial vectors after  GS orthonormalization '
           do i = 1, workspace_t%L
              write (6,*) 'debug >>> trial vector ',i
              call print_tensor(workspace_t%trial1(i),1.0d-18, "b1")
              call print_tensor(workspace_t%trial2(i),1.0d-18, "b2")
              if (rank.eq.1) then
                print *, "norm of new trial", sqrt(tensor_norm2(workspace_t%trial1(i)) + tensor_norm2(workspace_t%trial2(i))/2)
              else
                print *, "norm of new trial", sqrt(tensor_norm2(workspace_t%trial1(i)) + tensor_norm2(workspace_t%trial2(i))/4)
              end if
              write (6,*) ''
           end do
        end if
    end subroutine


    subroutine double_check_trial_vector(L,Lnew,trial1,trial2,rank)
        integer, intent(inout) :: L, Lnew
        type(talsh_tens_t), intent(inout) :: trial1(:), trial2(:)
        integer :: rank

        integer :: Lnew_tem, L_tem
        integer :: i, j, vi
        real(8) :: norm_r
        integer, allocatable :: valid_index(:)
        type(C_PTR) :: body_p, body_p1
        ! real(8) :: davidson_antisymme_threshold = 1.0d-6
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
            if (norm_r.gt.davidson_antisymme_threshold) then 
                L_tem = L_tem + 1
                valid_index(L_tem-L) = i
            end if 
        end do
        Lnew = L_tem
       
        if (debug_double_check) then
            print *, "valid index of vectors after checking", valid_index
        end if 

        j = 1
        ierr = talsh_tensor_dimensions(trial2(1),t2_rank,t2_dims)
        ierr = talsh_tensor_dimensions(trial1(1),t1_rank,t1_dims)

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

            ierr=talsh_tensor_get_body_access(trial2(vi),body_p,C8,int(0,C_INT),DEV_HOST)
            ierr=talsh_tensor_get_body_access(trial2(i),body_p1,C8,int(0,C_INT),DEV_HOST)
            ! ierr=talsh_init(tensor_position_you_waht_to_overwrite,ZERO)
            ! ierr=talsh_add(replacement+=vector_you_keep)
            call c_f_pointer(body_p,trial2_tens,t2_dims)
            call c_f_pointer(body_p1,trial2_tens_tem,t2_dims)
            trial2_tens_tem = trial2_tens
            
            ierr=talsh_tensor_get_body_access(trial1(vi),body_p,C8,int(0,C_INT),DEV_HOST)
            ierr=talsh_tensor_get_body_access(trial1(i),body_p1,C8,int(0,C_INT),DEV_HOST)
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
        deallocate(valid_index)

    end subroutine

    subroutine renormalize_vectors(v1,v2,rank)

       type(talsh_tens_t), intent(inout) :: v1, v2
       integer :: rank
       real :: norm_r
       complex(8) :: norm
    !    real(8) :: davidson_antisymme_threshold = 1.0d-6
       integer :: ierr
  

    !    rank=talsh_tensor_rank(v1)

       if (rank.eq.1) then
            norm_r = sqrt(tensor_norm2(v1) + tensor_norm2(v2)/2)
       else 
            norm_r = sqrt(tensor_norm2(v1) + tensor_norm2(v2)/4)
       end if 

    !    norm_r = sqrt(tensor_norm2(v1) + tensor_norm2(v2)/4) !XY 4 is for EE(ab.ij) but we need to take care for IP
       norm = dcmplx(norm_r, 0.0d0)

    !    print *, "norm in renormalize", norm_r
       if (norm_r.lt.davidson_antisymme_threshold) then ! to avoid numerical noisy 
            ierr=talsh_tensor_init(v1, ZERO)
            ierr=talsh_tensor_init(v2, ZERO)
       else
            call divide_tensor_elements_by_value(v1, norm)
            call divide_tensor_elements_by_value(v2, norm)
       end if

        ! if (norm_r .ge. safe_minimum ) then
        !     call divide_tensor_elements_by_value(v1, norm)
        !     call divide_tensor_elements_by_value(v2, norm)
        ! else
        !     call quit ('norm of tensors too small for division, quitting')
        ! end if

    end subroutine


    subroutine gram_schmidt_projection_tensors(factor, u1, u2, v1, v2, rank)
       complex(8), intent(inout) :: factor
       type(talsh_tens_t), intent(inout) :: u1, u2
       type(talsh_tens_t), intent(inout) :: v1, v2
       integer :: rank

       type(talsh_tens_t) :: result_tensor_num, result_tensor_den
       complex(8), pointer :: result_tens_num(:), result_tens_den(:)
       type(C_PTR):: body_p_num, body_p_den
       type(talsh_tens_t) :: one_tensor

       integer(C_INT)     :: one_dims(1)
       integer  :: ierr
       logical :: debug_GS_proj = .false.
       complex(8)    :: contrac_coeff

       one_dims(1) = 1

       ierr=talsh_tensor_construct(result_tensor_den,C8,one_dims(1:0),init_val=ZERO)
       ierr=talsh_tensor_get_body_access(result_tensor_den,body_p_den,C8,int(0,C_INT),DEV_HOST)
       call c_f_pointer(body_p_den,result_tens_den,one_dims)
       ierr=talsh_tensor_construct(result_tensor_num,C8,one_dims(1:0),init_val=ZERO)
       ierr=talsh_tensor_get_body_access(result_tensor_num,body_p_num,C8,int(0,C_INT),DEV_HOST)
       call c_f_pointer(body_p_num,result_tens_num,one_dims)

       ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

    !    rank=talsh_tensor_rank(v1)
       
       if (rank == 1) then
        contrac_coeff = ONE_HALF
       else
        contrac_coeff = ONE_QUARTER_C
       end if 

    !    if (rank.eq.1) then
    !       ierr=talsh_tensor_contract("num()+=u+(a)*v(a)",result_tensor_num,u1,v1)
    !       ierr=talsh_tensor_contract("den()+=u+(a)*u(a)",result_tensor_den,u1,u1)

    !       ierr=talsh_tensor_contract("num()+=u+(a,b,c)*v(a,b,c)",result_tensor_num,u2,v2,scale=ONE_HALF)
    !       ierr=talsh_tensor_contract("den()+=u+(a,b,c)*u(a,b,c)",result_tensor_den,u2,u2,scale=ONE_HALF)

    !    else if (rank.eq.2) then
    !    end if

       if (debug_GS_proj) then
           call print_tensor(u1, 1.0d-18, "u1 tensor")
           call print_tensor(v1, 1.0d-18, "v1 tensor")
           print *, "rank in projection:", rank
       end if 
       ierr=talsh_tensor_contract("num()+=u+(a,b)*v(a,b)",result_tensor_num,u1,v1)
       ierr=talsh_tensor_contract("den()+=u+(a,b)*u(a,b)",result_tensor_den,u1,u1)
       ierr=talsh_tensor_contract("num()+=u+(a,b,c,d)*v(a,b,c,d)",result_tensor_num,u2,v2,scale=contrac_coeff)
       ierr=talsh_tensor_contract("den()+=u+(a,b,c,d)*u(a,b,c,d)",result_tensor_den,u2,u2,scale=contrac_coeff)

       
       factor = result_tens_num(1) / result_tens_den(1)

       if (debug_GS_proj) then
          print *, "result_tens_den():", result_tens_den(1) 
          print *, "result_tens_num():", result_tens_num(1)
          print *, "factor: ", factor
       end if 

       ierr=talsh_tensor_destruct(result_tensor_num)
       ierr=talsh_tensor_destruct(result_tensor_den)
       ierr=talsh_tensor_destruct(one_tensor)

    end subroutine


    subroutine gram_schmidt_orthogonalization_tensors(L,Lnew,trial1,trial2,rank)
        integer, intent(in) :: L, Lnew
        type(talsh_tens_t), intent(inout) :: trial1(:), trial2(:)
        integer(INTD) :: rank

        real :: norm_r
        integer :: ierr, i, j
        complex(8) :: factor
        ! type(talsh_tens_t) :: one_tensor
        ! integer(C_INT)     :: one_dims(1)
        ! real(8) :: davidson_antisymme_threshold = 1.0d-6
        ! real(8) :: davidson_antisymme_threshold_num = 1.0d-12

        integer(C_INT)     :: one_dims(1)
        logical :: debug_GS_orthg = .false.



        ! one_dims(1) = 1
        ! ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
        ! rank=talsh_tensor_rank(trial1(1))
        ! if (rank.ne.2) call quit('rank not implemented in form_r_tensors')


        do j = L+1, Lnew
            if (debug_GS_orthg) then
                call print_tensor(trial1(j), 1.0d-18, "trial1(j) before ")
                call print_tensor(trial2(j), 1.0d-18, "trial2(j) before ")
            end if
            do i = 1, j-1
                if (debug_GS_orthg) then
                    print *, "I: ", i, "J: ", j
                    call print_tensor(trial1(i), 1.0d-10, "trial1(i)")
                    call print_tensor(trial2(i), 1.0d-10, "trial2(i)")
                end if 

                if (i.le.L) then    !The previous L vectors should be ok.
                    norm_r = 1.0
                else
                    norm_r = sqrt(tensor_norm2(trial1(i)) +  tensor_norm2(trial2(i))/4)
                end if
 
                if (norm_r.gt.davidson_antisymme_threshold) then                       
                    call gram_schmidt_projection_tensors(factor,trial1(i), trial2(i), trial1(j), trial2(j),rank)
                else
                    factor = ZERO ! complex(0.0D0,0.0D0)
                end if

                ! if (rank.eq.1) then
                !    ierr=talsh_tensor_add("t(a)+=p(a)",trial1(j),trial1(i),scale=-factor)
                !    ierr=talsh_tensor_add("t(a,b,c)+=p(a,b,c)",trial2(j),trial2(i),scale=-factor)
                ! else if (rank.eq.2) then
                ierr=talsh_tensor_add("t(a,b)+=p(a,b)",trial1(j),trial1(i),scale=-factor)
                ierr=talsh_tensor_add("t(a,b,c,d)+=p(a,b,c,d)",trial2(j),trial2(i),scale=-factor)
                ! end if

            end do

            call renormalize_vectors(trial1(j),trial2(j),rank)

         end do

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
                print '(4x,1e14.6,2x,SP,1e14.6,"i",$)', cmplx(vr,vi)
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

        if (davidson_complex_mode_with_reals) then
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
                write (davidson_output_file_unit,'(4x,a,i6,a,1e20.12,a,1e20.12,a)')  &
                      "Eigenvalue ",i," : ", vr," -",abs(vi)," i"
            else if (vi .gt. 0.0d0) then
                write (davidson_output_file_unit,'(4x,a,i6,a,1e20.12,a,1e20.12,a)')  &
                      "Eigenvalue ",i," : ", vr," +",abs(vi)," i"
            else
                write (davidson_output_file_unit,'(4x,a,i6,a,1e20.12)')  &
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
        if (davidson_debug) print *, "debug, evalues inside print_values_vectors_R"
        if (davidson_debug) print *, eValues(1,:)
        if (davidson_debug) print *, eValues(2,:)

        do i = 1, nroots_local
            vr = eValues(1,i)
            vi = eValues(2,i)
            if (abs(vi).lt.print_values_threshold) vi = 0.0d0
            if (abs(vr).lt.print_values_threshold) vr = 0.0d0

            if (vi .lt. 0.0d0) then
                write (davidson_output_file_unit,'(4x,a,i6,a,1e14.6,a,1e14.6,a)')  &
                      "Eigenvalue ",i," : ", vr," -",abs(vi)," i"
            else if (vi .gt. 0.0d0) then
                write (davidson_output_file_unit,'(4x,a,i6,a,1e14.6,a,1e14.6,a)')  &
                      "Eigenvalue ",i," : ", vr," +",abs(vi)," i"
            else
                write (davidson_output_file_unit,'(4x,a,i6,a,1e14.6)')  &
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
        if (davidson_complex_mode_with_reals) then
            vector_N_complex = vector_N / 2
            vector_norm = dznrm2(vector_N_complex, v, 1)
        end if
    end function

    subroutine form_r_tensors(k,eigenvectors_subspace,eigenvalues_subspace,workspace_t,prop_grad1, prop_grad2)
! we use delta{1,2,...} to store the values of r^k
        integer, intent(in) :: k
        complex(kind=8), intent(in)  :: eigenvectors_subspace(:,:)
        complex(kind=8), intent(in)  :: eigenvalues_subspace(:)
        type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
! aspg, added from (q)rsp
        ! complex(kind=8), pointer, contiguous, optional, intent(in) :: C_subspace(:)
        type(talsh_tens_t), intent(inout), optional :: prop_grad1, prop_grad2
! aspg, added from (q)rsp
        complex(kind=8) :: factor1, factor2, factor3
        integer(INTD) :: r1_rank 
        integer :: ierr, i
        logical :: debug_form_r = .false.
        type(talsh_tens_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
! we only need to know the rank of delta1 because the rank of r2 will always be r1_rank+2, that is: 
! r1 = 1 (1h for IP, 1p for EA) then r2 = 3 (2h1p for IP, 1h2p for EA)
! r1 = 2 (1h1p for EE, 2h for DIP, 2p for DEA) then r2 = 4 (2h2p for EE, 3h1p for DIP, 1h3p for DEA)
! etc
        r1_rank=talsh_tensor_rank(workspace_t%delta1(k))
        if (r1_rank.gt.2) call quit('rank not implemented in form_r_tensors')

!       ierr=talsh_tensor_init(workspace_t%delta1(k), ZERO)
!       ierr=talsh_tensor_init(workspace_t%delta2(k), ZERO)

!       if (debug_form_r) then
!          call davidson_print_vectors(eigenvectors_subspace,   label="eigenvectors_subspace in form_r")
!       end if

        do i = 1, workspace_t%L

           factor1 =  eigenvectors_subspace(i,k)
        !    factor1 = conjg(factor1)
           factor2 = -eigenvalues_subspace(k)*factor1
! aspg, added from (q)rsp
        !    if (present(C_subspace)) factor2 = -C_subspace(i)
        !    if (present(C_subspace)) factor3 = -C_subspace(i)
! aspg, added from (q)rsp

           if (debug_form_r) write(6,*) 'c(',i,',',k,')',factor1,'lambda*c',factor2

           if (r1_rank.eq.1) then
              
            !   call rotate_eigenvector(workspace_t%sigma1(i),workspace_t%sigma2(i),workspace_t)
              ierr=talsh_tensor_contract("r(a)+=s(a)*K()",workspace_t%delta1(k),workspace_t%sigma1(i),one_tensor,scale=factor1)
              ierr=talsh_tensor_contract("r(a)+=t(a)*K()",workspace_t%delta1(k),workspace_t%trial1(i),one_tensor,scale=factor2)
              ierr=talsh_tensor_contract("r(a,b,c)+=s(a,b,c)*K()",workspace_t%delta2(k),workspace_t%sigma2(i),&
                                & one_tensor,scale=factor1)
              ierr=talsh_tensor_contract("r(a,b,c)+=t(a,b,c)*K()",workspace_t%delta2(k),workspace_t%trial2(i),&
                                & one_tensor,scale=factor2)

            else if (r1_rank.eq.2) then
              ierr=talsh_tensor_add("r(a,b)+=s(a,b)",workspace_t%delta1(k),workspace_t%sigma1(i),scale=factor1)
              ierr=talsh_tensor_add("r(a,b)+=t(a,b)",workspace_t%delta1(k),workspace_t%trial1(i),scale=factor2)
              ierr=talsh_tensor_add("r(a,b,c,d)+=s(a,b,c,d)",workspace_t%delta2(k),workspace_t%sigma2(i),scale=factor1)
              ierr=talsh_tensor_add("r(a,b,c,d)+=t(a,b,c,d)",workspace_t%delta2(k),workspace_t%trial2(i),scale=factor2)

            !   call print_tensor(workspace_t%delta1(k), 1.0d-10, "delta2")

            !   if (present(C_subspace)) then
                ! ierr=talsh_tensor_add("r(a,b)+=t(a,b)",workspace_t%delta1(k),workspace_t%trial1(i),scale=factor3)
                ! ierr=talsh_tensor_add("r(a,b,c,d)+=t(a,b,c,d)",workspace_t%delta2(k),workspace_t%trial2(i),scale=factor3)
            !   end if 

           end if

!          if (debug_form_r) then
!             write (6,*) 'debug >>> in form_r_tensors, accumulated r1 root',k,' for L',i,' >>>'
!             call print_tensor_elements(workspace_t%delta1(k))
!             write (6,*) 'debug <<< in form_r_tensors, accumulated r1 for L',i,' <<<'
!             write (6,*) 'debug >>> in form_r_tensors, accumulated r2 root',k,' for L',i,' >>>'
!             call print_tensor_elements(workspace_t%delta2(k))
!             write (6,*) 'debug <<< in form_r_tensors, accumulated r2 for L',i,' <<<'
!          end if
       end do

       if (present(prop_grad1)) then
        ierr=talsh_tensor_contract("r(a,b)+=p(a,b)*K()",workspace_t%delta1(k),prop_grad1,&
        & one_tensor,scale=ONE)
        ierr=talsh_tensor_contract("r(a,b,c,d)+=p(a,b,c,d)*K()",workspace_t%delta2(k),prop_grad2,&
        & one_tensor,scale=ONE)
      end if 

! this should in the end work as what is below, but now we're doing a k at a time 
!      vectors_L = size(b,2)
!      vectors_N = size(b,1)
!      allocate(residue_a(vectors_N,vectors_L))
!      allocate(residue_b(vectors_N,vectors_L))
!      residue_a = 0.0d0
!      residue_b = 0.0d0
!!! this is residue_a(N,L) = sigma_vector(N,L)*c_small_M(L,L)
!      call zgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
!                 (1.0d0, 0.0d0), sigma_vector, vectors_N, c_small_M, vectors_L, &
!                 (0.0d0, 0.0d0), residue_a, vectors_N)
!!! this is residue_b(N,L) = sigma_vector(N,L)*b(L,L)
!      call zgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
!                 (1.0d0, 0.0d0), b, vectors_N, c_small_M, vectors_L,            &
!                 (0.0d0, 0.0d0), residue_b, vectors_N)
!      do k = 1, nroots
!         call zscal(vectors_N,e_small_M(:,k),residue_b(:,k),1)
!         r(:,k) = residue_a(:,k) - residue_b(:,k)
!      end do
!
! so
! r_a(:,j) = sigma(:,j)*c_small_M(j,j) for j = 1, L
! r_b(:,j) = b(:,j)*c_small_M(j,j) for j = 1, L
! r(:k) =  r_a(:,k) - lambda(k)*r_b(:,k) for k = 1, nroot 


       if (debug_form_r) then
          write (6,*) 'debug >>> in form_r_tensors, r1 root',k,' >>>'
          call print_tensor_elements(workspace_t%delta1(k))
          write (6,*) 'debug <<< in form_r_tensors, r1 <<<'
          write (6,*) 'debug >>> in form_r_tensors, r2 root',k,' >>>'
          call print_tensor_elements(workspace_t%delta2(k))
          write (6,*) 'debug <<< in form_r_tensors, r1 <<<'
       end if

       ierr=talsh_tensor_destruct(one_tensor)

    end subroutine    



    function correction_to_vectors(output_file_unit, e_small_M, c_small_M, workspace_t, prop_grad1, prop_grad2, correct_all)
        integer, intent(in)       :: output_file_unit
        type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
        complex(kind=8), intent(in)  :: e_small_M(:)
        complex(kind=8), intent(in)  :: c_small_M(:,:)
! aspg, added from (q)rsp
        ! complex(kind=8), pointer, contiguous, optional, intent(in) ::  C_small(:)
        type(talsh_tens_t), intent(inout), optional :: prop_grad1, prop_grad2
! aspg, added from (q)rsp
        logical, optional, intent(in) :: correct_all

        logical, pointer :: correction_to_vectors(:)

        integer :: ierr, k, unconverged_roots, nroots
        integer :: rank
        real(kind=8) :: norm_r = 0.0d0
        logical :: correct_all_local
        logical :: debug_correction_vectors = .false.

        if (present(correct_all)) then
            correct_all_local = correct_all
        else
            correct_all_local = .false. 
        end if

        ! if (present(C_small)) then
        if (present(prop_grad1)) then
            nroots = 1
        else
            nroots = workspace_t%nroots
        end if

        allocate(correction_to_vectors(nroots))
        correction_to_vectors = .true.

        unconverged_roots = 0
        ! rank = talsh_tensor_rank(workspace_t%delta1(1))
        rank = workspace_t%tens1_dims(1) 
        do k = 1, nroots
           ierr=talsh_tensor_init(workspace_t%delta1(k), ZERO)
           if (debug_correction_vectors) print *,'debug, init workspace_t%delta1(k) to zero, ierr',ierr
           ierr=talsh_tensor_init(workspace_t%delta2(k), ZERO)
           if (debug_correction_vectors) print *,'debug, init workspace_t%delta2(k) to zero, ierr',ierr
        end do
        do k = 1, nroots
! instead of having a new variable to store r(:,k), we will do collect the results on workspace_t%delta{1,2,...} 
        !    if (present(C_small)) then
           if (present(prop_grad1)) then
            !   call form_r_tensors(k,c_small_M,e_small_M,workspace_t,C_small)
              call form_r_tensors(k,c_small_M,e_small_M,workspace_t,prop_grad1,prop_grad2)
           else
              call form_r_tensors(k,c_small_M,e_small_M,workspace_t)
           end if
            
           if (rank .eq.1) then
              norm_r = sqrt(tensor_norm2(workspace_t%delta1(k)) +  tensor_norm2(workspace_t%delta2(k))/2)
           else 
              norm_r = sqrt(tensor_norm2(workspace_t%delta1(k)) +  tensor_norm2(workspace_t%delta2(k))/4)
           end if

! get the norm of the whole intermediate r(:,k) tensor, which means to calculate it over the different excitation ranks 
           if (debug_correction_vectors) then
              print *,'debug, norm of delta1(',k,'):',tensor_norm2(workspace_t%delta1(k))
              print *,'debug, norm of delta2(',k,'):',tensor_norm2(workspace_t%delta2(k))
           end if

           if (norm_r .lt. davidson_convergence_threshold) then 
              write(output_file_unit, "(4x,a,i6,a,1e14.6,4x,a)") "Root ",k,"  ||q_k,M|| ", norm_r, " <<< converged ! "
              correction_to_vectors(k) = .false.
           else
              write(output_file_unit, "(4x,a,i6,a,1e14.6)")  "Root ",k,"  ||q_k,M|| ", norm_r
              correction_to_vectors(k) = .true.
              unconverged_roots = unconverged_roots + 1
           end if
        end do
        write (davidson_output_file_unit,*) " "

        if (any(correction_to_vectors)) then
           if (correct_all_local) then
               unconverged_roots = nroots
               correction_to_vectors = .true. 
           end if
! forms delta (the correction vectors which will serve to create new trial vectors
! note here that, as seen above, workspace_t%delta{1,2,..} have been used to accumulate r(i,k) 
           do k = 1, nroots
              if (correction_to_vectors(k)) then
                 if (debug_correction_vectors) then
                    write (6,*) 'debug >>> in correction_to_vectors, precond1 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%precond1(k))
                    call print_tensor(workspace_t%precond1(k), 1.0d-10, "workspace_t%precond1(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, precond1 <<<'
                    write (6,*) 'debug >>> in correction_to_vectors, delta1 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta1(k))
                    call print_tensor(workspace_t%delta1(k), 1.0d-10, "workspace_t%delta1(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, delta1 <<<'
                 end if

                 call multiply_tensor_elements(workspace_t%delta1(k), workspace_t%precond1(k))

                 if (debug_correction_vectors) then
                    write (6,*) 'debug >>> in correction_to_vectors, result delta1 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta1(k))
                    call print_tensor(workspace_t%delta1(k), 1.0d-10, "workspace_t%delta1(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, result delta1 <<<'

                    write (6,*) 'debug >>> in correction_to_vectors, precond2 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%precond2(k)) 
                    call print_tensor(workspace_t%precond2(k), 1.0d-10, "workspace_t%precond2(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, precond2 <<<'
                    write (6,*) 'debug >>> in correction_to_vectors, delta2 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta2(k))
                    call print_tensor(workspace_t%delta2(k), 1.0d-10, "workspace_t%delta2(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, delta2 <<<'
                 end if

                !  call print_tensor(workspace_t%delta2(1), 1.0d-10, "delta2 before precond")
                !  call print_tensor(workspace_t%precond2(1), 1.0d-10, "precond2 before precond")
                 call multiply_tensor_elements(workspace_t%delta2(k), workspace_t%precond2(k))
                !  call print_tensor(workspace_t%delta2(1), 1.0d-10, "delta2 after precond")
                !  call print_tensor(workspace_t%precond2(1), 1.0d-10, "precond2 after precond")

                 if (debug_correction_vectors) then
                    write (6,*) 'debug >>> in correction_to_vectors, result delta2 root',k,' >>>'
                    ! call print_tensor_elements(workspace_t%delta2(k))
                    call print_tensor(workspace_t%delta2(k), 1.0d-10, "workspace_t%delta2(k)")
                    write (6,*) 'debug <<< in correction_to_vectors, result delta2 <<<'
                 end if

                !  norm_r = tensor_norm2(workspace_t%delta1(k)) +  tensor_norm2(workspace_t%delta2(k))
                !  print *,'debug, norm of delta+delta2():',norm_r

              else
                 ierr=talsh_tensor_init(workspace_t%delta1(k), ZERO)
                 ierr=talsh_tensor_init(workspace_t%delta2(k), ZERO)
              end if

              if (debug_correction_vectors) then
                 write (6,*) 'debug >>> in correction_to_vectors, delta1 root',k,' >>>'
                !  call print_tensor_elements(workspace_t%delta1(k))
                 call print_tensor(workspace_t%delta1(k), 1.0d-10, "workspace_t%delta1(k)")
                 write (6,*) 'debug <<< in correction_to_vectors, delta1 <<<'
                 write (6,*) 'debug >>> in correction_to_vectors, delta2 root',k,' >>>'
                !  call print_tensor_elements(workspace_t%delta2(k))
                 call print_tensor(workspace_t%delta2(k), 1.0d-10, "workspace_t%delta2(k)")
                 write (6,*) 'debug <<< in correction_to_vectors, delta2 <<<'
              end if

           end do
        else
            deallocate(correction_to_vectors)
            correction_to_vectors => NULL()
        end if

    end function

    subroutine form_preconditioner(workspace_t, eigenvalues)
        type(talsh_davidson_workspace_tens_t) :: workspace_t
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
        real(kind=8) :: small_number = 1.0d-99

! careful with dimension 1 of eigenvalues, 
        nroots    = workspace_t%nroots

! whether or not we want to precondition, we start by setting up the tensor to one
        do iroot = 1, nroots
           ierr=talsh_tensor_init(workspace_t%precond1(iroot), ONE)   
           ierr=talsh_tensor_init(workspace_t%precond2(iroot), ONE)
        end do

        if (davidson_use_preconditioner) then
! we want to precondition
! first prepare to access the rank-1 and rank-2 diagonal of Hbar 
           ierr = talsh_tensor_get_body_access(workspace_t%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
           ierr = talsh_tensor_get_body_access(workspace_t%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

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

              ierr = talsh_tensor_get_body_access(workspace_t%precond1(iroot),preco1_body_p,C8,0,DEV_HOST)
              ierr = talsh_tensor_get_body_access(workspace_t%precond2(iroot),preco2_body_p,C8,0,DEV_HOST)

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
                    if (dabs(x2py2).gt.safe_minimum) then
                       proj_r =  proj_r / x2py2
                       proj_i = -proj_i / x2py2
                       preco_1d_p(ja) = cmplx(proj_r,proj_i)
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
                          if (dabs(x2py2).gt.safe_minimum) then
                             proj_r =  proj_r / x2py2
                             proj_i = -proj_i / x2py2
                             preco_3d_p(ja,jb,jc) = cmplx(proj_r,proj_i)
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
                          write (6,*) 'debug safe minimum',safe_minimum
                          write (6,*) 'debug small number',small_number
                          write (6,*) 'debug for root ',iroot
                          write (6,*) 'debug E(',eigenvalue_r,',',eigenvalue_i,')'
                          write (6,*) 'debug Hii',real(hdiag_2d_p(ja,jb)),aimag(hdiag_2d_p(ja,jb))
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
                       if (dabs(x2py2).gt.safe_minimum) then
                          proj_r =  proj_r / x2py2
                          proj_i = -proj_i / x2py2
                          preco_2d_p(ja,jb) = cmplx(proj_r,proj_i)
                       end if
                    end do
                 end do

                 do jd = 1, workspace_t%tens2_dims(4)
                    do jc = 1, workspace_t%tens2_dims(3)
                       do jb = 1, workspace_t%tens2_dims(2)
                          do ja = 1, workspace_t%tens2_dims(1)
                             proj_r = eigenvalue_r - real(hdiag_4d_p(ja,jb,jc,jd))
                             proj_i = eigenvalue_i - aimag(hdiag_4d_p(ja,jb,jc,jd))
                       if (preconditioner_debug) then
                          write (6,*) 'debug safe minimum',safe_minimum
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
                             if (dabs(x2py2).gt.safe_minimum) then
                                proj_r =  proj_r / x2py2
                                proj_i = -proj_i / x2py2
                                preco_4d_p(ja,jb,jc,jd) = cmplx(proj_r,proj_i)
                             end if
                          end do
                       end do
                    end do
                 end do

              end if

              if (preconditioner_debug) then
                 write (6,*) 'debug >>> precon1_tens  root',iroot,' >>>'
                 call print_tensor(workspace_t%precond1(iroot), 1.0d-18, "precond1")
                 write (6,*) 'debug <<< precon1_tens <<<'
                 write (6,*) 'debug >>> precon2_tens  root',iroot,' >>>'
                 call print_tensor(workspace_t%precond2(iroot), 1.0d-18, "precond2")
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
                    write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                                 'i:',i,' j:',j,' kr',kr,' ki:',ki,'  vr:',vr,'  vi:',vi 

                    vt(j,kr) =  vr
                    vt(j,ki) = -vi

                     write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1e14.6,a,1e14.6,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           ' v(',i,',',j,')=(',vr,',',vi,') ->  vt(',j,',',i,')=(',vt(j,kr),',',vt(j,ki),')'
                else if (complex_dimension_v.eq.2) then

                    kr = L_step*(j-1) + 1
                    ki = kr + 1

                    vr = v(i, kr)
                    vi = v(i, ki)
                     write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           'i:',i,' j:',j,' kr',kr,' ki:',ki,'  vr:',vr,'  vi:',vi

                    vt(kr, i) =  vr
                    vt(ki, i) = -vi

                     write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1e14.6,a,1e14.6,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           ' v(',i,',',j,')=(',vr,',',vi,') ->  vt(',j,',',i,')=(',vt(kr,i),',',vt(ki,i),')'

                end if
            end do
        end do

    end subroutine

    subroutine talsh_solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left)
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
    subroutine talsh_solve_Linear_systems(left, A, B_input, Solution)
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
        Solution(:,NRHS) = X(:,NRHS)

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

    end subroutine talsh_solve_Linear_systems
! aspg, added from (q)rsp

    subroutine get_tensor_linear_systems(t_tensor,root_vect,workspace_t)
      
        type(talsh_comm_tens_t)              :: t_tensor
        complex(kind=8), pointer, contiguous :: root_vect(:,:)
        type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
        

        type(talsh_tens_t)          :: one_tensor
        integer(C_INT)              :: ierr
        integer                     :: iroot, N
        integer(C_INT)              :: one_dims(1)
        complex(kind=8)             :: factor

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

        N = size(root_vect,1)
        ! print *, "N in get_tensor_linear_systems", N


        ! print *, "dim1", size(root_vect,1)
        ! print *, "dim2", size(root_vect,2)
        ! print *, "roor_vect", root_vect(:,1)
        ! call print_tensor(t_tensor%t1, 1.0d-10, "t_tensor%t1")
        

         do iroot=1, N
            factor = root_vect(iroot,1)
            ierr=talsh_tensor_contract("t(a,i)+=B1(a,i)*X()",t_tensor%t1,workspace_t%trial1(iroot),one_tensor,scale=factor)
            ierr=talsh_tensor_contract("t(a,b,i,j)+=B2(a,b,i,j)*X()",t_tensor%t2,workspace_t%trial2(iroot),one_tensor,scale=factor)
            
         end do
        !  call print_tensor(t_tensor%t1, 1.0d-10, "t_tensor%t1")
        !  call print_tensor(t_tensor%t2, 1.0d-10, "t_tensor%t2")


    end subroutine get_tensor_linear_systems



    subroutine talsh_davidson_workspace_destroy(workspace_t)
        type(talsh_davidson_workspace_tens_t),  intent(inout) :: workspace_t
        integer :: ierr, i

        do i = 1, workspace_t%Lmax
           ierr = talsh_tensor_destruct(workspace_t%trial1(i))
           ierr = talsh_tensor_destruct(workspace_t%trial2(i))
           ierr = talsh_tensor_destruct(workspace_t%sigma1(i))
           ierr = talsh_tensor_destruct(workspace_t%sigma2(i))
        end do
        deallocate(workspace_t%trial1)
        deallocate(workspace_t%trial2)
        deallocate(workspace_t%sigma1)
        deallocate(workspace_t%sigma2)

        do i = 1, workspace_t%nroots
           ierr = talsh_tensor_destruct(workspace_t%delta1(i))
           ierr = talsh_tensor_destruct(workspace_t%delta2(i))
           ierr = talsh_tensor_destruct(workspace_t%precond1(i))
           ierr = talsh_tensor_destruct(workspace_t%precond2(i))
        end do
        deallocate(workspace_t%delta1)
        deallocate(workspace_t%delta2)
        deallocate(workspace_t%precond1)
        deallocate(workspace_t%precond2)

        deallocate(workspace_t%tens1_dims)
        deallocate(workspace_t%tens2_dims)

        ierr = talsh_tensor_destruct(workspace_t%Hdiag1)
        ierr = talsh_tensor_destruct(workspace_t%Hdiag2)

        if (associated(workspace_t%N)) deallocate(workspace_t%N)
        if (associated(workspace_t%M)) deallocate(workspace_t%M)

        write(davidson_output_file_unit, "(4x,a)") 'leaving talsh_davidson_workspace_destroy'
    end subroutine 

    subroutine talsh_davidson_workspace_create(workspace_t, nroots, nocc, nvirt, flavor, excitation_level, left_right)
        type(talsh_davidson_workspace_tens_t),  intent(inout) :: workspace_t
        character*2 :: flavor
        integer, intent(in) :: nocc, nvirt, nroots, excitation_level
        character*1, OPTIONAL                              :: left_right
        logical :: workspace_create_debug = .false.

        integer :: M_dims(2) 
        integer :: ierr, i
        integer :: n_continuum_occ, n_continuum_virt
        integer :: excitation_manifold_size

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
              workspace_t%N(2) = nvirt*n_continuum_virt*nocc*nocc
            !   workspace_t%N(2) = nvirt*n_continuum_virt*(nocc*(nocc-1))/2

           case("EA")
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
                print *, "come"
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

        if (nroots.lt.0) then
           print *, 'INFO: number of roots reset to the dimension of Hbar'
           workspace_t%nroots = excitation_manifold_size 
           workspace_t%Lmax   = excitation_manifold_size
        else
           if (nroots.gt.excitation_manifold_size) then
              print *, 'INFO: number of roots > dimension of Hbar, resetting it'
           end if
           workspace_t%nroots = min(nroots,excitation_manifold_size)
           print *, "excitation_manifold_size:", excitation_manifold_size
           print *, "davidson_max_subspace_size:", davidson_max_subspace_size
           workspace_t%Lmax   = min(davidson_max_subspace_size,excitation_manifold_size)
        end if
        workspace_t%L      = workspace_t%nroots
        workspace_t%Lprev  = workspace_t%nroots

        allocate(workspace_t%trial1(workspace_t%Lmax))
        allocate(workspace_t%trial2(workspace_t%Lmax))

        allocate(workspace_t%sigma1(workspace_t%Lmax))
        allocate(workspace_t%sigma2(workspace_t%Lmax))

        allocate(workspace_t%delta1(workspace_t%nroots))
        allocate(workspace_t%delta2(workspace_t%nroots))

        allocate(workspace_t%precond1(workspace_t%nroots))
        allocate(workspace_t%precond2(workspace_t%nroots))

!       M_dims     = (/ workspace_t%Lmax, workspace_t%Lmax /)
!       ierr=talsh_tensor_construct(workspace_t%M,C8,M_dims,init_val=ZERO)
        allocate(workspace_t%M(workspace_t%Lmax,workspace_t%Lmax))
        workspace_t%M = ZERO

        allocate(workspace_t%C(workspace_t%Lmax))
        workspace_t%C = ZERO
! singles
        do i = 1, workspace_t%Lmax
           ierr = talsh_tensor_construct(workspace_t%trial1(i),C8,workspace_t%tens1_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace_t%sigma1(i),C8,workspace_t%tens1_dims,init_val=ZERO)
        end do
        do i = 1, workspace_t%nroots
           ierr = talsh_tensor_construct(workspace_t%delta1(i),C8,workspace_t%tens1_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace_t%precond1(i),C8,workspace_t%tens1_dims,init_val=ZERO)
        end do
        ierr=talsh_tensor_construct(workspace_t%Hdiag1,C8,workspace_t%tens1_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(workspace_t%Norm_prop_grad,C8,workspace_t%tens1_dims,init_val=ZERO)

! doubles
        do i = 1, workspace_t%Lmax
           ierr=talsh_tensor_construct(workspace_t%trial2(i),C8,workspace_t%tens2_dims,init_val=ZERO)
           ierr=talsh_tensor_construct(workspace_t%sigma2(i),C8,workspace_t%tens2_dims,init_val=ZERO)
        end do
        do i = 1, workspace_t%nroots
           ierr = talsh_tensor_construct(workspace_t%delta2(i),C8,workspace_t%tens2_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace_t%precond2(i),C8,workspace_t%tens2_dims,init_val=ZERO)
        end do
        ierr=talsh_tensor_construct(workspace_t%Hdiag2,C8,workspace_t%tens2_dims,init_val=ZERO)

        workspace_t%Nmax = workspace_t%N(1) + workspace_t%N(2)


! debug statements
!       if (workspace_create_debug) then
           write(davidson_output_file_unit, "(4x,a)")  'in talsh_davidson_workspace_create'

           write(davidson_output_file_unit, "(4x,a,i6,2x,a,2x,i6)") &
     &  'created tensor for trial vectors, dimensions :',workspace_t%Nmax,"x",workspace_t%Lmax


        if (workspace_create_debug) then
           do i = 1, workspace_t%Lmax
              write(davidson_output_file_unit, "(4x,a,i4)") 'info on workspace_t%trial1, element',i
              call talsh_tensor_print_info(workspace_t%trial1(i))
              write(davidson_output_file_unit, "(4x,a,i4)") 'info on workspace_t%trial2, element',i
              call talsh_tensor_print_info(workspace_t%trial2(i))
           end do

           write(davidson_output_file_unit, "(4x,a)") 'info on workspace_t%Hdiag1' 
           call talsh_tensor_print_info(workspace_t%Hdiag1)
           write(davidson_output_file_unit, "(4x,a)") 'info on workspace_t%Hdiag2' 
           call talsh_tensor_print_info(workspace_t%Hdiag2)

           write(davidson_output_file_unit, "(4x,a,i6,2x,a,2x,i6)") &
     &  'created tensor for sigma vectors, dimensions :',workspace_t%Nmax,"x",workspace_t%Lmax

           write(davidson_output_file_unit, "(4x,a,i6,2x,a,2x,i6)") &
     &  'created tensor for correction vectors, dimensions :',workspace_t%Nmax,"x",workspace_t%nroots
      end if

        
   end subroutine 


   subroutine create_sigma(workspace_t,int_t,interm_t,right,left)
      use talsh_sigma_eomccsd
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      logical, intent(in) :: left,right
      integer :: ierr, iroot
      
      logical :: debug_sigma = .false., debug_print_before_calculating=.false.

    !   do iroot = (workspace_t%L - workspace_t%nroots + 1), workspace_t%L
      do iroot = workspace_t%Lprev, workspace_t%L

         if (debug_sigma) then
            write (6,*) 'root ',iroot
            if (debug_print_before_calculating) then
            write (6,*) 'debug >>> t_tens before sigma, root',iroot,' >>>'
!           call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
            print "(2x,a6,2x,a6,2x,a6,2x,a62x,a6)", "idx", "p  ", "q  ","r  ","s  "
            call print_tensor(workspace_t%trial1(iroot), 0.0d-18, "trial1")
            call print_tensor(workspace_t%trial2(iroot), 0.0d-18, "trial2")
            write (6,*) 'debug <<< t_tens before sigma <<<'

            write (6,*) 'debug >>> s_tens before sigma, root',iroot,' >>>'
!           call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
            print "(2x,a6,2x,a6,2x,a6,2x,a62x,a6)", "idx", "p  ", "q  ","r  ","s  "
            call print_tensor(workspace_t%sigma1(iroot), 0.0d-18, "sigma1")
            call print_tensor(workspace_t%sigma2(iroot), 0.0d-18, "sigma2")
            write (6,*) 'debug <<< s_tens before sigma <<<'
            else
               write (6,*) 'debug <<< silenced printout prior to calculation >>>'  
            end if ! silencing printout prior to calculated vectors
         end if
! debug info     
! aspg
         if (right.and..not.left) then
            if (associated(interm_t%sigma_rhs_p)) then
               call interm_t%sigma_rhs_p(int_t,interm_t, &
   &                         workspace_t%trial1(iroot),workspace_t%trial2(iroot),&
   &                         workspace_t%sigma1(iroot),workspace_t%sigma2(iroot))
            else
               stop 'RHS sigma vector routine has not been defined prior to calling the solver'
            end if
         else if (left.and..not.right) then
            if (associated(interm_t%sigma_lhs_p)) then
               call interm_t%sigma_lhs_p(int_t,interm_t, &
   &                         workspace_t%trial1(iroot),workspace_t%trial2(iroot),&
   &                         workspace_t%sigma1(iroot),workspace_t%sigma2(iroot))
            else
               stop 'LHS sigma vector routine has not been defined prior to calling the solver'
            end if
         else
               stop  'Calculation of both left and right sigma not implemented'
         end if
! aspg
! debug info     
         if (debug_sigma) then
            write (6,*) 'debug >>> t_tens, root',iroot,' >>>'
            print "(2x,a6,2x,a6,2x,a6,2x,a62x,a6)", "idx", "p  ", "q  ","r  ","s  "
            call print_tensor(workspace_t%trial1(iroot), 1.0d-18, "trial1")
            call print_tensor(workspace_t%trial2(iroot), 1.0d-18, "trial2")
            write (6,*) 'debug <<< t_tens <<<'

            write (6,*) 'debug >>> s_tens, root',iroot,' >>>'
            print "(2x,a6,2x,a6,2x,a6,2x,a62x,a6)", "idx", "p  ", "q  ","r  ","s  "
            call print_tensor(workspace_t%sigma1(iroot), 1.0d-18, "sigma1")
            call print_tensor(workspace_t%sigma2(iroot), 1.0d-18, "sigma2")
            write (6,*) 'debug <<< s_tens <<<'
         end if
! debug info     
! aspg

      end do

!     call quit('bb')

   end subroutine

! aspg, added from (q)rsp
   subroutine create_subspace_property_gradient(C_small,workspace_t,prop_grad1,prop_grad2,right,left)
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t ! here we will keep the M(Lmax,Lmax)
      complex(kind=8), pointer, contiguous, intent(inout) :: C_small(:)   ! this will be a copy of B, B_small(L)
      type(talsh_tens_t), intent(inout) :: prop_grad1, prop_grad2
      logical, intent(in) :: right, left
       
      integer :: mij_dims(2), mij_rank
      integer(INTD)      :: one_dims(1)
      type(talsh_tens_t) :: result_tensor
      complex(kind=8), pointer ::  ci_result_tens(:)
      type(C_PTR):: ci_body_p

      integer :: m_dims(2), m_rank
      complex(kind=8), pointer, contiguous:: m_tens(:,:)
      type(C_PTR):: m_body_p
      integer :: ierr, i, j
      logical :: debug_subspace_matrix = .false.

      if (associated(C_small)) deallocate(C_small)
      allocate(C_small(workspace_t%L))
      C_small = ZERO

!     Initialize scalars that are to be used as tensors in contractions
      one_dims(1) = 1
!     ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
      ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
      ierr=talsh_tensor_get_body_access(result_tensor,ci_body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(ci_body_p,ci_result_tens,one_dims)
      

! capital letters here are meant to be the roots
      if (debug_subspace_matrix) then
         write (6,*) 'foobar, before setup M'
         write (6,*) 'foobar, workspace_t%L',workspace_t%L
         write (6,*) 'foobar, workspace_t%Lprev',workspace_t%Lprev
      end if
      do i = workspace_t%Lprev, workspace_t%L
        if (debug_subspace_matrix) then
            call print_tensor(workspace_t%trial1(i), 1.0d-18, "trial1")
            call print_tensor(prop_grad1, 1.0d-18, "prop_grad1")
            call print_tensor(workspace_t%trial2(i), 1.0d-18, "trial2")
            call print_tensor(prop_grad2, 1.0d-18, "prop_grad2")
        end if 
! this will be a reduction to a scalar at each contractions, and we then will build the 
         if (workspace_t%N(1).ne.0) then
            if (right.and..not.left) then
               ierr=talsh_tensor_contract("CI()+=Bdag+(p,q)*C(p,q)",             &
   &                                                       result_tensor,        &
   &                                                       workspace_t%trial1(i),&
   &                                                       prop_grad1,scale=ONE)
            else if (left.and..not.right) then
               ierr=talsh_tensor_contract("CI()+=C+(p,q)*B(p,q)",                &
   &                                                       result_tensor,        &
&                                                          prop_grad1,           &
   &                                                       workspace_t%trial1(i))
            end if
         end if ! workspace_t%N(1).ne.0

         if (workspace_t%N(2).ne.0) then
            if (right.and..not.left) then
               ierr=talsh_tensor_contract("CI()+=Bdag+(p,q,r,s)*C(p,q,r,s)",     &
   &                                                       result_tensor,        &
   &                                                       workspace_t%trial2(i),&
   &                                                       prop_grad2,           &
                                                           scale=ONE_QUARTER_C)
            else if (left.and..not.right) then
               ierr=talsh_tensor_contract("CI()+=C+(p,q,r,s)*B(p,q,r,s)",       &
   &                                                       result_tensor,       &
   &                                                       prop_grad2,          &
                                                           workspace_t%trial2(i),&
                                                           scale=ONE_QUARTER_C)
            end if
         end if

         if (debug_subspace_matrix) write (6,*) 'i, c_i:', i, ci_result_tens(1)
         workspace_t%C(i) = ci_result_tens(1)*(-1)
        !  workspace_t%C(i) = ci_result_tens(1)
         ci_result_tens(1) = ZERO
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

   subroutine create_subspace_matrix(M_small,workspace_t,right,left,frequency)
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t ! here we will keep the M(Lmax,Lmax)
      complex(kind=8), pointer, contiguous, intent(inout) :: M_small(:,:)    ! this will be a copy of M, M_small(L,L)
      logical, intent(in) :: right, left
      complex(kind=8), intent(in), optional :: frequency  

      integer :: mij_dims(2), mij_rank
      integer(INTD)      :: one_dims(1)
      type(talsh_tens_t) :: result_tensor
      complex(8), pointer :: mij_result_tens(:)
      type(C_PTR):: mij_body_p

      integer :: m_dims(2), m_rank
      complex(8), pointer, contiguous:: m_tens(:,:)
      type(C_PTR):: m_body_p
      integer :: ierr, i, j
      complex(8) :: contrac_coeff
      logical :: debug_subspace_matrix = .false.

!! we need to fish out the small matrix from workspace_t%M, which is the slice LxL of the full tensor
!! we also need to conjugate transpose the trial vector 
      if (associated(M_small)) deallocate(M_small)
      allocate(M_small(workspace_t%L,workspace_t%L))
      M_small = ZERO 

!     Initialize scalars that are to be used as tensors in contractions
      one_dims(1) = 1
      ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
    !   print *, "ierr construct small M:", ierr 
      ierr=talsh_tensor_get_body_access(result_tensor,mij_body_p,C8,int(0,C_INT),DEV_HOST)
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
      print *, "workspce_tensor:", workspace_t%tens1_dims
      if (workspace_t%tens1_dims(1)==1) then  !To distinguish IP and EE
        contrac_coeff = ONE_HALF   !It's IP
      else 
        contrac_coeff = ONE_QUARTER_C  !It's EE
      end if

      do i=1, workspace_t%L
        !  do j = workspace_t%Lprev, workspace_t%L
        do j =1, workspace_t%L
! this will be a reduction to a scalar at each contractions, and we then will build the 
            if (workspace_t%N(1).ne.0) then
               if (right.and..not.left) then 
                  ierr=talsh_tensor_contract("MIJ()+=Bdag+(p,q)*S(p,q)",          &
   &                                                        result_tensor,        &
   &                                                        workspace_t%trial1(i),&
   &                                                        workspace_t%sigma1(j))
               else if (left.and..not.right) then 
                  ierr=talsh_tensor_contract("MIJ()+=S+(p,q)*B(p,q)",              &
   &                                                        result_tensor,        &
   &                                                        workspace_t%sigma1(i),&
   &                                                        workspace_t%trial1(j))
               end if
            end if

            if (workspace_t%N(2).ne.0) then
               if (right.and..not.left) then
                  ierr=talsh_tensor_contract("MIJ()+=Bdag+(p,q,r,s)*S(p,q,r,s)",  &
   &                                                        result_tensor,        &
   &                                                        workspace_t%trial2(i),&
   &                                                        workspace_t%sigma2(j),&
   &                                                        scale=contrac_coeff)
               else if (left.and..not.right) then
                  ierr=talsh_tensor_contract("MIJ()+=S+(p,q,r,s)*B(p,q,r,s)",      &
   &                                                        result_tensor,        &
   &                                                        workspace_t%sigma2(i),&
   &                                                        workspace_t%trial2(j),&
   &                                                        scale=contrac_coeff)
               end if
            end if

            if (debug_subspace_matrix) write (6,*) 'i,j, m_ij:', i,j, mij_result_tens(1)

            if (present(frequency)) then 
                if (i==j) then 
                    mij_result_tens(1) =mij_result_tens(1) - frequency
                end if 
            end if

            workspace_t%M(i,j) = mij_result_tens(1) 
            mij_result_tens(1) = ZERO
         end do
      end do

      M_small = workspace_t%M(1:workspace_t%L,1:workspace_t%L)

      if (debug_subspace_matrix) then
         write (6,*) ' >>> debug, subspace Hamiltonian matrix M'
         do i = 1, workspace_t%L
           print "(2x,i6,$)", i
           do j = 1, workspace_t%L
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

   subroutine talsh_solve_linear_system_davidson(results_t,int_t,interm_t, &
 &            workspace_t, projector_t, frequency, prop_grad1, prop_grad2)
      type(talsh_davidson_results_tens_t) :: results_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(talsh_projectors_tens_t), intent(inout) :: projector_t
      type(talsh_intg_tens_t), intent(inout) :: int_t
      type(talsh_tens_t), intent(inout) :: prop_grad1, prop_grad2
      complex(kind=8), intent(in) :: frequency

!     do something to populate the property gradient, within workspace_t, before calling the routine below
      call talsh_get_values_vectors_davidson(results_t,int_t,interm_t, &
 &          workspace_t, projector_t, frequency, prop_grad1, prop_grad2)

   end subroutine



! aspg, 16/04/2019
! for the moment we have a long argument list to get it going, later we try to shorten in
!   subroutine talsh_get_values_vectors_davidson(A, workspace, projector)
! removed property gradient since that will be avaliable from within workspace_t
   subroutine talsh_get_values_vectors_davidson(results_t,int_t,interm_t, &
 &               workspace_t, projector_t, frequency, prop_grad1, prop_grad2, tprop_tensor)
      type(talsh_davidson_results_tens_t) :: results_t
      type(talsh_intermediates_tens_t),  intent(inout) :: interm_t
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
      type(talsh_projectors_tens_t), intent(inout) :: projector_t
      type(talsh_intg_tens_t), intent(inout) :: int_t

      complex(kind=8), intent(in), optional :: frequency  
      type(talsh_tens_t), intent(inout), optional :: prop_grad1, prop_grad2
      type(talsh_comm_tens_t), intent(inout), optional    :: tprop_tensor
      complex(kind=8), pointer, contiguous :: x_small(:,:)        => NULL()
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal

!       type(talsh_lambda_tens_t), intent(in)   :: A
!       type(talsh_davidson_workspace_tens_t),  intent(inout) :: workspace
!       type(talsh_projectors), intent(inout) :: projector
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

        complex(kind=8) :: frequency_local
        logical :: solve_linear = .false.
        logical :: solve_eigen  = .true.

        integer :: niter
        integer :: i, j, k
        integer :: nroots, ierr

        ! if (davidson_solve_linear_system) then
        !     nroots = 1
        ! else
            nroots = workspace_t%nroots
        ! end if

! aspg need to think of a consistency check for whether the property gradient has been set
!      to replace this
!       if (davidson_solve_linear_system.and..not.present(prop_gradient_t1)) then
!           call quit('bla')
!       end if

        if (present(frequency)) then
           frequency_local = frequency
        else
           frequency_local = (0.0d0,0.0d0) 
        end if
!
! davidson_create_results() will allocate storage for results from davidson diagonalization when called a first time. 
! 
!!      call davidson_create_results(results,vectors_N,nroots,&
!!                                   davidson_solve_values,   &
!!                                   davidson_solve_lhs,      &
!!                                   davidson_solve_rhs)

        workspace_t%L = nroots
        workspace_t%Lprev = 1

        write(davidson_output_file_unit,"(4x,a,i6,a)") "Requested ",nroots," eigenvales"
        if ( nroots .gt. workspace_t%Nmax ) then 
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested larger than matrix order !"
           return 
        else if ( nroots .eq. workspace_t%Nmax ) then
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested equal to matrix order !"
        endif

       if (davidson_overlap_sorting) &
           call setup_overlap_sorting(c_small_r_old, c_small_l_old, nroots, davidson_solve_rhs, davidson_solve_lhs)


        write(davidson_output_file_unit,*) ""
        do niter = 1, davidson_max_iterations
            write(davidson_output_file_unit,*) "Iteration ", niter
            write(davidson_output_file_unit,*) ""

!           call apply_projector(workspace_t%L, projector_t, workspace_t%trial1,  workspace_t%trial2)
            call create_sigma(workspace_t,int_t,interm_t,davidson_solve_rhs,davidson_solve_lhs)
!           call apply_projector(workspace_t%L, projector_t, workspace_t%sigma1,  workspace_t%sigma2)
            call create_subspace_matrix(M_small,workspace_t,davidson_solve_rhs,davidson_solve_lhs,frequency_local)

! aspg, added from (q)rsp
            if (davidson_solve_linear_system) then
               call create_subspace_property_gradient(C_small,workspace_t, &
 &                  prop_grad1,prop_grad2,davidson_solve_rhs,davidson_solve_lhs)
               call calculate_Linear_systems(davidson_solve_lhs, M_small,C_small,x_small)
! aspg, added from (q)rsp
            else
               call calculate_eigenvalues_and_eigenvectors(davidson_output_file_unit,    &
                                                           M_small,                      &
                                                           e_small,                      &
                                                           c_small_r, c_small_l,         &
                                                           c_small_r_old, c_small_l_old, &
                                                           davidson_eigenvalues_shift,   &
                                                           follow_overlap=(davidson_overlap_sorting.and.(niter.gt.1)),  &
                                                           follow_overlap_scheme=davidson_overlap_sorting_scheme, &
                                                           apply_evalues_shift=davidson_apply_eigenvalues_shift, &
                                                           left=davidson_solve_lhs,      &
                                                           right=davidson_solve_rhs,     &
                                                           nroots=nroots,                &
                                                           verbose=davidson_verbose)
!                                                          verbose=.true.)
            end if

            if (davidson_solve_rhs.and.davidson_overlap_sorting) then
                c_small_r_old = c_small_r
            end if
            if (davidson_solve_lhs.and.davidson_overlap_sorting) then
                c_small_l_old = c_small_l
            end if
 
! aspg, added from (q)rsp
            if (davidson_solve_linear_system) then
                allocate(e_small(1))
                e_small(1) = frequency_local
            end if
! aspg, added from (q)rsp
            call form_preconditioner(workspace_t, e_small)

            if (associated(delta)) deallocate(delta)
            ! if (davidson_solve_rhs) then
                if (davidson_solve_linear_system) then
!                     delta => correction_to_vectors(davidson_output_file_unit, &
!    &                                        e_small, x_small, workspace_t, C_small)
                    delta => correction_to_vectors(davidson_output_file_unit, &
   &                                        e_small, x_small, workspace_t, prop_grad1,prop_grad2)
                else
                    if (davidson_solve_rhs) then
                        delta => correction_to_vectors(davidson_output_file_unit, &
   &                                               e_small, c_small_r, workspace_t)
                    else 
                        delta => correction_to_vectors(davidson_output_file_unit, &
                        &                          e_small, c_small_l, workspace_t)
                    end if
               end if
            ! end if
            if (.not.associated(delta)) then
                if (davidson_solve_linear_system) then 
                    write(davidson_output_file_unit,*) "Root converged."
                    call get_tensor_linear_systems(tprop_tensor,x_small,workspace_t)
                else
                    write(davidson_output_file_unit,*) "All roots converged, stopping."
                end if
                exit
            end if

! not converged yet, if asked to refresh the trial vectors do it now, otherwise create
! new trial vectors
            if (workspace_t%L .ge. workspace_t%Lmax) then
               write(davidson_output_file_unit, *) "Reached maximum subspace size, stopping."
               exit
            end if

            if (davidson_refresh_trial_vectors.and.(mod(niter,davidson_trial_vector_refresh_rate).eq.0)) then
                print *,'debug : refreshing trial vectors not implemented in this version'
            else
               call update_trial_vectors(delta,workspace_t)
            end if
 
        end do ! davidson iterations

        if (davidson_solve_rhs.and.davidson_overlap_sorting) deallocate(c_small_r_old)
        if (davidson_solve_lhs.and.davidson_overlap_sorting) deallocate(c_small_l_old)


        write(davidson_output_file_unit, *) "Final size of reduced subspace : ",workspace_t%L
        write(davidson_output_file_unit, *), "*****************EOM-CC Eigenvalules******************"

        if (davidson_solve_linear_system) then 
            write(davidson_output_file_unit, *) "No Eigenvales since Solving linear systems"
        else 
            call print_values_vectors(davidson_output_file_unit,e_small, c_small_r, c_small_l, &
                .false.,.false., nroots)
        end if 
        write(davidson_output_file_unit, *), "******************************************************"

        
!!      call form_eigenvalues_full_problem(results%eValues, nroots, e_small) 

!!      if (davidson_solve_lhs) then 
!!          call form_eigenvector_full_problem(results%eVectorsL, nroots,b,c_small_l,davidson_solve_rhs,davidson_solve_lhs)
!!      end if
!!      if (davidson_solve_rhs) then 
!!          call form_eigenvector_full_problem(results%eVectorsR, nroots,b,c_small_r,davidson_solve_rhs,davidson_solve_lhs)
!!      end if

!!      if (associated(b)) deallocate(b)
!!      deallocate(M_small)
!!      if (associated(sigma_vector)) deallocate(sigma_vector)
!!      deallocate(e_small)
!!      if (associated(c_small_r)) deallocate(c_small_r)
!!      if (associated(c_small_l)) deallocate(c_small_l)

    end subroutine




    ! subroutine rotate_eigenvector(sigma1, sigma2, workspace_t)
    !     type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace_t
    !     type(talsh_tens_t) sigma1, sigma2


    !     type(C_PTR):: body_p
    !     integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
    !     integer(C_INT):: host_arg_max

    !     integer(C_INT)              :: ierr
    !     integer(INTD), dimension(1) :: one_dims

    !     integer :: i,j,a 
    !     integer :: N_trial1, N_trail2

    !     type(talsh_tens_t)   :: result_tensor
    !     complex(8), pointer  :: result_tens(:), sigma1_tens(:), sigma2_tens(:)


    !     one_dims(1) = 1
    !     N_trail1 = size(workspace_t%trial1,1)
    !     N_trail2 = size(workspace_t%trial2,1)


    !     ierr=talsh_tensor_construct(result_tensor,C8,one_dims,init_val=ZERO)

    !     ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    !     call c_f_pointer(body_p,result_tens,one_dims)


    !     do i=1, N_trail1
    !         ierr=talsh_tensor_contract("V1()+=B(p,q)*S(p,q)",result_tensor,workspace_t%trial1(i),sigma1) ! 
    !         sigma1_tens(i) = result_tens(1)
    !     end do 


    !     do i=1, N_trail2
    !         ierr=talsh_tensor_contract("V2()+=B(a,b,c)*S(a,b,c)",result_tensor,workspace_t%trial2(i),sigma2) ! 
    !     end do 



    !     ierr=talsh_tensor_contract("V1(a,b)+=B(a,m)*V(m,b)",eigenvec_tensor_new1,workspace_t%trial1,eigenvec_tensor_old) ! 
    !     ierr=talsh_tensor_contract("V2(a,b)+=B(a,m)*V(m,b)",eigenvec_tensor_new2,workspace_t%trial2,eigenvec_tensor_old) ! 


    !     ierr=talsh_tensor_get_body_access(eigenvec_tensor_new1,body_p,C8,int(0,C_INT),DEV_HOST)
    !     call c_f_pointer(body_p,eigenvec_tens1,eigenvec_dims_s)

    !     ierr=talsh_tensor_get_body_access(eigenvec_tensor_new2,body_p,C8,int(0,C_INT),DEV_HOST)
    !     call c_f_pointer(body_p,eigenvec_tens2,eigenvec_dims_d)




    !     ierr=talsh_tensor_destruct(eigenvec_tensor_new1)
    !     ierr=talsh_tensor_destruct(eigenvec_tensor_new2)

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
        call talsh_solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left_local)
 

        ! print *, "eVectorsR(1) before zgemm", eVectorsR(:,1)
        ! print *, "eVectorsR(2) before zgemm", eVectorsR(:,2)
        ! print *, "eVectorsR(3) before zgemm", eVectorsR(:,3)
        ! print *, "eVectorsR(4) before zgemm", eVectorsR(:,4)
        ! print *, "eVectorsR(5) before zgemm", eVectorsR(:,5)
        ! print *, "eVectorsR(6) before zgemm", eVectorsR(:,6)
        ! print *, "eVectorsR(7) before zgemm", eVectorsR(:,7)
        ! print *, "eVectorsR(8) before zgemm", eVectorsR(:,8)
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
            call talsh_sort_values_vectors(sorting_order_local, eValues, eVectorsR, eVectorsL, &
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
        call print_values_vectors(davidson_output_file_unit,eValues, eVectorsR, eVectorsL, &
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

        call talsh_solve_Linear_systems(left, M, B, Solution)
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

        vectors_N = size(b,1)
        vectors_L = size(b,2)
        vectors_L_old = size(old_sigma_vector,2)

        allocate(update_sigma_vector(vectors_N,vectors_L))
        do i = 1, vectors_L_old 
           update_sigma_vector(:, i) = old_sigma_vector(:, i)  
        end do
        deallocate(old_sigma_vector)

        new_L = vectors_L - vectors_L_old
        if (davidson_complex_mode_with_reals) then
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

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_sigma_vector_full_matrix(vectors_N,vectors_L))

        if (davidson_complex_mode_with_reals) then
            vectors_N = vectors_N / 2
            call zgemm('N', 'N', vectors_N, vectors_L, vectors_N, (1.0d0, 0.0d0), M, vectors_N, &
                       b, vectors_N, (0.0d0, 0.0d0), create_sigma_vector_full_matrix, vectors_N)
        end if
    end function

    subroutine talsh_sort_values_vectors(order, e_small, c_small_r, c_small_l, right, left, apply_shift, energy_shift) 
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
        temp_vector = cmplx(0.0d0, 0.0d0)

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
       type(talsh_projectors_tens_t), intent(inout) :: projector_t
       type(talsh_tens_t), intent(inout) :: target1(:), target2(:)

       integer :: p, q, r, s, ierr, i
       complex(8), pointer, contiguous:: target1_tens(:,:), target2_tens(:,:,:,:)
       type(C_PTR):: target1_body_p, target2_body_p
       complex(8), pointer, contiguous:: proj1_tens(:,:), proj2_tens(:,:,:,:)
       type(C_PTR):: proj1_body_p, proj2_body_p

       logical :: debug = .true.

       ierr = talsh_tensor_get_body_access(projector_t%p1,proj1_body_p,C8,0,DEV_HOST)
       ierr = talsh_tensor_get_body_access(projector_t%p2,proj2_body_p,C8,0,DEV_HOST)
       call c_f_pointer(proj1_body_p,proj1_tens,projector_t%p1_dims(1:projector_t%p1_rank))
       call c_f_pointer(proj2_body_p,proj2_tens,projector_t%p2_dims(1:projector_t%p2_rank))

       do i = 1, nvectors

          ierr = talsh_tensor_get_body_access(target1(i),target1_body_p,C8,0,DEV_HOST)
          ierr = talsh_tensor_get_body_access(target2(i),target2_body_p,C8,0,DEV_HOST)

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
!         ierr=talsh_tensor_contract("T(p,q)+=T(p,q)*P(p,q)",target1(i),target1(i),p1)
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
!         ierr=talsh_tensor_contract("T(p,q,r,s)+=T(p,q,r,s)*P(p,q,r,s)",target2(i),target2(i),p2)
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

    subroutine talsh_davidson_get_eigenvalue_shift(apply_shift,shift)
       logical, intent(inout)      :: apply_shift
       real(kind=8), intent(inout) :: shift(2)

       shift = davidson_eigenvalues_shift
       apply_shift = davidson_apply_eigenvalues_shift 

    end subroutine

    subroutine divide_tensor_elements_by_value(a_tensor, denominator)
       type(talsh_tens_t), intent(inout)    :: a_tensor
       complex(8), intent(in) :: denominator

       integer(INTD)  ::  a_rank
       integer(INTD)  ::  a_dims1(1)
       integer(INTD)  ::  a_dims2(1:2)
       integer(INTD)  ::  a_dims3(1:3)
       integer(INTD)  ::  a_dims4(1:4)
       integer(INTD)  ::  ierr
       integer  ::  i, j, k, l, m, n, o, p
       type(C_PTR)         :: a_body
       complex(8), pointer :: a_tens1(:)
       complex(8), pointer :: a_tens2(:,:)
       complex(8), pointer :: a_tens3(:,:,:)
       complex(8), pointer :: a_tens4(:,:,:,:)

       a_rank=talsh_tensor_rank(a_tensor)

       if (a_rank.eq.1) then 
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims1)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens1,a_dims1)

         do i = 1, a_dims2(1)
            a_tens1(i)=a_tens1(i) / denominator 
         end do

       else if (a_rank.eq.2) then
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims2)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens2,a_dims2)

         do i = 1, a_dims2(1)
           do j = 1, a_dims2(2)
             a_tens2(i,j)=a_tens2(i,j) / denominator
           end do
         end do

       else if (a_rank.eq.3) then 
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims3)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens3,a_dims3)

         do i = 1, a_dims3(1)
            do j = 1, a_dims3(2)
              do k = 1, a_dims3(3)
                 a_tens3(i,j,k)=a_tens3(i,j,k) / denominator
             end do
           end do
         end do

       else if (a_rank.eq.4) then 
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims4)
         if (ierr.ne.0) call quit('error in divide_tensor_elements_by_value')

         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
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

       else
          call quit('unsupported tensor rank in divide_tensor_elements_by_value')
       end if



    end subroutine

    subroutine multiply_tensor_elements(a_tensor, b_tensor)
! this routine performs element-wise multiplication of two tensors
! and stores the outcome in the first argument
!
! a(i,j,k,...) = a(i,j,k,..)*b(i,j,k,...)
!
       implicit none
    
       type(talsh_tens_t), intent(inout)    :: a_tensor,b_tensor
    
       integer(INTD)  ::  a_rank, b_rank
       integer(INTD)  ::  a_dims1(1),b_dims1(1)
       integer(INTD)  ::  a_dims2(1:2),b_dims2(1:2)
       integer(INTD)  ::  a_dims3(1:3),b_dims3(1:3)
       integer(INTD)  ::  a_dims4(1:4),b_dims4(1:4)
       integer(INTD)  ::  ierr
       integer  ::  i, j, k, l
       type(C_PTR)         :: a_body,b_body
       complex(8), pointer :: a_tens1(:)
       complex(8), pointer :: b_tens1(:)
       complex(8), pointer :: a_tens2(:,:)
       complex(8), pointer :: b_tens2(:,:)
       complex(8), pointer :: a_tens3(:,:,:)
       complex(8), pointer :: b_tens3(:,:,:)
       complex(8), pointer :: a_tens4(:,:,:,:)
       complex(8), pointer :: b_tens4(:,:,:,:)
       complex(8) :: mult_val 
    
       a_rank=talsh_tensor_rank(a_tensor)
       b_rank=talsh_tensor_rank(b_tensor)
    
       if ((a_rank.eq.1).and.(b_rank.eq.1)) then
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims1)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
         ierr = talsh_tensor_dimensions(b_tensor,b_rank,b_dims1)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
    
         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens1,a_dims1)
         ierr=talsh_tensor_get_body_access(b_tensor,b_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(b_body,b_tens1,b_dims1)
    
! add test to check dimensions
!        a_tens1(:) = a_tens1(:)*b_tens1(:)
         do i = 1, a_dims1(1)
            a_tens1(i)=a_tens1(i)*b_tens1(i)
         end do
    
       else if ((a_rank.eq.2).and.(b_rank.eq.2)) then
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims2)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
         ierr = talsh_tensor_dimensions(b_tensor,b_rank,b_dims2)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
    
         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens2,a_dims2)
         ierr=talsh_tensor_get_body_access(b_tensor,b_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(b_body,b_tens2,b_dims2)
    
! add test to check dimensions
!        a_tens2(:,:) = a_tens2(:,:)*b_tens2(:,:)
         a_tens2 = a_tens2*b_tens2
!        do i = 1, a_dims2(1)
!          do j = 1, a_dims2(2)
!            a_tens2(i,j)=a_tens2(i,j)*b_tens2(i,j)
!          end do
!        end do
    
       else if ((a_rank.eq.3).and.(b_rank.eq.3)) then
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims3)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
         ierr = talsh_tensor_dimensions(b_tensor,b_rank,b_dims3)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
    
         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens3,a_dims3)
         ierr=talsh_tensor_get_body_access(b_tensor,b_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(b_body,b_tens3,b_dims3)
    
! add test to check dimensions
!        a_tens3(:,:,:) = a_tens3(:,:,:)*b_tens3(:,:,:)
         do i = 1, a_dims3(1)
            do j = 1, a_dims3(2)
              do k = 1, a_dims3(3)
                 a_tens3(i,j,k)=a_tens3(i,j,k)*b_tens3(i,j,k)
             end do
           end do
         end do
    
       else if ((a_rank.eq.4).and.(b_rank.eq.4)) then
         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims4)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
         ierr = talsh_tensor_dimensions(b_tensor,b_rank,b_dims4)
         if (ierr.ne.0) call quit('error in multiply_tensor_elements')
       
         ierr=talsh_tensor_get_body_access(a_tensor,a_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(a_body,a_tens4,a_dims4)
         ierr=talsh_tensor_get_body_access(b_tensor,b_body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(b_body,b_tens4,b_dims4)
    
! add test to check dimensions
!        a_tens4(:,:,:,:) = a_tens4(:,:,:,:)*b_tens4(:,:,:,:)
         a_tens4 =  a_tens4 * b_tens4
!        do i = 1, a_dims4(1)
!           do j = 1, a_dims4(2)
!             do k = 1, a_dims4(3)
!               do l = 1, a_dims4(4)
!                  mult_val = a_tens4(i,j,k,l)*b_tens4(i,j,k,l)
!                  a_tens4(i,j,k,l)=a_tens4(i,j,k,l)*b_tens4(i,j,k,l)
!               end do
!            end do
!          end do
!        end do
    
       else
          call quit('mismatch in tensor ranks in multiply_tensor')
       end if
    
    end subroutine 

    subroutine print_tensor_elements(tensor)
       implicit none
    
       type(talsh_tens_t), intent(inout)    :: tensor
    
       integer(INTD)  ::  rank
       integer(INTD)  ::  dims1(1)
       integer(INTD)  ::  dims2(1:2)
       integer(INTD)  ::  dims3(1:3)
       integer(INTD)  ::  dims4(1:4)
       integer(INTD)  ::  ierr
       integer  ::  i, j, k, l, b
       type(C_PTR)         :: body
       complex(8), pointer :: tens1(:)
       complex(8), pointer :: tens2(:,:)
       complex(8), pointer :: tens3(:,:,:)
       complex(8), pointer :: tens4(:,:,:,:)
       real(8)    :: threshold = 1.0e-10
       complex(8) :: value_to_print
    
       rank=talsh_tensor_rank(tensor)
    
       if (rank.eq.1) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims1)
         ierr=talsh_tensor_get_body_access(tensor,body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body,tens1,dims1)
    
         b = 0
         do i = 1, dims1(1)
            b = b + 1
            value_to_print = tens1(i)
            if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = ZERO
            else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = cmplx(real(value_to_print),0.0d0) 
            else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = cmplx(0.0d0,aimag(value_to_print)) 
            end if
            write(6,'(4x,i6,1e14.6,2x,SP,1e14.6,"i")') b, value_to_print
         end do
    
       else if (rank.eq.2) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims2)
         ierr=talsh_tensor_get_body_access(tensor,body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body,tens2,dims2)
    
         b = 0
         do j = 1, dims2(2)
           do i = 1, dims2(1)
            b = b + 1
            value_to_print = tens2(i,j)
            if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = ZERO
            else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = cmplx(real(value_to_print),0.0d0)
            else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
               value_to_print = cmplx(0.0d0,aimag(value_to_print))
            end if
            write(6,'(4x,i6,2x,i6,2x,i6,2x,1e14.6,2x,SP,1e14.6,"i")') b, i, j, value_to_print

           end do
         end do
    
       else if (rank.eq.3) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims3)
         ierr=talsh_tensor_get_body_access(tensor,body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body,tens3,dims3)

         do k = 1, dims3(3)
            do j = 1, dims3(2)
              do i = 1, dims3(1)
                 b = b + 1
                 value_to_print = tens3(i,j,k)
                 if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                    value_to_print = ZERO
                 else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
                    value_to_print = cmplx(real(value_to_print),0.0d0)
                 else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                    value_to_print = cmplx(0.0d0,aimag(value_to_print))
                 end if
                 write(6,'(4x,i6,2x,i6,2x,i6,2x,i6,2x,1e14.6,2x,SP,1e14.6,"i")') b, i, j, k, value_to_print

             end do
           end do
         end do
    
       else if (rank.eq.4) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims4)
         ierr=talsh_tensor_get_body_access(tensor,body,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body,tens4,dims4)

         b = 0
         do l = 1, dims4(4)
            do k = 1, dims4(3)
              do j = 1, dims4(2)
                do i = 1, dims4(1)
                   b = b + 1
                   value_to_print = tens4(i,j,k,l)
                   if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                      value_to_print = ZERO
                   else if ((real(value_to_print).gt.threshold).and.(aimag(value_to_print).le.threshold)) then
                      value_to_print = cmplx(real(value_to_print),0.0d0)
                   else if ((real(value_to_print).le.threshold).and.(aimag(value_to_print).le.threshold)) then
                      value_to_print = cmplx(0.0d0,aimag(value_to_print))
                   end if
                   write(6,'(4x,i6,2x,i6,2x,i6,2x,i6,2x,i6,2x,1e14.6,2x,SP,1e14.6,"i")') b, i, j, k, l, value_to_print

                end do
             end do
           end do
         end do
    
       else
          call quit('unsupported tensor rank in print_tensor_elements')
       end if
    
    end subroutine 

end module

