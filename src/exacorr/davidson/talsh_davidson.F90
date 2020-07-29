module talsh_davidson

use talsh
use tensor_algebra
use talsh_common_routines
use exacorr_datatypes
use talsh_ao_to_mo, only : tensor_norm2


implicit none
complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                         ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                         MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                         MINUS_TWO=(-2.D0,0.D0)


public talsh_davidson_get_eigenvalue_shift
public talsh_davidson_setup
public talsh_davidson_destroy
public talsh_davidson_driver
!public davidson_print_vectors
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

public talsh_davidson_results

private

! derived types

type talsh_davidson_results
    real(kind=8), pointer :: eValues(:,:) 
    real(kind=8), pointer :: eVectorsR(:,:)
    real(kind=8), pointer :: eVectorsL(:,:) 
end type


! module parameters

integer, parameter :: max_label_length = 30
integer, parameter :: nr_max_irreps = 32
real(kind=8), parameter :: davidson_eigenvalues_shift_threshold = 1.0e-3

! module variables, should be set with a call to davidson_setup

real(kind=8), save :: safe_minimum = 1.0d-36

real(kind=8), save :: davidson_convergence_threshold 
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
        type(talsh_davidson_results), intent(inout) :: results
        integer, intent(in)    :: N, nroots 
        integer                :: N_local
        logical, intent(in)    :: values, left, right

! we don't use complex variables but rather use twice the amount of real variables
! for symmetric matrices we only need one of the vectors, and we choose it to be the right one
! left vectors, if asked, are ignored
        N_local = N

        results%eValues => NULL()
        results%eVectorsR => NULL()
        results%eVectorsL => NULL()


        if (davidson_symmetric) then
            if (values) then 
                allocate(results%eValues(1,nroots))
                results%eValues   = 0.0d0
            end if
            if (right) then
                allocate(results%eVectorsR(N_local,nroots))
                results%eVectorsR = 0.0d0
            end if
        else
            if (values) then
                allocate(results%eValues(2,nroots))
                results%eValues   = 0.0d0
            end if
            if (right) then
                allocate(results%eVectorsR(N_local,nroots))
                results%eVectorsR = 0.0d0
            else
                allocate(results%eVectorsR(1,nroots))
                results%eVectorsR = 0.0d0
            end if
            if (left) then
                allocate(results%eVectorsL(N_local,nroots))
                results%eVectorsL = 0.0d0
            else
                allocate(results%eVectorsL(1,nroots))
                results%eVectorsL = 0.0d0
            end if
        end if
    end subroutine


    subroutine free_davidson_results(results)
        type(talsh_davidson_results), intent(inout) :: results

        if (associated(results%eValues)) deallocate(results%eValues)
        if (associated(results%eVectorsR)) deallocate(results%eVectorsR)
        if (associated(results%eVectorsL)) deallocate(results%eVectorsL)
    end subroutine


    subroutine talsh_davidson_setup(convergence_threshold, &
                              max_subspace_size, &
                              max_iterations,    & 
                              unit_output,       & 
                              unit_sigmas,       & 
                              unit_subspaces,    &
                              unit_results,      &
                              solve_values,      &
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
            davidson_max_subspace_size = 8 ! 1024
        end if

        if (present(max_iterations)) then
            davidson_max_iterations = max_iterations
        else
            davidson_max_iterations = 50
        end if

        if (present(solve_values)) then
            davidson_solve_values = solve_values
        else
            davidson_solve_values = .true. 
        end if

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
        type(talsh_davidson_results), intent(inout) :: results
        call free_davidson_results(results)
    end subroutine


    subroutine update_trial_vectors(delta,workspace)
        type(talsh_davidson_workspace) :: workspace
        logical, pointer, intent(inout) :: delta(:)

        integer :: ierr, i
        integer(INTD) :: rank

        rank=talsh_tensor_rank(workspace%trial1(1))
        if (rank.gt.2) call quit('rank not implemented in update_trial_vectors')

        workspace%Lprev = workspace%L
        do i = 1, size(delta,1)
           if (delta(i)) then
              workspace%L = workspace%L + 1

              ierr=talsh_tensor_init(workspace%trial1(workspace%L), ZERO)
              ierr=talsh_tensor_init(workspace%trial2(workspace%L), ZERO)
              if (rank.eq.1) then
                 ierr=talsh_tensor_add("t(a)+=d(a)",workspace%trial1(workspace%L),workspace%delta1(i))
                 ierr=talsh_tensor_add("t(a,b,c)+=d(a,b,c)",workspace%trial2(workspace%L),workspace%delta2(i))
              else if (rank.eq.2) then
                 ierr=talsh_tensor_add("t(a,b)+=d(a,b)",workspace%trial1(workspace%L),workspace%delta1(i))
                 ierr=talsh_tensor_add("t(a,b,c,d)+=d(a,b,c,d)",workspace%trial2(workspace%L),workspace%delta2(i))
              end if
           end if 
        end do
        call orthonormalize_trial_vectors(workspace%Lprev, workspace%L, rank, workspace%trial1, workspace%trial2)

    end subroutine

    subroutine orthonormalize_trial_vectors(L,Lnew,rank1,trial1,trial2)
        integer, intent(in) :: L, Lnew, rank1
        type(talsh_tens_t), pointer, intent(inout) :: trial1(:), trial2(:) 

        call gram_schmidt_orthogonalization_tensors(L,Lnew,rank1,trial1,trial2)
    end subroutine


    subroutine gram_schmidt_projection_R(projection, u, v)
        real(kind=8), intent(inout) :: projection(:)
        real(kind=8), intent(in) :: u(:)
        real(kind=8), intent(in) :: v(:)
        real(kind=8) :: numerator, denominator, factor
        real(kind=8) :: norm_denominator
        integer :: N, N_complex
        double precision ddot, dznrm2
        real(kind=8) :: numerator_complex(2), denominator_complex(2), factor_complex(2)

        N = size(u,1)

        projection = 0.0d0 

        if (davidson_complex_mode_with_reals) then
            N_complex = N / 2

            call xdotc(denominator_complex,N_complex, u, 1, u, 1) 
            call xdotc(numerator_complex,N_complex, u, 1, v, 1) 

            norm_denominator = dznrm2(N_complex, u, 1)

            if (norm_denominator .gt. safe_minimum) then
                call dladiv(numerator_complex(1),numerator_complex(2), &
                            denominator_complex(1),denominator_complex(2), &
                            factor_complex(1), factor_complex(2))
            else
                factor_complex(1) = 0.0d0
                factor_complex(2) = 0.0d0
            end if
            projection = u
            call zscal(N_complex, factor_complex, projection, 1)
        end if

    end subroutine 

    subroutine gram_schmidt_orthogonalization(vectors, start_from, verbose)
        real(kind=8), intent(inout) :: vectors(:,:)
        integer, intent(in) :: start_from
        logical, optional, intent(in) :: verbose

        real(kind=8), allocatable   :: projection(:)
        integer :: vectors_L, vectors_N
        integer :: i, j, k
        logical :: verbose_local
        real(kind=8) :: norm, vr, vi, pr, pi

        vectors_L = size(vectors,2)
        vectors_N = size(vectors,1)

        allocate(projection(vectors_N))

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        do j = start_from, vectors_L
            norm = vector_norm(vectors(:,j))
            if (norm .ge. safe_minimum ) then
                vectors(:,j) = vectors(:,j) / norm
            else
                stop ('division by zero in orthogonalization') 
            end if
            do i = 1, j-1
                call gram_schmidt_projection_R(projection, vectors(:,i),vectors(:,j))
                vectors(:,j) = vectors(:,j) - projection(:)
            end do
            norm = vector_norm(vectors(:,j))
            if (norm .ge. safe_minimum ) then
                vectors(:,j) = vectors(:,j) / norm
            else
                stop ('division by zero in orthogonalization')
            end if
        end do

        deallocate(projection)
    end subroutine 


    subroutine renormalize_vectors(v1,v2)
       type(talsh_tens_t), intent(inout) :: v1, v2
       real :: norm_r
       complex(8) :: norm

       norm_r = tensor_norm2(v1) +  tensor_norm2(v2)
       norm = CMPLX(norm_r, 0.0d0)
       norm_r = 1.0d0 / norm_r

       if (norm_r .ge. safe_minimum ) then
          call divide_tensor_elements_by_value(v1, norm)
          call divide_tensor_elements_by_value(v2, norm)
       else
          call quit ('norm of tensors too small for division, quitting')
       end if

    end subroutine

    subroutine gram_schmidt_orthogonalization_tensors(L,Lnew,rank1,trial1,trial2)
        integer, intent(in) :: L, Lnew, rank1
        type(talsh_tens_t), intent(inout) :: trial1(:), trial2(:)
        real :: norm_r
        integer(INTD) :: r1_rank
        integer :: ierr, i, j
        type(talsh_tens_t) :: projection 

        if (rank1.gt.2) call quit('rank not implemented in form_r_tensors')

        do j = L, Lnew
            call renormalize_vectors(trial1(j),trial2(j))
!           do i = 1, j-1
!               call gram_schmidt_projection_tensors(projection1, projectionvectors(:,i),vectors(:,j))

!               if (rank1.eq.1) then
!                  ierr=projection1
!                  ierr=projection2

!                  ierr=talsh_tensor_add("t(a)+=p(a)",trial1(j),projection1,scale=MINUS_ONE)
!                  ierr=talsh_tensor_add("t(a,b,c)+=p(a,b,c)",trial2(j),projection2,scale=MINUS_ONE)
!               else if (rank1.eq.2) then
!!add me 
!               end if
!           end do
!           call renormalize_vectors(trial1(j),trial2(j))
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
            vr = realpart(eValues(i))
            vi = imagpart(eValues(i))
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

    subroutine form_r_tensors(k,eigenvectors_subspace,eigenvalues_subspace,workspace)
! we use delta{1,2,...} to store the values of r^k
        integer, intent(in) :: k
        complex(kind=8), intent(in)  :: eigenvectors_subspace(:,:)
        complex(kind=8), intent(in)  :: eigenvalues_subspace(:)
        type(talsh_davidson_workspace), intent(inout) :: workspace
        complex(kind=8) :: factor1, factor2
        integer(INTD) :: r1_rank 
        integer :: ierr, i
        logical :: debug_form_r = .true.
        type(talsh_tens_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
! we only need to know the rank of delta1 because the rank of r2 will always be r1_rank+2, that is: 
! r1 = 1 (1h for IP, 1p for EA) then r2 = 3 (2h1p for IP, 1h2p for EA)
! r1 = 2 (1h1p for EE, 2h for DIP, 2p for DEA) then r2 = 4 (2h2p for EE, 3h1p for DIP, 1h3p for DEA)
! etc
        r1_rank=talsh_tensor_rank(workspace%delta1(k))
        if (r1_rank.gt.2) call quit('rank not implemented in form_r_tensors')

!       ierr=talsh_tensor_init(workspace%delta1(k), ZERO)
!       ierr=talsh_tensor_init(workspace%delta2(k), ZERO)

!       if (debug_form_r) then
!          call davidson_print_vectors(eigenvectors_subspace,   label="eigenvectors_subspace in form_r")
!       end if

        do i = 1, workspace%L

           factor1 =  eigenvectors_subspace(i,k)
           factor2 = -eigenvalues_subspace(k)*factor1
           if (debug_form_r) write(6,*) 'c(',i,',',k,')',factor1,'lambda*c',factor2

           if (r1_rank.eq.1) then
              ierr=talsh_tensor_contract("r(a)+=s(a)*K()",workspace%delta1(k),workspace%sigma1(i),one_tensor,scale=factor1)
              ierr=talsh_tensor_contract("r(a)+=t(a)*K()",workspace%delta1(k),workspace%trial1(i),one_tensor,scale=factor2)

              ierr=talsh_tensor_contract("r(a,b,c)+=s(a,b,c)*K()",workspace%delta2(k),workspace%sigma2(i),one_tensor,scale=factor1)
              ierr=talsh_tensor_contract("r(a,b,c)+=t(a,b,c)*K()",workspace%delta2(k),workspace%trial2(i),one_tensor,scale=factor2)

            else if (r1_rank.eq.2) then
              ierr=talsh_tensor_add("r(a,b)+=s(a,b)",workspace%delta1(k),workspace%sigma1(i),scale=factor1)
              ierr=talsh_tensor_add("r(a,b)+=t(a,b)",workspace%delta1(k),workspace%trial1(i),scale=factor2)

              ierr=talsh_tensor_add("r(a,b,c,d)+=s(a,b,c,d)",workspace%delta2(k),workspace%sigma2(i),scale=factor1)
              ierr=talsh_tensor_add("r(a,b,c,d)+=t(a,b,c,d)",workspace%delta2(k),workspace%trial2(i),scale=factor2)

           end if

!          if (debug_form_r) then
!             write (6,*) 'debug >>> in form_r_tensors, accumulated r1 root',k,' for L',i,' >>>'
!             call print_tensor_elements(workspace%delta1(k))
!             write (6,*) 'debug <<< in form_r_tensors, accumulated r1 for L',i,' <<<'
!             write (6,*) 'debug >>> in form_r_tensors, accumulated r2 root',k,' for L',i,' >>>'
!             call print_tensor_elements(workspace%delta2(k))
!             write (6,*) 'debug <<< in form_r_tensors, accumulated r2 for L',i,' <<<'
!          end if
       
       end do

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
          call print_tensor_elements(workspace%delta1(k))
          write (6,*) 'debug <<< in form_r_tensors, r1 <<<'
          write (6,*) 'debug >>> in form_r_tensors, r2 root',k,' >>>'
          call print_tensor_elements(workspace%delta2(k))
          write (6,*) 'debug <<< in form_r_tensors, r1 <<<'
       end if

       ierr=talsh_tensor_destruct(one_tensor)

    end subroutine    



    function correction_to_vectors(output_file_unit, nroots, e_small_M, c_small_M, workspace, correct_all)
        integer, intent(in)       :: nroots, output_file_unit
        type(talsh_davidson_workspace), intent(inout) :: workspace
        complex(kind=8), intent(in)  :: e_small_M(:)
        complex(kind=8), intent(in)  :: c_small_M(:,:)
        logical, optional, intent(in) :: correct_all

        logical, pointer :: correction_to_vectors(:)

        integer :: ierr, k, unconverged_roots
        real(kind=8) :: norm_r = 0.0d0
        logical :: correct_all_local
        logical :: debug_correction_vectors = .true.

        if (present(correct_all)) then
            correct_all_local = correct_all
        else
            correct_all_local = .false. 
        end if

        allocate(correction_to_vectors(nroots))
        correction_to_vectors = .true.

        unconverged_roots = 0
        do k = 1, nroots
           ierr=talsh_tensor_init(workspace%delta1(k), ZERO)
           print *,'debug, init workspace%delta1(k) to zero, ierr',ierr
           ierr=talsh_tensor_init(workspace%delta2(k), ZERO)
           print *,'debug, init workspace%delta2(k) to zero, ierr',ierr
        end do
        do k = 1, nroots
! instead of having a new variable to store r(:,k), we will do collect the results on workspace%delta{1,2,...} 
           call form_r_tensors(k,c_small_M,e_small_M,workspace)
            
! get the norm of the whole intermediate r(:,k) tensor, which means to calculate it over the different excitation ranks 
           print *,'debug, norm of delta1(',k,'):',tensor_norm2(workspace%delta1(k))
           print *,'debug, norm of delta2(',k,'):',tensor_norm2(workspace%delta2(k))

           norm_r = tensor_norm2(workspace%delta1(k)) +  tensor_norm2(workspace%delta2(k))

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
! note here that, as seen above, workspace%delta{1,2,..} have been used to accumulate r(i,k) 
           do k = 1, nroots
              if (correction_to_vectors(k)) then
                 call multiply_tensor_elements(workspace%delta1(k), workspace%precond1(k))
                 call multiply_tensor_elements(workspace%delta2(k), workspace%precond2(k))
              else
                 ierr=talsh_tensor_init(workspace%delta1(k), ZERO)
                 ierr=talsh_tensor_init(workspace%delta2(k), ZERO)
              end if

              if (debug_correction_vectors) then
                 write (6,*) 'debug >>> in correction_to_vectors, delta1 root',k,' >>>'
                 call print_tensor_elements(workspace%delta1(k))
                 write (6,*) 'debug <<< in correction_to_vectors, delta1 <<<'
                 write (6,*) 'debug >>> in correction_to_vectors, delta2 root',k,' >>>'
                 call print_tensor_elements(workspace%delta2(k))
                 write (6,*) 'debug <<< in correction_to_vectors, delta2 <<<'
              end if

           end do
        else
            deallocate(correction_to_vectors)
            correction_to_vectors => NULL()
        end if

    end function

    subroutine form_preconditioner(workspace, eigenvalues)
        type(talsh_davidson_workspace) :: workspace
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
        logical :: preconditioner_debug = .true.
        real(kind=8) :: small_number = 1.0d-10

! careful with dimension 1 of eigenvalues, 
        nroots    = size(eigenvalues,1)

! whether or not we want to precondition, we start by setting up the tensor to one
        do iroot = 1, nroots
           ierr=talsh_tensor_init(workspace%precond1(iroot), ONE)   
           ierr=talsh_tensor_init(workspace%precond2(iroot), ONE)
        end do

        if (davidson_use_preconditioner) then
! we want to precondition
! first prepare to access the rank-1 and rank-2 diagonal of Hbar 
           ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
           ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

           if (workspace%tens1_rank .eq. 1) then
              call c_f_pointer(hdiag1_body_p,hdiag_1d_p,workspace%tens1_dims(1:workspace%tens1_rank))
              call c_f_pointer(hdiag2_body_p,hdiag_3d_p,workspace%tens2_dims(1:workspace%tens2_rank))
           else if (workspace%tens1_rank .eq. 2) then
              call c_f_pointer(hdiag1_body_p,hdiag_2d_p,workspace%tens1_dims(1:workspace%tens1_rank))
              call c_f_pointer(hdiag2_body_p,hdiag_4d_p,workspace%tens2_dims(1:workspace%tens2_rank))
           end if

           write (6,*) 'debug eigenvalues:',eigenvalues
! now we get go root by root, and construct the preconditioner for each
           do iroot = 1, nroots

              eigenvalue_r = real(eigenvalues(iroot))
              eigenvalue_i = aimag(eigenvalues(iroot))

              ierr = talsh_tensor_get_body_access(workspace%precond1(iroot),preco1_body_p,C8,0,DEV_HOST)
              ierr = talsh_tensor_get_body_access(workspace%precond2(iroot),preco2_body_p,C8,0,DEV_HOST)

              if (workspace%tens1_rank .eq. 1) then
                 call c_f_pointer(preco1_body_p,preco_1d_p,workspace%tens1_dims(1:workspace%tens1_rank))
                 call c_f_pointer(preco2_body_p,preco_3d_p,workspace%tens2_dims(1:workspace%tens2_rank))

                 do ja = 1, workspace%tens1_dims(1)
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

                 do jc = 1, workspace%tens2_dims(3)
                    do jb = 1, workspace%tens2_dims(2)
                       do ja = 1, workspace%tens2_dims(1)
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

                 if (preconditioner_debug) then
                    write (6,*) 'debug >>> precon1_tens  root',iroot,' >>>'
                    b = 0
                    do ja = 1, workspace%tens1_dims(1)
                       b = b + 1
                       write (6,*) b, preco_1d_p(ja)
                    end do
                    write (6,*) 'debug <<< precon1_tens <<<'

                    write (6,*) 'debug >>> precon2_tens  root',iroot,' >>>'
                    do jc = 1, workspace%tens2_dims(3)
                       do jb = 1, workspace%tens2_dims(2)
                          do ja = 1, workspace%tens2_dims(1)
                              b = b + 1
                              write (6,*) b, ja, jb, jc, preco_3d_p(ja,jb,jc)   
                          end do
                       end do
                    end do
                    write (6,*) 'debug <<< precon2_tens <<<'
                 end if

              else if (workspace%tens1_rank .eq. 2) then
                 call c_f_pointer(preco1_body_p,preco_2d_p,workspace%tens1_dims(1:workspace%tens1_rank))
                 call c_f_pointer(preco2_body_p,preco_4d_p,workspace%tens2_dims(1:workspace%tens2_rank))

                 do jb = 1, workspace%tens1_dims(2)
                    do ja = 1, workspace%tens1_dims(1)
                       proj_r = eigenvalue_r - real(hdiag_2d_p(ja,jb))
                       proj_i = eigenvalue_i - aimag(hdiag_2d_p(ja,jb))
! aspg: added this check to filter out real or imaginary parts that are too small/noisy and would explode as denominators 
!       worst case this will set the preconditioner of this element to one (due to the initialization, see above)
                       if (dabs(proj_r).le.small_number) proj_r = 0.0d0
                       if (dabs(proj_i).le.small_number) proj_i = 0.0d0

                       x2py2 = proj_r*proj_r + proj_i*proj_i
                       if (dabs(x2py2).gt.safe_minimum) then
                          proj_r =  proj_r / x2py2
                          proj_i = -proj_i / x2py2
                          preco_2d_p(ja,jb) = cmplx(proj_r,proj_i)
                       end if
                    end do
                 end do

                 do jd = 1, workspace%tens2_dims(4)
                    do jc = 1, workspace%tens2_dims(3)
                       do jb = 1, workspace%tens2_dims(2)
                          do ja = 1, workspace%tens2_dims(1)
                             proj_r = eigenvalue_r - real(hdiag_4d_p(ja,jb,jc,jd))
                             proj_i = eigenvalue_i - aimag(hdiag_4d_p(ja,jb,jc,jd))
! aspg: added this check to filter out real or imaginary parts that are too small/noisy and would explode as denominators 
!       worst case this will set the preconditioner of this element to one (due to the initialization, see above)
                             if (dabs(proj_r).le.small_number) proj_r = 0.0d0
                             if (dabs(proj_i).le.small_number) proj_i = 0.0d0

                             x2py2 = proj_r*proj_r + proj_i*proj_i
                             if (dabs(x2py2).gt.safe_minimum) then
                                proj_r =  proj_r / x2py2
                                proj_i = -proj_i / x2py2
                                preco_4d_p(ja,jb,jc,jd) = cmplx(proj_r,proj_i)
                             end if
                          end do
                       end do
                    end do
                 end do

                 if (preconditioner_debug) then
                    write (6,*) 'debug >>> precon1_tens  root',iroot,' >>>'
                    b = 0
                    do jb = 1, workspace%tens1_dims(2)
                       do ja = 1, workspace%tens1_dims(1)
                          b = b + 1
                          write (6,*) b, ja, jb, preco_2d_p(ja,jb)
                       end do
                    end do
                    write (6,*) 'debug <<< precon1_tens <<<'

                    write (6,*) 'debug >>> precon2_tens  root',iroot,' >>>'
                    do jd = 1, workspace%tens2_dims(4)
                       do jc = 1, workspace%tens2_dims(3)
                          do jb = 1, workspace%tens2_dims(2)
                             do ja = 1, workspace%tens2_dims(1)
                                 b = b + 1
                                 write (6,*) b, ja, jb, jc, jd, preco_4d_p(ja,jb,jc,jd) 
                             end do
                          end do
                       end do
                    end do
                    write (6,*) 'debug <<< precon2_tens <<<'
                 end if

              end if
           end do ! iroot
       end if ! create preconditioner

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
!           do_right = 'N'
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
 

    subroutine talsh_davidson_workspace_destroy(workspace)
        type(talsh_davidson_workspace),  intent(inout) :: workspace
        integer :: ierr, i

      
        do i = 1, workspace%Lmax
           ierr = talsh_tensor_destruct(workspace%trial1(i))
           ierr = talsh_tensor_destruct(workspace%trial2(i))
           ierr = talsh_tensor_destruct(workspace%sigma1(i))
           ierr = talsh_tensor_destruct(workspace%sigma2(i))
        end do
        deallocate(workspace%trial1)
        deallocate(workspace%trial2)
        deallocate(workspace%sigma1)
        deallocate(workspace%sigma2)

        do i = 1, workspace%nroots
           ierr = talsh_tensor_destruct(workspace%delta1(i))
           ierr = talsh_tensor_destruct(workspace%delta2(i))
           ierr = talsh_tensor_destruct(workspace%precond1(i))
           ierr = talsh_tensor_destruct(workspace%precond2(i))
           ierr = talsh_tensor_destruct(workspace%CTtrial1(i))
           ierr = talsh_tensor_destruct(workspace%CTtrial2(i))
        end do
        deallocate(workspace%delta1)
        deallocate(workspace%delta2)
        deallocate(workspace%precond1)
        deallocate(workspace%precond2)
        deallocate(workspace%CTtrial1)
        deallocate(workspace%CTtrial2)

        deallocate(workspace%tens1_dims)
        deallocate(workspace%tens2_dims)

        ierr = talsh_tensor_destruct(workspace%Hdiag1)
        ierr = talsh_tensor_destruct(workspace%Hdiag2)

        if (associated(workspace%N)) deallocate(workspace%N)
        if (associated(workspace%M)) deallocate(workspace%M)

        write(davidson_output_file_unit, "(4x,a)") 'leaving talsh_davidson_workspace_destroy'
    end subroutine 

    subroutine talsh_davidson_workspace_create(workspace, nocc, nvirt, flavor, excitation_level)
        type(talsh_davidson_workspace),  intent(inout) :: workspace
        character*2 :: flavor
        integer, intent(in) :: nocc, nvirt, excitation_level

        integer :: M_dims(2) 
        integer :: ierr, i
!       
!   for the moment we force the number of roots to be 4
!
        workspace%nroots = 2
        workspace%Lmax   = davidson_max_subspace_size
        workspace%L      = workspace%nroots 
        workspace%Lprev  = 0

        allocate(workspace%N(excitation_level))

        allocate(workspace%trial1(workspace%Lmax))
        allocate(workspace%trial2(workspace%Lmax))

        allocate(workspace%CTtrial1(workspace%nroots))
        allocate(workspace%CTtrial2(workspace%nroots))

        allocate(workspace%sigma1(workspace%Lmax))
        allocate(workspace%sigma2(workspace%Lmax))

        allocate(workspace%delta1(workspace%nroots))
        allocate(workspace%delta2(workspace%nroots))

        allocate(workspace%precond1(workspace%nroots))
        allocate(workspace%precond2(workspace%nroots))

!       M_dims     = (/ workspace%Lmax, workspace%Lmax /)
!       ierr=talsh_tensor_construct(workspace%M,C8,M_dims,init_val=ZERO)
        allocate(workspace%M(workspace%Lmax,workspace%Lmax))
        workspace%M = ZERO
        
        select case(flavor)
           case("IP")
! singles
              workspace%tens1_rank = 1
              allocate(workspace%tens1_dims(workspace%tens1_rank))
              workspace%tens1_dims = (/ nocc /)
              workspace%N(1) = nocc 
! doubles
              workspace%tens2_rank = 3
              allocate(workspace%tens2_dims(workspace%tens2_rank))
              workspace%tens2_dims = (/ nvirt, nocc, nocc /)
              workspace%N(2) = nvirt*nocc*nocc

           case("EA")
! singles
              workspace%tens1_rank = 1
              allocate(workspace%tens1_dims(workspace%tens1_rank))
              workspace%tens1_dims = (/ nvirt /)
              workspace%N(1) = nvirt 
! doubles
              workspace%tens2_rank = 3
              allocate(workspace%tens2_dims(workspace%tens2_rank))
              workspace%tens2_dims = (/ nvirt, nvirt, nocc /)
              workspace%N(2) = nvirt*nvirt*nocc

           case("EE")
! singles
              workspace%tens1_rank = 2
              allocate(workspace%tens1_dims(workspace%tens1_rank))
              workspace%tens1_dims = (/ nvirt, nocc /)
              workspace%N(1) = nvirt*nocc
! doubles
              workspace%tens1_rank = 4
              allocate(workspace%tens2_dims(workspace%tens1_rank))
              workspace%tens2_dims = (/ nvirt, nvirt, nocc, nocc /)
              workspace%N(2) = nvirt*nvirt*nocc*nocc

        end select

! singles
        do i = 1, workspace%Lmax
           ierr = talsh_tensor_construct(workspace%trial1(i),C8,workspace%tens1_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace%sigma1(i),C8,workspace%tens1_dims,init_val=ZERO)
        end do
        do i = 1, workspace%nroots
           ierr = talsh_tensor_construct(workspace%delta1(i),C8,workspace%tens1_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace%precond1(i),C8,workspace%tens1_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace%CTtrial1(i),C8,workspace%tens1_dims,init_val=ZERO)
        end do
        ierr = talsh_tensor_construct(workspace%Hdiag1,C8,workspace%tens1_dims,init_val=ZERO)

! doubles
        do i = 1, workspace%Lmax
           ierr=talsh_tensor_construct(workspace%trial2(i),C8,workspace%tens2_dims,init_val=ZERO)
           ierr=talsh_tensor_construct(workspace%sigma2(i),C8,workspace%tens2_dims,init_val=ZERO)
        end do
        do i = 1, workspace%nroots
           ierr = talsh_tensor_construct(workspace%delta2(i),C8,workspace%tens2_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace%precond2(i),C8,workspace%tens2_dims,init_val=ZERO)
           ierr = talsh_tensor_construct(workspace%CTtrial2(i),C8,workspace%tens2_dims,init_val=ZERO)
        end do
        ierr=talsh_tensor_construct(workspace%Hdiag2,C8,workspace%tens2_dims,init_val=ZERO)

        workspace%Nmax = workspace%N(1) + workspace%N(2)


! debug statements
        write(davidson_output_file_unit, "(4x,a)")  'in talsh_davidson_workspace_create'

        write(davidson_output_file_unit, "(4x,a,i6,a,i6)") &
     &  'created tensor for trial vectors, dimensions :',workspace%Nmax,"x",workspace%Lmax


        do i = 1, workspace%Lmax
           write(davidson_output_file_unit, "(4x,a,i4)") 'info on workspace%trial1, element',i
           call talsh_tensor_print_info(workspace%trial1(i))
           write(davidson_output_file_unit, "(4x,a,i4)") 'info on workspace%trial2, element',i
           call talsh_tensor_print_info(workspace%trial2(i))
        end do

        write(davidson_output_file_unit, "(4x,a)") 'info on workspace%Hdiag1' 
        call talsh_tensor_print_info(workspace%Hdiag1)
        write(davidson_output_file_unit, "(4x,a)") 'info on workspace%Hdiag2' 
        call talsh_tensor_print_info(workspace%Hdiag2)

        write(davidson_output_file_unit, "(4x,a,i6,a,i6)") &
     &  'created tensor for sigma vectors, dimensions :',workspace%Nmax,"x",workspace%Lmax

        write(davidson_output_file_unit, "(4x,a,i6,a,i6)") &
     &  'created tensor for correction vectors, dimensions :',workspace%Nmax,"x",workspace%nroots

        
   end subroutine 


   subroutine create_sigma_right(exa_input,int_t,eom_t,t2_tensor,workspace)
      use talsh_eom_sigma, only : sigma_right_ip
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_davidson_workspace), intent(inout) :: workspace

      complex(8), pointer, contiguous:: s1_h_p_p(:)
      type(C_PTR):: s1_body_p

      complex(8), pointer, contiguous:: t1_h_p_p(:)
      type(C_PTR):: t1_body_p

      complex(8), pointer, contiguous:: t2_h_p_p(:,:,:)
      type(C_PTR):: t2_body_p

      complex(8), pointer, contiguous:: s2_h_p_p(:,:,:)
      type(C_PTR):: s2_body_p

      integer :: ierr, iroot, i, j, a, b

      do iroot = (workspace%L - workspace%nroots + 1), workspace%L

         write (6,*) 'root ',iroot
! aspg
! debug info
         ierr = talsh_tensor_get_body_access(workspace%sigma1(iroot),s1_body_p,C8,0,DEV_HOST)
         write (6,*) 's1_body_p',s1_body_p
         ! associating the pointer to sigma1 to the local, fortran s1_p pointer variable
         call c_f_pointer(s1_body_p,s1_h_p_p,workspace%tens1_dims(1:workspace%tens1_rank))

         ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),t1_body_p,C8,0,DEV_HOST)
         write (6,*) 't1_body_p',t1_body_p
         call c_f_pointer(t1_body_p,t1_h_p_p,workspace%tens1_dims(1:workspace%tens1_rank))

         ierr = talsh_tensor_get_body_access(workspace%sigma2(iroot),s2_body_p,C8,0,DEV_HOST)
         write (6,*) 's2_body_p',s2_body_p
         ! associating the pointer to sigma1 to the local, fortran s1_p pointer variable
         call c_f_pointer(s2_body_p,s2_h_p_p,workspace%tens2_dims(1:workspace%tens2_rank))

         ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),t2_body_p,C8,0,DEV_HOST)
         write (6,*) 't2_body_p',t2_body_p
         call c_f_pointer(t2_body_p,t2_h_p_p,workspace%tens2_dims(1:workspace%tens2_rank))
! debug info
! aspg
! aspg
! debug info     
         write (6,*) 'debug >>> t1_tens before sigma, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         do i = 1, workspace%tens1_dims(1)
            write (6,*) i, t1_h_p_p(i)
         end do
         write (6,*) 'debug <<< t1_tens before sigma <<<'
         write (6,*) 'debug >>> t2_tens before sigma, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         do i = 1, workspace%tens2_dims(3)
            do j = 1, workspace%tens2_dims(2)
               do a = 1, workspace%tens2_dims(1)
                  write (6,*) a, j, i, t2_h_p_p(a,j,i)
               end do
            end do
         end do
         write (6,*) 'debug <<< t2_tens before sigma <<<'

         write (6,*) 'debug >>> s1_tens before sigma, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         b = 0
         do i = 1, workspace%tens1_dims(1)
            b = b + 1
            write (6,*) b, s1_h_p_p(i)
         end do
         write (6,*) 'debug <<< s1_tens before sigma <<<'
         write (6,*) 'debug >>> s2_tens before sigma, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         do i = 1, workspace%tens2_dims(3)
            do j = 1, workspace%tens2_dims(2)
               do a = 1, workspace%tens2_dims(1)
                  b = b + 1
                  write (6,*) b, a, j, i, s2_h_p_p(a,j,i)
               end do
            end do
         end do
         write (6,*) 'debug <<< s2_tens before sigma <<<'
! debug info     
! aspg


         call sigma_right_ip(exa_input,int_t,eom_t,t2_tensor, &
   &                         workspace%trial1(iroot),workspace%trial2(iroot),&
   &                         workspace%sigma1(iroot),workspace%sigma2(iroot))

! aspg
! debug info     
         write (6,*) 'debug >>> t1_tens, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         b = 0
         do i = 1, workspace%tens1_dims(1)
            b = b + 1
            write (6,*) b, t1_h_p_p(i)
         end do
         write (6,*) 'debug <<< t1_tens <<<'
         write (6,*) 'debug >>> t2_tens, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         do i = 1, workspace%tens2_dims(3)
            do j = 1, workspace%tens2_dims(2) 
               do a = 1, workspace%tens2_dims(1)
                  b = b + 1
                  write (6,*) b, a, j, i, t2_h_p_p(a,j,i)
               end do
            end do
         end do
         write (6,*) 'debug <<< t2_tens <<<'

         write (6,*) 'debug >>> s1_tens, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         b = 0
         do i = 1, workspace%tens1_dims(1)
            b = b + 1
            write (6,*) b, s1_h_p_p(i)
         end do
         write (6,*) 'debug <<< s1_tens <<<'
         write (6,*) 'debug >>> s2_tens, root',iroot,' >>>'
!        call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
         do i = 1, workspace%tens2_dims(3)
            do j = 1, workspace%tens2_dims(2)
               do a = 1, workspace%tens2_dims(1)
                  b = b + 1
                  write (6,*) b, a, j, i, s2_h_p_p(a,j,i)
               end do
            end do
         end do
         write (6,*) 'debug <<< s2_tens <<<'
! debug info     
! aspg

      end do

!     call quit('bb')

   end subroutine


   subroutine create_subspace_matrix_right(M_small,workspace)
      type(talsh_davidson_workspace), intent(inout) :: workspace ! here we will keep the M(Lmax,Lmax)
      complex(kind=8), pointer, contiguous, intent(inout) :: M_small(:,:)    ! this will be a copy of M, M_small(L,L)

      integer :: mij_dims(2), mij_rank
      integer(INTD)      :: one_dims(1)
      type(talsh_tens_t) :: result_tensor
      complex(8), pointer :: mij_result_tens(:)
      type(C_PTR):: mij_body_p

      integer :: m_dims(2), m_rank
      complex(8), pointer, contiguous:: m_tens(:,:)
      type(C_PTR):: m_body_p
      integer :: ierr, i, j

!! we need to fish out the small matrix from workspace%M, which is the slice LxL of the full tensor
!! we also need to conjugate transpose the trial vector 
      write (6,*) 'workspace%L',workspace%L

      if (associated(M_small)) deallocate(M_small)
      allocate(M_small(workspace%L,workspace%L))
      M_small = ZERO 

!     Initialize scalars that are to be used as tensors in contractions
      one_dims(1) = 1
!     ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
      ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
      ierr=talsh_tensor_get_body_access(result_tensor,mij_body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(mij_body_p,mij_result_tens,one_dims)


! capital letters here are meant to be the roots
      write (6,*) 'foobar, before setup M'
      write (6,*) 'foobar, workspace%L',workspace%L
      write (6,*) 'foobar, workspace%Lprev',workspace%Lprev
      do i = workspace%Lprev + 1, workspace%L
         do j = workspace%Lprev + 1, workspace%L
! this will be a reduction to a scalar at each contractions, and we then will build the 
            ierr=talsh_tensor_contract("MIJ()+=Bdag+(e)*S(e)",               &
   &                                                    result_tensor,      &
   &                                                    workspace%trial1(i),&
   &                                                    workspace%sigma1(j))
            ierr=talsh_tensor_contract("MIJ()+=Bdag+(a,b,c)*S(a,b,c)",       &
   &                                                    result_tensor,      &
   &                                                    workspace%trial2(i),&
   &                                                    workspace%sigma2(j))
            write (6,*) 'i,j, m_ij:', i,j, mij_result_tens(1)
            workspace%M(i,j) = mij_result_tens(1) 
            mij_result_tens(1) = ZERO
         end do
      end do

      M_small = workspace%M(1:workspace%L,1:workspace%L)

      write (6,*) 'M small'
      do i = 1, workspace%L
         do j = 1, workspace%L
            write (6,*) M_small(i,j)
         end do
      end do

!      call quit('bb')

    end


! aspg, 16/04/2019
! for the moment we have a long argument list to get it going, later we try to shorten in
!   subroutine talsh_get_values_vectors_davidson(A, workspace, projector)
   subroutine talsh_get_values_vectors_davidson(exa_input,int_t,eom_t,t2_tensor, workspace, projector)
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_davidson_workspace), intent(inout) :: workspace
      type(talsh_projectors), intent(inout) :: projector
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal

!       type(talsh_lambda_tens), intent(in)   :: A
!       type(talsh_davidson_workspace),  intent(inout) :: workspace
!       type(talsh_projectors), intent(inout) :: projector
!
! nullifying all pointers at the start to give them a state, see discussion on http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html
! "[...] When a pointer is declared its status is undefined, and cannot be safely queried with the associated intrinsic. [...]"
! 
        complex(kind=8), pointer, contiguous :: M_small(:,:)        => NULL()
        complex(kind=8), pointer, contiguous :: c_small_r(:,:)      => NULL()
        complex(kind=8), pointer, contiguous :: c_small_l(:,:)      => NULL()
        complex(kind=8), pointer, contiguous :: e_small(:)          => NULL()

        complex(kind=8), pointer, contiguous :: c_small_r_old(:,:)  => NULL()
        complex(kind=8), pointer, contiguous :: c_small_l_old(:,:)  => NULL()
        logical, pointer :: delta(:) => NULL()

        integer :: niter
        integer :: i, j, k
        integer :: nroots, ierr

        nroots    = workspace%nroots
!
! davidson_create_results() will allocate storage for results from davidson diagonalization when called a first time. 
! 
!!      call davidson_create_results(results,vectors_N,nroots,&
!!                                   davidson_solve_values,   &
!!                                   davidson_solve_lhs,      &
!!                                   davidson_solve_rhs)

        workspace%L = nroots

        write(davidson_output_file_unit,"(4x,a,i6,a)") "Requested ",nroots," eigenvales"
        if ( nroots .gt. workspace%Nmax ) then 
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested larger than matrix order !"
           return 
        else if ( nroots .eq. workspace%Nmax ) then
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested equal to matrix order !"
        endif

        
        write(davidson_output_file_unit,*) ""
        do niter = 1, davidson_max_iterations
            write(davidson_output_file_unit,*) "Iteration ", niter
            write(davidson_output_file_unit,*) ""

            if (davidson_solve_rhs) then
               call create_sigma_right(exa_input,int_t,eom_t,t2_tensor,workspace)
               call create_subspace_matrix_right(M_small,workspace) 
            else
               stop 'unknown choice of side to solve for in davidson'  
            end if

!!          if (associated(projector)) call apply_projector(projector,sigma_vector)

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
                                                        verbose=.true.)
!                                                       verbose=davidson_verbose)

            if (davidson_solve_rhs.and.davidson_overlap_sorting) then
                c_small_r_old = c_small_r
            end if
            if (davidson_solve_lhs.and.davidson_overlap_sorting) then
                c_small_l_old = c_small_l
            end if
 
            call form_preconditioner(workspace, e_small)

            if (associated(delta)) deallocate(delta)
           if (davidson_solve_rhs) delta => correction_to_vectors(davidson_output_file_unit, nroots,  &
                                                                  e_small, c_small_r, workspace)

            if (.not.associated(delta)) then
                write(davidson_output_file_unit,*) "All roots converged."
            end if

! not converged yet, if asked to refresh the trial vectors do it now, otherwise create
! new trial vectors
            if (davidson_refresh_trial_vectors.and.(mod(niter,davidson_trial_vector_refresh_rate).eq.0)) then
            else
               call update_trial_vectors(delta,workspace)
            end if
 
            if (workspace%L .ge. workspace%Lmax) then
               write(davidson_output_file_unit, *) "Reached maximum subspace size, stopping."
               exit
            end if
        end do ! davidson iterations

        if (davidson_solve_rhs.and.davidson_overlap_sorting) deallocate(c_small_r_old)
        if (davidson_solve_lhs.and.davidson_overlap_sorting) deallocate(c_small_l_old)

        write(davidson_output_file_unit, *) "Final size of reduced subspace : ",workspace%L

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


    subroutine setup_overlap_sorting(c_r_old, c_l_old, nroots, solve_rhs, solve_lhs)
        real(kind=8), intent(inout), pointer :: c_r_old(:,:) 
        real(kind=8), intent(inout), pointer :: c_l_old(:,:)  
        integer, intent(in) :: nroots
        logical, intent(in) :: solve_rhs, solve_lhs

        integer :: i, factor 
!       real(kind=8), pointer :: c_small_r_old_dum(:,:)  => NULL()

        if (davidson_complex_mode_with_reals) then
            factor = 2 
        else
            factor = 1 
        end if

        if (solve_rhs) then
            allocate(c_r_old(factor*nroots,nroots))
            c_r_old = 0.0d0
            do i = 1, nroots*factor, factor
                c_r_old(i,(i+1)/2) = 1.0d0
            enddo
        endif
        if (solve_lhs) then
            allocate(c_l_old(factor*nroots,nroots))
            c_l_old = 0.0d0
            do i = 1, nroots*factor, factor
                c_l_old(i,(i+1)/2) = 1.0d0
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
            follow_overlap_scheme_local = 1 ! 1 = scheme calculating overlap matrix ; 2 = scheme not explicitly calculating overlap matrix
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

        call talsh_solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left_local)

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
                    call sort_by_overlap_of_subspace_left_or_right(eValues, eVectorsR)
                else if (left_local) then
                    call sort_by_overlap_of_subspace_left_or_right(eValues, eVectorsL)
                end if
            end if ! selection of overlap sorting scheme 

        else
            call talsh_sort_values_vectors(sorting_order_local, eValues, eVectorsR, eVectorsL, &
                                     right_local, left_local, apply_evalues_shift, eValuesShift)
        end if

        call print_values_vectors(davidson_output_file_unit,eValues, eVectorsR, eVectorsL, &
                                    .true., .true., nroots=nroots_to_print)

    end subroutine

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

        real(kind=8), allocatable :: temp_vector(:,:)
        real(kind=8) :: temp
        integer, allocatable :: indices_sorted_array(:)
        real(kind=8), allocatable :: e_small_real(:,:)
        integer :: vectors_N, vectors_L, i_sorted, i

        vectors_L = size(e_small)
        allocate(e_small_real(2,vectors_L))
        allocate(indices_sorted_array(vectors_L))
        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
            e_small_real(1,i) = realpart(e_small(i))
            e_small_real(2,i) = imagpart(e_small(i))
        end do

        if (apply_shift) then
           call quicksort_vector_complex_as_real_shifted(energy_shift, e_small_real, indices_sorted_array, 1, vectors_L)
        else
           call quicksort_vector_complex_as_real(e_small_real, indices_sorted_array, 1, vectors_L)
        end if

        do i = 1, vectors_L
            e_small(i) = complex(e_small_real(1,i), e_small_real(2,i))
        end do

        allocate(temp_vector(vectors_L,vectors_L))
        temp_vector = 0.0d0

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
              overlap_copy(1,i) = realpart(overlap(iroot,i))
              overlap_copy(2,i) = imagpart(overlap(iroot,i))
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

    subroutine sort_by_overlap_of_subspace_left_or_right(e_small, c_small)
!
! this routine sorts the eigenvectors with respect to the overlap between them
! so that the diagonal of the sorted eigenvector matrix contains the largest
! element of each row.
!
! this should be equivalent to sorting the eigenpairs in descending order of overlap
! with a set of initial trial vectors V = {v_1, v_2, v_3, ...} since, following the
! reasoning by Zuev et al, J. Comput. Chem 2014 (DIO:10.1002/jcc23800):
!
! a. considering the set of new vectors for the full problem 
!      U = { u_1, u_2, u_3, ...} = { V c_small_1, V c_small_2, V c_small_3, ...}  
!    
!    the overlap of the initial and new vectors is given as:
!      S_{1;k} = < v_1 | u_k > = < v_1 | V c_small_k >
!      S_{2;k} = < v_2 | u_k > = < v_2 | V c_small_k >
!      [...]
!
!    where c_small_k is the k-th element of c_small(=the eigenvectors for the reduced subspace)
!
! b. we have that the overlap matrix elements are given by
!      S_{1;k} = < v_1 |  V c_small_k > = \sum_j < v_1 |  v_j > c_small_{j;k} = \delta_{1j} c_small_{j;k} = c_small_{1;k}
!      S_{2;k} = < v_2 |  V c_small_k > = \sum_j < v_2 |  v_j > c_small_{j;k} = \delta_{2j} c_small_{j;k} = c_small_{2;k}
!      [...]
!
!    as, due to the orthonormality of the trial vectors in V, the overlap for vectors in the full bases
!    can be obtained from the vectors in the reduced subspace. 
!   
! c. now, in order to find the maximum overlap with each of the trial vectors, we have to identify
!    the largest element of each corresponding line in the reduced vectors.
! 
!    therefore, in the routine below, by reshuffling the columns of c_small so the maximum elements of 
!    each row of c_small is placed in the diagonal, we achive a mapping that places the maximum overlap for v_1 
!    in column 1, for v_2 in column 2 and so on.
!
! NOTE: at this point nothing is done to ensure that no vector is placed twice, or removed from c_small altogether;
! something that may be an issue if, for instance, the largest element on two or more lines are found for the
! same column (so we'd end up repeating the same column a number of times) - though even in the case of using the
! same c_small_k, one would ! still have linearly independent u vectors since 
!    u_1 = v_1 c_small_k 
!    u_2 = v_2 c_small_k
!
! that said, a possible wat forward if this becomes an issue is would be to use the next to largest elements that haven't
! been used before.
!
        complex(kind=8), intent(inout) :: c_small(:,:) ! left or right-side reduced subspace eigenvectors
        complex(kind=8), intent(inout) :: e_small(:) ! reduced subspace eigenvalues

        complex(kind=8),allocatable :: temp(:), temp_vector(:,:)
        integer, allocatable :: indices_sorted_array(:), already_sorted(:)
        integer :: vectors_L, vectors_N, i_sorted, i, max_overlap(1), iroot, jroot
        real(kind=8),allocatable :: norms(:)
        double precision dnrm2, dznrm2


        vectors_L = size(c_small,2)

        allocate(indices_sorted_array(vectors_L))
        allocate(already_sorted(vectors_L))
        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
        end do
        already_sorted = -1

        allocate(norms(vectors_L))
        allocate(temp(vectors_L))
        allocate(temp_vector(vectors_L,vectors_L))

        do iroot = 1, vectors_L
! the loop over jroot will create a temporary array with the norm of the iroot-th roow of the eigenvectors
          do jroot = 1, vectors_L
              norms(jroot) = dznrm2(vectors_L, c_small(:,jroot), 1) 
          end do
          max_overlap = maxloc(norms)
          i = max_overlap(1)
          indices_sorted_array(iroot) = i
        enddo

        do i = 1, vectors_L
            i_sorted = indices_sorted_array(i)
            temp_vector(:,i) = c_small(:,i_sorted)
        end do
        c_small = temp_vector

        do i = 1, vectors_L
            i_sorted = indices_sorted_array(i)
            temp(i) = e_small(i_sorted)
        end do
        e_small = temp

        deallocate(norms)
        deallocate(temp)
        deallocate(temp_vector)
        deallocate(indices_sorted_array)

    end subroutine

    subroutine apply_projector(projector,b)
        logical, intent(in) :: projector(:)
        real(kind=8), intent(inout)  :: b(:,:)
        integer :: vectors_L, vectors_N, i, j, k 

! the values of projector mean 
!  .true. : the entry will remain as is
! .false. : the entre will be set to zero

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (davidson_complex_mode_with_reals) vectors_N = vectors_N / 2

        do i = 1, vectors_N
           if (projector(i).eqv..false.) then
              if (davidson_complex_mode_with_reals) then
                 j = 2*(i-1) + 1
                 k = j + 1
                 b(j,:) = 0.0d0
                 b(k,:) = 0.0d0
              else 
                 b(i,:) = 0.0d0
              end if
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
         do i = 1, a_dims2(1)
           do j = 1, a_dims2(2)
             a_tens2(i,j)=a_tens2(i,j)*b_tens2(i,j)
           end do
         end do
    
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
         do i = 1, a_dims4(1)
            do j = 1, a_dims4(2)
              do k = 1, a_dims4(3)
                do l = 1, a_dims4(4)
                   a_tens4(i,j,k,l)=a_tens4(i,j,k,l)*b_tens4(i,j,k,l)
                end do
             end do
           end do
         end do
    
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

