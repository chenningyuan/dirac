module mfsolver_cfg

implicit none

private 

public mfsolver_setup
public mfsolver_config_t

type mfsolver_config_t

   real(kind=8) :: safe_minimum

   real(kind=8) :: convergence_threshold 
   real(kind=8) :: antisymme_threshold
   real(kind=8) :: eigenvalues_shift(2)

   integer      :: max_subspace_size 
   integer      :: max_iterations
   integer      :: output_file_unit 
   integer      :: subspaces_file_unit 
   integer      :: sigmas_file_unit 
   integer      :: results_file_unit 
   integer      :: overlap_sorting_scheme
   integer      :: trial_vector_refresh_rate
   integer      :: maximum_trial_vector_refresh_rate 

   logical      :: apply_eigenvalues_shift
   logical      :: solve_rhs 
   logical      :: solve_lhs 
   logical      :: solve_values
   logical      :: solve_linear_system
   logical      :: verbose 
   logical      :: symmetric 
   logical      :: complex_mode_with_reals 
   logical      :: save_subspaces 
   logical      :: save_sigmas 
   logical      :: save_results
   logical      :: overlap_sorting
   logical      :: debug 
   logical      :: use_preconditioner
   logical      :: refresh_trial_vectors 

! 
! variables for controlling i/o
!
   character(len=16) :: left_sigma_file_name 
   character(len=16) :: left_vectors_file_name
   character(len=16) :: left_subspace_file_name

   character(len=16) :: right_sigma_file_name
   character(len=16) :: right_vectors_file_name
   character(len=16) :: right_subspace_file_name

   integer :: right_sigma_file_unit
   integer :: right_vectors_file_unit
   integer :: right_subspace_file_unit

   integer :: left_sigma_file_unit
   integer :: left_vectors_file_unit
   integer :: left_subspace_file_unit

   real(kind=8) :: eigenvalues_shift_threshold
   integer :: max_label_length
   integer :: nr_max_irreps

end type mfsolver_config_t

contains 

    subroutine mfsolver_setup(config,&
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

        type(mfsolver_config_t), intent(inout) :: config
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

! defaults for variables

        config%eigenvalues_shift_threshold = 1.0e-3
        config%max_label_length = 30
        config%nr_max_irreps = 32
        config%safe_minimum = 1.0d-36


        config%debug = .false.
        config%use_preconditioner = .true.
        config%refresh_trial_vectors = .false.
        config%maximum_trial_vector_refresh_rate = 10

        config%left_sigma_file_name     = 'LSIGMAV.sym'
        config%left_vectors_file_name   = 'LEVECTO.sym'
        config%left_subspace_file_name  = 'LSUBSPA.sym'

        config%right_sigma_file_name    = 'RSIGMAV.sym'
        config%right_vectors_file_name  = 'REVECTO.sym'
        config%right_subspace_file_name = 'RSUBSPA.sym'

        config%right_sigma_file_unit    = 140
        config%right_vectors_file_unit  = 141
        config%right_subspace_file_unit = 142

        config%left_sigma_file_unit     = 240
        config%left_vectors_file_unit   = 241
        config%left_subspace_file_unit  = 242

! now handling variables

        call set_safe_minimum(config%safe_minimum,config%debug) 

        if (present(convergence_threshold)) then
            config%convergence_threshold = convergence_threshold
        else
            config%convergence_threshold = 1.0d-8 
        end if

        if (config%convergence_threshold.lt.1.0d-8) then
            config%antisymme_threshold = config%convergence_threshold * 0.1
        else
            config%antisymme_threshold = 1.0d-9
        end if 


        if (present(energy_shift)) then
            config%eigenvalues_shift = energy_shift
            config%apply_eigenvalues_shift = .true.

            if (abs(config%eigenvalues_shift(1)).lt.config%eigenvalues_shift_threshold) &
               config%eigenvalues_shift(1) = 0.0d0 
 
            if (abs(config%eigenvalues_shift(2)).lt.config%eigenvalues_shift_threshold) &
               config%eigenvalues_shift(2) = 0.0d0 

            if ((abs(config%eigenvalues_shift(1)).lt.config%eigenvalues_shift_threshold).and. &
                (abs(config%eigenvalues_shift(2)).lt.config%eigenvalues_shift_threshold)) then
               config%apply_eigenvalues_shift = .false.
            end if
        else
            config%eigenvalues_shift = (/0.0d0, 0.0d0/) 
            config%apply_eigenvalues_shift = .false.
        end if

        if (present(max_subspace_size)) then
            config%max_subspace_size = max_subspace_size  
        else
            config%max_subspace_size = 2048
        end if

        if (present(max_iterations)) then
            config%max_iterations = max_iterations
        else
            config%max_iterations = 100
        end if

        if (present(solve_values)) then
            config%solve_values = solve_values
        else
            config%solve_values = .true. 
        end if

        if (present(solve_linear_sys)) then
            config%solve_linear_system = solve_linear_sys
            if (config%solve_linear_system) config%solve_values = .false.
        else
            config%solve_linear_system = .false. 
        end if
        if (config%solve_linear_system.and.config%solve_values) then
           call quit('inconsistent setup. cannot calculate eigenvalues and solve linear system')
        endif
 

        if (present(solve_left)) then
            config%solve_lhs = solve_left 
        else
            config%solve_lhs = .false. 
        end if

        if (present(solve_right)) then
            config%solve_rhs = solve_right
        else
            config%solve_rhs = .true. 
        end if

        if (present(overlap_sorting_scheme)) then
            config%overlap_sorting_scheme = overlap_sorting_scheme
        else
            config%overlap_sorting_scheme = 1 ! 1 = scheme calculating overlap matrix ; 2 = scheme not explicitly calculating overlap matrix
        end if

        if (present(overlap_sorting)) then
            config%overlap_sorting = overlap_sorting
        else
            config%overlap_sorting = .false.
        end if

        if (present(use_preconditioner)) then
            config%use_preconditioner = use_preconditioner
        else
            config%use_preconditioner = .true.
        end if

        if (present(verbose)) then
            config%verbose = verbose
            config%debug   = verbose
        else
            config%verbose = .false.
        end if

        if (present(symmetric)) then
            config%symmetric = symmetric
        else
            config%symmetric = .false. 
        end if

        if (present(complex_mode)) then
            config%complex_mode_with_reals = complex_mode 
        else
            config%complex_mode_with_reals = .false. 
        end if

        if (present(refresh_trial_rate)) then
            if (refresh_trial_rate.le.0) then
                config%refresh_trial_vectors = .false.
            else 
                config%refresh_trial_vectors = .true.
                if (refresh_trial_rate.le.config%maximum_trial_vector_refresh_rate) then
                    config%trial_vector_refresh_rate = config%maximum_trial_vector_refresh_rate 
                else
                    config%trial_vector_refresh_rate = refresh_trial_rate
                end if
            end if
        else
            config%refresh_trial_vectors = .false.
        end if

        if (present(unit_output)) then
            config%output_file_unit = unit_output
        end if

        if (present(save_subspaces)) then
           config%save_subspaces = save_subspaces
        else
           config%save_subspaces = .false. 
        end if

        if (present(save_sigmas)) then
           config%save_sigmas = save_sigmas
        else
           config%save_sigmas = .false. 
        end if

        if (present(save_results)) then
           config%save_results = save_results
        else
           config%save_results = .false. 
        end if

        if (present(unit_subspaces)) then
            config%subspaces_file_unit = unit_subspaces
        else
            config%subspaces_file_unit = 137 
        end if

        if (present(unit_sigmas)) then
            config%sigmas_file_unit = unit_sigmas
        else
            config%sigmas_file_unit = 138
        end if

        if (present(unit_results)) then
            config%results_file_unit = unit_results
        else
            config%results_file_unit = 138
        end if

        if (config%complex_mode_with_reals.and.config%overlap_sorting) then
            config%overlap_sorting_scheme = 2
        end if

        call show_configurations(config)
    end subroutine


    subroutine show_configurations(config)
        type(mfsolver_config_t), intent(inout) :: config

        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) "Configuration variables for matrix-free diagonalizer"
        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) "  output written to unit        : ", config%output_file_unit
        write(config%output_file_unit, *) "  verbose output                : ", config%verbose

        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) " convergence control"
        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) "  convergence threshold         : ", config%convergence_threshold
        write(config%output_file_unit, *) "  maximum subspace size         : ", config%max_subspace_size
        write(config%output_file_unit, *) "  maximum number of iterations  : ", config%max_iterations
        write(config%output_file_unit, *) "  refresh trial vectors         : ", config%refresh_trial_vectors 
        if (config%refresh_trial_vectors) then
            write(config%output_file_unit, '(A,I4,A)') &
                         "     set to refresh vectors every ",config%trial_vector_refresh_rate," interations"
        end if

        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) " restart/data storage"
        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) "  save subspaces                :", config%save_subspaces
        if (config%save_subspaces) &
        write(config%output_file_unit, *) "  subspaces written to unit     : ", config%subspaces_file_unit
        write(config%output_file_unit, *) "  save sigma vectors            : ", config%save_sigmas
        if (config%save_sigmas) &
        write(config%output_file_unit, *) "  sigma vectors written to unit : ", config%sigmas_file_unit
        write(config%output_file_unit, *) "  save results                  :", config%save_results
        if (config%save_results) &
        write(config%output_file_unit, *) "  results written to unit       : ", config%results_file_unit

        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) " diagonalization characteristics"
        write(config%output_file_unit, *) ""
        write(config%output_file_unit, *) "  solve for right eigenvectors  : ", config%solve_rhs
        write(config%output_file_unit, *) "  solve for left eigenvectors   : ", config%solve_lhs
        write(config%output_file_unit, *) "  symmetric eigenproblem        : ", config%symmetric
        write(config%output_file_unit, *) "  complex mode (real variables) : ", config%complex_mode_with_reals
        write(config%output_file_unit, *) "  root following via overlap    : ", config%overlap_sorting
        if (config%overlap_sorting) then
            write(config%output_file_unit, *) "     overlap following scheme   : ", config%overlap_sorting_scheme
        end if
        write(config%output_file_unit, *) "  energy shift in eigval sorting: ",config%apply_eigenvalues_shift
!       if (config%apply_eigenvalues_shift) then
           write(config%output_file_unit, *) "     shift value, real  part    : ",config%eigenvalues_shift(1)
           write(config%output_file_unit, *) "     shift value, imag. part    : ",config%eigenvalues_shift(2)
!       end if

        write(config%output_file_unit, *) ""

    end subroutine


    subroutine set_safe_minimum(safe_minimum, debug)
        double precision dlamch
        real(kind=8) :: safe_minimum 
        logical :: debug
        safe_minimum = dlamch('S')
        if (debug) print *, "debug, safe minimum (LAPACK) is ",safe_minimum
    end subroutine 

end module

