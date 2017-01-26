module davidson


implicit none

public davidson_setup
public davidson_cleanup
public davidson_driver

!public sort_values_vectors
!public set_safe_minimum
!public davidson_create_results
!public free_davidson_results

interface davidson_driver
    module procedure get_values_vectors_davidson
end interface

public davidson_results

private

! derived types

type davidson_results
    real(kind=8), pointer :: eValuesRe(:) => NULL()
    real(kind=8), pointer :: eValuesIm(:) => NULL()
    real(kind=8), pointer :: eVectorsR(:,:) => NULL()
    real(kind=8), pointer :: eVectorsL(:,:) => NULL()
end type


! module parameters

integer, parameter :: max_label_length = 30
integer, parameter :: nr_max_irreps = 32

real(kind=8), parameter :: zero_r8 = 0.0d0
real(kind=8), parameter :: one_r8  = 1.0d0

! module variables, should be set with a call to davidson_setup

real(kind=8), save :: safe_minimum = 1.0d-36

real(kind=8), save :: davidson_convergence_threshold 

integer, save      :: davidson_max_subspace_size 
integer, save      :: davidson_max_iterations
integer, save      :: davidson_output_file_unit 
integer, save      :: davidson_subspaces_file_unit 
integer, save      :: davidson_sigmas_file_unit 
integer, save      :: davidson_results_file_unit 

logical, save      :: davidson_solve_rhs 
logical, save      :: davidson_solve_lhs 
logical, save      :: davidson_solve_values
logical, save      :: davidson_verbose 
logical, save      :: davidson_symmetric 
logical, save      :: davidson_complex_matrix 
logical, save      :: davidson_save_subspaces 
logical, save      :: davidson_save_sigmas 
logical, save      :: davidson_save_results
logical, save      :: davidson_overlap_sorting

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

!
! create_sigma_vector is the interface used with the diagonalizer. 
! inside it we will handle left and right sigma vector cases, and the sigma_eom module 
! will have to be initialized so that it knows whether to do left or right vectors
!  
   function create_sigma_vector(A, r, old_sigma_vector)
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal
      use sigma_eom, only : intermediates, create_sigma_vector_left, create_sigma_vector_right

      real(kind=8), pointer :: create_sigma_vector(:,:)
      real(kind=8), pointer :: old_sigma_vector(:,:)
      type(intermediates) :: A
      real(kind=8), intent(in) :: r(:,:)

      if (davidson_solve_rhs.and..not.davidson_solve_lhs) then 
!     do right vectors, this is the default
          create_sigma_vector => create_sigma_vector_right(A, r, old_sigma_vector)

      else if (davidson_solve_lhs.and..not.davidson_solve_rhs) then
!     do left vectors vectors
          create_sigma_vector => create_sigma_vector_left(A, r, old_sigma_vector)

      else
! do both at the same time; this will require some thinking, one option is to return
! a matrix of (N,(Lr+Lf)), L(r/l) the number of rhs/lhs vecturs, containing Lr rhs vectors 
! in the first Lr columns and the transpose of sigma_L in the remaining columns. it is easy
! if there's always the same number of columns flr left and right, otherwise there's no 
! simple way to keep track of the dimensions.
          create_sigma_vector => NULL()
      endif 
   end function


! this should never be used to update an instance, only to create new ones. it's better to have results
! for left and right separately if we determine these in two steps. but we must check whether an 
! element of the derived type has been allocated already... perhaps this should then be changed to a
! function returning a results type, and have the davidson driver return that to its caller
    subroutine davidson_create_results(results, N, nroots, values, left, right)
        type(davidson_results), intent(inout) :: results
        integer, intent(in)    :: N, nroots 
        integer                :: nroots_local 
        logical, intent(in)    :: values, left, right

! we don't use complex variables but rather use twice the amount of real variables
! for symmetric matrices we only need one of the vectors, and we choose it to be the right one
! left vectors, if asked, are ignored
        nroots_local = nroots 
        if (davidson_symmetric) then
            if (davidson_complex_matrix) nroots_local = 2*nroots
            if (values) allocate(results%eValuesRe(nroots_local))
            if (right)  allocate(results%eVectorsR(N,nroots_local))
        else
            if (values) then
                allocate(results%eValuesRe(nroots_local))
                allocate(results%eValuesIm(nroots_local))
            end if

            if (davidson_complex_matrix) nroots_local = 2*nroots_local
            if (right) allocate(results%eVectorsR(N,nroots_local))
            if (left)  allocate(results%eVectorsL(N,nroots_local))
        end if
    end subroutine


    subroutine free_davidson_results(results)
        type(davidson_results), intent(inout) :: results

!        if (associated(results%eValuesRe)) deallocate(results%eValuesRe)
        if (allocated(results%eValuesRe)) deallocate(results%eValuesRe)
        if (associated(results%eValuesIm)) deallocate(results%eValuesIm)
        if (associated(results%eVectorsR)) deallocate(results%eVectorsR)
        if (associated(results%eVectorsL)) deallocate(results%eVectorsL)

    end subroutine


    subroutine davidson_setup(convergence_threshold, &
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
                              complex_matrix,    &
                              verbose,           & 
                              overlap_sorting,   &
                              save_results,      &
                              save_sigmas,       &
                              save_subspaces)

        real(kind=8), intent(in), optional :: convergence_threshold
        integer, intent(in), optional :: max_subspace_size 
        integer, intent(in), optional :: max_iterations
        integer, intent(in), optional :: unit_output
        integer, intent(in), optional :: unit_subspaces
        integer, intent(in), optional :: unit_sigmas
        integer, intent(in), optional :: unit_results
        logical, optional :: solve_values
        logical, optional :: solve_right
        logical, optional :: solve_left
        logical, optional :: verbose
        logical, optional :: symmetric
        logical, optional :: overlap_sorting
        logical, optional :: save_sigmas 
        logical, optional :: save_subspaces
        logical, optional :: save_results
        logical, optional :: complex_matrix

        call set_safe_minimum() 

        if (present(convergence_threshold)) then
            davidson_convergence_threshold = convergence_threshold
        else
            davidson_convergence_threshold = 1.0d-5 
        end if

        if (present(max_subspace_size)) then
            davidson_max_subspace_size = max_subspace_size  
        else
            davidson_max_subspace_size = 1024
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

        if (present(overlap_sorting)) then
            davidson_overlap_sorting = overlap_sorting
        else
            davidson_overlap_sorting = .true.
        end if

        if (present(verbose)) then
            davidson_verbose = verbose
        else
            davidson_verbose = .false.
        end if

        if (present(symmetric)) then
            davidson_symmetric = symmetric
        else
            davidson_symmetric = .false. 
        end if

        if (present(complex_matrix)) then
            davidson_complex_matrix = complex_matrix 
        else
            davidson_complex_matrix = .false. 
        end if

        if (present(unit_output)) then
            davidson_output_file_unit = unit_output
        end if

        if (present(save_subspaces)) then
           davidson_save_subspaces = save_subspaces
        else
           davidson_save_subspaces = .true. 
        end if

        if (present(save_sigmas)) then
           davidson_save_sigmas = save_sigmas
        else
           davidson_save_sigmas = .true. 
        end if

        if (present(save_results)) then
           davidson_save_results = save_results
        else
           davidson_save_results = .true. 
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
        write(davidson_output_file_unit, *) "  complex matrix                : ", davidson_complex_matrix
        write(davidson_output_file_unit, *) "  root following via overlap    : ", davidson_overlap_sorting

        write(davidson_output_file_unit, *) ""

    end subroutine


    subroutine set_safe_minimum()
        double precision dlamch
        safe_minimum = dlamch('S')
        print *, "debug, safe minimum (LAPACK) is ",safe_minimum
    end subroutine 


    subroutine davidson_cleanup(results)
        type(davidson_results), intent(inout) :: results
        call free_davidson_results(results)
    end subroutine


    function new_subspace_vectors_R(new_L, N, old_subspace_vectors_R, verbose)
!   new_L is the dimension of the subspace to be created 
!   N is the dimension of the original problem
!   new_subspace_vectors_R  is the subspace we create or update
!   old_subspace_vectors_R (if present) is the old subspace
        real(kind=8), pointer :: new_subspace_vectors_R(:,:)

        integer, intent (in) :: new_L, N
        real(kind=8), intent(inout), optional, pointer :: old_subspace_vectors_R(:,:)
        logical, optional, intent(in) :: verbose

        real(kind=8), pointer :: old_subspace_vectors_local(:,:)
        integer :: i, j, old_L
        logical :: verbose_local

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        allocate(new_subspace_vectors_R(N, new_L))
        new_subspace_vectors_R = zero_r8

        if (.not.present(old_subspace_vectors_R)) then
            do i = 1, N
                do j = 1, new_L
                    if (i .eq. j) then
                        new_subspace_vectors_R(i,j) = one_r8
                    end if
                end do
            end do
        else
            old_subspace_vectors_local => old_subspace_vectors_R
            old_L = size(old_subspace_vectors_local,2)
            do i = 1, old_L 
                new_subspace_vectors_R(:,i) = old_subspace_vectors_local(:,i)
            end do

            nullify(old_subspace_vectors_local)
            deallocate(old_subspace_vectors_R)
        end if
    end function

    function updated_subspace_vectors_R(old_vectors_R, delta, verbose)
        real(kind=8), pointer:: updated_subspace_vectors_R(:,:)

        real(kind=8), intent(inout), pointer :: old_vectors_R(:,:)
        real(kind=8), intent(in)      :: delta(:,:)
        logical, optional, intent(in) :: verbose

        integer :: vectors_L, vectors_N
        integer :: vectors_Ld, vectors_Nd
        integer :: new_L
        integer :: i
        logical :: verbose_local
        real(kind=8), allocatable :: c(:,:)

        vectors_N = size(old_vectors_R,1)
        vectors_L = size(old_vectors_R,2)

        vectors_Nd= size(delta,1)
        vectors_Ld= size(delta,2)

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        new_L = vectors_L + vectors_Ld
        updated_subspace_vectors_R => new_subspace_vectors_R(new_L,vectors_N,old_subspace_vectors_R=old_vectors_R)
        do i = 1, vectors_Ld  
           updated_subspace_vectors_R(:, vectors_L + i) = delta(:,i)
        end do

        if (verbose_local) &
            call print_vectors_R(updated_subspace_vectors_R, label="new subspace vectors before GS")

        call orthonormalize_subspace(updated_subspace_vectors_R, verbose=verbose_local)

        if (associated(old_vectors_R)) deallocate(old_vectors_R)
    end function

    subroutine orthonormalize_subspace(b, method, verbose)
        real(kind=8), intent(inout)   :: b(:,:)
        logical, optional, intent(in) :: verbose
        integer, optional, intent(in) :: method

        integer :: method_local
        logical :: verbose_local

        if (present(method)) then
            method_local = method 
        else
            method_local = 1 ! 1 = GS, 2 = QR
        end if

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        if (method_local .eq. 1) then
            call gram_schmidt_orthogonalization(b, verbose=verbose_local)
        else if (method_local .eq. 2) then
!           cal QR(b, verbose=verbose_local)
        end if


    end subroutine

    subroutine gram_schmidt_projection_R(projection, u, v)
        real(kind=8), intent(inout) :: projection(:)
        real(kind=8), intent(in) :: u(:)
        real(kind=8), intent(in) :: v(:)
        real(kind=8) :: numerator, denominator, factor
        integer :: N
        double precision ddot

        N = size(u,1)
        denominator = ddot(N, u, 1, u, 1) 
        numerator   = ddot(N, v, 1, u, 1) 

        projection = zero_r8
        if (denominator .gt. safe_minimum) then
            factor = numerator/denominator
        else
            factor = zero_r8
        end if
        if (.false.) then
            print *, "GS projection, factor: ",factor
        end if
        projection = factor * u
    end subroutine 

    subroutine gram_schmidt_orthogonalization(vectors, verbose)
        real(kind=8), intent(inout) :: vectors(:,:)
        logical, optional, intent(in) :: verbose

        real(kind=8), allocatable   :: projection(:)
        integer :: vectors_L, vectors_N
        integer :: i, j, k
        logical :: verbose_local
        real(kind=8) :: norm

        vectors_L = size(vectors,2)
        vectors_N = size(vectors,1)

        allocate(projection(vectors_N))

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        do i = 1, vectors_L
            norm = vector_norm(vectors(:,i))
            if (verbose_local) then
                print *, "debug, Gram-Schmidt orthogonalization"
                print *, "  renormalizing vector v(",i,") : (norm = ",norm,")"
                do k = 1, vectors_N
                    print "(4x,1e24.16)", vectors(k,i)
                end do
                print *, ""
            end if 
            if (norm .ge. safe_minimum ) then
                vectors(:,i) = vectors(:,i) / norm
            else
                stop ('division by zero in orthogonalization') 
            end if
            if (verbose_local) then
                print *, "renormalized  vector v(",i,") : " 
                do k = 1, vectors_N
                    print "(4x,1e24.16)", vectors(k,i)
                end do
                print *, ""
            end if
            do j = i+1, vectors_L
                call gram_schmidt_projection_R(projection, vectors(:,i),vectors(:,j))
                if (verbose_local) then
                    print *," projection of v(",i,") onto v(",j,"):"
                    do k = 1, vectors_N
                        print "(4x,1e24.16,4x,1e24.16)", vectors(k,j), projection(k)
                    end do
                    print *, ""
                end if
                vectors(:,j) = vectors(:,j) - projection(:)
            end do
            if (verbose_local) then
                call print_vectors_R(vectors, label="all vectors                   ")
            end if
        end do

        deallocate(projection)
    end subroutine 

    subroutine print_vectors_R(vectors, label, nvectors)
        real(kind=8), intent(in) :: vectors(:,:)
        integer, intent(in), optional :: nvectors
        character(len=max_label_length), intent(in), optional :: label 

        integer :: L, N, i, j
        character(len=max_label_length) :: label_local

        N = size(vectors,1)
        L = size(vectors,2)
        if (present(nvectors)) L = nvectors

        if (present(label)) then
            label_local = label
            print *, ""
            print "(a)", label_local(1:max_label_length)
            print *, ""
        end if
        do i = 1, N
            print "(2x,i6,$)", i
            do j = 1, L
                print "(4x,1e24.16,$)", vectors(i,j)
            end do
            print *, ""
        end do
    end subroutine


    subroutine print_values_vectors_R(output_file_unit, eValuesRe, eValuesIm, eVectorsR, eVectorsL, left, right, nroots)
        real(kind=8), intent(in) :: eVectorsL(:,:)
        real(kind=8), intent(in) :: eVectorsR(:,:)
        real(kind=8), intent(in) :: eValuesRe(:), eValuesIm(:)
        integer, optional, intent(in) :: nroots
        integer, intent(in) :: output_file_unit

        logical :: left, right
        integer :: vectors_N
        integer :: nroots_local
        integer :: i
        character(len=max_label_length) :: label

        vectors_N = size(eValuesRe,1)
        nroots_local = vectors_N
        if (present(nroots)) nroots_local = nroots

        label = "Eigenvalues"
        do i = 1, nroots_local 
            if (eValuesIm(i) .eq. zero_r8) then
                write(output_file_unit,"(4x,a,i6,a,1e24.16)") "Eigenvalue ",i," : ",eValuesRe(i)
            else
                write(output_file_unit,"(4x,a,i6,a,1e24.16,a,1e24.16,a)") &
                      "Eigenvalue ",i," : ",eValuesRe(i)," + ",eValuesIm(i), " i "
            end if
        end do
        if (left) then
            label = "Left Eigenvectors            "
            call print_vectors_R(eVectorsL, label=label) 
        end if
        if (right) then
            label = "Right Eigenvectors           "
            call print_vectors_R(eVectorsR, label=label) 
        end if

    end subroutine 

    function vector_norm(v)
        real(kind=8) :: vector_norm
        real(kind=8), intent(in) :: v(:)
        integer :: vector_N, i
        double precision dnrm2

        vector_N = size(v,1)
        vector_norm = dnrm2(vector_N, v, 1)
    end function

    function correction_to_vectors(output_file_unit, nroots, e_small_M, c_small_M, sigma_vector, b, preconditioner, correct_all)
        real(kind=8), pointer     :: correction_to_vectors(:,:)
        integer, intent(in)       :: nroots, output_file_unit
        real(kind=8), intent(in)  :: e_small_M(:)
        real(kind=8), intent(in)  :: c_small_M(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        real(kind=8), pointer :: preconditioner(:,:)
        logical, optional, intent(in) :: correct_all

        integer :: vectors_L, vectors_N
        integer :: k, unconverged_roots
        real(kind=8) :: norm_r = zero_r8
        real(kind=8), allocatable :: r(:,:)
        logical, allocatable :: vectors_to_update(:)
        logical :: correct_all_local

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (present(correct_all)) then
            correct_all_local = correct_all
        else
            correct_all_local = .false. 
        end if

        allocate(vectors_to_update(nroots))
        vectors_to_update = .true.

        allocate(r(vectors_N,nroots))
        call form_r(r, nroots, e_small_M, c_small_M, sigma_vector, b)

        unconverged_roots = 0
        do k = 1, nroots
            norm_r = vector_norm(r(:,k))
            if (norm_r .lt. davidson_convergence_threshold) then 
               write(output_file_unit, "(4x,a,i6,a,1e24.16,a,1e24.16,a)") "Root ",k,"  ||q_k,M|| ", norm_r,  &
                                                      " threshold ", davidson_convergence_threshold, &
                                                      " <<< converged ! "
               vectors_to_update(k) = .false.
            else
               write(output_file_unit, "(4x,a,i6,a,1e24.16,a,1e24.16)")  "Root ",k,"  ||q_k,M|| ", norm_r, &
                                                     " threshold ", davidson_convergence_threshold
               vectors_to_update(k) = .true.
               unconverged_roots = unconverged_roots + 1
            end if
        end do
        print *, ""

        if (any(vectors_to_update)) then
            if (correct_all_local) then
                unconverged_roots = nroots
                vectors_to_update = .true.
            end if
            allocate(correction_to_vectors(vectors_N,unconverged_roots))
            correction_to_vectors = zero_r8
            call form_delta(correction_to_vectors, vectors_to_update, r, preconditioner)
        else
            correction_to_vectors => NULL()
        end if

        deallocate(vectors_to_update)
        deallocate(r)
    end function

    function form_preconditioner(M_diag, eigenvectors)
        real(kind=8), pointer     :: form_preconditioner(:,:)
        real(kind=8), intent(in)  :: eigenvectors(:)
        real(kind=8), intent(in)  :: M_diag(:)

        integer :: i, k, l
        integer :: vectors_N
        integer :: nroots

        nroots    = size(eigenvectors,1)
        vectors_N = size(M_diag,1)

        allocate(form_preconditioner(vectors_N,nroots))
        form_preconditioner = 0.0d0
        do k = 1, nroots
            form_preconditioner(:,k) = eigenvectors(k) - M_diag(:)
            do l = 1, vectors_N
                if (dabs(form_preconditioner(l,k)) .gt. safe_minimum) then
                    form_preconditioner(l,k) = one_r8 / form_preconditioner(l,k)
                else
!                    stop ('detected division by zero')
                    form_preconditioner(l,k) = one_r8 
                end if
            end do
        end do
    end function

    subroutine form_delta(delta, vectors_to_update, r, preconditioner)
        real(kind=8), intent(inout) :: delta(:,:)
        logical, intent(in)         :: vectors_to_update(:)
        real(kind=8), intent(in)    :: preconditioner(:,:)
        real(kind=8), intent(in)    :: r(:,:)

        integer :: nroots
        integer :: i, k

        nroots    = size(vectors_to_update,1)
 
        i = 1
        do k = 1, nroots
            if (vectors_to_update(k)) then
               delta(:,i) = preconditioner(:,k)*r(:,k)
               delta(:,i) = delta(:,i) / vector_norm(delta(:,i))
               i = i + 1
           end if
        end do
    end subroutine 

    subroutine form_r(r, nroots,  e_small_M, c_small_M, sigma_vector,  b)
        real(kind=8), intent(inout) :: r(:,:)
        integer, intent(in)       :: nroots
        real(kind=8), intent(in)  :: e_small_M(:)
        real(kind=8), intent(in)  :: c_small_M(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        integer :: k
        real(kind=8) :: alpha_lapack = one_r8
        real(kind=8) :: beta_lapack  = zero_r8
        real(kind=8), allocatable :: residue_a(:,:)
        real(kind=8), allocatable :: residue_b(:,:)

        vectors_L = size(b,2)
        vectors_N = size(b,1)

        allocate(residue_a(vectors_N,vectors_L))
        allocate(residue_b(vectors_N,vectors_L))

        residue_a = zero_r8
        residue_b = zero_r8

        call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                   alpha_lapack, sigma_vector, vectors_N, c_small_M, vectors_L, &
                   beta_lapack, residue_a, vectors_N)
        call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                   alpha_lapack, b, vectors_N, c_small_M, vectors_L,            &
                   beta_lapack, residue_b, vectors_N)
        if (.false.) then
            print *, "debug in form_r, residue_a :",residue_a
            print *, "debug in fomr_r, residue_b :",residue_b
        end if
        do k = 1, nroots
            r(:,k) = residue_a(:,k) - e_small_M(k)*residue_b(:,k)
        end do

        deallocate(residue_a)
        deallocate(residue_b)
    end subroutine

    function create_reduced_problem(b, sigma_vector)
        real(kind=8), pointer     :: create_reduced_problem(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        real(kind=8) :: alpha_lapack = 1.0d0
        real(kind=8) :: beta_lapack = 1.0d0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_reduced_problem(vectors_L,vectors_L))
        create_reduced_problem = zero_r8
        call dgemm('T', 'N', vectors_L, vectors_L, vectors_N, alpha_lapack, b, vectors_N, &
                   sigma_vector, vectors_N, beta_lapack, create_reduced_problem, vectors_L)

    end function

    subroutine solve_values_vectors(M, eValuesRe, eValuesIm, eVectorsR, eVectorsL, left)
        real(kind=8), intent(inout), pointer :: M(:,:)
        real(kind=8), intent(inout), pointer :: eVectorsL(:,:)
        real(kind=8), intent(inout), pointer :: eVectorsR(:,:)
        real(kind=8), intent(inout), pointer :: eValuesRe(:), eValuesIm(:)
        logical, intent(in) :: left 

        integer :: vectors_L, i

        integer :: lwork, info
        character(len=1) :: do_left = 'N', do_right = 'V'
        real(kind=8), allocatable :: work(:)

        vectors_L = size(M,1)

        if (associated(eVectorsL)) deallocate(eVectorsL)
        if (associated(eVectorsR)) deallocate(eVectorsR)
        if (associated(eValuesRe)) deallocate(eValuesRe)
        if (associated(eValuesIm)) deallocate(eValuesIm)

        if (left) then
            do_left = 'V'
            allocate(eVectorsL(vectors_L,vectors_L))
            eVectorsL = zero_r8
        end if
        allocate(eVectorsR(vectors_L,vectors_L))
        allocate(eValuesRe(vectors_L))
        allocate(eValuesIm(vectors_L))

        eVectorsR = zero_r8
        eValuesRe = zero_r8
        eValuesIm = zero_r8

! first query optimal lwork value
        lwork = -1
        allocate(work(4*vectors_L))
        call dgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                    eValuesRe, eValuesIm, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, INFO)
        lwork = work(1)
        deallocate(work)
        allocate(work(lwork))
        call dgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                    eValuesRe, eValuesIm, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, INFO)
        deallocate(work)
        if (info .gt. 0) then
            print *, 'in small_solve_eigenvalue: not all eigenvalues converged' 
            print *, 'converged'
            do i = info+1, vectors_L 
                print *, " real part : ",eValuesRe(i), " imaginary part : ",eValuesIm(i)
            end do
        end if

    end subroutine
 
    function get_diagonal_elements(M)
        real(kind=8), intent(in)    :: M(:,:)
        real(kind=8), pointer :: get_diagonal_elements(:)

        integer :: vectors_N
        integer :: i

        vectors_N = size(M,1)
        allocate(get_diagonal_elements(vectors_N))

        do i = 1, vectors_N
            get_diagonal_elements(i) = M(i,i)
        end do
    end function

    subroutine get_values_vectors_davidson(nroots, A, results, b_guess)
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal
        use sigma_eom, only : intermediates, get_matrix_size, create_diagonal
        type(intermediates), intent(in)   :: A
! back to general stuff
        type(davidson_results), intent(inout) :: results
        real(kind=8), intent(in), pointer :: b_guess(:,:)
        integer, intent(in) :: nroots
!
! nullifying all pointers at the start to give them a state, see discussion on http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html
! "[...] When a pointer is declared its status is undefined, and cannot be safely queried with the associated intrinsic. [...]"
! 
        real(kind=8), pointer :: M_diag(:)           => NULL()
        real(kind=8), pointer :: delta(:,:)          => NULL()
        real(kind=8), pointer :: preconditioner(:,:) => NULL()
        real(kind=8), pointer :: sigma_vector(:,:)   => NULL()
        real(kind=8), pointer :: b(:,:)              => NULL()
        real(kind=8), pointer :: M_small(:,:)        => NULL()
        real(kind=8), pointer :: c_small_r(:,:)      => NULL()
        real(kind=8), pointer :: c_small_l(:,:)      => NULL()
        real(kind=8), pointer :: e_small_r(:)        => NULL()
        real(kind=8), pointer :: e_small_i(:)        => NULL()

        real(kind=8), pointer :: c_small_r_old(:,:)  => NULL()
        real(kind=8), pointer :: c_small_l_old(:,:)  => NULL()

        integer :: vectors_N, L, min_subspace_size, max_subspace_size, niter
        integer :: i, j

        vectors_N = get_matrix_size(A)
!
! davidson_create_results() will allocate storage for results from davidson diagonalization when called a first time. 
! 
        call davidson_create_results(results,vectors_N,nroots,&
                                     davidson_solve_values,   &
                                     davidson_solve_lhs,      &
                                     davidson_solve_rhs)

! todo: replace lines below with a call to a function which also performs the allocation, like
!        M_diag => create_diagonal(A)
        allocate(M_diag(size(A%H_II))) 
        M_diag = A%H_II 

        min_subspace_size = nroots 
        max_subspace_size = max(nroots, davidson_max_subspace_size)
        max_subspace_size = min(vectors_N, max_subspace_size)

        write(davidson_output_file_unit,"(4x,a,i6,a)") "Requested ",nroots," eigenvales"
        if ( nroots .gt. vectors_N ) then 
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested larger than matrix order !"
           return 
        endif

        b => b_guess

        if (davidson_verbose) then
            write(davidson_output_file_unit,"(4x,a,i4,a,i8)") &
                  "Start by mapping M onto a subspace of size ",min_subspace_size, &
                  ", will extend it up to ",max_subspace_size
            call print_vectors_R(b, label="Initial subspace vectors      ")
        end if

       if (davidson_solve_rhs.and.davidson_overlap_sorting) then
           allocate(c_small_r_old(nroots,nroots))
           c_small_r_old = 0.0d0
           do i = 1, nroots
               c_small_r_old(i,i) = 1.0d0
           enddo
       endif
       if (davidson_solve_lhs.and.davidson_overlap_sorting) then
           allocate(c_small_l_old(nroots,nroots))
           c_small_l_old = 0.0d0
           do i = 1, nroots
               c_small_l_old(i,i) = 1.0d0
           enddo
       endif
 

        do niter = 1, davidson_max_iterations
            write(davidson_output_file_unit,*) ""
            write(davidson_output_file_unit,*) "Iteration ", niter
            write(davidson_output_file_unit,*) ""

            if (associated(M_small)) deallocate(M_small)

            sigma_vector => create_sigma_vector(A, b, sigma_vector)


            if (davidson_solve_rhs) M_small => create_reduced_problem(b, sigma_vector)
            if (davidson_solve_lhs) M_small => create_reduced_problem(sigma_vector, b)

            call calculate_eigenvalues_and_eigenvectors(davidson_output_file_unit,    &
                                                        M_small,                      &
                                                        e_small_r, e_small_i,         &
                                                        c_small_r, c_small_l,         &
                                                        c_small_r_old, c_small_l_old, &
                                                        left=davidson_solve_lhs,      &
                                                        right=davidson_solve_rhs,     &
                                                        verbose=davidson_verbose)

! store old c_small vectors, that is, subspace vectors for the overlap selection in the next iterative step   
            if (davidson_solve_rhs.and.davidson_overlap_sorting) c_small_r_old = c_small_r(1:nroots,1:nroots)
            if (davidson_solve_lhs.and.davidson_overlap_sorting) c_small_l_old = c_small_l(1:nroots,1:nroots)

            if (associated(preconditioner)) deallocate(preconditioner)
            if (associated(delta)) deallocate(delta)

            preconditioner =>  form_preconditioner(M_diag, e_small_r)
           if (davidson_solve_rhs) delta => correction_to_vectors(davidson_output_file_unit, nroots,  &
                                           e_small_r, c_small_r, sigma_vector, & 
                                           b, preconditioner)

           if (davidson_solve_lhs) delta => correction_to_vectors(davidson_output_file_unit, nroots,  &
                                           e_small_r, c_small_l, sigma_vector, & 
                                           b, preconditioner)

            if (.not.associated(delta)) then
                write(davidson_output_file_unit,*) "All roots converged."
                 if (associated(preconditioner)) deallocate(preconditioner)
                exit
            end if
            b => updated_subspace_vectors_R(b, delta, verbose=davidson_verbose)

            if (davidson_verbose) then
                call print_vectors_R(delta, label="delta                         ")
                call print_vectors_R(b, label="Updated subspace vectors      ")
            end if
  
            L = size(b,2)
            if (L .eq. max_subspace_size) then
               if (associated(preconditioner)) deallocate(preconditioner)
               if (associated(delta)) deallocate(delta)
                write(davidson_output_file_unit, *) "Reached maximum subspace size, stopping."
                exit
            end if
        end do

        if (davidson_solve_rhs.and.davidson_overlap_sorting) deallocate(c_small_r_old)
        if (davidson_solve_lhs.and.davidson_overlap_sorting) deallocate(c_small_l_old)

        write(davidson_output_file_unit, *) "Final size of reduced subspace : ",L

        call form_eigenvalues_full_problem(results%eValuesRe, results%eValuesIm, nroots, e_small_r, e_small_i) 

        if (davidson_solve_lhs) &
            call form_eigenvector_full_problem(results%eVectorsL, nroots,b,c_small_l,davidson_solve_lhs)
        if (davidson_solve_rhs) &
            call form_eigenvector_full_problem(results%eVectorsR, nroots,b,c_small_r,davidson_solve_rhs)

        if(associated(b)) deallocate(b)
        deallocate(M_small)
        deallocate(sigma_vector)
        deallocate(M_diag)
        deallocate(e_small_r)
        deallocate(e_small_i)
        deallocate(c_small_r)
        if (associated(c_small_l)) deallocate(c_small_l)

    end subroutine

    subroutine calculate_eigenvalues_and_eigenvectors(output_file_unit, M, eValuesRe, eValuesIm, eVectorsR, eVectorsL, &
                                eVectorsR_old, eVectorsL_old, left, right, sorting_order, nroots, verbose)
        real(kind=8), intent(inout), pointer   :: M(:,:)
        real(kind=8), intent(inout), pointer   :: eValuesRe(:), eValuesIm(:)
        real(kind=8), intent(inout), pointer   :: eVectorsR(:,:), eVectorsL(:,:)
        real(kind=8), intent(in)   :: eVectorsR_old(:,:), eVectorsL_old(:,:) !old vectors
        logical, optional, intent(in) :: left, right, verbose
        integer, optional, intent(in) :: sorting_order, nroots
        integer, intent(in) :: output_file_unit 

        logical :: left_local, right_local, verbose_local
        integer :: sorting_order_local, nroots_local, nroot
        real(kind=8), allocatable :: overlap(:,:)
        real(kind=8) :: alpha_lapack, beta_lapack

        alpha_lapack = 1.0d0
        beta_lapack  = 1.0d0

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

        if (present(nroots)) then
            nroots_local = nroots
        else
            nroots_local = size(M,1) 
        end if

        call solve_values_vectors(M, eValuesRe, eValuesIm, eVectorsR, eVectorsL, left_local)

        if (davidson_overlap_sorting) then
! calculate overlap between old and new vectors. Overlap(nroot,mroot) = Old(nroot,nroot)*New(1:nroot,mroot)   
! Max(Overlap(iroot,mroot)) for all iroot = 1,nroot 
            nroot = size(eVectorsR_old,1)
            allocate(overlap(nroot,nroots_local))
            overlap = 0.0d0 

            if(right) then    
                 call dgemm('N', 'N', nroot, nroots_local, nroot, alpha_lapack, eVectorsR_old, &
                            nroot, eVectorsR(1:nroot,:), nroot, beta_lapack, overlap, nroot)
            else
                 call dgemm('N', 'N', nroot, nroots_local, nroot, alpha_lapack, eVectorsL_old, &
                            nroot, eVectorsL(1:nroot,:), nroot, beta_lapack, overlap, nroot)
            endif

            call sort_by_overlap(-1, eValuesRe, eValuesIm, overlap, eVectorsR, eVectorsL, &
                                 right_local, left_local) 
            deallocate(overlap)
        else
            call sort_values_vectors(sorting_order_local, eValuesRe, eValuesIm, eVectorsR, eVectorsL, &
                                     right_local, left_local)
        end if

        if (davidson_verbose) &
        call print_values_vectors_R(davidson_output_file_unit,eValuesRe, eValuesIm, eVectorsR, eVectorsL, &
                                    left_local, right_local, nroots=nroots_local)

    end subroutine

    subroutine form_eigenvalues_full_problem(eValuesRe, eValuesIm, nroots, e_small_r, e_small_i)
        real(kind=8), intent(inout) :: eValuesRe(:), eValuesIm(:)
        real(kind=8), intent(in)    :: e_small_r(:), e_small_i(:)
        integer, intent(in) :: nroots

        integer :: i

        do i = 1, nroots
            eValuesRe(i) = e_small_r(i)
            eValuesIm(i) = e_small_i(i)
        end do
    end subroutine

    subroutine form_eigenvector_full_problem(eVectorsFull,nroots,b,eVectorsSmall,left_vector) 
        real(kind=8), intent(inout) :: eVectorsFull(:,:)
        real(kind=8), intent(in)    :: eVectorsSmall(:,:)
        real(kind=8), intent(in)    :: b(:,:)
        logical, intent(in) :: left_vector
        integer, intent(in) :: nroots

        integer :: i
        real(kind=8), allocatable :: bc(:,:)
        integer :: vectors_L, vectors_N
        real(kind=8) :: alpha_lapack = 1.0d0
        real(kind=8) :: beta_lapack = 1.0d0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (left_vector) then
            allocate(bc(vectors_L,vectors_N))
            bc = zero_r8
!           call print_vectors_R(eVectorsSmall, label="L eigenvectors small problem  ")
            call dgemm('T', 'T', vectors_L, vectors_N, vectors_L, alpha_lapack, eVectorsSmall, &
                       vectors_L, b, vectors_N, beta_lapack, bc, vectors_L)
!           call print_vectors_R(bc, label="bc                            ")
            do i = 1, nroots
                eVectorsFull(:,i) = bc(i,:) 
            end do

        else
            allocate(bc(vectors_N,vectors_L))
            bc = zero_r8
!           call print_vectors_R(eVectorsSmall, label="R eigenvectors small problem  ")
            call dgemm('N', 'N', vectors_N, vectors_L, vectors_L, alpha_lapack, b, &
                       vectors_N, eVectorsSmall, vectors_L, beta_lapack, bc, vectors_N)

!           call print_vectors_R(bc, label="bc                            ")
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

        real(kind=8) :: alpha_lapack = 1.0d0
        real(kind=8) :: beta_lapack = 0.0d0

        vectors_N = size(b,1)
        vectors_L = size(b,2)
        vectors_L_old = size(old_sigma_vector,2)

        allocate(update_sigma_vector(vectors_N,vectors_L))
        do i = 1, vectors_L_old 
           update_sigma_vector(:, i) = old_sigma_vector(:, i)  
        end do
        deallocate(old_sigma_vector)
!ashee        nullify(old_sigma_vector)

        new_L = vectors_L - vectors_L_old
        call dgemm('N', 'N', vectors_N, new_L, vectors_N, alpha_lapack, M, vectors_N, &
       &           b(:, vectors_L_old:vectors_L), vectors_N, beta_lapack,             &
       &           update_sigma_vector(:,vectors_L_old:vectors_L), vectors_N)
    end function


    function create_sigma_vector_full_matrix(M, b)
        real(kind=8), pointer :: create_sigma_vector_full_matrix(:,:)
        real(kind=8), intent(in)  :: M(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        real(kind=8) :: alpha_lapack = 1.0
        real(kind=8) :: beta_lapack = 0.0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_sigma_vector_full_matrix(vectors_N,vectors_L))
        call dgemm('N', 'N', vectors_N, vectors_L, vectors_N, alpha_lapack, M, vectors_N, &
                   b, vectors_N, beta_lapack, create_sigma_vector, vectors_N)

    end function

    subroutine sort_values_vectors(order, e_small_r, e_small_i, c_small_r, c_small_l, right, left) 
        integer, intent(in) :: order ! 1: ascending, -1: descending
        logical, intent(in) :: right ! if true, sort also the right eigenvectors
        logical, intent(in) :: left  ! if true, sort also the left eigenvectors
        real(kind=8), intent(inout) :: e_small_r(:)
        real(kind=8), intent(inout) :: e_small_i(:)
        real(kind=8), intent(inout) :: c_small_r(:,:)
        real(kind=8), intent(inout) :: c_small_l(:,:)

        real(kind=8), allocatable :: temp_vector(:,:)
        real(kind=8) :: temp
        integer, allocatable :: indices_sorted_array(:)
        integer :: vectors_L, i_sorted, i

        vectors_L = size(e_small_r,1)
        allocate(indices_sorted_array(vectors_L))
        allocate(temp_vector(vectors_L,vectors_L))
        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
        end do

        call quicksort_vector(e_small_r, indices_sorted_array, 1, vectors_L)
        ! we sort the indices that would sort the eigenvalues array, and apply the change to all other arrays
        do i = 1, vectors_L
         i_sorted = indices_sorted_array(i)

         temp = e_small_i(i)
         e_small_i(i) = e_small_i(i_sorted)
         e_small_i(i_sorted) = temp
        end do

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

    end subroutine


    subroutine sort_by_overlap(order, e_small_r, e_small_i, overlap, c_small_r, c_small_l, right, left) 
        integer, intent(in) :: order ! 1: ascending, -1: descending
        logical, intent(in) :: right ! if true, sort also the right eigenvectors
        logical, intent(in) :: left  ! if true, sort also the left eigenvectors
        real(kind=8), intent(inout) :: overlap(:,:)
        real(kind=8), intent(inout) :: c_small_r(:,:)
        real(kind=8), intent(inout) :: c_small_l(:,:)
        real(kind=8), intent(inout) :: e_small_r(:)
        real(kind=8), intent(inout) :: e_small_i(:)

        real(kind=8), allocatable :: temp_vector(:,:)
        real(kind=8),allocatable :: temp(:)
        integer, allocatable :: indices_sorted_array(:)
        integer :: vectors_L, i_sorted, i, max_overlap(1), iroot, nroot

        vectors_L = size(overlap,2)
        nroot     = size(overlap,1)
        allocate(indices_sorted_array(vectors_L))
        allocate(temp_vector(vectors_L,vectors_L))
        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
        end do

        do iroot = 1, nroot

          max_overlap = maxloc(abs(overlap(iroot:iroot,:)),2)
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
             temp(i) = e_small_r(i_sorted)
            end do
            e_small_r = temp

        deallocate(temp)
        deallocate(indices_sorted_array)

    end subroutine


    recursive subroutine quicksort_vector(v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:)
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, t
        integer :: i, j, ti

        x = v( (first + last) / 2 )
        i = first
        j = last
        do
           do while (v(i) < x)
              i = i + 1
           end do
           do while (x < v(j))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(i)
           v(i) = v(j)
           v(j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) call quicksort_vector(v, indices_sorted, first, i-1)
        if (j+1 < last)  call quicksort_vector(v, indices_sorted, j+1, last)
    end subroutine 


! I/O for davidson

    subroutine io_davidson_setup(N, L, nroot, complex_mode, left, right, symmetry, node)
        integer, intent(in) :: N, L, nroot, node, symmetry
        logical, intent(in) :: left, right, complex_mode
        integer :: reclen_real, reclen_vector, factor
        real(kind=8) :: dummy

        factor = 1
        if (complex_mode) factor = 2

        inquire (iolength = reclen_real) dummy

        left_sigma_file_name     = 'LSIGMAV.sym'//trim(str(symmetry))//"."//trim(str(node))
        left_vectors_file_name   = 'LEVECTO.sym'//trim(str(symmetry))//"."//trim(str(node))
        left_subspace_file_name  = 'LSUBSPA.sym'//trim(str(symmetry))//"."//trim(str(node))

        right_sigma_file_name    = 'RSIGMAV.sym'//trim(str(symmetry))//"."//trim(str(node))
        right_vectors_file_name  = 'REVECTO.sym'//trim(str(symmetry))//"."//trim(str(node))
        right_subspace_file_name = 'RSUBSPA.sym'//trim(str(symmetry))//"."//trim(str(node))

        if (right) then
            reclen_vector = factor*N*L*reclen_real
            open(unit=right_sigma_file_unit,file=right_sigma_file_name,&
                 form="unformatted", access="direct", recl=reclen_vector)
            open(unit=right_subspace_file_unit, file=right_subspace_file_name, &
                 form="unformatted", access="direct", recl=reclen_vector)

            reclen_vector = factor*N*nroot*reclen_real + 2 
            open(unit=right_vectors_file_unit, file=right_vectors_file_name, &
                 form="unformatted", access="direct", recl=reclen_vector)
        end if

        if (left) then
            reclen_vector = factor*N*L*reclen_real
            open(unit=left_sigma_file_unit, file=left_sigma_file_name, &
                 form='unformatted', access='direct', recl=reclen_vector)
            open(unit=left_subspace_file_unit, file=left_subspace_file_name, &
                 form='unformatted', access='direct', recl=reclen_vector)

            reclen_vector = factor*N*nroot*reclen_real + 2 
            open(unit=left_vectors_file_unit, file=left_vectors_file_name, &
                 form='unformatted', access='direct', recl=reclen_vector)
        end if

    end subroutine

    subroutine io_davidson_finalize(left,right)
        logical, intent(in) :: left, right

        if (right) then
            close(unit=right_sigma_file_unit)
            close(unit=right_vectors_file_unit)
            close(unit=right_subspace_file_unit)
        end if

        if (left) then
            close(unit=left_sigma_file_unit)
            close(unit=left_vectors_file_unit)
            close(unit=left_subspace_file_unit)
        end if

    end subroutine

    subroutine write_davidson_values_vectors(file_unit,valuesRe,valuesIm,vectors,column_start,column_end)
        real(kind=8), intent(in) :: vectors(:,:), valuesRe(:), valuesIm(:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            write(file_unit, rec=i) valuesRe(i), valuesIm(i), vectors(:,i)
        end do

    end subroutine

    subroutine read_davidson_values_vectors(file_unit,valuesRe,valuesIm,vectors,column_start,column_end)
        real(kind=8), intent(inout) :: vectors(:,:), valuesRe(:), valuesIm(:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            read(file_unit, rec=i) valuesRe(i), valuesIm(i), vectors(:,i)
        end do

    end subroutine

    subroutine write_davidson_data_1d(file_unit,data1d)
        real(kind=8), intent(inout) :: data1d(:)
        integer, intent(in) :: file_unit

        write(file_unit, rec=1) data1d(:)
    end subroutine

    subroutine read_davidson_data_1d(file_unit,data1d)
        real(kind=8), intent(inout) :: data1d(:)
        integer, intent(in) :: file_unit

        read(file_unit, rec=1) data1d(:)
    end subroutine

    subroutine write_davidson_data_2d(file_unit,data2d,column_start,column_end)
        real(kind=8), intent(in) :: data2d(:,:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            write(file_unit, rec=i) data2d(:,i)
        end do
    end subroutine

    subroutine read_davidson_data_2d(file_unit,data2d,column_start,column_end)
        real(kind=8), intent(inout) :: data2d(:,:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            read(file_unit, rec=i) data2d(:,i)
        end do
    end subroutine

! taken from http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    character(len=20) function str(k)
!       "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end module

