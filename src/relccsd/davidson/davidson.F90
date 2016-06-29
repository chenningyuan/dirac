module davidson


implicit none

!public davidson_solve_nonsymmetric_R
!public davidson_solve_symmetric_R
public davidson_setup
public davidson_driver
public sort_values_vectors
public set_safe_minimum

private

! module parameters

integer, parameter :: max_label_length = 30
integer, parameter :: nr_max_irreps = 32

real(kind=8), parameter :: zero_r8 = 0.0d0
real(kind=8), parameter :: one_r8  = 1.0d0

! module variables

real(kind=8), save :: safe_minimum = 1.0d-36
real(kind=8), save :: convergence_threshold = 1.0d-2
integer, save      :: max_subspace_size = 2048
integer, save      :: roots_per_symmetry(nr_max_irreps)
integer, save      :: max_iterations = 100

contains 

    subroutine davidson_setup(threshold, subspace_size, nr_roots)
        real(kind=8), intent(in) :: threshold
        integer, intent(in) :: subspace_size 
        integer, intent(in) :: nr_roots(nr_max_irreps)

        integer :: i

        call set_safe_minimum()
        convergence_threshold = threshold
        max_subspace_size = subspace_size  
        do i = 1, nr_max_irreps
            roots_per_symmetry(i) = nr_roots(i)
        end do
    end subroutine

    subroutine set_threshold(threshold)
        real(kind=8), intent(in) :: threshold
        convergence_threshold = threshold
    end

    subroutine set_safe_minimum()
        double precision dlamch
        safe_minimum = dlamch('S')
        print *, " safe minimum is ",safe_minimum
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
                    print "(4x,1e14.6)", vectors(k,i)
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
                    print "(4x,1e14.6)", vectors(k,i)
                end do
                print *, ""
            end if
            do j = i+1, vectors_L
                call gram_schmidt_projection_R(projection, vectors(:,i),vectors(:,j))
                if (verbose_local) then
                    print *," projection of v(",i,") onto v(",j,"):"
                    do k = 1, vectors_N
                        print "(4x,1e14.6,4x,1e14.6)", vectors(k,j), projection(k)
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
            do j = 1, L
                print "(4x,1e14.6,$)", vectors(i,j)
            end do
            print *, ""
        end do
    end subroutine


    subroutine print_values_vectors_R(eValuesRe, eValuesIm, eVectorsR, eVectorsL, left, right, nroots)
        real(kind=8), intent(in) :: eVectorsL(:,:)
        real(kind=8), intent(in) :: eVectorsR(:,:)
        real(kind=8), intent(in) :: eValuesRe(:), eValuesIm(:)
        integer, optional, intent(in) :: nroots

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
                print *, "Eigenvalue ",i," : ",eValuesRe(i)
            else
                print *, "Eigenvalue ",i," : ",eValuesRe(i)," + ",eValuesIm(i), " i "
            end if
        end do
        if (left) then
            label = "Left Eigenvalues            "
            call print_vectors_R(eVectorsL, label=label) 
        end if
        if (right) then
            label = "Right Eigenvalues            "
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

    function correction_to_vectors(nroots, e_small_M, c_small_M, sigma_vector, b, preconditioner, correct_all)
        real(kind=8), pointer     :: correction_to_vectors(:,:)
        integer, intent(in)       :: nroots
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
            if (norm_r .lt. convergence_threshold) then 
               print "(4x,a,i6,a,1e14.6,a,1e14.6,a)", "Root ",k,"  ||q_k,M|| ", norm_r,  &
                                                      " threshold ", convergence_threshold, &
                                                      " <<< converged ! "
               vectors_to_update(k) = .false.
            else
               print "(4x,a,i6,a,1e14.6,a,1e14.6)",  "Root ",k,"  ||q_k,M|| ", norm_r, &
                                                     " threshold ", convergence_threshold
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
                    stop ('detected division by zero')
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
        real(kind=8) :: alpha_lapack = 1.0
        real(kind=8) :: beta_lapack = 1.0

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


    subroutine davidson_driver(nroots, Intm, ptr_to_sigma, ptr_to_diagonal, symmetric, eValuesRe, eValuesIm, eVectorsL, eVectorsR, verbose, b_guess)
        integer, intent(in) :: nroots 
        real(kind=8), intent(in)    :: M(:,:)
        real(kind=8), intent(inout), optional :: eVectorsL(:,:)
        real(kind=8), intent(inout), optional :: eVectorsR(:,:)
        real(kind=8), intent(inout) :: eValuesRe(:), eValuesIm(:)
        real(kind=8), optional, pointer :: b_guess(:,:)
        logical, intent(in), optional :: verbose, symmetric
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
        integer :: vectors_N, L, min_subspace_size, niter
        integer :: i, j
        logical :: left, right, verbose_local, symmetric_local, eVectorsL_local, eVectorsR_local


        if (present(eVectorsL)) then
            left = .true.
        else
            left = .false.
        end if
            
        if (present(eVectorsR)) then
            right = .true.
        else
            right = .false.
        end if

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false. 
        end if

        if (present(symmetric)) then
            symmetric_local = symmetric
        else
            symmetric_local = .false. 
        end if

        vectors_N = size(M,1)
        M_diag => get_diagonal_elements(M)

        min_subspace_size = nroots 
        max_subspace_size = max(nroots, max_subspace_size)
        max_subspace_size = min(vectors_N, max_subspace_size)

        print "(4x,a,i6,a)", "Requested ",nroots," eigenvales"
        if ( nroots .gt. vectors_N ) then 
           print "(4x,a)", "Number of roots requested larger than matrix order !"
           return 
        endif

        if (present(b_guess)) then
            b => b_guess
        else
            b => new_subspace_vectors_R(min_subspace_size,vectors_N)
        end if

        if (verbose) then
            print "(4x,a,i4,a,i8)", "Start by mapping M onto a subspace of size ",min_subspace_size, &
                     ", will extend it up to ",max_subspace_size
            call print_vectors_R(b, label="Initial subspace vectors      ")
        end if

        do niter = 1, max_iterations
            print *, ""
            print *, "Iteration ", niter
            print *, ""

            if (associated(sigma_vector)) deallocate(sigma_vector)
            if (associated(M_small)) deallocate(M_small)

            sigma_vector => create_sigma_vector(M, b)
            M_small => create_reduced_problem(b, sigma_vector)

            call calculate_eigenvalues_and_eigenvectors(M_small, e_small_r, e_small_i, &
                                                c_small_r, c_small_l, left=left, verbose=verbose)

            if (associated(preconditioner)) deallocate(preconditioner)
            if (associated(delta)) deallocate(delta)

            preconditioner =>  form_preconditioner(M_diag, e_small_r)
            delta => correction_to_vectors(nroots, e_small_r, c_small_r, sigma_vector, b, preconditioner)
            if (.not.associated(delta)) then
                print *, "All roots converged, stopping ..."
                exit
            end if
            b => updated_subspace_vectors_R(b, delta, verbose=verbose)
            if (verbose) then
                call print_vectors_R(delta, label="delta                         ")
                call print_vectors_R(b, label="Updated subspace vectors      ")
            end if
            L = size(b,2)
            if (L .eq. max_subspace_size) then
                print *, "Reached maximum subspace size, stopping ..."
                exit
            end if
        end do

        print *, "Final size of reduced subspace : ",L

        if (associated(sigma_vector)) deallocate(sigma_vector)
        if (associated(M_small)) deallocate(M_small)

        sigma_vector => create_sigma_vector(M, b)
        M_small => create_reduced_problem(b, sigma_vector)

        if (verbose) then
            call print_vectors_R(b, label="Final subspace vectors        ")
            print "(4x,a,i6,a)", "After iteration ", niter,",  Eigenvalues for the smaller problem"
        end if
        call calculate_eigenvalues_and_eigenvectors(M_small, e_small_r, e_small_i, &
                                            c_small_r, c_small_l, left=left, verbose=verbose)

        call form_eigenvalues_full_problem(eValuesRe, eValuesIm, nroots, e_small_r, e_small_i) 
        if (left)  call form_eigenvector_full_problem(eVectorsL,nroots,b,c_small_l, .true.)
        if (right) call form_eigenvector_full_problem(eVectorsR,nroots,b,c_small_r, .false.)

        deallocate(b)
        deallocate(M_small)
        deallocate(sigma_vector)
        deallocate(M_diag)

    end subroutine

    subroutine calculate_eigenvalues_and_eigenvectors(M, eValuesRe, eValuesIm, eVectorsR, eVectorsL, &
                                                      left, right, sorting_order, nroots, verbose)
        real(kind=8), intent(inout), pointer   :: M(:,:)
        real(kind=8), intent(inout), pointer   :: eValuesRe(:), eValuesIm(:)
        real(kind=8), intent(inout), pointer   :: eVectorsR(:,:), eVectorsL(:,:)
        logical, optional, intent(in) :: left, right, verbose
        integer, optional, intent(in) :: sorting_order, nroots

        logical :: left_local, right_local, verbose_local
        integer :: sorting_order_local, nroots_local

        if (present(sorting_order)) then
            sorting_order_local = sorting_order
        else
            sorting_order_local = 1 
        end if

        if (present(right)) then
            right_local = right
        else
            right_local = .true.
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
        call sort_values_vectors(sorting_order_local, eValuesRe, eValuesIm, eVectorsR, eVectorsL, left_local)
        if (verbose_local) &
        call print_values_vectors_R(eValuesRe, eValuesIm, eVectorsR, eVectorsL, left_local, right_local, nroots=nroots_local)

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
        real(kind=8) :: alpha_lapack = 1.0
        real(kind=8) :: beta_lapack = 0.0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (left_vector) then
            allocate(bc(vectors_L,vectors_N))
            bc = zero_r8
            call print_vectors_R(eVectorsSmall, label="L eigenvectors small problem  ")
            call dgemm('T', 'T', vectors_L, vectors_N, vectors_L, alpha_lapack, eVectorsSmall, &
                       vectors_N, b, vectors_N, beta_lapack, bc, vectors_L)
            call print_vectors_R(bc, label="bc                            ")
            do i = 1, nroots
                eVectorsFull(:,i) = bc(i,:) 
            end do

        else
            allocate(bc(vectors_N,vectors_L))
            bc = zero_r8
            call print_vectors_R(eVectorsSmall, label="R eigenvectors small problem  ")
            call dgemm('N', 'N', vectors_N, vectors_L, vectors_L, alpha_lapack, b, &
                       vectors_N, eVectorsSmall, vectors_L, beta_lapack, bc, vectors_N)

            call print_vectors_R(bc, label="bc                            ")
            do i = 1, nroots
                eVectorsFull(:,i) = bc(:,i) 
            end do
        end if

        deallocate (bc)
    end subroutine

    function create_sigma_vector(M, b)
        real(kind=8), pointer :: create_sigma_vector(:,:)
        real(kind=8), intent(in)  :: M(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        real(kind=8) :: alpha_lapack = 1.0
        real(kind=8) :: beta_lapack = 0.0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_sigma_vector(vectors_N,vectors_L))
        call dgemm('N', 'N', vectors_N, vectors_L, vectors_N, alpha_lapack, M, vectors_N, &
                   b, vectors_N, beta_lapack, create_sigma_vector, vectors_N)

    end function

    subroutine sort_values_vectors(order, e_small_r, e_small_i, c_small_r, c_small_l, left) 
        integer, intent(in) :: order ! 1: ascending, -1: descending
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

        do i = 1, vectors_L
         i_sorted = indices_sorted_array(i)
         temp_vector(:,i) = c_small_r(:,i_sorted)
        end do
        c_small_r = temp_vector

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

end module

