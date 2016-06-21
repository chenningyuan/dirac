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

    function updated_subspace_vectors_R(old_vectors_R, delta, nroots, verbose)
        real(kind=8), pointer:: updated_subspace_vectors_R(:,:)
        real(kind=8), allocatable :: c(:,:)

        real(kind=8), intent(inout), pointer :: old_vectors_R(:,:)
        real(kind=8), intent(in)      :: delta(:,:)
        integer, intent(in)           :: nroots
        logical, optional, intent(in) :: verbose

        integer :: vectors_L, vectors_N, new_L
        integer :: i
        logical :: verbose_local

        vectors_N = size(old_vectors_R,1)
        vectors_L = size(old_vectors_R,2)

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        new_L = vectors_L + nroots
        updated_subspace_vectors_R => new_subspace_vectors_R(new_L,vectors_N,old_subspace_vectors_R=old_vectors_R)
        do i = 1, nroots
           updated_subspace_vectors_R(:, vectors_L + i) = delta(:,i)
        end do

        if (verbose_local) call print_vectors_R(updated_subspace_vectors_R, "new subspace vectors before GS")
        call gram_schmidt_orthogonalization(updated_subspace_vectors_R, verbose=verbose_local)
        if (associated(old_vectors_R)) deallocate(old_vectors_R)
    end function

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
        if (.true.) then
            print *, "GS projection, factor: ",factor
        end if
        projection = factor * u
    end subroutine 

    subroutine gram_schmidt_orthogonalization(vectors, verbose)
        real(kind=8), intent(inout) :: vectors(:,:)
        logical, optional, intent(in) :: verbose

        real(kind=8), allocatable   :: projection(:)
        real(kind=8), allocatable   :: norms(:)

        integer :: vectors_L, vectors_N
        integer :: i, j, k
        logical :: verbose_local
        real(kind=8) :: norm

        vectors_L = size(vectors,2)
        vectors_N = size(vectors,1)

!       allocate(norms(vectors_L))
        allocate(projection(vectors_N))

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

!       norms = zero_r8 

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
            vectors(:,i) = vectors(:,i) / norm
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

!       do i = 1, vectors_L
!           norms(i) = vectors(:,i) / snrm2(vectors_N, vectors(:,i), 1)
!       end do

        deallocate(projection)
!       deallocate(norms)
    end subroutine 

    subroutine print_vectors_R(vectors, label)
        real(kind=8), intent(in) :: vectors(:,:)
        character(len=max_label_length), intent(in), optional :: label 

        integer :: L, N, i, j
        character(len=max_label_length) :: label_local

        N = size(vectors,1)
        L = size(vectors,2)

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


    subroutine print_eigenvalues_eigenvectors_R(eValuesRe, eValuesIm, eVectorsR, eVectorsL, left, right, nroots)
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
        do i = 1, vectors_N
            if (eValuesIm(i) .eq. zero_r8) then
                print *, "Eigenvalue ",i," : ",eValuesRe(i)
            else
                print *, "Eigenvalue ",i," : ",eValuesRe(i)," + ",eValuesIm(i), " i "
            end if
        end do
        if (right) then
            label = "Right Eigenvalues            "
            call print_vectors_R(eVectorsR, label) 
        end if
        if (left) then
            label = "Left Eigenvalues            "
            call print_vectors_R(eVectorsL, label) 
        end if


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

!       print *, "debug, before first dgemm call"
        call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                   alpha_lapack, sigma_vector, vectors_N, c_small_M, vectors_L, &
                   beta_lapack, residue_a, vectors_N) 
!       print *, "debug, before second dgemm call"
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

    function vector_norm(v)
        real(kind=8) :: vector_norm
        real(kind=8), intent(in) :: v(:)
        integer :: vector_N, i
        double precision ddot
        double precision dnrm2

        vector_N = size(v,1)
        vector_norm = dnrm2(vector_N, v, 1)
!       vector_norm = ddot(vector_N, v, 1, v, 1)
!       vector_norm = dsqrt(vector_norm)

    end function

    function correction_vectors(nroots, M_diag, e_small_M, c_small_M, sigma_vector, b)
        real(kind=8), pointer     :: correction_vectors(:,:)
        integer, intent(in)       :: nroots
        real(kind=8), intent(in)  :: M_diag(:)
        real(kind=8), intent(in)  :: e_small_M(:)
        real(kind=8), intent(in)  :: c_small_M(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        integer :: k, converged_roots
        real(kind=8) :: norm_r = zero_r8
        real(kind=8), allocatable :: r(:,:)

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(r(vectors_N,nroots))
!       print *, "debug, calling form_r"
        call form_r(r, nroots, e_small_M, c_small_M, sigma_vector, b)

        converged_roots = 0
        do k = 1, nroots
            norm_r = vector_norm(r(:,k))
            if (norm_r .lt. convergence_threshold) then 
               print "(4x,a,i6,a,1e14.6,a,1e14.6,a)", "root ",k,"  ||q_k,M|| ", norm_r,  &
                                                      " threshold ", convergence_threshold, &
                                                      " <<< converged ! "
               converged_roots = converged_roots + 1
            else
               print "(4x,a,i6,a,1e14.6,a,1e14.6)",  "root ",k,"  ||q_k,M|| ", norm_r, &
                                                     " threshold ", convergence_threshold
            end if
        end do
        print *, ""

        allocate(correction_vectors(vectors_N,nroots))
        correction_vectors = zero_r8
        if (converged_roots .ne. nroots) then
!           print *, "debug, calling form_delta"
            call form_delta(correction_vectors, nroots, r, M_diag, e_small_M)
        end if

        deallocate(r)
    end function

    function is_converged(correction_vectors)
        logical :: is_converged
        real(kind=8) :: norm_d, correction_vectors(:,:)
        integer :: i, vectors_L

        vectors_L = size(correction_vectors,2)

        norm_d = zero_r8
        do i = 1, vectors_L
            norm_d = vector_norm(correction_vectors(:,i))
            if (norm_d .gt. convergence_threshold) then
                is_converged = .false.
            else
                is_converged = .true.
            end if
        end do
    
    end function

    subroutine form_delta(delta, nroots, r, M_diag, e_small_M)
        real(kind=8), intent(inout) :: delta(:,:)
        integer, intent(in)       :: nroots
        real(kind=8), intent(in)  :: e_small_M(:)
        real(kind=8), intent(in)  :: M_diag(:)
        real(kind=8), intent(in)  :: r(:,:)

        integer :: vectors_N
        integer :: k
        real(kind=8) :: denominator

        vectors_N = size(delta,1)

        do k = 1, nroots
            denominator = e_small_M(k) - M_diag(k)
            if (dabs(denominator) .lt. safe_minimum) then
                print *,"error !"
                denominator = 1000
            else
                denominator = one_r8 / denominator
            end if
           delta(:,k) = denominator*r(:,k)
           delta(:,k) = delta(:,k) / vector_norm(delta(:,k))
       end do
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
!       call print_vectors_R(create_reduced_problem, "M_small, before construction  ")
        call dgemm('T', 'N', vectors_L, vectors_L, vectors_N, alpha_lapack, b, vectors_N, &
                   sigma_vector, vectors_N, beta_lapack, create_reduced_problem, vectors_L)
!       call print_vectors_R(create_reduced_problem, "M_small, after  construction  ")

    end function

    subroutine solve_small_eigenvalue(A_small, eValuesRe, eValuesIm, eVectorsR, eVectorsL, left)
        real(kind=8), intent(inout), pointer :: A_small(:,:)
        real(kind=8), intent(inout), pointer :: eVectorsL(:,:)
        real(kind=8), intent(inout), pointer :: eVectorsR(:,:)
        real(kind=8), intent(inout), pointer :: eValuesRe(:), eValuesIm(:)
        logical, intent(in) :: left 

        integer :: vectors_L, i

        integer :: lwork, info
        character(len=1) :: do_left = 'N', do_right = 'V'
        real(kind=8), allocatable :: work(:)

!       print *, " > debug, at the start of solve_small_eigenvalue()"
        vectors_L = size(A_small,1)
!       print *, "    size of A_small :",vectors_L

        if (associated(eVectorsL)) deallocate(eVectorsL)
        if (associated(eVectorsR)) deallocate(eVectorsR)
        if (associated(eValuesRe)) deallocate(eValuesRe)
        if (associated(eValuesIm)) deallocate(eValuesIm)

        if (left) then
            do_left = 'V'
            allocate(eVectorsL(vectors_L,vectors_L))
            eVectorsL = zero_r8
!           print *,"after alloc VL"
        end if
        allocate(eVectorsR(vectors_L,vectors_L))
!           print *,"after alloc VR"
        allocate(eValuesRe(vectors_L))
!           print *,"after alloc er"
        allocate(eValuesIm(vectors_L))
!           print *,"aftoer alloc ei"

        eVectorsR = zero_r8
        eValuesRe = zero_r8
        eValuesIm = zero_r8

! first query optimal lwork value
        lwork = -1
!       print *,"before alloc work #1 size work = ",(4*vectors_L) 
        allocate(work(4*vectors_L))
!       print *,"after  alloc work #1"
        call dgeev (do_left, do_right, vectors_L, A_small, vectors_L,      &
                    eValuesRe, eValuesIm, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, INFO)
        lwork = work(1)
!       print *,"before dealloc work #1, size work = ",lwork
        deallocate(work)
!       print *,"after  dealloc work #1"
!       print *,"before alloc work #2"
        allocate(work(lwork))
!       print *,"after  alloc work #2"
        call dgeev (do_left, do_right, vectors_L, A_small, vectors_L,      &
                    eValuesRe, eValuesIm, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, INFO)
!       print *,"before dealloc work #2"
        deallocate(work)
!       print *,"after  dealloc work #2"
        if (info .gt. 0) then
            print *, 'in small_solve_eigenvalue: not all eigenvalues converged' 
            print *, 'converged'
            do i = info+1, vectors_L 
                print *, " real part : ",eValuesRe(i), " imaginary part : ",eValuesIm(i)
            end do
        end if
        if (.false.) then
            do i = 1, vectors_L 
                print *, " real      part : ",eValuesRe(i)," imaginary part : ",eValuesIm(i)
            end do
        end if
!       print *, " > debug, at the end   of solve_small_eigenvalue()"

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


    subroutine davidson_driver(nroots, M, symmetric, eValuesRe, eValuesIm, eVectorsL, eVectorsR, left, right, verbose)
        real(kind=8), intent(in)    :: M(:,:)
        real(kind=8), intent(inout) :: eVectorsL(:,:)
        real(kind=8), intent(inout) :: eVectorsR(:,:)
        real(kind=8), intent(inout) :: eValuesRe(:), eValuesIm(:)
        integer, intent(in) :: nroots 
        logical, intent(in) :: left, right, verbose, symmetric
!
! nullifying all pointers at the start to give them a state, see discussion on http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html
! "[...] When a pointer is declared its status is undefined, and cannot be safely queried with the associated intrinsic. [...]"
! 
        real(kind=8), pointer :: M_diag(:)           => NULL()
        real(kind=8), pointer :: delta(:,:)          => NULL()
        real(kind=8), pointer :: sigma_vector(:,:)   => NULL()
        real(kind=8), pointer :: b(:,:)              => NULL()
        real(kind=8), pointer :: M_small(:,:)        => NULL()
        real(kind=8), pointer :: c_small_r(:,:)      => NULL()
        real(kind=8), pointer :: c_small_l(:,:)      => NULL()
        real(kind=8), pointer :: e_small_r(:)        => NULL()
        real(kind=8), pointer :: e_small_i(:)        => NULL()
        integer :: vectors_N, L, min_subspace_size, niter
        integer :: i, j
        integer :: sorting_order = 1


        vectors_N = size(M,1)
        M_diag => get_diagonal_elements(M)

        min_subspace_size = nroots 
        max_subspace_size = max(nroots, max_subspace_size)
        max_subspace_size = min(vectors_N, max_subspace_size)
        max_subspace_size = max_subspace_size - nroots

        print "(4x,a,i6,a)", "Requested ",nroots," eigenvales"
        if ( nroots .gt. vectors_N ) then 
           print "(4x,a)", "Number of roots requested larger than matrix order !"
           return 
        endif

        b => new_subspace_vectors_R(min_subspace_size,vectors_N)
        if (verbose) then
            print "(4x,a,i4,a,i8)", "Start by mapping M onto a subspace of size ",min_subspace_size, &
                     ", will extend it up to ",max_subspace_size
            call print_vectors_R(b, "Initial subspace vectors      ")
        end if

        niter = 1
        if (min_subspace_size.le.max_subspace_size) L = min_subspace_size 

        do L = min_subspace_size, max_subspace_size, nroots 
            print *, ""
            print *, "Iteration ", niter
            print *, ""

            if (associated(sigma_vector)) deallocate(sigma_vector)
            if (associated(M_small)) deallocate(M_small)

            sigma_vector => create_sigma_vector(M, b)
            M_small => create_reduced_problem(b, sigma_vector)

            call solve_small_eigenvalue(M_small, e_small_r, e_small_i, c_small_r, c_small_l, left)
            call sort_values_vectors(sorting_order,e_small_r, e_small_i, c_small_r, c_small_l, left) 

!           if (verbose) &
            if (.true.) &
            call print_eigenvalues_eigenvectors_R(e_small_r, e_small_i, c_small_r, c_small_l, left, right, nroots=nroots)

            delta => correction_vectors(nroots, M_diag, e_small_r, c_small_r, sigma_vector, b)
            if (verbose) call print_vectors_R(delta, "delta                         ")

            if (is_converged(delta)) then
                print *, "All roots converged, stopping ..."
                exit 
            else
                print *, "not all roots converged, continuing ..."
            end if

            b => updated_subspace_vectors_R(b, delta, nroots, verbose=.true.)
            if (verbose) call print_vectors_R(b, "Updated subspace vectors      ")
            niter = niter + 1
        end do

        print *, "Final size of reduced subspace : ",L

        if (associated(M_small)) deallocate(M_small)
        if (associated(sigma_vector)) deallocate(sigma_vector)

        sigma_vector => create_sigma_vector(M, b)
        M_small => create_reduced_problem(b, sigma_vector)

        call solve_small_eigenvalue(M_small, e_small_r, e_small_i, c_small_r, c_small_l, .false.)
        call sort_values_vectors(sorting_order,e_small_r, e_small_i, c_small_r, c_small_l, .false.) 

        print "(4x,a,i6,a)", "After iteration ", (niter-1),",  Eigenvalues for the smaller problem"
        call print_eigenvalues_eigenvectors_R(e_small_r, e_small_i, c_small_r, c_small_l, left, right, nroots=nroots)

        call print_vectors_R(b, "Final subspace vectors        ")
        call form_eigenvalues_full_problem(eValuesRe, eValuesIm, nroots, e_small_r, e_small_i) 
        if (right) call form_eigenvector_full_problem(eVectorsR,nroots,b,c_small_r, .false.)
        if (left)  call form_eigenvector_full_problem(eVectorsL,nroots,b,c_small_l, .true.)

        deallocate(b)
        deallocate(M_small)
        deallocate(sigma_vector)
        deallocate(M_diag)

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
            call dgemm('N', 'T', vectors_L, vectors_N, vectors_L, alpha_lapack, eVectorsSmall, &
                       vectors_N, b, vectors_N, beta_lapack, bc, vectors_L)

            do i = 1, nroots
                eVectorsFull(:,i) = bc(i,:) 
            end do

        else
            allocate(bc(vectors_N,vectors_L))
            bc = zero_r8
            call print_vectors_R(eVectorsSmall, label="eigenvectors small problem    ")
            call dgemm('N', 'N', vectors_N, vectors_L, vectors_L, alpha_lapack, b, &
                       vectors_N, eVectorsSmall, vectors_L, beta_lapack, bc, vectors_N)

            call print_vectors_R(bc, label="bc                            ")
            do i = 1, nroots
                eVectorsFull(:,i) = bc(:,i) 
            end do
        end if

        deallocate (bc)
    end subroutine

    function create_sigma_vector(A, b)
        real(kind=8), pointer :: create_sigma_vector(:,:)
        real(kind=8), intent(in)  :: A(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        real(kind=8) :: alpha_lapack = 1.0
        real(kind=8) :: beta_lapack = 0.0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_sigma_vector(vectors_N,vectors_L))
        call dgemm('N', 'N', vectors_N, vectors_L, vectors_N, alpha_lapack, A, vectors_N, &
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
!       print *, "blah!"
!       print *,"original order, ",e_small_r
        call quicksort_vector(e_small_r, indices_sorted_array, 1, vectors_L)
!       print *,"sorted eigenvalues",e_small_r
!       print *,"mapping between original and sorted indices",indices_sorted_array
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

recursive subroutine quicksort_vector(a, indices_sorted_array, first, last)
! routine taken from https://gist.github.com/t-nissie/479f0f16966925fa29ea; this is GPL so this has to be replaced by a clean version
  real(kind=8)  a(:), x, t
  integer indices_sorted_array(:)
  integer first, last
  integer i, j, ti

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     ti = indices_sorted_array(i)
     indices_sorted_array(i) = indices_sorted_array(j)
     indices_sorted_array(j) = ti
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort_vector(a, indices_sorted_array, first, i-1)
  if (j+1 < last)  call quicksort_vector(a, indices_sorted_array, j+1, last)
end subroutine 


!   subroutine davidson_solve_nonsymmetric_R
!   
!   end subroutine

end module

