module davidson


implicit none

!public davidson_solve_nonsymmetric_R
!public davidson_solve_symmetric_R
public davidson_setup
public davidson_driver

private

! module parameters

integer, parameter :: max_label_length = 10
integer, parameter :: nr_max_irreps = 32

real(kind=8), parameter :: zero_r8 = 0.0d0
real(kind=8), parameter :: one_r8  = 1.0d0

! module variables

real(kind=8), save :: safe_minimum = 1.0d-12
real(kind=8), save :: convergence_threshold = 1.0d-4
integer, save      :: max_subspace_size = 128
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
    end subroutine 


    function new_subspace_vectors_R(new_L, N, old_subspace_vectors_R, verbose)
!   new_L is the dimension of the subspace to be created 
!   N is the dimension of the original problem
!   new_subspace_vectors_R  is the subspace we create or update
!   old_subspace_vectors_R (if present) is the old subspace
        real(kind=8), allocatable :: new_subspace_vectors_R(:,:)
        integer, intent (in) :: new_L, N
        real(kind=8), intent(inout), optional, allocatable, target :: old_subspace_vectors_R(:,:)
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
        new_subspace_vectors_R = 0.0d0

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
        real(kind=8), allocatable     :: updated_subspace_vectors_R(:,:), c(:,:)

        real(kind=8), intent(inout), allocatable :: old_vectors_R(:,:)
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
        updated_subspace_vectors_R = new_subspace_vectors_R(new_L,vectors_N,old_subspace_vectors_R=old_vectors_R)

        do i = 1, nroots
           updated_subspace_vectors_R(:, vectors_L + i) = delta(:,i)
        end do

        if (allocated(old_vectors_R)) deallocate(old_vectors_R)
    end function

    subroutine gram_schmidt_projection_R(projection, u, v)
        real(kind=8), intent(out)   :: projection(:)
        real(kind=8), intent(inout) :: u(:)
        real(kind=8), intent(in)    :: v(:)
        real(kind=8) :: numerator, denominator, factor
        integer :: N
        double precision ddot

        N = size(u,1)

        denominator = ddot(N, u, 1, u, 1) 
        numerator   = ddot(N, v, 1, u, 1) 

        projection = 0.0d0
        if (denominator .gt. safe_minimum) then
            factor = numerator/denominator
        else
            factor = 0.0
        end if
        projection = factor * u
    end subroutine 

    subroutine gram_schmidt_orthogonalization(vectors, verbose)
        real(kind=8), intent(inout) :: vectors(:,:)
        logical, optional, intent(in) :: verbose

        real(kind=8), allocatable   :: projection(:)
        real(kind=8), allocatable   :: norms(:)

        integer :: vectors_L, vectors_N
        integer :: i, j
        logical :: verbose_local

        vectors_L = size(vectors,2)
        vectors_N = size(vectors,1)

        allocate(norms(vectors_L))
        allocate(projection(vectors_L))

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        norms = 0.0d0

        do i = 1, vectors_L
            if (verbose_local) then
                print *, "debug, Gram-Schmidt orthogonalization"
                print *, "  renormalizing vector",vectors(:,i)
            end if 
            vectors(:,i) = vectors(:,i) / vector_norm(vectors(:,i))
            if (verbose_local) then
                print *, "renormalized  vector v(",i,") :",vectors(:,i)
            end if
            do j = i+1, vectors_L
                if (verbose_local) then
                    print *," projection of v(",i,") onto v(",j,"):",vectors(:,j)
                end if
                call gram_schmidt_projection_R(projection, vectors(:,i),vectors(:,j))
                vectors(:,j) = vectors(:,j) - projection(:)
            end do
            if (verbose_local) then
                call print_vectors_R(vectors, label="all vectors")
            end if
        end do

!       do i = 1, vectors_L
!           norms(i) = vectors(:,i) / snrm2(vectors_N, vectors(:,i), 1)
!       end do

        deallocate(projection)
        deallocate(norms)
    end subroutine 

    subroutine print_vectors_R(vectors, label)
        real(kind=8), intent(in) :: vectors(:,:)
        character(len=max_label_length), intent(in), optional :: label 

        integer :: L, N, i, j
        character(len=max_label_length) :: label_local

        L = size(vectors,2)
        N = size(vectors,1)

        if (present(label)) then
            label_local = label
            print "(a)", label_local(1:max_label_length)
        end if
        do i = 1, N
            print "(4x,1f14.6,a,$)" , (vectors(i,j), j= 1,L)
            print *, ""
        end do
    end subroutine


    subroutine print_eigenvalues_eigenvectors_R(eValuesRe, eValuesIm, eVectorsR, eVectorsL, left, right)
        real(kind=8), intent(in) :: eVectorsL(:,:)
        real(kind=8), intent(in) :: eVectorsR(:,:)
        real(kind=8), intent(in) :: eValuesRe(:), eValuesIm(:)

        logical :: left, right
        integer :: vectors_N
        integer :: i
        character(len=max_label_length) :: label

        vectors_N = size(eValuesRe,1)
        label = "Eigenvalues"
        do i = 1, vectors_N
            if (eValuesIm(i) .eq. zero_r8) then
                print *, "Eigenvalue ",i," : ",eValuesRe(i)
            else
                print *, "Eigenvalue ",i," : ",eValuesRe(i)," + ",eValuesIm(i), " i "
            end if
        end do
        if (right) then
            label = "Right Eigenvalues"
            call print_vectors_R(eVectorsR, label) 
        end if
        if (left) then
            label = "Left Eigenvalues"
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

        call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                   alpha_lapack, sigma_vector, vectors_N, c_small_M, vectors_L, &
                   beta_lapack, residue_a, vectors_L) 
        call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                   alpha_lapack, b, vectors_L, c_small_M, vectors_L,            &
                   beta_lapack, residue_b, vectors_L)
        if (.true.) then
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

        vector_N = size(v,1)
        vector_norm = ddot(vector_N, v, 1, v, 1)
        vector_norm = dsqrt(vector_norm)

    end function

    function correction_vectors(nroots, M_diag, e_small_M, c_small_M, sigma_vector, b)
        real(kind=8), allocatable :: correction_vectors(:,:)
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
        allocate(correction_vectors(vectors_N,nroots))

        call form_r(r, nroots, e_small_M, c_small_M, sigma_vector, b)

        converged_roots = 0
        do k = 1, nroots
            norm_r = vector_norm(r(:,k))
            if (norm_r .lt. convergence_threshold) then 
               print *, "root ",k,"  ||q_k,M|| ",norm_r, " threshold ", convergence_threshold, " <<< converged ! "
               converged_roots = converged_roots + 1
            else
               print *, "root ",k,"  ||q_k,M|| ", norm_r, " threshold ", convergence_threshold
            end if
       end do
       print *, ""

       if (converged_roots .ne. nroots) then
           call form_delta(correction_vectors, nroots, r, M_diag, e_small_M)
       end if

       deallocate(r)
    end function

    function is_converged(correction_vectors)
        logical :: is_converged
        real(kind=8) :: norm_d, correction_vectors(:,:)
        integer :: i, vectors_N

        vectors_N = size(correction_vectors,1)

        norm_d = zero_r8
        do i = 1, vectors_N
            norm_d = norm_d + vector_norm(correction_vectors(:,i))
        end do
        if (norm_d .gt. zero_r8) then
            is_converged = .false.
        else
            is_converged = .true.
        end if 
    
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
        real(kind=8), allocatable :: create_reduced_problem(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)

        integer :: vectors_L, vectors_N
        real(kind=8) :: alpha_lapack = 1.0
        real(kind=8) :: beta_lapack = 1.0

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        print *, "debug, in create_reduced_problem, dimensions ",vectors_L," by ", vectors_L
        allocate(create_reduced_problem(vectors_L,vectors_L))
        call print_vectors_R(create_reduced_problem, "M_small, before construction")
        call dgemm('T', 'N', vectors_L, vectors_L, vectors_N, alpha_lapack, b, vectors_N, &
                   sigma_vector, vectors_N, beta_lapack, create_reduced_problem, vectors_L)
        call print_vectors_R(create_reduced_problem, "M_small, after  construction")

    end function

    subroutine solve_small_eigenvalue(A_small, eValuesRe, eValuesIm, eVectorsR, eVectorsL, left)
        real(kind=8), intent(inout) :: A_small(:,:)
        real(kind=8), intent(inout), allocatable :: eVectorsL(:,:)
        real(kind=8), intent(inout), allocatable :: eVectorsR(:,:)
        real(kind=8), intent(inout), allocatable :: eValuesRe(:), eValuesIm(:)
        logical, intent(in) :: left 

        integer :: vectors_L, i

        integer :: lwork, info
        character(len=1) :: do_left = 'N', do_right = 'V'
        real(kind=8), allocatable :: work(:)

        vectors_L = size(A_small,1)
        if (left) then
            do_left = 'V'
            if (.not.allocated(eVectorsL)) allocate(eVectorsL(vectors_L,vectors_L))
        end if
        if (.not.allocated(eVectorsR)) allocate(eVectorsR(vectors_L,vectors_L))
        if (.not.allocated(eValuesRe)) allocate(eValuesRe(vectors_L))
        if (.not.allocated(eValuesIm)) allocate(eValuesIm(vectors_L))

! first query optimal lwork value
        lwork = -1
        allocate(work(4*vectors_L))
        call dgeev (do_left, do_right, vectors_L, A_small, vectors_L,      &
                    eValuesRe, eValuesIm, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, INFO)
        lwork = work(1)
        deallocate(work)
        allocate(work(lwork))
        call dgeev (do_left, do_right, vectors_L, A_small, vectors_L,      &
                    eValuesRe, eValuesIm, eVectorsL, vectors_L, eVectorsR, &
                    vectors_L, WORK, LWORK, INFO)
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
        real(kind=8), allocatable :: get_diagonal_elements(:)

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

        real(kind=8), allocatable :: M_diag(:)
        real(kind=8), allocatable :: delta(:,:)
        real(kind=8), allocatable :: sigma_vector(:,:)
        real(kind=8), allocatable :: b(:,:)
        real(kind=8), allocatable :: M_small(:,:)
        real(kind=8), allocatable :: c_small_r(:,:)
        real(kind=8), allocatable :: c_small_l(:,:)
        real(kind=8), allocatable :: e_small_r(:)
        real(kind=8), allocatable :: e_small_i(:)
        integer :: vectors_N, L, min_subspace_size, niter
        integer :: i, j

        vectors_N = size(M,1)
        M_diag = get_diagonal_elements(M)

        min_subspace_size = nroots 
        max_subspace_size = max(nroots, max_subspace_size)

        print *, "Requested ",nroots," eigenvales"
        if ( nroots .gt. vectors_N ) then 
           print *, "Number of roots requested larger than matrix order !"
           return 
        endif

        b  = new_subspace_vectors_R(min_subspace_size,vectors_N)
        if (verbose) then
            print "(4x,a,i4,a,i8)", "Start by mapping M onto a subspace of size ",min_subspace_size, &
                     ", will extend it up to ",max_subspace_size
            print *, "Initial subspace vectors : "
            do i = 1, vectors_N
                print "(4x,1f14.6,a,$)", (b(i,j), j = 1, min_subspace_size)
                print *, ""
            end do
            print *, ""
        end if

        niter = 1
        do L = min_subspace_size, max_subspace_size, nroots 
            if (allocated(M_small)) deallocate(M_small)
            if (allocated(sigma_vector)) deallocate(sigma_vector)

            sigma_vector = create_sigma_vector(M, b)
            M_small = create_reduced_problem(b, sigma_vector)

            call solve_small_eigenvalue(M_small, e_small_r, e_small_i, c_small_r, c_small_l, .false.)
!           call sort_values_vectors(e_small_r, e_small_i, c_small_r, c_small_l, .false.) 

            print *, "Iteration ", niter,",  projection onto the smaller problem"
            call print_eigenvalues_eigenvectors_R(e_small_r, e_small_i, c_small_r, c_small_l, .false., .false.)
            print *, "Calculating correction to subspace"

            delta = correction_vectors(nroots, M_diag, e_small_r, c_small_r, sigma_vector, b)
            if (is_converged(delta)) exit 

            b = updated_subspace_vectors_R(b, delta, nroots)
            niter = niter + 1
        end do

        print *, "Final size of reduced subspace : ",L

        if (allocated(M_small)) deallocate(M_small)
        if (allocated(sigma_vector)) deallocate(sigma_vector)

        sigma_vector = create_sigma_vector(M, b)
        M_small = create_reduced_problem(b, sigma_vector)

        call solve_small_eigenvalue(M_small, e_small_r, e_small_i, c_small_r, c_small_l, .false.)
!       call sort_values_vectors(e_small_r, e_small_i, c_small_r, c_small_l, .false.) 

        print *, "After iteration ", (niter-1),",  Eigenvalues for the smaller problem"
        call print_eigenvalues_eigenvectors_R(e_small_r, e_small_i, c_small_r, c_small_l, .false., .false.)

        deallocate(b)
        deallocate(M_small)
        deallocate(sigma_vector)
        deallocate(M_diag)
    end subroutine

    function create_sigma_vector(A, b)
        real(kind=8), allocatable :: create_sigma_vector(:,:)
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


!   subroutine davidson_solve_nonsymmetric_R
!   
!   end subroutine

end module

