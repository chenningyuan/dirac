! Unit test of the Foldy-Wouthuysen (FW) coefficient construction in foldy_wouthuysen.F90.
!
! The FW coefficients are constructed from a set of orthonormal 4-component spinors with large-component and
! small-component coefficients CL and CS and the overlap matrices SLL and SSS of the two basis sets.  Checks:
!   1. the Gram-matrix formulation (construct_fw_coefficients) agrees with the former algorithm, based on
!      SLL^(+-1/2) and SSS^(+-1/2) (reference_fw below, taken from the previous version of foldy_wouthuysen.F90),
!      when the latter is well defined
!   2. the FW coefficients are orthonormal in the large-component metric:  CFW^+ SLL CFW = CL^+ SLL CL + CS^+ SSS CS = 1
!   3. the construction stays accurate when SSS is numerically singular (nearly linearly dependent small-component
!      functions as in uncontracted Cartesian Dyall sets), where the former algorithm stops
! All three algebras of linear_algebra.F90 are tested: real (rcq=1), complex (rcq=2) and quaternion (rcq=4).
! The test is self-checking: the exit status is non-zero if any check fails.
program test_fw
  use foldy_wouthuysen
  use linear_algebra
  implicit none
  integer, parameter :: dp = kind(1.D0)
  real(dp), parameter :: tol_ref  = 1.D-10  ! new vs reference implementation
  real(dp), parameter :: tol_orth = 1.D-10  ! orthonormality of the FW coefficients (also for a singular SSS)
  integer :: rcq, nfail
  character(len=16) :: mode

  ! With an argument the program runs one deliberately invalid input and is expected to STOP in a guard of
  ! construct_fw_coefficients (used by ctest with WILL_FAIL and the expected message):
  !   norm    spinors not orthonormal in the four-component metric
  !   count   more spinors than large-component basis functions
  !   nolarge a spinor without large component
  call get_command_argument(1, mode)
  if (len_trim(mode) > 0) then
     call bad_input(trim(mode))
     write(*,'(A)') 'test_fw: the guard did not stop the program'   ! not reached if the guard works
     stop 0
  end if

  nfail = 0
  do rcq = 1, 4
     if (rcq == 3) cycle
     call check(rcq, 30, 50, 6, .false., nfail)
     call check(rcq, 60, 90, 14, .false., nfail)
     call check(rcq, 60, 90, 14, .true., nfail)
  end do
  if (nfail > 0) then
     write(*,'(A,I0,A)') 'test_fw: ', nfail, ' check(s) FAILED'
     error stop 1
  end if
  write(*,'(A)') 'test_fw: all checks passed'

contains

  subroutine random_overlap(n, shift, S, singular)
    ! symmetric positive (semi)definite overlap in component 1, other components zero (as in ROSE)
    integer, intent(in) :: n
    real(dp), intent(in) :: shift
    real(dp), intent(out) :: S(:,:,:)
    logical, intent(in) :: singular
    real(dp) :: A(n,n)
    integer :: i
    call random_number(A); A = A - 0.5D0
    if (singular) then
       ! three functions nearly linearly dependent on others (differences 1e-9) -> eigenvalues of S ~ 1e-18..1e-17
       A(:,n)   = A(:,n-1)   + 1.D-9*A(:,n)
       A(:,n-2) = A(:,n-3)   + 1.D-9*A(:,n-2)
       A(:,n-4) = A(:,n-5)   + 1.D-9*A(:,n-4)
    end if
    S = 0.D0
    S(:,:,1) = matmul(A, transpose(A))
    if (.not. singular) then
       do i = 1, n
          S(i,i,1) = S(i,i,1) + shift
       end do
    end if
  end subroutine

  subroutine reference_fw(CL, CS, SLL, SSS, CFW)
    ! Former algorithm (construct_fw_coefficients of foldy_wouthuysen.F90 before the Gram formulation): transform to
    ! orthonormal large and small bases with SLL^(1/2), SSS^(1/2), construct X and R there, transform back.
    real(dp), intent(in) :: CL(:,:,:), CS(:,:,:), SLL(:,:,:), SSS(:,:,:)
    real(dp), allocatable, intent(out) :: CFW(:,:,:)
    real(dp), allocatable :: X(:,:,:), R(:,:,:), SLLs(:,:,:), SLLi(:,:,:), SSSs(:,:,:), SSSi(:,:,:)
    real(dp), allocatable :: CLo(:,:,:), CSo(:,:,:), CFWo(:,:,:), IL(:,:,:), IS(:,:,:)
    integer :: nL, nS, np, rc, i
    nL = size(CL,1); np = size(CL,2); rc = size(CL,3); nS = size(CS,1)
    allocate(X(nS,nL,rc), R(nL,nL,rc), CFW(nL,np,rc), SLLs(nL,nL,rc), SLLi(nL,nL,rc), SSSs(nS,nS,rc), SSSi(nS,nS,rc))
    allocate(CLo(nL,np,rc), CSo(nS,np,rc), CFWo(nL,np,rc), IL(nL,nL,rc), IS(nS,nS,rc))
    call hermitian_matrix_sqrt_pair(SLL, SLLs, SLLi, 'reference_fw')
    call hermitian_matrix_sqrt_pair(SSS, SSSs, SSSi, 'reference_fw')
    call matrix_mul(SLLs, CL, CLo)
    call matrix_mul(SSSs, CS, CSo)
    IL = 0.D0; IS = 0.D0
    do i = 1, nL; IL(i,i,1) = 1.D0; end do
    do i = 1, nS; IS(i,i,1) = 1.D0; end do
    call construct_fw_x_matrix(CLo, CSo, IL, X)
    call construct_fw_r_matrix(X, IL, IS, R)
    call matrix_mul(R, CLo, CFWo)
    call matrix_mul(SLLi, CFWo, CFW)
  end subroutine

  subroutine bad_input(kind_)
    character(len=*), intent(in) :: kind_
    integer, parameter :: nL = 30, nS = 50
    integer :: np, i
    real(dp), allocatable :: SLL(:,:,:), SSS(:,:,:), CL(:,:,:), CS(:,:,:), CbL(:,:,:), CbS(:,:,:)
    real(dp), allocatable :: GL(:,:,:), GS(:,:,:), G(:,:,:), Gsq(:,:,:), Gisq(:,:,:), CLn(:,:,:), CSn(:,:,:), CFW(:,:,:)
    np = 6
    if (kind_ == 'count') np = nL + 4
    allocate(SLL(nL,nL,1), SSS(nS,nS,1), CL(nL,np,1), CS(nS,np,1), CbL(np,nL,1), CbS(np,nS,1))
    allocate(GL(np,np,1), GS(np,np,1), G(np,np,1), Gsq(np,np,1), Gisq(np,np,1), CLn(nL,np,1), CSn(nS,np,1))
    call random_number(CL); call random_number(CS); CL = CL - 0.5D0; CS = CS - 0.5D0
    call random_overlap(nL, 0.5D0, SLL, .false.)
    call random_overlap(nS, 0.5D0, SSS, .false.)
    if (kind_ == 'nolarge') CL(:,1,:) = 1.D-7*CL(:,1,:)
    call transpose_conjg(CL, CbL); call matrix_mul(CbL, SLL, CL, GL)
    call transpose_conjg(CS, CbS); call matrix_mul(CbS, SSS, CS, GS)
    G = GL + GS
    if (kind_ /= 'count') then
       call hermitian_matrix_sqrt_pair(G, Gsq, Gisq, 'test_fw: normalisation')
       call matrix_mul(CL, Gisq, CLn); call matrix_mul(CS, Gisq, CSn)
    else
       CLn = CL; CSn = CS
    end if
    if (kind_ == 'norm') CSn = 1.1D0*CSn          ! orthonormality destroyed by 10% in the small component
    call construct_fw_coefficients(CLn, CSn, SLL, SSS, CFW)
  end subroutine

  subroutine check(rcq_, nL, nS, np, singular, nfail_)
    integer, intent(in) :: rcq_, nL, nS, np
    logical, intent(in) :: singular
    integer, intent(inout) :: nfail_
    real(dp) :: SLL(nL,nL,rcq_), SSS(nS,nS,rcq_), CL(nL,np,rcq_), CS(nS,np,rcq_), CLn(nL,np,rcq_), CSn(nS,np,rcq_)
    real(dp) :: GL(np,np,rcq_), GS(np,np,rcq_), G(np,np,rcq_), Gsq(np,np,rcq_), Gisq(np,np,rcq_), Id(np,np,rcq_)
    real(dp) :: CbL(np,nL,rcq_), CbS(np,nS,rcq_), T(np,np,rcq_)
    real(dp), allocatable :: C1(:,:,:), C2(:,:,:), CbF(:,:,:)
    real(dp) :: d_ref, d_orth
    integer :: i
    logical :: ok
    character(len=24) :: label

    call random_number(CL); call random_number(CS)
    CL = CL - 0.5D0; CS = CS - 0.5D0
    call random_overlap(nL, 0.5D0, SLL, .false.)
    call random_overlap(nS, 0.5D0, SSS, singular)
    Id = 0.D0
    do i = 1, np
       Id(i,i,1) = 1.D0
    end do

    ! make the spinors orthonormal in the four-component metric: CL <- CL G^(-1/2), CS <- CS G^(-1/2), G = GL + GS
    call transpose_conjg(CL, CbL); call matrix_mul(CbL, SLL, CL, GL)
    call transpose_conjg(CS, CbS); call matrix_mul(CbS, SSS, CS, GS)
    G = GL + GS
    call hermitian_matrix_sqrt_pair(G, Gsq, Gisq, 'test_fw: normalisation')
    call matrix_mul(CL, Gisq, CLn)
    call matrix_mul(CS, Gisq, CSn)

    call construct_fw_coefficients(CLn, CSn, SLL, SSS, C1)

    ! check 2 (and 3 for a singular SSS): CFW^+ SLL CFW = 1
    allocate(CbF(np,nL,rcq_))
    call transpose_conjg(C1, CbF); call matrix_mul(CbF, SLL, C1, T)
    d_orth = maxval(abs(T - Id))

    ! check 1: agreement with the reference implementation (only meaningful if SSS is well conditioned)
    d_ref = 0.D0
    if (.not. singular) then
       call reference_fw(CLn, CSn, SLL, SSS, C2)
       d_ref = maxval(abs(C1 - C2))/maxval(abs(C2))
    end if

    ok = (d_orth < tol_orth) .and. (d_ref < tol_ref) .and. all(C1 == C1)
    write(label,'(A,I0,A,L1)') 'rcq=', rcq_, ' singular=', singular
    write(*,'(A,A,I4,I4,I4,A,ES9.2,A,ES9.2,A)') trim(label), ' nL,nS,nspinor =', nL, nS, np, &
         '  |CFW+SLL CFW - 1| = ', d_orth, '  rel.diff to reference = ', d_ref, merge('  ok    ', '  FAILED', ok)
    if (.not. ok) nfail_ = nfail_ + 1
  end subroutine

end program
