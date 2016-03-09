subroutine dgemm_mkl95(a,b,c,transa,transb,alpha,beta)
!
!***********************************************************************
! fortran77 call:
! call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!***********************************************************************
!
!     ----- use statements -----
!   use mkl77_blas, only: mkl77_gemm
!     ----- implicit statement -----
    implicit none
!     ----- kind parameter -----
    integer, parameter :: WP=kind(1.0D+00)
!     ----- arguments -----
    character(len=1), intent(in), optional :: transa,transb
    real(WP), intent(in), optional :: alpha,beta
    real(WP), intent(in), dimension (:,:) :: a,b
    real(WP), intent(inout), dimension (:,:) :: c
!     ----- local variables -----
    character(len=4), parameter :: srname='gemm'
    character(len=1) :: o_transa
    character(len=1) :: o_transb
    integer :: m,n,k,lda,ldb,ldc
    real(WP) :: o_alpha,o_beta
!     ----- intrinsic functions -----
!   intrinsic max, present, size
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
    if(present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1
    endif
    if(present(beta)) then
        o_beta = beta
    else
        o_beta = 1
    endif
    if(present(transa)) then
        o_transa = transa
    else
        o_transa = 'n'
    endif
    if(present(transb)) then
        o_transb = transb
    else
        o_transb = 'n'
    endif
    if((o_transa.eq.'N'.or.o_transa.eq.'n')) then
        k = size(a,2)
    else
        k = size(a,1)
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call dgemm(o_transa,o_transb,m,n,k,o_alpha,a,lda,b,ldb,o_beta,c,ldc)
end subroutine dgemm_mkl95
