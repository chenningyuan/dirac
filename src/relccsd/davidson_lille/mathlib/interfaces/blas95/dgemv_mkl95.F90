subroutine dgemv_mkl95(a,x,y,alpha,beta,trans)
!
!***********************************************************************
! fortran77 call:
! call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
!***********************************************************************
!
    implicit none
!     ----- parameters -----
    integer, parameter :: WP=kind(1.0D+00)
    character(len=4), parameter :: srname = 'gemv'
!     ----- arguments -----
    character(len=1), intent(in), optional :: trans
    real(WP), intent(in), optional :: alpha
    real(WP), intent(in), optional :: beta
    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: a(:,:)
    real(WP), intent(inout) :: y(:)
!     ----- local variables -----
    character(len=1) :: o_trans
    integer :: incx,incy,m,n,lda
    real(WP) :: o_alpha,o_beta
!
    intrinsic max, present, size
!
!     -----
!     start
!     -----
!
    if(present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1
    endif
    if(present(beta)) then
        o_beta = beta
    else
        o_beta = 0
    endif
    if(present(trans)) then
        o_trans = trans
    else
        o_trans = 'n'
    endif
    incx = 1
    incy = 1
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    call dgemv(o_trans,m,n,o_alpha,a,lda,x,incx,o_beta,y,incy)
end subroutine dgemv_mkl95
