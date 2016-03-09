subroutine sscal_mkl95(x,a)
!
!*******************************************************************************
! fortran77 call:
! sscal(n,x,incx)
!*******************************************************************************
!
!   use mkl77_blas, only: mkl77_scal
    implicit none
!     ----- kind parameter -----
    integer, parameter :: WP = kind(1.0E0)
!     ----- arguments -----
    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: a
!     ----- local variables -----
    character(len=4), parameter :: srname = 'SCAL'
    integer :: incx,n
!     ----- intrinsic functions -----
!   intrinsic size
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
    incx = 1
    n = size(x)
!     ----- call blas77 routine -----
    call sscal(n,a,x,incx)
end subroutine sscal_mkl95
