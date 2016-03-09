subroutine zscal_mkl95(x,a)
!
!*******************************************************************************
! fortran77 call:
! dzscal(n,x,incx)
!*******************************************************************************
!
!   use mkl77_blas, only: mkl77_scal
    implicit none
!     ----- kind parameter -----
    integer, parameter :: WP = kind(1.0D0)
!     ----- arguments -----
    complex(WP), intent(in) :: x(:)
    complex(WP), intent(in) :: a
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
    call zscal(n,a,x,incx)
end subroutine zscal_mkl95
