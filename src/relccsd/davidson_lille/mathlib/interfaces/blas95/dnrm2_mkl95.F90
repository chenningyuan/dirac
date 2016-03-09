function dnrm2_mkl95(x)
!
!*******************************************************************************
! fortran77 call:
! dnrm2(n,x,incx)
!*******************************************************************************
!
!   use mkl77_blas, only: mkl77_nrm2
    implicit none
!     ----- kind parameter -----
    integer, parameter :: WP = kind(1.0D0)
    real(WP) :: dnrm2_mkl95
!     ----- arguments -----
    real(WP), intent(in) :: x(:)
!     ----- local variables -----
    character(len=4), parameter :: srname = 'NRM2'
    integer :: incx,n
    real(WP) :: dnrm2
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
    dnrm2_mkl95 = dnrm2(n,x,incx)
end function dnrm2_mkl95
