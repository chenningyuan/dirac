function snrm2_mkl95(x)
!
!*******************************************************************************
! fortran77 call:
! snrm2(n,x,incx)
!*******************************************************************************
!
!   use mkl77_blas, only: mkl77_nrm2
    implicit none
!     ----- kind parameter -----
    integer, parameter :: WP = kind(1.0E0)
    real(WP) :: snrm2_mkl95
!     ----- arguments -----
    real(WP), intent(in) :: x(:)
!     ----- local variables -----
    character(len=4), parameter :: srname = 'NRM2'
    integer :: incx,n
    real(WP) :: snrm2
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
    snrm2_mkl95 = snrm2(n,x,incx)
end function snrm2_mkl95
