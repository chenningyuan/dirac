function scnrm2_mkl95(x)
!
!*******************************************************************************
! fortran77 call:
! scnrm2(n,x,incx)
!*******************************************************************************
!
!   use mkl77_blas, only: mkl77_nrm2
    implicit none
!     ----- kind parameter -----
    integer, parameter :: WP = kind(1.0E0)
    real(WP) :: scnrm2_mkl95
!     ----- arguments -----
    complex(WP), intent(in) :: x(:)
!     ----- local variables -----
    character(len=4), parameter :: srname = 'NRM2'
    integer :: incx,n
    real(WP) :: scnrm2
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
    scnrm2_mkl95 = scnrm2(n,x,incx)
end function scnrm2_mkl95
