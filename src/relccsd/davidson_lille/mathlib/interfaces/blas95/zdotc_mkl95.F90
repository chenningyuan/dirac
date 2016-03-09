function zdotc_mkl95(x,y,ix,iy)
!
!***********************************************************************
! fortran77 call:
! val = sdot(n,x,incx,y,incy)
!***********************************************************************
!
!     ----- use statements -----
!  use mkl77_blas, only: mkl77_sdot
   implicit none
!     ----- kind parameter -----
   integer, parameter :: WP = kind(1.0D+00)
   complex(WP) :: zdotc_mkl95
!     ----- arguments -----
   integer, intent(in), optional :: ix,iy
   complex(WP), intent(in), dimension (:) :: x,y
!     ----- local variables -----
   character(len=3), parameter :: srname='dot'
   integer :: incx,incy,n
   complex(WP) :: zdotc
!     ----- intrinsic functions -----
!  intrinsic size
!
!     +++++
!     start
!     +++++
!
   if(present(ix)) then
      incx = ix
   else
      incx = 1
   endif
   if(present(iy)) then
      incy = iy
   else
      incy = 1
   endif
   n = size(x)
   zdotc_mkl95 = zdotc(n,x,incx,y,incy)
end function zdotc_mkl95
