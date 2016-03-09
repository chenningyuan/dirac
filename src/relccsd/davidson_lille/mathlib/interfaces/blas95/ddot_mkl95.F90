function ddot_mkl95(x,y,ix,iy)
!
!***********************************************************************
! fortran77 call:
! val = ddot(n,x,incx,y,incy)
!***********************************************************************
!
!     ----- use statements -----
!  use mkl77_blas, only: mkl77_ddot
   implicit none
!     ----- kind parameter -----
   integer, parameter :: WP = kind(1.0D+00)
!     ----- arguments -----
   integer, intent(in), optional :: ix,iy
   real(WP), intent(in), dimension (:) :: x,y
   real(WP) :: ddot_mkl95
!     ----- local variables -----
   character(len=3), parameter :: srname='dot'
   integer :: incx,incy,n
   real(WP) :: ddot
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
   ddot_mkl95 = ddot(n,x,incx,y,incy)
end function ddot_mkl95
