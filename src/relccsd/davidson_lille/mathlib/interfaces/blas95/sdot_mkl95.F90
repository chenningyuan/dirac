function sdot_mkl95(x,y,ix,iy)
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
   integer, parameter :: WP = kind(1.0E+000)
   real(WP) :: sdot_mkl95
!     ----- arguments -----
   integer, intent(in), optional :: ix,iy
   real(WP), intent(in), dimension (:) :: x,y
!     ----- local variables -----
   character(len=3), parameter :: srname='dot'
   integer :: incx,incy,n
   real(WP) :: sdot
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
   sdot_mkl95 = sdot(n,x,incx,y,incy)
end function sdot_mkl95
