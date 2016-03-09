subroutine zaxpy_mkl95(x,y,a,ix,iy)
!
!***********************************************************************
! fortran77 call:
! zaxpy(n,a,x,incx,y,incy)
!***********************************************************************
!
!     ----- use statements -----
!  use mkl77_blas, only: mkl77_axpy
!     ----- implicit statement -----
   implicit none
!     ----- kind parameter -----
   integer, parameter :: WP=kind(1.0D+00)
!     ----- arguments -----
   integer, intent(in), optional :: ix,iy
   complex(WP), intent(in), optional :: a
   complex(WP), intent(in) :: x(:)
   complex(WP), intent(inout) :: y(:)
!
!     ----- local variables -----
   character(len=4), parameter :: srname='axpy'
   complex(WP) :: o_a
   integer :: incx,incy,n
!     ----- intrinsic functions -----
!  intrinsic present, size
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
   if(present(a)) then
       o_a = a
   else
       o_a = dcmplx(1.0D+00,0.0D+00)
   endif
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
   call zaxpy(n,o_a,x,incx,y,incy)
end subroutine zaxpy_mkl95
