module blas95
!
!***********************************************************************
!   ----- this module contains interface to call blas routines -----
!   ----- in the MKL-95 style for those who use lapack or essl -----
!-----------------------------------------------------------------------
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use f95_precision
interface axpy
   subroutine caxpy_mkl95(x,y,a,ix,iy)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(in), optional :: ix,iy
      complex(WP), intent(in), optional :: a
      complex(WP), intent(in) :: x(:)
      complex(WP), intent(inout) :: y(:)
   end subroutine caxpy_mkl95
   subroutine daxpy_mkl95(x,y,a,ix,iy)
      integer, parameter :: WP = kind(1.0D+00)
      integer, intent(in), optional :: ix,iy
      real(WP), intent(in), optional :: a
      real(WP), intent(in) :: x(:)
      real(WP), intent(inout) :: y(:)
   end subroutine daxpy_mkl95
   subroutine saxpy_mkl95(x,y,a,ix,iy)
      integer, parameter :: WP = kind(1.0E+00)
      integer, intent(in), optional :: ix,iy
      real(WP), intent(in), optional :: a
      real(WP), intent(in) :: x(:)
      real(WP), intent(inout) :: y(:)
   end subroutine saxpy_mkl95
   subroutine zaxpy_mkl95(x,y,a,ix,iy)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(in), optional :: ix,iy
      complex(WP), intent(in), optional :: a
      complex(WP), intent(in) :: x(:)
      complex(WP), intent(inout) :: y(:)
   end subroutine zaxpy_mkl95
end interface axpy

interface dot
   function ddot_mkl95(x,y,ix,iy)
      integer, parameter :: WP = kind(1.0D+00)
      integer, intent(in), optional :: ix,iy
      real(WP), intent(in), dimension (:) :: x,y
      real(WP) :: ddot_mkl95
   end function ddot_mkl95
   function sdot_mkl95(x,y,ix,iy)
      integer, parameter :: WP = kind(1.0E+000)
      real(WP) :: sdot_mkl95
      integer, intent(in), optional :: ix,iy
      real(WP), intent(in), dimension (:) :: x,y
   end function sdot_mkl95
end interface dot

interface dotc
   function zdotc_mkl95(x,y,ix,iy)
      integer, parameter :: WP = kind(1.0D+00)
      complex(WP) :: zdotc_mkl95
      integer, intent(in), optional :: ix,iy
      complex(WP), intent(in), dimension (:) :: x,y
   end function zdotc_mkl95
end interface dotc

interface gemm
   subroutine cgemm_mkl95(a,b,c,transa,transb,alpha,beta)
       integer, parameter :: WP=kind(1.0E+00)
       character(len=1), intent(in), optional :: transa,transb
       complex(WP), intent(in), optional :: alpha,beta
       complex(WP), intent(in), dimension (:,:) :: a,b
       complex(WP), intent(inout), dimension (:,:) :: c
   end subroutine cgemm_mkl95
   subroutine dgemm_mkl95(a,b,c,transa,transb,alpha,beta)
       integer, parameter :: WP=kind(1.0D+00)
       character(len=1), intent(in), optional :: transa,transb
       real(WP), intent(in), optional :: alpha,beta
       real(WP), intent(in), dimension (:,:) :: a,b
       real(WP), intent(inout), dimension (:,:) :: c
   end subroutine dgemm_mkl95
   subroutine sgemm_mkl95(a,b,c,transa,transb,alpha,beta)
       integer, parameter :: WP=kind(1.0E+00)
       character(len=1), intent(in), optional :: transa,transb
       real(WP), intent(in), optional :: alpha,beta
       real(WP), intent(in), dimension (:,:) :: a,b
       real(WP), intent(inout), dimension (:,:) :: c
   end subroutine sgemm_mkl95
   subroutine zgemm_mkl95(a,b,c,transa,transb,alpha,beta)
       integer, parameter :: WP=kind(1.0D+00)
       character(len=1), intent(in), optional :: transa,transb
       complex(WP), intent(in), optional :: alpha,beta
       complex(WP), intent(in), dimension (:,:) :: a,b
       complex(WP), intent(inout), dimension (:,:) :: c
   end subroutine zgemm_mkl95
end interface gemm

interface gemv
   subroutine dgemv_mkl95(a,x,y,alpha,beta,trans)
       integer, parameter :: WP=kind(1.0D+00)
       character(len=1), intent(in), optional :: trans
       real(WP), intent(in), optional :: alpha
       real(WP), intent(in), optional :: beta
       real(WP), intent(in) :: x(:)
       real(WP), intent(in) :: a(:,:)
       real(WP), intent(inout) :: y(:)
   end subroutine dgemv_mkl95
   subroutine sgemv_mkl95(a,x,y,alpha,beta,trans)
       integer, parameter :: WP=kind(1.0E+00)
       character(len=1), intent(in), optional :: trans
       real(WP), intent(in), optional :: alpha
       real(WP), intent(in), optional :: beta
       real(WP), intent(in) :: x(:)
       real(WP), intent(in) :: a(:,:)
       real(WP), intent(inout) :: y(:)
   end subroutine sgemv_mkl95
   subroutine cgemv_mkl95(a,x,y,alpha,beta,trans)
       integer, parameter :: WP=kind(1.0E+00)
       character(len=1), intent(in), optional :: trans
       complex(WP), intent(in), optional :: alpha
       complex(WP), intent(in), optional :: beta
       complex(WP), intent(in) :: a(:,:)
       complex(WP), intent(in) :: x(:)
       complex(WP), intent(inout) :: y(:)
   end subroutine cgemv_mkl95
   subroutine zgemv_mkl95(a,x,y,alpha,beta,trans)
       integer, parameter :: WP=kind(1.0D+00)
       character(len=1), intent(in), optional :: trans
       complex(WP), intent(in), optional :: alpha
       complex(WP), intent(in), optional :: beta
       complex(WP), intent(in) :: a(:,:)
       complex(WP), intent(in) :: x(:)
       complex(WP), intent(inout) :: y(:)
   end subroutine zgemv_mkl95
end interface gemv

interface nrm2
   function dnrm2_mkl95(x)
       integer, parameter :: WP = kind(1.0D+00)
       real(WP) :: dnrm2_mkl95
       real(WP), intent(in) :: x(:)
   end function dnrm2_mkl95
   function dznrm2_mkl95(x)
       integer, parameter :: WP = kind(1.0D+00)
       real(WP) :: dznrm2_mkl95
       complex(WP), intent(in) :: x(:)
   end function dznrm2_mkl95
   function scnrm2_mkl95(x)
       integer, parameter :: WP = kind(1.0E+00)
       real(WP) :: scnrm2_mkl95
       complex(WP), intent(in) :: x(:)
   end function scnrm2_mkl95
   function snrm2_mkl95(x)
       integer, parameter :: WP = kind(1.0E+00)
       real(WP) :: snrm2_mkl95
       real(WP), intent(in) :: x(:)
   end function snrm2_mkl95
end interface nrm2

interface scal
   subroutine dscal_mkl95(x,a)
       integer, parameter :: WP = kind(1.0D+00)
       real(WP), intent(in) :: a
       real(WP), intent(inout) :: x(:)
   end subroutine dscal_mkl95
   subroutine cscal_mkl95(x,a)
       integer, parameter :: WP = kind(1.0E+00)
       complex(WP), intent(in) :: a
       complex(WP), intent(inout) :: x(:)
   end subroutine cscal_mkl95
   subroutine csscal_mkl95(x,a)
       integer, parameter :: WP = kind(1.0E+00)
       real(WP), intent(in) :: a
       complex(WP), intent(inout) :: x(:)
   end subroutine csscal_mkl95
   subroutine sscal_mkl95(x,a)
       integer, parameter :: WP = kind(1.0E+00)
       real(WP), intent(in) :: a
       real(WP), intent(inout) :: x(:)
   end subroutine sscal_mkl95
   subroutine zdscal_mkl95(x,a)
       integer, parameter :: WP = kind(1.0D+00)
       real(WP), intent(in) :: a
       complex(WP), intent(inout) :: x(:)
   end subroutine zdscal_mkl95
   subroutine zscal_mkl95(x,a)
       integer, parameter :: WP = kind(1.0D+00)
       complex(WP), intent(in) :: a
       complex(WP), intent(inout) :: x(:)
   end subroutine zscal_mkl95
end interface scal
end module blas95
module mkl95_precision
    integer, parameter :: SP=kind(1.0E0)
    integer, parameter :: DP=kind(1.0D0)
end module mkl95_precision
module mkl95_blas_inter
   interface cgemm_mkl95
      subroutine cgemm_mkl95(a,b,c,transa,transb,alpha,beta)
         use mkl95_precision, only: wp => sp
         character(len=1), intent(in), optional :: transa,transb
         complex(WP), intent(in), optional :: alpha,beta
         complex(WP), intent(in), dimension (:,:) :: a,b
         complex(WP), intent(inout), dimension (:,:) :: c
      end  subroutine cgemm_mkl95
   end interface cgemm_mkl95
   interface dgemv_mkl95
      subroutine dgemv_mkl95(a,x,y,alpha,beta,trans)
          integer, parameter :: WP=kind(1.0D+00)
          character(len=1), intent(in), optional :: trans
          real(WP), intent(in), optional :: alpha
          real(WP), intent(in), optional :: beta
          real(WP), intent(in) :: x(:)
          real(WP), intent(in) :: a(:,:)
          real(WP), intent(inout) :: y(:)
      end subroutine dgemv_mkl95
   end interface dgemv_mkl95
   interface sgemv_mkl95
      subroutine sgemv_mkl95(a,x,y,alpha,beta,trans)
          integer, parameter :: WP=kind(1.0E+00)
          character(len=1), intent(in), optional :: trans
          real(WP), intent(in), optional :: alpha
          real(WP), intent(in), optional :: beta
          real(WP), intent(in) :: x(:)
          real(WP), intent(in) :: a(:,:)
          real(WP), intent(inout) :: y(:)
      end subroutine sgemv_mkl95
   end interface sgemv_mkl95
   interface cgemv_mkl95
      subroutine cgemv_mkl95(a,x,y,alpha,beta,trans)
          integer, parameter :: WP=kind(1.0E+00)
          character(len=1), intent(in), optional :: trans
          complex(WP), intent(in), optional :: alpha
          complex(WP), intent(in), optional :: beta
          complex(WP), intent(in) :: a(:,:)
          complex(WP), intent(in) :: x(:)
          complex(WP), intent(inout) :: y(:)
      end subroutine cgemv_mkl95
   end interface cgemv_mkl95
   interface zgemv_mkl95
      subroutine zgemv_mkl95(a,x,y,alpha,beta,trans)
          integer, parameter :: WP=kind(1.0D+00)
          character(len=1), intent(in), optional :: trans
          complex(WP), intent(in), optional :: alpha
          complex(WP), intent(in), optional :: beta
          complex(WP), intent(in) :: a(:,:)
          complex(WP), intent(in) :: x(:)
          complex(WP), intent(inout) :: y(:)
      end subroutine zgemv_mkl95
   end interface zgemv_mkl95
   interface cscal_mkl95
      subroutine cscal_mkl95(x)
         use mkl95_precision, only: wp=> sp
         complex(WP), intent(inout) :: x(:)
      end subroutine cscal_mkl95
   end interface cscal_mkl95
   interface daxpy_mkl95
      subroutine daxpy_mkl95(x,y,a,ix,iy)
         use mkl95_precision, only: wp => dp
         integer, intent(in), optional :: ix,iy
         real(WP), intent(in), optional :: a
         real(WP), intent(in) :: x(:)
         real(WP), intent(inout) :: y(:)
      end subroutine daxpy_mkl95
   end interface daxpy_mkl95
   interface ddot_mkl95
      function ddot_mkl95(x,y,ix,iy)
         use mkl95_precision, only: wp => dp
         integer, intent(in), optional :: ix,iy
         real(WP), intent(in), dimension (:) :: x,y
         real(WP) :: ddot_mkl95
      end function ddot_mkl95
   end interface ddot_mkl95
   interface dgemm_mkl95
      subroutine dgemm_mkl95(a,b,c,transa,transb,alpha,beta)
         use mkl95_precision, only: wp => dp
         character(len=1), intent(in), optional :: transa,transb
         real(WP), intent(in), optional :: alpha,beta
         real(WP), intent(in), dimension (:,:) :: a,b
         real(WP), intent(inout), dimension (:,:) :: c
      end subroutine dgemm_mkl95
   end interface dgemm_mkl95
   interface dnrm2_mkl95
      function dnrm2_mkl95(x)
         use mkl95_precision, only: wp => dp
         real(WP), intent(in) :: x(:)
         real(WP) :: dnrm2_mkl95
      end function dnrm2_mkl95
   end interface dnrm2_mkl95
   interface dscal_mkl95
      subroutine dscal_mkl95(x)
         use mkl95_precision, only: wp => dp
         real(WP), intent(inout) :: x(:)
      end subroutine dscal_mkl95
   end interface dscal_mkl95
   interface dznrm2_mkl95
      function dznrm2_mkl95(x)
         use mkl95_precision, only: wp => dp
         complex(WP), intent(in) :: x(:)
         real(WP) :: dznrm2_mkl95
      end function dznrm2_mkl95
   end interface dznrm2_mkl95
   interface saxpy_mkl95
      subroutine saxpy_mkl95(x,y,a,ix,iy)
         use mkl95_precision, only: wp => sp
         integer, intent(in), optional :: ix,iy
         real(WP), intent(in), optional :: a
         real(WP), intent(in) :: x(:)
         real(WP), intent(inout) :: y(:)
      end subroutine saxpy_mkl95
   end interface saxpy_mkl95
   interface scnrm2_mkl95
      function scnrm2_mkl95(x)
         use mkl95_precision, only: wp => sp
         complex(WP), intent(in) :: x(:)
         real(WP) :: scnrm2_mkl95
      end function scnrm2_mkl95
   end interface scnrm2_mkl95
   interface sdot_mkl95
      function sdot_mkl95(x,y,ix,iy)
         use mkl95_precision, only: wp => sp
         integer, intent(in), optional :: ix,iy
         real(WP), intent(in), dimension (:) :: x,y
         real(WP) :: sdot_mkl95
      end function sdot_mkl95
   end interface sdot_mkl95
   interface sgemm_mkl95
      subroutine sgemm_mkl95(a,b,c,transa,transb,alpha,beta)
         use mkl95_precision, only: wp => sp
         character(len=1), intent(in), optional :: transa,transb
         real(WP), intent(in), optional :: alpha,beta
         real(WP), intent(in), dimension (:,:) :: a,b
         real(WP), intent(inout), dimension (:,:) :: c
      end subroutine sgemm_mkl95
   end interface sgemm_mkl95
   interface snrm2_mkl95
      function snrm2_mkl95(x)
         use mkl95_precision, only: wp => sp
         real(WP), intent(in) :: x(:)
         real(WP) :: snrm2_mkl95
      end function snrm2_mkl95
   end interface snrm2_mkl95
   interface sscal_mkl95
      subroutine sscal_mkl95(x)
         use mkl95_precision, only: wp => sp
         real(WP), intent(inout) :: x(:)
      end subroutine sscal_mkl95
   end interface sscal_mkl95
   interface zaxpy_mkl95
      subroutine zaxpy_mkl95(x,y,a,ix,iy)
         use mkl95_precision, only: wp => dp
         integer, intent(in), optional :: ix,iy
         complex(WP), intent(in), optional :: a
         complex(WP), intent(in) :: x(:)
         complex(WP), intent(inout) :: y(:)
         end subroutine zaxpy_mkl95
   end interface zaxpy_mkl95
   interface zdotc_mkl95
      function zdotc_mkl95(x,y,ix,iy)
         use mkl95_precision, only: wp => dp
         integer, intent(in), optional :: ix,iy
         complex(WP), intent(in), dimension (:) :: x,y
         complex(WP) :: zdotc_mkl95
      end function zdotc_mkl95
   end interface zdotc_mkl95
   interface zgemm_mkl95
      subroutine zgemm_mkl95(a,b,c,transa,transb,alpha,beta)
         use mkl95_precision, only: wp => dp
         character(len=1), intent(in), optional :: transa,transb
         complex(WP), intent(in), optional :: alpha,beta
         complex(WP), intent(in), dimension (:,:) :: a,b
         complex(WP), intent(inout), dimension (:,:) :: c
      end subroutine zgemm_mkl95
   end interface zgemm_mkl95
   interface zscal_mkl95
      subroutine zscal_mkl95(x)
         use mkl95_precision, only: wp => dp
         complex(WP), intent(inout) :: x(:)
      end subroutine zscal_mkl95
   end interface zscal_mkl95
end module mkl95_blas_inter
