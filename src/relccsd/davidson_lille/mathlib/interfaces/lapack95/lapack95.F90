module lapack95
!
!***********************************************************************
!   ----- this module contains interface to call lapack routines -----
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
interface getrf
   subroutine cgetrf_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(out), optional :: info
      complex(WP), intent(inout), dimension (:,:) :: a
      integer, intent(out), optional, target, dimension (:) :: ipiv
   end subroutine cgetrf_mkl95
   subroutine dgetrf_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(out), optional :: info
      real(WP), intent(inout), dimension (:,:) :: a
      integer, intent(out), optional, target, dimension (:) :: ipiv
   end subroutine dgetrf_mkl95
   subroutine sgetrf_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(out), optional :: info
      real(WP), intent(inout), dimension (:,:) :: a
      integer, intent(out), optional, target, dimension (:) :: ipiv
   end subroutine sgetrf_mkl95
   subroutine zgetrf_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(out), optional :: info
      complex(WP), intent(inout), dimension (:,:) :: a
      integer, intent(out), optional, target, dimension (:) :: ipiv
   end subroutine zgetrf_mkl95
end interface getrf

interface getri
   subroutine cgetri_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(out), optional :: info
      integer, intent(in), dimension (:) :: ipiv
      complex(WP), intent(inout), dimension (:,:) :: a
   end subroutine cgetri_mkl95
   subroutine dgetri_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(out), optional :: info
      integer, intent(in), dimension (:) :: ipiv
      real(WP), intent(inout), dimension (:,:) :: a
   end subroutine dgetri_mkl95
   subroutine sgetri_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(out), optional :: info
      integer, intent(in), dimension (:) :: ipiv
      real(WP), intent(inout), dimension (:,:) :: a
   end subroutine sgetri_mkl95
   subroutine zgetri_mkl95(a,ipiv,info)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(out), optional :: info
      integer, intent(in), dimension (:) :: ipiv
      complex(WP), intent(inout), dimension (:,:) :: a
   end subroutine zgetri_mkl95
end interface getri

interface geev
   subroutine cgeev_mkl95(a,w,vl,vr,info)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(out), optional :: info
      complex(WP), intent(out), dimension (:) :: w
      complex(WP), intent(inout), dimension (:,:) :: a
      complex(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
   end subroutine cgeev_mkl95
   subroutine dgeev_mkl95(a,wr,wi,vl,vr,info)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(out), optional :: info
      real(WP), intent(inout), dimension (:,:) :: a
      real(WP), intent(out), dimension (:) :: wr,wi
      real(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
   end subroutine dgeev_mkl95
   subroutine sgeev_mkl95(a,wr,wi,vl,vr,info)
      integer, parameter :: WP=kind(1.0E+00)
      integer, intent(out), optional :: info
      real(WP), intent(inout), dimension (:,:) :: a
      real(WP), intent(out), dimension (:) :: wr,wi
      real(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
   end subroutine sgeev_mkl95
   subroutine zgeev_mkl95(a,w,vl,vr,info)
      integer, parameter :: WP=kind(1.0D+00)
      integer, intent(out), optional :: info
      complex(WP), intent(out), dimension (:) :: w
      complex(WP), intent(inout), dimension (:,:) :: a
      complex(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
   end subroutine zgeev_mkl95
end interface geev

interface hpev
   subroutine chpev_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0E+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      complex(WP), intent(inout), dimension (:) :: a
      complex(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine chpev_mkl95
   subroutine zhpev_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0D+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      complex(WP), intent(inout), dimension (:) :: a
      complex(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine zhpev_mkl95
end interface hpev

interface hpevd
   subroutine chpevd_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0E+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      complex(WP), intent(inout), dimension (:) :: a
      complex(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine chpevd_mkl95
   subroutine zhpevd_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0D+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      complex(WP), intent(inout), dimension (:) :: a
      complex(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine zhpevd_mkl95
end interface hpevd


interface spev
   subroutine dspev_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0D+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      real(WP), intent(inout), dimension (:) :: a
      real(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine dspev_mkl95
   subroutine sspev_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0E+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      real(WP), intent(inout), dimension (:) :: a
      real(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine sspev_mkl95
end interface spev

interface spevd
   subroutine dspevd_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0D+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      real(WP), intent(inout), dimension (:) :: a
      real(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine dspevd_mkl95
   subroutine sspevd_mkl95(a,w,uplo,z,info)
      integer, parameter :: WP=kind(1.0E+00)
      character(len=1), intent(in), optional :: uplo
      integer, intent(out), optional :: info
      real(WP), intent(out), dimension (:) :: w
      real(WP), intent(inout), dimension (:) :: a
      real(WP), intent(out), optional, target, dimension (:,:) :: z
   end subroutine sspevd_mkl95
end interface spevd

interface gesvd
   subroutine cgesvd_mkl95(a,s,u,vt,ww,job,info)
       integer, parameter :: WP=kind(1.0E+00)
       character(len=1), intent(in), optional :: job
       integer, intent(out), optional :: info
       complex(WP), intent(inout), dimension (:,:) :: a
       complex(WP), intent(out), dimension (:) :: s
       complex(WP), intent(out), optional, target, dimension (:,:) :: u,vt
       complex(WP), intent(out), optional, dimension (:) :: ww
   end subroutine cgesvd_mkl95
   subroutine dgesvd_mkl95(a,s,u,vt,ww,job,info)
       integer, parameter :: WP=kind(1.0D+00)
       character(len=1), intent(in), optional :: job
       integer, intent(out), optional :: info
       real(WP), intent(inout), dimension (:,:) :: a
       real(WP), intent(out), dimension (:) :: s
       real(WP), intent(out), optional, target, dimension (:,:) :: u,vt
       real(WP), intent(out), optional, dimension (:) :: ww
   end subroutine dgesvd_mkl95
   subroutine zgesvd_mkl95(a,s,u,vt,ww,job,info)
       integer, parameter :: WP=kind(1.0D+00)
       character(len=1), intent(in), optional :: job
       integer, intent(out), optional :: info
       complex(WP), intent(inout), dimension (:,:) :: a
       complex(WP), intent(out), dimension (:) :: s
       complex(WP), intent(out), optional, target, dimension (:,:) :: u,vt
       complex(WP), intent(out), optional, dimension (:) :: ww
   end subroutine zgesvd_mkl95
end interface gesvd

#ifdef _ESSL_
interface xerbla
   subroutine xerbla_ibm(srname,info)
      character (len=*), intent(in) :: srname
      integer, intent(in) :: info
   end subroutine xerbla_ibm
end interface xerbla
#endif
end module lapack95
module mkl95_precision
    integer, parameter :: SP=kind(1.0E0)
    integer, parameter :: DP=kind(1.0D0)
end module mkl95_precision
module mkl95_lapack_inter
   interface cgeev_mkl95
      subroutine cgeev_mkl95(a,w,vl,vr,info)
         use mkl95_precision, only: wp => sp
         integer, intent(out), optional :: info
         complex(WP), intent(out), dimension (:) :: w
         complex(WP), intent(inout), dimension (:,:) :: a
         complex(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
      end subroutine cgeev_mkl95
   end interface cgeev_mkl95
   interface cgesvd_mkl95
      subroutine cgesvd_mkl95(a,s,u,vt,ww,job,info)
         use mkl95_precision, only: wp => sp
         character(len=1), intent(in), optional :: job
         integer, intent(out), optional :: info
         complex(WP), intent(inout), dimension (:,:) :: a
         complex(WP), intent(out), dimension (:) :: s
         complex(WP), intent(out), optional, target, dimension (:,:) :: u,vt
         complex(WP), intent(out), optional, dimension (:) :: ww
         end subroutine cgesvd_mkl95
   end interface cgesvd_mkl95
   interface cgetrf_mkl95
      subroutine cgetrf_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => sp
         integer, intent(out), optional :: info
         complex(WP), intent(inout), dimension (:,:) :: a
         integer, intent(out), optional, target, dimension (:) :: ipiv
      end subroutine cgetrf_mkl95
   end interface cgetrf_mkl95
   interface cgetri_mkl95
      subroutine cgetri_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => sp
         integer, intent(out), optional :: info
         integer, intent(in), dimension (:) :: ipiv
         complex(WP), intent(inout), dimension (:,:) :: a
      end subroutine cgetri_mkl95
   end interface cgetri_mkl95
   interface chpev_mkl95
      subroutine chpev_mkl95(a,w,uplo,z,info)
         use mkl95_precision, only: wp => sp
         character(len=1), intent(in), optional :: uplo
         integer, intent(out), optional :: info
         real(WP), intent(out), dimension (:) :: w
         complex(WP), intent(inout), dimension (:) :: a
         complex(WP), intent(out), optional, target, dimension (:,:) :: z
      end subroutine chpev_mkl95
   end interface chpev_mkl95
   interface dgeev_mkl95
      subroutine dgeev_mkl95(a,wr,wi,vl,vr,info)
         use mkl95_precision, only: wp => dp
         integer, intent(out), optional :: info
         real(WP), intent(inout), dimension (:,:) :: a
         real(WP), intent(out), dimension (:) :: wr,wi
         real(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
      end subroutine dgeev_mkl95
   end interface dgeev_mkl95
   interface dgesvd_mkl95
      subroutine dgesvd_mkl95(a,s,u,vt,ww,job,info)
         use mkl95_precision, only: wp => dp
         character(len=1), intent(in), optional :: job
         integer, intent(out), optional :: info
         real(WP), intent(inout), dimension (:,:) :: a
         real(WP), intent(out), dimension (:) :: s
         real(WP), intent(out), optional, target, dimension (:,:) :: u,vt
         real(WP), intent(out), optional, dimension (:) :: ww
         end subroutine dgesvd_mkl95
   end interface dgesvd_mkl95
   interface dgetrf_mkl95
      subroutine dgetrf_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => dp
         integer, intent(out), optional :: info
         real(WP), intent(inout), dimension (:,:) :: a
         integer, intent(out), optional, target, dimension (:) :: ipiv
      end subroutine dgetrf_mkl95
   end interface dgetrf_mkl95
   interface dgetri_mkl95
      subroutine dgetri_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => dp
         integer, intent(out), optional :: info
         integer, intent(in), dimension (:) :: ipiv
         real(WP), intent(inout), dimension (:,:) :: a
      end subroutine dgetri_mkl95
   end interface dgetri_mkl95
   interface dspev_mkl95
      subroutine dspev_mkl95(a,w,uplo,z,info)
         use mkl95_precision, only: wp => dp
         character(len=1), intent(in), optional :: uplo
         integer, intent(out), optional :: info
         real(WP), intent(out), dimension (:) :: w
         real(WP), intent(inout), dimension (:) :: a
         real(WP), intent(out), optional, target, dimension (:,:) :: z
      end subroutine dspev_mkl95
   end interface dspev_mkl95
   interface sgeev_mkl95
      subroutine sgeev_mkl95(a,wr,wi,vl,vr,info)
         use mkl95_precision, only: wp => sp
         integer, intent(out), optional :: info
         real(WP), intent(inout), dimension (:,:) :: a
         real(WP), intent(out), dimension (:) :: wr,wi
         real(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
      end subroutine sgeev_mkl95
   end interface sgeev_mkl95
   interface sgetrf_mkl95
      subroutine sgetrf_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => sp
         integer, intent(out), optional :: info
         real(WP), intent(inout), dimension (:,:) :: a
         integer, intent(out), optional, target, dimension (:) :: ipiv
      end subroutine sgetrf_mkl95
   end interface sgetrf_mkl95
   interface sgetri_mkl95
      subroutine sgetri_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => sp
         integer, intent(out), optional :: info
         integer, intent(in), dimension (:) :: ipiv
         real(WP), intent(inout), dimension (:,:) :: a
      end subroutine sgetri_mkl95
   end interface sgetri_mkl95
   interface sspev_mkl95
      subroutine sspev_mkl95(a,w,uplo,z,info)
         use mkl95_precision, only: wp => sp
         character(len=1), intent(in), optional :: uplo
         integer, intent(out), optional :: info
         real(WP), intent(out), dimension (:) :: w
         real(WP), intent(inout), dimension (:) :: a
         real(WP), intent(out), optional, target, dimension (:,:) :: z
      end subroutine sspev_mkl95
   end interface sspev_mkl95
   interface zgeev_mkl95
      subroutine zgeev_mkl95(a,w,vl,vr,info)
         use mkl95_precision, only: wp => dp
         integer, intent(out), optional :: info
         complex(WP), intent(out), dimension (:) :: w
         complex(WP), intent(inout), dimension (:,:) :: a
         complex(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
      end subroutine zgeev_mkl95
   end interface zgeev_mkl95
   interface zgesvd_mkl95
      subroutine zgesvd_mkl95(a,s,u,vt,ww,job,info)
         use mkl95_precision, only: wp => dp
         character(len=1), intent(in), optional :: job
         integer, intent(out), optional :: info
         complex(WP), intent(inout), dimension (:,:) :: a
         complex(WP), intent(out), dimension (:) :: s
         complex(WP), intent(out), optional, target, dimension (:,:) :: u,vt
         complex(WP), intent(out), optional, dimension (:) :: ww
         end subroutine zgesvd_mkl95
   end interface zgesvd_mkl95
   interface zgetrf_mkl95
      subroutine zgetrf_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => dp
         integer, intent(out), optional :: info
         complex(WP), intent(inout), dimension (:,:) :: a
         integer, intent(out), optional, target, dimension (:) :: ipiv
      end subroutine zgetrf_mkl95
   end interface zgetrf_mkl95
   interface zgetri_mkl95
      subroutine zgetri_mkl95(a,ipiv,info)
         use mkl95_precision, only: wp => dp
         integer, intent(out), optional :: info
         integer, intent(in), dimension (:) :: ipiv
         complex(WP), intent(inout), dimension (:,:) :: a
      end subroutine zgetri_mkl95
   end interface zgetri_mkl95
   interface zhpev_mkl95
      subroutine zhpev_mkl95(a,w,uplo,z,info)
         use mkl95_precision, only: wp => dp
         character(len=1), intent(in), optional :: uplo
         integer, intent(out), optional :: info
         real(WP), intent(out), dimension (:) :: w
         complex(WP), intent(inout), dimension (:) :: a
         complex(WP), intent(out), optional, target, dimension (:,:) :: z
      end subroutine zhpev_mkl95
   end interface zhpev_mkl95

end module mkl95_lapack_inter
