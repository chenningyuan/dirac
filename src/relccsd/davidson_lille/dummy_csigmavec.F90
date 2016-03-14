module dummy_csigmavec
! module to calculate sigma vectors for complex variables.
! it can work in two ways: with complex variables themselves, or splitting 
!
use kindvar
use mconst
use mdebug
!use mercod
!use matpri
use f95_precision, only: wp => dp
implicit none

private

public cdavhv


interface cdavhv
   module procedure cdavhv1,cdavhv2
end interface cdavhv

contains

   subroutine cdavhv1(hvc,vc,n,nvec)
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: n,nvec
      complex (kind=K_C8), intent(in), dimension (:,:) :: vc
      complex (kind=K_C8), intent(inout), dimension (:,:) :: hvc
!
!     ----- local variables -----
!
      integer :: citype = 1
!
!     ----- allocatable variables -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('CDAVHV1')
#endif
      select case (citype)
         case (1)
            print *, 'dummy call going thru CDAVHV1'
            call c1_dummy_sigma(hvc,vc,n,nvec)
      end select
#ifdef _MYDEBUG_
      call subret('CDAVHV1',1)
#endif
      return
   end subroutine cdavhv1
   subroutine cdavhv2(hvr,hvi,vr,vi,n,nvec)
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr,vi
      double precision, intent(inout), dimension (:,:) :: hvr,hvi
!
!     ----- local variables -----
!
      integer :: citype = 1
!
!     ----- allocatable variables -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('DCLCHV2')
#endif
      select case (citype)
         case (1)
            call c2_cashv(hvr,hvi,vr,vi,n,nvec)
      end select
#ifdef _MYDEBUG_
      call subret('DCLCHV2',1)
#endif
      return
   end subroutine cdavhv2

   subroutine c1_dummy_sigma(hvc,vc,n,nvec)
      integer (kind=K_I4), intent(in) :: n,nvec
      complex (kind=K_C8), intent(in), dimension (:,:) :: vc
      complex (kind=K_C8), intent(inout), dimension (:,:) :: hvc

      print *, 'doing nothing inside c1_dummy_sigma'

   end subroutine c1_dummy_sigma

   subroutine c2_dummy_sigma(hvr,hvi,vr,vi,n,nvec)
      integer (kind=K_I4), intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr,vi
      double precision, intent(inout), dimension (:,:) :: hvr,hvi

      print *, 'doing nothing inside c2_dummy_sigma'

   end subroutine c2_dummy_sigma

end module dummy_csigmavec
