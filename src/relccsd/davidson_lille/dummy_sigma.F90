module dummy_sigmavec
! module to calculate sigma vectors for complex variables.
! it can work in two ways: with complex variables themselves, or splitting 
!
use kindvar
use mconst
use mdebug
use f95_precision, only: wp => dp
implicit none

private

public davhv

interface davhv
   module procedure cdavhv
   module procedure cdavhv_vrvi
   module procedure rdavhv
end interface davhv

contains

!
! driver routines according to type of arguments
!
! in these routines the kind of sigma vector to be constructed is selected (e.g.
! eom-ee/eom-ip/eom-ea, or different 
! 

   subroutine cdavhv(hvc,vc,n,nvec)
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer , intent(in) :: n,nvec
      complex (kind=K_C8), intent(in), dimension (:,:) :: vc
      complex (kind=K_C8), intent(inout), dimension (:,:) :: hvc
!
!     ----- local variables -----
!
      integer :: citype = 1
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
   end subroutine cdavhv

   subroutine cdavhv_vrvi(hvr,hvi,vr,vi,n,nvec)
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer , intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr,vi
      double precision, intent(inout), dimension (:,:) :: hvr,hvi
!
!     ----- local variables -----
!
      integer :: citype = 1
!
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
   end subroutine cdavhv_vrvi

   subroutine rdavhv(hv,vr,n,nvec)
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer , intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr
      double precision, intent(out), dimension (:,:) :: hv
!
!     ----- local variables -----
!
      integer :: citype = 1
!
#ifdef _MYDEBUG_
      call subent('DAVHV')
#endif
      select case (citype)
         case (1)
            call r_dummy_sigma(hv,vr,n,nvec)
      end select
#ifdef _MYDEBUG_
      call subret('DAVHV',1)
#endif
      return
   end subroutine rdavhv


!
! routines to compute sigma vectors according to type of arguments
!

   subroutine c1_dummy_sigma(hvc,vc,n,nvec)
      integer , intent(in) :: n,nvec
      complex (kind=K_C8), intent(in), dimension (:,:) :: vc
      complex (kind=K_C8), intent(inout), dimension (:,:) :: hvc

      print *, 'doing nothing inside c1_dummy_sigma'

   end subroutine c1_dummy_sigma

   subroutine c2_dummy_sigma(hvr,hvi,vr,vi,n,nvec)
      integer , intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr,vi
      double precision, intent(inout), dimension (:,:) :: hvr,hvi

      print *, 'doing nothing inside c2_dummy_sigma'

   end subroutine c2_dummy_sigma

   subroutine r_dummy_sigma(hv,vr,n,nvec)
      implicit none
      integer , intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr
      double precision, intent(out), dimension (:,:) :: hv

      print *, 'doing nothing inside r_dummy_sigma'

   end subroutine r_dummy_sigma

end module dummy_sigmavec
