module dummy_sigmavec
!
!*******************************************************************************
!   ----- this module contains general data for Davidson -----
!   -----  and routines to compute H*|v> in real algebra -----
!===============================================================================
!   The hamiltonian may be
!   - stored on disk
!   - a full CI (CAS)
!   - a RAS CI
!-------------------------------------------------------------------------------
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
   use kindvar
   use mconst
   use mdebug
   use mercod
   implicit none

   private

   public davhv

contains

   subroutine davhv(hv,vr,n,nvec)
!
!*******************************************************************************
!              ----- driver routine  to compute H*|V> -----
!-------------------------------------------------------------------------------
!     call the right routine according to the value of citype/cimeth 
!     in module mcityp
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!
!     ----- interfaces -----
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer (kind=K_I4), intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr
      double precision, intent(out), dimension (:,:) :: hv
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
      call subent('DAVHV')
#endif
      select case (citype)
         case (1)
            call c_dummy_sigma(hv,vr,n,nvec)
      end select
#ifdef _MYDEBUG_
      call subret('DAVHV',1)
#endif
      return
   end subroutine

   subroutine c_dummy_sigma(hv,vr,n,nvec)
      implicit none
      integer (kind=K_I4), intent(in) :: n,nvec
      double precision, intent(in), dimension (:,:) :: vr
      double precision, intent(out), dimension (:,:) :: hv
   
      print *, 'doing nothing inside c_dummy_sigma'
   
   end subroutine c_dummy_sigma

end module dummy_sigmavec 
