subroutine calvhv(vhv,hvr,hvi,vr,vi,ndet,nvec)
!
!***********************************************************************
!     ----- compute the <v|H|v> matrix -----
!=======================================================================
!       VHV (R8) = lower triangle <v|H|v> (OUTPUT)
!        HV (R8) = H*|v matrix (INPUT)
!         V (R8) = vectors v(ndet,ndet) (INPUT)
!        NP (I4) = dimension of hvr, hvi, vr, vi (INPUT)
!      NDET (I4) = dimension of the original space (INPUT)
!      NVEC (I4) = number of vectors |v> (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use kindvar
   use mconst
#ifdef _MYDEBUG_
   use mdebug
#endif
   use f95_precision, only: wp => dp
   use blas95, only: dot
   implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
   integer, intent(in) :: ndet,nvec
   complex (kind=K_C8), intent(out), dimension (:) :: vhv
   double precision, intent(in), dimension (:,:) :: hvr,hvi,vr,vi
!
!     ----- local variables -----
!
   integer :: j,ij
   double precision :: dumr,dumi
!
!     +++++
!     start
!     +++++
!
!     ----- ATTENTION -----
!     ----- we have |HV> = H*|v>, but we compute -----
!     ----- <v|H|w> = <HV|w> (for ZHPEV)         -----
!     ----- ZHPEV uses the upper triangle of vhv -----
!     ----- vhv(j,nvec) = vhv(j + (nvec*(nvec-1))/2)
!
#ifdef _MYDEBUG_
   call subent('CALVHV')
#endif
   ij = nvec*(nvec-1)/2
   do j=1,nvec
      ij = ij+1
      dumr = dot(vr(1:ndet,j),hvr(1:ndet,nvec))+  &
             dot(vi(1:ndet,j),hvi(1:ndet,nvec))
      dumi = dot(vr(1:ndet,j),hvi(1:ndet,nvec))-  &
             dot(vi(1:ndet,j),hvr(1:ndet,nvec))
      vhv(ij) = dcmplx(dumr,dumi)
   enddo
#ifdef _MYDEBUG_
   call subret('CALVHV',1)
#endif
   return
end subroutine calvhv
