subroutine dclvhv(vhv,hv,vr,ndet,nvec)
!
!***********************************************************************
!     ----- Fullfill the <v|H|v> matrix -----
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
   integer :: ndet,nvec
   double precision, intent(inout), dimension (:) :: vhv
   double precision, intent(in), dimension (:,:) :: hv,vr
!
!     ----- local variables -----
!
   integer :: j,ij
!  double precision :: ddot
!
!     +++++
!     start
!     +++++
!
!     ----- ATTENTION -----
!     ----- on calcule |HV> = H*|v>, mais on calcule   -----
!     ----- <v|H|v> = <HV|v>                           -----
!
!     ----- version avec parties reelles et imaginaires separees -----
!
#ifdef _MYDEBUG_
   call subent('DCLVHV')
#endif
   ij = nvec*(nvec-1)/2
   do j=1,nvec
      ij = ij+1
!        dumr = ZERO
!        do k=1,ndet
!           dumr = hv(k,nvec)*vr(k,j)+dumr
!        enddo
!        vhv(ij) = dumr
!     vhv(ij) = ddot(ndet,hv(1,nvec),1,vr(1,j),1)
      vhv(ij) = dot(hv(1:ndet,nvec),vr(1:ndet,j))
   enddo
#ifdef _MYDEBUG_
   call subret('DCLVHV',2)
#endif
   return
end subroutine dclvhv
