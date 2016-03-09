subroutine dtrf(cv,mv,pr,np,ndet,nvec,wr,wkdim)
!
!***********************************************************************
!     ----- Expand the eigenvectors on the basis of determinants -----
!=======================================================================
!     CV = vectors in the small basis (INPUT)
!      P = vectors in the big basis (INPUT/OUTPUT)
!     NP = dimension of pr
!    NCF = dimension in the big basis (INPUT)
!   NVEC = number of vectors |v>  (INPUT)
!      W = working array (dimension nvec)
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
   use blas95, only: gemm
   implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
   integer, intent(in) :: mv,np,ndet,nvec,wkdim
   double precision, intent(in), dimension (:,:) :: cv
   double precision, intent(inout), dimension (:,:) :: pr,wr
!
!     ----- local variables -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('DTRF')
#endif
   wr(1:ndet,1:nvec) = ZERO
!  call dgemm('N','N',ndet,nvec,nvec,ONE,pr,np,cv,mv,ZERO,wr,wkdim)
   call gemm(pr(1:ndet,1:nvec),cv(1:nvec,1:nvec),wr(1:ndet,1:nvec),  &
             'N','N',ONE,ZERO)
   pr(1:ndet,1:nvec) = wr(1:ndet,1:nvec)
#ifdef _MYDEBUG_
   call subret('DTRF',1)
#endif
   return
end subroutine dtrf
