subroutine trf(cv,vr,vi,ndet,nvec,win,wout,wkdim)
!
!***********************************************************************
!     ----- Develop the eigenvectors on the basis of determinants -----
!=======================================================================
!     CV (R8) = vectors in the small basis (INPUT)
!  VR,VI (R8) = vectors in the large basis (INPUT/OUTPUT)
!   NDET (i4) = dimension of the large basis (INPUT)
!   NVEC (I4) = number of eigenvectors (INPUT)
!    WIN (R8) = working array INPUT of size HBUF (SCRATCH)
!   WOUT (R8) = working array INPUT of size HBUF (SCRATCH)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
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
   integer, intent(in) :: ndet,nvec,wkdim
   double precision, intent(inout), dimension (:,:) :: vr,vi
   complex(kind=K_C8), intent(in), dimension (:,:) :: cv
   complex(kind=K_C8), intent(inout), dimension (:,:) :: win,wout
!
!     ----- local variables -----
!
   integer :: i1,i2,m
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('TRF')
#endif
   do i1=1,ndet,wkdim
      i2 = min(ndet,i1+wkdim-1)
      m = i2-i1+1
      win(1:m,1:nvec) = dcmplx(vr(i1:i2,1:nvec),vi(i1:i2,1:nvec))
      call gemm(win(1:m,1:nvec),cv(1:nvec,1:nvec),wout(1:m,1:nvec),  &
                       'N','N',CONE,CZERO)
!     call gemm('N','N',m,nvec,nvec,CONE,win,wkdim,cv,mv,CZERO,wout,wkdim)
      vr(i1:i2,1:nvec) = dble(wout(1:m,1:nvec))
      vi(i1:i2,1:nvec) = dimag(wout(1:m,1:nvec))
   enddo
#ifdef _MYDEBUG_
   call subret('TRF',1)
#endif
   return
end subroutine trf
