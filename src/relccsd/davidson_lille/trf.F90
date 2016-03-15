module davidson_trf

use kindvar
use mconst
use blas95, only: gemm

implicit none

private

public trf

interface trf
   module procedure ctrf_vrvi
   module procedure ctrf
end interface 

contains

   subroutine ctrf_vrvi(cv,vr,vi,ndet,nvec,win,wout,wkdim)
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
      use f95_precision, only: wp => dp
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
      do i1=1,ndet,wkdim
         i2 = min(ndet,i1+wkdim-1)
         m = i2-i1+1
         win(1:m,1:nvec) = dcmplx(vr(i1:i2,1:nvec),vi(i1:i2,1:nvec))
         call gemm(win(1:m,1:nvec),cv(1:nvec,1:nvec),wout(1:m,1:nvec),  &
                          'N','N',CONE,CZERO)
!        call gemm('N','N',m,nvec,nvec,CONE,win,wkdim,cv,mv,CZERO,wout,wkdim)
         vr(i1:i2,1:nvec) = dble(wout(1:m,1:nvec))
         vi(i1:i2,1:nvec) = dimag(wout(1:m,1:nvec))
      enddo
   end subroutine ctrf_vrvi

   subroutine ctrf(c,v,ndet,nvec,wc)
!
!*******************************************************************************
!     ----- Develop the eigenvectors on the basis of determinants -----
!===============================================================================
!      C (R8) = vectors in the small basis (INPUT)
!      V (R8) = vectors in the large basis (INPUT/OUTPUT)
!   NDET (i4) = dimension of the large basis (INPUT)
!   NVEC (I4) = number of eigenvectors (INPUT)
!     WC (R8) = working array of same shape as VC (SCRATCH)
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer , intent(in) :: ndet,nvec
      complex (kind=K_C8), intent(in), dimension (:,:) :: c
      complex (kind=K_C8), intent(inout), dimension (:,:) :: v,wc
!
      call gemm(v(1:ndet,1:nvec),c(1:nvec,1:nvec),wc(1:ndet,1:nvec),  &
                   'N','N',CONE,CZERO)
      v(1:ndet,1:nvec) = wc(1:ndet,1:nvec)
   end subroutine ctrf

end module davidson_trf

