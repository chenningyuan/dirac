module davidson_vhv_symmetric 

use kindvar
use mconst

implicit none

private

public calvhv_symmetric

interface calvhv_symmetric
   module procedure calvhv_c
   module procedure calvhv_vrvi
   module procedure calvhv_c_fullmatrix
end interface calvhv_symmetric

contains
subroutine calvhv_vrvi(vhv,hvr,hvi,vr,vi,ndet,nvec)
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
#ifdef _MYDEBUG_
   use mdebug
#endif
   use f95_precision, only: wp => dp
   use blas95, only: dot
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
end subroutine calvhv_vrvi

   subroutine calvhv_c(vhv,hvc,vc,ndet,nvec)
!
!*******************************************************************************
!     ----- compute the <v|H|v> matrix -----
!===============================================================================
!       VHV (R8) = lower triangle <v|H|v> (OUTPUT)
!        HV (R8) = H*|v matrix (INPUT)
!         V (R8) = vectors v(ndet,ndet) (INPUT)
!        NP (I4) = dimension of hvr, hvi, vr, vi (INPUT)
!      NDET (I4) = dimension of the original space (INPUT)
!      NVEC (I4) = number of vectors |v> (INPUT)
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
     use blas95, only: dotc
     implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
     integer, intent(in) :: ndet,nvec
     complex (kind=K_C8), intent(inout), dimension (:) :: vhv
     complex (kind=K_C8), intent(in), dimension (:,:) :: hvc,vc
!
!     ----- local variables -----
!
     integer :: j,ij
!
!     ----- statement functions -----
!
!
!     +++++
!     start
!     +++++
!
!     ----- ATTENTION -----
!     ----- we have |HV> = H|v>, but we compute -----
!     ----- <v|H|w> = <HV|w> (for ZHPEV)        -----
!     ----- ZHPEV uses the upper triangle of vhv -----
!     ----- vhv(j,nvec) = vhv(j + (nvec*(nvec-1))/2)
!
     ij = nvec*(nvec-1)/2
     do j=1,nvec
        ij = ij+1
        vhv(ij) = dotc(vc(1:ndet,j),hvc(1:ndet,nvec))
     enddo
   end subroutine calvhv_c

   subroutine calvhv_c_fullmatrix(vhv_full,hvc,vc,ndet,nvec)
!
!*******************************************************************************
!     ----- compute the <v|H|v> matrix -----
!===============================================================================
!       VHV (C8) = full matrix <v|H|v> for reduced problem (OUTPUT)
!        HV (C8) = sigma vectors (H*|v matrix) (INPUT)
!         V (C8) = vectors v(ndet,ndet) (INPUT)
!        NP (I4) = dimension of hvr, hvi, vr, vi (INPUT)
!      NDET (I4) = dimension of the original space (INPUT)
!      NVEC (I4) = number of vectors |v> (INPUT)
!-------------------------------------------------------------------------------
!     written by: Andre Gomes 
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
     use blas95, only: dotc
     implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
     integer, intent(in) :: ndet,nvec
     complex (kind=K_C8), intent(inout), dimension (:,:) :: vhv_full
     complex (kind=K_C8), intent(in), dimension (:,:) :: hvc,vc
!
!     ----- local variables -----
!
     integer :: i, j
!
     do i = 1, nvec
         do j = 1, nvec
            vhv_full(i,j) = dotc(vc(1:ndet,j),hvc(1:ndet,nvec))
         enddo
     enddo

   end subroutine calvhv_c_fullmatrix

end module
