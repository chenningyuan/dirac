module davidson_schmidt

use kindvar
use mconst
!use mdebug
!use mercod
!use matpri
use f95_precision, only: wp => dp
implicit none

public dav_schmi

interface dav_schmi
   module procedure cschmi
   module procedure cschmi_vrvi
   module procedure rschmi
end interface dav_schmi

private

contains

   subroutine cschmi(vc,n,ivec,ni,stol)
!
!*******************************************************************************
!  ----- vectors |ivec+1> => |ivec+ni> are first made orthogonal -----
!  ----- to vectors |1> => |ivec> an then among themselves       -----
!===============================================================================
!     PR = array of the vectors (real part)   (R8)
!     PI = array of the vectors (imaginary part) (R8)
!     ND = leading dimension of the arrays PR et PI  (I4)
!      N = size of the basis (dimension of the problem) (I4)
!   IVEC = number of vectors already orthonormalized (I4)
!     NI = number of vectors to be orthonormalized (I4). these vectors
!          are stored from ivec+1 to ivec+ni
!          At return NI is set to the number of kept  corrections
!   STOL = convergence criterion for the correction (R8)(input)
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      use blas95, only: dotc,nrm2
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n,ivec
      integer, intent(inout) :: ni
      double precision, intent(in) :: stol
      complex (kind=K_C8), intent(inout), dimension (:,:) :: vc
!
!     ----- local variables -----
!
      integer :: i,j,nj
      double precision :: dum,gg
      complex (kind=K_C8) :: sc
!
!     +++++
!     start
!     +++++
!
!     ----- Schmidt orthogonalisation        -----
!     ----- keep the large norms (>0.1*stol) -----
!
      nj = 0
      do i=ivec+1,ivec+ni
         do j=1,ivec+nj
            sc = dotc(vc(1:n,j),vc(1:n,i))
            vc(1:n,i) = -sc*vc(1:n,j)+vc(1:n,i)
         enddo
         gg = nrm2(vc(1:n,i))
         if(gg.gt.stol*0.1D+00) then
            dum = ONE/gg
            nj = nj+1
            vc(1:n,i) = vc(1:n,i)*dum
            if(ivec+nj.ne.i) vc(1:n,ivec+nj) = vc(1:n,i)
         endif
      enddo
      ni = nj
      return
   end subroutine cschmi

   subroutine cschmi_vrvi(vr,vi,n,ivec,ni,stol)
!
!***********************************************************************
!  ----- vectors |ivec+1> => |ivec+ni> are first made orthogonal -----
!  ----- to vectors |1> => |ivec> an then among themselves       -----
!=======================================================================
!     VR = array of the vectors (real part)   (R8)
!     VI = array of the vectors (imaginary part) (R8)
!     ND = leading dimension of the arrays PR et PI  (I4)
!      N = size of the basis (dimension of the problem) (I4)
!   IVEC = number of vectors already orthonormalized (I4)
!     NI = number of vectors to be orthonormalized (I4). these vectors
!          are stored from ivec+1 to ivec+ni
!          At return NI is set to the number of kept  corrections
!   STOL = convergence criterion for the correction (R8)(input)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
     use blas95, only: dot
     implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
     integer, intent(in) :: n,ivec
     integer, intent(inout) :: ni
     double precision, intent(in) :: stol
     double precision, intent(inout), dimension (:,:) :: vr,vi
!
!     ----- local variables -----
!
     integer :: i,j,nj
     double precision :: sr,si,dum,gg
!
!     +++++
!     start
!     +++++
!
!     ----- Schmidt orthogonalisation -----
!     ----- keep the large norms (>stol)    -----
!
     nj = 0
     do i=ivec+1,ivec+ni
        do j=1,ivec+nj
           sr = dot(vr(1:n,j),vr(1:n,i))+dot(vi(1:n,j),vi(1:n,i))
           si = dot(vr(1:n,j),vi(1:n,i))-dot(vi(1:n,j),vr(1:n,i))
           vr(1:n,i) = -sr*vr(1:n,j)+vr(1:n,i)
           vi(1:n,i) = -sr*vi(1:n,j)+vi(1:n,i)
           vr(1:n,i) = si*vi(1:n,j)+vr(1:n,i)
           vi(1:n,i) = -si*vr(1:n,j)+vi(1:n,i)
        enddo
        gg = sqrt(sum(vr(1:n,i)**2+vi(1:n,i)**2))
        if(gg.gt.stol*0.1D+00) then
           dum = ONE/gg
           vr(1:n,i) = vr(1:n,i)*dum
           vi(1:n,i) = vi(1:n,i)*dum
           nj = nj+1
           if(ivec+nj.ne.i) then
              vr(1:n,ivec+nj) = vr(1:n,i)
              vi(1:n,ivec+nj) = vi(1:n,i)
           endif
        endif
     enddo
     ni = nj
     return
   end subroutine cschmi_vrvi

   subroutine rschmi(vr,n,ivec,ni,stol)
!
!*******************************************************************************
!  ----- vectors |ivec+1> => |ivec+ni> are first made orthogonal -----
!  ----- to vectors |1> => |ivec> an then among themselves       -----
!===============================================================================
!     PR = array of the vectors (real part)   (R8)
!     ND = leading dimension of the array PR  (I4)
!      N = size of the basis (dimension of the problem) (I4)
!   IVEC = number of vectors already orthonormalized (I4)
!     NI = number of vectors to be orthonormalized (I4). these vectors
!          are stored from ivec+1 to ivec+ni
!          At return NI is set to the number of kept  corrections
!   STOL = convergence criterion for the correction (R8)(input)
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      use blas95, only: dot,nrm2
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n,ivec
      integer, intent(inout) :: ni
      double precision, intent(in) :: stol
      double precision, intent(inout), dimension (:,:) :: vr
!
!     ----- local variables -----
!
      integer :: i,j,nj
      double precision :: dum,gg
      double precision :: sr
!
!     +++++
!     start
!     +++++
!
!     ----- Schmidt orthogonalisation        -----
!     ----- keep the large norms (>0.1*stol) -----
!
      nj = 0
      do i=ivec+1,ivec+ni
         do j=1,ivec+nj
            sr = dot(vr(1:n,j),vr(1:n,i))
            vr(1:n,i) = -sr*vr(1:n,j)+vr(1:n,i)
         enddo
         gg = nrm2(vr(1:n,i))
         if(gg.gt.stol*0.1D+00) then
            dum = ONE/gg
            nj = nj+1
            vr(1:n,i) = vr(1:n,i)*dum
            if(ivec+nj.ne.i) vr(1:n,ivec+nj) = vr(1:n,i)
         endif
      enddo
      ni = nj
      return
   end subroutine rschmi


end module davidson_schmidt
