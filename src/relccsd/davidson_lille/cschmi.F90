subroutine cschmi(vr,vi,n,ivec,ni,stol)
!
!***********************************************************************
!  ----- vectors |ivec+1> => |ivec+ni> are first made orthogonal -----
!  ----- to vectors |1> => |ivec> an then among themselves       -----
!=======================================================================
!     PR = array of the vectors (real part)   (R8)
!     PI = array of the vectors (imaginary part) (R8)
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
  use kindvar
  use mparam
  use mconst
   use f95_precision, only: wp => dp
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
interface
subroutine orthon(vr,vi,n,nvec)
  integer, intent(in) :: n,nvec
  double precision, intent(in), dimension (:,:) :: vr,vi
end subroutine orthon
end interface
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
end subroutine cschmi
