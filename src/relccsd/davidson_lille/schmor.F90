subroutine schmor(vr,vi,n,ivec,ni,gg,stol)
!
!***********************************************************************
!  ----- vectors |ivec+1> => |ivec+ni> are first made orthogonal -----
!  ----- to vectors |1> => |ivec> an then among themselves       -----
!=======================================================================
!     VR = array of the vectors (real part)   (R8)
!     VI = array of the vectors (imaginary part) (R8)
!      N = size of the basis (dimension of the problem) (I4)
!   IVEC = number of vectors already orthonormalized (I4)
!     NI = number of vectors to be orthonormalized (I4). these vectors
!          are stored from ivec+1 to ivec+ni
!          At return NI is set to the number of kept  corrections
!     GG = norm of the vector
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
   use mconst
#ifdef _MYDEBUG_
   use mdebug
#endif
   use matpri
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
   double precision, intent(out) :: gg
   double precision, intent(inout), dimension (:,:) :: vr,vi
!
!     ----- local variables -----
!
   integer :: i,j,nj,ncor
   double precision :: sr,si,dum
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('SCHMOR')
#endif
!
!     ----- test the norm of the correction -----
!     ----- keep the large norms (>stol)    -----
!
   nj = 0
   do i=ivec+1,ivec+ni
      gg = sqrt(sum(vr(1:n,i)**2+vi(1:n,i)**2))
      if(gg.gt.stol) then
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
!
!     ----- Schmidt orthogonalisation -----
!     ----- keep the large norms (>stol)    -----
!
   ncor = nj
   nj = 0
   do i=ivec+1,ivec+ncor
      do j=1,i-1
!        sr = ddot(n,vr(1,j),1,vr(1,i),1)+ddot(n,vi(1,j),1,vi(1,i),1)
!        si = ddot(n,vr(1,j),1,vi(1,i),1)-ddot(n,vi(1,j),1,vr(1,i),1)
         sr = dot(vr(1:n,j),vr(1:n,i))+dot(vi(1:n,j),vi(1:n,i))
         si = dot(vr(1:n,j),vi(1:n,i))-dot(vi(1:n,j),vr(1:n,i))
!           if(abs(sr).gt.SMALL) then
         vr(1:n,i) = -sr*vr(1:n,j)+vr(1:n,i)
         vi(1:n,i) = -sr*vi(1:n,j)+vi(1:n,i)
!           endif
!           if(abs(si).gt.SMALL) then
         vr(1:n,i) = si*vi(1:n,j)+vr(1:n,i)
         vi(1:n,i) = -si*vr(1:n,j)+vi(1:n,i)
!           endif
      enddo
      gg = sqrt(sum(vr(1:n,i)**2+vi(1:n,i)**2))
      if(gg.gt.stol) then
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
#ifdef _MYDEBUG_
   call subret('SCHMOR',1)
#endif
   return
end subroutine schmor
