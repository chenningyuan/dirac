subroutine dschmi(pr,n,ivec,ni,stol)
!
!***********************************************************************
!     ----- Schmidt orthogonalisation of the vectors ----
!=======================================================================
!     PR (R8) = real part of the vectors (INPUT/OUTPUT)
!     ND (I4) = row dimension of pr (INPUT)
!      N (I4) = dimension of the basis (INPUT)
!   IVEC (I4) = number of vectors already orthonormalised (INPUT)
!     NI (I4) = number of vectors to be orthonormalised. These
!               vectors placed from ivec+1 to ivec+ni (INPUT)
!     GG (R8) = norm of the vector (OUTUT)
!   STOL (R8) = convergence criterion on the gradient (INPUT)
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
  use mdebug
  use f95_precision, only: wp => dp
  use blas95, only: dot,axpy,nrm2
  implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
  integer, intent(in) :: n,ivec
  integer, intent(inout) :: ni
  double precision, intent(in) :: stol
  double precision, intent(inout), dimension (:,:) :: pr
!
!     ----- local variables -----
!
  integer :: i,j,nj,K
  double precision :: sij,dum,gg
  save sij
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
  call subent('DSCHMI')
#endif
  nj = 0
  do i=ivec+1,ivec+ni
     do j=1,ivec+nj
!       sij = ddot(n,pr(1,j),1,pr(1,i),1)
!       sij = dot(pr(1:n,j),pr(1:n,i))
        sij = ZERO
        do K=1,n
           sij = sij+pr(k,j)*pr(k,i)
        enddo
!       call daxpy(n,-sij,pr(1,j),1,pr(1,i),1)
        call axpy(pr(1:n,j),pr(1:n,i),-sij)
     enddo
!    gg = dnrm2(n,pr(1,i),1)
     gg = nrm2(pr(1:n,i))
     if(gg.gt.stol*0.1D+00) then
        dum = ONE/gg
        pr(1:n,i) = pr(1:n,i)*dum
        nj = nj+1
        if(ivec+nj.ne.i) then
           pr(1:n,ivec+nj) = pr(1:n,i)
        endif
     endif
  enddo
  ni = nj
#ifdef _MYDEBUG_
  call subret('DSCHMI',1)
#endif
  return
end subroutine dschmi
