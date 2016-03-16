module fragh_mod

  use kindvar
  use mparam
#ifdef _MYDEBUG_
  use mdebug
#endif
  implicit none

   public fragh

contains

subroutine fragh(e,mvec,mvmax,nseg,nfrag)
!
!***********************************************************************
!     ----- try to segment the set of vectors in Davidson -----
!     ----- take care of degeneracies                     -----
!=======================================================================
!      E (R8) = eigenvalues (INPUT)
!   MVEC (I4) = number of eigenvectors (INPUT)
!  MVMAX (I4) = maximum dimension of the model space (INPUT)
!   NSEG (I4) = an array containing the number of vectors in
!               each segment (OUPUT)
!  NFRAG (I4) = number of segments created (OUTPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
!     ----- parameters -----
!
  double precision, parameter :: SMALL=1.0D-03
!
!     ----- arguments -----
!
  integer, intent(in) :: mvec
  integer, intent(inout) :: mvmax
  integer, intent(out) :: nfrag
  integer, intent(out), dimension (:) :: nseg
  double precision, intent(in), dimension (:) :: e
!
!     ----- local variables -----
!
  integer :: nmx,ncor,iseg,nb,i,j,ndeg,maxdeg,ntot
  double precision :: val
  integer, dimension(MAXDAV) :: ideg
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
  call subent('FRAGH')
#endif
  mvmax = min(MAXDAV,mvmax)
  nmx = mvmax
  ncor = MAXCOR
  iseg = (nmx-mvec)/ncor
  nb = iseg
!
!     ----- Special case for space of dimension 1 -----
!
  if(mvec.eq.1) then
     nfrag = 1
     nseg(nfrag) = 1
#ifdef _MYDEBUG_
     call subret('FRAGH',1)
#endif
     return
  endif
!
!     ----- look for a border between 2 sets of degenerate vectors -----
!
  ndeg = 1
  ideg(ndeg) = 1
  val = e(1)
  do i=2,mvec
     if(abs(val-e(i)).lt.SMALL) then
        ideg(ndeg) = ideg(ndeg)+1
     else
        ndeg = ndeg+1
        ideg(ndeg) = 1
        val = e(i)
     endif
  enddo
  maxdeg = 0
  do i=1,ndeg
     maxdeg = max(maxdeg,ideg(i))
  enddo
  if(iseg.lt.maxdeg) then
     ncor = (nmx-mvec)/maxdeg
     iseg = maxdeg
  endif
!
!     ---- Determine the blocks of Davidson's diagonalization ----
!
  nb = ideg(1)
  if(iseg.ge.mvec) then
     nfrag = 1
     nseg(nfrag) = mvec
  else
     nfrag = 0
     ntot = 0
     do i=2,ndeg
        nb = nb+ideg(i)
        if(nb.ge.iseg) then
           nfrag = nfrag+1
           nseg(nfrag) = nb-ideg(i)
           ntot = ntot+nseg(nfrag)
           nb = ideg(i)
        endif
     enddo
     if(ntot.lt.mvec) then
        nfrag = nfrag+1
        nseg(nfrag) = mvec-ntot
     endif
  endif
  ntot = 0
  do i=1,nfrag
     do j=1,nseg(i)
        ntot = ntot+1
     enddo
  enddo
#ifdef _MYDEBUG_
  call subret('FRAGH',2)
#endif
  return
end subroutine fragh

end module fragh_mod
