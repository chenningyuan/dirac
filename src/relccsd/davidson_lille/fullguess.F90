subroutine fullguess(pr,e,isel,n,msize,mvec,hbuf,elfile)
!
!***********************************************************************
! ----- Generate a guess by diagonalization within the space of -----
! ----- degenerate configurations                               -----
!=======================================================================
!       PR (R8) = real part of the eigenvectors. (OUTPUT)
!                 dimension pr(np,*)
!                 In input, pr contains the result of the guess
!        E (R8) = eigenvalues (OUTPUT)
!     ISEL (I4) = indix of the selected determinants (INPUT)
!                 isel(i) = 0 if det. |i> not selected
!                 isel(i) = j if det. |i> is the jth in the selected space
!        N (I4) = size of the full determinant space (INPUT)
!    MSIZE (I4) = size of the space to diagonalize (selected determinants)
!                 (INPUT)
!     MVEC (I4) = number of eigenvectors/eigenvalues to keep (INPUT)
!     HBUF (I4) = dimension of the buffers (INPUT)
!   ELFILE (I4) = unit number of the non relativistic part of H (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use kindvar
   use mconst
   use mdebug
   use mercod
   use jpsort
   use f95_precision, only: wp => dp
!  use lapack95, only: spev
   use lapack95, only: spevd
   USE MATPRI
!
!     ----- interfaces -----
!
   use hred_mod
   implicit none
!
!     ----- parameters -----
!
   double precision, parameter :: EPS=1.0D-20
!
!     ----- arguments -----
!
   integer, intent(in) :: n,msize,mvec,hbuf,elfile
   integer, intent(in), dimension(:) :: isel
   double precision, intent(out), dimension (:) :: e
   double precision, intent(out), dimension (:,:) :: pr
!
!     ----- local variables -----
!
   integer :: mvm,i,info,ios,nx
!
!     ----- allocatable variables -----
!
   integer, allocatable, dimension (:) :: icol,elstar
   double precision, allocatable, dimension (:) :: hmr,etmp,vhv
   double precision, allocatable, dimension (:,:) :: cv
!
!     ----- formats -----
!
  9999 format(/,5x,'The lowest electronic configurations generate more ',  &
              'states than required:',/,5x,  &
              'the guess is generate for the complete configuration space ',  &
              'with',i8,' determinants.',/,5x,'The lowest ',i4,  &
              ' states are selected automatically.')
  9998 format(/,5x,'FULLGUESS error: allocation failed for HMR HMI ICOL')
  9997 format(/,5x,'FULLGUESS error: allocation failed for VHV and CV')
  9996 format(//,5x,'ERROR in DSPEV called by FULLGUESS; info = ', i4)
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
   call subent('FULGUES')
#endif
   write(6,9999) msize,mvec
   mvm = msize
   nx = (mvm*(mvm+1))/2
   allocate(hmr(hbuf),icol(hbuf),elstar(n+1),stat=ios)
   if(ios.ne.0) then
      write(6,9998)
      call pricod(ios)
      stop
   endif
   allocate(etmp(mvm),cv(mvm,mvm),vhv(nx),stat=ios)
   if(ios.ne.0) then
      write(6,9997)
      call pricod(ios)
      stop
   endif
!#ifndef _ESSL_
!   allocate(w(3*mvm),stat=ios)
!   if(ios.ne.0) then
!      write(6,9997)
!      call pricod(ios)
!      stop
!   endif
!#else
!   allocate(w(2*mvm),stat=ios)
!   if(ios.ne.0) then
!      write(6,9997)
!      call pricod(ios)
!      stop
!   endif
!#endif
!
!     ----- Get the relevant Hamiltonian submatrix on the space -----
!     ----- of selected determinants                            -----
!
   rewind(93)
   read(93) elstar
   call hred(vhv,hmr,isel,icol,n,hbuf,elfile,elstar)
!
!     ----- Diagonalize the relevant submatrix -----
!
!   call spev(vhv(1:nx),etmp,'U',cv,info)
   call spevd(vhv(1:nx),etmp,'U',cv,info)
   if(info.ne.0) then
      write(6,9996) info
      stop
   endif
!
!     ----- Construct the eigenvectors on the full space -----
!
   e(1:mvec) = etmp(1:mvec)
   do i=1,n
      if(isel(i).ne.0) then
         pr(i,1:mvec) = cv(isel(i),1:mvec)
      endif
   enddo
!  deallocate(w)
   deallocate(hmr,icol,etmp,cv,vhv,elstar)
#ifdef _MYDEBUG_
   call subret('FULGUES',1)
#endif
   return
end subroutine fullguess
