module davidson_lille
!
!*******************************************************************************
!   ----- this module contains general data for Davidson -----
!   -----  and routines to compute H*|v> in real algebra -----
!===============================================================================
!   The hamiltonian may be
!   - stored on disk
!   - a full CI (CAS)
!   - a RAS CI
!-------------------------------------------------------------------------------
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
   use kindvar
   use mconst
   use mdebug
   use mercod
   use davidson_trf
   use davidson_vhv_symmetric
   use davidson_schmidt
   use davidson_trial_vectors_utils
   use dummy_sigmavec
   USE MATPRI

   implicit none

   public davidson_solver

   interface davidson_solver
      module procedure fdeigen
   end interface

   private
!
!     ----- parameters -----
!
      integer , parameter :: local_maxcor=10,local_maxdav=400
      integer , parameter :: hbuf = 4096
!
!     ----- global variables -----
!
   integer , save :: itmax,mvm
   double precision, save :: ftol,gtol,stol
!
!     ----- allocatable variables -----
!
      integer , save, allocatable, dimension (:) :: icol,elstar,  &
                                                               mseg,indix
!  
      double precision, save, allocatable, dimension (:) :: hmr,vhv,eold,ff,  &
                                                            gg,dgg
      double precision, save, allocatable, dimension (:,:) :: work,cv,hvr
!
!    a mettre dans le module specifique de DA_DSK
!     integer , intent(in) :: hbuf,ifile  
!     private :: icol,elstar,hvr,work,indix
contains
   subroutine davini(np,nvec)
!
!*******************************************************************************
!              ----- initialize Davidson stuff -----
!-------------------------------------------------------------------------------
!     allocate the required arrays
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!
!     ----- interfaces -----
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer , intent(in) :: np,nvec
!
!     ----- local variables -----
!
      integer  :: ios
!
!     ----- allocatable variables -----
!
!     ----- formats -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('DAVINI')
#endif
!
!     ----- mandatory arrays of fdeigen -----
!
      mvm = min((local_maxcor+1)*nvec,local_maxdav)
      mvm = 426
      allocate(cv(mvm,mvm),vhv((mvm*(mvm+1))/2),eold(mvm),ff(mvm),gg(mvm),  &
               dgg(mvm),work(np,mvm),  &
               stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'DAVINI error: allocation 1 failed with ios=',ios
         call pricod(ios)
         stop
      endif
      allocate(hvr(np,mvm),mseg(mvm),stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'DAVINI error: allocation 2 failed with ios=',ios
         call pricod(ios)
         stop
      endif
!
!     ----- arrays specific to the CI method (if needed) -----
!
#ifdef _MYDEBUG_
      call subret('DAVINI',1)
#endif
   end subroutine davini
   subroutine davend
!
!*******************************************************************************
!              ----- clear Davidson stuff -----
!-------------------------------------------------------------------------------
!     deallocate the required arrays
!===============================================================================
!     written by: J.-P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
!     ----- interfaces -----
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
!     ----- local variables -----
!
!     ----- allocatable variables -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('DAVEND')
#endif
      deallocate(mseg,hvr,work,dgg,gg,ff,eold,vhv,cv)
#ifdef _MYDEBUG_
      call subret('DAVEND',1)
#endif
   end subroutine davend

   subroutine fdeigen(e,hd,vr,n,mvec,ier,vecfil)
!
!*******************************************************************************
! ----- Davidson's diagonalisation for a real symmetric matrix -----
!===============================================================================
!       PR (R8) = real part of the eigenvectors. (INPUT/OUTPUT)
!                 dimension pr(np,*)
!                 In input, pr contains the result of the guess
!       NP (I4) = physical dimension of pr (INPUT). NP >=N
!                 for efficiency , NP = N is recommanded
!        N (I4) = size of the matrix to diagonalise (INPUT)
!     MVEC (I4) = number of eigenvectors/eigenvalues to find (INPUT)
!     EOLD (R8) = array of size mvec (SCRATCH)
!       FF (R8) = convergence test on energy:
!                 FF=(e(n+1)-e(n))/e(n). Convergence is reached as
!                 soon as ff < ftol. (OUTPUT)
!       GG (R8) = convergence test on the gradient. Convergen is
!                 reached as soon as norm(gg) < gtol (OUTPUT)
!      DGG (R8) = convergence test on the direction of the correction:
!                 |dgg> = H**(-1)|G>. Convergence is reached as
!                 soon as norm(dgg) < stol (OUTPUT)
!      IER (I4) = succes index
!                 =  -1  no convergence
!                 =   1  convergence on G
!                 =  10  convergence on DGG
!                 = 100  convergence on energy
!                 if several criteria are fullfiled simultaneously,
!                 IER is the sum of the individual IER (ex: 101 = E+G)
!      HVR (R8) = working array. Size =  (np,*). (SCRATCH)
!      VHV (C16)= working array. Size at least equal to mv*(mv+1)/2
!                  (SCRATCH)
!       CV (C16)= array cv(mv,mv) (SCRATCH)
!      HMR (R8) = buffer to read the Hamitonian matrix on disk, of size
!                 HBUF (SCRATCH)
!     ICOL (I8) = buffer of size HBUF to read the indices of the non zero
!                 elements of the Hamitonian matrix saved on disk (SCRATCH)
!       HD (R8) = array containing the diagonal of H. Minimal size N (INPUT)
!        E (R8) = working array of size at least equal to mv. (SCRATCH)
!       MV (I4) = dimension of CV, maximum size of the reduced space (INPUT)
!        W (R8) = working arrays of size at least equal to 6*Mv (SCRATCH).
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      use kindvar
      use mconst
      use mdebug
!     use mdavi
!     use mgenrl
      use jpsort
      use mphase
      use mercod
      use f95_precision, only: wp => dp
      use lapack95, only: spev
      USE MATPRI
!
!     ----- interfaces -----
!
!  use dclvhv_mod
!  use savvls_mod
   use fiflsh_mod
   implicit none
!
!     ----- parameters -----
!
      double precision, parameter :: EPS=epsilon(ONE)
!
!     ----- arguments -----
!
      integer , intent(in) :: n,mvec
      integer , intent(out) :: ier
      integer , intent(in), optional :: vecfil
      double precision, intent(out), dimension (:) :: e
      double precision, intent(in), dimension (:) :: hd
      double precision, intent(inout), dimension (:,:) :: vr
!
!     ----- local variables -----
!
      character (len=1) :: reset
      integer  :: i,ivec,ni,its,info,nseg,n0,n1,n2,nfrag,ifrag,ndiag,  &
                             nprev,ncor,ier0,wkdim,np,mv
      double precision :: dumr,denom,ffm,ggm,dggm
      logical (kind=K_L4) :: permut,per
!
!     ----- allocatable variables -----
!
!     ----- formats -----
!
    9999 format(2x,a1,2x,'ITS=',i4,(t15,' IVEC=',i3,' E=',f12.8,'  DE=',d12.4,  &
                '  DG=',d12.4,'  DV=',d12.4))
    9998 format(//,5x,'*** FDEIGEN WARNING: permutation(s) of ',  &
                      'eigenvectors occurred during the iterations.',  &
                 /,5x,'*** SOME EIGENVALUES MAY BE MISSING')
    9997 format(//,5x,'ERROR in DSPEV called by FDEIGEN')
    9996 format(/,10x,'LDZ=',i4,' but N=',i4)
!
!     +++++
!     start
!     +++++
!
!     ----- initialisation:                -----
!     ----- Schmidt's orthonormalisation   -----
!     ----- initial diagonalisation        -----
!     ----- computation of H*|v>, <v|H|v>  -----
!     -----           and |g> = (H-e*I)|v> -----
!
#ifdef _MYDEBUG_
      call subent('FDEIGEN')
#endif
      if(n.eq.1) then
         e(1) = hd(1)
         vr(1,1) = ONE
         ier = 111
#ifdef _MYDEBUG_
         call subret('FDEIGEN',1)
#endif
         return
      endif
      np = size(vr,1)
      call davini(np,mvec)
      mv = mvm
!  rewind 93
!  read(93) elstar
      reset = ' '
      ier0 = 1000
      ier = 0
      ivec = 0
      ni = mvec
!     CALL DBGPRT('VC 0')
!     CALL CARMAT(VR,N,MVEC)
      call dav_schmi(vr,n,0,ni,stol)
!     CALL DBGPRT('VC 00')
!     CALL CARMAT(VR,N,MVEC)
      ndiag = mvec
      call davhv(hvr,vr,n,mvec)
      do ivec=1,mvec
         call dclvhv(vhv,hvr,vr,n,ivec)
      enddo
!     CALL DBGPRT('INIT VHV=')
!     CALL CARMAT(VHV,mvec)
      ff = ZERO
      gg = ZERO
      dgg = ZERO
      call spev(vhv(1:mvec*(mvec+1)/2),e,'U',cv,info)
      call cphase(cv,mvec,mvec)
!     CALL DBGPRT('CV=')
!     CALL CARMAT(CV,MVEC,MVEC)
      if(info.ne.0) then
         write(6,9997)
         write(6,9996) mv,mvec
         stop
      endif
      wkdim = size(work,1)
      call trf(cv,vr,n,mvec,work)
      call trf(cv,hvr,n,mvec,work)
!     DO IVEC=1,2
!        CALL DBGPRT('TRF: IVEC=',IVEC)
!        CALL DBGPRT('V')
!        DO IER=1,N
!        IF(ABS(VR(IER,IVEC)).GT.TENM10) WRITE(6,"(I15,1pd23.12)") IER+127,VR(IER,IVEC)
!        ENDDO
!     ENDDO
!     DO IVEC=1,2
!        CALL DBGPRT('TRF: IVEC=',IVEC)
!        CALL DBGPRT('HV')
!        DO IER=1,N
!        IF(ABS(HVR(IER,IVEC)).GT.TENM10) WRITE(6,"(I15,1pd23.12)") IER+127,HVR(IER,IVEC)
!        ENDDO
!     ENDDO
!
!    ------------------------------
!    ----- loop over segments -----
!    ------------------------------
!
      call find_trial_vectors_fragmentation(e,mvec,mvm,mseg,nfrag)
      vhv(1:mvec*(mvec+1)/2) = ZERO
      do ivec=1,mvec
         vhv(ivec*(ivec+1)/2) = e(ivec)
      enddo
      n0 = 0
      do ifrag=1,nfrag
         nseg = mseg(ifrag)
!
!    ----- n1, n2: limits of the current fragment              -----
!    ----- n0    : number of converged vectors in the previous -----
!    -----         fragments
!    ----- ndiag: size of the diagonalisation for the fragment -----
!    -----        at iteration its                             -----
!    ----- nprev: size of the diagonalisation for the fragment -----
!    -----        at the previous iteration (its-1)            -----
!    ----- ncor : nb of corrections for the current            -----
!
         n1 = n0+1
         n2 = n0+nseg
         ndiag = mvec
!
!     ----- vr(ni) contains (H-e(ivec)*I)|ivec> -----
!
         nprev = ndiag
         ni = ndiag
         do ivec=n1,n2
            ni = ni+1
            ggm = ZERO
            do i=1,n
               dumr = hvr(i,ivec)-e(ivec)*vr(i,ivec)
               denom = ONE/(hd(i)-e(ivec)+EPS)
               vr(i,ni) = dumr*denom
               ggm = ggm+dumr**2
            enddo
            ggm = sqrt(ggm)
            if(ggm.lt.gtol) then
               ni = ni-1
            endif
         enddo
         ncor = ni-ndiag
!
!     ----------------------------
!     ----- start iterations -----
!     ----------------------------
!
         permut = .TRUE.
         do its=1,itmax
            call dav_schmi(vr,n,nprev,ncor,stol)
            ndiag = nprev+ncor
!           if(ncor.eq.0) exit
            if(ncor.gt.0) call davhv(hvr(1:n,nprev+1:ndiag),  &
                                     vr(1:n,nprev+1:ndiag),n,ncor)
            ni = nprev
            do ivec=1,ncor
               ni = nprev+ivec
               call dclvhv(vhv,hvr,vr(1:n,1:ni),n,ni)
            enddo
            eold(1:nprev) = e(1:nprev)
            call spev(vhv(1:ndiag*(ndiag+1)/2),e,'U',cv,info)
            if(info.ne.0) then
               write(6,9997)
               write(6,9996) mv,ndiag
               stop
            endif
            per = check_trial_vectors(cv,e,ndiag,mvec)
            call trf(cv,vr,n,ndiag,work)
            call trf(cv,hvr,n,ndiag,work)
            if(ndiag+nseg.gt.mvm) then
               reset = '*'
               ndiag = mvec
            endif
            nprev = ndiag
            vhv(1:ndiag*(ndiag+1)/2) = ZERO
            do i=1,ndiag
               vhv(i*(i+1)/2) = e(i)
            enddo
!
!     ----- compute gradient and new direction -----
!
            ggm = ZERO
            dggm = ZERO
            ffm = ZERO
            ni = ndiag
            do ivec=n1,n2
               gg(ivec) = ZERO
               dgg(ivec) = ZERO
               ni = ni+1
               do i=1,n
                  dumr = hvr(i,ivec)-e(ivec)*vr(i,ivec)
                  denom = ONE/(hd(i)-e(ivec)+EPS)
                  vr(i,ni) = dumr*denom
                  gg(ivec) = dumr**2+gg(ivec)
                  dgg(ivec) = vr(i,ni)**2+dgg(ivec)
               enddo
               gg(ivec) = sqrt(gg(ivec))
               ggm = max(ggm,gg(ivec))
               dgg(ivec) = sqrt(dgg(ivec))
               dggm = max(dggm,dgg(ivec))
               ff(ivec) = abs((e(ivec)-eold(ivec))/e(ivec))
               ffm = max(ffm,ff(ivec))
               if(ff(ivec).lt.ftol.and.gg(ivec).lt.gtol) ni = ni-1
            enddo
            ncor = ni-ndiag
            write(6,9999) reset,its,(ivec,e(ivec),ff(ivec),gg(ivec),dgg(ivec),  &
                           ivec=n1,n2)
            reset = ' '
            call fiflsh(6)
!
!     ----- if converged, return -----
!
            ier = 0
            if(ggm.le.gtol) ier = 1
            if(dggm.le.stol) ier = 10+ier
            if(ffm.le.ftol) ier = 100+ier
            if(ier.gt.100) then
               if(permut.and.per) then
                  write(6,9998)
                  permut = .FALSE.
               endif
               ier0 = min(ier0,ier)
!              goto 200
               exit
            else
               ier = -1
            endif
!
!     ----- Save eigenvectors on vecfil file for subsquent restart -----
!
            if(present(vecfil).and.mod(its,5).eq.0) then
!              call savvls(vr,e,n,mvec,vecfil)
               write(6,*) 'dans fdeigen pas encore de suavegarde'
            endif
         enddo
!
!               --------------------------
!               ----- end iterations -----
!               --------------------------
!
!    200 continue
         n0 = n0+nseg
         if(ifrag.lt.nfrag) then
!           call dclchv(hvr,vr,n,mvec,hmr,icol,hbuf,elfile,elstar,work)
            call davhv(hvr(1:n,1:mvec),vr(1:n,1:mvec),n,mvec)
            do ivec=1,mvec
               call dclvhv(vhv,hvr,vr,n,ivec)
            enddo
            call spev(vhv(1:mvec*(mvec+1)/2),e,'U',cv,info)
            if(info.ne.0) then
               write(6,9997)
               write(6,9996) mv,mvec
               stop
            endif
            per = check_trial_vectors(cv,e,mvec,mvec)
            call trf(cv,vr,n,mvec,work)
            call trf(cv,hvr,n,mvec,work)
            vhv(1:mvec*(mvec+1)/2) = ZERO
            do i=1,mvec
               vhv(i*(i+1)/2) = e(i)
            enddo
         endif
      enddo
!          ----- end of loop over slices of energy -----
!
!     ----- sort final eigenvectors in increasing eigenvalue order -----
!
      allocate(indix(mvec))
      indix(1:mvec) = (/ (i,i=1,mvec) /)
      call gensort(mvec,e,indix)
      hvr(1:n,1:mvec) = vr(1:n,1:mvec)
      vr(1:n,1:mvec) = hvr(1:n,indix(1:mvec))
      ier = ier0
      deallocate(indix)
      call davend
#ifdef _MYDEBUG_
      call subret('FDEIGEN',2)
#endif
      return
   end subroutine fdeigen
end module davidson_lille
