module cdavidson_lille
!*******************************************************************************
!   ----- this module contains general data for Davidson    -----
!   -----  and routines to compute H*|v> in complex algebra -----
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
   use matpri
   use davidson_vhv_symmetric
   use davidson_schmidt
   use dummy_sigmavec
   use f95_precision, only: wp => dp
   implicit none

!  private
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
                                                            sostar,mseg,indix
!  
   double precision, save, allocatable, dimension (:) :: eold,ff,gg,dgg
   double precision, save, allocatable, dimension (:) :: hr,hi
   double precision, save, allocatable, dimension (:,:) :: work,hvr,hvi
   complex (kind=K_C8), save, allocatable, dimension (:) :: vhv
   complex (kind=K_C8), save, allocatable, dimension (:,:) :: cv,hvc,cwork
!

contains
   subroutine cdavini(np,nvec)
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
      call subent('CDAVINI')
#endif
!
!     ----- mandatory arrays of eiprmn -----
!
      mvm = min((local_maxcor+1)*nvec,local_maxdav)
      allocate(cv(mvm,mvm),vhv((mvm*(mvm+1))/2),eold(mvm),ff(mvm),gg(mvm),  &
               dgg(mvm),work(np,mvm),  &
               stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'DAVINI error: allocation 1 failed with ios=',ios
         call pricod(ios)
         stop
      endif
      allocate(hvr(np,mvm),hvi(np,mvm),mseg(mvm),stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'DAVINI error: allocation 2 failed with ios=',ios
         call pricod(ios)
         stop
      endif
!
#ifdef _MYDEBUG_
      call subret('CDAVINI',1)
#endif
   end subroutine cdavini
   subroutine cdavend
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
      call subent('CDAVEND')
#endif
      deallocate(mseg,hvi,hvr,work,dgg,gg,ff,eold,vhv,cv)
#ifdef _MYDEBUG_
      call subret('CDAVEND',1)
#endif
   end subroutine cdavend

   subroutine eiprmn(eso,hd,vr,vi,n,mvec,ier,vecfil)
!
!*******************************************************************************
! ----- Davidson's diagonalisation for a complex hermitian matrix -----
!===============================================================================
!       VR (R8) = real part of the eigenvectors. (INPUT/OUTPUT)
!                 dimension vr(np,*)
!       VI (R8) = imaginary part of the eigenvectors. (INPUT/OUTPUT)
!                 dimension vi(np,*)
!                 In input, vr and pi contain the result of the guess
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
!       HD (R8) = array containing the diagonal of H. Minimal size N (INPUT)
!      ESO (R8) = working array of size at least equal to mv. (SCRATCH)
!
!     written by: J. P. Flament
!                 PhLAM
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      use jpsort
      use mphase
      use mercod
      use lapack95, only: hpev
!
!     ----- interfaces -----
!
      use chkprm_mod
      use fragh_mod
!     use savvso_mod
      use trf_mod
      use fiflsh_mod
      implicit none
!
!     ----- parameters -----
!
      double precision, parameter :: EPS=epsilon(ONE),TOCM=219474.0D+00
!
!     ----- arguments -----
!
      integer , intent(in) :: n,mvec
      integer , intent(out) :: ier
      integer , intent(in), optional :: vecfil
      double precision, intent(out), dimension (:) :: eso
      double precision, intent(in), dimension (:) :: hd
      double precision, intent(inout), dimension (:,:) :: vr,vi
!
!     ----- local variables -----
!
      character (len=1) :: reset
      integer :: i,ivec,ni,its,info,wkdim,nseg,n0,n1,n2,nfrag, &
                             ifrag,ios,ndiag,nprev,ncor,ier0,wk0,np,mv
      double precision :: dumi,dumr,denom,ffm,ggm,dggm
      logical (kind=K_L4) :: permut,per
!
!     ----- allocatable variables -----
!
      complex (kind=K_C8), allocatable, dimension (:,:) :: itemp,otemp
!
!     ----- formats -----
!
  9999 format(2x,a1,2x,'ITS=',i4,(t15,' IVEC=',i3,' E=',f12.8,'  DE=',d12.4,  &
              '  GG=',d12.4,'  DV=',d12.4))
  9998 format(//,5x,'*** EIPRMN WARNING: permutation(s) of ',  &
              'eigenvectors occurred during the iterations.',/,5x,  &
              '*** SOME EIGENVALUES MAY BE MISSING')
  9997 format(//,5x,'ERROR in ZHPEV called by EIPRMN')
  9996 format(/,10x,'LDZ=',i4,' but N=',i4)
  9995 format(/,5x,'EIPRMN ERROR: wkdim = 0. Increase parameter HBUF ',  &
              'in file spinor.f')
  9989 format(/,5x,i5,' elements of the triadonal matrix failed to converge')
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
      call subent('EIPRMN')
#endif
      if(n.eq.1) then
         eso(1) = hd(1)
         vr(1,1) = ONE
         vi(1,1) = ZERO
         ier = 111
#ifdef _MYDEBUG_
         call subret('EIPRMN',1)
#endif
         return
      endif
      np = size(vr,1)
      call cdavini(np,mvec)
      mv = mvm
      wk0 = 0
      reset = ' '
      ier0 = 1000
      ier = 0
      ivec = 0
      ni = mvec
      call dav_schmi(vr,vi,n,ivec,ni,stol)
      call davhv(hvr,hvi,vr,vi,n,mvec)
      do ivec=1,mvec
         call calvhv_symmetric(vhv,hvr,hvi,vr,vi,n,ivec)
      enddo
      where(abs(vhv(1:mvec*(mvec+1)/2)).lt.TENM10) vhv(1:mvec*(mvec+1)/2) = CZERO
      call hpev(vhv(1:mvec*(mvec+1)/2),eso,'U',cv,info)
      if(info.lt.0) then
         write(6,9997)
         write(6,9996) mv,mvec
         stop
      else if(info.gt.0) then
         write(6,9997)
         write(6,9989) info
         stop
      endif
      wkdim = 4*hbuf/mvec
!     if(wkdim.ne.wk0) then
!        wk0 = wkdim
         if(allocated(itemp)) deallocate(itemp,otemp)
         allocate(itemp(wkdim,mvec),otemp(wkdim,mvec),stat=ios)
         if(ios.ne.0) then
            write(6,"(5x,a)") 'EIPRMN error: allocation failed for itemp'
            call pricod(ios)
            stop
         endif
!     endif
      call trf(cv,vr,vi,n,mvec,itemp,otemp,wkdim)
      call trf(cv,hvr,hvi,n,mvec,itemp,otemp,wkdim)
!
!    ------------------------------
!    ----- loop over segments -----
!    ------------------------------
!
      call fragh(eso,mvec,mvm,mseg,nfrag)
      vhv(1:mvec*(mvec+1)/2) = CZERO
      do ivec=1,mvec
         vhv(ivec*(ivec+1)/2) = cmplx(eso(ivec),ZERO,kind=K_C8)
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
!     ----- vr(ni)/vi(ni) contain (H-e(ivec)*I)|ivec> -----
!
         ggm = ZERO
         dggm = ZERO
         ffm = ZERO
         nprev = ndiag
         ni = ndiag
         do ivec=n1,n2
            ni = ni+1
            gg(ivec) = ZERO
            dgg(ivec) = ZERO
            do i=1,n
               denom = ONE/(hd(i)-eso(ivec)+EPS)
               dumr = hvr(i,ivec)-eso(ivec)*vr(i,ivec)
               dumi = hvi(i,ivec)-eso(ivec)*vi(i,ivec)
               vr(i,ni) = dumr*denom
               vi(i,ni) = dumi*denom
               gg(ivec) = gg(ivec)+dumr**2+dumi**2
            enddo
            gg(ivec) = sqrt(gg(ivec))
            ggm = gg(ivec)
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
            call dav_schmi(vr,vi,n,nprev,ncor,stol)
            ndiag = nprev+ncor
            if(ncor.gt.0) then
               call davhv(hvr(1:n,nprev+1:ndiag),hvi(1:n,nprev+1:ndiag),  &
                           vr(1:n,nprev+1:ndiag),vi(1:n,nprev+1:ndiag),  &
                           n,ncor)
            endif
            ni = nprev
            do ivec=1,ncor
               ni = nprev+ivec
               call calvhv_symmetric(vhv,hvr,hvi,vr(1:n,1:ni),vi(1:n,1:ni),n,ni)
            enddo
            where(abs(vhv(1:ndiag*(ndiag+1)/2)).lt.TENM10)  &
                      vhv(1:ndiag*(ndiag+1)/2) = CZERO
            eold(1:nprev) = eso(1:nprev)
            call hpev(vhv(1:ndiag*(ndiag+1)/2),eso,'U',cv,info)
            if(info.lt.0) then
               write(6,9997)
               write(6,9996) mv,ndiag
               stop
            else if(info.gt.0) then
               write(6,9997)
               write(6,9989) info
               stop
            endif
            per = chkprm(cv,eso,ndiag,mvec)
            wkdim = 4*hbuf/ndiag
            if(wkdim.eq.0) then
               write(6,9995)
               stop
            endif
!           if(wkdim.ne.wk0) then
!              wk0 = wkdim
               if(allocated(itemp)) deallocate(otemp,itemp)
               allocate(itemp(wkdim,ndiag),otemp(wkdim,ndiag),stat=ios)
               if(ios.ne.0) then
                  write(6,"(5x,a)") 'EIPRMN error: allocation failed for itemp'
                  call pricod(ios)
                  stop
               endif
!           endif
            call trf(cv,vr,vi,n,ndiag,itemp,otemp,wkdim)
            call trf(cv,hvr,hvi,n,ndiag,itemp,otemp,wkdim)
            if(ndiag+nseg.gt.mvm) then
               reset = '*'
               ndiag = mvec
            endif
            nprev = ndiag
            vhv(1:ndiag*(ndiag+1)/2) = CZERO
            do i=1,ndiag
               vhv(i*(i+1)/2) = cmplx(eso(i),ZERO,kind=K_C8)
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
                  denom = ONE/(hd(i)-eso(ivec)+EPS)
                  dumr = hvr(i,ivec)-eso(ivec)*vr(i,ivec)
                  dumi = hvi(i,ivec)-eso(ivec)*vi(i,ivec)
                  vr(i,ni) = dumr*denom
                  vi(i,ni) = dumi*denom
                  gg(ivec) = dumr**2+dumi**2+gg(ivec)
                  dgg(ivec) = vr(i,ni)**2+vi(i,ni)**2+dgg(ivec)
               enddo
               gg(ivec) = sqrt(gg(ivec))
               ggm = max(ggm,gg(ivec))
               dgg(ivec) = sqrt(dgg(ivec))
               dggm = max(dggm,dgg(ivec))
               ff(ivec) = abs((eso(ivec)-eold(ivec))/eso(ivec))
               ffm = max(ffm,ff(ivec))
               if(ff(ivec).lt.ftol.and.gg(ivec).lt.gtol) ni = ni-1
            enddo
            ncor = ni-ndiag
            write(6,9999) reset,its,(ivec,eso(ivec),ff(ivec),gg(ivec),dgg(ivec),  &
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
               exit
            else
               ier = -1
            endif
!
!     ----- Save eigenvectors on vecfil file for subsquent restart -----
!
            if(present(vecfil).and.mod(its,5).eq.0) then
!              call savvso(vr,vi,eso,eold,n,mvec,vecfil)
               write(6,*) 'dans eiprmn pas encore de sauvegarde'
            endif
         enddo
!
!               --------------------------
!               ----- end iterations -----
!               --------------------------
!
         n0 = n0+nseg
         if(ifrag.lt.nfrag) then
            call davhv(hvr(1:n,1:mvec),hvi(1:n,1:mvec),  &
                        vr(1:n,1:mvec),vi(1:n,1:mvec),n,mvec)
            do ivec=1,mvec
               call calvhv_symmetric(vhv,hvr,hvi,vr,vi,n,ivec)
            enddo
            call hpev(vhv(1:mvec*(mvec+1)/2),eso,'U',cv,info)
            if(info.lt.0) then
               write(6,9997)
               write(6,9996) mv,mvec
               stop
            else if(info.gt.0) then
               write(6,9997)
               write(6,9989) info
               stop
            endif
            per = chkprm(cv,eso,mvec,mvec)
            wkdim = 4*hbuf/mvec
            if(wkdim.eq.0) then
               write(6,9995)
               stop
            endif
!           if(wkdim.ne.wk0) then
!              wk0 = wkdim
               if(allocated(itemp)) deallocate(otemp,itemp)
               allocate(itemp(wkdim,mvec),otemp(wkdim,mvec),stat=ios)
               if(ios.ne.0) then
                  write(6,"(5x,a)") 'EIPRMN error: allocation failed for itemp'
                  call pricod(ios)
                  stop
               endif
!           endif
            call trf(cv,vr,vi,n,mvec,itemp,otemp,wkdim)
            call trf(cv,hvr,hvi,n,mvec,itemp,otemp,wkdim)
            vhv(1:mvec*(mvec+1)/2) = CZERO
            do i=1,mvec
               vhv(i*(i+1)/2) = cmplx(eso(i),ZERO,kind=K_C8)
            enddo
         endif
      enddo
!          ----- end of loop over slices of energy -----
!
!     ----- sort final eigenvectors in increasing eigenvalue order -----
!
      allocate(indix(mvec))
      indix(1:mvec) = (/ (i,i=1,mvec) /)
      call gensort(mvec,eso,indix)
      hvr(1:n,1:mvec) = vr(1:n,1:mvec)
      hvi(1:n,1:mvec) = vi(1:n,1:mvec)
      vr(1:n,1:mvec) = hvr(1:n,indix(1:mvec))
      vi(1:n,1:mvec) = hvi(1:n,indix(1:mvec))
      ier = ier0
      call cphase(vr,vi,n,mvec)
      deallocate(indix)
      if(allocated(itemp)) deallocate(otemp,itemp)
      call cdavend
#ifdef _MYDEBUG_
      call subret('EIPRMN',2)
#endif
      return
   end subroutine eiprmn


   subroutine cdavini1(np,nvec)
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
      integer  :: ios,mvmp
!
!     ----- allocatable variables -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('DAVINI1')
#endif
!
!     ----- mandatory arrays of eiprmn or cdeigen -----
!
      mvm = min((local_maxcor+1)*nvec,local_maxdav)
      mvmp = mvm+iand(mvm,1)
      allocate(cv(mvmp,mvm),vhv((mvm*(mvm+1))/2),eold(mvm),ff(mvm),gg(mvm),  &
               dgg(mvm),cwork(np,mvm),  &
               stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'DAVINI1 error: allocation 1 failed with ios=',ios
         call pricod(ios)
         stop
      endif
      allocate(hvc(np,mvm),mseg(mvm),stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'DAVINI1 error: allocation 2 failed with ios=',ios
         call pricod(ios)
         stop
      endif
!
!     ----- arrays specific to the CI method (if needed) -----
!
#ifdef _MYDEBUG_
      call subret('DAVINI1',1)
#endif
   end subroutine cdavini1
   subroutine cdavend1
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
      call subent('DAVEND1')
#endif
      deallocate(mseg,hvc,cwork,dgg,gg,ff,eold,vhv,cv)
#ifdef _MYDEBUG_
      call subret('DAVEND1',1)
#endif
   end subroutine cdavend1

   subroutine cdeigen(eso,hd,vc,n,mvec,ier,vecfil)
!
!*******************************************************************************
! ----- Davidson's diagonalisation for a complex hermitian matrix -----
!===============================================================================
!       VC (C8) = complex eigenvectors. (INPUT/OUTPUT)
!                 dimension vc(np,*)
!                 In input, vr and pi contain the result of the guess
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
!       HD (R8) = array containing the diagonal of H. Minimal size N (INPUT)
!      ESO (R8) = working array of size at least equal to mv. (SCRATCH)
!
!     written by: J. P. Flament
!                 PhLAM
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
      use jpsort
      use mphase
!     use mlinalg
      use lapack95, only: hpev
      use blas95, only:gemm
!
!     ----- interfaces -----
!
!     use jpferr_mod
      use chkprm_mod
      use fragh_mod
!     use savvso_mod
      use fiflsh_mod
      implicit none
!
!     ----- parameters -----
!
      double precision, parameter :: EPS=epsilon(ONE),TOCM=219474.0D+00
!
!     ----- arguments -----
!
      integer , intent(in) :: n,mvec
      integer , intent(out) :: ier
      integer , intent(in), optional :: vecfil
      double precision, intent(out), dimension (:) :: eso
      double precision, intent(in), dimension (:) :: hd
      complex (kind=K_C8), intent(inout), dimension (:,:) :: vc
!
!     ----- local variables -----
!
      character (len=1) :: reset
      integer  :: i,ivec,ni,its,info,nseg,n0,n1,n2,nfrag,ifrag,  &
                             ndiag,nprev,ncor,ier0,wk0,np,mv,ios
      double precision :: denom,ffm,ggm,dggm,dumg,dumgg,xeso,extra,dter
      complex (kind=K_C8) :: dumc
      logical (kind=K_L4) :: permut,per
      double precision, dimension (3,4) :: amat
!
!     ----- allocatable variables -----
!
      double precision, allocatable, dimension (:,:) :: xdg,xe
      complex (kind=K_C8), allocatable, dimension (:,:) :: itemp,otemp
!
!     ----- formats -----
!
  9999 format(2x,a1,2x,'ITS=',i4,(t15,' IVEC=',i3,' E=',f12.8,'  DE=',d12.4,  &
              '  GG=',d12.4,'  DV=',d12.4))
  9998 format(//,5x,'*** CDEIGEN WARNING: permutation(s) of ',  &
              'eigenvectors occurred during the iterations.',/,5x,  &
              '*** SOME EIGENVALUES MAY BE MISSING')
  9997 format(//,5x,'ERROR in ZHPEV called by CDEIGEN')
  9996 format(/,10x,'LDZ=',i4,' but N=',i4)
  9989 format(/,5x,i5,' elements of the triadonal matrix failed to converge')
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
      call subent('CDEIGEN')
#endif
      if(n.eq.1) then
         eso(1) = hd(1)
         vc(1,1) = ONE
         ier = 111
#ifdef _MYDEBUG_
         call subret('CDEIGEN',1)
#endif
         return
      endif
      allocate(xe(3,mvec),xdg(3,mvec),stat=ios)
      if(ios.ne.0) then
         write(6,"(5x,a,i4)") 'TESGAS error: allocation 0 failed with ios=',ios
         call pricod(ios)
         call jpferr('Progam stop in -CDEIGEN-')
      endif
      np = size(vc,1)
      call cdavini1(np,mvec)
      mv = mvm
      wk0 = 0
      reset = ' '
      ier0 = 1000
      ier = 0
      ivec = 0
      ni = mvec
      call dav_schmi(vc,n,ivec,ni,stol)
      call davhv(hvc,vc,n,mvec)
      ITS = 0
      do ivec=1,mvec
         call calvhv_symmetric(vhv,hvc,vc,n,ivec)
      enddo
      where(abs(vhv(1:mvec*(mvec+1)/2)).lt.TENM10)  &
                vhv(1:mvec*(mvec+1)/2) = CZERO
      call hpev(vhv(1:(mvec*(mvec+1))/2),eso,'U',cv,info)
      if(info.lt.0) then
         write(6,9997)
         write(6,9996) mv,mvec
         stop
      else if(info.gt.0) then
         write(6,9997)
         write(6,9989) info
         stop
      endif   
      call ctrf(cv,vc,n,mvec,cwork)
      call ctrf(cv,hvc,n,mvec,cwork)
!
!    ------------------------------
!    ----- loop over segments -----
!    ------------------------------
!
      call fragh(eso,mvec,mvm,mseg,nfrag)
      vhv(1:mvec*(mvec+1)/2) = CZERO
      do ivec=1,mvec
         vhv(ivec*(ivec+1)/2) = cmplx(eso(ivec),ZERO,kind=K_C8)
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
!     ----- vc(ni) contains (H-e(ivec)*I)|ivec> -----
!
         ggm = ZERO
         dggm = ZERO
         ffm = ZERO
         nprev = ndiag
         ni = ndiag
         do ivec=n1,n2
            ni = ni+1
            gg(ivec) = ZERO
            dgg(ivec) = ZERO
            do i=1,n
               denom = ONE/(hd(i)-eso(ivec)+EPS)
               dumc = hvc(i,ivec)-eso(ivec)*vc(i,ivec)
               vc(i,ni) = dumc*denom
               gg(ivec) = dumc*conjg(dumc)+gg(ivec)
            enddo
            gg(ivec) = sqrt(gg(ivec))
            ggm = gg(ivec)
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
            call dav_schmi(vc,n,nprev,ncor,stol)
            ndiag = nprev+ncor
            if(ncor.gt.0) call davhv(hvc(:,nprev+1:ndiag),  &
                                      vc(:,nprev+1:ndiag),n,ncor)
            ni = nprev
            do ivec=1,ncor
               ni = nprev+ivec
               call calvhv_symmetric(vhv,hvc,vc(:,1:ni),n,ni)
            enddo
            eold(1:nprev) = eso(1:nprev)
            call hpev(vhv(1:(ndiag*(ndiag+1))/2),eso,'U',cv,info)
            if(info.lt.0) then
               write(6,9997)
               write(6,9996) mv,ndiag
               stop
            else if(info.gt.0) then
               write(6,9997)
               write(6,9989) info
               stop
            endif
            per = chkprm(cv,eso,ndiag,mvec)
            call ctrf(cv,vc,n,ndiag,cwork)
            call ctrf(cv,hvc,n,ndiag,cwork)
            if(ndiag+nseg.gt.mvm) then
               reset = '*'
               ndiag = mvec
            endif
            nprev = ndiag
            vhv(1:ndiag*(ndiag+1)/2) = CZERO
            do i=1,ndiag
               vhv(i*(i+1)/2) = cmplx(eso(i),ZERO,kind=K_C8)
            enddo
!
!     ----- compute gradient and new direction -----
!
            ggm = ZERO
            dggm = ZERO
            ffm = ZERO
            ni = ndiag
            do ivec=n1,n2
               extra = eso(ivec)
               if(its.eq.1) then
                  xe(1,ivec) = eso(ivec)
                  xdg(1,ivec) = gg(ivec)
               else if(its.eq.2) then
                  xe(2,ivec) = eso(ivec)
                  xdg(2,ivec) = gg(ivec)
               else if(its.eq.3) then
                  xe(3,ivec) = eso(ivec)
                  xdg(3,ivec) = gg(ivec)
               else
                  xe(1,ivec) = xe(2,ivec)
                  xdg(1,ivec) = xdg(2,ivec)
                  xe(2,ivec) = xe(3,ivec)
                  xdg(2,ivec) = xdg(3,ivec)
                  xe(3,ivec) = eso(ivec)
                  xdg(3,ivec) = gg(ivec)
               endif
               if(its.ge.3) then
                  do i=1,3
                     xeso = xdg(i,ivec)
                     amat(i,1) = xeso**2
                     amat(i,2) = xeso
                     amat(i,3) = ONE
                     amat(i,4) = xe(i,ivec)
                  enddo
                  call resyeq(amat,3,dter)
                  extra = amat(3,4)
               endif
               dumg = ZERO
               dumgg = ZERO
               ni = ni+1
               do i=1,n
                  denom = ONE/(hd(i)-extra+EPS)
                  dumc = hvc(i,ivec)-eso(ivec)*vc(i,ivec)
                  vc(i,ni) = dumc*denom
                  dumg = dumc*conjg(dumc)+dumg
                  dumgg = vc(i,ni)*conjg(vc(i,ni))+dumgg
               enddo
               gg(ivec) = sqrt(dumg)
               ggm = max(ggm,gg(ivec))
               dgg(ivec) = sqrt(dumgg)
               dggm = max(dggm,dgg(ivec))
               ff(ivec) = abs((eso(ivec)-eold(ivec))/eso(ivec))
               ffm = max(ffm,ff(ivec))
               if(ff(ivec).lt.ftol.and.gg(ivec).lt.gtol) ni = ni-1
            enddo
            ncor = ni-ndiag
            write(6,9999) reset,its,(ivec,eso(ivec),ff(ivec),gg(ivec),dgg(ivec),  &
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
               goto 200
!              exit
            else
               ier = -1
            endif
!
!     ----- Save eigenvectors on vecfil file for subsquent restart -----
!
            if(present(vecfil).and.mod(its,5).eq.0) then
!              call savvso(vr,vi,eso,eold,n,mvec,vecfil)
               write(6,*) 'dans CDEIGEN pas encore de sauvegarde'
            endif
         enddo
!
!               --------------------------
!               ----- end iterations -----
!               --------------------------
!
     200 continue
         n0 = n0+nseg
         if(ifrag.lt.nfrag) then
            call davhv(hvc(:,1:mvec),vc(:,1:mvec),n,mvec)
            do ivec=1,mvec
               call calvhv_symmetric(vhv,hvc,vc,n,ivec)
            enddo
            call hpev(vhv(1:(mvec*(mvec+1))/2),eso,'U',cv,info)
            if(info.lt.0) then
               write(6,9997)
               write(6,9996) mv,mvec
               stop
            else if(info.gt.0) then
               write(6,9997)
               write(6,9989) info
               stop
            endif
            per = chkprm(cv,eso,mvec,mvec)
            call ctrf(cv,vc,n,mvec,cwork)
            call ctrf(cv,hvc,n,mvec,cwork)
            vhv(1:mvec*(mvec+1)/2) = CZERO
            do i=1,mvec
               vhv(i*(i+1)/2) = cmplx(eso(i),ZERO,kind=K_C8)
            enddo
         endif
      enddo
!          ----- end of loop over slices of energy -----
!
!     ----- sort final eigenvectors in increasing eigenvalue order -----
!
      allocate(indix(mvec))
      indix(1:mvec) = (/ (i,i=1,mvec) /)
      call gensort(mvec,eso,indix)
      hvc(:,1:mvec) = vc(:,1:mvec)
      vc(:,1:mvec) = hvc(:,indix(1:mvec))
      ier = ier0
      call cphase(vc,n,mvec)
      deallocate(indix)
      if(allocated(itemp)) deallocate(otemp,itemp)
      deallocate(xe,xdg)
      call cdavend1
#ifdef _MYDEBUG_
      call subret('CDEIGEN',2)
#endif
      return
   end subroutine cdeigen
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
      use blas95, only: gemm
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer , intent(in) :: ndet,nvec
      complex (kind=K_C8), intent(in), dimension (:,:) :: c
      complex (kind=K_C8), intent(inout), dimension (:,:) :: v,wc
!
!     ----- local variables -----
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
      call subent('CTRF')
#endif
      call gemm(v(1:ndet,1:nvec),c(1:nvec,1:nvec),wc(1:ndet,1:nvec),  &
                   'N','N',CONE,CZERO)
      v(1:ndet,1:nvec) = wc(1:ndet,1:nvec)
#ifdef _MYDEBUG_
      call subret('CTRF',1)
#endif
      return
   end subroutine ctrf
end module cdavidson_lille
