subroutine eiprmn(vr,vi,np,n,mvec,eold,ff,gg,dgg,ier,hbuf,hd,eso,mv,delta,d,  &
                   mdress,heff,sfile,elfile,vecfil)
!
!***********************************************************************
! ----- Davidson's diagonalisation for a complex hermitian matrix -----
!=======================================================================
!       PR (R8) = real part of the eigenvectors. (INPUT/OUTPUT)
!                 dimension vr(np,*)
!       PI (R8) = imaginary part of the eigenvectors. (INPUT/OUTPUT)
!                 dimension vi(np,*)
!                 In input, vr and pi contain the result of the guess
!       NP (I4) = physical dimension of vr and pi (INPUT). NP >=N
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
!  HVR,HVI (R8) = working array. Size =  (np,*). (SCRATCH)
!      VHV (C16)= working array. Size at least equal to mv*(mv+1)/2
!                  (SCRATCH)
!       CV (C16)= array cv(mv,mv) (SCRATCH)
!  HMR,HMI (R8) = buffer to read the Hamitonian matrix on disk, of size
!                 HBUF (SCRATCH)
!     ICOL (I8) = buffer of size HBUF to read the indices of the non zero
!                 elements of the Hamitonian matrix saved on disk (SCRATCH)
!     HBUF (I4) = dimension of the buffers (INPUT)
!       HD (R8) = array containing the diagonal of H. Minimal size N (INPUT)
!      ESO (R8) = working array of size at least equal to mv. (SCRATCH)
!       MV (I4) = dimension of CV, maximum size of the reduced space (INPUT)
!    DELTA (R8) = array of size mvec, which contains the energetic shift
!                 used for the contruction of the Effective Hamiltonian
!        D (R8) = array containing the eigenvectors of the diagonalisation
!                 of H, without SO and shift (INPUT)
!   MDRESS (I4) = number of vectors to be dressed for the effective
!                 Hamiltonian
!     HEFF (L4) = flag for the calculation of the effective Hamiltonian
!     W,WR (R8) = working arrays of size at least equal to 6*Mv (SCRATCH).
!    SFILE (I4) = unit number of the SO part of H (INPUT)
!   ELFILE (I4) = unit number of the non relativistic part of H (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use kindvar
   use mconst
   use mdavi
   use matpri
   use jpsort
   use mgenrl
   use mdebug
   use mphase
   use mercod
   use mmulsy
   use f95_precision, only: wp => dp
!  use lapack95, only: hpev
   use lapack95, only: hpevd
!
!     ----- interfaces -----
!
   use cschmi_mod
   use calchv_mod
   use calvhv_mod
   use chkprm_mod
   use savvso_mod
   use trf_mod
   use fragh_mod
   use fiflsh_mod
   implicit none
!
!     ----- parameters -----
!
!  double precision, parameter :: EPS=1.0D-20,TOCM=219474.0D+00
   double precision, parameter :: EPS=epsilon(1.0D0),TOCM=219474.0D+00
!
!     ----- arguments -----
!
   integer, intent(in) :: np,n,mvec,hbuf,mv,mdress,sfile,elfile,vecfil
   integer, intent(out) :: ier
   logical, intent(in) :: heff
   double precision, intent(out), dimension (:) :: eold,ff,gg,dgg,eso
   double precision, intent(in), dimension (:) :: hd,delta
   double precision, intent(in), dimension (:,:) :: d
   double precision, intent(inout), dimension (:,:) :: vr,vi
!
!     ----- local variables -----
!
   character (len=1) :: reset
   integer :: i,ivec,ni,its,info,wkdim,nseg,n0,n1,n2,nfrag,mvm,ifrag,ios,  &
              ndiag,nprev,ncor,ier0,wk0
   double precision :: dumi,dumr,denom,ffm,ggm,dggm
!  COMPLEX*16 :: CVAL,DVAL
   logical :: permut,per
!
!     ----- allocatable variables -----
!
   integer, allocatable, dimension (:) :: icol,mseg,indix,elstar,sostar
   double precision, allocatable, dimension (:) :: hmr,hmi
   double precision, allocatable, dimension (:,:) :: hvr,hvi,wdres
   complex (kind=K_C8), allocatable, dimension (:) :: vhv
   complex (kind=K_C8), allocatable, dimension (:,:) :: itemp,otemp,cv
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
  9994 format(/,5x,'EIPRMN error: allocation failed for HMR HMI ICOL')
  9993 format(/,5x,'EIPRMN error: allocation failed for VHV and CV')
! 9992 format(/,5x,'EIPRMN error: allocation failed for W and WR')
  9991 format(/,5x,'Spin-orbit matrix in the LS-basis (cm-1)',/,  &
                5x,'----------------------------------------')
  9990 format(//)
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
   mvm = mv
   allocate(hmr(hbuf),hmi(hbuf),icol(hbuf),elstar(n+1),sostar(n+1), &
            stat=ios)
   if(ios.ne.0) then
      write(6,9994)
      call pricod(ios)
      stop
   endif
   allocate(cv(mvm,mvm),vhv((mvm*(mvm+1))/2),stat=ios)
   if(ios.ne.0) then
      write(6,9993)
      call pricod(ios)
      stop
   endif
   allocate(hvr(np,mvm),hvi(np,mvm),wdres(hbuf,mvec),mseg(mvm),stat=ios)
   if(ios.ne.0) then
      write(6,9994)
      call pricod(ios)
      stop
   endif
!
   wk0 = 0
   rewind 93
   read(93) elstar,sostar
   reset = ' '
   ier0 = 1000
   ier = 0
   ivec = 0
   ni = mvec
   call cschmi(vr,vi,n,ivec,ni,stol)
   call calchv(hvr,hvi,vr,vi,n,mvec,hmr,hmi,icol,hbuf,delta,d,mdress,  &
               heff,wdres,sfile,elfile,elstar,sostar)
!  do ivec=1,mvec
!     call calvhv(vhv,hvr,hvi,vr,vi,n,ivec)
!  enddo
   call mmsy(vr(1:n,1:mvec),vi(1:n,1:mvec),hvr(1:n,1:mvec),hvi(1:n,1:mvec),  &
             vhv,n,mvec)
   where(abs(vhv(1:mvec*(mvec+1)/2)).lt.TENM10) vhv(1:mvec*(mvec+1)/2) = CZERO
   if(index(check,'SO/ST').ne.0) then
      write(6,9991)
      vhv(1:mvec*(mvec+1)/2) = vhv(1:mvec*(mvec+1)/2)*TOCM
      call carmat(vhv,mvec)
      vhv(1:mvec*(mvec+1)/2) = vhv(1:mvec*(mvec+1)/2)/TOCM
      write(6,9990)
   endif
!  call hpev(vhv(1:mvec*(mvec+1)/2),eso,'U',cv,info)
   call hpevd(vhv(1:mvec*(mvec+1)/2),eso,'U',cv,info)
   if(info.ne.0) then
      write(6,9997)
      write(6,9996) mv,mvec
      stop
   endif
   wkdim = 4*hbuf/mvec
   if(allocated(itemp)) deallocate(itemp,otemp)
   allocate(itemp(wkdim,mvec),otemp(wkdim,mvec),stat=ios)
   if(ios.ne.0) then
      write(6,"(5x,a)") 'EIPRMN error: allocation failed for itemp'
      call pricod(ios)
      stop
   endif
   call trf(cv,vr,vi,n,mvec,itemp,otemp,wkdim)
   call trf(cv,hvr,hvi,n,mvec,itemp,otemp,wkdim)
!
!    ------------------------------
!    ----- loop over segments -----
!    ------------------------------
!
   call fragh(eso,mvec,mvm,mseg,nfrag)
   vhv(1:mvec*(mvec+1)/2) = dcmplx(ZERO,ZERO)
   do ivec=1,mvec
      vhv(ivec*(ivec+1)/2) = dcmplx(eso(ivec),ZERO)
   enddo
   n0 = 0
   do ifrag=1,nfrag
      nseg = mseg(ifrag)
!
!    ----- n1, n2: limits of the current fragment ----
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
      nprev = ndiag
      ni = ndiag
      do ivec=n1,n2
         ni = ni+1
         ggm = ZERO
         do i=1,n
            denom = ONE/(hd(i)-eso(ivec)+EPS)
            dumr = hvr(i,ivec)-eso(ivec)*vr(i,ivec)
            dumi = hvi(i,ivec)-eso(ivec)*vi(i,ivec)
            vr(i,ni) = dumr*denom
            vi(i,ni) = dumi*denom
            ggm = ggm+dumr**2+dumi**2
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
         call cschmi(vr,vi,n,nprev,ncor,stol)
         ndiag = nprev+ncor
         if(ncor.gt.0) then
            call calchv(hvr(1:np,nprev+1:ndiag),hvi(1:np,nprev+1:ndiag),  &
                        vr(1:np,nprev+1:ndiag),vi(1:np,nprev+1:ndiag),n,  &
                        ncor,hmr,hmi,icol,hbuf,delta,d,mdress,heff,wdres,  &
                        sfile,elfile,elstar,sostar)
         endif
!        ni = nprev
!        do ivec=1,ncor
!           ni = nprev+ivec
!           call calvhv(vhv,hvr,hvi,vr(1:np,1:ni),vi(1:np,1:ni),n,ni)
!        enddo
         call mmsy(vr(1:n,1:ndiag),vi(1:n,1:ndiag),hvr(1:n,1:ndiag), &
                   hvi(1:n,1:ndiag),vhv,n,ndiag)
         do ivec=1,nprev
            eold(ivec) = eso(ivec)
         enddo
!        call hpev(vhv(1:ndiag*(ndiag+1)/2),eso,'U',cv,info)
         call hpevd(vhv(1:ndiag*(ndiag+1)/2),eso,'U',cv,info)
         if(info.ne.0) then
            write(6,9997)
            write(6,9996) mv,ndiag
            stop
         endif
         per = chkprm(cv,eso,ndiag,mvec)
         wkdim = 4*hbuf/ndiag
         if(wkdim.eq.0) then
            write(6,9995)
            stop
         endif
         if(allocated(itemp)) deallocate(itemp,otemp)
         allocate(itemp(wkdim,ndiag),otemp(wkdim,ndiag),stat=ios)
         if(ios.ne.0) then
            write(6,"(5x,a)") 'EIPRMN error: allocation failed for itemp'
            call pricod(ios)
            stop
         endif
         call trf(cv,vr,vi,n,ndiag,itemp,otemp,wkdim)
         call trf(cv,hvr,hvi,n,ndiag,itemp,otemp,wkdim)
         if(ndiag+nseg.gt.mvm) then
            reset = '*'
            ndiag = mvec
         endif
         nprev = ndiag
         vhv(1:ndiag*(ndiag+1)/2) = dcmplx(ZERO,ZERO)
         do i=1,ndiag
            vhv(i*(i+1)/2) = dcmplx(eso(i),ZERO)
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
            goto 200
         else
            ier = -1
         endif
         if(mod(its,5).eq.0) then
            call savvso(vr,vi,eso,eold,n,mvec,vecfil)
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
         call calchv(hvr,hvi,vr,vi,n,mvec,hmr,hmi,icol,hbuf,delta,d,  &
                     mdress,heff,wdres,sfile,elfile,elstar,sostar)
!        do ivec=1,mvec
!           call calvhv(vhv,hvr,hvi,vr,vi,n,ivec)
!        enddo
         call mmsy(vr(1:n,1:mvec),vi(1:n,1:mvec),hvr(1:n,1:mvec),  &
                   hvi(1:n,1:mvec),vhv,n,mvec)
!        call hpev(vhv(1:mvec*(mvec+1)/2),eso,'U',cv,info)
         call hpevd(vhv(1:mvec*(mvec+1)/2),eso,'U',cv,info)
         if(info.ne.0) then
            write(6,9997)
            write(6,9996) mv,mvec
            stop
         endif
         per = chkprm(cv,eso,mvec,mvec)
         wkdim = 4*hbuf/mvec
         if(wkdim.eq.0) then
            write(6,9995)
            stop
         endif
         if(allocated(itemp)) deallocate(itemp,otemp)
         allocate(itemp(wkdim,mvec),otemp(wkdim,mvec),stat=ios)
         if(ios.ne.0) then
            write(6,"(5x,a)") 'EIPRMN error: allocation failed for itemp'
            call pricod(ios)
            stop
         endif
         call trf(cv,vr,vi,n,mvec,itemp,otemp,wkdim)
         call trf(cv,hvr,hvi,n,mvec,itemp,otemp,wkdim)
         vhv(1:mvec*(mvec+1)/2) = dcmplx(ZERO,ZERO)
         do i=1,mvec
            vhv(i*(i+1)/2) = dcmplx(eso(i),ZERO)
         enddo
      endif
   enddo
!
!     ----- sort final eigenvectors in increasing eigenvalue order -----
!
   deallocate(hmr,hmi,icol,wdres,cv,vhv)
   allocate(indix(mvec))
   indix(1:mvec) = (/ (i,i=1,mvec) /)
   call gensort(mvec,eso,indix)
   hvr(1:n,1:mvec) = vr(1:n,1:mvec)
   hvi(1:n,1:mvec) = vi(1:n,1:mvec)
   vr(1:n,1:mvec) = hvr(1:n,indix(1:mvec))
   vi(1:n,1:mvec) = hvi(1:n,indix(1:mvec))
   call cphase(vr,vi,n,mvec)
   deallocate(hvr,hvi,indix)
   if(allocated(itemp)) deallocate(itemp,otemp)
#ifdef _MYDEBUG_
   call subret('EIPRMN',2)
#endif
   return
end subroutine eiprmn
