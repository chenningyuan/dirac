      subroutine bummer( text, ierr, errtyp )
c
c  process a program error.
c  input:
c  text  = character string to be printed.
c  ierr  = internal program error to be printed.
c  errtyp = 0 for warning.  traceback may be generated. execution
c             continues.
c         = 1 for nonfatal error.  traceback may be generated.
c             execution is stopped. jcl condition code is set to allow
c             subsequent program steps to continue if possible.
c         = 2 for fatal error.  traceback may be generated.
c             execution is stopped. jcl condition code is set to abort
c             subsequent program steps if possible.
c
c  entry ibummr must be called prior to bummer() to set the output
c  unit and to perform any additional initialization.
c
c  version log:
c  24-apr-92 %val() blocks added for ibm rs6000. -rls
c  11-sep-91 parerr() calls added. -rjh
c  13-mar-91 posix code added. -rls
c  01-mar-89 write(stderr,*) and call exit() for unix machines (rls).
c  05-jul-88 unicos version (rls).
c  03-nov-87 ibm version (dcc).
c  10-sep-87 written by ron shepard.
c
      implicit integer(a-z)
      character*(*) text
      integer iunit, errtyp, ierr
      integer f77err
      integer nlist,        stderr
      save    nlist,        stderr
      data    nlist / 6 /,  stderr / 0 /
 
      if ( errtyp .eq. 0 ) then
c        # print a warning message and continue execution.
         write(nlist,6010) 'bummer (warning):', text, ierr
         write(stderr,6010) 'bummer (warning):', text, ierr
         return
      elseif ( errtyp .eq. 1 ) then
c        # print a warning message, stop execution.
         write(nlist,6010) 'bummer (nonfatal):', text, ierr
         write(stderr,6010) 'bummer (nonfatal):', text, ierr
         call exit( 0 )
      elseif ( errtyp .eq. 2 ) then
c        # print an error message, stop execution, and abort job
c        # sequence.
         write(nlist,6010) 'bummer (fatal):', text, ierr
         write(stderr,6010) 'bummer (fatal):', text, ierr
         call exit( 1 )
      else
c        # unknown error level.  treat as a fatal error.
         write(nlist,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
         write(stderr,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
         call exit( 1 )
      endif
c
c     # this statement is not executed, it is included
c     # just to avoid compiler warnings. -rls
      stop 'bummer error'
c
c     # initialization...
c
      entry ibummr( iunit )
c
c     # save the listing unit for use later.
c
      nlist = iunit
c
      return
6010  format(1x,a,a,i10)
6020  format(1x,a,i10,a,i10)
      end

*deck sifcfg
      subroutine sifcfg(
     & nbft,  ifmt,    lrec,    nmax,    ierr )
c
c  return a set of consistent configuration parameters for a standard
c  integral file structure.  ( for 1-e integrals )
c
c  input:
c  nbft = total number of basis functions.
c
c  notes: (1) ifmt, lrec, and nmax may eventually have meaning on input
c         in future versions of this routine.  for now the input values
c         are ignored.
c         (2) for extensibility, the input variables may be passed
c         into this routine as array elements, instead of scalars, in
c         future versions of this routine.
c
c  output:
c  ifmt = ifmt parameter.
c  lrec = actual buffer length to be written.
c  nmax = maximum number of values in each record.
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 for itype error.
c       = -2 for nbft error.
c       = -3 for lrec error.
c
      implicit none
      integer nbft, ifmt, lrec, nmax, ierr
      integer nbig, chunk, nword, space
c
c     # nbftmx = largest nbft consistent with label packing methods.
c     # lrecmn = minimum allowed record length.
c     # lrecmx = maximum allowed record length.  this should be
c     #          consistent with the record "dword" encoding.
c     # lrcinc = record length block size increments.
c     # ndeflt = default number of block increments. such that
c     #          the default record length = (lrecmn + ndeflt * lrcinc)
c  
      integer    nbftmx
      parameter( nbftmx=2**16-1 )
      integer    lrecmn,     lrecmx,         lrcinc,      ndeflt
      parameter( lrecmn=512, lrecmx=2**16,   lrcinc=512,  ndeflt=3 )
c
c     # check the input parameters...
c
      ierr = 0
c
      if ( nbft .le. 0    .or.    nbft .gt. nbftmx ) then
         ierr = -2
         return
      endif
c
c     # adjust lrec for efficient I/O.
c     # machine-dependent quirks, magic record lengths, etc. should
c     # be localized here and in the above parameter definitions..
c
      lrec = lrecmn + (ndeflt * lrcinc)

      lrec = min( lrecmx, lrec )
      if ( lrec .lt. lrecmn ) then
         ierr = -3
         return
      endif
c
c     # for this record length, compute the number of values.
c
c     # 1-e records: dword // values // fcore // packed_labels // ibitv
c     # 2-e records: dword // values //          packed_labels // ibitv
c
      if ( nbft .lt. 2**8 ) then
c        # use standard 8-bit orbital label packing.
         ifmt=0
      else
c        # use standard 16-bit orbital label packing.
         ifmt=1
      endif
c
c     # compute a reasonable upper bound to nmax.
c     # chunk is used to constrain nmax so that ulab8() or ulab16() do
c     # not overwrite the labels(*) array.
c     #
c     # another choice would be to require the programmer to account
c     # for these overruns when allocating labels(*).  this would
c     # result in more efficient buffer(*) use, but the chunk choice
c     # is simpler.
c
c     # note that the programmer must explicitly account for ibitv(*)
c     # overruns; constraining chunk to 64 is too wasteful. -rls
c
      if ( ifmt .eq. 0 ) then
c        # n + n/4 <= lrec
         nbig  = ( 4 * lrec ) / 5
         chunk = 4
      elseif ( ifmt .eq. 1 ) then
c        # n + n/2 <= lrec
         nbig  = ( 2 * lrec ) / 3
         chunk = 2
      endif
c
c     # round up to a higher multiple of chunk.
      nbig = ( ( nbig + 2*chunk ) / chunk ) * chunk
c
c     # account for dword and fcore.
      nword = 2
c
c     # loop backwards from the upper bound until a valid nmax is found.
      do 10 nmax = nbig, 0, -chunk
         space = nword + nmax + ( ( nmax + chunk - 1) / chunk ) 
         if ( space .le. lrec ) goto 11
 10   continue
 11   continue
      if ( nmax .le. 0 ) then
         ierr = -4
         return
      endif
c
c     # ifmt, lrec, and nmax are all ok.
c
      return
      end

*deck sifpre
      subroutine sifpre( nlist, nenrgy, energy, ietype )
c
c  print out the energy(*) array with character labels.
c
c  24-apr-92 cnvginf values added. -rls
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
      integer  nlist,  nenrgy, ietype(*)
      real*8   energy(*)
      integer  i,      itypea, itypeb
      character*8      chrtyp
c
      do 10 i = 1, nenrgy
         itypea = ietype(i) / 1024
         itypeb = mod( ietype(i), 1024 )
         call siftyp( itypea, itypeb, chrtyp )
c
         if ( ietype(i) .ge. 0 ) then
c           # fcore from a 1-e hamiltonian array.
            write(nlist,6010) i, energy(i), ietype(i), 'fcore', chrtyp
         elseif ( itypea .eq. 0 ) then
c           # core energy value.
            write(nlist,6010) i, energy(i), ietype(i), 'core', chrtyp
         elseif ( itypea .eq. -1 ) then
c           # total energy value.
            write(nlist,6010) i, energy(i), ietype(i), 'total', chrtyp
         elseif ( itypea .eq. -2 ) then
c           # energy or wave function convergence value.
            write(nlist,6010) i, energy(i), ietype(i), 'cnvginf', chrtyp
         else
c           # undefined energy type.
            write(nlist,6010) i, energy(i), ietype(i), 'unknown', chrtyp
         endif
 10   continue
      return
c
 6010 format(' energy(', i2, ')=', 1pe20.12, ', ietype=', i5,
     & ', ', a7, ' energy of type: ', a )
      end

 

*deck sifsce
      function sifsce( nenrgy, energy, ietype )
c
c  select and sum the core energies in the energy(*) array.
c
c  usage: real*8   sifsce
c         external sifsce
c         total_core = sifsce(...)
c
c  input: nenrgy = number of energy(*) values.
c         ietype(1:nenrgy) = energy types.
c         energy(1:nenrgy) = energy array.
c
c  output: sifsce = total core energy such that
c                       total_potential = total_electronic + sifsce(...)
c                   is the total clamped-nucleus, born-oppenheimer
c                   potential.
c
      implicit none
      integer  nenrgy, ietype(nenrgy)
      real*8   sifsce, energy(nenrgy)
      integer  i,      itypea
      real*8   ecore
c
      ecore = (0)
      do i = 1, nenrgy
c        # 0 <= ietype elements are 1-e hamiltonian frozen core terms.
c        # ietype < 0  elements with itypea=0 are other core terms to
c                      be added.  e.g. ietype=-1 = nuclear repulsion.
         itypea = ietype(i) / 1024
         if ( (ietype(i) .ge. 0 ) .or.
     &    ( (ietype(i) .lt. 0) .and. (itypea .eq. 0) ) ) then
            ecore = ecore + energy(i)
         endif
      enddo   
      sifsce = ecore
      return
      end


*deck siftyp
      subroutine siftyp( itypea, itypeb, chrtyp )
c
c  return a character description of the integral type.
c
c  input:
c  itypea, itypeb = generic and specific integral or energy(*) types.
c
c  output:
c  chrtyp = character description. (this should be at least character*8
c           in the calling program.)
c
c  08-oct-90 total energy and convergence types added. -rls
c  01-jun-90 table version. -rls
c  11-aug-89 written by ron shepard.
c
      implicit none
c
c     # dummy...
      integer          itypea,  itypeb
      character*(*)    chrtyp
c
c     # local...
      integer  i
c
      integer  typ1e,  st1e,   end1e
      integer  typ2e,  st2e,   end2e
      integer  typc,   stc,    endc
      integer  typte,  stte,   endte
      integer  typcv,  stcv,   endcv
c
c     # typ1e = number of defined 1-e array types.
c     # typ2e = number of defined 2-e array types.
c     # typc  = number of core energy types.
c     # typte = number of total energy types.
c     # typcv = number of convergence types.
c
      parameter( typ1e = 16, st1e = 1,        end1e = typ1e        )
      parameter( typ2e =  2, st2e = end1e +1, end2e = end1e +typ2e )
      parameter( typc  =  1, stc  = end2e +1, endc  = end2e +typc  )
      parameter( typte = 14, stte = endc  +1, endte = endc  +typte )
      parameter( typcv = 10, stcv = endte +1, endcv = endte +typcv )
c
      integer    ntype
      parameter( ntype = endcv )
c
      integer ltypea(ntype), ltypeb(ntype)
      character*8 lctype(ntype)
c
c     # these array types can be in any order, so new ones can be added
c     # to the end or inserted into the middle as appropriate. -rls
c
c     # warning: case dependent code:  do not change the case of the
c     # following character strings.
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st1e, end1e ) /
     & 0,    0, 'S1(*)',     0,    1, 'T1(*)',
     & 0,    2, 'V1(*)',     0,    3, 'Veff(*)',
     & 0,    4, 'VFC(*)',    0,    5, 'Vref(*)',
     & 0,    6, 'H1(*)',     0,    7, 'D1(*)',
     & 0,    8, 'F(*)',      0,    9, 'Q(*)',
     & 1,    0, 'X(*)',      1,    1, 'Y(*)',
     & 1,    2, 'Z(*)',
     & 2,    0, 'Im(SO:x)',  2,    1, 'Im(SO:y)',
     & 2,    2, 'Im(SO:z)'  /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st2e, end2e ) /
     & 3,    0, '1/r12',     3,    1, 'd2(*)'     /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stc, endc ) /
     & 0,   -1, 'Nuc.Rep.' /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stte, endte ) /
     & -1,   0, 'SCF',      -1,   -1, 'MCSCF',
     & -1,  -2, 'MRSDCI',   -1,   -3, 'CPF',
     & -1,  -4, 'ACPF',     -1,   -5, 'LCC-SD',
     & -1,  -6, 'MRPT',     -1,   -7, 'Bk',
     & -1,  -8, 'DV1',      -1,   -9, 'DV2',
     & -1, -10, 'EPOPLE',   -1,  -11, 'S.O. CI',
     & -1, -12, 'SR-SDCI',  -1,  -13, 'UCEPA'  /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stcv, endcv ) /
     & -2,   0, 'SCF-D.E.', -2,   -1, 'SCF-D.D1',
     & -2,  -2, 'MC-D.E.',  -2,   -3, 'MC-Wnorm',
     & -2,  -4, 'MC-Knorm', -2,   -5, 'MC-ApxDE',
     & -2,  -6, 'Bk-Resid', -2,   -7, 'CI-Resid',
     & -2,  -8, 'CI-D.E.',  -2,   -9, 'CI-ApxDE' /
c
      do 10 i = 1, ntype
         if ( itypea .eq. ltypea(i) ) then
            if ( itypeb .eq. ltypeb(i) ) then
               chrtyp = lctype(i)
               return
            endif
         endif
 10   continue
c     # loop exit means unrecognized type.
      chrtyp = 'Unknown'
c
      return
      end



      subroutine setiarray(
     &   itypea, itypax, itypbx, btypmx, btypes, iarray)
!
!     # set iarray to point into array(*,iarray).
!
      integer itypea, itypax, itypbx, btypmx
      integer btypes(0:btypmx), iarray

         if ( itypea .eq. -1 ) then
!           # single array mode. accumulate into array(*,1) always.
            iarray = 1
         elseif ( itypax.eq.0 ) then
!           # correct generic array type for s, t, v, c.
!           # mask btypmx with btypes(*) to see if this array should
!           # be accumulated somewhere.
            if ( itypbx .le. btypmx) then
               iarray = btypes(itypbx)
            else
               iarray = -1
            endif
         elseif ( itypax.eq.2 ) then
!           # correct generic array type for soc.
!           # mask btypmx with btypes(*) to see if this array should
!           # be accumulated somewhere.
            if ( itypbx .le. btypmx) then
               iarray = 4 + btypes(itypbx)
            else
               iarray = -1
            endif
         else
!           # ignore this array.
            iarray = -1
         endif
      return
      end

!----------------------- end colib --------------------

      subroutine gcentr(ica,nc)
c  generate group of center-interchange operators
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      dimension ica(mcu,msu,*), nc(*)

      if (ng .eq. 1) go to 56
c
c  set up identity operator
c
      do 16 is=1,ns
        icu=nc(is)
        if(icu.eq.1) go to 16
        do 12 ic=1,icu
          ica(ic,is,1)=ic
   12   continue
   16 continue
c
c        do 52(+3) ig=2,ng
c
      ig=2
   20 do 52 jg=2,ig
      igp=ig
      jgp=jg
c
c  generate new operator candidate
c
   24 do 32 is=1,ns
        icu=nc(is)
        if(icu.eq.1) go to 32
        do 28 ic=1,icu
          ica(ic,is,ng+1)=ica(ica(ic,is,igp),is,jgp)
   28   continue
   32 continue
c
c  check against list of operators
c
      do 44 kg=1,ng
        do 40 is=1,ns
          icu=nc(is)
          if(icu.eq.1) go to 40
          do 36 ic=1,icu
            if(ica(ic,is,ng+1).ne.ica(ic,is,kg)) go to 44
   36     continue
   40   continue
c       # not a new operator
        go to 48
   44 continue
c     # new operator found
      ng=ng+1
   48 if(jgp.eq.ig) go to 52
c
c  take product in other order
c
      igp=jg
      jgp=ig
      go to 24
   52 continue
      if(ig.eq.ng) go to 56
      ig=ig+1
      go to 20
c
c  set up operators for unique nuclei
c
   56 do 64 is=1,ns
      if(nc(is).eq.1) then
        do 60 ig=1,ng
          ica(1,is,ig)=1
   60   continue
      endif
   64 continue
      return
      end

      subroutine stvint(hpt,hwt)
c
c  gauss-hermite quadrature using minimum-point formula
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0)
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension hpt(*), hwt(*)
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj)/2
      imax=(npts*(npts+1))/2
      imin=imax-npts+1
      do 40 i=imin,imax
        dum=hpt(i)/t
        ptx=dum+x0
        pty=dum+y0
        ptz=dum+z0
        ax=ptx-xi
        ay=pty-yi
        az=ptz-zi
        bx=ptx-xj
        by=pty-yj
        bz=ptz-zj
        dum=hwt(i)
c       # take care to avoid (0.0d0)**0
        dumx = dum
        dumy = dum
        dumz = dum
        do 20 ipwr = 1,ni-1
          dumx = dumx*ax
          dumy = dumy*ay
          dumz = dumz*az
   20   continue
        do 30 ipwr = 1,nj-1
          dumx = dumx*bx
          dumy = dumy*by
          dumz = dumz*bz
   30   continue
        xint = xint + dumx
        yint = yint + dumy
        zint = zint + dumz
c       xint=xint+dum*ax**(ni-1)*bx**(nj-1)
c       yint=yint+dum*ay**(ni-1)*by**(nj-1)
c       zint=zint+dum*az**(ni-1)*bz**(nj-1)
   40 continue
      return
      end


      function qcomp(alpha,dfac,n,l,t,xk)
c
c  compute q(n,l) scaled by sqrt(pi)*exp(-t) to prevent overflows
c  arguments are alpha, xk, and t=xk**2/(4*alpha)
c  no restriction on the magnitude of t
c  increase dfac array to raise n, l restrictions
c
      implicit real*8 (a-h,o-z)
      parameter (am1=-1.0d0, a0=0.0d0, accpow=1.0d-14,
     1  accasy=1.0d-10, a1rtpi=0.56418958354775629d0, a1=1.0d0,
     2  a2=2.0d0, a4=4.0d0)
      dimension dfac(*)
      dimension tmin(0:8)
      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1  18.0d0,15.0d0/
c
      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c
c  use alternating series (n+l.le.22.and.l.le.10)
c
      if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=xkp * dfac(n+l+1)/
     1  ((alpha+alpha)**((n-l)/2) * sqrt(a4*alpha) * dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
   10 if(num.ne.0) then
        fnum=num
        term=term*fnum*t/(xden*xc)
        xc=xc+am1
        sum=sum+term
        num=num+2
        xden=xden+a2
        go to 10
      endif
c
      qcomp=prefac*sum
      return
c
   30 if(t.lt.tmin(min(n,8))) go to 60
c
c  use asymptotic series (arbitrary n,l)
c
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp/((alpha+alpha)*sqrt(a4*alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(abs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
c
c  use power series (n+l.le.22.and.l.le.10)
c
   60 if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) then
        prefac=prefac/sqrt(a4*alpha)
      else
        prefac=a1rtpi*prefac
      endif
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end

*deck droot
      subroutine droot(nroots,n1,xx,ff,c,s,a,rt,r,w,uf,wf)
c
c  ith root of the jth rys polynomial is returned in r(i,j)
c  with the corresponding weight factor in w(i,j).  j=1,2,...,n
c  this version uses christoffel formula for weights.
c
*mdc*if real16 vax ibm cray
*      implicit real*16 (a-h,o-z)
*      real*8 xx,uf,wf
*mdc*else
      implicit real*8 (a-h,o-z)
      real*8 xx,uf,wf
*mdc*endif
*mdc*if real16q ibm vax rs6000 fujitsu
*c *** this mdc block is not necessary since the constants will be ***
*c *** correctly converted at compile time.  should it be removed? ***
*      parameter (a0=0.0q0, a1s2=0.5q0, a1=1.0q0, a4=4.0q0)
*mdc*else
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
*mdc*endif
c     common /ffm/ x,nn,n1,k
c     dimension uf(*), wf(*)
c     dimension ff(19), c(10,10), s(10,10), a(10), rt(10),
c    &  r(9,9), w(9,9)
      dimension ff(*), c(n1,*), s(n1,*), a(*), rt(*),
     &  r(nroots,*), w(nroots,*), uf(*), wf(*)
c
      n=max(nroots,2)
      x=xx
c     n1=n+1
      nn=n+n
      call dfunc(nn,x,ff,c)
      do 15 i=1,n1
        do 10 j=1,n1
          s(i,j)=ff(i+j-1)
   10   continue
   15 continue
      call dsmit(n1,c,s,a,rt)
      do 25 i=1,n
        do 20 j=1,i
          w(i,j)=a0
          r(i,j)=a0
   20   continue
   25 continue
      wsum=ff(1)
      w(1,1)=wsum
      r(1,1)=ff(2)/wsum
      dum=sqrt(c(2,3)**2-a4*c(1,3)*c(3,3))
      r(1,2)=a1s2*(-c(2,3)-dum)/c(3,3)
      r(2,2)=a1s2*(-c(2,3)+dum)/c(3,3)
      if(n.eq.2) go to 50
      do 30 i=3,n1
        rt(i)=a1
   30 continue
      rt(1)=r(1,2)
      rt(2)=r(2,2)
      do 45 k=3,n
        k1=k+1
        do 35 i=1,k1
          a(i)=c(i,k1)
   35   continue
        call dnode(k,a,rt)
        do 40 i=1,k
          r(i,k)=rt(i)
   40   continue
   45 continue
   50 continue
      do 70 k=2,n
        do 65 i=1,k
          root=r(i,k)
          dum=a1/ff(1)
          do 60 j=1,k-1
            j1=j+1
            poly=c(j1,j1)
            do 55 m=1,j
              poly=poly*root+c(j1-m,j1)
   55       continue
            dum=dum+poly*poly
   60     continue
          w(i,k)=a1/dum
   65   continue
   70 continue
      do 75 k=1,nroots
        dum=r(k,nroots)
        uf(k)=dum/(a1-dum)
        wf(k)=w(k,nroots)
   75 continue
      return
      end
*deck dfunc
      subroutine dfunc(nn,x,ff,term)
*mdc*if real16 vax cray ibm
*      implicit real*16 (a-h,o-z)
*mdc*else
      implicit real*8 (a-h,o-z)
*mdc*endif
*mdc*if real16q ibm vax rs6000 fujitsu
*      parameter (a0=0.0q0, a1s2=0.5q0, a1=1.0q0 )
*      parameter (pie4=0.7853981633974483096156608 q0)
*mdc*else
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0)
      parameter (pie4=0.7853981633974483096156608 d0)
*mdc*endif
*mdc*if fps164
*      parameter (xsw=19.0d0)
*mdc*elseif cray
*      parameter (xsw=33.0d0)
*mdc*elseif real16q vax ibm rs6000 fujitsu
*      parameter (xsw=39.0q0)
*mdc*else
      parameter (xsw=15.0d0)
*mdc*endif
c     common /ffm/ x,nn,n1,k1
      dimension ff(*), term(*)
c
      e=a1s2*exp(-x)
      fac0=nn
      fac0=fac0+a1s2
      if(x.gt.xsw) go to 50
c
c  use power series
c
   10 fac=fac0
      term0=e/fac
      sum=term0
      ku=(x-fac0)
c     # sum increasing terms forwards
      do 20 k=1,ku
        fac=fac+a1
        term0=term0*x/fac
        sum=sum+term0
   20 continue
      i=1
      fac=fac+a1
      term(1)=term0*x/fac
      suma=sum+term(1)
      if(sum.eq.suma) go to 85
   40 i=i+1
      fac=fac+a1
      term(i)=term(i-1)*x/fac
      sum1=suma
      suma=suma+term(i)
      if(sum1-suma) 40,85,40
c
c  use asymptotic series
c
   50 sum=sqrt(pie4/x)
      if(nn.eq.0) go to 70
      fac=-a1s2
      do 60 k=1,nn
      fac=fac+a1
   60 sum=sum*fac/x
   70 i=1
      term(1)=-e/x
      suma=sum+term(1)
      if(sum.eq.suma) go to 85
      fac=fac0
      ku=(x+fac0-a1)
      do 80 i=2,ku
        fac=fac-a1
        term(i)=term(i-1)*fac/x
        sum1=suma
        suma=suma+term(i)
        if(sum1.eq.suma) go to 85
   80 continue
c     # xsw set too low. use power series.
      go to 10
c
c  sum decreasing terms backwards
c
   85 sum1=a0
c
c  11-feb-90 cray vectorization suppressed to avoid a numeric exception
c            which occured for x=5670. - rick ross/rls
cdir$ novector
      do 90 k=1,i
        sum1=sum1+term(i+1-k)
   90 continue
cdir$ vector
      ff(nn+1)=sum+sum1
c
c  use recurrence relation
c
      if(nn.eq.0) return
      do 95 k=1,nn
        fac0=fac0-a1
        ff(nn+1-k)=(e+x*ff(nn+2-k))/fac0
   95 continue
      return
      end
*deck dsmit
      subroutine dsmit(n1,c,s,v,y)
c
c  return an n1 by n1 triangular matrix c such that
c  c(t)sc=i, where i is an n1 by n1 identity matrix.
c
*mdc*if real16 vax ibm cray
*      implicit real*16 (a-h,o-z)
*mdc*else
      implicit real*8 (a-h,o-z)
*mdc*endif
*mdc*if real16q ibm vax
*c *** this mdc block is not necessary since the constants will be ***
*c *** correctly converted at compile time.  should it be removed? ***
*      parameter (a0=0.0q0, a1=1.0q0)
*mdc*else
      parameter (a0=0.0d0, a1=1.0d0)
*mdc*endif
c     common /ffm/ x,nn,n1,k1
      dimension c(n1,*), s(n1,*), v(*), y(*)
c
      do 20 i=1,n1
        do 10 j=1,i
          c(i,j)=a0
   10   continue
   20 continue
      do 80 j=1,n1
        kmax=j-1
        fac=s(j,j)
        do 30 k=1,kmax
          v(k)=a0
          y(k)=s(k,j)
   30   continue
        do 60 k=1,kmax
          dot=a0
          do 40 i=1,k
            dot=c(i,k)*y(i)+dot
   40     continue
          do 50 i=1,k
            v(i)=v(i)-dot*c(i,k)
   50     continue
          fac = fac - dot * dot
   60   continue
        fac=a1/sqrt(fac)
        c(j,j)=fac
        do 70 k=1,kmax
          c(k,j)=fac*v(k)
   70   continue
   80 continue
      return
      end
*deck dnode
      subroutine dnode(k,a,rt)
c  return in rt(i) the ith root of a polynomial
c  of order k whose mth coefficient is stored in a(m+1).
c  the initial values in rt must bracket the final values.
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1s16=6.25d-2, a1s4=0.25d0, a3s4=0.75d0)
      parameter (accrt=1.0d-12)
 
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c     common /ffm/ x,nn,n1,k
      dimension a(*), rt(*)
 
      k1=k+1
      r2=a0
      p2=a(1)
      do 95 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
        p2=p2*r2+a(k1-i)
   10 continue
      prod=p1*p2
      if(prod.ge.a0) then
        write (nlist,15) m,k
   15   format(/' dnode: root number',i4,
     +    ' was not found for polynomial of order',i4//)
        call bummer('dnode: root not found, m=',m,2)
      endif
c
      r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      dr = a1s16 * dr
      r5=r-dr
      r6=r+dr
      if(abs(delta).lt.accrt.or.r5.eq.r.or.r6.eq.r) go to 90
      if(r5.lt.r3) r5=r3
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
        p5=p5*r5+a(k1-i)
        p6=p6*r6+a(k1-i)
   40 continue
   45 prod=p5*p6
      if(prod.lt.a0) go to 30
      prod=p3*p5
      if(prod.gt.a0) go to 60
      r5=a1s4*r3+a3s4*r5
      p5=a(k1)
      do 50 i=1,k
        p5=p5*r5+a(k1-i)
   50 continue
      go to 45
   60 r6=a1s4*r4+a3s4*r6
      p6=a(k1)
      do 70 i=1,k
        p6=p6*r6+a(k1-i)
   70 continue
      go to 45
   90 rt(m)=r
   95 continue
      return
      end
