      subroutine argos_1e( WORK, LWORK, a, eta, ica, icb, icxst,
     &   iprst,  ipt,    lmnp1,  lmnv,      mcons,  nc,     ncon,
     &   nf,     nfct,   ngw,    npair,     nprir,  nrcr,   nt,
     &   ntl,    ntu,    x,      y,         z,      zet,    lcore )

c
c  read the user input and compute the 1-e integrals.
c
      implicit logical(a-z)

      DIMENSION WORK(*)
      INTEGER   LWORK

#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h" 
#include "argos_print.h" 

      integer INT_AREP
 
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
 
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
 
      integer        ipq
      common /dim21/ ipq(1000) !256->1000
 
!     # dummy:
      integer lcore
      integer
     & ica(mcu,msu,*), icb(4,24,*), icxst(2,*), iprst(*),
     & ipt(*),         lmnp1(*),    lmnv(3,*),
     & mcons(*),       nc(*),       ncon(*),    nf(*),
     & nfct(*),        ngw(*),      npair(2,*), nprir(2,mstu,*),
     & nrcr(*),        nt(*),       ntl(*),     ntu(*)
 
      real*8
     & a(*),     eta(mrcru,mconu,*), x(mcu,*), y(mcu,*),
     & z(mcu,*), zet(mconu,*)
 
!     # local:
      integer i, n, narray
      integer
     & lcr(msup),       lls(msup),      nkcrl(6,msup), nkcru(6,msup),
     & nklsl(4,msup),   nklsu(4,msup),  mcrs(msup),    nir(kaordp),
     & maords(mgcsup),  icxsv1(mgcsup), nd(mstup),     nso(mstup),
     & nsopr(mstup),    nblpr(mstup),   nopir(mstup),  ilxyz(3,mstup),
     & la(mrup,kaordp), ncr(mcrup),     nct(msfup),    mgcs(msfup),
     & lb(msfrup),      ms(msfrup),     mnl(msfrup),   mau(mnrup),
     & mics(mnrup),     mjcs(mnrup),    icxast(mnsfup),
     & icxsv2(mgcsup,mgcsup),           idp(mstup,mstup,mstup),
     & ihwsp,           igh,            lmn1,          nn,
     & lmn1u,           ncru,           lproju,        ndfac,
     & lmax,            l1max2,         lmnpwr,        l2m1
 
      real*8 chg(msup), zcr(mcrup), ccr(mcrup), eps
 
      character*3 ityp(mstup), mtype(msup)
      integer  atebyt, forbyt
      external atebyt, forbyt
 
      kaords = kaordp
      mconsu = mconsp
      mgcsu  = mgcsup
      mru    = mrup
      mcsu   = mcsup
      mctu   = mctup
      mcru   = mcrup
      msfu   = msfup
      mgu    = mgup
      msfru  = msfrup
      mnru   = mnrup
 
!     # lcore = total available workspace on entry to this routine.
      lcore  = mblu
 
      do i = 1, 3
        inam(i) = i
      enddo
 
!     # set up ipq(*).
      do i = 1, 1000  ! For temporary 255-> 1000 yc
        ipq(i) = (i * (i - 1) ) / 2
      enddo
 
!     # reserve space for c(mctu,mcsu,mgcsu) at the end of a(1:mblu).
!     # afterwards, a(1:mcxu) is the remaining available workspace.
!     # a(*)--> 1:...55:c(1:csize)
 
      ipt(55) = mblu + 1 - atebyt( mctu * mcsu * mgcsu )
 
      mcxu = ipt(55) - ipt(1)
 
      if ( mcxu .lt. 0 ) then
        call bummer('seg1mn: before syminp() mcxu=',mcxu,faterr)
      endif
 
!     # input processing.
      call RECP_INPUT( a,      a(ipt(55)), ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,       la,     lb,     lcr,
     &   lls,    lmnp1,  lmnv,   lmn1u,      lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,         mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,        nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,      nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,          y,      z,      zcr,
     &   zet,    lcore,  N1MAX,  L1REC)
 
!     # allocate space for il(*), and compute the
!     # symmetry coefficient products.
!     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(*)...,
!     #          55:c(*)
 
      ipt(3) = ipt(2) + forbyt( nnbft )
 
!     # mcxu is the remaining space.
      mcxu = ipt(55) - ipt(3)
 
      call argos_socfpd( a(ipt(55)), a(ipt(3)), icxast, icxst,
     &   icxsv1, icxsv2, idp,        a(ipt(2)), iprst,  la,
     &   lb,     maords, mau,        mcons,     mgcs,   mics,
     &   mjcs,   nblpr,  nc,         nd,        nf,     nir,
     &   npair,  nprir,  nrcr,       nsopr,     nt,     ntl,
     &   ntu  )
 
!     # compute the one-electron integrals.

!     # on return from socfpd(*), cx(1:mcxu) was computed.
!     # allocate space.
!     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
!     #          4:array(1:mpru,1:narray), 5:iarray(1:iamax,1:narray),
!     #          6:buffer(1:l1rec), 7:values(1:n1max),
!     #          8:labels(1:2,1:n1max), 9:hpt(1:ihwsp), 10:hwt(1:ihwsp),
!     #          11:dfac(1:ndfac), 12:binom(1:xyzdim), 13:lmf(1:l1max2),
!     #          14:lml(1:l1max2), 15:lmx(1:lmnpwr), 16:lmy(1:lmnpwr),
!     #          17:lmz(1:lmnpwr), 18:zlm(1:lmnpwr),
!     #          19:flmtx(1:3,1:lmnu2), 20:mc(1:3,1:l2m1),
!     #          21:mr(1:3,1:l2m1), 22:h2(1:nblu)
 
      if ( inam(nu) .eq. 5 ) then
        narray = 3 !    spin-orbit integral (x,y,z components)
      else
        narray = 1 ! no spin-orbit integral
      endif
 
      ihwsp = (lmn1u*(lmn1u+1))/2
      ipt(4)  = ipt(3)  + atebyt( mcxu )
      ipt(5)  = ipt(4)  + atebyt( narray * mpru ) !           :array
      ipt(9)  = ipt(5)  + forbyt( narray * mpru ) !           :iarray
!     ipt(7)  = ipt(6)  + atebyt( l1rec )         !omit ipt(6):buffer
!     ipt(8)  = ipt(7)  + atebyt( n1max )         !omit ipt(7):values
!     ipt(9)  = ipt(8)  + forbyt( 2 * n1max )     !omit ipt(8):labels
      ipt(10) = ipt(9)  + atebyt( ihwsp )
      ipt(22) = ipt(10) + atebyt( ihwsp )
 
      mblu = lcore - ipt(22) + 1
 
!     # compute gauss-hermite points and weights for s, t, v integrals.
!     # 10**(-12) needed to get 6-point values to sufficient accuracy.
      eps = 1.0d-12
      igh = 0
      do 40 lmn1 = 1, lmn1u
        call hermit(lmn1,a(ipt(9)+igh),a(ipt(10)+igh),eps)
        igh = igh + lmn1
   40 continue
 
      write(nlist,*)

      INT_AREP = 0
      IF (NU.eq.4) INT_AREP = 1
 
      do 60 n = 1, nu
 
        if( n.eq.4 ) then
!         # allocate space.
          ipt(10) = ipt(9)  + atebyt( 35 )
          ipt(11) = ipt(10) + atebyt( 35 )
          ndfac = max(4*lmn1u+2*lproju-3,6*lproju+3,4*lmn1u-1,
     &      2*lmn1u+2*lproju+1,4,ncru+4*lmn1u+2*lproju-1)
          ipt(12) = ipt(11) + atebyt( ndfac )
          ipt(13) = ipt(12) + atebyt( (lmn1u*(lmn1u+1))/2 )
          lmax =  max(1,lmn1u-1 + max(lmn1u-1,lproju))
          l1max2 = (lmax+1)**2
          ipt(14) = ipt(13) + forbyt( l1max2 )
          ipt(15) = ipt(14) + forbyt( l1max2 )
          lmnpwr = (((lmax*(lmax+2)*(lmax+4))/3)*(lmax+3)+
     &               (lmax+2)**2*(lmax+4))/16
          ipt(16) = ipt(15) + forbyt( lmnpwr )
          ipt(17) = ipt(16) + forbyt( lmnpwr )
          ipt(18) = ipt(17) + forbyt( lmnpwr )
          ipt(19) = ipt(18) + atebyt( lmnpwr )
          ipt(20) = ipt(19) + atebyt( 3*lproju**2 )
          l2m1 = 2*lproju-1
          ipt(21) = ipt(20) + forbyt( 3*l2m1 )
          ipt(22) = ipt(21) + forbyt( 3*l2m1 )
          call cortab(a(ipt(12)),a(ipt(11)),eps,a(ipt(19)),a(ipt(9)),
     &      a(ipt(10)),a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),
     &      a(ipt(17)),lmax,lmn1u,lproju,a(ipt(20)),a(ipt(21)),ndfac,
     &      a(ipt(18)))
        endif
 
        nnam  = inam(n)
        call argos_oneint( 
     &        a,       ccr,     chg,        a(ipt(3)),
     &        eta,       a(ipt(22)), ica,
     &        icb,       icxast,     icxst,     a(ipt(2)),  ilxyz,
     &        iprst,     ipt,        lcr,       lls,        lmnp1,
     &        lproju,    mcons,      mcrs,      nc,         ncon,
     &        ncr,       nf,         nfct,      ngw,        nkcrl,
     &        nkcru,     nklsl,      nklsu,     nopir,      npair,
     &        nprir,     nrcr,       nt,        ntl,        ntu,
     &        x,         y,          z,         zcr,        zet,
     &        a(ipt(4)), a(ipt(5)),   
     &        nso, INT_AREP, N1MAX,  L1REC )
 
   60 continue
      write(nlist,*)
      return
      end


*deck l2opxv
      subroutine l2opxv( lmn, v, itl, itu, lmnv, lenscr, scr, eval,
     & leig )
c
c  compute the matrix vector product, (L^2) * v, in an unnormalized
c  cartesian basis, and determine if v(*) is an eigenvector of the
c  total angular momentum operator.
c
c  input:
c  lmn = (l + m + n) where l, m, and n are the exponents of x, y, and z
c        respectively.  lmn=0 for cartesian "s" functions, 1 for
c        cartesian "p" functions, 2 for cartesian "d" functions, etc.
c  v(1:xyzdim) = input vector.  the elements are coefficients of the
c              cartesian basis functions (x**l * y**m * z**n).
c              where xyzdim = ((lmn + 1) * (lmn + 2))/2
c  lmnv(1:3,1:lmnvmx) = cartesian exponents of each basis function.
c  lenscr = scratch vector length.
c  scr(1:lenscr) = scratch work vector.  This must be at least as large
c                  as the cartesian subspace dimension.
c                  lenscr >= xyzdim
c
c  output:
c  eval = vector expectation value = <v| L^2 |v> / <v|v>
c  leig = eigenvector return code.
c       = -4 for vector norm error.
c       = -3 for lmnv(*,*) inconsistency.
c       = -2 for lenscr error.
c       = -1 if v(*) is not an eigenvector of L^2
c       = principal quantum number if v(*) is an eigenvector.
c         0 for s-type vectors, 1 for p-type vectors, etc.
c         the eigenvalue is given by (leig*(leig+1)) to within numerical
c         precision.
c  28-oct-96 account for nonorthogonality of cartesian gaussians (rmp)
c  11-jun-91 written by ron shepard with suggestions from r. m. pitzer.
c
      implicit logical(a-z)
c
c     # dummy:
c     integer lmn, lenv, lenscr, leig
      integer itl, itu, lmn, lenscr, leig
      integer lmnv(3,*)
      real*8 eval
      real*8 v(*), scr(*)
c
c     # local:
      integer l, m, n, i, j, it, jt, xyzdim, ll, mm, nn, p, q
      real*8 vi, rnorm, vnorm2
c
c     # small is used to determine if a vector is close enough to be
c     #       called an eigenvector.  this should be about 1000x the
c     #       machine epsilon to avoid sqrt() precision problems.
c
      real*8    zero,     one,     two,     fourth,        small
      parameter(zero=0d0, one=1d0, two=2d0, fourth=0.25d0, small=1d-10)
 
#include "argos_error.h"
 
      integer  iodfac
      external iodfac
c
c     # the [l,m,n] basis function is assigned the internal index:
c     #  (lx*(lx+1))/2 + n + 1   with lx=(lmn-l)=(m+n)
c     # this convention is inconsistent with the rest of this program,
c     # but it is used anyway since it is convenient and since the
c     # computed results are local.
c
      integer llp1
      llp1(l) = (l * (l + 1)) / 2 + 1
c
      leig = -1
      eval = -one
c
c     # check to make sure scr(*) is large enough.
      xyzdim = llp1( lmn) + lmn
      if ( xyzdim .gt. lenscr ) then
         call bummer('l2opxv: need larger scr(*), (xyzdim-lenscr)=',
     &    (xyzdim-lenscr), wrnerr )
         leig = -2
         return
      endif
c
c     # initialize the scratch vector.
      call wzero(xyzdim,scr,1)
c
c     # compute the matrix-vector product and the vector norm.
c     # matrix elements of lop(*,*) are computed from the operator
c     # definition in cartesian space.
c
      vnorm2 = zero
      do 30 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         vi = v(i)
         do 20 jt = itl, it-1
            j = jt - (itl - 1)
            vnorm2 = vnorm2 + (two * vi) * v(j) *
     &       iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
   20    continue
         vnorm2 = vnorm2 + vi * vi * iodfac(l,l,m,m,n,n)
c
         if ( (l+m+n) .ne. lmn ) then
c           # inconsistent exponents in the basis function.
            call bummer('l2opxv: (l+m+n-lmn)=', (l+m+n-lmn), wrnerr )
            leig = -3
            return
         endif
c
c        # L^2 is sparse in this representation.  each vi contributes
c        # to, at most, only 7 elements in the matrix-vector product.
c
         if ( l .ge. 2 ) then
c           # include -l*(l-1) * ( [l-2,m,n+2] + [l-2,m+2,n] ) terms.
            ll = l * (l - 1)
            p = llp1( lmn - l + 2 ) + n
            scr(p) = scr(p) - vi * (ll)
            p = p + 2
            scr(p) = scr(p) - vi * (ll)
         endif
c
         if ( m .ge. 2 ) then
c           # include -m*(m-1) * ( [l,m-2,n+2] + [l+2,m-2,n] ) terms.
            mm = m * (m - 1)
            p = llp1( lmn - l ) + n + 2
            scr(p) = scr(p) - vi * (mm)
            p = llp1( lmn - l - 2 ) + n
            scr(p) = scr(p) - vi * (mm)
         endif
c
         if ( n .ge. 2 ) then
c           # include -n*(n-1) * ( [l,m+2,n-2] + [l+2,m,n-2] ) terms.
            nn = n * (n - 1)
            p = llp1( lmn - l ) + n - 2
            scr(p) = scr(p) - vi * (nn)
            p = llp1( lmn - l - 2 ) + n - 2
            scr(p) = scr(p) - vi * (nn)
         endif
c
c        # include the 2*(l*m+l*n+m*n+l+m+n)*[l,m,n] diagonal term.
c
         p = llp1( lmn - l ) + n
         scr(p) = scr(p) + vi * (2 * (l * (m + n) + m * n + l + m + n))
   30 continue
c
      if ( vnorm2 .le. small ) then
         call bummer('l2mxv: small vector norm, lnm=', lmn, wrnerr)
         leig = -4
         return
      endif
c
c     # compute the expectation value.
c
      eval  = zero
      do 50 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         p = llp1( lmn - l ) + n
         do 40 jt = itl, itu
            j = jt - (itl - 1)
            eval = eval + v(j) * scr(p) *
     &                  iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
   40    continue
   50 continue
      eval = eval / vnorm2
c
c     # compute the residual norm.
c     #
      rnorm = zero
c
      do 70 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         p = llp1( lmn - l ) + n
         do 60 jt = itl, it-1
            j = jt - (itl - 1)
            q = llp1( lmn - lmnv(1,jt) ) + lmnv(3,jt)
            rnorm = rnorm + two * ( scr(p) - eval * v(i) ) *
     &                            ( scr(q) - eval * v(j) ) *
     &                iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
   60    continue
        rnorm = rnorm + (scr(p) - eval * v(i))**2 * iodfac(l,l,m,m,n,n)
   70 continue
c
c     # normalize w.r.t. |v|=1.
      rnorm = rnorm / vnorm2
c
      rnorm = sqrt( rnorm )
      if ( rnorm .gt. small ) then
c
c        # v(*) is not an eigenvector.
c
         leig = -1
      else
c
c        # v(*) is an eigenvector.
c
c        # determine leig such that eval = leig*(leig+1)
c
c        # the following assignment should truncate to the
c        # next smaller integer value...
         leig = ( sqrt( eval + fourth ) )
      endif
c
      return
      end


      subroutine argos_socfpd( c,     cx,         icxast, icxst,
     &   icxsv1, icxsv2, idp,         il,         iprst,  la,
     &   lb,     maords, mau,         mcons,      mgcs,   mics,
     &   mjcs,   nblpr,  nc,          nd,         nf,     nir,
     &   npair,  nprir,  nrcr,        nsopr,      nt,     ntl,
     &   ntu  )

      implicit logical(a-z)
 
!     # il(*) packing factor. see also oneint() and wtint2().
      integer    ilfact
      parameter( ilfact=1024 )

#include "argos_files.h"
#include "argos_error.h"
#include "argos_print.h"

      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
 
      integer        ipq
      common /dim21/ ipq(1000) !256->1000
 
      integer
     & icxast(*),        icxst(2,*),   icxsv1(*),    icxsv2(mgcsu,*),
     & idp(mstu,mstu,*), il(nnbft),    iprst(*),     la(mru,*),
     & lb(*),            maords(*),    mau(*),       mcons(*),
     & mgcs(*),          mics(*),      mjcs(*),      nblpr(*),
     & nc(*),            nd(*),        nf(*),        nir(*),
     & npair(2,*),       nrcr(*),      nsopr(*),     nprir(2,mstu,*),
     & nt(*),            ntl(*),       ntu(*)
      real*8 c(mctu,mcsu,*), cx(*)
c
c     # local:
      integer mcxu1, icx, ipair, ijsf, isf, isfr, igcs, jgcs, is,
     & ifu, icu, isfis, isfris, jsf, jsfr, js, jsfjs, jsfrjs, jfu,
     & jcu, if, isfrif, iaords, iru, ircru, itl, itu, jf, jaords,
     & jtl, jsfrjf, jtu, jru, ist, jrcru, ircr, isfrib, jrcr,
     & jsfrjb, iesfb, npr, ipr, ics, ir, lai, iau, jcs, jr, laj,
     & isv1, icxu, icxa, isv2, ict, ic, ictic, jct, jc, jctjc, it,
     & jt, ia, icxl, mnsf, npru
      logical ea, easfb, es, esf, esfb, esfbc, msfbct
c
      real*8     a0
      parameter (a0=0.0d0)

      IF (PRINT_LEVEL .GE. 5)
     &   PRINT *,' * ECPINT : entering... ARGOS_SOCFPD'

c     # check the il(*) packing factor.
      if ( nbft .ge. ilfact ) then
        call bummer('socfpd: larger ilfact required, nbft=',
     &   nbft, faterr )
      endif
c
      ea = inam(nu).eq.5
c
c     # mcxu1 is the original mcxu offset.
      mcxu1 = mcxu
      icx   = 0
      ipair = 0
c
c     # ijsf = one dimensionalization of is, js, if, jf
c     # isfr = one dimensionalization of is, if, ir
      ijsf = 0
      isf  = 0
      isfr = 0
      do 12 igcs = 1, ngcs
        icxsv1(igcs) = 0
        do 10 jgcs = 1, ngcs
          icxsv2(jgcs,igcs) = 0
10      continue
12    continue
      do 340 is = 1, ns
        ifu = nf(is)
        if(ifu.eq.0) go to 340
        icu = nc(is)
        isfis = isf
        isfris = isfr
        jsf = 0
        jsfr = 0
        do 330 js = 1,is
          jsfjs = jsf
          jsfrjs = jsfr
          es = is.eq.js
          if(.not.es) then
            jfu = nf(js)
            if(jfu.eq.0) go to 330
            jcu = nc(js)
          endif
          isf = isfis
          isfr = isfris
          do 324 if = 1,ifu
            isf = isf+1
            isfrif = isfr
            igcs = mgcs(isf)
            iaords = maords(igcs)
            iru = nir(iaords)
            ircru = nrcr(mcons(isf))
            itl = ntl(isf)
            itu = ntu(isf)
            jsf = jsfjs
            jsfr = jsfrjs
            if(es) jfu = if
            do 322 jf = 1,jfu
              jsf = jsf+1
              jgcs = mgcs(jsf)
              jaords = maords(jgcs)
              jtl = ntl(jsf)
              jsfrjf = jsfr
              ijsf = ijsf+1
              esf = isf.eq.jsf
              if(esf) then
                if(ircru.ge.2) then
                  jtu = ntu(jsf)
                  jru = nir(jaords)
                else
                  npair(2,ijsf) = 0
                  icxst(2,ijsf) = 0
                  do 20 ist = 1,nst
                    nprir(2,ist,ijsf) = 0
20                continue
                endif
              else
                jrcru = nrcr(mcons(jsf))
                jtu   = ntu(jsf)
                jru   = nir(jaords)
              endif
              iprst(ijsf) = ipair
              isfr = isfrif
              do 320 ircr = 1,ircru
                isfrib = isfr
                jsfr = jsfrjf
                if(esf) jrcru = ircr
                do 310 jrcr = 1,jrcru
                  jsfrjb = jsfr
                  esfb = esf.and.ircr.eq.jrcr
                  easfb = ea.and.esfb
                  if(esfb) then
                    iesfb = 1
                    npru = (iru * (iru + 1)) / 2
                  else
                    iesfb = 2
                    npru = iru * jru
                  endif
c
                  if ( npru .gt. mnru ) then
                    call bummer('change mnrup (one place) to ',
     &               npru,faterr)
                  endif
c
                  npr = 0
                  do 75 ist = 1, nst
                    ipr = npr
                    ics = 0
                    isfr = isfrib
                    do 70 ir = 1, iru
                      isfr = isfr + 1
                      lai = la(ir,iaords)
                      iau = nd(lai)
                      if(esfb) jru = ir
                      jsfr = jsfrjb
                      jcs = 0
                      do 65 jr = 1, jru
                        jsfr = jsfr + 1
                        laj = la(jr,jaords)
                        if ( idp(lai,laj,ist) .ne. 0 ) then
                          ipair = ipair + 1
                          il(ipair) = (nsopr(lai) + lb(isfr))
     &                     * ilfact + (nsopr(laj) + lb(jsfr))
c
                          npr = npr + 1
                          mau(npr)  = iau
                          mics(npr) = ics
                          mjcs(npr) = jcs
                        endif
                        jcs = jcs + nd(laj)
65                    continue
                      ics = ics + iau
70                  continue
                    nprir(iesfb,ist,ijsf) = npr - ipr
75                continue
                  npair(iesfb,ijsf) = npr
c
                  if ( npr .eq. 0 ) go to 310
c
                  if ( esfb ) then
                    isv1 = icxsv1(igcs)
                    if ( isv1 .ne. 0 ) then
                      icxst(1,ijsf) = icxst(1,isv1)
                      if(ea) icxast(ijsf) = icxast(isv1)
                      go to 310
                    endif
                    icxsv1(igcs) = ijsf
                    icxst(1,ijsf) = icx
                    icxu = icx+(npr*ipq(icu*nt(isf)+1))
                    if ( ea ) then
                      mcxu = mcxu-(npr*ipq(icu*nt(isf)+1))
                      icxa = mcxu
                      icxast(ijsf) = icxa
                    endif
                  else
                    isv2 = icxsv2(jgcs,igcs)
                    if ( isv2 .ne. 0 ) then
                      icxst(2,ijsf) = icxst(2,isv2)
                      go to 310
                    endif
                    icxsv2(jgcs,igcs) = ijsf
                    icxst(2,ijsf) = icx
                    icxu = icx+npr*icu*nt(isf)*jcu*nt(jsf)
                  endif
                  if ( icxu .gt. mcxu ) then
                    call bummer('socfpd: (icxu-mcxu)=',
     &               (icxu-mcxu), faterr )
                  endif
c
c                 # compute new cx block
c
                  ict = 0
                  do 298 ic = 1, icu
                    ictic = ict
                    if(esfb) jcu = ic
                    jct = 0
                    do 296 jc = 1, jcu
                      jctjc = jct
                      esfbc = esfb .and. ic .eq. jc
                      ict = ictic
                      do 294 it = itl, itu
                        ict = ict + 1
                        if(esfbc) jtu = it
                        jct = jctjc
                        do 292 jt = jtl, jtu
                          jct = jct + 1
                          msfbct =  .not. esfb .or.
     &                     (ic.eq.jc .and. it.eq.jt)
                          do 290 ipr = 1, npr
                            icx = icx + 1
                            cx(icx) = a0
                            if(easfb) then
                              icxa = icxa + 1
                              cx(icxa) = a0
                            endif
                            iau = mau(ipr)
                            ics = mics(ipr)
                            jcs = mjcs(ipr)
                            do 280 ia = 1, iau
                              ics = ics + 1
                              jcs = jcs + 1
                              cx(icx) = cx(icx)
     &                         + (c(ict,ics,igcs) * c(jct,jcs,jgcs))
                              if(easfb) cx(icxa) = cx(icxa)
     &                         + (c(ict,ics,igcs) * c(jct,jcs,jgcs))
                              if(msfbct) go to 280
                              cx(icx) = cx(icx)
     &                         + (c(jct,ics,igcs) * c(ict,jcs,jgcs))
                              if(easfb) cx(icxa) = cx(icxa)
     &                         - (c(jct,ics,igcs) * c(ict,jcs,jgcs))
280                         continue
290                       continue
292                     continue
294                   continue
296                 continue
298               continue
310             continue
320           continue
322         continue
324       continue
330     continue
340   continue
c
      if ( ea ) then
c
c       # move antisymmetric coefficient products down.
c
c       #       offsets: 0           mcxu1  mblu
c       #                |            |     |
c       # before: a(*)-->.............c(*)
c       #
c       #                0      mcxu mcxu1  mblu
c       #                |       |    |     |
c       # now:    a(*)-->cx(*)...x(*),c(*)
c       #
c       #                0   mcxu2 mcxu  (mblu+mcxu)
c       #                |     |   |        |
c       # after:  a(*)-->cx(*),x(*),........
c
c       # on return, mcxu2 is the offset for the s.o. coefficients.
        mcxu2 = icx
c
        icxl  = icx + 1
        icxu  = mcxu - icx
        mcxu  = mcxu1 - icxu
        do icx = icxl, mcxu
          cx(icx) = cx(icx+icxu)
        enddo

c       # make pointers relative to mcxu2 instead of 0.
        do mnsf = 1, ijsf
          icxast(mnsf) = icxast(mnsf) - icxu
        enddo

c       # on return, mcxu and mblu divide the workspace as:
c       # a(*)-->[cx//x](1:mcxu),w(1:mblu)
        mblu = mblu - mcxu
        IF (PRINT_LEVEL .GE. 3) write (nlist,6010) mcxu, mcxu2, mblu
6010    format(/' socfpd: mcxu=',i9,' mcxu2=',i9,' left=',i9)

      else
c       # on return, mcxu and mblu divide the workspace as:
c       # a(*)-->cx(1:mcxu),w(1:mblu)
        mcxu = icx
        mblu = mblu - mcxu
        IF (PRINT_LEVEL .GE. 3) write (nlist,6020) mcxu, mblu
6020    format(/' socfpd: mcxu=',i9,' left=',i9)
      endif
 

      IF (PRINT_LEVEL .GE. 5)
     &   PRINT *,' * ECPINT : exiting...  ARGOS_SOCFPD'

      return
      end

*deck hermit
      subroutine hermit(nn,x,a,eps)
c          calculates the zeros  x(i)  of the nn-th order
c        hermite polynomial.  the largest zero will be
c        stored in x(1).  also calculates the corresponding
c        coefficients  a(i)  of the nn-th order gauss-hermite
c        quadrature formula of degree 2*nn-1.  the factor of
c        sqrt(pi) has been removed from the a(i).
c
c  a. h. stroud & d. secrest, gaussian quadrature formulas,
c  prentice-hall, 1966
c
      implicit real*8 (a-h,o-z)
      parameter (sixth = 1.0d0/6.0d0)
      dimension x(*), a(*)

      IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : HERMIT ENTER * '
c
      fn = (nn)
      n1 = nn - 1
      n2 = (nn+1)/2
      cc = 1.0d0
      s  = 0.0d0
      do 15 i=1,n1
        s  = s + 0.5d0
        cc = s*cc
   15 continue
      s  = (2.0d0*fn+1.0d0)**sixth
      do 30 i=1,n2
        if( i.eq.1 ) then
c         # largest zero
          xt = s**3 - 1.85575d0/s
        elseif( i.eq.2 ) then
c         # second zero
          xt = xt - 1.14d0*fn**0.426d0/xt
        elseif( i.eq.3 ) then
c         # third zero
          xt = 1.86d0*xt - 0.86d0*x(1)
        elseif( i.eq.4 ) then
c         # fourth zero
          xt = 1.91d0*xt - 0.91d0*x(2)
        else
c         # all other zeros
          xt = 2.0d0*xt - x(i-2)
        endif
c
        call hroot(xt,nn,dpn,pn1,eps)
        x(i) = xt
        a(i) = cc/dpn/pn1
c       write (6,'(2i4,2d25.17)') nn, i, xt, a(i)
        ni = nn-i+1
        x(ni) = -xt
        a(ni) = a(i)
   30 continue

      IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : HERMIT EXIT  * '

      return
      end
*deck hroot
      subroutine hroot(x,nn,dpn,pn1,eps)
c     improves the approximate root  x
c     in addition we also obtain
c         dpn = derivative of h(n) at x
c         pn1 = value of h(n-1) at x

      implicit real*8 (a-h,o-z)

c     # iter = 5 sufficient for 8-byte accuracy up to nn = 7
      do 14 iter=1,10
        call hrecur_argos(p,dp,pn1,x,nn)
        d  = p/dp
        x  = x - d
        if( abs(d).le.eps ) go to 16
   14 continue
   16 dpn = dp
      return
      end

*deck hrecur_argos
      subroutine hrecur_argos(pn,dpn,pn1,x,nn)
c
      implicit real*8 (a-h,o-z)
c
      p1 = 1.0d0
      p  = x
      dp1 = 0.0d0
      dp  = 1.0d0
      do 20 j=2,nn
        fj = (j)
        fj2 = 0.5d0*(fj-1.0d0)
        q  = x*p - fj2*p1
        dq = x*dp + p - fj2*dp1
        p1 = p
        p  = q
        dp1 = dp
        dp  = dq
   20 continue
      pn  = p
      dpn = dp
      pn1 = p1
      return
      end

      subroutine cortab(binom,dfac,eps,flmtx,hpt,hwt,lmf,lml,lmx,lmy,
     &  lmz,lmax,lmn1u,lproju,mc,mr,ndfac,zlm)

c     # tables for core potential and spin-orbit integrals.

      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0, a3=3.0d0)
      dimension binom(*),dfac(*),flmtx(3,*),hpt(*),hwt(*),lmf(*),
     &  lml(*),lmx(*),lmy(*),lmz(*),mc(3,*),mr(3,*),zlm(*)
c
c     # compute gauss-hermite points and weights for c, z integrals.
      igh =1
      nn = 5
      do 120 i = 1, 3
        call hermit(nn,hpt(igh),hwt(igh),eps)
        igh = igh + nn
        nn = 2*nn
  120 continue
c     # compute double factorials.
      dfac(1)=a1
      dfac(2)=a1
      fi=a0
      do 130 i=1,ndfac-2
        fi=fi+a1
        dfac(i+2)=fi * dfac(i)
  130 continue
c     # compute binomial coefficients.
      inew=1
      binom(1)=a1
      do 150 j=1,lmn1u-1
        inew=inew+1
        binom(inew)=a1
        do 140 i=1,j-1
          inew=inew+1
          binom(inew)=((j-i+1)*binom(inew-1))/(i)
  140   continue
        inew=inew+1
        binom(inew)=a1
  150 continue
c     # compute tables by recursion for real spherical harmonics.  they
c     # are indexed by l, m and sigma.  the sequence number of the
c     # harmonic with quantum numbers l, m and sigma is given by
c     #            l**2+2*m+1-sigma
c     # lmf(index) and lml(index) hold the positions of the first and
c     # last terms of the harmonic in the lmx, lmy, lmz, and zlm arrays.
c     # the harmonics with angular momentum l are generated from those
c     # with angular momenta l-1 and l-2.
c     # for m = 0,1,2,...,l-1, the recursion relation
c     z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
c                  sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
c     # is used.
c     # for m = l, the recursion relation
c     x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) =
c                  sqrt((2*l))/((2*l+1)))*Z(l,l,s)
c     # is used.
c     # l=0
      lmf(1) = 1
      lml(1) = 1
      lmx(1) = 0
      lmy(1) = 0
      lmz(1) = 0
      zlm(1) = a1
c     # l=1
      lmf(2) = 2
      lml(2) = 2
      lmx(2) = 0
      lmy(2) = 0
      lmz(2) = 1
      zlm(2) = sqrt(a3)
      lmf(3) = 3
      lml(3) = 3
      lmx(3) = 0
      lmy(3) = 1
      lmz(3) = 0
      zlm(3) = zlm(2)
      lmf(4) = 4
      lml(4) = 4
      lmx(4) = 1
      lmy(4) = 0
      lmz(4) = 0
      zlm(4) = zlm(2)
      nterm=4
      do 270 lang=2,lmax
        do 240 mang=0,lang-1
          anum = ((2*lang-1)*(2*lang+1))
          aden = ((lang-mang)*(lang+mang))
          coef1 = sqrt(anum/aden)
          anum = ((lang+mang-1)*(lang-mang-1)*(2*lang+1))
          aden = (2*lang-3)*aden
          coef2 = sqrt(anum/aden)
          nsigma=min(1,mang)
          do 230 isigma=nsigma,0,-1
            indexh=lang**2+2*mang+1-isigma
            lone=lang-1
            ltwo=lang-2
            ione=lone**2+2*mang+1-isigma
            itwo=ltwo**2+2*mang+1-isigma
            lmf(indexh)=lml(indexh-1)+1
            lml(indexh)=lml(indexh-1)
            nxy=(mang-isigma+2)/2
            iu=lmf(ione)+nxy-1
            do 200 i=lmf(ione),iu
              lml(indexh)=lml(indexh)+1
              j=lml(indexh)
              lmx(j)=lmx(i)
              lmy(j)=lmy(i)
              lmz(j)=lmz(i)+1
              zlm(j)=zlm(i)*coef1
              nterm=nterm+1
  200       end do
            if(ltwo.ge.mang) then
              il=iu+1
              do 210 i=il,lml(ione)
                lml(indexh)=lml(indexh)+1
                j=lml(indexh)
                k=lmf(itwo)+i-il
                lmx(j)=lmx(k)
                lmy(j)=lmy(k)
                lmz(j)=lmz(k)
                zlm(j)=zlm(i)*coef1-zlm(k)*coef2
                nterm=nterm+1
  210         enddo
              il=lml(itwo)-nxy+1
              if(mod(lang-mang,2).eq.0) then
                do 220 i=il,lml(itwo)
                  lml(indexh)=lml(indexh)+1
                  j=lml(indexh)
                  lmx(j)=lmx(i)
                  lmy(j)=lmy(i)
                  lmz(j)=lmz(i)
                  zlm(j)=-zlm(i)*coef2
                  nterm=nterm+1
  220           end do
              endif
            endif
  230     end do
  240   end do
        anum = (2*lang+1)
        aden = (2*lang)
        coef = sqrt(anum/aden)
        mang=lang
        isigma=1
        indexh=lang**2+2*mang+1-isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lml(indexh-1)
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonic (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        k=lmf(isigm)
        do 250 i=lmf(isig),lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)+zlm(k))*coef
          k=k+1
          nterm=nterm+1
  250   enddo
        if(mod(mang,2).eq.1) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=zlm(k)*coef
          nterm=nterm+1
        endif
        isigma=0
        indexh=lang**2+2*mang+1-isigma
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonuc (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lmf(indexh)
        j=lml(indexh)
        i=lmf(isig)
        lmx(j)=lmx(i)+1
        lmy(j)=lmy(i)
        lmz(j)=lmz(i)
        zlm(j)=zlm(i)*coef
        nterm=nterm+1
        k=lmf(isigm)
        do 260 i=lmf(isig)+1,lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)-zlm(k))*coef
          k=k+1
          nterm=nterm+1
  260   enddo
        if(mod(mang,2).eq.0) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          k=lml(isigm)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=-zlm(k)*coef
          nterm=nterm+1
        endif
  270 end do
      ixy = 0
      iz = 0
      do 300 lang=1,lproju
        do 290 mang=0,lang-1
          nsigma=min(1,mang)
          ndelta=max(0,1-mang)
          anum = ((lang-mang)*(lang+mang+1))
          aden = (2*(2-ndelta))
          coef=sqrt(anum/aden)
          do 280 isigma=nsigma,0,-1
            isign=2*isigma-1
            ixy = ixy+1
            flmtx(1,ixy) = (isign)*coef
            flmtx(2,ixy) = coef
            if(mang.ne.0) then
              iz=iz+1
              flmtx(3,iz) = -(mang*isigma)
            endif
  280     enddo
  290   enddo
        iz=iz+1
        flmtx(3,iz) = -(lang)
  300 enddo
c     # column and row indices for angular momentum matrix elements.
      iadd = 1
      do 310 i=1,2*lproju-1
        mc(1,i) = i
        mc(2,i) = i
        mc(3,i) = i+1
        mr(1,i) = i+iadd
        mr(2,i) = i+2
        mr(3,i) = i+2
        iadd = 4-iadd
  310 continue
      return
      end

! iodfac -> FN_IODFAC 
      function iodfac(l1,l2,m1,m2,n1,n2)
c
c  compute the product of odd number factorials which gives the
c  one-center overlaps of cartesian gaussian aos.
c
      if(mod(l1+l2,2).eq.1 .or.
     &   mod(m1+m2,2).eq.1 .or.
     &   mod(n1+n2,2).eq.1) then
        iodfac = 0
      else
        iprd = 1
        do 20 i=3,l1+l2-1,2
          iprd = i*iprd
   20   continue
        do 30 i=3,m1+m2-1,2
          iprd = i*iprd
   30   continue
        do 40 i=3,n1+n2-1,2
          iprd = i*iprd
   40   continue
        iodfac = iprd
      endif
      return
      end
