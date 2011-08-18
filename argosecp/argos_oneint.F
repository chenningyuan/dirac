      subroutine argos_oneint( 
     &       a,          ccr,        chg,        cx,
     &       eta,        h2,         ica,
     &       icb,        icxast,     icxst,      il,         ilxyz,
     &       iprst,      ipt,        lcr,        lls,        lmnp1,
     &       lproju,     mcons,      mcrs,       nc,         ncon,
     &       ncr,        nf,         nfct,       ngw,        nkcrl,
     &       nkcru,      nklsl,      nklsu,      nopir,      npair,
     &       nprir,      nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zcr,        zet,
     &       array,      iarray,           
     &       n_bfn_sym,  INT_AREP,   N1MAX,  L1REC   )

      USE ARGOS_ONEINT_NUC
      USE ARGOS_WRITE

      implicit logical(a-z)
 
#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_print.h"
#include "argos_error.h"

      integer    INT_AREP
      integer    N1MAX,  L1REC
 
      real*8         tol
      common /parmr/ tol
 
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
c
      real*8
     & dxij,   dyij,   dzij,   fnfct,  rr,     xij,    xijm,
     & yij,    yijm,   zij,    zijm
      integer  ibl1,   ibl2,   icxi1,  icxi2,  ij,     ijsf,
     & ic,     icons,  igu,
     & ircru,  is,     isf,    itl,    itu,    jc,     jcons,
     & jgu,    jrcru,  js,     jsf,    jtl,    jtu,    lit,
     & ljt,    nblt1,  nc2,    nc1,    nop,    ntij1,  ntij2
      logical
     & esf,    esfc,   igueq1, jgueq1
      common /one/     dxij,   dyij,   dzij,   fnfct,  rr,
     & xij,    xijm,   yij,    yijm,   zij,    zijm,   ibl1,
     & ibl2,   icxi1,  icxi2,  ij,     ijsf,   ic,     icons,
     & igu,    ircru,  is,     isf,    itl,    itu,    jc,
     & jcons,  jgu,    jrcru,  js,     jsf,    jtl,    jtu,
     & lit,    ljt,    nblt1,  nc2,    nc1,    nop,    ntij1,
     & ntij2,  esf,    esfc,   igueq1, jgueq1
c
      integer        ipq
      common /dim21/ ipq(1000) !256->1000
 
!     # /bufout/ holds some output integral file parameters.
      integer
     & ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
 
      integer    nipv
      parameter( nipv=2 )
c
c     # dummy:
      integer           ica(mcu,msu,*), icb(4,24,*), icxast(*),
     & icxst(2,*),      il(nnbft),      ilxyz(3,*),  iprst(*),
     & ipt(*),          lcr(*),         lls(*),      lmnp1(*),
     & lproju,          mcons(*),       mcrs(*),     nc(*),
     & ncon(*),         ncr(*),         nf(*),       nfct(*),
     & ngw(*),          nkcrl(6,*),     nkcru(6,*),  nklsl(4,*),
     & nklsu(4,*),      nopir(*),       npair(2,*),  nprir(2,mstu,*),
     & nrcr(*),         nt(*),          ntl(*),      ntu(*),
     & iarray(mpru,*)
      real*8
     & a(*),         ccr(*),       chg(*),        cx(*),
     & h2(*),
     & x(mcu,*),     y(mcu,*),     z(mcu,*),      eta(mrcru,mconu,*),
     & zcr(*),       zet(mconu,*), array(mpru,*)
 
c     # local:
      integer ixyz, ist, ifu, icu, isfis, jsfjs, ncc, igwu,
     & SYM_I, ICA_MIN, icc, igw, jfu, if, jf, iop,
     & ipr, itst, iprm, inx, jnx, nwt, 
     & ipr1, ipr2, iju, iiju
      logical  es
 
      integer  atebyt, forbyt
      external atebyt, forbyt
 
      real*8     zero
      parameter( zero=0d0 )
 
      integer last
      integer iapt(3)

      integer NUC_J
      integer n_bfn_sym(mstup)  ! (nbpsy or nso)
      integer N_INTEGRAL(3)     ! number of integral element of each type
      integer MN_INTEGRAL       ! maximum of N_INTEGRAL(i,i=1-3)
 
      IF (PRINT_LEVEL.GE.5) THEN 
         PRINT *,' * ECPINT : entering... ARGOS_ONEINT'
         PRINT *,'   MPRU = ',MPRU  ! array(MPRU,*)
!        PRINT *,'   NBLU = ',NBLU  ! h2(NBLU)
      ENDIF
 
!     --------------------------------------
!     # operator symmetry
!     --------------------------------------
 
      call izero_argos( nst,   nopir, 1 )
      call izero_argos( 3*nst, ilxyz, 1 )

      N_INTEGRAL     = 0
      MN_INTEGRAL    = 0
 
      if ( nnam .eq. 5 ) then
c       # x, y, and z components of spin-orbit integrals are all
c       # computed together.
        do ixyz = 1, 3
          ist                   = lxyzir(ixyz)
          nopir(ist)            = nopir(ist) + 1
          ilxyz(nopir(ist),ist) = ixyz
        enddo
        nop     = 3
        iapt(1) = 0
        iapt(2) = 0
        iapt(3) = 0
      else
        nopir(1)   = 1
        ilxyz(1,1) = 1
        nop        = 1
        iapt(1)    = 0
      endif
 
      ijsf = 0
      isf  = 0

      do 712 is = 1,ns
      ifu = nf(is)
      if(ifu.eq.0) go to 712
      icu=nc(is)
      isfis=isf
      jsf = 0
      do 708 js = 1,is
      jsfjs=jsf
      es = is.eq.js
      ic=icu
      ncc = 0


!     --------------------------------------
!     center combinations for is=js or is>js
!     --------------------------------------
      IF(es) THEN 
!        center combinations for is=js
         CALL ONEINT_NUC_ISEQJS(is,js,ic,ng,ncc,mcu,msu,
     &                          ica,icb,nfct,ngw,PRINT_LEVEL)
      ELSE
!        center combinations for is>js
         if(nf(js).eq.0) goto 708
         CALL ONEINT_NUC_ISGTJS(is,js,nc,ng,ic,mcu,msu,
     &                          ica,icb,ncc,nfct,PRINT_LEVEL)
      ENDIF

      if(ncc.gt.mccu) then
        call bummer('change mccup (one place) to ',ncc,faterr)
      endif

!     --------------------------------------
!     evaluate integrals
!     --------------------------------------
      isf=isfis
      do 704 if = 1, ifu
      isf = isf + 1
      icons=mcons(isf)
      ircru=nrcr(icons)
      lit=lmnp1(icons)
      igu=ncon(icons)
      itl=ntl(isf)
      itu=ntu(isf)
      jsf=jsfjs
      if(es) jfu=if
      do 700 jf = 1, jfu
      jsf = jsf + 1
      esf=isf.eq.jsf
      jcons=mcons(jsf)
      jrcru=nrcr(jcons)
      ljt=lmnp1(jcons)
      ijsf = ijsf + 1
      ibl1=0
      ibl2=0
      do ist=1,nst
        iop=nopir(ist)
        if(iop.ne.0) then
          if(esf) ibl1=ibl1+iop*nprir(1,ist,ijsf)
          ibl2=ibl2+iop*nprir(2,ist,ijsf)
        endif
      enddo
      iju=nt(isf)*nt(jsf)
      ntij2=iju*npair(2,ijsf)
      icxi2=icxst(2,ijsf)
      if(esf) then
        nblt1=npair(1,ijsf)
        ntij1=iju*nblt1
        nblu=ircru*ibl1+ipq(ircru)*ibl2
        if(nnam.eq.5) then
          icxi1=icxast(ijsf)
        else
          icxi1=icxst(1,ijsf)
        endif
      else
        nblu=ircru*jrcru*ibl2
      endif

      if(nblu.eq.0) go to 696

c     # allocate space.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
c     #          4:array(1:mpru,1:narray), 5:iarray(1:iamax,1:narray),
c     #          6:buffer(1:l1rec), 7:values(1:n1max),
c     #          8:labels(1:2,1:n1max), 9:hpt(1:ihwsp), 10:hwt(1:ihwsp),
c     #          11:dfac(1:ndfac), 12:binom(1:xyzdim), 13:lmf(1,l1max2),
c     #          14:lml(1:l1max2), 15:lmx(1:lmnpwr), 16:lmy(1:lmnpwr),
c     #          17:lmz(1:lmnpwr), 18:zlm(1:lmnpwr),
c     #          19:flmtx(1:3,1:lmnu2), 20:mc(1:3,1:l2m1),
c     #          21:mr(1:3,1:l2m1), 22:h2(1:nblu), 23:ijx(1:iiju),
c     #          24:ijy(1,iiju), 25:ijz(1:iiju), 26:g2(1:ic1)
      ipt(23) = ipt(22) + atebyt( nblu )
      if( nnam.le.3 ) then
        iiju = iju
      else
        iiju = 0
      endif
      ipt(24) = ipt(23) + forbyt( iiju )
      ipt(25) = ipt(24) + forbyt( iiju )
      ipt(26) = ipt(25) + forbyt( iiju )

      ipr=iprst(ijsf)
 
      if ( (ipt(26)-1) .gt. mblu ) then
        call bummer('oneint: mblu too small ',(ipt(26)-1-mblu),faterr)
      endif
 
      jgu=ncon(jcons)
      jtl=ntl(jsf)
      jtu=ntu(jsf)
      call wzero(nblu,h2,1)
      if ( (.not. es) .or. esf ) then
        do icc = 1, ncc
          ic = icb(1,1,icc)
          jc = icb(2,1,icc)
          fnfct = nfct(icc)
            call argos_stvcz( a,          ccr,        chg,        cx,
     &      eta,        a(ipt(26)), h2,         a(ipt(23)), a(ipt(24)),
     &      a(ipt(25)), ilxyz,      ipt,        lcr,        lls,
     &      a(ipt(1)),  lproju,     mcrs,       nc,         ncr,
     &      nkcrl,      nkcru,      nklsl,      nklsu,      nopir,
     &      nprir,      nt,         x,          y,          z,
     &      zcr,        zet )
        enddo
      else
        do icc = 1, ncc
          igwu = ngw(icc)
          itst = 0
          do 520 iprm = 1, 2
            if ( iprm .eq. 1 ) then
              inx=1
              jnx=2
            else
              inx=2
              jnx=1
            endif
            jc = icb(jnx,1,icc)
            if ( jc .eq. itst ) go to 520
            ic = icb(inx,1,icc)
            nwt = 1
            do igw = 2, igwu
              if( icb(jnx,igw,icc) .lt. jc ) go to 520
              nwt = nwt + 1
            enddo
            fnfct = (nwt * nfct(icc))
            call argos_stvcz( a,          ccr,        chg,        cx,
     &      eta,        a(ipt(26)), h2,         a(ipt(23)), a(ipt(24)),
     &      a(ipt(25)), ilxyz,      ipt,        lcr,        lls,
     &      a(ipt(1)),  lproju,     mcrs,       nc,         ncr,
     &      nkcrl,      nkcru,      nklsl,      nklsu,      nopir,
     &      nprir,      nt,         x,          y,          z,
     &      zcr,    zet )
            itst = jc
520       continue
        enddo
      endif
 
!     # put integrals and labels into the intermediate storage arrays.
      CALL ARGOS_WRITE_INTERMEDIATE(ircru,jrcru,nst,mstu,nprir,ijsf,
     &                              nopir,ilxyz,iapt,mpru,iarray,ipr,
     &                              N_INTEGRAL,array,h2,esf)
696   continue
700   continue
704   continue
708   continue
712   continue

!     # write 1-e integral arrays.
      MN_INTEGRAL = N_INTEGRAL(1) + N_INTEGRAL(2) + N_INTEGRAL(3)
!     print *,'stop here'
!     stop
      CALL ARGOS_WRITE_MAIN(nop,iapt,ipr,ipr1,ipr2, N1MAX,  L1REC,
     &                      il,array,iarray,n_bfn_sym,
     &                      N_INTEGRAL,MN_INTEGRAL,INT_AREP,nnam,
     &                      mpru,nbft,nnbft,nst,PRINT_LEVEL)


      if ( nnam .eq. 5 ) then
!       # all done with the antisymmetric coefficients.
!       # release the workspace, and adjust mcxu and mblu.
!       # this should probably be done in the higher-level
!       # calling program, but for now it's done down here.
        mblu = mblu + mcxu - mcxu2
        mcxu = mcxu2
      endif

      IF (PRINT_LEVEL.GE.5) THEN 
         PRINT *,' * ECPINT : exiting...  ARGOS_ONEINT'
      ENDIF

      end



