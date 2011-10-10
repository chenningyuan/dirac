      subroutine RECP_ONE_MAIN( 
     &       a,          ccr,        chg,        eta,        ica,
     &       icb,        icxast,     icxst,      ilxyz,
     &       iprst,      ipt,        lcr,        lls,        lmnp1,
     &       lproju,     mcons,      mcrs,       nc,         ncon,
     &       ncr,        nf,         nfct,       ngw,        nkcrl,
     &       nkcru,      nklsl,      nklsu,      nopir,      npair,
     &       nprir,      nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zcr,        zet,
     &       n_bfn_sym,  INT_AREP,N1MAX,L1REC,mblu,mcxu,mcxu2 )

      USE ARGOS_ONEINT_NUC
      USE RECP_IPT
      USE RECP_OUTPUT
      USE RECP_WRITE

      implicit logical(a-z)
 
#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"
#include "inc_print.h"

      integer    INT_AREP
      integer    N1MAX,  L1REC
 
      real*8         tol
      common /parmr/ tol
 
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
c
      integer
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir(3),      inam(5),     nnam, mdum(32)
      integer mblu,mcxu,mcxu2
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
!     # /bufout/ holds some output integral file parameters.
      integer
     & ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
 
      integer    nipv
      parameter( nipv=2 )
 
!     # dummy:
      integer           ica(mcu,msu,*), icb(4,24,*), icxast(*),
     & icxst(2,*),      ilxyz(3,*),  iprst(*),
     & ipt(*),          lcr(*),         lls(*),      lmnp1(*),
     & lproju,          mcons(*),       mcrs(*),     nc(*),
     & ncon(*),         ncr(*),         nf(*),       nfct(*),
     & ngw(*),          nkcrl(6,*),     nkcru(6,*),  nklsl(4,*),
     & nklsu(4,*),      nopir(*),       npair(2,*),  nprir(2,mstu,*),
     & nrcr(*),         nt(*),          ntl(*),      ntu(*)
      real*8
     & a(*),         ccr(*),       chg(*),       
     & x(mcu,*),     y(mcu,*),     z(mcu,*),      eta(mrcru,mconu,*),
     & zcr(*),       zet(mconu,*)
 
!     # local:
      integer ixyz, ist, ifu, icu, isfis, jsfjs, ncc, igwu,
     & SYM_I, ICA_MIN, icc, igw, jfu, if, jf, iop,
     & ipr, itst, iprm, inx, jnx, nwt, 
     & ipr1, ipr2, iju
      logical  es
 
      integer  atebyt, forbyt
      external atebyt, forbyt
      real*8     zero
      parameter( zero=0d0 )
      integer i,last,iapt(3),nblu,NUC_J
      integer n_bfn_sym(mstup)  ! (nbpsy or nso)
      integer N_INTEGRAL(3)     ! number of integral element of each type
      integer MN_INTEGRAL       ! maximum of N_INTEGRAL(i,i=1-3)

      CALL OUTPUT_LOC('RECP_ONE_MAIN','E') 
 
!     --------------------------------------
!     # operator symmetry
!     --------------------------------------
 
      call izero_argos( nst,   nopir, 1 )
      call izero_argos( 3*nst, ilxyz, 1 )

      N_INTEGRAL(1)  = 0
      N_INTEGRAL(2)  = 0
      N_INTEGRAL(3)  = 0
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
        OPEN (UNIT=64,FILE='RECP_ISA1') 
        OPEN (UNIT=65,FILE='RECP_ISA2')
        OPEN (UNIT=66,FILE='RECP_ISA3')
        OPEN (UNIT=67,FILE='RECP_ISB1')
        OPEN (UNIT=68,FILE='RECP_ISB2')
        OPEN (UNIT=69,FILE='RECP_ISB3')
      else
        nopir(1)   = 1
        ilxyz(1,1) = 1
        nop        = 1
        iapt(1)    = 0
        OPEN (UNIT=64,FILE='RECP_ISA1')
        OPEN (UNIT=67,FILE='RECP_ISB1')
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
        nblu=ircru*ibl1+( ircru*(ircru-1)/2 )*ibl2
        if(nnam.eq.5) then
          icxi1=icxast(ijsf)
        else
          icxi1=icxst(1,ijsf)
        endif
      else
        nblu=ircru*jrcru*ibl2
      endif

      if(nblu.eq.0) go to 696

!     # allocate space.
      CALL RECP_IPTA_ONEINT(nblu)
!     call wzero(nblu,h2,1)
      do i = 1, nblu
         IPT_ONE_H2(i) = 0.0d0
      enddo

      ipr=iprst(ijsf)
 
 
      jgu=ncon(jcons)
      jtl=ntl(jsf)
      jtu=ntu(jsf)

      if ( (.not. es) .or. esf ) then
        do icc = 1, ncc
           ic = icb(1,1,icc)
           jc = icb(2,1,icc)
           fnfct = nfct(icc)
           call RECP_STVCZ( a,ccr,chg,eta,
     &     ilxyz,ipt,lcr,lls,IPT_LMNV,lproju,mcrs,nc,ncr,
     &     nkcrl,nkcru,nklsl,nklsu,nopir,nprir,nt,x,y,z,zcr,zet )
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
            call RECP_STVCZ( a,ccr,chg,eta,ilxyz,ipt,lcr,lls,IPT_LMNV,
     &           lproju,mcrs,nc,ncr,nkcrl,nkcru,nklsl,nklsu,nopir,nprir,
     &           nt,x,y,z,zcr,zet )
            itst = jc
520       continue
        enddo
      endif
 
!     # put integrals and labels into the intermediate storage arrays.
      CALL RECP_WRITE_ISA(ircru,jrcru,nst,mstu,nprir,ijsf,
     &        nopir,ilxyz,iapt,mpru,ipr,N_INTEGRAL,IPT_ONE_H2,esf)

!     Deallocate IPT_ONE_H2
      CALL RECP_IPTD_ONEINT
696   continue
700   continue
704   continue
708   continue
712   continue

      IF (nnam.EQ.5) THEN
         CLOSE (UNIT=64)
         CLOSE (UNIT=65)
         CLOSE (UNIT=66)
         CLOSE (UNIT=67)
         CLOSE (UNIT=68)
         CLOSE (UNIT=69)
      ELSE
         CLOSE (UNIT=64)
         CLOSE (UNIT=67)
      ENDIF
!     # write 1-e integral arrays.
      MN_INTEGRAL = N_INTEGRAL(1) + N_INTEGRAL(2) + N_INTEGRAL(3)

      CALL RECP_WRITE_MAIN(nop,iapt,ipr,ipr1,ipr2, N1MAX,  L1REC,
     &        n_bfn_sym,N_INTEGRAL,MN_INTEGRAL,INT_AREP,nnam,
     &        mpru,nbft,nnbft,nst)

      if ( nnam .eq. 5 ) then
!       # all done with the antisymmetric coefficients.
!       # release the workspace, and adjust mcxu and mblu.
!       # this should probably be done in the higher-level
!       # calling program, but for now it's done down here.
!       mblu = mblu + mcxu - mcxu2
!       mcxu = mcxu2
      endif

      CALL OUTPUT_LOC('RECP_ONE_MAIN','X') 
      end



