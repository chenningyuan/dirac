*deck cints
      subroutine cints( a,              ccr,            crda,
     &                  crdb,           etai,           etaj,
     &                  g2,             g1,             gout,
     &                  ipt,            lcr,            lambu,
     &                  ltot1,          mcrs,           mproju,
     &                  nc,             ncr,            nkcrl,
     &                  nkcru,          x,              y,
     &                  z,              zcr,            zet )
c
c  lcru is the max l value + 1 for the potential.
c  ncr contains the value of n for each term.
c  zcr contains the value of alpha for each term.
c  ccr contains the coefficient of each term.
c
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
      common /parmr/ tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      common /dim21/ ipq(1000) !256->1000
      dimension a(*),ccr(*),crda(lit,3),crdb(ljt,3),etai(mrcru,*),
     1  etaj(mrcru,*),g2(*),g1(*),gout(*),ipt(*),lcr(*),
     2  mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*),x(mcu,*),
     3  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
 
      do i=1,3
        crda(1,i)=a1
        crdb(1,i)=a1
      enddo
      if(esfc) iu=ipq(lit+1)
      if(igueq1) then
        fctr1=a4*etai(1,1)
      else
        fctr1=a4
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
c
c  i primitive
c
      do 490 ig=1,igu
      if(.not.igueq1) call wzero(nc1,g1,1)
      ai=zet(ig,icons)
      tai=ai+ai
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) call wzero(ij,gout,1)
      aj=zet(jg,jcons)
      taj=aj+aj
      aa=ai+aj
      taa=aa+aa
      aaa=(ai-aj)/aa
      apr=ai*aj/aa
      aarr1=apr*rr
c
c  form core potential integrals
c
      xp=xij+aaa*xijm
      yp=yij+aaa*yijm
      zp=zij+aaa*zijm
      fctr2=fctr1
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      do 350 ks=1,ns
        kcrs=mcrs(ks)
        if(kcrs.eq.0) go to 350
        lcru=lcr(kcrs)
        if(lcru.lt.0) go to 340
        do 330 kc=1,nc(ks)
          xc=x(kc,ks)
          yc=y(kc,ks)
          zc=z(kc,ks)
          xka=xc-xi
          yka=yc-yi
          zka=zc-zi
          ca=sqrt(xka*xka+yka*yka+zka*zka)
          if(lit.ge.2) then
            crda(2,1)=xka
            crda(2,2)=yka
            crda(2,3)=zka
            if(lit.ge.3) then
              do i=1,3
                do j=3,lit
                  crda(j,i)=crda(2,i)*crda(j-1,i)
                enddo
              enddo
            endif
          endif
          xkb=xc-xj
          ykb=yc-yj
          zkb=zc-zj
          cb=sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
          if(ljt.ge.2) then
            crdb(2,1)=xkb
            crdb(2,2)=ykb
            crdb(2,3)=zkb
            if(ljt.ge.3) then
              do 240 i=1,3
                do 230 j=3,ljt
                  crdb(j,i)=crdb(2,i)*crdb(j-1,i)
  230           continue
  240         continue
            endif
          endif

!         type1 pseudo
          if(aarr1.le.tol) then
            xk=xp-xc
            yk=yp-yc
            zk=zp-zc
            call pseud1(a,a(ipt(31)),ccr,crda,crdb,gout,ipt,a(ipt(1)),
     &        ltot1,ncr,nkcrl,nkcru,a(ipt(32)),a(ipt(33)),a(ipt(34)),
     &        a(ipt(35)),a(ipt(36)),zcr)
          endif

!         type2 pseudo
          aarr2=apr*(ca-cb)**2
          if(aarr2.le.tol.and.lcru.ne.0) then
!           call pseud2( ccr, gout, lmnv, ncr, nkcrl, nkcru, zcr )
            call pseud2( a, a(ipt(37)), a(ipt(38)), ccr, crda, crdb,
     &        gout, ipt, lambu, ltot1, mproju, ncr, nkcrl,
     &        nkcru, a(ipt(39)), zcr )
          endif
  330   continue
  340   continue
  350 continue
!
!  j transformation
!
      if(jgueq1) go to 400
 
      j1=0
      do 380 jrcr=1,jrcru
        do 370 i=1,ij
          g1(j1+i)=g1(j1+i)+gout(i)*etaj(jrcr,jg)
  370   continue
        j1=j1+ij
  380 continue
  390 continue
!
!  i transformation
!
  400 if(igueq1) return
 
      if(esfc) i1=0
      ij2=0
      do 470 ircr=1,ircru
        j1=0
        if(esfc) then
          do 440 jrcr=1,ircr
            do 430 i=1,iu
              i1n=i1+i
              do 420 j=1,iu
                g2(ij2+j)=g2(ij2+j)+g1(j1+j)*etai(ircr,ig)+
     1                              g1(i1n)*etai(jrcr,ig)
                i1n=i1n+iu
  420         continue
              ij2=ij2+iu
              j1=j1+iu
  430       continue
  440     continue
          i1=i1+ij
        else
          do 460 jrcr=1,jrcru
            do 450 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      return
      end

      subroutine pseud1(a,ang,ccr,crda,crdb,gout,ipt,lmnv,ltot1,ncr,
     &  nkcrl,nkcru,q,qsum,xab,yab,zab,zcr)
c
c  compute type 1 core potential integrals
c
      USE RECP_FUNCTION
      USE RECP_INT_C2

      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0)
      common /parmr/ tol
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),ang(ltot1,*),ccr(*),crda(lit,3),crdb(ljt,3),
     &  gout(*),ipt(*),lmnv(3,*),ncr(*),nkcrl(6,*),nkcru(6,*),
     &  q(ltot1,*),qsum(ltot1,*),xab(*),yab(*),zab(*),zcr(*)
c     # should include kcr loop in recur1 and use recurrence relations
c     # on qsum
 
      rp2=xk*xk+yk*yk+zk*zk

      if(rp2.eq.0.0d0) then
        rp=0.0d0
        arp2=0.0d0
        alpt=0.0d0
        rk=0.0d0
        lamu=1
      else
        rp=sqrt(rp2)
        xk=xk/rp
        yk=yk/rp
        zk=zk/rp
        arp2=aa*rp2
        alpt=aa*arp2
        rk=taa*rp
        lamu=ltot1
      endif
!
!  compute radial integrals and sum over potential terms
!
      call wzero((ltot1*lamu),qsum,1)
      kcrl=nkcrl(1,kcrs)
      kcru=nkcru(1,kcrs)
c     do 40 kcr=kcrl,kcru
c       alpha=aa+zcr(kcr)
c       # exponential factor from q functions included in dum
c       dum=aarr1+zcr(kcr)*arp2/alpha
c       if(dum.gt.tol) go to 40
c       prd=fctr2*ccr(kcr)*exp(-dum)
c       t=alpt/alpha
c       call recur1(alpha,a(ipt(11)),ncr(kcr),ltot1,q,t,rk)
c       do 30 lam=1,lamu
c         nhi=ltot1-mod(ltot1-lam,2)
c         do 20 n=lam,nhi,2
c           qsum(n,lam)=qsum(n,lam)+prd*q(n,lam)
c  20     continue
c  30   continue
c  40 continue

      call RECP_INT_C2_RAD1
     &    (aa,aarr1,alpt,arp2,ccr,a(ipt(11)),fctr2,kcrl,kcru,
     &     lamu,ltot1,ncr,qsum,rk,tol,zcr)

      ijt=0
      do 90 it=itl,itu
        na1=lmnv(1,it)+1
        la1=lmnv(2,it)+1
        ma1=lmnv(3,it)+1
        do 80 jt=jtl,jtu
          ijt=ijt+1
          s=a0
          nb1=lmnv(1,jt)+1
          lb1=lmnv(2,jt)+1
          mb1=lmnv(3,jt)+1
!
!  compute angular integrals
!
          call RECP_FUNCTION_FACAB
     &         (a(ipt(12)),na1,nb1,crda(1,1),crdb(1,1),xab)
          call RECP_FUNCTION_FACAB
     &         (a(ipt(12)),la1,lb1,crda(1,2),crdb(1,2),yab)
          call RECP_FUNCTION_FACAB
     &         (a(ipt(12)),ma1,mb1,crda(1,3),crdb(1,3),zab)
          call ang1(ang,a(ipt(11)),na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,
     &      a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &      ltot1,xab,yab,zab,xk,yk,zk,a(ipt(18)))
!         call RECP_INT_C2_ANG1 
!    &        (ang,a(ipt(11)),na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,
!    &         a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
!    &         ltot1,xab,yab,zab,xk,yk,zk,a(ipt(18)))

!
!  combine angular and radial integrals
!
          do 70 lam=1,lamu
            nhi=ltot1-mod(ltot1-lam,2)
            do 60 n=lam,nhi,2
              s=s+ang(n,lam)*qsum(n,lam)
   60       continue
   70     continue
          gout(ijt)=gout(ijt)+s
   80   continue
   90 continue
      return
      end

*deck pseud2
      subroutine pseud2( a, anga, angb, ccr, crda, crdb, gout, ipt,
     &  lambu, ltot1, mproju, ncr, nkcrl, nkcru, qsum, zcr )
c
c  compute type 2 core potential integrals
c
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0)
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),anga(lit,mproju,*),angb(ljt,mproju,*),ccr(*),
     &  crda(lit,3),crdb(ljt,3),gout(*),ipt(*),ncr(*),
     &  nkcrl(6,*),nkcru(6,*),qsum(ltot1,lambu,*),zcr(*)
c
      if(ca.eq.a0) then
        rka=a0
        lmau=1
      else
        xka=-xka/ca
        yka=-yka/ca
        zka=-zka/ca
        rka=tai*ca
        lmau=lcru+(lit-1)
      endif
      if(cb.eq.a0) then
        rkb=a0
        lmbu=1
      else
        xkb=-xkb/cb
        ykb=-ykb/cb
        zkb=-zkb/cb
        rkb=taj*cb
        lmbu=lcru+(ljt-1)
      endif
      if((ca.eq.a0).and.(cb.eq.a0)) then
        lhi=min(lcru,lit,ljt)
        llo=mod((lit-1),2)+1
        if(llo.ne.mod((ljt-1),2)+1.or.llo.gt.lhi) return
        inc=2
      elseif(ca.eq.a0) then
        lhi=min(lcru,lit)
        llo=mod((lit-1),2)+1
        if(llo.gt.lhi) return
        inc=2
      elseif(cb.eq.a0) then
        lhi=min(lcru,ljt)
        llo=mod((ljt-1),2)+1
        if(llo.gt.lhi) return
        inc=2
      else
        lhi=lcru
        llo=1
        inc=1
      endif
      do 88 l=llo,lhi,inc
      mhi=l+l-1
      lmalo=max(l-(lit-1),1)
      lmahi=min(lmau,l+(lit-1))
      lmblo=max(l-(ljt-1),1)
      lmbhi=min(lmbu,l+(ljt-1))
c
c  compute radial integrals
c
      kcrl=nkcrl(l+1,kcrs)
      kcru=nkcru(l+1,kcrs)
            call rad2( a,          ccr,        ipt,        kcrl,
     &                 kcru,       l,          lambu,      lmahi,
     &                 lmalo,      lmbhi,      lmblo,      ltot1,
     &                 ncr,        qsum,       rka,        rkb,
     &                 zcr )
c
c  compute angular integrals and combine with radial integrals
c
      ijt=0
      do 84 it=itl,itu
      call ang2(anga,a(ipt(12)),crda,a(ipt(11)),it,l,lit,lmalo,lmahi,
     &  a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &  a(ipt(1)),mproju,xka,yka,zka,a(ipt(18)))
      do 80 jt=jtl,jtu
      ijt=ijt+1
      s=a0
      call ang2(angb,a(ipt(12)),crdb,a(ipt(11)),jt,l,ljt,lmblo,lmbhi,
     &  a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &  a(ipt(1)),mproju,xkb,ykb,zkb,a(ipt(18)))
      do 76 lama=lmalo,lmahi
        ldifa1=abs(l-lama)+1
        nlmau=lit-mod(lit-ldifa1,2)
        do 72 lamb=lmblo,lmbhi
          ldifb=abs(l-lamb)
          nlmbu=(ljt-1)-mod((ljt-1)-ldifb,2)
          nlo=ldifa1+ldifb
          nhi=nlmau+nlmbu
          do 68 n=nlo,nhi,2
            nlmalo=max(ldifa1,n-nlmbu)
            nlmahi=min(nlmau,n-ldifb)
            angp=a0
            do 60 m=1,mhi
              do 56 nlma=nlmalo,nlmahi,2
                angp=angp+anga(nlma,m,lama)*angb((n+1)-nlma,m,lamb)
   56         continue
   60       continue
            s=s+angp*qsum(n,lamb,lama)
   68     continue
   72   continue
   76 continue
      gout(ijt)=gout(ijt)+s
   80 continue
   84 continue
   88 continue
      return
      end

      subroutine lsints( a,              ccr,            crda,
     &                   crdb,           etai,           etaj,
     &                   g2,             g1,             gout,
     &                   ipt,            lls,            lambu,
     &                   ltot1,          mcrs,           mproju,
     &                   nc,             ncr,            nklsl,
     &                   nklsu,          x,              y,
     &                   z,              zcr,            zet )
c
c  lcru is the max l value + 1 for the potential.
c  ncr contains the value of n for each term.
c  zcr contains the value of alpha for each term.
c  ccr contains the coefficient of each term.
c
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
      common /parmr/ tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      common /dim21/ ipq(1000) !256->1000
      dimension a(*),ccr(*),crda(lit,3),crdb(ljt,3),etai(mrcru,*),
     1  etaj(mrcru,*),g2(*),g1(*),gout(*),ipt(*),lls(*),
     2  mcrs(*),nc(*),ncr(*),nklsl(4,*),nklsu(4,*),x(mcu,*),
     3  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
c
      do 120 i=1,3
        crda(1,i)=a1
  120 crdb(1,i)=a1
      if(esfc) iu=ipq(lit+1)
      if(igueq1) then
        fctr1=a4*etai(1,1)
      else
        fctr1=a4
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
c
c  i primitive
c
      do 490 ig=1,igu
      if(.not.igueq1) call wzero(nc1,g1,1)
      ai=zet(ig,icons)
      tai=ai+ai
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) call wzero(ij,gout,1)
      aj=zet(jg,jcons)
      taj=aj+aj
      aa=ai+aj
      taa=aa+aa
      apr=ai*aj/aa
c
c  form spin-orbit potential integrals
c
      fctr2=fctr1
      if((esfc.and..not.igueq1).and.ig.eq.jg) fctr2=a1s2*fctr2
      do 350 ks=1,ns
      kcrs=mcrs(ks)
      if(kcrs.eq.0) go to 350
      lcru=lls(kcrs)+1
      if(lcru.le.1) go to 340
      do 330 kc=1,nc(ks)
      xc=x(kc,ks)
      yc=y(kc,ks)
      zc=z(kc,ks)
      xka=xc-xi
      yka=yc-yi
      zka=zc-zi
      ca=sqrt(xka*xka+yka*yka+zka*zka)
      if(lit.ge.2) then
        crda(2,1)=xka
        crda(2,2)=yka
        crda(2,3)=zka
        if(lit.ge.3) then
          do i=1,3
            do j=3,lit
              crda(j,i)=crda(2,i)*crda(j-1,i)
            enddo
          enddo
        endif
      endif
      xkb=xc-xj
      ykb=yc-yj
      zkb=zc-zj
      cb=sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
      if(ljt.ge.2) then
        crdb(2,1)=xkb
        crdb(2,2)=ykb
        crdb(2,3)=zkb
        if(ljt.ge.3) then
          do 240 i=1,3
            do 230 j=3,ljt
              crdb(j,i)=crdb(2,i)*crdb(j-1,i)
  230       continue
  240     continue
        endif
      endif
      aarr2=apr*(ca-cb)**2
      if(aarr2.le.tol) then
        call pseud3( a, a(ipt(37)), a(ipt(38)), ccr, crda, crdb,
     &    a(ipt(19)), gout, ipt, lambu, ltot1, a(ipt(20)), a(ipt(21)),
     &    mproju, ncr, nklsl, nklsu, a(ipt(39)), zcr )
      endif
  330 continue
  340 continue
  350 continue
c
c  j transformation
c
      if(jgueq1) go to 400
c
      j1=0
      do 380 jrcr=1,jrcru
        do 370 i=1,ij
          g1(j1+i)=g1(j1+i)+gout(i)*etaj(jrcr,jg)
  370   continue
        j1=j1+ij
  380 continue
  390 continue
  400 continue
c
c  i transformation
c
c 400 if(igueq1) return
      if(igueq1) return
c
      if(esfc) i1=0
      ij2=0
      do 470 ircr=1,ircru
        j1=0
        if(esfc) then
          do 440 jrcr=1,ircr
            i1m=i1
            do 430 i=1,iu
              i1n=i1m
              do 420 j=1,iu
                do 410 k=1,3
                  g2(ij2+k)=g2(ij2+k)+g1(j1+k)*etai(ircr,ig)-
     &                                g1(i1n+k)*etai(jrcr,ig)
  410           continue
                ij2=ij2+3
                j1=j1+3
                i1n=i1n+(3*iu)
  420         continue
              i1m=i1m+3
  430       continue
  440     continue
          i1=i1+ij
        else
          do 460 jrcr=1,jrcru
            do 450 i=1,ij
              g2(ij2+i)=g2(ij2+i)+g1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      return
      end
*deck pseud3
      subroutine pseud3( a, anga, angb, ccr, crda, crdb, flmtx, gout,
     &  ipt, lambu, ltot1, mc, mr, mproju, ncr, nklsl, nklsu, qsum,
     &  zcr )
c
c  computes spin-orbit potential integrals
c
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0)
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),anga(lit,mproju,*),angb(ljt,mproju,*),ccr(*),
     &  crda(lit,3),crdb(ljt,3),flmtx(3,*),gout(*),ipt(*),mc(3,*),
     &  mr(3,*),ncr(*),nklsl(4,*),nklsu(4,*),qsum(ltot1,lambu,*),
     &  zcr(*)
c
      if(ca.eq.a0) then
        rka=a0
        lmau=1
      else
        xka=-xka/ca
        yka=-yka/ca
        zka=-zka/ca
        rka=tai*ca
        lmau=lcru+(lit-1)
      endif
      if(cb.eq.a0) then
        rkb=a0
        lmbu=1
      else
        xkb=-xkb/cb
        ykb=-ykb/cb
        zkb=-zkb/cb
        rkb=taj*cb
        lmbu=lcru+(ljt-1)
      endif
      if((ca.eq.a0).and.(cb.eq.a0)) then
        lhi=min(lcru,lit,ljt)
        llo=mod(lit,2)+2
        if(llo.ne.(mod(ljt,2)+2).or.llo.gt.lhi) return
        inc=2
      elseif(ca.eq.a0) then
        lhi=min(lcru,lit)
        llo=mod(lit,2)+2
        if(llo.gt.lhi) return
        inc=2
      elseif(cb.eq.a0) then
        lhi=min(lcru,ljt)
        llo=mod(ljt,2)+2
        if(llo.gt.lhi) return
        inc=2
      else
        lhi=lcru
        llo=2
        inc=1
      endif
      do 88 l=llo,lhi,inc
      mhi=l+l-3
      lmalo=max(l-(lit-1),1)
      lmahi=min(lmau,l+(lit-1))
      lmblo=max(l-(ljt-1),1)
      lmbhi=min(lmbu,l+(ljt-1))
c
c  compute radial integrals
c
      kcrl=nklsl(l-1,kcrs)
      kcru=nklsu(l-1,kcrs)
            call rad2( a,          ccr,        ipt,        kcrl,
     &                 kcru,       l,          lambu,      lmahi,
     &                 lmalo,      lmbhi,      lmblo,      ltot1,
     &                 ncr,        qsum,       rka,        rkb,
     &                 zcr )
c
c  compute angular integrals and combine with radial integrals
c
      ijt=0
      do 84 it=itl,itu
      call ang2(anga,a(ipt(12)),crda,a(ipt(11)),it,l,lit,lmalo,lmahi,
     &  a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &  a(ipt(1)),mproju,xka,yka,zka,a(ipt(18)))
      do 80 jt=jtl,jtu
      call ang2(angb,a(ipt(12)),crdb,a(ipt(11)),jt,l,ljt,lmblo,lmbhi,
     &  a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &  a(ipt(1)),mproju,xkb,ykb,zkb,a(ipt(18)))
      do 76 lama=lmalo,lmahi
      ldifa1=abs(l-lama)+1
      nlmau=lit-mod(lit-ldifa1,2)
      do 72 lamb=lmblo,lmbhi
      ldifb=abs(l-lamb)
      nlmbu=(ljt-1)-mod((ljt-1)-ldifb,2)
      nlo=ldifa1+ldifb
      nhi=nlmau+nlmbu
      do 68 n=nlo,nhi,2
      nlmalo=max(ldifa1,n-nlmbu)
      nlmahi=min(nlmau,n-ldifb)
      do 64 kt=1,3
        s=a0
        do 60 m=1,mhi
          angp=a0
          do 56 nlma=nlmalo,nlmahi,2
            angp=angp+
     1        anga(nlma,mr(kt,m),lama)*angb((n+1)-nlma,mc(kt,m),lamb)-
     2        anga(nlma,mc(kt,m),lama)*angb((n+1)-nlma,mr(kt,m),lamb)
   56     continue
          s=s+angp*flmtx(kt,(l-2)**2+m)
   60   continue
        gout(ijt+kt)=gout(ijt+kt)+s*qsum(n,lamb,lama)
   64 continue
   68 continue
   72 continue
   76 continue
      ijt=ijt+3
   80 continue
   84 continue
   88 continue
      return
      end
*deck rad2
      subroutine rad2( a,          ccr,        ipt,        kcrl,
     &                 kcru,       l,          lambu,      lmahi,
     &                 lmalo,      lmbhi,      lmblo,      ltot1,
     &                 ncr,        qsum,       rka,        rkb,
     &                 zcr )
c
c  compute type 2 radial integrals.
c
c  28-nov-90 new version replaced old version. -rmp
c  19-jan-97 put bessel formula into a separate subroutine, qbess. -rmp
c
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0, eps1=1.0d-15, a1=1.0d0, a2=2.0d0,
     1  a4=4.0d0, a50=50.0d0)
      common /parmr/ tol
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*), ccr(*), ipt(*), ncr(*), qsum(ltot1,lambu,*),
     &  zcr(*)
c
      call wzero(ltot1*lambu*lmahi,qsum,1)
c     # sum over potential terms
      do 68 kcr=kcrl,kcru
      npi=ncr(kcr)
      alpha=aa+zcr(kcr)
      rc=(rka+rkb)/(alpha+alpha)
      arc2=alpha*rc*rc
      dum=aarr2+zcr(kcr)*arc2/aa
      if(dum.gt.tol) go to 68
      prd=fctr2*ccr(kcr)*exp(-dum)
c
      if(rka.eq.a0.and.rkb.eq.a0) then
c       # rka=0 and rkb=0
        rk=a0
        t=a0
        qsum(ltot1,1,1)=qsum(ltot1,1,1)+prd*
     &             qcomp(alpha,a(ipt(11)),npi+ltot1-1,0,t,rk)
      elseif(rka.eq.a0) then
c       # rka=0 and rkb>0
        rk=rkb
        t=arc2
        do 12 lamb=l,lmbhi
          qsum(lamb-l+lit,lamb,1)=qsum(lamb-l+lit,lamb,1)+prd*
     1             qcomp(alpha,a(ipt(11)),npi+lamb-l+lit-1,lamb-1,t,rk)
   12   continue
      elseif(rkb.eq.a0) then
c       # rka>0 and rkb=0
        rk=rka
        t=arc2
        do 16 lama=l,lmahi
          qsum(lama-l+ljt,1,lama)=qsum(lama-l+ljt,1,lama)+prd*
     1             qcomp(alpha,a(ipt(11)),npi+lama-l+ljt-1,lama-1,t,rk)
   16   continue
      elseif(npi.eq.2) then
c       # rka>0 and rkb>0; use bessel function formula.
c       # to be applicable for a set of integrals, must have nu.le.l and
c       # nu.eq.(integer), where nu=l+1-npi/2, so it is used here only
c       # for the npi=2 case.  It can't be used at all for npi=(odd) and
c       # only for partial sets for npi=0
        nu=l
            call qbess( alpha,      a(ipt(40)), a(ipt(41)), a(ipt(42)),
     &      a(ipt(12)), a(ipt(43)), a(ipt(44)), a(ipt(45)), a(ipt(11)),
     &      l,          lambu,      lmahi,      lmbhi,      ltot1,
     &      nu,         prd,        qsum,       rka,        rkb,
     &      a(ipt(46)) )
      elseif(arc2.ge.a50) then
c       # rka>0 and rkb>0; use pts and wts method
c       # estimate radial integrals and compare to threshold
        qlim=abs(prd)/(max(a1,(rc+rc)*rka)*max(a1,(rc+rc)*rkb))*
     1    sqrt(a4*(tai+tai)**lit*(taj+taj)**ljt*sqrt(tai*taj)/alpha)
        if(rc.lt.ca) then
          nlim=npi
          qlim=qlim*ca**(lit-1)
        else
          nlim=npi+(lit-1)
        endif
        if(rc.lt.cb) then
          qlim=qlim*cb**(ljt-1)
        else
          nlim=nlim+(ljt-1)
        endif
        if(qlim*rc**nlim.ge.eps1) then
          call ptwt(a(ipt(47)),arc2,a(ipt(48)),a(ipt(11)),npi,l,lambu,
     &      ltot1,lmahi,lmbhi,alpha,a(ipt(49)),a(ipt(50)),rc,rka,rkb,
     &      prd,a(ipt(9)),a(ipt(10)),qsum)
        endif
      else
c       # rka>0 and rkb>0; use partially asymptotic method
        call qpasy(alpha,a(ipt(11)),npi,l,lambu,lmahi,lmbhi,ltot1,rka,
     &    rkb,fctr2*ccr(kcr),dum+arc2,qsum)
      endif
   68 continue
c
      if(rka.eq.a0.and.rkb.ne.a0) then
c       # rka=0 and rkb>0
        f2lmb3=(2*lmbhi-3)
        do 76 lamb=lmbhi-2,lmblo,-1
          nlo=abs(lamb-l+1)+lit+1
          nhi=ljt-mod((ljt-1)-abs(lamb-l),2)+lit-1
          do 72 n=nlo,nhi,2
            qsum(n,lamb,1)=qsum(n,lamb+2,1)+
     &                     (f2lmb3/rkb)*qsum(n-1,lamb+1,1)
   72     continue
          f2lmb3=f2lmb3-a2
   76   continue
      elseif(rka.ne.a0.and.rkb.eq.a0) then
c       # rka>0 and rkb=0
        f2lma3=(2*lmahi-3)
        do 84 lama=lmahi-2,lmalo,-1
          nlo=abs(lama-l+1)+ljt+1
          nhi=lit-mod((lit-1)-abs(lama-l),2)+ljt-1
          do 80 n=nlo,nhi,2
            qsum(n,1,lama)=qsum(n,1,lama+2)+
     &                     (f2lma3/rka)*qsum(n-1,1,lama+1)
   80     continue
          f2lma3=f2lma3-a2
   84   continue
      elseif(rka.ne.a0.and.rkb.ne.a0) then
c       # rka>0 and rkb>0
        f2lma3=(lmahi+lmahi+1)
        do 96 lama=lmahi,lmalo,-1
          ldifa1=abs(l-lama)+1
          f2lmb3=(2*lmbhi+1)
          do 92 lamb=lmbhi,lmblo,-1
            ldifb=abs(l-lamb)
            nlo=ldifa1+ldifb
            nhi=(ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
            do 88 n=nlo,nhi,2
              if(n-(lama+lamb).eq.(1-l-l)) go to 88
              if(lama.gt.(lmahi-2).or.n.le.(abs(l-lama-2)+ldifb)) then
c               # lamb recursion
                qsum(n,lamb,lama)=qsum(n,lamb+2,lama)+
     1                            (f2lmb3/rkb)*qsum(n-1,lamb+1,lama)
              else
c               # lama recursion
                qsum(n,lamb,lama)=qsum(n,lamb,lama+2)+
     1                            (f2lma3/rka)*qsum(n-1,lamb,lama+1)
              endif
   88       continue
            f2lmb3=f2lmb3-a2
   92     continue
          f2lma3=f2lma3-a2
   96   continue
      endif
      return
      end
*deck qbess
      subroutine qbess( alpha,      apwr,       aterm1,     aterm2,
     &      binom,      bpref,      bpwr,       bterm1,     dfac,
     &      l,          lambu,      lmahi,      lmbhi,      ltot1,
     &      nu,         prd,        qsum,       rka,        rkb,
     &      ssi )
c
c  compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
c  using the bessel function formula for
c  lama=max(l,nu) to lmahi, lamb=max(l,nu) to lmbhi, n=lama+lamb-l-l
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1=1.0d0, a2=2.0d0, a4=4.0d0)
      dimension apwr(*), aterm1(*), aterm2(*), binom(*), bpref(*),
     &  bpwr(*), bterm1(*), dfac(*), qsum(ltot1,lambu,*), ssi(*)
c
c     nu=l+1-npi/2
c     # bessel function formula applies to all npi=2 cases, no npi=1
c     # cases, and some npi=0 cases.
      num1=nu-1
      lmlo=max(l,nu)
      lmaphi=lmahi-num1
      lmbphi=lmbhi-num1
      fcta=rka/(alpha+alpha)
      fct=rka*fcta
      fctb=rkb/(alpha+alpha)
      bpref(1)=fctb**num1
      do 24 lambp=2,lmbphi
        bpref(lambp)=fctb*bpref(lambp-1)
   24 continue
      apwr(1)=a1
      apwr(2)=fct
      do 28 lambp=3,lmbphi
        apwr(lambp)=fct*apwr(lambp-1)
   28 continue
      fct=rkb*fctb
      bpwr(1)=a1
      bpwr(2)=fct
      do 32 lamap=3,lmaphi
        bpwr(lamap)=fct*bpwr(lamap-1)
   32 continue
      lmihi=lmaphi+lmbphi+(nu-2)
      call ssibfn(lmihi-1,rka*fctb,ssi)
      do 36 lami=1,lmihi
        ssi(lami)=ssi(lami)/dfac(lami+lami-1)
   36 continue
      lmplo=lmlo-num1
      fctra=(alpha+alpha)**(nu-2)*fcta**(lmlo-1)*prd/sqrt(a4*alpha)*
     &  ((dfac(2*(2*lmplo+num1)-1)/dfac(2*(lmplo+num1)+1))*
     &  dfac(2*num1+1))/dfac(2*(lmplo+num1)+1)
      fctran=(2*(2*lmplo+num1)-1)
      fctrad=(2*(lmplo+num1)+1)
      do 64 lamap=lmplo,lmaphi
        fctru=a1
        fctrun=(nu+num1)
        fctrud=(2*(lamap+num1)+1)
        do 40 iu=1,lmbphi
          aterm1(iu)=fctru*apwr(iu)
          fctru=fctru*fctrun/fctrud
          fctrun=fctrun+a2
          fctrud=fctrud+a2
   40   continue
        do 44 it=1,lamap
          bterm1(it)=binom((lamap*(lamap-1))/2+it)*bpwr(it)
   44   continue
        fctrb=fctra
        fctrbn=(2*(lamap+lmplo+num1)-1)
        fctrbd=(2*(lmplo+num1)+1)
        do 60 lambp=lmplo,lmbphi
          n=((2*(nu-l)-1)+lamap)+lambp
          do 48 iu=1,lambp
            aterm2(iu)=binom((lambp*(lambp-1))/2+iu)*aterm1(iu)
   48     continue
          sum=a0
          fctrt=a1
          fctrtn=(nu+num1)
          fctrtd=(2*(lambp+num1)+1)
          do 56 it=1,lamap
            do 52 iu=1,lambp
              sum=sum+aterm2(iu)*(fctrt*bterm1(it))*
     &          ssi((it+(nu-2))+iu)
   52       continue
            fctrt=fctrt*fctrtn/fctrtd
            fctrtn=fctrtn+a2
            fctrtd=fctrtd+a2
   56     continue
          qsum(n,lambp+num1,lamap+num1)=qsum(n,lambp+num1,lamap+num1)+
     &      fctrb * bpref(lambp) * sum
          fctrb=fctrb*fctrbn/fctrbd
          fctrbn=fctrbn+a2
          fctrbd=fctrbd+a2
   60   continue
        fctra=fcta*fctra*fctran/fctrad
        fctran=fctran+a2
        fctrad=fctrad+a2
   64 continue
      return
      end
