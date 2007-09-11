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
      common /parmr/ cutoff, tol
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
      common /dim21/ ipq(256)
      dimension a(*),ccr(*),crda(lit,3),crdb(ljt,3),etai(mrcru,*),
     1  etaj(mrcru,*),g2(*),g1(*),gout(*),ipt(*),lcr(*),
     2  mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*),x(mcu,*),
     3  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
c
      do 120 i=1,3
        crda(1,i)=a1
        crdb(1,i)=a1
  120 continue
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
              do 220 i=1,3
                do 210 j=3,lit
                  crda(j,i)=crda(2,i)*crda(j-1,i)
  210           continue
  220         continue
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
          if(aarr1.le.tol) then
            xk=xp-xc
            yk=yp-yc
            zk=zp-zc
            call pseud1(a,a(ipt(31)),ccr,crda,crdb,gout,ipt,a(ipt(1)),
     &        ltot1,ncr,nkcrl,nkcru,a(ipt(32)),a(ipt(33)),a(ipt(34)),
     &        a(ipt(35)),a(ipt(36)),zcr)
          endif
          aarr2=apr*(ca-cb)**2
          if(aarr2.le.tol.and.lcru.ne.0) then
c           call pseud2( ccr, gout, lmnv, ncr, nkcrl, nkcru, zcr )
            call pseud2( a, a(ipt(37)), a(ipt(38)), ccr, crda, crdb,
     &        gout, ipt, lambu, ltot1, mproju, ncr, nkcrl,
     &        nkcru, a(ipt(39)), zcr )
          endif
  330   continue
  340   continue
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
c
c  i transformation
c
  400 if(igueq1) return
c
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
c
c  RCS $Revision$  $Date$
c
*deck pseud1
      subroutine pseud1(a,ang,ccr,crda,crdb,gout,ipt,lmnv,ltot1,ncr,
     &  nkcrl,nkcru,q,qsum,xab,yab,zab,zcr)
c
c  compute type 1 core potential integrals
c
c  03-jun-92 fujitsu compiler directive added to workaround
c            a compiler optimization bug (fix by Roger Edberg). -rls
c
      implicit real*8 (a-h,o-z)
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0)
      common /parmr/ cutoff, tol
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
c
      rp2=xk*xk+yk*yk+zk*zk
      if(rp2.eq.a0) then
        rp=a0
        arp2=a0
        alpt=a0
        rk=a0
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
c
c  compute radial integrals and sum over potential terms
c
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
      call rad1(aa,aarr1,alpt,arp2,ccr,a(ipt(11)),fctr2,kcrl,kcru,
     &  lamu,ltot1,ncr,qsum,rk,tol,zcr)
c
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
c
c  compute angular integrals
c
          call facab(a(ipt(12)),na1,nb1,crda(1,1),crdb(1,1),xab)
          call facab(a(ipt(12)),la1,lb1,crda(1,2),crdb(1,2),yab)
          call facab(a(ipt(12)),ma1,mb1,crda(1,3),crdb(1,3),zab)
          call ang1(ang,a(ipt(11)),na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,
     &      a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &      ltot1,xab,yab,zab,xk,yk,zk,a(ipt(18)))
c
c  combine angular and radial integrals
c
          do 70 lam=1,lamu
            nhi=ltot1-mod(ltot1-lam,2)
*mdc*if fujitsu
*c           # need vocl here to force scalar execution - compiler bug??
*c           # 03-jun-92 Roger Edberg
**vocl loop,scalar
*mdc*endif
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
*deck lsints
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
      common /parmr/ cutoff, tol
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
      common /dim21/ ipq(256)
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
          do 220 i=1,3
            do 210 j=3,lit
              crda(j,i)=crda(2,i)*crda(j-1,i)
  210       continue
  220     continue
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
      common /parmr/ cutoff, tol
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
*deck twoint_argos
      subroutine twoint_argos( a,          cx,         eta,        h4,
     &      ica,         icb,        icxst,      il,         iprst,
     &      ipt,         lmnp1,      lmnv,       mcons,      nc,
     &      ncon,        nf,         nfct,       ngw,        npair,
     &      nprir,       nrcr,       nt,         ntl,        ntu,
     &      x,           y,          z,          zet,        buffer,
     &      values,      labels,     ibitv )
c
c  11-sep-91 final parallel version. -rls/rjh
c  01-dec-90 SIFS version. -rls
c
      implicit logical(a-z)
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      logical       eprsf,    esfi,    esfk
      common /lgcl/ eprsf(3), esfi(3), esfk(3)
c
      integer
     & itl,    itu,    jtl,    jtu,    ktl,    ktu,    ltl,    ltu,
     & jcu,    kcu,    lcu,    inx,    jnx,    knx,    lnx,    nwt,
     & nwt1,           ijsf,           klsf,           icxs,
     & kcxs,           npri,           nprk,           iesfb,  kesfb,
     & ircru,  jrcru,  krcru,  lrcru
      common /ntgr/
     & itl,    itu,    jtl,    jtu,    ktl,    ktu,    ltl,    ltu,
     & jcu,    kcu,    lcu,    inx,    jnx,    knx,    lnx,    nwt,
     & nwt1(3),        ijsf(3),        klsf(3),        icxs(2,3),
     & kcxs(2,3),      npri(2,3),      nprk(2,3),      iesfb,  kesfb,
     & ircru,  jrcru,  krcru,  lrcru
c
      integer
     & ntij,   ntik,   ntil,   ntkl,   ntjl,   ntjk,   ntij1,
     & ntjk1,  ntkl1,  ijcxst,         ikcxst,         ilcxst,
     & klcxst,         jlcxst,         jkcxst
      common /cxindx/
     & ntij,   ntik,   ntil,   ntkl,   ntjl,   ntjk,   ntij1,
     & ntjk1,  ntkl1,  ijcxst(2),      ikcxst(2),      ilcxst(2),
     & klcxst(2),      jlcxst(2),      jkcxst(2)
c
c     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
c
      integer       l2rec
      equivalence ( l2rec, info(4) )
      integer       n2max
      equivalence ( n2max, info(5) )
c
      integer        ipq
      common /dim21/ ipq(256)
c
      integer       npriri,        nprirk
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
c
      integer    nipv
      parameter( nipv=4 )
c
c     # dummy:
       integer         ica(mcu,msu,*), icb(4,24,*),     icxst(2,*),
     &  il(*),         iprst(*),       ipt(*),          lmnp1(*),
     &  lmnv(3,*),     mcons(*),       nc(*),           ncon(*),
     &  nf(*),         nfct(*),        ngw(*),          npair(2,*),
     &  nrcr(*),       nt(*),          nprir(2,mstu,*), ntl(*),
     &  ntu(*),        labels(nipv,*), ibitv(*)
c
       real*8
     &  a(*),      cx(*),     eta(mrcru,mconu,*),  h4(*),
     &  x(mcu,*),  y(mcu,*),  z(mcu,*),            zet(mconu,*),
     &  buffer(*), values(*)
c     # ibm+convex f77 bug; actual:
c     # integer labels(nipv,n2max), ibitv( ((n2max+63)/64)*64 )
c     # real*8 buffer(l2rec), values(n2max)
c
c     # local:
       integer npkflg, isf, is, ifu, isfis, iksfis, ilsfis, ic,
     &  icu, ieqs1, jsf, js, jfu, ieqs2, jsfjs, ijsfjs, jlsfjs,
     &  ksf, ks, kfu, ieqs3, ksfks, iksfks, jksfks, lsf, ls, lfu,
     &  ncc, lsfls, klsfls, jlsfls, if, iksfif, ilsfif, ieqsf1,
     &  jf, jlsfjf, ieqsf2, kf, ieqsf3, lf, ibl, iscm, ist, nbli,
     &  nblk, jesfb, iijgt, iijx, iijy, iijz, kklgt, kklx, kkly,
     &  kklz, kblu, icc, itst, igwu, iprm, lc1, ngw1, igw, kc1,
     &  itst1, kc2, lc2, ic1, jc1, mxblu
c
      integer    iju(3),      ijb(2,3), klu(3),  klb(2,3), niqd(3),
     & nqd(3)
      integer    iperm(4,24), ip21(4),  ip31(6), ip41(12), ip51(4),
     & ip61(12), ip71(12),    ip75(3),  ip83(4), ip85(6)
c
      logical eijs, ejks, ekls
c
      character*8 chrtyp
c
      integer    nomore
      parameter( nomore = 2 )
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
*mdc*if parallel
*      integer me, nproc, ijkl
*      integer  nodeid, nnodes
*      external nodeid, nnodes
*mdc*endif
c
      data iperm /
     & 1,2,3,4, 2,1,3,4, 1,3,2,4, 3,1,2,4, 2,3,1,4, 3,2,1,4,
     & 1,2,4,3, 2,1,4,3, 1,4,2,3, 4,1,2,3, 2,4,1,3, 4,2,1,3,
     & 1,3,4,2, 3,1,4,2, 1,4,3,2, 4,1,3,2, 3,4,1,2, 4,3,1,2,
     & 2,3,4,1, 3,2,4,1, 2,4,3,1, 4,2,3,1, 3,4,2,1, 4,3,2,1 /
c
      data ip21 /1,7,13,19/
      data ip31 /1,3,9,5,11,17/
      data ip41 /1,3,5,7,9,11,13,15,17,19,21,23/
      data ip51 /10,4,2,1/
      data ip61 /4,2,1,10,8,7,16,14,13,22,20,19/
      data ip71 /17,11,9,18,5,3,12,6,1,10,4,2/
      data ip75 /9,3,1/
      data ip83 /1,2,7,8/
      data ip85 /1,3,7,9,13,15/
c
c     # initialize the integral record counters and the pk bit-vector.
c
      call w2stup( ibitv )
c
      itypea = 3
      itypeb = 0
c
      npkflg = 0
c
      nwt1(3) = 1
      isf = 0
c     # ijsf(1) is ijsf
      ijsf(1) = 0
c     # ijsf(2) is iksf
      ijsf(2) = 0
c     # ijsf(3) is ilsf
      ijsf(3) = 0
c
*mdc*if parallel
*c     # These are used to associate shell blocks to processes.
*      me = nodeid()
*      nproc = nnodes()
*      ijkl = me
*mdc*endif
c
      do 990 is = 1,ns
      ifu = nf(is)
      if(ifu.eq.0) go to 990
      isfis= isf
      iksfis = ijsf(2)
      ilsfis = ijsf(3)
      ic = nc(is)
      icu = nc(is)
      ieqs1 = 8
      jsf = 0
c     # klsf(2) is jlsf
      klsf(2) = 0
c     # klsf(3) is jksf
      klsf(3) = 0
      do 980 js = 1,is
      eijs=is.eq.js
      if(eijs) then
        ieqs1=4
      else
        jfu = nf(js)
        if(jfu.eq.0) go to 980
      endif
      ieqs2 = ieqs1
      jcu = nc(js)
      jsfjs = jsf
      ijsfjs = ijsf(1)
      jlsfjs = klsf(2)
      ksf = 0
      ijsf(2) = iksfis
c     # klsf(1) is klsf
      klsf(1) = 0
      do 970 ks = 1,js
      ejks=js.eq.ks
      if(ejks) then
        ieqs2=ieqs1-2
      else
        kfu = nf(ks)
        if(kfu.eq.0) go to 970
      endif
      ieqs3 = ieqs2
      kcu = nc(ks)
      ksfks= ksf
      iksfks = ijsf(2)
      jksfks = klsf(3)
      lsf = 0
      ijsf(3) = ilsfis
      klsf(2) = jlsfjs
      do 960 ls = 1,ks
      ekls=ks.eq.ls
      if(ekls) then
        ieqs3=ieqs2-1
      else
        lfu = nf(ls)
        if(lfu.eq.0) go to 960
      endif
      ncc = 0
      lcu = nc(ls)
      lsfls = lsf
      klsfls = klsf(1)
      jlsfls = klsf(2)
      call cecob( ic, ica, icb, ieqs3, is, js, ks, ls, ncc, nfct, ngw )
      isf = isfis
      ijsf(1) = ijsfjs
      ijsf(2) = iksfks
      do 930 if = 1, ifu
      isf = isf + 1
      ircru=nrcr(mcons(isf))
      itl = ntl(isf)
      itu = ntu(isf)
      iksfif = ijsf(2)
      ilsfif = ijsf(3)
      ieqsf1 = 8
      jsf = jsfjs
      klsf(2) = jlsfls
      klsf(3) = jksfks
      if(eijs) jfu=if
      do 930 jf = 1, jfu
      jlsfjf = klsf(2)
      jsf = jsf + 1
      jrcru=nrcr(mcons(jsf))
      jtl = ntl(jsf)
      jtu = ntu(jsf)
      ijsf(1) = ijsf(1) + 1
      iju(1) = nt(isf)*nt(jsf)
      ntij=iju(1)*npair(2,ijsf(1))
      esfi(1)=isf.eq.jsf
      if(esfi(1)) then
        ieqsf1=4
        ntij1=iju(1)*npair(1,ijsf(1))
        ijb(1,1)=ircru
      endif
      ijb(2,1)=ircru*jrcru
      ieqsf2 = ieqsf1
      ksf = ksfks
      ijsf(2) = iksfif
      klsf(1) = klsfls
      if(ejks) kfu=jf
      do 930 kf = 1, kfu
      ksf = ksf + 1
      krcru=nrcr(mcons(ksf))
      ktl = ntl(ksf)
      ktu = ntu(ksf)
      ijsf(2) = ijsf(2) + 1
      iju(2) = nt(isf)*nt(ksf)
      ntik=iju(2)*npair(2,ijsf(2))
      klsf(3) = klsf(3) + 1
      klu(3) = nt(jsf)*nt(ksf)
      ntjk=klu(3)*npair(2,klsf(3))
      esfi(2)=isf.eq.ksf
      if(esfi(2)) ijb(1,2)=ircru
      ijb(2,2)=ircru*krcru
      esfk(3)=jsf.eq.ksf
      if(esfk(3)) then
        ieqsf2=ieqsf1-2
        ntjk1=klu(3)*npair(1,klsf(3))
        klb(1,3)=jrcru
      endif
      klb(2,3)=jrcru*krcru
      ieqsf3 = ieqsf2
      lsf = lsfls
      ijsf(3) = ilsfif
      klsf(2) = jlsfjf
      if(ekls) lfu=kf
      do 920 lf = 1, lfu
      lsf = lsf + 1
      lrcru=nrcr(mcons(lsf))
      ltl = ntl(lsf)
      ltu = ntu(lsf)
      ijsf(3) = ijsf(3) + 1
      iju(3) = nt(isf)*nt(lsf)
      ntil=iju(3)*npair(2,ijsf(3))
      klsf(1) = klsf(1) + 1
      klu(1) = nt(ksf)*nt(lsf)
      ntkl=klu(1)*npair(2,klsf(1))
      klsf(2) = klsf(2) + 1
      klu(2) = nt(jsf)*nt(lsf)
      ntjl=klu(2)*npair(2,klsf(2))
      esfi(3)=isf.eq.lsf
      if(esfi(3)) ijb(1,3)=ircru
      ijb(2,3)=ircru*lrcru
      esfk(2)=jsf.eq.lsf
      if(esfk(2)) klb(1,2)=jrcru
      klb(2,2)=jrcru*lrcru
      esfk(1)=ksf.eq.lsf
      if(esfk(1)) then
        ieqsf3=ieqsf2-1
        ntkl1=klu(1)*npair(1,klsf(1))
        klb(1,1)=krcru
      endif
      klb(2,1)=krcru*lrcru
c
*mdc*if parallel
*      ijkl = ijkl + 1
*      if ( mod(ijkl,nproc) .ne. 0 ) goto 920
*mdc*endif
c
c     # set up arrays by permutation
      ibl = 0
      do 40 iscm = 1, 3
        if( (iscm .ne. 1) .and.
     &   (max( ijsf(iscm), klsf(iscm)) .eq. ijsf(iscm-1)) ) then
c         # (same block of integrals)
          nqd(iscm) = nqd(iscm-1)
          if(esfi(iscm-1)) then
            npri(1,iscm)=npri(1,iscm-1)
            do 22 ist=1,nst
              npriri(1,iscm,ist)=npriri(1,iscm-1,ist)
22          continue
          endif
          npri(2,iscm)=npri(2,iscm-1)
          do 24 ist=1,nst
            npriri(2,iscm,ist)=npriri(2,iscm-1,ist)
24        continue
          if ( esfk(iscm-1) ) then
            nprk(1,iscm)=nprk(1,iscm-1)
            do 26 ist=1,nst
              nprirk(1,iscm,ist)=nprirk(1,iscm-1,ist)
26          continue
          endif
          nprk(2,iscm)=nprk(2,iscm-1)
          do 28 ist=1,nst
            nprirk(2,iscm,ist)=nprirk(2,iscm-1,ist)
28        continue
          niqd(iscm) = 0
        else
c         # (new block of integrals)
          if ( esfi(iscm) ) then
            npri(1,iscm) = npair(1,ijsf(iscm))
            do 30 ist = 1, nst
              npriri(1,iscm,ist)=nprir(1,ist,ijsf(iscm))
30          continue
          endif
          npri(2,iscm) = npair(2,ijsf(iscm))
          do 32 ist = 1, nst
            npriri(2,iscm,ist) = nprir(2,ist,ijsf(iscm))
32        continue
          if ( esfk(iscm) ) then
            nprk(1,iscm) = npair(1,klsf(iscm))
            do 34 ist = 1, nst
              nprirk(1,iscm,ist) = nprir(1,ist,klsf(iscm))
34          continue
          endif
          nprk(2,iscm) = npair(2,klsf(iscm))
          do 36 ist = 1, nst
            nprirk(2,iscm,ist) = nprir(2,ist,klsf(iscm))
36        continue
          nblu = 0
          eprsf(iscm) = ijsf(iscm).eq.klsf(iscm)
          do 38 ist = 1, nst
            if ( esfi(iscm) ) then
              nbli = ijb(1,iscm) * npriri(1,iscm,ist)
     &         + ipq(ijb(1,iscm)) * npriri(2,iscm,ist)
            else
              nbli = ijb(2,iscm)*npriri(2,iscm,ist)
            endif
            if ( eprsf(iscm) ) then
              nblu = nblu + (nbli * (nbli + 1)) / 2
            else
              if ( esfk(iscm) ) then
                nblk = klb(1,iscm) * nprirk(1,iscm,ist)
     &           + ipq(klb(1,iscm)) * nprirk(2,iscm,ist)
              else
                nblk = klb(2,iscm) * nprirk(2,iscm,ist)
              endif
              nblu = nblu + nbli * nblk
            endif
38        continue
          nqd(iscm)  = nblu
          niqd(iscm) = nblu
        endif
40    continue
c
      do 60 jesfb = 1, 2
        ijcxst(jesfb) = icxst(jesfb,ijsf(1))
        ikcxst(jesfb) = icxst(jesfb,ijsf(2))
        ilcxst(jesfb) = icxst(jesfb,ijsf(3))
        klcxst(jesfb) = icxst(jesfb,klsf(1))
        jlcxst(jesfb) = icxst(jesfb,klsf(2))
        jkcxst(jesfb) = icxst(jesfb,klsf(3))
60    continue
c
      do 910 iscm = 1, 3
      if ( niqd(iscm) .eq. 0 ) go to 910
      iesfb = 2
      kesfb = 2
c
c     # allocate nblu workspace for this shell-triple component.
      nblu   = nqd(iscm)
c     # allocate space.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
c     #          4:buffer(1:l2rec), 5:values(1:n2max),
c     #          6:labels(1:4,1:n2max), 7:ibitv(1:n2max), 8:h4(1:nblu)
c     #          9:ijgt(1:iju(iscm)), 10:ijx(1:iju(iscm)),
c     #          11:ijy(1:iju(iscm)), 12:ijz(1:iju(iscm)),
c     #          13:klgt(1:klu(iscm)), 14:klx(1:klu(iscm)),
c     #          15:kly(1:klu(iscm)), 16:klz(1:klu(iscm)),
c     #          17:g4(1:ic3)
      ipt(9)  = ipt(8)  + atebyt( nblu )
      ipt(10) = ipt(9)  + forbyt( iju(iscm) )
      ipt(11) = ipt(10) + forbyt( iju(iscm) )
      ipt(12) = ipt(11) + forbyt( iju(iscm) )
      ipt(13) = ipt(12) + forbyt( iju(iscm) )
      ipt(14) = ipt(13) + forbyt( klu(iscm) )
      ipt(15) = ipt(14) + forbyt( klu(iscm) )
      ipt(16) = ipt(15) + forbyt( klu(iscm) )
      ipt(17) = ipt(16) + forbyt( klu(iscm) )
        if ( (ipt(17)-1) .gt. mblu ) then
          call bummer('twoint: mblu too small ',(ipt(17)-1-mblu),faterr)
        endif
c
c     # push the current alloction onto the high-water mark stack.
      call h2opsh( (ipt(17)-ipt(8)) )
c
      call wzero(nblu,h4,1)
c
      nwt = 1
      inx = 1
      jnx = 2
      knx = 3
      lnx = 4
      go to (100,200,300,400,500,600,700,800), ieqsf3
c
c  isf = jsf = ksf = lsf section
c
  100 do 110 icc = 1, ncc
c     # ieqsf3=1, ieqs3=1
  110       call sf3eq1( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        is,
     &       isf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
      go to 900
c
c  isf = jsf = ksf > lsf section
c
  200 do 299 icc=1,ncc
      if(ieqs3.ne.1) go to 220
c     # ieqsf3=2, ieqs3=1
      itst = 0
      igwu = ngw(icc)
      do 219 iprm = 1, 4
      inx = iperm(1,ip21(iprm))
      jnx = iperm(2,ip21(iprm))
      knx = iperm(3,ip21(iprm))
      lnx = iperm(4,ip21(iprm))
      lc1 = icb(lnx,1,icc)
      if(lc1.eq.itst) go to 219
      nwt = 1
      if(igwu.eq.1) go to 217
      ngw1 = 1
      do 216 igw=2,igwu
        if(icb(lnx,igw,icc)-lc1) 218, 214, 215
  214   ngw1 = ngw1 + 1
  215   nwt = nwt + 1
  216 continue
      nwt = nwt/ngw1
  217       call sf3eq2( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        is,
     &       isf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       ngw,        nrcr,       nt,
     &       ntl,        ntu,        x,          y,          z,
     &       zet )
  218 itst= lc1
  219 continue
      go to 299
c     # ieqsf3=2, ieqs3=2
  220       call sf3eq2( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        is,
     &       isf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       ngw,        nrcr,       nt,
     &       ntl,        ntu,        x,          y,          z,
     &       zet )
  299 continue
      go to 900
c
c  isf = jsf > ksf = lsf section
c
  300 do 399 icc = 1, ncc
      if(ieqs3.ne.1) go to 330
c     # ieqsf3=3, ieqs3=1
      itst = 0
      igwu = ngw(icc)
      do 319 iprm = 1, 6
      inx = iperm(1,ip31(iprm))
      jnx = iperm(2,ip31(iprm))
      knx = iperm(3,ip31(iprm))
      lnx = iperm(4,ip31(iprm))
      kc1 = icb(knx,1,icc)
      lc1 = icb(lnx,1,icc)
      itst1=ipq(kc1)+lc1
      if (itst1 .le. itst) go to 319
      nwt = 1
      if (igwu .eq. 1) go to 317
      ngw1 = 1
      do 316 igw = 2, igwu
        kc2 = icb(knx,igw,icc)
        lc2 = icb(lnx,igw,icc)
        if(ipq(max(kc2,lc2))+min(kc2,lc2)-itst1) 318, 314, 315
  314   ngw1 = ngw1 + 1
  315   nwt = nwt + 1
  316 continue
      nwt = nwt/ngw1
  317       call sf3eq3( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        ks,         ksf,        a(ipt(9)),
     &       a(ipt(10)), a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)),
     &       a(ipt(15)), a(ipt(16)), lmnp1,      lmnv,       mcons,
     &       nc,         ncon,       nfct,       ngw,        nrcr,
     &       nt,         ntl,        ntu,        x,          y,
     &       z,          zet )
  318 itst= itst1
  319 continue
      go to 399
c     # ieqsf3=3, ieqs3=3
  330       call sf3eq3( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        ks,         ksf,        a(ipt(9)),
     &       a(ipt(10)), a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)),
     &       a(ipt(15)), a(ipt(16)), lmnp1,      lmnv,       mcons,
     &       nc,         ncon,       nfct,       ngw,        nrcr,
     &       nt,         ntl,        ntu,        x,          y,
     &       z,          zet )
  399 continue
      go to 900
c
c  isf = jsf > ksf > lsf section
c
  400 do 499 icc = 1, ncc
      go to (410,420,430,440), ieqs3
c     # ieqsf3=4, ieqs3=1
  410 itst = 0
      igwu = ngw(icc)
      do 419 iprm = 1, 12
      inx = iperm(1,ip41(iprm))
      jnx = iperm(2,ip41(iprm))
      knx = iperm(3,ip41(iprm))
      lnx = iperm(4,ip41(iprm))
      lc1 = icb(lnx,1,icc)
      kc1 = icb(knx,1,icc)
      itst1 = kc1 + lc1 * icu
      if (itst1 .le. itst) go to 419
      nwt = 1
      if (igwu .eq. 1) go to 417
      ngw1 = 1
      do 416 igw = 2, igwu
        if(icb(knx,igw,icc)+icb(lnx,igw,icc)*icu-itst1) 418, 414, 415
  414   ngw1 = ngw1 + 1
  415   nwt = nwt + 1
  416 continue
      nwt = nwt/ngw1
  417       call sf3eq4( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        ks,         ksf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  418 itst= itst1
  419 continue
      go to 499
c     # ieqsf3=4, ieqs3=2; use first part of ip41 in place of ip42
  420 itst = 0
      igwu = ngw(icc)
      do 429 iprm = 1, 3
      inx = iperm(1,ip41(iprm))
      jnx = iperm(2,ip41(iprm))
      knx = iperm(3,ip41(iprm))
      kc1 = icb(knx,1,icc)
      if(kc1.eq.itst) go to 429
      nwt = 1
      if(igwu.eq.1) go to 427
      ngw1 = 1
      do 426 igw = 2, igwu
        if(icb(knx,igw,icc)-kc1) 428, 424, 425
  424   ngw1 = ngw1 + 1
  425   nwt = nwt + 1
  426 continue
      nwt = nwt /ngw1
  427       call sf3eq4( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        ks,         ksf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  428 itst = kc1
  429 continue
      go to 499
c     # ieqsf3=4, ieqs3=3; use first part of ip21 in place of ip43
  430 itst = 0
      igwu = ngw(icc)
      do 439 iprm = 1, 2
      knx = iperm(3,ip21(iprm))
      lnx = iperm(4,ip21(iprm))
      lc1 = icb (lnx,1,icc)
      if(lc1.eq.itst) go to 439
      nwt = 1
      if (igwu .eq. 1) go to 437
      ngw1 = 1
      do 436 igw = 2, igwu
        if(icb(lnx,igw,icc)-lc1) 438, 434, 435
  434   ngw1 = ngw1 + 1
  435   nwt = nwt + 1
  436 continue
      nwt = nwt /ngw1
  437       call sf3eq4( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        ks,         ksf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  438 itst = lc1
  439 continue
      go to 499
c     # ieqsf3=4, ieqs3=4
  440       call sf3eq4( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        ks,         ksf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  499 continue
      go to 900
c
c  isf > jsf = ksf = lsf section
c
  500 do 599 icc = 1, ncc
      if(ieqs3.ne.1) go to 550
c     # ieqsf3=5, ieqs3=1
      itst = 0
      igwu = ngw(icc)
      do 519 iprm = 1, 4
      inx = iperm(1,ip51(iprm))
      jnx = iperm(2,ip51(iprm))
      knx = iperm(3,ip51(iprm))
      lnx = iperm(4,ip51(iprm))
      ic1 = icb(inx,1,icc)
      if (ic1 .eq. itst) go to 519
      nwt = 1
      if (igwu .eq. 1) go to 517
      ngw1 = 1
      do 516 igw = 2, igwu
        if (icb(inx,igw,icc) - ic1) 518, 514, 515
  514   ngw1 = ngw1 + 1
  515   nwt = nwt + 1
  516 continue
      nwt = nwt/ngw1
  517       call sf3eq5( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        is,
     &       isf,        js,         jsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       ngw,        nrcr,       nt,
     &       ntl,        ntu,        x,          y,          z,
     &       zet )
  518 itst = ic1
  519 continue
      go to 599
c     # ieqsf3=5, ieqs3=5
  550       call sf3eq5( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        is,
     &       isf,        js,         jsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       ngw,        nrcr,       nt,
     &       ntl,        ntu,        x,          y,          z,
     &       zet )
  599 continue
      go to 900
c
c  isf > jsf = ksf > lsf section
c
  600 do 699 icc=1,ncc
      go to (610,620,601,601,650,660), ieqs3
  601 call bummer('stop at 601; ieqs3 = ',ieqs3,faterr)
c     # ieqsf3=6, ieqs3=1
  610 itst = 0
      igwu = ngw(icc)
      do 619 iprm = 1, 12
      inx = iperm(1,ip61(iprm))
      jnx = iperm(2,ip61(iprm))
      knx = iperm(3,ip61(iprm))
      lnx = iperm(4,ip61(iprm))
      lc1 = icb(lnx,1,icc)
      ic1 = icb(inx,1,icc)
      itst1 = ic1 + lc1 * icu
      if (itst1 .le. itst) go to 619
      nwt = 1
      if (igwu .eq. 1) go to 617
      ngw1 = 1
      do 616 igw = 2, igwu
        if (icb(inx,igw,icc)+icb(lnx,igw,icc)*icu-itst1) 618, 614, 615
  614   ngw1 = ngw1+ 1
  615   nwt = nwt + 1
  616 continue
      nwt = nwt/ngw1
  617       call sf3eq6( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  618 itst = itst1
  619 continue
      go to 699
c     # ieqsf3=6, ieqs3=2; use first part of ip61 in place of ip62
  620 itst = 0
      igwu = ngw(icc)
      do 629 iprm = 1, 3
      inx = iperm(1,ip61(iprm))
      jnx = iperm(2,ip61(iprm))
      knx = iperm(3,ip61(iprm))
      ic1 = icb(inx,1,icc)
      if (ic1 .eq. itst) go to 629
      nwt = 1
      if (igwu .eq. 1) go to 627
      ngw1 = 1
      do 626 igw = 2, igwu
        if (icb(inx,igw,icc) - ic1) 628, 624, 625
  624   ngw1 = ngw1 + 1
  625   nwt = nwt + 1
  626 continue
      nwt = nwt/ngw1
  627       call sf3eq6( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  628 itst = ic1
  629 continue
      go to 699
c     # ieqsf3=6, ieqs3=5; use first part of ip21 in place of ip65
  650 itst = 0
      igwu = ngw(icc)
      do 659 iprm = 1, 3
      jnx = iperm(2,ip21(iprm))
      knx = iperm(3,ip21(iprm))
      lnx = iperm(4,ip21(iprm))
      lc1 = icb(lnx,1,icc)
      if (lc1 .eq. itst) go to 659
      nwt = 1
      if (igwu .eq. 1) go to 657
      ngw1 = 1
      do 656 igw = 2, igwu
        if (icb(lnx,igw,icc) - lc1) 658, 654, 655
  654   ngw1 = ngw1 + 1
  655   nwt = nwt + 1
  656 continue
      nwt = nwt/ngw1
  657       call sf3eq6( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  658 itst = lc1
  659 continue
      go to 699
c     # ieqsf3=6, ieqs3=6
  660       call sf3eq6( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ls,
     &       lsf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  699 continue
      go to 900
c
c  isf > jsf > ksf = lsf section
c
  700 do 799 icc=1,ncc
      go to (710,701,730,701,750,701,770), ieqs3
  701 call bummer('stop at 701; ieqs3 = ',ieqs3,faterr)
c     # ieqsf3=7, ieqs3=1
  710 itst = 0
      igwu = ngw(icc)
      do 719 iprm = 1, 12
      inx = iperm(1,ip71(iprm))
      jnx = iperm(2,ip71(iprm))
      knx = iperm(3,ip71(iprm))
      lnx = iperm(4,ip71(iprm))
      jc1 = icb(jnx,1,icc)
      ic1 = icb(inx,1,icc)
      itst1 = ic1 + jc1 * icu
      if (itst1 .le. itst) go to 719
      nwt = 1
      if (igwu .eq. 1) go to 717
      ngw1 = 1
      do 716 igw = 2, igwu
        if(icb(inx,igw,icc)+icb(jnx,igw,icc)*icu-itst1) 718, 714, 715
  714   ngw1 = ngw1 + 1
  715   nwt = nwt + 1
  716 continue
      nwt = nwt/ngw1
  717       call sf3eq7( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  718 itst= itst1
  719 continue
      go to 799
c     # ieqsf3=7, ieqs3=3
  730 itst = 0
      igwu = ngw(icc)
      do 739 iprm = 1, 2
      inx = iperm(1,iprm)
      jnx = iperm(2,iprm)
      jc1 = icb(jnx,1,icc)
      if (jc1 .eq. itst) go to 739
      nwt = 1
      if (igwu .eq. 1) go to 737
      ngw1 = 1
      do 736 igw = 2, igwu
        if(icb(jnx,igw,icc) - jc1) 738, 734, 735
  734   ngw1 = ngw1 + 1
  735   nwt = nwt + 1
  736 continue
      nwt = nwt /ngw1
  737       call sf3eq7( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  738 itst = jc1
  739 continue
      go to 799
c     # ieqsf3=7, ieqs3=5
  750 itst = 0
      igwu = ngw(icc)
      do 759 iprm = 1, 3
      jnx = iperm(2,ip75(iprm))
      knx = iperm(3,ip75(iprm))
      lnx = iperm(4,ip75(iprm))
      jc1 = icb(jnx,1,icc)
      if (jc1 .eq. itst) go to 759
      nwt = 1
      if (igwu .eq. 1) go to 757
      ngw1 = 1
      do 756  igw = 2, igwu
        if (icb(jnx,igw,icc) - jc1) 758, 754, 755
  754   ngw1 = ngw1 + 1
  755   nwt = nwt + 1
  756 continue
      nwt = nwt/ngw1
  757       call sf3eq7( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  758 itst = jc1
  759 continue
      go to 799
c     # ieqsf3=7, ieqs3=7
  770       call sf3eq7( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        a(ipt(9)),  a(ipt(10)), a(ipt(11)), a(ipt(12)),
     &       a(ipt(13)), a(ipt(14)), a(ipt(15)), a(ipt(16)), lmnp1,
     &       lmnv,       mcons,      nc,         ncon,       nfct,
     &       ngw,        nrcr,       nt,         ntl,        ntu,
     &       x,          y,          z,          zet )
  799 continue
      go to 900
c
c  isf > jsf > ksf > lsf section
c
  800 do 899 icc=1,ncc
      go to (810,820,830,840,850,860,870,880), ieqs3
c     # ieqsf3=8, ieqs3=1
  810 itst = 0
      igwu = ngw(icc)
      do 819 iprm = 1, 24
      inx = iperm(1,iprm)
      jnx = iperm(2,iprm)
      knx = iperm(3,iprm)
      lnx = iperm(4,iprm)
      lc1 = icb(lnx,1,icc)
      kc1 = icb(knx,1,icc)
      jc1 = icb(jnx,1,icc)
      itst1 = jc1 + (kc1 + lc1 * icu) * icu
      if (itst1 .le. itst) go to 819
      nwt = 1
      do 816 igw = 2, igwu
        if(icb(jnx,igw,icc)+(icb(knx,igw,icc)+icb(lnx,igw,icc)*icu)*icu
     1    .lt.itst1) go to 818
        nwt = nwt + 1
  816 continue
            call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
  818 itst = itst1
  819 continue
      go to 899
c     # ieqsf3=8, ieqs3=2
  820 itst = 0
      igwu = ngw(icc)
      do 829 iprm = 1, 6
      inx = iperm(1,iprm)
      jnx = iperm(2,iprm)
      knx = iperm(3,iprm)
      kc1 = icb(knx,1,icc)
      jc1 = icb(jnx,1,icc)
      itst1 = jc1 + kc1 * icu
      if (itst1 .le. itst) go to 829
      nwt = 1
      do 826 igw = 2, igwu
        if(icb(jnx,igw,icc)+icb(knx,igw,icc)*icu.lt.itst1) go to 828
        nwt = nwt + 1
  826 continue
            call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
  828 itst = itst1
  829 continue
      go to 899
c     # ieqsf3=8, ieqs3=3
  830 itst = 0
      igwu = ngw(icc)
      do 839 iprm = 1, 4
      inx = iperm(1,ip83(iprm))
      jnx = iperm(2,ip83(iprm))
      knx = iperm(3,ip83(iprm))
      lnx = iperm(4,ip83(iprm))
      lc1 = icb(lnx,1,icc)
      jc1 = icb(jnx,1,icc)
      itst1 = jc1 + lc1 * icu
      if (itst1 .le. itst) go to 839
      nwt = 1
      do 836 igw = 2, igwu
        if(icb(jnx,igw,icc)+icb(lnx,igw,icc)*icu.lt.itst1) go to 838
        nwt = nwt + 1
  836 continue
            call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
  838 itst = itst1
  839 continue
      go to 899
c     # ieqsf3=8, ieqs3=4
  840 itst = 0
      igwu = ngw(icc)
      do 849 iprm = 1, 2
      inx = iperm(1,iprm)
      jnx = iperm(2,iprm)
      jc1 = icb(jnx,1,icc)
      if (jc1 .eq. itst) go to 849
      nwt = 1
      if (igwu .eq. 1) go to 847
      if(icb(jnx,2,icc).lt.jc1) go to 849
      nwt = 2
  847       call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
      itst = jc1
  849 continue
      go to 899
c     # ieqsf3=8, ieqs3=5
  850 itst = 0
      igwu = ngw(icc)
      do 859 iprm = 1, 6
      jnx = iperm(2,ip85(iprm))
      knx = iperm(3,ip85(iprm))
      lnx = iperm(4,ip85(iprm))
      lc1 = icb(lnx,1,icc)
      kc1 = icb(knx,1,icc)
      itst1 = kc1 + lc1 * jcu
      if (itst1 .le. itst) go to 859
      nwt = 1
      if (igwu .eq. 1) go to 857
      do 856 igw = 2, igwu
        if(icb(knx,igw,icc)+icb(lnx,igw,icc)*jcu.lt.itst1) go to 858
        nwt = nwt + 1
  856 continue
  857       call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
  858 itst = itst1
  859 continue
      go to 899
c     # ieqsf3=8, ieqs3=6; use first part of ip85 in place of ip86
  860 itst = 0
      igwu = ngw(icc)
      do 869 iprm = 1, 2
      jnx = iperm(2,ip85(iprm))
      knx = iperm(3,ip85(iprm))
      kc1 = icb(knx,1,icc)
      if (kc1 .eq. itst) go to 869
      nwt = 1
      if (igwu .eq. 1) go to 867
      if(icb(knx,2,icc).lt.kc1) go to 869
      nwt = 2
  867       call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
      itst = kc1
  869 continue
      go to 899
c     # ieqsf3=8, ieqs3=7; use first part of ip21 in place of ip87
  870 itst = 0
      igwu = ngw(icc)
      do 879 iprm = 1, 2
      knx = iperm(3,ip21(iprm))
      lnx = iperm(4,ip21(iprm))
      lc1 = icb (lnx,1,icc)
      if (lc1 .eq. itst) go to 879
      nwt = 1
      if (igwu .eq. 1) go to 877
      if(icb(lnx,2,icc).lt.lc1) go to 879
      nwt = 2
  877       call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
      itst = lc1
  879 continue
      go to 899
c     # ieqsf3=8, ieqs3=8
  880       call sf3eq8( a,          cx,         eta,        a(ipt(17)),
     &       h4,         icb,        icc,        ipt,        iscm,
     &       is,         isf,        js,         jsf,        ks,
     &       ksf,        ls,         lsf,        a(ipt(9)),  a(ipt(10)),
     &       a(ipt(11)), a(ipt(12)), a(ipt(13)), a(ipt(14)), a(ipt(15)),
     &       a(ipt(16)), lmnp1,      lmnv,       mcons,      nc,
     &       ncon,       nfct,       nrcr,       nt,         ntl,
     &       ntu,        x,          y,          z,          zet )
  899 continue
  900 continue
c
c     # write out the integrals for this component of the shell-triple.
            call wtint2( h4,     il,     iprst,  iscm,   buffer, values,
     &           labels, ibitv )
c
c     # reclaim kblu workspace, and pop back to the intitial level.
c     mblu = mblu + kblu
      call h2opop
c
  910 continue
c     # the complete set of integrals has been calculated for the
c     # construction of a disjoint subset of p and k elements.
c     # set the pk flag.
      if ( ibuf .ne. 0 ) then
        npkflg = npkflg - ibitv(ibuf) + 1
        ibitv(ibuf) = +1
      endif
  920 continue
  930 continue
  960 continue
  970 continue
  980 continue
  990 continue
c
c     # write out the last 2-e integral record.
c
      call wtlab2( nomore, buffer, values, labels, ibitv )
c
      call siftyp( itypea, itypeb, chrtyp )
      write( nlist, 6010 ) numout, chrtyp, npkflg, nrec
c
c     # get the high-water mark for a(*) workspace usage.
      call h2omtr( mxblu )
      write( nlist, '(a,i9)' ) ' twoint: maximum mblu needed = ',mxblu
c
      return
6010  format(/' twoint:',i12,1x,a,' integrals and',i9,
     & ' pk flags were written in',i6,' records.'/)
      end
*deck cecob
      subroutine cecob(
     & ic, ica, icb, ieqs3, is, js, ks, ls, ncc, nfct, ngw )
c
      implicit integer(a-z)
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1  jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2  npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      dimension ica(mcu,msu,*), icb(4,24,*), nfct(*), ngw(*)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      go to (100,200,300,400,500,600,700,800), ieqs3
c
c  center combinations for is=js=ks=ls
c
  100 do 170 jc = 1, ic
      do 165 kc = 1, jc
      do 160 lc = 1, kc
      igw = 0
      do 140 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      kc1 = ica(kc,is,ig)
      lc1 = ica(lc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic.and.kc1.ne.ic.and.lc1.ne.ic) go to 140
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      kc2 = max(kc1,lc1)
      lc2 = min(kc1,lc1)
      ic3 = max(ic2,kc2)
      kc2 = min(ic2,kc2)
      lc3 = min(jc2,lc2)
      jc2 = max(jc2,lc2)
      jc3 = max(jc2,kc2)
      kc3 = min(jc2,kc2)
      do 115 icc = ncc, 1, -1
        if(jc3 - icb(2,1,icc)) 115, 105, 120
  105   if(kc3 - icb(3,1,icc)) 115, 110, 120
  110   if(lc3 - icb(4,1,icc)) 115, 145, 120
  115 continue
  120 if(jc3 .ne. jc .or. kc3 .ne. kc .or. lc3 .ne. lc) go to 140
      do 130 igw1 = igw, 1, -1
        if (ic1.eq.icb(1,igw1,ncc+1).and.
     &      jc1.eq.icb(2,igw1,ncc+1).and.
     &      kc1.eq.icb(3,igw1,ncc+1).and.
     &      lc1.eq.icb(4,igw1,ncc+1)) go to 140
  130 continue
      igw = igw + 1
      icb(1,igw,ncc+1) = ic1
      icb(2,igw,ncc+1) = jc1
      icb(3,igw,ncc+1) = kc1
      icb(4,igw,ncc+1) = lc1
  140 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 150
  145 nfct(icc) = nfct(icc) + ic
  150 continue
  160 continue
  165 continue
  170 continue
      do 190  icc = 1, ncc
        jc1  = icb(2,1,icc)
        kc1  = icb(3,1,icc)
        lc1  = icb(4,1,icc)
        idc = 1
        if (jc1 .ne. ic) idc = 2
        if (kc1 .ne. jc1) idc = idc + 1
        if (lc1 .ne. kc1) idc = idc + 1
        nfct(icc) = nfct(icc)/idc
  190 continue
      go to 900
c
c  center combinations for is=js=ks>ls
c
  200 do 270 jc = 1, ic
      do 265 kc = 1, jc
      do 260 lc = 1, lcu
      igw = 0
      do 240 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      kc1 = ica(kc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic.and.kc1.ne.ic) go to 240
      lc1 = ica(lc,ls,ig)
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      ic3 = max(ic2,kc1)
      kc2 = min(ic2,kc1)
      jc3 = max(jc2,kc2)
      kc3 = min(jc2,kc2)
      do 215 icc = ncc, 1, -1
        if(jc3 - icb(2,1,icc)) 215, 205, 220
  205   if(kc3 - icb(3,1,icc)) 215, 210, 220
  210   if(lc1 - icb(4,1,icc)) 215, 245, 220
  215 continue
  220 if(jc3 .ne. jc .or. kc3 .ne. kc .or. lc1 .ne. lc) go to 240
      do 230 igw1 = igw, 1, -1
        if (ic1.eq.icb(1,igw1,ncc+1).and.
     &      jc1.eq.icb(2,igw1,ncc+1).and.
     &      kc1.eq.icb(3,igw1,ncc+1)) go to 240
  230 continue
      igw = igw + 1
      icb(1,igw,ncc+1) = ic1
      icb(2,igw,ncc+1) = jc1
      icb(3,igw,ncc+1) = kc1
      icb(4,igw,ncc+1) = lc1
  240 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 260
  245 nfct(icc) = nfct(icc) + ic
  260 continue
  265 continue
  270 continue
      do 290 icc = 1, ncc
        jc1  = icb(2,1,icc)
        kc1  = icb(3,1,icc)
        idc = 1
        if ( jc1 .ne. ic) idc = 2
        if (kc1 .ne. jc1) idc = idc + 1
        nfct(icc) = nfct(icc)/idc
  290 continue
      go to 900
c
c  center combinations for is=js>ks=ls
c
  300 do 370 jc = 1, ic
      do 365 kc = 1, kcu
      do 360 lc = 1, kc
      igw = 0
      do 340 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic) go to 340
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ks,ig)
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      kc2 = max(kc1,lc1)
      lc2 = min(kc1,lc1)
      do 315 icc = ncc, 1, -1
        if (jc2 - icb(2,1,icc)) 315, 305, 320
  305   if (kc2 - icb(3,1,icc)) 315, 310, 320
  310   if (lc2 - icb(4,1,icc)) 315, 345, 320
  315 continue
  320 if (jc2 .ne. jc .or. kc2 .ne. kc .or. lc2 .ne. lc) go to 340
      do 330 igw1 = igw, 1, -1
        if( ic1.eq.icb(1,igw1,ncc+1).and.
     &      jc1.eq.icb(2,igw1,ncc+1).and.
     &      kc1.eq.icb(3,igw1,ncc+1).and.
     &      lc1.eq.icb(4,igw1,ncc+1)) go to 340
  330 continue
      igw = igw + 1
      icb (1,igw,ncc+1) = ic1
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
  340 continue
      ncc = ncc + 1
      nfct(ncc) =  ic
      ngw(ncc) = igw
      go to 360
  345 nfct(icc) = nfct(icc) + ic
  360 continue
  365 continue
  370 continue
      go to 480
c
c  center combinations for is=js>ks>ls
c
  400 do 470 jc = 1, ic
      do 465 kc = 1, kcu
      do 460 lc = 1, lcu
      igw = 0
      do 440 ig = 1, ng
      ic1 = ica(ic,is,ig)
      jc1 = ica(jc,is,ig)
      if(ic1.ne.ic.and.jc1.ne.ic) go to 440
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ls,ig)
      ic2 = max(ic1,jc1)
      jc2 = min(ic1,jc1)
      do 415 icc = ncc, 1, -1
        if (jc2 - icb(2,1,icc)) 415, 405, 420
  405   if (kc1 - icb(3,1,icc)) 415, 410, 420
  410   if (lc1 - icb(4,1,icc)) 415, 445, 420
  415 continue
  420 if(jc2 .ne. jc .or. kc1 .ne. kc .or. lc1 .ne. lc) go to 440
      do 430 igw1 = igw, 1, -1
        if (ic1.eq.icb(1,igw1,ncc+1).and.
     &      jc1.eq.icb(2,igw1,ncc+1)) go to 440
  430 continue
      igw = igw + 1
      icb(1,igw,ncc+1) = ic1
      icb(2,igw,ncc+1) = jc1
      icb(3,igw,ncc+1) = kc1
      icb(4,igw,ncc+1) = lc1
  440 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 460
  445 nfct(icc) = nfct(icc) + ic
  460 continue
  465 continue
  470 continue
  480 do 490 icc = 1, ncc
        if(icb(2,1,icc).ne.ic) nfct(icc)=nfct(icc)/2
  490 continue
      go to 900
c
c  center combinations for is>js=ks=ls
c
  500 do 570 jc = 1, jcu
      do 565 kc = 1, jc
      do 560 lc = 1, kc
      igw =0
      do 540 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 540
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,js,ig)
      lc1 = ica(lc,js,ig)
      jc2 = max (jc1,kc1)
      kc2 = min (jc1,kc1)
      jc3 = max (jc2,lc1)
      lc2 = min (jc2,lc1)
      kc3 = max (kc2,lc2)
      lc3 = min (kc2,lc2)
      do 515 icc = ncc, 1, -1
      if (jc3 - icb(2,1,icc)) 515, 505, 520
  505 if (kc3 - icb(3,1,icc)) 515, 510, 520
  510 if (lc3 - icb(4,1,icc)) 515, 545, 520
  515 continue
  520 if (jc3.ne.jc.or.kc3.ne.kc.or.lc3.ne.lc) go to 540
      do 530 igw1 = igw, 1, -1
        if (jc1.eq.icb(2,igw1,ncc+1).and.
     &      kc1.eq.icb(3,igw1,ncc+1).and.
     &      lc1.eq.icb(4,igw1,ncc+1)) go to 540
  530 continue
      igw = igw + 1
      icb (1,igw,ncc+1) = ic
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
  540 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 560
  545 nfct(icc) = nfct(icc) + ic
  560 continue
  565 continue
  570 continue
      go to 900
c
c  center combinations for is>js=ks>ls
c
  600 do 670 jc = 1, jcu
      do 665 kc = 1, jc
      do 660 lc = 1, lcu
      igw = 0
      do 640 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 640
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,js,ig)
      lc1 = ica(lc,ls,ig)
      jc2 = max (jc1, kc1)
      kc2 = min (jc1, kc1)
      do 615 icc= ncc, 1, -1
        if (jc2 - icb(2,1,icc)) 615, 605, 620
  605   if (kc2 - icb(3,1,icc)) 615, 610, 620
  610   if (lc1 - icb(4,1,icc)) 615, 645, 620
  615 continue
  620 if (jc2.ne.jc.or.kc2.ne.kc.or.lc1.ne.lc) go to 640
      do 630 igw1 = igw, 1, -1
        if (jc1.eq.icb(2,igw1,ncc+1).and.
     &      kc1.eq.icb(3,igw1,ncc+1)) go to 640
  630 continue
      igw = igw + 1
      icb (1,igw,ncc+1) = ic
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
  640 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw (ncc) = igw
      go to 660
  645 nfct(icc) = nfct(icc) + ic
  660 continue
  665 continue
  670 continue
      go to 900
c
c  center combinations for is>js>ks=ls
c
  700 do 770 jc = 1, jcu
      do 765 kc = 1, kcu
      do 760 lc = 1, kc
      igw = 0
      do 740 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 740
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ks,ig)
      kc2 = max (kc1, lc1)
      lc2 = min (kc1, lc1)
      do 715 icc = ncc, 1, -1
        if (jc1 - icb(2,1,icc)) 715, 705, 720
  705   if (kc2 - icb(3,1,icc)) 715, 710, 720
  710   if (lc2 - icb(4,1,icc)) 715, 745, 720
  715 continue
  720 if (jc1.ne.jc.or.kc2.ne.kc.or.lc2.ne.lc) go to 740
      do 730 igw1 = igw, 1, -1
        if(kc1.eq.icb(3,igw1,ncc+1).and.
     &     lc1.eq.icb(4,igw1,ncc+1)) go to 740
  730 continue
      igw = igw + 1
      icb (1,igw,ncc+1) = ic
      icb (2,igw,ncc+1) = jc1
      icb (3,igw,ncc+1) = kc1
      icb (4,igw,ncc+1) = lc1
  740 continue
      ncc = ncc + 1
      nfct(ncc) = ic
      ngw(ncc) = igw
      go to 760
  745 nfct(icc) = nfct(icc) + ic
  760 continue
  765 continue
  770 continue
      go to 900
c
c  center combinations for is>js>ks>ls
c
  800 do 870 jc = 1, jcu
      do 865 kc = 1, kcu
      do 860 lc = 1, lcu
      do 840 ig = 1, ng
      if (ica(ic,is,ig) .ne. ic) go to 840
      jc1 = ica(jc,js,ig)
      kc1 = ica(kc,ks,ig)
      lc1 = ica(lc,ls,ig)
      do 815 icc = ncc, 1, -1
        if (jc1 - icb(2,1,icc)) 815, 805, 840
  805   if (kc1 - icb(3,1,icc)) 815, 810, 840
  810   if (lc1 - icb(4,1,icc)) 815, 845, 840
  815 continue
  840 continue
      ncc = ncc + 1
      icb (1, 1, ncc) = ic
      icb (2, 1, ncc) = jc
      icb (3, 1, ncc) = kc
      icb (4, 1, ncc) = lc
      nfct(ncc) = ic
      go to 860
  845 nfct(icc) = nfct(icc) + ic
  860 continue
  865 continue
  870 continue
  900 if(ncc.gt.mccu) then
        write (nlist,*) 'cecob: mccup too small', ncc
        call bummer('cecob: ncc=',ncc,faterr)
      endif
      return
      end
*deck sf3eq1
      subroutine sf3eq1( a,        cx,       eta,       g4,
     &         h4,       icb,      icc,      ipt,       is,
     &         isf,      ijgt,     ijx,      ijy,       ijz,
     &         klgt,     klx,      kly,      klz,       lmnp1,
     &         lmnv,     mcons,    nc,       ncon,      nfct,
     &         ngw,      nrcr,     nt,       ntl,       ntu,
     &         x,        y,        z,        zet )
      implicit real*8 (a-h,o-z)
      logical eij, eij1, ekl, esfc, esfcij, esfckl
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1  jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2  npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /cxindx/ ntij,ntik,ntil,ntkl,ntjl,ntjk,ntij1,ntjk1,ntkl1,
     1  ijcxst(2),ikcxst(2),ilcxst(2),klcxst(2),jlcxst(2),jkcxst(2)
      common /dim21/ ipq(256)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*),cx(*),eta(mrcru,mconu,*),g4(*),h4(*),icb(4,24,*),
     1  ipt(*),ijgt(*),ijx(*),ijy(*),ijz(*),klgt(*),klx(*),kly(*),
     2  klz(*),lmnp1(*),lmnv(3,*),mcons(*),nc(*),ncon(*),nfct(*),ngw(*),
     3  nrcr(*),nt(*),ntl(*),ntu(*),x(mcu,*),y(mcu,*),z(mcu,*),
     4  zet(mconu,*)
      dimension esfc(3),esfcij(3),esfckl(3),icxbg(3,3),jbl(6),kcxbg(3,3)
c
      ic=icb(1,1,icc)
      jc=icb(2,1,icc)
      kc=icb(3,1,icc)
      lc=icb(4,1,icc)
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(1,2)=icxbg(1,1)+(kc-jc)*ntij1
      icxbg(1,3)=icxbg(1,2)+(lc-kc)*ntij1
      kcxbg(1,1)=ijcxst(1)+ipq((kc-1)*nt(isf)+1)*npri(1,1)+(lc-1)*ntij1
      kcxbg(1,2)=ijcxst(1)+ipq((jc-1)*nt(isf)+1)*npri(1,1)+(lc-1)*ntij1
      kcxbg(1,3)=kcxbg(1,2)+(kc-lc)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*((nc(is)-1)*ntij)
      icxbg(2,2)=icxbg(2,1)+(kc-jc)*ntij
      icxbg(3,2)=icxbg(2,2)-(ic-kc)*((nc(is)-1)*ntij)
      icxbg(2,3)=icxbg(2,2)+(lc-kc)*ntij
      icxbg(3,3)=icxbg(2,3)-(ic-lc)*((nc(is)-1)*ntij)
      kcxbg(2,1)=icxbg(2,3)+(kc-ic)*(nc(is)*ntij)
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*((nc(is)-1)*ntij)
      kcxbg(2,2)=kcxbg(2,1)+(jc-kc)*(nc(is)*ntij)
      kcxbg(3,2)=kcxbg(2,2)-(jc-lc)*((nc(is)-1)*ntij)
      kcxbg(2,3)=kcxbg(2,2)+(kc-lc)*ntij
      kcxbg(3,3)=kcxbg(2,3)-(jc-kc)*((nc(is)-1)*ntij)
      esfcij(1)=ic.eq.jc
      esfckl(3)=jc.eq.kc
      esfckl(1)=kc.eq.lc
      esfcij(2)=ic.eq.kc
      esfckl(2)=jc.eq.lc
      esfcij(3)=ic.eq.lc
      esfc(1)=esfcij(2).and.esfckl(2)
      esfc(2)=esfcij(1).and.esfckl(1)
      esfc(3)=esfc(2)
      call izero_argos(6,jbl,1)
      do 10 ist=1,nst
        jbl(1)=jbl(1)+ipq(npriri(1,1,ist)+1)
        jbl(3)=jbl(3)+ipq(npriri(2,1,ist)+1)
        jbl(4)=jbl(4)+npriri(1,1,ist)**2
        jbl(5)=jbl(5)+npriri(2,1,ist)*npriri(1,1,ist)
        jbl(6)=jbl(6)+npriri(2,1,ist)**2
   10 continue
      if(esfcij(1).or.esfckl(1)) then
        nwt1(1) = 1
        if(esfckl(3)) then
          iscmu = 1
          go to 30
        else
          iscmu = 2
          go to 28
        endif
      endif
      if(esfckl(3)) then
        iscml = 2
        iscmu = 3
        nwt1(2) = 1
        go to 32
      endif
      if(ngw(icc).eq.1) go to 26
      isw2 = 2
      isw3 = 3
      do 18 igw=2,ngw(icc)
        icb1 = icb(inx,igw,icc)
        icb2 = icb(knx,igw,icc)
        itst1 = max(icb1, icb2)
        itst2 = min(icb1, icb2)
        if((itst1.eq.ic.and.itst2.eq.jc).or.
     &     (itst1.eq.kc.and.itst2.eq.lc)) isw2 = 1
        icb2 = icb(lnx,igw,icc)
        itst1 = max(icb1, icb2)
        itst2 = min(icb1, icb2)
        if((itst1.eq.ic.and.itst2.eq.kc).or.
     &     (itst1.eq.jc.and.itst2.eq.lc)) isw3 = min(isw3,2)
        if((itst1.eq.ic.and.itst2.eq.jc).or.
     &     (itst1.eq.kc.and.itst2.eq.lc)) isw3 = 1
   18 continue
      if(isw3.eq.1) then
        if(isw2.eq.1) then
          iscmu = 1
          nwt1(1) = 3
          go to 30
        else
          iscmu = 2
          nwt1(1) = 2
          go to 28
        endif
      elseif(isw3.eq.2) then
        iscmu = 2
        nwt1(1) = 1
        nwt1(2) = 2
        go to 30
      elseif(isw2.eq.1) then
        iscml = 2
        iscmu = 3
        nwt1(2) = 2
        go to 32
      endif
   26 iscmu = 3
      nwt1(1) = 1
   28 nwt1(2) = 1
   30 iscml = 1
   32 ngt=nt(isf)*(nt(isf)*(nt(isf)*nt(isf)))
      do 98 iscm=iscml,iscmu
      fnfct=nwt1(iscm)*nfct(icc)
      if(iscm-2) 34, 35, 36
   34 call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,jc,is,isf,kc,is,isf,lc,
     &  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     &  nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 38
   35 call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,kc,is,isf,jc,is,isf,lc,
     &  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     &  nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 38
   36 call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,lc,is,isf,jc,is,isf,kc,
     &  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     &  nrcr,nt,ntl,ntu,x,y,z,zet)
   38 ndx=1
      do 96 ircr=1,ircru
      if(esfc(iscm)) krcru=ircr
      if(esfcij(iscm)) jrcru=ircr
      do 94 jrcr=1,jrcru
      if(ircr-jrcr) 50, 52, 54
   50 ijrcr=ipq(jrcr)+ircr
      ije=jrcr-1
      iesfb1=2
      icxbg1=icxbg(3,iscm)
      ngti1=(nt(isf)*nt(isf))
      ngtj1=(nt(isf)*(nt(isf)*nt(isf)))
      go to 58
   52 ijrcr=ipq(ircr+1)
      ije=ircr
      iesfb1=1
      go to 56
   54 ijrcr=ipq(ircr)+jrcr
      ije=ircr-1
      iesfb1=2
   56 icxbg1=icxbg(iesfb1,iscm)
      ngti1=(nt(isf)*(nt(isf)*nt(isf)))
      ngtj1=(nt(isf)*nt(isf))
   58 ijn=ijrcr-ije
      eij1=esfcij(iscm).and.ircr.eq.jrcr
      do 92 krcr=1,krcru
      if(esfc(iscm).and.ircr.eq.krcr) then
        lrcru = jrcr
      elseif(esfckl(iscm)) then
        lrcru = krcr
      elseif(esfc(iscm)) then
        lrcru = ircru
      endif
      do 90 lrcr=1,lrcru
      if(krcr-lrcr) 60, 62, 64
   60 klrcr=ipq(lrcr)+krcr
      kle=lrcr-1
      kesfb=2
      kcxbg1=kcxbg(3,iscm)
      ngtk1=1
      ngtl1=nt(isf)
      go to 68
   62 klrcr=ipq(krcr+1)
      kle=krcr
      kesfb=1
      go to 66
   64 klrcr=ipq(krcr)+lrcr
      kle=krcr-1
      kesfb=2
   66 kcxbg1=kcxbg(kesfb,iscm)
      ngtk1=nt(isf)
      ngtl1=1
   68 kln=klrcr-kle
      ekl=esfckl(iscm).and.krcr.eq.lrcr
      if(ijrcr-klrcr) 70, 72, 74
   70 if(kesfb.eq.1) then
        iblb=(kle-1)*jbl(1)+kln*jbl(3)+(ipq(kle-1)+ije)*jbl(4)+(kln*
     1    (kle-1)+ijn)*jbl(5)+ipq(kln)*jbl(6)-jbl(iesfb1+3)
      else
        iblb=kle*jbl(1)+(kln-1)*jbl(3)+ipq(kle)*jbl(4)+((kln-1)*kle+
     1    ije)*jbl(5)+(ipq(kln)-kln+1+ijn)*jbl(6)-jbl(iesfb1+4)
      endif
      iesfb=kesfb
      kesfb=iesfb1
      icxs(iesfb,iscm)=kcxbg1
      kcxbg1=icxbg1
      eij=ekl
      ekl=eij1
      ngti=ngtk1
      ngtj=ngtl1
      ngtk=ngti1
      ngtl=ngtj1
      go to 78
   72 iblb=ije*jbl(1)+ijn*jbl(3)+ipq(ije)*jbl(4)+ijn*ije*jbl(5)+
     1  (ipq(ijn+1)-ijn)*jbl(6)-jbl(iesfb1+kesfb-1)
      go to 76
   74 if(iesfb1.eq.1) then
        iblb=(ije-1)*jbl(1)+ijn*jbl(3)+(ipq(ije-1)+kle)*jbl(4)+(ijn*
     1    (ije-1)+kln)*jbl(5)+ipq(ijn)*jbl(6)-jbl(kesfb+3)
      else
        iblb=ije*jbl(1)+(ijn-1)*jbl(3)+ipq(ije)*jbl(4)+((ijn-1)*ije+
     1    kle)*jbl(5)+(ipq(ijn)-ijn+1+kln)*jbl(6)-jbl(kesfb+4)
      endif
   76 iesfb=iesfb1
      icxs(iesfb,iscm)=icxbg1
      eij=eij1
      ngti=ngti1
      ngtj=ngtj1
      ngtk=ngtk1
      ngtl=ngtl1
   78 ndxi=ndx
      do 88 it=itl,itu
        if(((esfc(iscm).and.ircr.eq.krcr).and.jrcr.eq.lrcr)) ktu = it
        ndxj=ndxi
        if(eij) jtu = it
        do 86 jt=itl,jtu
          ndxk=ndxj
          kcxs(kesfb,iscm)=kcxbg1
          do 84 kt=itl,ktu
            ndxl=ndxk
            if(esfc(iscm).and.ircr.eq.krcr.and.jrcr.eq.lrcr.and.
     &         it.eq.kt) then
              ltu = jt
            elseif(ekl) then
              ltu = kt
            elseif(esfc(iscm).and.ircr.eq.krcr.and.jrcr.eq.lrcr) then
              ltu = itu
            endif
            do 82 lt=itl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              if(ijrcr.eq.klrcr) then
                call aoso2e( h4, cx, iscm, val )
              else
                call aoso2( h4, cx, iscm, val )
              endif
   80         ndxl=ndxl+ngtl
              kcxs(kesfb,iscm)=kcxs(kesfb,iscm)+npri(kesfb,1)
   82       continue
            ndxk=ndxk+ngtk
   84     continue
          ndxj=ndxj+ngtj
          icxs(iesfb,iscm)=icxs(iesfb,iscm)+npri(iesfb,1)
   86   continue
        ndxi=ndxi+ngti
   88 continue
      ndx=ndx+ngt
   90 continue
   92 continue
   94 continue
   96 continue
   98 continue
      return
      end
*deck sf3eq2
      subroutine sf3eq2( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      is,
     &         isf,      ls,       lsf,      ijgt,     ijx,
     &         ijy,      ijz,      klgt,     klx,      kly,
     &         klz,      lmnp1,    lmnv,     mcons,    nc,
     &         ncon,     nfct,     ngw,      nrcr,     nt,
     &         ntl,      ntu,      x,        y,        z,
     &         zet )
      implicit real*8 (a-h,o-z)
      logical ec
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1  jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2  npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /cxindx/ ntij,ntik,ntil,ntkl,ntjl,ntjk,ntij1,ntjk1,ntkl1,
     1  ijcxst(2),ikcxst(2),ilcxst(2),klcxst(2),jlcxst(2),jkcxst(2)
      common /dim21/ ipq(256)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*),cx(*),eta(mrcru,mconu,*),g4(*),h4(*),icb(4,24,*),
     1  ipt(*),ijgt(*),ijx(*),ijy(*),ijz(*),klgt(*),klx(*),kly(*),
     2  klz(*),lmnp1(*),lmnv(3,*),mcons(*),nc(*),ncon(*),nfct(*),ngw(*),
     3  nrcr(*),nt(*),ntl(*),ntu(*),x(mcu,*),y(mcu,*),z(mcu,*),
     4  zet(mconu,*)
      dimension ec(3), icxbg(3,3), jbl(2), kcxbg(3)
c
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(1,2)=icxbg(1,1)+(kc-jc)*ntij1
      kcxbg(3)=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg(1)=kcxbg(3)+(kc-ic)*(nc(ls)*ntil)
      kcxbg(2)=kcxbg(1)+(jc-kc)*(nc(ls)*ntil)
      icxbg(1,3)=ijcxst(1)+ipq((jc-1)*nt(isf)+1)*npri(1,1)+(kc-1)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*((nc(is)-1)*ntij)
      icxbg(2,2)=icxbg(2,1)+(kc-jc)*ntij
      icxbg(3,2)=icxbg(2,2)-(ic-kc)*((nc(is)-1)*ntij)
      icxbg(2,3)=icxbg(2,2)+(jc-ic)*nc(is)*ntij
      icxbg(3,3)=icxbg(2,3)-(jc-kc)*((nc(is)-1)*ntij)
      ec(1)=ic.eq.jc
      ec(3)=jc.eq.kc
      ec(2)=ic.eq.kc
      do 12 jesfb=1,2
        jbl(jesfb)=0
        do 10 ist=1,nst
          jbl(jesfb)=jbl(jesfb)+npriri(jesfb,1,ist)*nprirk(2,1,ist)
   10   continue
   12 continue
      jbl1=(ircru*lrcru)*jbl(1)
      jbl2=(ircru*lrcru)*jbl(2)
      jbld=jbl2-jbl1
      if(ec(1)) then
        nwt1(1) = 1
        if(ec(3)) then
          iscmu = 1
          go to 48
        else
          iscmu = 2
          go to 44
        endif
      endif
      if(ec(3)) then
        iscml = 2
        iscmu = 3
        nwt1(2) = 1
        go to 52
      endif
      if(ngw(icc).eq.1) go to 40
      isw2 = 2
      isw3 = 3
      do 24 igw=2,ngw(icc)
        if(icb(lnx,igw,icc).ne.lc) go to 24
        if(icb(jnx,igw,icc).eq.kc) isw2=1
        if(icb(inx,igw,icc).eq.jc) isw3=min(isw3,2)
        if(icb(inx,igw,icc).eq.kc) isw3=1
   24 continue
      if(isw3.eq.1) then
        if(isw2.eq.1) then
          iscmu = 1
          nwt1 (1) = 3
          go to 48
        else
          iscmu = 2
          nwt1(1) = 2
          go to 44
        endif
      endif
      if(isw3.eq.2) then
        iscmu = 2
        nwt1(1) = 1
        nwt1(2) = 2
        go to 48
      endif
      if(isw2.eq.1) then
        iscml = 2
        iscmu = 3
        nwt1(2) = 2
        go to 52
      endif
   40 iscmu = 3
      nwt1(1) = 1
   44 nwt1(2) = 1
   48 iscml = 1
   52 nwt2 = nwt * nfct(icc)
      ngt=nt(isf)*(nt(isf)*(nt(isf)*nt(lsf)))
      do 98 iscm=iscml,iscmu
      fnfct=nwt1(iscm)*nwt2
      if(iscm-2) 56, 60, 64
   56 call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,jc,is,isf,kc,ls,lsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   60 call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,kc,is,isf,jc,ls,lsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   64 call jandk(a,eta,g4,ipt,is,isf,jc,is,isf,kc,is,isf,ic,ls,lsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
   68 ndx=1
      do 96 ircr=1,ircru
      if(ec(iscm)) jrcru=ircr
      do 94 jrcr=1,jrcru
      if(ircr-jrcr) 70, 72, 74
   70 iblb=(ipq(jrcr)+(ircr-1))*jbl2-(jrcr-1)*jbld
      iesfb=2
      icxbg1=icxbg(3,iscm)
      ngti=(nt(isf)*nt(lsf))
      ngtj=(nt(isf)*(nt(isf)*nt(lsf)))
      go to 78
   72 iblb=ipq(ircr)*jbl2+(ircr-1)*jbl1
      iesfb=1
      go to 76
   74 iblb=(ipq(ircr)+(jrcr-1))*jbl2-(ircr-1)*jbld
      iesfb=2
   76 icxbg1=icxbg(iesfb,iscm)
      ngti=(nt(isf)*(nt(isf)*nt(lsf)))
      ngtj=(nt(isf)*nt(lsf))
   78 do 92 krcr=1,ircru
      do 90 lrcr=1,lrcru
      icxs(iesfb,iscm)=icxbg1
      ndxi=ndx
      do 88 it=itl,itu
        if(ircr.eq.jrcr.and.ec(iscm)) jtu = it
        ndxj=ndxi
        do 86 jt=itl,jtu
          kcxs(2,iscm)=kcxbg(iscm)
          ndxl=ndxj
          do 84 kt=itl,itu
            do 82 lt=ltl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, iscm, val )
   80         ndxl=ndxl+1
              kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,1)
   82       continue
   84     continue
          ndxj=ndxj+ngtj
          icxs(iesfb,iscm)=icxs(iesfb,iscm)+npri(iesfb,1)
   86   continue
        ndxi=ndxi+ngti
   88 continue
      iblb=iblb+jbl(iesfb)
      ndx=ndx+ngt
   90 continue
   92 continue
   94 continue
   96 continue
   98 continue
      return
      end
*deck sf3eq3
      subroutine sf3eq3( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      iscmci,
     &         is,       isf,      ks,       ksf,      ijgt,
     &         ijx,      ijy,      ijz,      klgt,     klx,
     &         kly,      klz,      lmnp1,    lmnv,     mcons,
     &         nc,       ncon,     nfct,     ngw,      nrcr,
     &         nt,       ntl,      ntu,      x,        y,
     &         z,        zet )
      implicit real*8 (a-h,o-z)
      logical iceqjc, kceqlc
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1  jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2  npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /cxindx/ ntij,ntik,ntil,ntkl,ntjl,ntjk,ntij1,ntjk1,ntkl1,
     1  ijcxst(2),ikcxst(2),ilcxst(2),klcxst(2),jlcxst(2),jkcxst(2)
      common /dim21/ ipq(256)
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*),cx(*),eta(mrcru,mconu,*),g4(*),icb(4,24,*),h4(*),
     1  ipt(*),ijgt(*),ijx(*),ijy(*),ijz(*),klgt(*),klx(*),kly(*),
     2  klz(*),lmnp1(*),lmnv(3,*),mcons(*),nc(*),ncon(*),nfct(*),ngw(*),
     3  nrcr(*),nt(*),ntl(*),ntu(*),x(mcu,*),y(mcu,*),z(mcu,*),
     4  zet(mconu,*)
      dimension icxbg(3,3), jbl(2,2), jbld(2), kcxbg(3,3)
c
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      iceqjc=ic.eq.jc
      kceqlc=kc.eq.lc
      ngt=(nt(isf)*nt(ksf))**2
      if(iscmci.eq.2) go to 60
c
      do 16 jesfb=1,2
        do 14 lesfb=1,2
          jbl(jesfb,lesfb)=0
          do 12 ist=1,nst
            jbl(jesfb,lesfb)=jbl(jesfb,lesfb)+
     &       npriri(jesfb,1,ist)*nprirk(lesfb,1,ist)
   12     continue
          jbld(jesfb)=jbl(jesfb,2)-jbl(jesfb,1)
   14   continue
   16 continue
      jblii=ipq(krcru)*jbl(1,2)+krcru*jbl(1,1)
      jblij=ipq(krcru)*jbl(2,2)+krcru*jbl(2,1)
      jblid=jblij-jblii
      fnfct=(nwt*nfct(icc))
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*(nc(is)-1)*ntij
      kcxbg(1,1)=klcxst(1)+ipq((kc-1)*nt(ksf)+1)*nprk(1,1)+(lc-1)*ntkl1
      kcxbg(2,1)=klcxst(2)+((kc-1)*nc(ks)+(lc-1))*ntkl
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*(nc(ks)-1)*ntkl
      call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,jc,ks,ksf,kc,ks,ksf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      ndx=1
      do 54 ircr=1,ircru
      if(iceqjc) jrcru=ircr
      do 52 jrcr=1,jrcru
      if(ircr-jrcr) 18, 20, 22
   18 iblkl=(ipq(jrcr)+(ircr-1))*jblij-(jrcr-1)*jblid
      iesfb=2
      icxbg1=icxbg(3,1)
      ngti=(nt(ksf)*nt(ksf))
      ngtj=(nt(isf)*(nt(ksf)*nt(ksf)))
      go to 26
   20 iblkl=ipq(ircr)*jblij+(ircr-1)*jblii
      iesfb=1
      go to 24
   22 iblkl=(ipq(ircr)+(jrcr-1))*jblij-(ircr-1)*jblid
      iesfb=2
   24 icxbg1=icxbg(iesfb,1)
      ngti=(nt(isf)*(nt(ksf)*nt(ksf)))
      ngtj=(nt(ksf)*nt(ksf))
   26 do 50 krcr=1,krcru
      if(kceqlc) lrcru=krcr
      do 48 lrcr=1,lrcru
      if(krcr-lrcr) 28, 30, 32
   28 iblb=iblkl+(ipq(lrcr)+(krcr-1))*jbl(iesfb,2)-(lrcr-1)*jbld(iesfb)
      kesfb=2
      kcxbg1=kcxbg(3,1)
      ngtk=1
      ngtl=nt(ksf)
      go to 36
   30 iblb=iblkl+ipq(krcr)*jbl(iesfb,2)+(krcr-1)*jbl(iesfb,1)
      kesfb=1
      go to 34
   32 iblb=iblkl+(ipq(krcr)+(lrcr-1))*jbl(iesfb,2)-(krcr-1)*jbld(iesfb)
      kesfb=2
   34 kcxbg1=kcxbg(kesfb,1)
      ngtk=nt(ksf)
      ngtl=1
   36 icxs(iesfb,1)=icxbg1
      ndxi=ndx
      do 46 it=itl,itu
        if(iceqjc.and.ircr.eq.jrcr) jtu = it
        ndxj=ndxi
        do 44 jt=itl,jtu
          kcxs(kesfb,1)=kcxbg1
          ndxk=ndxj
          do 42 kt=ktl,ktu
            if(kceqlc.and.krcr.eq.lrcr) ltu = kt
            ndxl=ndxk
            do 40 lt=ktl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 38
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 1, val )
   38         ndxl=ndxl+ngtl
              kcxs(kesfb,1)=kcxs(kesfb,1)+nprk(kesfb,1)
   40       continue
            ndxk=ndxk+ngtk
   42     continue
          ndxj=ndxj+ngtj
        icxs(iesfb,1)=icxs(iesfb,1)+npri(iesfb,1)
   44   continue
        ndxi=ndxi+ngti
   46 continue
      ndx=ndx+ngt
   48 continue
   50 continue
   52 continue
   54 continue
      return
c
   60 icxbg(2,2)=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      icxbg(2,3)=icxbg(2,2)+(lc-kc)*ntik
      kcxbg(2,2)=icxbg(2,3)+(jc-ic)*nc(ks)*ntik
      kcxbg(2,3)=kcxbg(2,2)+(kc-lc)*ntik
      jbl1=0
      jbl2=0
      do 62 ist=1,nst
        jbl1=jbl1+ipq(npriri(2,2,ist)+1)
        jbl2=jbl2+npriri(2,2,ist)**2
   62 continue
      jbldi=jbl2-jbl1
      if(iceqjc.or.kceqlc) then
        fnfct=(nwt*nfct(icc))
      else
        do 64 igw=2,ngw(icc)
          if((icb(inx,igw,icc).eq.ic.and.
     &        icb(jnx,igw,icc).eq.jc).or.
     &       (icb(knx,igw,icc).eq.kc.and.
     &        icb(lnx,igw,icc).eq.lc)) go to 66
   64   continue
        iscmu=3
        fnfct=(nwt*nfct(icc))
        go to 68
   66   fnfct=(2*nwt*nfct(icc))
      endif
      iscmu=2
   68 do 98 iscm=2,iscmu
      if(iscm.eq.2) then
        call jandk(a,eta,g4,ipt,is,isf,ic,ks,ksf,kc,is,isf,jc,ks,ksf,lc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      else
        call jandk(a,eta,g4,ipt,is,isf,ic,ks,ksf,lc,is,isf,jc,ks,ksf,kc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      endif
      ndx=1
      ikrcr=0
      do 96 ircr=1,ircru
      if((iceqjc.and.kceqlc)) jrcru=ircr
      do 94 krcr=1,krcru
      ikrcr=ikrcr+1
      jlrcr=0
      do 92 jrcr=1,jrcru
      if(iceqjc.and.kceqlc) then
        if(ircr.eq.jrcr) then
          lrcru = krcr
        else
          lrcru = krcru
        endif
      endif
      do 90 lrcr=1,lrcru
      jlrcr=jlrcr+1
      if(ikrcr-jlrcr) 70, 72, 74
   70 iblb=(ipq(jlrcr)+(ikrcr-1))*jbl2-(jlrcr-1)*jbldi
      icxs(2,iscm)=kcxbg(2,iscm)
      kcxbg1=icxbg(2,iscm)
      ngtk=1
      ngtl=(nt(isf)*nt(ksf))
      go to 78
   72 iblb=(ipq(ikrcr+1)-1)*jbl2-(ikrcr-1)*jbldi
      go to 76
   74 iblb=(ipq(ikrcr)+(jlrcr-1))*jbl2-(ikrcr-1)*jbldi
   76 icxs(2,iscm)=icxbg(2,iscm)
      kcxbg1=kcxbg(2,iscm)
      ngtk=(nt(isf)*nt(ksf))
      ngtl=1
   78 ndxk=ndx
      do 88 it=itl,itu
        if(((iceqjc.and.kceqlc).and.ikrcr.eq.jlrcr)) jtu=it
        do 86 kt=ktl,ktu
          kcxs(2,iscm)=kcxbg1
          ndxl=ndxk
          do 84 jt=itl,jtu
            if((iceqjc.and.kceqlc).and.ikrcr.eq.jlrcr) then
              if(it.eq.jt) then
                ltu = kt
              else
                ltu = ktu
              endif
            endif
            do 82 lt=ktl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              if(ikrcr.eq.jlrcr) then
                call aoso2e( h4, cx, iscm, val )
              else
                call aoso2( h4, cx, iscm, val )
              endif
   80         ndxl=ndxl+ngtl
              kcxs(2,iscm)=kcxs(2,iscm)+npri(2,2)
   82       continue
   84     continue
          ndxk=ndxk+ngtk
          icxs(2,iscm)=icxs(2,iscm)+npri(2,2)
   86   continue
   88 continue
      ndx=ndx+ngt
   90 continue
   92 continue
   94 continue
   96 continue
   98 continue
      return
      end
