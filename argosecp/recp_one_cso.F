      subroutine cints(a,ccr,etai,etaj,g2,lcr,lambu, 
     &           ltot1,mcrs,mproju,nc,ncr,nkcrl,nkcru,x,y,z,zcr,zet )
!  lcru is the max l value + 1 for the potential.
!  ncr contains the value of n for each term.
!  zcr contains the value of alpha for each term.
!  ccr contains the coefficient of each term.
      USE RECP_IPT 
      USE RECP_OUTPUT 
      implicit real*8 (a-h,o-z)
#include "inc_print.h"
      logical esf, esfc, igueq1, jgueq1
      parameter (a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
      common /parmr/ tol
c
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),ccr(*),etai(mrcru,*),
     1  etaj(mrcru,*),g2(*),lcr(*),
     2  mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*),x(mcu,*),
     3  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
      integer i 

      CALL OUTPUT_LOC('cint','E')
      do i=1,3
        IPT_CZ_CRDA(1,i)=a1
        IPT_CZ_CRDB(1,i)=a1
      enddo
      if(esfc) iu=(lit+1)*((lit+1)-1)/2
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
      if(.not.igueq1) then
        do i = 1, nc1
           IPT_CZ_G1(i) = 0.0d0
        enddo
      endif
      ai=zet(ig,icons)
      tai=ai+ai
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) then
        do i = 1, ij
           IPT_CZ_GOUT(i) = 0.0d0
        enddo
      endif
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
            IPT_CZ_CRDA(2,1)=xka
            IPT_CZ_CRDA(2,2)=yka
            IPT_CZ_CRDA(2,3)=zka
            if(lit.ge.3) then
              do i=1,3
                do j=3,lit
                  IPT_CZ_CRDA(j,i)=IPT_CZ_CRDA(2,i)*IPT_CZ_CRDA(j-1,i)
                enddo
              enddo
            endif
          endif
          xkb=xc-xj
          ykb=yc-yj
          zkb=zc-zj
          cb=sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
          if(ljt.ge.2) then
            IPT_CZ_CRDB(2,1)=xkb
            IPT_CZ_CRDB(2,2)=ykb
            IPT_CZ_CRDB(2,3)=zkb
            if(ljt.ge.3) then
              do 240 i=1,3
                do 230 j=3,ljt
                  IPT_CZ_CRDB(j,i)=IPT_CZ_CRDB(2,i)*IPT_CZ_CRDB(j-1,i)
  230           continue
  240         continue
            endif
          endif

!         type1 pseudo
          if(aarr1.le.tol) then
            xk=xp-xc
            yk=yp-yc
            zk=zp-zc
            call pseud1(a,ccr,IPT_CZ_CRDA,IPT_CZ_CRDB,IPT_CZ_GOUT,
     &           IPT_LMNV,ltot1,ncr,nkcrl,nkcru,zcr)
          endif

!         type2 pseudo
          aarr2=apr*(ca-cb)**2
          if(aarr2.le.tol.and.lcru.ne.0) then
            call pseud2( a, ccr, IPT_CZ_CRDA, IPT_CZ_CRDB,
     &        IPT_CZ_GOUT,lambu,ltot1,mproju,ncr,nkcrl,nkcru, zcr )
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
          IPT_CZ_G1(j1+i)=IPT_CZ_G1(j1+i)+IPT_CZ_GOUT(i)*etaj(jrcr,jg)
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
                g2(ij2+j)=g2(ij2+j)+IPT_CZ_G1(j1+j)*etai(ircr,ig)+
     1                              IPT_CZ_G1(i1n)*etai(jrcr,ig)
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
              g2(ij2+i)=g2(ij2+i)+IPT_CZ_G1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      CALL OUTPUT_LOC('cint','X')
      return
      end


      subroutine lsints( a,ccr,etai,etaj,g2,lls,lambu,ltot1,mcrs,
     &                   mproju,nc,ncr,nklsl,nklsu,x,y,z,zcr,zet )
c
c  lcru is the max l value + 1 for the potential.
c  ncr contains the value of n for each term.
c  zcr contains the value of alpha for each term.
c  ccr contains the coefficient of each term.

      USE RECP_IPT
      USE RECP_OUTPUT
      implicit real*8 (a-h,o-z)
#include "inc_print.h"
      logical esf, esfc, igueq1, jgueq1
      parameter (a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
      common /parmr/ tol
c
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
c
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),ccr(*),etai(mrcru,*),
     1  etaj(mrcru,*),g2(*),lls(*),
     2  mcrs(*),nc(*),ncr(*),nklsl(4,*),nklsu(4,*),x(mcu,*),
     3  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)
      integer i
      CALL OUTPUT_LOC('LSINTS','E') 
      do 120 i=1,3
        IPT_CZ_CRDA(1,i)=a1
  120   IPT_CZ_CRDB(1,i)=a1
      if(esfc) iu=(lit+1)*((lit+1)-1)/2
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
      if(.not.igueq1) then
        do i = 1, nc1
           IPT_CZ_G1(i) = 0.0d0
        enddo
      endif
      ai=zet(ig,icons)
      tai=ai+ai
c
c  j primitive
c
      if(esfc) jgu=ig
      do 390 jg=1,jgu
      if(.not.jgueq1) then
        do i = 1, ij
           IPT_CZ_GOUT(i) = 0.0d0 
        enddo
      endif
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
        IPT_CZ_CRDA(2,1)=xka
        IPT_CZ_CRDA(2,2)=yka
        IPT_CZ_CRDA(2,3)=zka
        if(lit.ge.3) then
          do i=1,3
            do j=3,lit
              IPT_CZ_CRDA(j,i)=IPT_CZ_CRDA(2,i)*IPT_CZ_CRDA(j-1,i)
            enddo
          enddo
        endif
      endif
      xkb=xc-xj
      ykb=yc-yj
      zkb=zc-zj
      cb=sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
      if(ljt.ge.2) then
        IPT_CZ_CRDB(2,1)=xkb
        IPT_CZ_CRDB(2,2)=ykb
        IPT_CZ_CRDB(2,3)=zkb
        if(ljt.ge.3) then
          do 240 i=1,3
            do 230 j=3,ljt
              IPT_CZ_CRDB(j,i)=IPT_CZ_CRDB(2,i)*IPT_CZ_CRDB(j-1,i)
  230       continue
  240     continue
        endif
      endif
      aarr2=apr*(ca-cb)**2
      if(aarr2.le.tol) then
        call pseud3( a, ccr, IPT_CZ_CRDA, IPT_CZ_CRDB,
     &    IPT_CZ_GOUT,lambu,ltot1,
     &    mproju, ncr, nklsl, nklsu, zcr )
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
          IPT_CZ_G1(j1+i)=IPT_CZ_G1(j1+i)+IPT_CZ_GOUT(i)*etaj(jrcr,jg)
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
                  g2(ij2+k)=g2(ij2+k)+IPT_CZ_G1(j1+k)*etai(ircr,ig)-
     &                                IPT_CZ_G1(i1n+k)*etai(jrcr,ig)
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
              g2(ij2+i)=g2(ij2+i)+IPT_CZ_G1(j1+i)*etai(ircr,ig)
  450       continue
            ij2=ij2+ij
            j1=j1+ij
  460     continue
        endif
  470 continue
  490 continue
      CALL OUTPUT_LOC('LSINTS','X') 
      return
      end
 
