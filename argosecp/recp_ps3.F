      subroutine pseud3( a,ccr,crda,crdb,gout,
     &  ipt,lambu,ltot1,mproju,ncr,nklsl,nklsu,zcr )
!  computes spin-orbit potential integrals

      USE RECP_IPT
      USE RECP_OUTPUT
      implicit real*8 (a-h,o-z)
#include "inc_print.h"
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0)
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
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),ccr(*),
     &  crda(lit,3),crdb(ljt,3),gout(*),ipt(*),
     &  ncr(*),nklsl(4,*),nklsu(4,*),zcr(*)
      CALL OUTPUT_LOC('RECP_PSEUD3','E')
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
            call rad2( a,    ccr,        ipt,        kcrl,
     &                 kcru, l,          lambu,      lmahi,
     &                 lmalo,lmbhi,      lmblo,      ltot1,
     &                 ncr,  rka,        rkb,
     &                 zcr )
c
c  compute angular integrals and combine with radial integrals
c
      ijt=0
      do 84 it=itl,itu
      call ang2(IPT_P23_ANGA,crda,it,l,lit,lmalo,lmahi,
     &  IPT_LMNV,mproju,xka,yka,zka)
      do 80 jt=jtl,jtu
      call ang2(IPT_P23_ANGB,crdb,jt,l,ljt,lmblo,lmbhi,
     &  IPT_LMNV,mproju,xkb,ykb,zkb)
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
     1        IPT_P23_ANGA(nlma,IPT_SH_MR(kt,m),lama)
     &          *IPT_P23_ANGB((n+1)-nlma,IPT_SH_MC(kt,m),lamb)
     2       -IPT_P23_ANGA(nlma,IPT_SH_MC(kt,m),lama)
     &          *IPT_P23_ANGB((n+1)-nlma,IPT_SH_MR(kt,m),lamb)
   56     continue
          s=s+angp*IPT_SH_FLMTX(kt,(l-2)**2+m)
   60   continue
        gout(ijt+kt)=gout(ijt+kt)+s*IPT_P23_QSUM(n,lamb,lama)
   64 continue
   68 continue
   72 continue
   76 continue
      ijt=ijt+3
   80 continue
   84 continue
   88 continue
      CALL OUTPUT_LOC('RECP_PSEUD3','X')
      return
      end
