      subroutine pseud1(a,ccr,crda,crdb,gout,lmnv,ltot1,ncr,
     &  nkcrl,nkcru,zcr)
!  compute type 1 core potential integrals

      USE RECP_FUNCTION
      USE RECP_IPT
      USE RECP_OUTPUT
      USE RECP_RAD1
      implicit real*8 (a-h,o-z)
#include "inc_print.h"
      logical esf, esfc, igueq1, jgueq1
      parameter (a0=0.0d0)
      common /parmr/ tol
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /callin/ xka,yka,zka,ca,xkb,ykb,zkb,cb,tai,taj,aa,taa,
     1  aarr1,aarr2,xk,yk,zk,fctr2,kcrs,lcru
      dimension a(*),ccr(*),crda(lit,3),crdb(ljt,3),
     &  gout(*),lmnv(3,*),ncr(*),nkcrl(6,*),nkcru(6,*),zcr(*)
      integer i,j
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
      do i = 1, ltot1
!        do j = 1, lamu
         do j = 1, ltot1
            IPT_P1_QSUM(I,J) = 0.0d0
         enddo
      enddo

      kcrl=nkcrl(1,kcrs)
      kcru=nkcru(1,kcrs)
!     do 40 kcr=kcrl,kcru
!       alpha=aa+zcr(kcr)
!       # exponential factor from q functions included in dum
!       dum=aarr1+zcr(kcr)*arp2/alpha
!       if(dum.gt.tol) go to 40
!       prd=fctr2*ccr(kcr)*exp(-dum)
!       t=alpt/alpha
!       call recur1(alpha,a(ipt(11)),ncr(kcr),ltot1,q,t,rk)
!       do 30 lam=1,lamu
!         nhi=ltot1-mod(ltot1-lam,2)
!         do 20 n=lam,nhi,2
!           IPT_P1_QSUM(n,lam)=IPT_P1_QSUM(n,lam)+prd*IPT_P1_Q(n,lam)
!  20     continue
!  30   continue
!  40 continue

      call RECP_RAD1_MAIN
     &    (aa,aarr1,alpt,arp2,ccr,fctr2,kcrl,kcru,
     &     lamu,ltot1,ncr,rk,tol,zcr)

      ijt=0
      do 90 it=itl,itu
!       print *,'pseud1-it',it 
        na1=lmnv(1,it)+1
        la1=lmnv(2,it)+1
        ma1=lmnv(3,it)+1
        do 80 jt=jtl,jtu
          ijt=ijt+1
          s=a0
!         print *,'pseud1-jt',jt 
          nb1=lmnv(1,jt)+1
          lb1=lmnv(2,jt)+1
          mb1=lmnv(3,jt)+1
!
!  compute angular integrals
!
          call RECP_FUNCTION_FACAB
     &         (IPT_SH_BINOM,na1,nb1,crda(1,1),crdb(1,1),IPT_P1_XAB)
          call RECP_FUNCTION_FACAB
     &         (IPT_SH_BINOM,la1,lb1,crda(1,2),crdb(1,2),IPT_P1_YAB)
          call RECP_FUNCTION_FACAB
     &         (IPT_SH_BINOM,ma1,mb1,crda(1,3),crdb(1,3),IPT_P1_ZAB)
          call ang1(na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,
     &      ltot1,xk,yk,zk)
!
!  combine angular and radial integrals
!
          do 70 lam=1,lamu
            nhi=ltot1-mod(ltot1-lam,2)
            do 60 n=lam,nhi,2
              s=s+IPT_P1_ANG(n,lam)*IPT_P1_QSUM(n,lam)
   60       continue
   70     continue
          gout(ijt)=gout(ijt)+s
   80   continue
   90 continue
      return
      end


