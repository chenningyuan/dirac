      subroutine recp_stvcz(a,ccr,chg,eta,
     &           ilxyz,lcr,lls,lmnv,lproju,mcrs,nc,ncr,
     &           nkcrl,nkcru,nklsl,nklsu,nopir,nprir,nt,x,y,z,zcr,zet)
 
      USE RECP_IPT
      USE RECP_OUTPUT
      implicit real*8 (a-h,o-z)
#include "inc_print.h"
      logical chsign,ec,esf,esfb,esfc,esfbc,igueq1,jgueq1
      parameter (a1s2=0.5d0)
 
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
 
      common /parmr/ tol
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
 
      integer
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir(3),      inam(5),     nnam, mdum(32)
 
      common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
     1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ic,
     2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
     3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension a(*),ccr(*),chg(*),eta(mrcru,mconu,*),
     1  ilxyz(3,*),lcr(*),
     2  lls(*),lmnv(3,*),mcrs(*),nc(*),ncr(*),nkcrl(6,*),nkcru(6,*),
     3  nklsl(4,*),nklsu(4,*),nopir(*),nprir(2,mstu,*),nt(*),x(mcu,*),
     4  y(mcu,*),z(mcu,*),zcr(*),zet(mconu,*)

      integer ic0,ic1,ltot1,mproju,lproju,lamau,lambu,i,IPQ_TEMP 

      CALL OUTPUT_LOC('ARGOS_STVCZ','E')

      lijt = lit*ljt
      ec   = ic.eq.jc
      esfc = esf.and.ec
      icx2 = icxi2+((ic-1)*nc(js)+jc-1)*ntij2
      if(.not.esfc.and.ircru.ne.1.and.esf) then
        ibld=ibl2-ibl1
        icx3=icx2-(ic-jc)*(nc(is)-1)*ntij2
      endif
      if(esf) then
        IPQ_TEMP = ((ic-1)*nt(isf)+1)*(((ic-1)*nt(isf)+1)-1)/2
        icx1=icxi1+IPQ_TEMP*nblt1+(jc-1)*ntij1
      endif
!
!  ishell
!
      xi=x(ic,is)
      yi=y(ic,is)
      zi=z(ic,is)
!
!  jshell
!
      xj=x(jc,js)
      yj=y(jc,js)
      zj=z(jc,js)
      xij=a1s2*(xi+xj)
      dxij=xi-xj
      xijm=a1s2*dxij
      yij=a1s2*(yi+yj)
      dyij=yi-yj
      yijm=a1s2*dyij
      zij=a1s2*(zi+zj)
      dzij=zi-zj
      zijm=a1s2*dzij
      rr=dxij*dxij+dyij*dyij+dzij*dzij
!
!  prepare items for pairs of (i,j) functions
!
      ij=0
      do it=itl,itu
!       print *,'stvcz-it',it
        nx=ljt*lmnv(1,it)+1
        ny=ljt*lmnv(2,it)+1
        nz=ljt*lmnv(3,it)+1
        do jt=jtl,jtu
          ij=ij+1
        enddo
      enddo
      ij=nop*ij
      nc1=jrcru*ij
      if(esfc) then
        IPQ_TEMP = (ircru+1)*((ircru+1)-1)/2
        nc2=IPQ_TEMP*ij
      else
        nc2=ircru*nc1
      endif
      igueq1=igu.eq.1
      jgueq1=jgu.eq.1
      if(igueq1) then
        ic1=0
      else
        ic1=nc2
      endif
      if(jgueq1) then
        ic0=0
      else
        ic0=nc1
      endif

!     Allocate variables
      CALL RECP_IPTA_STVCZ(ic0,ic1,ij,lit,ljt,ltot1,mproju,lproju,
     &     lamau,lambu)
!!     do i = 1, nc2
!     do i = 1, ic1
!        IPT_CZ_G2(i) = 0.0d0
!     enddo

!     ------------------
!     Calculate integral 
!     ------------------
      if     (nnam.eq.1) then
!       PRINT *, ' *skip Overlap(S) integral'
      elseif (nnam.eq.2) then
!       PRINT *, ' *skip Kinetic(T) integral'
      elseif (nnam.eq.3) then
!       PRINT *, ' *skip Coulomb(V) integral'
      elseif (nnam.eq.4) then
            call cints( a,ccr,eta(1,1,icons),eta(1,1,jcons),IPT_CZ_G2,
     &                  lcr,lambu,ltot1,mcrs,mproju,nc,ncr,nkcrl,
     &                  nkcru,x,y,z,zcr,zet )
      elseif (nnam.eq.5) then
            call lsints( a,ccr,eta(1,1,icons),eta(1,1,jcons),IPT_CZ_G2, 
     &                   lls,lambu,ltot1,mcrs,mproju,nc,ncr,nklsl,
     &                   nklsu,x,y,z,zcr,zet )
      endif

!     -------------------------------

      ndx=0
      ndxa=0
      do 92 ircr=1,ircru
      if(esfc) jrcru=ircr
      do 90 jrcr=1,jrcru
      esfb=esf.and.ircr.eq.jrcr
      esfbc=esfb.and.ec
      chsign=nnam.ge.5.and.esf.and.ircr.lt.jrcr
      if(.not.esfb) go to 60
      icx=icx1
      go to 66
   60 if((.not.esf).or.esfc) go to 64
      if(ircr.gt.jrcr) go to 62
      IPQ_TEMP = (jrcr)*((jrcr)-1)/2
      ndxa=(IPQ_TEMP+(ircr-1))*ibl2-(jrcr-1)*ibld
      icx=icx3
      ngti=nop
      ngtj=(nop*nt(isf))
      go to 70
   62 if(jrcr.eq.1) then
        IPQ_TEMP = (ircr)*((ircr)-1)/2
        ndxa=(IPQ_TEMP+(jrcr-1))*ibl2-(ircr-1)*ibld
      endif
   64 icx=icx2
   66 ngti=(nop*nt(jsf))
      ngtj=nop
   70 ndxi=ndx

      do 88 it=itl,itu
      ndxj=ndxi
      if(esfbc) jtu=it
      ndxi=ndxi+ngti
      do 86 jt=jtl,jtu
        index=ndxj
        ndxj=ndxj+ngtj
        ndxb=ndxa
        do 84 ist=1,nst
          iopu=nopir(ist)
          if(esfb) then
            npr=nprir(1,ist,ijsf)
          else
            npr=nprir(2,ist,ijsf)
          endif
c          if(iopu.eq.0) go to 82
c          if(npr.eq.0) go to 84
          do 80 iop=1,iopu
            IF (ic1.NE.0) THEN
               IF ((index+ilxyz(iop,ist)).GT.ic1) THEN
                  print *,'IPT_CZ_G2 index exceeded maximum'
                  print *,(index+ilxyz(iop,ist)),'>',ic1,'(ic1)'
                  stop
               ENDIF
               val = IPT_CZ_G2(index+ilxyz(iop,ist))
            ELSEIF (ic0.NE.0) THEN
               IF ((index+ilxyz(iop,ist)).GT.ic0) THEN
                  print *,'IPT_CZ_G1 index exceeded maximum'
                  print *,(index+ilxyz(iop,ist)),'>',ic0,'(ic0)'
                  stop
               ENDIF
               val = IPT_CZ_G1(index+ilxyz(iop,ist))
!              print *,'g2-1',val
            ELSE
               IF ((index+ilxyz(iop,ist)).GT.ij) THEN
                  print *,'IPT_CZ_GOUT index exceeded maximum'
                  print *,(index+ilxyz(iop,ist)),'>',ij,'(ij)'
                  stop
               ENDIF
               val = IPT_CZ_GOUT(index+ilxyz(iop,ist))
!              print *,'g2-2',val
            ENDIF
            if(chsign) val=-val
c      print *, 'fnfct = ',fnfct
            val=fnfct*val

            do 76 ibl=1,npr
              IPT_ONE_H2(ndxb+ibl)=IPT_ONE_H2(ndxb+ibl)
     &                             +IPT_CX(icx+ibl)*val
c      print *, 'h2(',ndxb,'+',ibl,')=',h2(ndxb+ibl)
   76       continue
   78       ndxb=ndxb+npr
   80     continue
   82     icx=icx+npr
   84   continue
   86 continue
   88 continue
      ndxa=ndxb
      ndx=ndx+ij
   90 continue
   92 continue

!     Deallocate variables
      CALL RECP_IPTD_STVCZ
      CALL OUTPUT_LOC('ARGOS_STVCZ','X')

      end

