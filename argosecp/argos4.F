c
c  RCS $Revision$  $Date$
c
*deck sf3eq4
      subroutine sf3eq4( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      iscmci,
     &         is,       isf,      ks,       ksf,      ls,
     &         lsf,      ijgt,     ijx,      ijy,      ijz,
     &         klgt,     klx,      kly,      klz,      lmnp1,
     &         lmnv,     mcons,    nc,       ncon,     nfct,
     &         ngw,      nrcr,     nt,       ntl,      ntu,
     &         x,        y,        z,        zet )
      implicit real*8 (a-h,o-z)
      logical iceqjc
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
      dimension icxbg(3,3), jbl(2), kcxbg(3)
c
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      iceqjc=ic.eq.jc
      do 12 jesfb=iscmci,2
        jbl(jesfb)=0
        do 10 ist=1,nst
          jbl(jesfb)=jbl(jesfb)+
     &      npriri(jesfb,iscmci,ist)*nprirk(2,iscmci,ist)
   10   continue
   12 continue
      if(iscmci.eq.2) go to 60
c
      jbl1=(krcru*lrcru)*jbl(1)
      jbl2=(krcru*lrcru)*jbl(2)
      jbld=jbl2-jbl1
      fnfct=(nwt*nfct(icc))
      icxbg(1,1)=ijcxst(1)+ipq((ic-1)*nt(isf)+1)*npri(1,1)+(jc-1)*ntij1
      icxbg(2,1)=ijcxst(2)+((ic-1)*nc(is)+(jc-1))*ntij
      icxbg(3,1)=icxbg(2,1)-(ic-jc)*(nc(is)-1)*ntij
      kcxbg(1)=klcxst(2)+((kc-1)*nc(ls)+(lc-1))*ntkl
      call jandk(a,eta,g4,ipt,is,isf,ic,is,isf,jc,ks,ksf,kc,ls,lsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      ndx=1
      ngt=nt(isf)*(nt(isf)*(nt(ksf)*nt(lsf)))
      do 48 ircr=1,ircru
      if(iceqjc) jrcru=ircr
      do 46 jrcr=1,jrcru
      if(ircr-jrcr) 18, 20, 22
   18 iblb=(ipq(jrcr)+(ircr-1))*jbl2-(jrcr-1)*jbld
      iesfb=2
      icxbg1=icxbg(3,1)
      ngti=(nt(ksf)*nt(lsf))
      ngtj=(nt(isf)*(nt(ksf)*nt(lsf)))
      go to 26
   20 iblb=ipq(ircr)*jbl2+(ircr-1)*jbl1
      iesfb=1
      go to 24
   22 iblb=(ipq(ircr)+(jrcr-1))*jbl2-(ircr-1)*jbld
      iesfb=2
   24 icxbg1=icxbg(iesfb,1)
      ngti=(nt(isf)*(nt(ksf)*nt(lsf)))
      ngtj=(nt(ksf)*nt(lsf))
   26 do 44 krcr=1,krcru
      do 42 lrcr=1,lrcru
      icxs(iesfb,1)=icxbg1
      ndxi=ndx
      do 40 it=itl,itu
        if(iceqjc.and.ircr.eq.jrcr) jtu=it
        ndxj=ndxi
        do 38 jt=itl,jtu
          kcxs(2,1)=kcxbg(1)
          ndxl=ndxj
          do 36 kt=ktl,ktu
            do 34 lt=ltl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 30
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 1, val )
   30         ndxl=ndxl+1
              kcxs(2,1)=kcxs(2,1)+nprk(2,1)
   34       continue
   36     continue
          ndxj=ndxj+ngtj
          icxs(iesfb,1)=icxs(iesfb,1)+npri(iesfb,1)
   38   continue
        ndxi=ndxi+ngti
   40 continue
      iblb=iblb+jbl(iesfb)
      ndx=ndx+ngt
   42 continue
   44 continue
   46 continue
   48 continue
      return
c
   60 icxbg(2,2)=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      icxbg(2,3)=icxbg(2,2)+(jc-ic)*nc(ks)*ntik
      kcxbg(3)=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg(2)=kcxbg(3)+(jc-ic)*nc(ls)*ntil
      if(iceqjc) then
        fnfct=(nwt*nfct(icc))
      else
        do 64 igw=2,ngw(icc)
          if(icb(knx,igw,icc).eq.kc.and.
     &       icb(lnx,igw,icc).eq.lc.and.
     &       icb(inx,igw,icc).eq.jc) go to 66
   64   continue
        iscmu=3
        fnfct=(nwt*nfct(icc))
        go to 68
   66   fnfct=(2*nwt*nfct(icc))
      endif
      iscmu=2
   68 do 98 iscm=2,iscmu
      if(iscm.eq.2) then
        call jandk(a,eta,g4,ipt,is,isf,ic,ks,ksf,kc,is,isf,jc,ls,lsf,lc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      else
        call jandk(a,eta,g4,ipt,is,isf,jc,ks,ksf,kc,is,isf,ic,ls,lsf,lc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      endif
      ndx=0
      iblb=0
      do 96 ircr=1,ircru
      do 94 krcr=1,krcru
      do 92 jrcr=1,ircru
      do 90 lrcr=1,lrcru
      icxs(2,iscm)=icxbg(2,iscm)
      do 88 it=itl,itu
        do 86 kt=ktl,ktu
          kcxs(2,iscm)=kcxbg(iscm)
          do 84 jt=itl,itu
            do 82 lt=ltl,ltu
              ndx=ndx+1
              val=g4(ndx)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, iscm, val )
   80         kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,2)
   82       continue
   84     continue
          icxs(2,iscm)=icxs(2,iscm)+npri(2,2)
   86   continue
   88 continue
      iblb=iblb+jbl(2)
   90 continue
   92 continue
   94 continue
   96 continue
   98 continue
      return
      end
*deck sf3eq5
      subroutine sf3eq5( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      is,
     &         isf,      js,       jsf,      ijgt,     ijx,
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
      dimension ec(3), icxbg(3), jbl(2), kcxbg(3,3)
c
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      icxbg(1)=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      icxbg(2)=icxbg(1)+(kc-jc)*ntij
      icxbg(3)=icxbg(2)+(lc-kc)*ntij
      kcxbg(1,1)=jkcxst(1)+ipq((kc-1)*nt(jsf)+1)*nprk(1,1)+(lc-1)*ntjk1
      kcxbg(1,3)=jkcxst(1)+ipq((jc-1)*nt(jsf)+1)*nprk(1,1)+(kc-1)*ntjk1
      kcxbg(1,2)=kcxbg(1,3)+(lc-kc)*ntjk1
      kcxbg(2,1)=jkcxst(2)+((kc-1)*nc(js)+(lc-1))*ntjk
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*((nc(js)-1)*ntjk)
      kcxbg(2,2)=kcxbg(2,1)+(jc-kc)*nc(js)*ntjk
      kcxbg(3,2)=kcxbg(2,2)-(jc-lc)*((nc(js)-1)*ntjk)
      kcxbg(2,3)=kcxbg(2,2)+(kc-lc)*ntjk
      kcxbg(3,3)=kcxbg(2,3)-(jc-kc)*((nc(js)-1)*ntjk)
      ec(1)=kc.eq.lc
      ec(2)=jc.eq.lc
      ec(3)=jc.eq.kc
      do 12 jesfb=1,2
        jbl(jesfb)=0
        do 10 ist=1,nst
          jbl(jesfb)=jbl(jesfb)+npriri(2,1,ist)*nprirk(jesfb,1,ist)
   10   continue
   12 continue
      jbld=jbl(2)-jbl(1)
      jblkl=ipq(jrcru)*jbl(2)+jrcru*jbl(1)
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
        if(icb(inx,igw,icc).ne.ic) go to 24
        if (icb(knx,igw,icc) .eq. jc) isw2 = 1
        if (icb(lnx,igw,icc) .eq. kc) isw3 = min(isw3,2)
        if (icb(lnx,igw,icc) .eq. jc) isw3 = 1
   24 continue
      if(isw3.eq.1) then
        if(isw2.eq.1) then
          iscmu = 1
          nwt1(1) = 3
          go to 48
        else
          iscmu = 2
          nwt1(1) = 2
          go to 44
        endif
      elseif(isw3.eq.2) then
        iscmu = 2
        nwt1(1) = 1
        nwt1(2) = 2
        go to 48
      elseif(isw2.eq.1) then
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
      ngtj=nt(jsf)**2
      do 98 iscm=iscml,iscmu
      fnfct=nwt1(iscm)*nwt2
      if(iscm-2) 56, 60, 64
   56 call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,jc,js,jsf,kc,js,jsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   60 call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,kc,js,jsf,jc,js,jsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      go to 68
   64 call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,lc,js,jsf,jc,js,jsf,kc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
   68 ndxj=1
      iblbij=0
      do 96 ircr=1,ircru
      do 94 jrcr=1,jrcru
      do 92 krcr=1,jrcru
      if(ec(iscm)) lrcru=krcr
      do 90 lrcr=1,lrcru
      if(krcr-lrcr) 70, 72, 74
   70 iblb=iblbij+(ipq(lrcr)+(krcr-1))*jbl(2)-(lrcr-1)*jbld
      kesfb=2
      kcxbg1=kcxbg(3,iscm)
      ngtk=1
      ngtl=nt(jsf)
      go to 78
   72 iblb=iblbij+ipq(krcr)*jbl(2)+(krcr-1)*jbl(1)
      kesfb=1
      go to 76
   74 iblb=iblbij+(ipq(krcr)+(lrcr-1))*jbl(2)-(krcr-1)*jbld
      kesfb=2
   76 kcxbg1=kcxbg(kesfb,iscm)
      ngtk=nt(jsf)
      ngtl=1
   78 icxs(2,iscm)=icxbg(iscm)
      do 88 it=itl,itu
        do 86 jt=jtl,jtu
          kcxs(kesfb,iscm)=kcxbg1
          ndxk=ndxj
          do 84 kt=jtl,jtu
            if(krcr.eq.lrcr.and.ec(iscm)) ltu=kt
            ndxl=ndxk
            do 82 lt=jtl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, iscm, val )
   80         ndxl=ndxl+ngtl
              kcxs(kesfb,iscm)=kcxs(kesfb,iscm)+nprk(kesfb,1)
   82       continue
            ndxk=ndxk+ngtk
   84     continue
          ndxj=ndxj+ngtj
          icxs(2,iscm)=icxs(2,iscm)+npri(2,1)
   86   continue
   88 continue
   90 continue
   92 continue
      iblbij=iblbij+jblkl
   94 continue
   96 continue
   98 continue
      return
      end
*deck sf3eq6
      subroutine sf3eq6( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      iscmci,
     &         is,       isf,      js,       jsf,      ls,
     &         lsf,      ijgt,     ijx,      ijy,      ijz,
     &         klgt,     klx,      kly,      klz,      lmnp1,
     &         lmnv,     mcons,    nc,       ncon,     nfct,
     &         ngw,      nrcr,     nt,       ntl,      ntu,
     &         x,        y,        z,        zet )
      implicit real*8 (a-h,o-z)
      logical jceqkc
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
      dimension icxbg(3),jbl(2), kcxbg(3,3)
c
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      jceqkc=jc.eq.kc
      jesfbl=(5-iscmci)/2
      do 12 jesfb=jesfbl,2
        jbl(jesfb)=0
        do 10 ist=1,nst
          jbl(jesfb)=jbl(jesfb)+
     &      npriri(2,iscmci,ist)*nprirk(jesfb,iscmci,ist)
   10   continue
   12 continue
      if(iscmci.eq.3) go to 60
c
      icxbg(1)=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      icxbg(2)=icxbg(1)+(kc-jc)*ntij
      kcxbg(2,1)=klcxst(2)+((kc-1)*nc(ls)+(lc-1))*ntkl
      kcxbg(2,2)=kcxbg(2,1)+(jc-kc)*nc(ls)*ntkl
      if(jceqkc) then
        fnfct=(nwt*nfct(icc))
      else
        do 16 igw=2,ngw(icc)
          if(icb(inx,igw,icc).eq.ic.and.
     &       icb(jnx,igw,icc).eq.kc.and.
     &       icb(lnx,igw,icc).eq.lc) go to 24
   16   continue
        iscml=1
        fnfct=(nwt*nfct(icc))
        go to 28
   24   fnfct=(2*nwt*nfct(icc))
      endif
      iscml=2
   28 do 48 iscm=iscml,2
      if(iscm.eq.1) then
        call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,jc,js,jsf,kc,ls,lsf,lc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      else
        call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,kc,js,jsf,jc,ls,lsf,lc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      endif
      ndx=0
      iblb=0
      do 46 ircr=1,ircru
      do 44 jrcr=1,jrcru
      do 42 krcr=1,krcru
      do 40 lrcr=1,lrcru
      icxs(2,iscm)=icxbg(iscm)
      do 38 it=itl,itu
        do 36 jt=jtl,jtu
          kcxs(2,iscm)=kcxbg(2,iscm)
          do 34 kt=ktl,ktu
            do 32 lt=ltl,ltu
              ndx=ndx+1
              val=g4(ndx)
              if(abs(val).le.cutoff) go to 30
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, iscm, val )
   30         kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,1)
   32       continue
   34     continue
          icxs(2,iscm)=icxs(2,iscm)+npri(2,1)
   36   continue
   38 continue
      iblb=iblb+jbl(2)
   40 continue
   42 continue
   44 continue
   46 continue
   48 continue
      return
c
   60 fnfct=(nwt*nfct(icc))
      jbld=jbl(2)-jbl(1)
      jbljk=ipq(jrcru)*jbl(2)+jrcru*jbl(1)
      icxbg(3)=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg(1,3)=jkcxst(1)+ipq((jc-1)*nt(jsf)+1)*nprk(1,3)+(kc-1)*ntjk1
      kcxbg(2,3)=jkcxst(2)+((jc-1)*nc(js)+(kc-1))*ntjk
      kcxbg(3,3)=kcxbg(2,3)-(jc-kc)*(nc(js)-1)*ntjk
      call jandk(a,eta,g4,ipt,is,isf,ic,ls,lsf,lc,js,jsf,jc,js,jsf,kc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      ndxl=1
      ngtl=nt(jsf)**2
      iblbil=0
      do 98 ircr=1,ircru
      do 96 lrcr=1,lrcru
      do 94 jrcr=1,jrcru
      if(jceqkc) krcru=jrcr
      do 92 krcr=1,krcru
      if(jrcr-krcr) 70, 72, 74
   70 iblb=iblbil+(ipq(krcr)+(jrcr-1))*jbl(2)-(krcr-1)*jbld
      kesfb=2
      kcxbg1=kcxbg(3,3)
      ngtj=1
      ngtk=nt(jsf)
      go to 78
   72 iblb=iblbil+ipq(krcr)*jbl(2)+(krcr-1)*jbl(1)
      kesfb=1
      go to 76
   74 iblb=iblbil+(ipq(jrcr)+(krcr-1))*jbl(2)-(jrcr-1)*jbld
      kesfb=2
   76 kcxbg1=kcxbg(kesfb,3)
      ngtj=nt(jsf)
      ngtk=1
   78 icxs(2,3)=icxbg(3)
      do 88 it=itl,itu
        do 86 lt=ltl,ltu
          kcxs(kesfb,3)=kcxbg1
          ndxj=ndxl
          do 84 jt=jtl,jtu
            if(jrcr.eq.krcr.and.jceqkc) ktu=jt
            ndxk=ndxj
            do 82 kt=ktl,ktu
              val=g4(ndxk)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 3, val )
   80         ndxk=ndxk+ngtk
              kcxs(kesfb,3)=kcxs(kesfb,3)+nprk(kesfb,3)
   82       continue
            ndxj=ndxj+ngtj
   84     continue
          ndxl=ndxl+ngtl
          icxs(2,3)=icxs(2,3)+npri(2,3)
   86   continue
   88 continue
   92 continue
   94 continue
      iblbil=iblbil+jbljk
   96 continue
   98 continue
      return
      end
*deck sf3eq7
      subroutine sf3eq7( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      iscmci,
     &         is,       isf,      js,       jsf,      ks,
     &         ksf,      ijgt,     ijx,      ijy,      ijz,
     &         klgt,     klx,      kly,      klz,      lmnp1,
     &         lmnv,     mcons,    nc,       ncon,     nfct,
     &         ngw,      nrcr,     nt,       ntl,      ntu,
     &         x,        y,        z,        zet )
      implicit real*8 (a-h,o-z)
      logical kceqlc
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
      dimension icxbg(3), jbl(2), kcxbg(3,3)
c
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      kceqlc=kc.eq.lc
      do 12 jesfb=iscmci,2
      jbl(jesfb)=0
      do 10 ist=1,nst
        jbl(jesfb)=jbl(jesfb)+
     &    npriri(2,iscmci,ist)*nprirk(jesfb,iscmci,ist)
   10   continue
   12 continue
      if(iscmci.eq.2) go to 60
c
      jbld=jbl(2)-jbl(1)
      jblkl=ipq(krcru)*jbl(2)+krcru*jbl(1)
      fnfct=(nwt*nfct(icc))
      icxbg(1)=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      kcxbg(1,1)=klcxst(1)+ipq((kc-1)*nt(ksf)+1)*nprk(1,1)+(lc-1)*ntkl1
      kcxbg(2,1)=klcxst(2)+((kc-1)*nc(ks)+(lc-1))*ntkl
      kcxbg(3,1)=kcxbg(2,1)-(kc-lc)*(nc(ks)-1)*ntkl
      call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,jc,ks,ksf,kc,ks,ksf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      ndxj=1
      ngtj=nt(ksf)**2
      iblbij=0
      do 48 ircr=1,ircru
      do 46 jrcr=1,jrcru
      do 44 krcr=1,krcru
      if(kceqlc) lrcru=krcr
      do 42 lrcr=1,lrcru
      if(krcr-lrcr) 18, 20, 22
   18 iblb=iblbij+(ipq(lrcr)+(krcr-1))*jbl(2)-(lrcr-1)*jbld
      kesfb=2
      kcxbg1=kcxbg(3,1)
      ngtk=1
      ngtl=nt(ksf)
      go to 26
   20 iblb=iblbij+ipq(krcr)*jbl(2)+(krcr-1)*jbl(1)
      kesfb=1
      go to 24
   22 iblb=iblbij+(ipq(krcr)+(lrcr-1))*jbl(2)-(krcr-1)*jbld
      kesfb=2
   24 kcxbg1=kcxbg(kesfb,1)
      ngtk=nt(ksf)
      ngtl=1
   26 icxs(2,1)=icxbg(1)
      do 40 it=itl,itu
        do 38 jt=jtl,jtu
          kcxs(kesfb,1)=kcxbg1
          ndxk=ndxj
          do 36 kt=ktl,ktu
            if(krcr.eq.lrcr.and.kceqlc) ltu=kt
            ndxl=ndxk
            do 34 lt=ltl,ltu
              val=g4(ndxl)
              if(abs(val).le.cutoff) go to 30
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 1, val )
   30         ndxl=ndxl+ngtl
              kcxs(kesfb,1)=kcxs(kesfb,1)+nprk(kesfb,1)
   34       continue
            ndxk=ndxk+ngtk
   36     continue
          ndxj=ndxj+ngtj
          icxs(2,1)=icxs(2,1)+npri(2,1)
   38   continue
   40 continue
      continue
   42 continue
   44 continue
      iblbij=iblbij+jblkl
   46 continue
   48 continue
      return
c
   60 icxbg(2)=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      icxbg(3)=icxbg(2)+(lc-kc)*ntik
      kcxbg(2,3)=jkcxst(2)+((jc-1)*nc(ks)+(kc-1))*ntjk
      kcxbg(2,2)=kcxbg(2,3)+(lc-kc)*ntjk
      if(kceqlc) then
        fnfct=(nwt*nfct(icc))
      else
        do 64 igw=2,ngw(icc)
          if(icb(inx,igw,icc).eq.ic.and.
     &       icb(jnx,igw,icc).eq.jc.and.
     &       icb(knx,igw,icc).eq.lc) go to 66
   64   continue
        iscmu = 3
        fnfct=(nwt*nfct(icc))
        go to 68
   66   fnfct=(2*nwt*nfct(icc))
      endif
      iscmu = 2
   68 do 98 iscm=2,iscmu
      if(iscm.eq.2) then
        call jandk(a,eta,g4,ipt,is,isf,ic,ks,ksf,kc,js,jsf,jc,ks,ksf,lc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      else
        call jandk(a,eta,g4,ipt,is,isf,ic,ks,ksf,lc,js,jsf,jc,ks,ksf,kc,
     1    ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2    nrcr,nt,ntl,ntu,x,y,z,zet)
      endif
      ndx=0
      iblb=0
      do 96 ircr=1,ircru
      do 94 krcr=1,krcru
      do 92 jrcr=1,jrcru
      do 90 lrcr=1,lrcru
      icxs(2,iscm)=icxbg(iscm)
      do 88 it=itl,itu
        do 86 kt=ktl,ktu
          kcxs(2,iscm)=kcxbg(2,iscm)
          do 84 jt=jtl,jtu
            do 82 lt=ltl,ltu
              ndx=ndx+1
              val=g4(ndx)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, iscm, val )
   80         kcxs(2,iscm)=kcxs(2,iscm)+nprk(2,2)
   82       continue
   84     continue
          icxs(2,iscm)=icxs(2,iscm)+npri(2,2)
   86   continue
   88 continue
      iblb=iblb+jbl(2)
   90 continue
   92 continue
   94 continue
   96 continue
   98 continue
      return
      end
*deck sf3eq8
      subroutine sf3eq8( a,        cx,       eta,      g4,
     &         h4,       icb,      icc,      ipt,      iscm,
     &         is,       isf,      js,       jsf,      ks,
     &         ksf,      ls,       lsf,      ijgt,     ijx,
     &         ijy,      ijz,      klgt,     klx,      kly,
     &         klz,      lmnp1,    lmnv,     mcons,    nc,
     &         ncon,     nfct,     nrcr,     nt,       ntl,
     &         ntu,      x,        y,        z,        zet )
      implicit real*8 (a-h,o-z)
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
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension a(*),cx(*),eta(mrcru,mconu,*),g4(*),h4(*),icb(4,24,*),
     1  ipt(*),ijgt(*),ijx(*),ijy(*),ijz(*),klgt(*),klx(*),kly(*),
     2  klz(*),lmnp1(*),lmnv(3,*),mcons(*),nc(*),ncon(*),nfct(*),
     3  nrcr(*),nt(*),ntl(*),ntu(*),x(mcu,*),y(mcu,*),z(mcu,*),
     4  zet(mconu,*)
c
      fnfct=nwt*nfct(icc)
      ic=icb(inx,1,icc)
      jc=icb(jnx,1,icc)
      kc=icb(knx,1,icc)
      lc=icb(lnx,1,icc)
      jbl2=0
      do 10 ist=1,nst
        jbl2=jbl2+npriri(2,iscm,ist)*nprirk(2,iscm,ist)
   10 continue
      ndx=0
      iblb=0
      if(iscm-2) 18, 48, 78
c
   18 icxbg=ijcxst(2)+((ic-1)*nc(js)+(jc-1))*ntij
      kcxbg=klcxst(2)+((kc-1)*nc(ls)+(lc-1))*ntkl
      call jandk(a,eta,g4,ipt,is,isf,ic,js,jsf,jc,ks,ksf,kc,ls,lsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      do 38 ircr=1,ircru
      do 36 jrcr=1,jrcru
      do 34 krcr=1,krcru
      do 32 lrcr=1,lrcru
      icxs(2,1)=icxbg
      do 30 it=itl,itu
        do 28 jt=jtl,jtu
          kcxs(2,1)=kcxbg
          do 26 kt=ktl,ktu
            do 24 lt=ltl,ltu
              ndx=ndx+1
              val=g4(ndx)
              if(abs(val).le.cutoff) go to 20
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 1, val )
   20         kcxs(2,1)=kcxs(2,1)+nprk(2,1)
   24       continue
   26     continue
          icxs(2,1)=icxs(2,1)+npri(2,1)
   28   continue
   30 continue
      iblb=iblb+jbl2
   32 continue
   34 continue
   36 continue
   38 continue
      return
c
   48 icxbg=ikcxst(2)+((ic-1)*nc(ks)+(kc-1))*ntik
      kcxbg=jlcxst(2)+((jc-1)*nc(ls)+(lc-1))*ntjl
      call jandk(a,eta,g4,ipt,is,isf,ic,ks,ksf,kc,js,jsf,jc,ls,lsf,lc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      do 68 ircr=1,ircru
      do 66 krcr=1,krcru
      do 64 jrcr=1,jrcru
      do 62 lrcr=1,lrcru
      icxs(2,2)=icxbg
      do 60 it=itl,itu
        do 58 kt=ktl,ktu
          kcxs(2,2)=kcxbg
          do 56 jt=jtl,jtu
            do 54 lt=ltl,ltu
              ndx=ndx+1
              val=g4(ndx)
              if(abs(val).le.cutoff) go to 50
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 2, val )
   50         kcxs(2,2)=kcxs(2,2)+nprk(2,2)
   54       continue
   56     continue
          icxs(2,2)=icxs(2,2)+npri(2,2)
   58   continue
   60 continue
      iblb=iblb+jbl2
   62 continue
   64 continue
   66 continue
   68 continue
      return
c
   78 icxbg=ilcxst(2)+((ic-1)*nc(ls)+(lc-1))*ntil
      kcxbg=jkcxst(2)+((jc-1)*nc(ks)+(kc-1))*ntjk
      call jandk(a,eta,g4,ipt,is,isf,ic,ls,lsf,lc,js,jsf,jc,ks,ksf,kc,
     1  ijgt,ijx,ijy,ijz,klgt,klx,kly,klz,lmnp1,lmnv,mcons,ncon,
     2  nrcr,nt,ntl,ntu,x,y,z,zet)
      do 98 ircr=1,ircru
      do 96 lrcr=1,lrcru
      do 94 jrcr=1,jrcru
      do 92 krcr=1,krcru
      icxs(2,3)=icxbg
      do 90 it=itl,itu
        do 88 lt=ltl,ltu
          kcxs(2,3)=kcxbg
          do 86 jt=jtl,jtu
            do 84 kt=ktl,ktu
              ndx=ndx+1
              val=g4(ndx)
              if(abs(val).le.cutoff) go to 80
              val=fnfct*val
              ibl=iblb
              call aoso2( h4, cx, 3, val )
   80         kcxs(2,3)=kcxs(2,3)+nprk(2,3)
   84       continue
   86     continue
          icxs(2,3)=icxs(2,3)+npri(2,3)
   88   continue
   90 continue
      iblb=iblb+jbl2
   92 continue
   94 continue
   96 continue
   98 continue
      return
      end
*deck jandk
      subroutine jandk( a,      eta,    g4,     ipt,    is,     isf,
     &  ic,     js,     jsf,    jc,     ks,     ksf,    kc,     ls,
     &  lsf,    lc,     ijgt,   ijx,    ijy,    ijz,    klgt,   klx,
     &  kly,    klz,    lmnp1,  lmnv,   mcons,  ncon,   nrcr,   nt,
     &  ntl,    ntu,    x,      y,      z,      zet )
      implicit real*8 (a-h,o-z)
      logical esfc, esfcij, esfckl, igueq1, jgueq1, kgueq1, lgueq1
      parameter (a1s2=0.5d0)
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
      common /shlnos/lit,ljt,lkt,llt,lklt,ljklt,ij,kl,ijkl
      common /shlinf/ xi,yi,zi,rij,xk,yk,zk,rkl,igu,jgu,kgu,lgu,nc4,
     1  nc3,nc2,nc1,igueq1,jgueq1,kgueq1,lgueq1
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      common /setinta/ dxij,dyij,dzij,dxkl,dykl,dzkl,bp01,b00,b10,xcp00,
     1  xc00,ycp00,yc00,zcp00,zc00,f00,ni,nj,nk,nl,nmax,mmax
      common /dim21/ ipq(256)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer  atebyt, forbyt
      external atebyt, forbyt
c
      dimension a(*), eta(mrcru,mconu,*), g4(*), ipt(*), ijgt(*),
     1  ijx(*), ijy(*), ijz(*), klgt(*), klx(*), kly(*), klz(*),
     2  lmnp1(*), lmnv(3,*), mcons(*), ncon(*), nrcr(*), nt(*), ntl(*),
     3  ntu(*), x(mcu,*), y(mcu,*), z(mcu,*), zet(mconu,*)
c
      icons=mcons(isf)
      lit=lmnp1(icons)
      jcons=mcons(jsf)
      ljt=lmnp1(jcons)
      kcons=mcons(ksf)
      lkt=lmnp1(kcons)
      lcons=mcons(lsf)
      llt=lmnp1(lcons)
      lklt=lkt*llt
      ljklt=ljt*lklt
      ngtl=1
      ngtk=nt(lsf)
      ngtj=ngtk*nt(ksf)
      ngti=ngtj*nt(jsf)
c
c  ishell
c
      xi=x(ic,is)
      yi=y(ic,is)
      zi=z(ic,is)
      igu=ncon(icons)
      mini=ntl(isf)
      maxi=ntu(isf)
      ircru=nrcr(icons)
c
c  jshell
c
      xj=x(jc,js)
      yj=y(jc,js)
      zj=z(jc,js)
      jgu=ncon(jcons)
      minj=ntl(jsf)
      maxj=ntu(jsf)
      jrcru=nrcr(jcons)
c
c  prepare items for pairs of (i,j) functions
c
      esfcij=isf.eq.jsf.and.ic.eq.jc
      xij=a1s2*(xi+xj)
      dxij=xi-xj
      xijm=a1s2*dxij
      yij=a1s2*(yi+yj)
      dyij=yi-yj
      yijm=a1s2*dyij
      zij=a1s2*(zi+zj)
      dzij=zi-zj
      zijm=a1s2*dzij
      rij=dxij*dxij+dyij*dyij+dzij*dzij
      ij=0
      do 150 i=mini,maxi
        nx=ljklt*lmnv(1,i) + 1
        ny=ljklt*lmnv(2,i) + 1
        nz=ljklt*lmnv(3,i) + 1
        do 140 j=minj,maxj
          ij=ij+1
          ijx(ij)=nx+lklt*lmnv(1,j)
          ijy(ij)=ny+lklt*lmnv(2,j)
          ijz(ij)=nz+lklt*lmnv(3,j)
          ijgt(ij)=ngti*(i-mini)+ngtj*(j-minj)+1
  140   continue
  150 continue
c
c  kshell
c
      xk=x(kc,ks)
      yk=y(kc,ks)
      zk=z(kc,ks)
      kgu=ncon(kcons)
      mink=ntl(ksf)
      maxk=ntu(ksf)
      krcru=nrcr(kcons)
c
c  lshell
c
      xl=x(lc,ls)
      yl=y(lc,ls)
      zl=z(lc,ls)
      lgu=ncon(lcons)
      minl=ntl(lsf)
      maxl=ntu(lsf)
      lrcru=nrcr(lcons)
c
c  prepare items for pairs of (k,l) functions
c
      esfckl=ksf.eq.lsf.and.kc.eq.lc
      xkl=a1s2*(xk+xl)
      dxkl=xk-xl
      xklm=a1s2*dxkl
      ykl=a1s2*(yk+yl)
      dykl=yk-yl
      yklm=a1s2*dykl
      zkl=a1s2*(zk+zl)
      dzkl=zk-zl
      zklm=a1s2*dzkl
      rkl=dxkl*dxkl+dykl*dykl+dzkl*dzkl
      kl=0
      do 350 k=mink,maxk
        nx=llt*lmnv(1,k)
        ny=llt*lmnv(2,k)
        nz=llt*lmnv(3,k)
        do 340 l=minl,maxl
          kl=kl+1
          klx(kl)=nx+lmnv(1,l)
          kly(kl)=ny+lmnv(2,l)
          klz(kl)=nz+lmnv(3,l)
          klgt(kl)=ngtk*(k-mink)+ngtl*(l-minl)
  340   continue
  350 continue
c
c  prepare miscellaneous items
c
      nroots=(lit+ljt+lkt+llt-2)/2
      esfc=isf.eq.ksf.and.ic.eq.kc.and.jsf.eq.lsf.and.jc.eq.lc
      ijkl=ij*kl
      ngout=ngti*nt(isf)
c
c  allocate space for transforming ao integrals
c
      igueq1=igu.eq.1
      jgueq1=jgu.eq.1
      kgueq1=kgu.eq.1
      lgueq1=lgu.eq.1
      nc1=lrcru*ngout
      if(esfckl) then
        nc2=ipq(krcru+1)*ngout
      else
        nc2=krcru*nc1
        if(igueq1.and.esfc.and.(.not.jgueq1)) then
          nc3=ipq(jrcru+1)*ngout
          nc4=nc3
          go to 660
        endif
      endif
      nc3=jrcru*nc2
      if(esfc) then
        if(esfcij) then
          nc4=ipq(ircru+1)
        else
          nc4=ircru*jrcru
        endif
        nc4=ipq(nc4+1)*ngout
      else
        if(esfcij) then
          nc4=ipq(ircru+1)*nc2
        else
          nc4=ircru*nc3
        endif
      endif
c
  660 if(igueq1) then
        ic3=0
      else
        ic3=nc4
      endif
      if(jgueq1) then
        ic2=0
      else
        ic2=nc3
      endif
      if(kgueq1) then
        ic1=0
      else
        ic1=nc2
      endif
      if(lgueq1) then
        ic0=0
      else
        ic0=nc1
      endif
c
      nxyzin = nroots*lit*ljklt
      n1 = nroots + 1
c     # allocate space.
c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:il(1:nnbft), 3:cx(1:mcxu),
c     #          4:buffer(1:l2rec), 5:values(1:n2max),
c     #          6:labels(1:4,1:n2max), 7:ibitv(1:n2max), 8:h4(1:nblu),
c     #          9:ijgt(1:iju(iscm)), 10:ijx(1:iju(iscm)),
c     #          11:ijy(1:iju(iscm)), 12:ijz(1:iju(iscm)),
c     #          13:klgt(1:klu(iscm)), 14:klx(1:klu(iscm)),
c     #          15:kly(1:klu(iscm)), 16:klz(1:klu(iscm)), 17:g4(1:ic3),
c     #          18:g3(1:ic2), 19:g2(1:ic1), 20:g1(1:ic0),
c     #          21:gout(1:ngout), 22:in1(1:lit+ljt-1),
c     #          23:in(1:lit+ljt-1), 24:kn(1:lkt+llt-1),
c     #          25:uf(1:nroots), 26:wf(1:nroots), 27:ff(1:nroots+n1),
c     #          (cf and sf must be kept together for droot)
c     #          28:cf(1:n1,1:n1), 29:sf(1:n1,1:n1), 30:af(1:n1),
c     #          31:rt(1:n1), 32:r(1:nroots,1:nroots),
c     #          (34-36 use same space as 27-33)
c     #          33:w(1:nroots,1:nroots), 34:xin(1:nxyzin),
c     #          35:yin(1:nxyzin), 36:zin(1:nxyzin)
      ipt(18) = ipt(17) + atebyt( ic3 )
      ipt(19) = ipt(18) + atebyt( ic2 )
      ipt(20) = ipt(19) + atebyt( ic1 )
      ipt(21) = ipt(20) + atebyt( ic0 )
      ipt(22) = ipt(21) + atebyt( ngout )
      ipt(23) = ipt(22) + forbyt( lit + ljt - 1 )
      ipt(24) = ipt(23) + forbyt( lit + ljt - 1 )
      ipt(25) = ipt(24) + forbyt( lkt + llt - 1 )
      ipt(26) = ipt(25) + atebyt( nroots )
      ipt(27) = ipt(26) + atebyt( nroots )
      ipt(28) = ipt(27) + 2*atebyt( nroots + n1 )
      ipt(29) = ipt(28) + 2*atebyt( n1**2 )
      ipt(30) = ipt(29) + 2*atebyt( n1**2 )
      ipt(31) = ipt(30) + 2*atebyt( n1 )
      ipt(32) = ipt(31) + 2*atebyt( n1 )
      ipt(33) = ipt(32) + 2*atebyt( nroots**2 )
      ipt(34) = ipt(27)
      ipt(35) = ipt(34) + atebyt( nxyzin )
      ipt(36) = ipt(35) + atebyt( nxyzin )
      ipt(37) = max( ipt(33) + 2*atebyt( nroots**2 ),
     &               ipt(36) + atebyt( nxyzin ) )
      if ( (ipt(37)-1) .gt. mblu ) then
        call bummer('jandk: mblu too small ',(ipt(37)-1-mblu),faterr)
      endif
c     # set the high-water mark for g4(*) usage.
      call h2oset( (ipt(37)-ipt(17)) )
c
c  compute two-electron integrals
c
      call wzero(nc4,g4,1)
      if(ijkl.eq.1) then
            call s0000( eta(1,1,icons), eta(1,1,jcons), eta(1,1,kcons),
     &  eta(1,1,lcons), g4,             a(ipt(18)),     a(ipt(19)),
     &  a(ipt(20)),     a(ipt(21)),     zet )
      else
            call genral( a,              eta(1,1,icons), eta(1,1,jcons),
     &   eta(1,1,kcons), eta(1,1,lcons), g4,             a(ipt(18)),
     &   a(ipt(19)),     a(ipt(20)),     a(ipt(21)),     ijgt,
     &   ijx,            ijy,            ijz,            klgt,
     &   klx,            kly,            klz,            a(ipt(22)),
     &   a(ipt(23)),     a(ipt(24)),     ipt,            nroots,
     &   a(ipt(25)),     a(ipt(26)),     n1,             a(ipt(34)),
     &   a(ipt(35)),     a(ipt(36)),     zet )
      endif
      return
      end
*deck s0000
      subroutine s0000( etai,           etaj,          etak,
     &  etal,           g4,             g3,            g2,
     &  g1,             gout,           zet )
      implicit real*8 (a-h,o-z)
      logical esfc, esfcij, esfckl, igueq1, jgueq1, kgueq1, lgueq1
      parameter (a1s2=0.5d0, asrtpi=1.1283791670955126d0,
     1  pie4=7.85398163397448d-01)
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /shlinf/ xi,yi,zi,rij,xk,yk,zk,rkl,igu,jgu,kgu,lgu,nc4,
     1  nc3,nc2,nc1,igueq1,jgueq1,kgueq1,lgueq1
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      dimension etai(mrcru,*),etaj(mrcru,*),etak(mrcru,*),
     1  etal(mrcru,*),g4(*),g3(*),g2(*),g1(*),gout(*),zet(mconu,*)
c
      abx=xij-xkl
      aby=yij-ykl
      abz=zij-zkl
      if(igueq1) then
        fctr1=asrtpi*etai(1,1)
      else
        fctr1=asrtpi
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
      if(kgueq1) fctr1=fctr1*etak(1,1)
      if(lgueq1) fctr1=fctr1*etal(1,1)
      kgl=1
c
c  i primitive
c
      do 550 ig=1,igu
      if(.not.igueq1) call wzero(nc3,g3,1)
      fctr2=fctr1
      ai=zet(ig,icons)
      arij=ai*rij
      if(esfc) kgl=ig
c
c  j primitive
c
      if(esfcij) jgu=ig
      do 350 jg=1,jgu
      if(.not.jgueq1) call wzero(nc2,g2,1)
      if((esfcij.and..not.igueq1).and.jg.eq.ig) fctr2=a1s2*fctr2
      aj=zet(jg,jcons)
      aa=ai+aj
      aarij=aj*arij/aa
      if(aarij.gt.tol) go to 350
      fctr3=fctr2/aa
      aaa=(ai-aj)/aa
      px=abx+aaa*xijm
      py=aby+aaa*yijm
      pz=abz+aaa*zijm
c
c  k primitive
c
      do 250 kg=kgl,kgu
      if(.not.kgueq1) call wzero(nc1,g1,1)
      bk=zet(kg,kcons)
      brkl=bk*rkl
c
c  l primitive
c
      if(esfc.and.(kg.eq.ig)) then
        lgl=jg
      else
        lgl=1
      endif
      if(esfckl) lgu=kg
      do 150 lg=lgl,lgu
      bl=zet(lg,lcons)
      bb=bk+bl
      dum=aarij+bl*brkl/bb
      if(dum.gt.tol) go to 150
      expe=fctr3*exp(-dum)/bb
      if(((esfc.and..not.(igueq1.and.jgueq1)).and.(kg.eq.ig)).and.
     1  lg.eq.jg) expe=a1s2*expe
      if((esfckl.and..not.kgueq1).and.lg.eq.kg) expe=a1s2*expe
      bbb=(bk-bl)/bb
      apb=aa+bb
      x=((px-bbb*xklm)**2+(py-bbb*yklm)**2+(pz-bbb*zklm)**2)*aa*bb/apb
      if(x.gt.5.0d+00) go to 50
      if(x.gt.1.0d+00) go to 30
      if(x.gt.3.0d-7) go to 20
      ww1=1.0d+00-x/3.0d+00
      go to 100
   20 f1=          ((((((((-8.36313918003957d-08*x+1.21222603512827d-06
     1)*x-1.15662609053481d-05 )*x+9.25197374512647d-05
     2)*x-6.40994113129432d-04 )*x+3.78787044215009d-03
     3)*x-1.85185172458485d-02 )*x+7.14285713298222d-02
     4)*x-1.99999999997023d-01 )*x+3.33333333333318d-01
      ww1=(x+x)*f1+exp(-x)
      go to 100
   30 if(x.gt.3.0d+00) go to 40
      y=x-2.0d+00
      f1=        ((((((((((-1.61702782425558d-10*y+1.96215250865776d-09
     1)*y-2.14234468198419d-08 )*y+2.17216556336318d-07
     2)*y-1.98850171329371d-06 )*y+1.62429321438911d-05
     3)*y-1.16740298039895d-04 )*y+7.24888732052332d-04
     4)*y-3.79490003707156d-03 )*y+1.61723488664661d-02
     5)*y-5.29428148329736d-02 )*y+1.15702180856167d-01
      ww1=(x+x)*f1+exp(-x)
      go to 100
   40 y=x-4.0d+00
      f1=        ((((((((((-2.62453564772299d-11*y+3.24031041623823d-10
     1)*y-3.614965656163d-09)*y+3.760256799971d-08)*y-3.553558319675d-07
     2)*y+3.022556449731d-06)*y-2.290098979647d-05)*y+1.526537461148d-04
     3)*y-8.81947375894379d-04 )*y+4.33207949514611d-03
     4)*y-1.75257821619926d-02 )*y+5.28406320615584d-02
      ww1=(x+x)*f1+exp(-x)
      go to 100
   50 if(x.gt.15.0d+00) go to 70
      e=exp(-x)
      if(x.gt.10.0d+00) go to 60
      ww1=    (((((( 4.6897511375022d-01/x-6.9955602298985d-01)/x
     1+5.3689283271887d-01)/x-3.2883030418398d-01)/x
     2+2.4645596956002d-01)/x-4.9984072848436d-01)/x
     3-3.1501078774085d-06)*e + sqrt(pie4/x)
      go to 100
   60 ww1=       (((-1.8784686463512d-01/x+2.2991849164985d-01)/x
     1-4.9893752514047d-01)/x-2.1916512131607d-05)*e + sqrt(pie4/x)
      go to 100
   70 if(x.gt.33.0d+00) go to 90
      e=exp(-x)
      ww1=        (( 1.9623264149430d-01/x-4.9695241464490d-01)/x
     1-6.0156581186481d-05)*e + sqrt(pie4/x)
      go to 100
   90 ww1=sqrt(pie4/x)
  100 gout(1)=ww1*expe/sqrt(apb)
c
c  l transformation
c
      if(lgueq1) go to 200
c
      do 140 lrcr=1,lrcru
        g1(lrcr)=g1(lrcr)+gout(1)*etal(lrcr,lg)
  140 continue
  150 continue
c
c  k transformation
c
  200 if(kgueq1) go to 280
c
      kl2=0
      do 240 krcr=1,krcru
        if(esfckl) then
          do 210 lrcr=1,krcr
            g2(kl2+lrcr)=g2(kl2+lrcr)+g1(lrcr)*etak(krcr,kg)+
     1                                g1(krcr)*etak(lrcr,kg)
  210     continue
          kl2=kl2+krcr
        else
          do 230 lrcr=1,lrcru
            g2(kl2+lrcr)=g2(kl2+lrcr)+g1(lrcr)*etak(krcr,kg)
  230     continue
          kl2=kl2+lrcru
        endif
  240 continue
  250 continue
c
c  j transformation
c
  280 if(jgueq1) go to 410
c
      jkl3=0
      if(esfc.and.igueq1) then
        do 310 jrcr=1,jrcru
          do 300 lrcr=1,jrcr
            g3(jkl3+lrcr)=g3(jkl3+lrcr)+g2(lrcr)*etaj(jrcr,jg)+
     1                                  g2(jrcr)*etaj(lrcr,jg)
  300     continue
          jkl3=jkl3+jrcr
  310   continue
      else
        do 340 jrcr=1,jrcru
          kl2=0
          do 330 krcr=1,krcru
            if(esfckl) lrcru=krcr
            do 320 lrcr=1,lrcru
              g3(jkl3+lrcr)=g3(jkl3+lrcr)+g2(kl2+lrcr)*etaj(jrcr,jg)
  320       continue
            jkl3=jkl3+lrcru
            kl2=kl2+lrcru
  330     continue
  340   continue
      endif
  350 continue
c
c  i transformation
c
  410 if(igueq1) return
c
      if(esfc) ij3=-nc2+1
      if(esfcij) ikl3i=0
      ijkl4=0
      do 550 ircr=1,ircru
        jkl3j=0
        if(esfcij) jrcru=ircr
        do 540 jrcr=1,jrcru
          jkl3=jkl3j
          if(esfcij) ikl3=ikl3i
          if(esfc) krcru=ircr
          do 520 krcr=1,krcru
            if(esfc) then
              if(esfckl) then
                if(krcr.eq.ircr) then
                  lrcru=jrcr
                else
                  lrcru=krcr
                endif
                do 450 lrcr=1,lrcru
                  g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+
     1                           g3(jkl3+lrcr)*etai(ircr,ig)+
     2                           g3(ikl3+lrcr)*etai(jrcr,ig)+
     3                           g3(ij3+nc2*lrcr)*etai(krcr,ig)+
     4                           g3(ij3+nc2*krcr)*etai(lrcr,ig)
  450           continue
                ijkl4=ijkl4+lrcru
                jkl3=jkl3+krcr
                ikl3=ikl3+krcr
              else
                if(krcr.eq.ircr) then
                  lrcru=jrcr
                else
                  lrcru=jrcru
                endif
                do 470 lrcr=1,lrcru
                  g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+
     1                           g3(jkl3+lrcr)*etai(ircr,ig)+
     2                           g3(ij3+nc2*lrcr)*etai(krcr,ig)
  470           continue
                ijkl4=ijkl4+lrcru
                jkl3=jkl3+jrcru
              endif
            else
              if(esfckl) lrcru=krcr
              if(esfcij) then
                do 490 lrcr=1,lrcru
                  g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+
     1                           g3(jkl3+lrcr)*etai(ircr,ig)+
     2                           g3(ikl3+lrcr)*etai(jrcr,ig)
  490           continue
                ijkl4=ijkl4+lrcru
                jkl3=jkl3+lrcru
                ikl3=ikl3+lrcru
              else
                do 510 lrcr=1,lrcru
                  g4(ijkl4+lrcr)=g4(ijkl4+lrcr)+
     1                           g3(jkl3+lrcr)*etai(ircr,ig)
  510           continue
                ijkl4=ijkl4+lrcru
                jkl3=jkl3+lrcru
              endif
            endif
  520     continue
          jkl3j=jkl3j+nc2
          if(esfc) ij3=ij3+1
  540   continue
        if(esfcij) ikl3i=ikl3i+nc2
  550 continue
      return
      end
*deck genral
      subroutine genral( a,              etai,           etaj,
     &   etak,           etal,           g4,             g3,
     &   g2,             g1,             gout,           ijgt,
     &   ijx,            ijy,            ijz,            klgt,
     &   klx,            kly,            klz,            in1,
     &   in,             kn,             ipt,            nroots,
     &   uf,             wf,             n1,             xin,
     &   yin,            zin,            zet )
      implicit real*8 (a-h,o-z)
      logical esfc, esfcij, esfckl, igueq1, jgueq1, kgueq1, lgueq1
      parameter (a1s2=0.5d0, asrtpi=1.1283791670955126d0)
      common /parmr/ cutoff, tol
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      common /shlinf/ xi,yi,zi,rij,xk,yk,zk,rkl,igu,jgu,kgu,lgu,nc4,
     1  nc3,nc2,nc1,igueq1,jgueq1,kgueq1,lgueq1
      common /shlnos/lit,ljt,lkt,llt,lklt,ljklt,ij,kl,ijkl
      common/misc/xij,xijm,yij,yijm,zij,zijm,xkl,xklm,ykl,yklm,zkl,zklm,
     1icons,ircru,jcons,jrcru,kcons,krcru,lcons,lrcru,esfc,esfcij,esfckl
      common /setinta/ dxij,dyij,dzij,dxkl,dykl,dzkl,bp01,b00,b10,xcp00,
     1  xc00,ycp00,yc00,zcp00,zc00,f00,ni,nj,nk,nl,nmax,mmax
      dimension a(*),etai(mrcru,*),etaj(mrcru,*),etak(mrcru,*),
     1  etal(mrcru,*),g4(*),g3(*),g2(*),g1(*),gout(*),ijgt(*),ijx(*),
     2  ijy(*),ijz(*),klgt(*),klx(*),kly(*),klz(*),in1(*),in(*),kn(*),
     3  ipt(*),uf(*),wf(*),xin(*),yin(*),zin(*),zet(mconu,*)
c
      lijklt=lit*ljklt
      mmu=nroots*lijklt
      ni=lit-1
      nj=ljt-1
      nk=lkt-1
      nl=llt-1
      nmax=ni+nj
      mmax=nk+nl
      ihi=nmax+1
      do 12 i=1,ihi
        if(i.le.lit) then
          in1(i)=ljklt*(i-1)+1
        else
          in1(i)=ljklt*(lit-1)+lklt*(i-lit)+1
        endif
   12 continue
      ihi=mmax+1
      do 22 i=1,ihi
        if(i.le.lkt) then
          kn(i)=llt*(i-1)
        else
          kn(i)=llt*(lkt-1)+i-lkt
        endif
   22 continue
      if(igueq1) then
        fctr1=asrtpi*etai(1,1)
      else
        fctr1=asrtpi
      endif
      if(jgueq1) fctr1=fctr1*etaj(1,1)
      if(kgueq1) fctr1=fctr1*etak(1,1)
      if(lgueq1) fctr1=fctr1*etal(1,1)
      kgl=1
c
c  i primitive
c
      do 550 ig=1,igu
      if(.not.igueq1) call wzero(nc3,g3,1)
      fctr2=fctr1
      ai=zet(ig,icons)
      arij=ai*rij
      if(esfc) kgl=ig
c
c  j primitive
c
      if(esfcij) jgu=ig
      do 350 jg=1,jgu
      if(.not.jgueq1) call wzero(nc2,g2,1)
      if((esfcij.and..not.igueq1).and.jg.eq.ig) fctr2=a1s2*fctr2
      aj=zet(jg,jcons)
      aa=ai+aj
      aarij=aj*arij/aa
      if(aarij.gt.tol) go to 350
      fctr3=fctr2/aa
      aaa=(ai-aj)/aa
      xa=xij+aaa*xijm
      ya=yij+aaa*yijm
      za=zij+aaa*zijm
      axai=aa*(xa-xi)
      ayai=aa*(ya-yi)
      azai=aa*(za-zi)
      axak=aa*(xa-xk)
      ayak=aa*(ya-yk)
      azak=aa*(za-zk)
c
c  k primitive
c
      do 250 kg=kgl,kgu
      if(.not.kgueq1) call wzero(nc1,g1,1)
      bk=zet(kg,kcons)
      brkl=bk*rkl
c
c  l primitive
c
      if(esfc.and.(kg.eq.ig)) then
        lgl=jg
      else
        lgl=1
      endif
      if(esfckl) lgu=kg
      do 150 lg=lgl,lgu
      if(.not.lgueq1) call wzero(ijkl,gout,1)
      bl=zet(lg,lcons)
      bb=bk+bl
      dum=aarij+bl*brkl/bb
      if(dum.gt.tol) go to 150
      apb=aa+bb
      expe=fctr3*exp(-dum)/(bb*sqrt(apb))
      if(((esfc.and..not.(igueq1.and.jgueq1)).and.(kg.eq.ig)).and.
     1  lg.eq.jg) expe=a1s2*expe
      if((esfckl.and..not.kgueq1).and.lg.eq.kg) expe=a1s2*expe
      bbb=(bk-bl)/bb
      xb=xkl+bbb*xklm
      yb=ykl+bbb*yklm
      zb=zkl+bbb*zklm
      atb=aa*bb
      rho=atb/apb
      xx=rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
      bxbi=bb*(xb-xi)
      bybi=bb*(yb-yi)
      bzbi=bb*(zb-zi)
      bxbk=bb*(xb-xk)
      bybk=bb*(yb-yk)
      bzbk=bb*(zb-zk)
      c1x=bxbk+axak
      c2x=aa*bxbk
      c3x=bxbi+axai
      c4x=bb*axai
      c1y=bybk+ayak
      c2y=aa*bybk
      c3y=bybi+ayai
      c4y=bb*ayai
      c1z=bzbk+azak
      c2z=aa*bzbk
      c3z=bzbi+azai
      c4z=bb*azai
c
c  roots and weights for quadrature
c
      if(nroots.le.3) then
        call rt123(nroots,xx,uf,wf)
      elseif(nroots.eq.4) then
        call root4(xx,uf,wf)
      elseif(nroots.eq.5) then
        call root5(xx,uf,wf)
      else
        call droot(nroots,n1,xx,a(ipt(27)),a(ipt(28)),a(ipt(29)),
     &    a(ipt(30)),a(ipt(31)),a(ipt(32)),a(ipt(33)),uf,wf)
      endif
c
c  compute two-electron integrals for each root
c
      mm=0
      ihi=nmax+1
      do 120 m=1,nroots
        u2=uf(m)*rho
        f00=expe*wf(m)
        do 100 i=1,ihi
          in(i)=in1(i)+mm
  100   continue
        dum=atb+u2*apb
        dum2=dum+dum
        bp01=(aa+u2)/dum2
        b00=u2/dum2
        b10=(bb+u2)/dum2
        xcp00=(u2*c1x+c2x)/dum
        xc00 =(u2*c3x+c4x)/dum
        ycp00=(u2*c1y+c2y)/dum
        yc00 =(u2*c3y+c4y)/dum
        zcp00=(u2*c1z+c2z)/dum
        zc00 =(u2*c3z+c4z)/dum
        call xyzint(in,kn,xin,yin,zin)
        mm=mm+lijklt
  120 continue
c
c  form (i,j//k,l) integrals over functions
c
      do 126 i=1,ij
        nx=ijx(i)-lijklt
        ny=ijy(i)-lijklt
        nz=ijz(i)-lijklt
        nijgt=ijgt(i)
        do 124 k=1,kl
          mx=nx+klx(k)
          my=ny+kly(k)
          mz=nz+klz(k)
          n=nijgt+klgt(k)
          do 122 mm=lijklt,mmu,lijklt
            gout(n)=gout(n)+xin(mx+mm)*yin(my+mm)*zin(mz+mm)
  122     continue
  124   continue
  126 continue
c
c  l transformation
c
      if(lgueq1) go to 200
c
      l1=0
      do 140 lrcr=1,lrcru
        do 135 n=1,ijkl
          g1(l1+n)=g1(l1+n)+gout(n)*etal(lrcr,lg)
  135   continue
        l1=l1+ijkl
  140 continue
  150 continue
c
c  k transformation
c
  200 if(kgueq1) go to 280
c
      if(esfckl) k1=0
      kl2=0
      do 240 krcr=1,krcru
        l1=0
        if(esfckl) then
          do 210 lrcr=1,krcr
            do 205 n=1,ijkl
              g2(kl2+n)=g2(kl2+n)+g1(l1+n)*etak(krcr,kg)+
     1                            g1(k1+n)*etak(lrcr,kg)
  205       continue
            kl2=kl2+ijkl
            l1=l1+ijkl
  210     continue
          k1=k1+ijkl
        else
          do 230 lrcr=1,lrcru
            do 225 n=1,ijkl
              g2(kl2+n)=g2(kl2+n)+g1(l1+n)*etak(krcr,kg)
  225       continue
            kl2=kl2+ijkl
            l1=l1+ijkl
  230     continue
        endif
  240 continue
  250 continue
c
c  j transformation
c
  280 if(jgueq1) go to 410
c
      jkl3=0
      if(esfc.and.igueq1) then
        ij2=0
        do 310 jrcr=1,jrcru
          kl2=0
          do 305 lrcr=1,jrcr
            do 300 nij=1,ij
              ij2n=ij2+nij
              do 290 nkl=1,ij
                g3(jkl3+nkl)=g3(jkl3+nkl)+g2(kl2+nkl)*etaj(jrcr,jg)+
     1                                    g2(ij2n)*etaj(lrcr,jg)
                ij2n=ij2n+ij
  290         continue
              jkl3=jkl3+ij
              kl2=kl2+ij
  300       continue
  305     continue
          ij2=ij2+ijkl
  310   continue
      else
        do 340 jrcr=1,jrcru
          kl2=0
          do 335 krcr=1,krcru
            if(esfckl) lrcru=krcr
            do 330 lrcr=1,lrcru
              do 325 n=1,ijkl
                g3(jkl3+n)=g3(jkl3+n)+g2(kl2+n)*etaj(jrcr,jg)
  325         continue
              jkl3=jkl3+ijkl
              kl2=kl2+ijkl
  330       continue
  335     continue
  340   continue
      endif
  350 continue
c
c  i transformation
c
  410 if(igueq1) return
c
      if(esfc) ij3=0
      if(esfcij) ikl3i=0
      ijkl4=0
      do 545 ircr=1,ircru
      jkl3j=0
      if(esfcij) jrcru=ircr
      do 540 jrcr=1,jrcru
      if(esfc) jkl3k=jkl3j
      if(esfcij) ikl3k=ikl3i
      if(esfc.and.esfcij) kij3=ij3
      if(esfc) krcru=ircr
      do 520 krcr=1,krcru
      if(esfc) then
        jkl3=jkl3k
        lij3=ij3
        if(esfckl) then
          ikl3=ikl3k
          if(krcr.eq.ircr) then
            lrcru=jrcr
          else
            lrcru=krcr
          endif
          do 450 lrcr=1,lrcru
            do 445 n=1,ijkl
              g4(ijkl4+n)=g4(ijkl4+n)+g3(jkl3+n)*etai(ircr,ig)+
     1                                g3(ikl3+n)*etai(jrcr,ig)+
     2                                g3(lij3+n)*etai(krcr,ig)+
     3                                g3(kij3+n)*etai(lrcr,ig)
  445       continue
            ijkl4=ijkl4+ijkl
            jkl3=jkl3+ijkl
            ikl3=ikl3+ijkl
            lij3=lij3+nc2
  450     continue
          kij3=kij3+nc2
          jkl3k=jkl3k+(krcr*ijkl)
          ikl3k=ikl3k+(krcr*ijkl)
        else
          if(krcr.eq.ircr) then
            lrcru=jrcr
          else
            lrcru=jrcru
          endif
          do 470 lrcr=1,lrcru
            do 465 nij=1,ij
              lij3n=lij3+nij
              do 463 nkl=1,ij
                g4(ijkl4+nkl)=g4(ijkl4+nkl)+g3(jkl3+nkl)*etai(ircr,ig)+
     1                                      g3(lij3n)*etai(krcr,ig)
                lij3n=lij3n+ij
  463         continue
              ijkl4=ijkl4+ij
              jkl3=jkl3+ij
  465       continue
            lij3=lij3+nc2
  470     continue
          jkl3k=jkl3k+nc1
        endif
      else
        if(esfckl) lrcru=krcr
        if(esfcij) then
          do 490 lrcr=1,lrcru
            do 485 n=1,ijkl
              g4(ijkl4+n)=g4(ijkl4+n)+g3(jkl3j+n)*etai(ircr,ig)+
     1                                g3(ikl3k+n)*etai(jrcr,ig)
  485       continue
            ijkl4=ijkl4+ijkl
            jkl3j=jkl3j+ijkl
            ikl3k=ikl3k+ijkl
  490     continue
        else
          do 510 lrcr=1,lrcru
            do 505 n=1,ijkl
              g4(ijkl4+n)=g4(ijkl4+n)+g3(jkl3j+n)*etai(ircr,ig)
  505       continue
            ijkl4=ijkl4+ijkl
            jkl3j=jkl3j+ijkl
  510     continue
        endif
      endif
  520 continue
      if(esfc) jkl3j=jkl3j+nc2
      if(esfc) ij3=ij3+ijkl
  540 continue
      if(esfcij) ikl3i=ikl3i+nc2
  545 continue
  550 continue
      return
      end
*deck xyzint
      subroutine xyzint(i,k,xint,yint,zint)
      implicit real*8 (a-h,o-z)
      logical n0, n1, m0, m1
      parameter (a0=0.0d0, a1=1.0d0)
      common /setinta/ dxij,dyij,dzij,dxkl,dykl,dzkl,bp01,b00,b10,xcp00,
     1  xc00,ycp00,yc00,zcp00,zc00,f00,nimax,njmax,nkmax,nlmax,nmax,mmax
      common /shlnos/lit,ljt,lkt,llt,lklt,ljklt,ij,kl,ijkl
      dimension i(*), k(*), xint(*), yint(*), zint(*)
c
      n0=nmax.eq.0
      n1=nmax.le.1
      m0=mmax.eq.0
      m1=mmax.le.1
c
c     # i(0,0)
c
      i1=i(1)
      xint(i1)=a1
      yint(i1)=a1
      zint(i1)=f00
      if(n0.and.m0) return
      if(n0) go to 5
c
c     # i(1,0)
c
      i2=i(2)
      xint(i2)=xc00
      yint(i2)=yc00
      zint(i2)=zc00*f00
c
    5 if(m0) go to 10
c
c     # i(0,1)
c
      k2=k(2)
      i3=i1+k2
      xint(i3)=xcp00
      yint(i3)=ycp00
      zint(i3)=zcp00*f00
      if(n0) go to 10
c
c     # i(1,1)
c
      i3=i2+k2
      cp10=b00
      xint(i3)=xcp00*xint(i2)+cp10
      yint(i3)=ycp00*yint(i2)+cp10
      zint(i3)=zcp00*zint(i2)+cp10*f00
c
   10 if(n1) go to 40
      c10=a0
      i3=i1
      i4=i2
      do 30 n=2,nmax
        c10=c10+b10
c
c       # i(n,0)
c
        i5=i(n+1)
        xint(i5)=c10*xint(i3)+xc00*xint(i4)
        yint(i5)=c10*yint(i3)+yc00*yint(i4)
        zint(i5)=c10*zint(i3)+zc00*zint(i4)
        if(m0) go to 20
        cp10=cp10+b00
c
c       # i(n,1)
c
        i3=i5+k2
        xint(i3)=xcp00*xint(i5)+cp10*xint(i4)
        yint(i3)=ycp00*yint(i5)+cp10*yint(i4)
        zint(i3)=zcp00*zint(i5)+cp10*zint(i4)
c
   20   i3=i4
        i4=i5
   30 continue
c
   40 if(m1) go to 70
      cp01=a0
      c01=b00
      i3=i1
      i4=i1+k2
      do 60 m=2,mmax
        cp01=cp01+bp01
c
c       # i(0,m)
c
        i5=i1+k(m+1)
        xint(i5)=cp01*xint(i3)+xcp00*xint(i4)
        yint(i5)=cp01*yint(i3)+ycp00*yint(i4)
        zint(i5)=cp01*zint(i3)+zcp00*zint(i4)
        if(n0) go to 50
        c01=c01+b00
c
c       # i(1,m)
c
        i3=i2+k(m+1)
        xint(i3)=xc00*xint(i5)+c01*xint(i4)
        yint(i3)=yc00*yint(i5)+c01*yint(i4)
        zint(i3)=zc00*zint(i5)+c01*zint(i4)
c
   50   i3=i4
        i4=i5
   60 continue
c
   70 if(n1.or.m1) go to 100
c
c     # i(n,m)
c
      c01=b00
      k3=k2
      do 90 m=2,mmax
        k4=k(m+1)
        c01=c01+b00
        i3=i1
        i4=i2
        c10=b10
        do 80 n=2,nmax
          i5=i(n+1)
          xint(i5+k4)=c10*xint(i3+k4)+xc00*xint(i4+k4)+c01*xint(i4+k3)
          yint(i5+k4)=c10*yint(i3+k4)+yc00*yint(i4+k4)+c01*yint(i4+k3)
          zint(i5+k4)=c10*zint(i3+k4)+zc00*zint(i4+k4)+c01*zint(i4+k3)
          c10=c10+b10
          i3=i4
          i4=i5
   80   continue
        k3=k4
   90 continue
c
  100 if(njmax.eq.0) go to 200
c
c     # i(ni,nj,m)
c
      m=0
      i5=i(nmax+1)
  110 low=nimax
      km=k(m+1)
  120 n=nmax
      i3=i5+km
  130 i4=i(n)+km
      xint(i3)=xint(i3)+dxij*xint(i4)
      yint(i3)=yint(i3)+dyij*yint(i4)
      zint(i3)=zint(i3)+dzij*zint(i4)
      i3=i4
      n=n-1
      if(n.gt.low) go to 130
      low=low+1
      if(low.lt.nmax) go to 120
      if(nimax.eq.0) go to 160
      i3=lklt+km+i1
      do 150 nj=1,njmax
        i4=i3
        do 140 ni=1,nimax
          xint(i4)=xint(i4+ljklt-lklt)+dxij*xint(i4-lklt)
          yint(i4)=yint(i4+ljklt-lklt)+dyij*yint(i4-lklt)
          zint(i4)=zint(i4+ljklt-lklt)+dzij*zint(i4-lklt)
          i4=i4+ljklt
  140   continue
        i3=i3+lklt
  150 continue
c
  160 m=m+1
      if(m.le.mmax) go to 110
c
  200 if(nlmax.eq.0) go to 300
c
c     # i(ni,nj,nk,nl)
c
      i5=k(mmax+1)
      ia=i1
      ni=0
  210 nj=0
      ib=ia
  220 low=nkmax
  230 m=mmax
      i3=ib+i5
  240 i4=ib+k(m)
      xint(i3)=xint(i3)+dxkl*xint(i4)
      yint(i3)=yint(i3)+dykl*yint(i4)
      zint(i3)=zint(i3)+dzkl*zint(i4)
      i3=i4
      m=m-1
      if(m.gt.low) go to 240
      low=low+1
      if(low.lt.mmax) go to 230
      if(nkmax.eq.0) go to 270
      i3=ib+1
      do 260 nl=1,nlmax
        i4=i3
        do 250 nk=1,nkmax
          xint(i4)=xint(i4+llt-1)+dxkl*xint(i4-1)
          yint(i4)=yint(i4+llt-1)+dykl*yint(i4-1)
          zint(i4)=zint(i4+llt-1)+dzkl*zint(i4-1)
          i4=i4+llt
  250   continue
        i3=i3+1
  260 continue
c
  270 nj=nj+1
      ib=ib+lklt
      if(nj.le.njmax) go to 220
      ni=ni+1
      ia=ia+ljklt
      if(ni.le.nimax) go to 210
  300 continue
      return
      end
*deck aoso2e
      subroutine aoso2e( h4, cx, iscm, val )
      implicit real*8 (a-h,o-z)
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      parameter (a0=0.0d0)
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1  jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2  npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
      dimension cx(*), h4(*)
c
      icx = icxs(iesfb,iscm)
      jcx = kcxs(iesfb,iscm)
      do 72 ist=1,nst
      kcxst = jcx
      lcxst = icx
      do 64 ijbl=1,npriri(iesfb,iscm,ist)
      icx=icx+1
      jcx=jcx+1
      if(cx(icx).eq.a0.and.cx(jcx).eq.a0) then
        ibl=ibl+ijbl
      else
        kcx = kcxst
        lcx = lcxst
        if(icx.eq.jcx) then
          do 48 klbl=1,ijbl
            kcx=kcx+1
            ibl=ibl+1
            h4(ibl)=h4(ibl)+(cx(icx)*val)*cx(kcx)
   48     continue
        else
          do 56 klbl=1,ijbl
            kcx=kcx+1
            lcx=lcx+1
            ibl=ibl+1
            h4(ibl)=h4(ibl)+(cx(icx)*val)*cx(kcx)+(cx(jcx)*val)*cx(lcx)
   56     continue
        endif
      endif
   64 continue
   72 continue
      return
      end
*deck aoso2
      subroutine aoso2( h4, cx, iscm, val )
      implicit real*8 (a-h,o-z)
c
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
c
      parameter (a0=0.0d0)
      common /ntgr/itl,itu,jtl,jtu,ktl,ktu,ltl,ltu,jcu,kcu,lcu,inx,
     1  jnx,knx,lnx,nwt,nwt1(3),ijsf(3),klsf(3),icxs(2,3),kcxs(2,3),
     2  npri(2,3),nprk(2,3),iesfb,kesfb,ircru,jrcru,krcru,lrcru
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
      common /sf3ao/ ibl
c
      dimension cx(*), h4(*)
      icx = icxs(iesfb,iscm)
      kcxst = kcxs(kesfb,iscm)
      do 56 ist=1,nst
      npairi = npriri(iesfb,iscm,ist)
      npairk = nprirk(kesfb,iscm,ist)
      if(npairi.eq.0) go to 52
      if(npairk.ne.0) go to 20
      icx=icx+npairi
      go to 56
   20 do 48 ijbl=1,npairi
        icx=icx+1
        if(cx(icx).eq.a0) go to 44
        do 40 klbl=1,npairk
          h4(ibl+klbl)=h4(ibl+klbl)+cx(kcxst+klbl)*(cx(icx)*val)
   40   continue
   44   ibl=ibl+npairk
   48 continue
   52 kcxst = kcxst + npairk
   56 continue
      return
      end
*deck wtint2
      subroutine wtint2( h4,      il,     iprst,  iscm,
     & buffer, values, labels, ibitv )
c
c  assign orbital labels and write records of 2-e integrals.
c
c  01-dec-90 SIFS version. ls(2,*) used as a 2d array. -rls
c
      implicit logical(a-z)
c
      integer    ilfact
      parameter( ilfact=1024 )
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
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
      real*8         cutoff, tol
      common /parmr/ cutoff, tol
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
     & kcxs,           npri,           nprk,           iesfbx, kesfbx,
     & ircrux, jrcrux, krcrux, lrcrux
      common /ntgr/
     & itl,    itu,    jtl,    jtu,    ktl,    ktu,    ltl,    ltu,
     & jcu,    kcu,    lcu,    inx,    jnx,    knx,    lnx,    nwt,
     & nwt1(3),        ijsf(3),        klsf(3),        icxs(2,3),
     & kcxs(2,3),      npri(2,3),      nprk(2,3),      iesfbx, kesfbx,
     & ircrux, jrcrux, krcrux, lrcrux
c
      integer       npriri,        nprirk
      common /ikpr/ npriri(2,3,8), nprirk(2,3,8)
c
      real*8
     & xij,    xijm,   yij,    yijm,   zij,    zijm,   xkl,    xklm,
     & ykl,    yklm,   zkl,    zklm
      integer
     & icons,  ircru,  jcons,  jrcru,  kcons,  krcru,  lcons,  lrcru
      logical
     & esfc,   esfcij, esfckl
      common /misc/
     & xij,    xijm,   yij,    yijm,   zij,    zijm,   xkl,    xklm,
     & ykl,    yklm,   zkl,    zklm,
     & icons,  ircru,  jcons,  jrcru,  kcons,  krcru,  lcons,  lrcru,
     & esfc,   esfcij, esfckl
c
      integer    nipv
      parameter( nipv=4 )
c
c     # dummy:
      integer iscm
      integer il(*), iprst(*), labels(nipv,*), ibitv(*)
      real*8 h4(*), buffer(*), values(*)
c     # ibm+convex f77 bug; actual:
c     # integer labels(nipv,n2max), ibitv( ((n2max+63)/64)*64 )
c     # real*8 buffer(l2rec), values(n2max)
c
c     # local:
      integer ibl, ipr, kprib, ircr, jrcr, iesfb, iprkb, kpr, krcr,
     & lrcr, kesfb, ist, npairi, npairk, kpri, ipair, kpair
      integer    msame
      parameter( msame=0 )
c
      ibl = 0
      ipr = iprst(ijsf(iscm))
      if ( eprsf(iscm) ) then
        kprib = ipr
      else
        kprib = iprst(klsf(iscm))
      endif
      do 176 ircr = 1, ircru
        if ( esfi(iscm) ) jrcru = ircr
        do 174 jrcr = 1, jrcru
          if ( esfi(iscm) .and. (ircr .eq. jrcr) ) then
            iesfb = 1
          else
            iesfb = 2
          endif
          iprkb = ipr
          kpr = kprib
          if ( eprsf(iscm) ) krcru = ircr
          do 172 krcr = 1, krcru
            if ( eprsf(iscm) ) then
              if ( ircr .eq. krcr ) then
                lrcru = jrcr
              elseif ( esfk(iscm) ) then
                lrcru = krcr
              else
                lrcru = jrcru
              endif
            else
              if ( esfk(iscm) ) lrcru = krcr
            endif
            do 170 lrcr = 1,lrcru
              if ( esfk(iscm) .and. (krcr .eq. lrcr) ) then
                kesfb = 1
              else
                kesfb = 2
              endif
              ipr = iprkb
              do 160 ist = 1, nst
                npairi = npriri(iesfb,iscm,ist)
                npairk = nprirk(kesfb,iscm,ist)
                if ( npairi .eq. 0 ) then
                  kpr = kpr + npairk
                  go to 160
                elseif ( npairk .eq. 0 ) then
                  ipr = ipr + npairi
                  go to 160
                endif
                kpri = kpr
                do 150 ipair = 1, npairi
                  ipr = ipr + 1
                  if ( eprsf(iscm) .and. (ircr .eq. krcr)
     &             .and. (jrcr .eq. lrcr) ) npairk = ipair
                  kpr = kpri
                  do 140 kpair = 1, npairk
                    kpr = kpr + 1
                    ibl = ibl + 1
                    if ( abs(h4(ibl)) .gt. cutoff ) then
                      if ( ibuf .eq. n2max ) then
                        call wtlab2( msame, buffer, values, labels,
     &                   ibitv )
                      endif
                      ibuf           = ibuf + 1
                      labels(1,ibuf) = il(ipr) / ilfact
                      labels(2,ibuf) = mod( il(ipr), ilfact )
                      labels(3,ibuf) = il(kpr) / ilfact
                      labels(4,ibuf) = mod( il(kpr), ilfact )
                      values(ibuf)   = h4(ibl)
                    endif
140               continue
150             continue
160           continue
170         continue
172       continue
174     continue
176   continue
c
      return
      end
*deck wthead
      subroutine wthead( nsym,   ns,     nbft,   ntitle, title,
     &   repnuc, ityp,   nbpsy,  mtype,  ms,     mnl,    info,
     &   aoints, nlist,  map )
c
c  write the header info to the output aoints file.
c
c  input:
c  nsym = number of symmetry blocks.
c  ns   = number of symmetry inequivalent atom labels.
c  nbft = total number of basis functions.
c  ntitle = number of titles.
c  title(1:ntitle) = titles.
c  repnuc = nuclear repulsion energy.
c  ityp(1:nsym) = input symmetry labels.
c  nbpsy(1:nsym) = number of basis functions in each symmetry block.
c  mtype(1:ns) = atom labels.
c  ms(1:nbft) = basis_function-to-center mapping vector.
c  mnl(1:nbft) = integer codes for basis function types.
c                1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g
c  info(*) = info(*) array for the output file.
c  aoints = integral output file unit number.
c  nlist = listing file output unit.
c  map(1:2,1:nbft) = temporary scratch array for writing the output
c                    map(*) vectors.
c
c  01-dec-90 SIFS version written by ron shepard.
c
      implicit logical(a-z)
c
      integer    nbfmx,     nmap,   ninfo,   nenrgy,   nrtype
      parameter( nbfmx=255, nmap=2, ninfo=5, nenrgy=1, nrtype=-1 )
c
c     # dummy:
      integer nsym, ns, nbft, ntitle, aoints, nlist
      character*80  title(1:ntitle)
      real*8        repnuc
      character*3   ityp(1:nsym)
      integer       nbpsy(1:nsym)
      character*3   mtype(1:ns)
      integer       ms(1:nbft)
      integer       mnl(1:nbft)
      integer       info(1:ninfo)
      integer       map(nbft,nmap)
c
c     # local:
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer    bmin,   bmax
      parameter( bmin=0, bmax=17 )
c
      integer      clip, ierr, i, j
      character*8  bfnlab(nbfmx)
      character*4  slabel(8)
      real*8       energy(nenrgy)
      integer      ietype(nenrgy)
      integer      imtype(nmap)
      character*2  btype(bmin:bmax)
c
      real*8   sifsce
      external sifsce
c
c     # map vector types.
      data imtype/ 3, 4 /
c
c     # basis function types.
      data btype/ '??', '_s', '_p', '3s', '_d', '4p', '_f', '5s', '5d',
     &      '_g', '6p', '6f', '_h', '7s', '7d', '7g', '_i', '??' /
c
c     # this routine is a bit of a hack.  the sifs data should
c     # be more intimately incorporated into argos. -rls
c
      if ( nbft .gt. nbfmx ) then
c       # this should eventually be handled more gracefully.
c       # the only problem is the bfnlab(*) array, which should
c       # really be a dummy argument anyway. -rls
        call bummer('wthead: (nbft-nbfmx)=', (nbft-nbfmx), faterr )
      endif
c
c     # set map(*,*) and create bfnlab(*) from the input arrays.
c
      do 20 i = 1, nbft
        map(i,1) = ms(i)
        map(i,2) = mnl(i)
        clip = min( max( bmin, mnl(i) ), bmax )
        write( bfnlab(i), fmt='(i3,a3,a2)' )
     &   i, mtype(ms(i)), btype(clip)
        do 10 j = 4, 8
c         # replace embedded spaces with a printing character to
c         # avoid confusion when printing labels later.
          if ( bfnlab(i)(j:j) .eq. ' ' ) then
            bfnlab(i)(j:j) = '_'
          endif
10      continue
20    continue
c
c     # create the output symmetry labels.
c
      do 40 i = 1, nsym
        slabel(i)(1:1) = ' '
        slabel(i)(2:4) = ityp(i)
40    continue
c
c     # the info(*) entries were passed in.
c
c     info(1) = fsplit
c     info(2) = l1rec
c     info(3) = n1max
c     info(4) = l2rec
c     info(5) = n2max
c
c     # setup for a single nuclear repulsion energy.
c
      energy(1) = repnuc
      ietype(1) = nrtype
c
c     # write the headers.
c
      call sifwh( aoints,  ntitle,  nsym,   nbft,   ninfo,   nenrgy,
     &   nmap,    title,   nbpsy,   slabel, info,   bfnlab,  ietype,
     &   energy,  imtype,  map,     ierr )
c
      if ( ierr .ne. 0 ) then
        call bummer('wthead: from sifwh(), ierr=', ierr, faterr )
      endif
c
      write(nlist,6100)'output SIFS file header information:'
      write(nlist,'(1x,a)') (title(i),i=1,ntitle)
c
      write(nlist,6100)'output energy(*) values:'
      call sifpre( 6, nenrgy, energy, ietype )
c
      write(nlist,6200)'total core energy =',
     & sifsce( nenrgy, energy, ietype )
c
      write(nlist,6040) nsym, nbft
c
      write(nlist,6050)'symmetry  =',(i,i=1,nsym)
      write(nlist,6060)'slabel(*) =',(slabel(i),i=1,nsym)
      write(nlist,6050)'nbpsy(*)  =',(nbpsy(i),i=1,nsym)
c
      write(nlist,6100)'info(*) =',(info(i),i=1,ninfo)
c
      write(nlist,6100)'output orbital labels, i:bfnlab(i)='
      write(nlist,6080)(i,bfnlab(i),i=1,nbft)
c
      write(nlist,6100)'bfn_to_center map(*), i:map(i)'
      write(nlist,6110) (i,map(i,1),i=1,nbft)
c
      write(nlist,6100)'bfn_to_orbital_type map(*), i:map(i)'
      write(nlist,6110) (i,map(i,2),i=1,nbft)
c
6040  format(/' nsym =',i2,' nbft=',i4/)
6050  format(1x,a,8i5)
6060  format(1x,a,8(1x,a4))
6080  format(10(i4,':',a8))
6100  format(/1x,a,(1x,10i10))
6110  format(10(i4,':',i3))
6200  format(/1x,a,(1x,1p4e20.12))
c
      return
      end
*deck w1zzzz
*deck w1stup
*deck wtlab1
      subroutine w1zzzz( last, fcore, buffer, values, labels )
c
c  01-dec-90 SIFS version. -rls
c  11-aug-87 written by ron shepard.
c
      implicit logical(a-z)
c
c     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
c
      integer       l1rec
      equivalence ( l1rec, info(2) )
      integer       n1max
      equivalence ( n1max, info(3) )
c
      integer    nipv
      parameter( nipv=2 )
c
c     # dummy:
      integer last
      integer labels(nipv,*)
      real*8 fcore
      real*8 buffer(*), values(*)
c     # ibm+convex f77 bug; actual:
c     # integer labels(nipv,n1max)
c     # real*8 buffer(l1rec), values(n1max)
c
c     # local:
      integer ierr
      integer ibitv(1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer    ibvtyp
      parameter( ibvtyp=0 )
c
c     # this routine is only to set up the entry points and should
c     # never be called.
      call bummer('w1zzzz called',0,faterr)
c
c**********************************************************************
      entry w1stup
c**********************************************************************
c
c     # perform initialization for the 1-e integral records.
c
      ibuf   = 0
      numout = 0
      nrec   = 0
      return
c
c**********************************************************************
      entry wtlab1( last, fcore, buffer, values, labels )
c**********************************************************************
c
c     # encode and write a 1-e integral buffer.
c
      numout = numout + ibuf
c
      call sifew1( ntape,  info,   nipv,   ibuf,   last,   itypea,
     &     itypeb, ifmt1,  ibvtyp, values, labels, fcore,  ibitv,
     &     buffer, nrec,   ierr )
c
      if ( ierr .ne. 0 ) then
        call bummer('wtlab1: from sifew1(), ierr=',ierr,faterr)
      endif
c
c     # allow for partial record writes.
      numout = numout - ibuf
c
      return
      end
*deck w2zzzz
*deck w2stup
*deck wtlab2
      subroutine w2zzzz( last, buffer, values, labels, ibitv )
c
c  01-dec-90 SIFS version. -rls
c  11-aug-87 written by ron shepard.
c
      implicit logical(a-z)
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
      integer    nipv
      parameter( nipv=4 )
c
c     # dummy:
      integer last, labels(nipv,*), ibitv(*)
      real*8 buffer(*), values(*)
c     # ibm+convex f77 bug; actual:
c     # integer labels(nipv,n2max), ibitv( ((n2max+63)/64)*64 )
c     # real*8 buffer(l2rec), values(n2max)
c
c     # local:
      integer ierr
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer    ibvtyp,   iwait
      parameter( ibvtyp=1, iwait=0 )
c
c     # reqnum is an asynchronous i/o handle.
c     # it must be saved between sifew2() and sif2w8() calls.
      integer reqnum
      save    reqnum
      data    reqnum / 0 /
c
c     # this routine is only to set up the entry points and should
c     # never be called.
      call bummer('w2zzzz called',0,faterr)
c
c**********************************************************************
      entry w2stup( ibitv )
c**********************************************************************
c
c     # perform initialization for the 2-e integral records.
c
      ibuf   = 0
      numout = 0
      nrec   = 0
c
      call izero_argos( (((n2max+63)/64)*64), ibitv, 1 )
      return
c
c**********************************************************************
      entry wtlab2( last, buffer, values, labels, ibitv )
c**********************************************************************
c
c     # encode and write a 2-e integral record.
c
      numout = numout + ibuf
c
c     # wait for i/o on the last buffer to complete.
      if ( reqnum .ne. 0 ) then
        call sif2w8( ntape, info, reqnum, ierr )
c
        if ( ierr .ne. 0 ) then
          call bummer('wtlab2: from sif2w8(), ierr=',ierr,faterr)
        endif
      endif
c
c     # encode and initiate the write of the record.
c     # [this call returns reqnum associated with the i/o request.]
c
      call sifew2( ntape,  info,   nipv,   ibuf,   last,   itypea,
     &     itypeb, ifmt2,  ibvtyp, values, labels, ibitv,  buffer,
     &     iwait,  nrec,   reqnum, ierr )
c
      if ( ierr .ne. 0 ) then
        call bummer('wtlab2: from sifew2(), ierr=',ierr,faterr)
      endif
c
c     # reset numout and ibitv(*), allowing for partial record writes.
      numout = numout - ibuf
c
      call izero_argos( (n2max-ibuf), ibitv(ibuf+1), 1 )
c
      return
      end
