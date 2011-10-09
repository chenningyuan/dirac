*deck gcentr
      subroutine gcentr(ica,nc)
c
c  generate group of center-interchange operators
c
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
c
      dimension ica(mcu,msu,*), nc(*)
c
      if (ng .eq. 1) go to 56
c
c  set up identity operator
c
      do 16 is=1,ns
        icu=nc(is)
        if(icu.eq.1) go to 16
        do 12 ic=1,icu
          ica(ic,is,1)=ic
   12   continue
   16 continue
c
c        do 52(+3) ig=2,ng
c
      ig=2
   20 do 52 jg=2,ig
      igp=ig
      jgp=jg
c
c  generate new operator candidate
c
   24 do 32 is=1,ns
        icu=nc(is)
        if(icu.eq.1) go to 32
        do 28 ic=1,icu
          ica(ic,is,ng+1)=ica(ica(ic,is,igp),is,jgp)
   28   continue
   32 continue
c
c  check against list of operators
c
      do 44 kg=1,ng
        do 40 is=1,ns
          icu=nc(is)
          if(icu.eq.1) go to 40
          do 36 ic=1,icu
            if(ica(ic,is,ng+1).ne.ica(ic,is,kg)) go to 44
   36     continue
   40   continue
c       # not a new operator
        go to 48
   44 continue
c     # new operator found
      ng=ng+1
   48 if(jgp.eq.ig) go to 52
c
c  take product in other order
c
      igp=jg
      jgp=ig
      go to 24
   52 continue
      if(ig.eq.ng) go to 56
      ig=ig+1
      go to 20
c
c  set up operators for unique nuclei
c
   56 do 64 is=1,ns
      if(nc(is).eq.1) then
        do 60 ig=1,ng
          ica(1,is,ig)=1
   60   continue
      endif
   64 continue
      return
      end
*deck stvint
      subroutine stvint(hpt,hwt)
c
c  gauss-hermite quadrature using minimum-point formula
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0)
      common /stv/ xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension hpt(*), hwt(*)
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj)/2
      imax=(npts*(npts+1))/2
      imin=imax-npts+1
      do 40 i=imin,imax
        dum=hpt(i)/t
        ptx=dum+x0
        pty=dum+y0
        ptz=dum+z0
        ax=ptx-xi
        ay=pty-yi
        az=ptz-zi
        bx=ptx-xj
        by=pty-yj
        bz=ptz-zj
        dum=hwt(i)
c       # take care to avoid (0.0d0)**0
        dumx = dum
        dumy = dum
        dumz = dum
        do 20 ipwr = 1,ni-1
          dumx = dumx*ax
          dumy = dumy*ay
          dumz = dumz*az
   20   continue
        do 30 ipwr = 1,nj-1
          dumx = dumx*bx
          dumy = dumy*by
          dumz = dumz*bz
   30   continue
        xint = xint + dumx
        yint = yint + dumy
        zint = zint + dumz
c       xint=xint+dum*ax**(ni-1)*bx**(nj-1)
c       yint=yint+dum*ay**(ni-1)*by**(nj-1)
c       zint=zint+dum*az**(ni-1)*bz**(nj-1)
   40 continue
      return
      end


*deck qcomp
      function qcomp(alpha,dfac,n,l,t,xk)
c
c  compute q(n,l) scaled by sqrt(pi)*exp(-t) to prevent overflows
c  arguments are alpha, xk, and t=xk**2/(4*alpha)
c  no restriction on the magnitude of t
c  increase dfac array to raise n, l restrictions
c
      implicit real*8 (a-h,o-z)
      parameter (am1=-1.0d0, a0=0.0d0, accpow=1.0d-14,
     1  accasy=1.0d-10, a1rtpi=0.56418958354775629d0, a1=1.0d0,
     2  a2=2.0d0, a4=4.0d0)
      dimension dfac(*)
      dimension tmin(0:8)
      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1  18.0d0,15.0d0/
c
      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c
c  use alternating series (n+l.le.22.and.l.le.10)
c
      if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=xkp * dfac(n+l+1)/
     1  ((alpha+alpha)**((n-l)/2) * sqrt(a4*alpha) * dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
   10 if(num.ne.0) then
        fnum=num
        term=term*fnum*t/(xden*xc)
        xc=xc+am1
        sum=sum+term
        num=num+2
        xden=xden+a2
        go to 10
      endif
c
      qcomp=prefac*sum
      return
c
   30 if(t.lt.tmin(min(n,8))) go to 60
c
c  use asymptotic series (arbitrary n,l)
c
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp/((alpha+alpha)*sqrt(a4*alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(abs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
c
c  use power series (n+l.le.22.and.l.le.10)
c
   60 if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) then
        prefac=prefac/sqrt(a4*alpha)
      else
        prefac=a1rtpi*prefac
      endif
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end
