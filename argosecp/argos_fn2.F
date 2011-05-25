*deck droot
      subroutine droot(nroots,n1,xx,ff,c,s,a,rt,r,w,uf,wf)
c
c  ith root of the jth rys polynomial is returned in r(i,j)
c  with the corresponding weight factor in w(i,j).  j=1,2,...,n
c  this version uses christoffel formula for weights.
c
*mdc*if real16 vax ibm cray
*      implicit real*16 (a-h,o-z)
*      real*8 xx,uf,wf
*mdc*else
      implicit real*8 (a-h,o-z)
      real*8 xx,uf,wf
*mdc*endif
*mdc*if real16q ibm vax rs6000 fujitsu
*c *** this mdc block is not necessary since the constants will be ***
*c *** correctly converted at compile time.  should it be removed? ***
*      parameter (a0=0.0q0, a1s2=0.5q0, a1=1.0q0, a4=4.0q0)
*mdc*else
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a4=4.0d0)
*mdc*endif
c     common /ffm/ x,nn,n1,k
c     dimension uf(*), wf(*)
c     dimension ff(19), c(10,10), s(10,10), a(10), rt(10),
c    &  r(9,9), w(9,9)
      dimension ff(*), c(n1,*), s(n1,*), a(*), rt(*),
     &  r(nroots,*), w(nroots,*), uf(*), wf(*)
c
      n=max(nroots,2)
      x=xx
c     n1=n+1
      nn=n+n
      call dfunc(nn,x,ff,c)
      do 15 i=1,n1
        do 10 j=1,n1
          s(i,j)=ff(i+j-1)
   10   continue
   15 continue
      call dsmit(n1,c,s,a,rt)
      do 25 i=1,n
        do 20 j=1,i
          w(i,j)=a0
          r(i,j)=a0
   20   continue
   25 continue
      wsum=ff(1)
      w(1,1)=wsum
      r(1,1)=ff(2)/wsum
      dum=sqrt(c(2,3)**2-a4*c(1,3)*c(3,3))
      r(1,2)=a1s2*(-c(2,3)-dum)/c(3,3)
      r(2,2)=a1s2*(-c(2,3)+dum)/c(3,3)
      if(n.eq.2) go to 50
      do 30 i=3,n1
        rt(i)=a1
   30 continue
      rt(1)=r(1,2)
      rt(2)=r(2,2)
      do 45 k=3,n
        k1=k+1
        do 35 i=1,k1
          a(i)=c(i,k1)
   35   continue
        call dnode(k,a,rt)
        do 40 i=1,k
          r(i,k)=rt(i)
   40   continue
   45 continue
   50 continue
      do 70 k=2,n
        do 65 i=1,k
          root=r(i,k)
          dum=a1/ff(1)
          do 60 j=1,k-1
            j1=j+1
            poly=c(j1,j1)
            do 55 m=1,j
              poly=poly*root+c(j1-m,j1)
   55       continue
            dum=dum+poly*poly
   60     continue
          w(i,k)=a1/dum
   65   continue
   70 continue
      do 75 k=1,nroots
        dum=r(k,nroots)
        uf(k)=dum/(a1-dum)
        wf(k)=w(k,nroots)
   75 continue
      return
      end
*deck dfunc
      subroutine dfunc(nn,x,ff,term)
*mdc*if real16 vax cray ibm
*      implicit real*16 (a-h,o-z)
*mdc*else
      implicit real*8 (a-h,o-z)
*mdc*endif
*mdc*if real16q ibm vax rs6000 fujitsu
*      parameter (a0=0.0q0, a1s2=0.5q0, a1=1.0q0 )
*      parameter (pie4=0.7853981633974483096156608 q0)
*mdc*else
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0)
      parameter (pie4=0.7853981633974483096156608 d0)
*mdc*endif
*mdc*if fps164
*      parameter (xsw=19.0d0)
*mdc*elseif cray
*      parameter (xsw=33.0d0)
*mdc*elseif real16q vax ibm rs6000 fujitsu
*      parameter (xsw=39.0q0)
*mdc*else
      parameter (xsw=15.0d0)
*mdc*endif
c     common /ffm/ x,nn,n1,k1
      dimension ff(*), term(*)
c
      e=a1s2*exp(-x)
      fac0=nn
      fac0=fac0+a1s2
      if(x.gt.xsw) go to 50
c
c  use power series
c
   10 fac=fac0
      term0=e/fac
      sum=term0
      ku=(x-fac0)
c     # sum increasing terms forwards
      do 20 k=1,ku
        fac=fac+a1
        term0=term0*x/fac
        sum=sum+term0
   20 continue
      i=1
      fac=fac+a1
      term(1)=term0*x/fac
      suma=sum+term(1)
      if(sum.eq.suma) go to 85
   40 i=i+1
      fac=fac+a1
      term(i)=term(i-1)*x/fac
      sum1=suma
      suma=suma+term(i)
      if(sum1-suma) 40,85,40
c
c  use asymptotic series
c
   50 sum=sqrt(pie4/x)
      if(nn.eq.0) go to 70
      fac=-a1s2
      do 60 k=1,nn
      fac=fac+a1
   60 sum=sum*fac/x
   70 i=1
      term(1)=-e/x
      suma=sum+term(1)
      if(sum.eq.suma) go to 85
      fac=fac0
      ku=(x+fac0-a1)
      do 80 i=2,ku
        fac=fac-a1
        term(i)=term(i-1)*fac/x
        sum1=suma
        suma=suma+term(i)
        if(sum1.eq.suma) go to 85
   80 continue
c     # xsw set too low. use power series.
      go to 10
c
c  sum decreasing terms backwards
c
   85 sum1=a0
c
c  11-feb-90 cray vectorization suppressed to avoid a numeric exception
c            which occured for x=5670. - rick ross/rls
cdir$ novector
      do 90 k=1,i
        sum1=sum1+term(i+1-k)
   90 continue
cdir$ vector
      ff(nn+1)=sum+sum1
c
c  use recurrence relation
c
      if(nn.eq.0) return
      do 95 k=1,nn
        fac0=fac0-a1
        ff(nn+1-k)=(e+x*ff(nn+2-k))/fac0
   95 continue
      return
      end
*deck dsmit
      subroutine dsmit(n1,c,s,v,y)
c
c  return an n1 by n1 triangular matrix c such that
c  c(t)sc=i, where i is an n1 by n1 identity matrix.
c
*mdc*if real16 vax ibm cray
*      implicit real*16 (a-h,o-z)
*mdc*else
      implicit real*8 (a-h,o-z)
*mdc*endif
*mdc*if real16q ibm vax
*c *** this mdc block is not necessary since the constants will be ***
*c *** correctly converted at compile time.  should it be removed? ***
*      parameter (a0=0.0q0, a1=1.0q0)
*mdc*else
      parameter (a0=0.0d0, a1=1.0d0)
*mdc*endif
c     common /ffm/ x,nn,n1,k1
      dimension c(n1,*), s(n1,*), v(*), y(*)
c
      do 20 i=1,n1
        do 10 j=1,i
          c(i,j)=a0
   10   continue
   20 continue
      do 80 j=1,n1
        kmax=j-1
        fac=s(j,j)
        do 30 k=1,kmax
          v(k)=a0
          y(k)=s(k,j)
   30   continue
        do 60 k=1,kmax
          dot=a0
          do 40 i=1,k
            dot=c(i,k)*y(i)+dot
   40     continue
          do 50 i=1,k
            v(i)=v(i)-dot*c(i,k)
   50     continue
          fac = fac - dot * dot
   60   continue
        fac=a1/sqrt(fac)
        c(j,j)=fac
        do 70 k=1,kmax
          c(k,j)=fac*v(k)
   70   continue
   80 continue
      return
      end
*deck dnode
      subroutine dnode(k,a,rt)
c
c  return in rt(i) the ith root of a polynomial
c  of order k whose mth coefficient is stored in a(m+1).
c  the initial values in rt must bracket the final values.
c
*mdc*if real16 vax cray ibm
*      implicit real*16 (a-h,o-z)
*mdc*else
      implicit real*8 (a-h,o-z)
*mdc*endif
*mdc*if real16q ibm vax
*c *** this mdc block is not necessary since the constants will be ***
*c *** correctly converted at compile time.  should it be removed? ***
*      parameter (a0=0.0q0, a1s16=6.25q-2, a1s4=0.25q0, a3s4=0.75q0)
*mdc*else
      parameter (a0=0.0d0, a1s16=6.25d-2, a1s4=0.25d0, a3s4=0.75d0)
*mdc*endif
*mdc*if cray
*      parameter (accrt=1.0d-21)
*mdc*elseif real16q vax ibm rs6000 fujitsu
*      parameter (accrt=1.0q-21)
*mdc*else
      parameter (accrt=1.0d-12)
*mdc*endif
c
      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
c
c     common /ffm/ x,nn,n1,k
      dimension a(*), rt(*)
c
c     # bummer error types.
      integer    wrnerr,  nfterr,  faterr
      parameter (wrnerr=0,nfterr=1,faterr=2)
c
      k1=k+1
      r2=a0
      p2=a(1)
      do 95 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
        p2=p2*r2+a(k1-i)
   10 continue
      prod=p1*p2
      if(prod.ge.a0) then
        write (nlist,15) m,k
   15   format(/' dnode: root number',i4,
     +    ' was not found for polynomial of order',i4//)
        call bummer('dnode: root not found, m=',m,faterr)
      endif
c
      r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      dr = a1s16 * dr
      r5=r-dr
      r6=r+dr
      if(abs(delta).lt.accrt.or.r5.eq.r.or.r6.eq.r) go to 90
      if(r5.lt.r3) r5=r3
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
        p5=p5*r5+a(k1-i)
        p6=p6*r6+a(k1-i)
   40 continue
   45 prod=p5*p6
      if(prod.lt.a0) go to 30
      prod=p3*p5
      if(prod.gt.a0) go to 60
      r5=a1s4*r3+a3s4*r5
      p5=a(k1)
      do 50 i=1,k
        p5=p5*r5+a(k1-i)
   50 continue
      go to 45
   60 r6=a1s4*r4+a3s4*r6
      p6=a(k1)
      do 70 i=1,k
        p6=p6*r6+a(k1-i)
   70 continue
      go to 45
   90 rt(m)=r
   95 continue
      return
      end
