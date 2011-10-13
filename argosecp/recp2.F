!     subroutine l2opxv(lmn,v,itl,itu,lmnv,eval,leig)
      subroutine l2opxv(lmn,itl,itu,lmnv,eval,leig,ics,igcs)
c  compute the matrix vector product, (L^2) * v, in an unnormalized
c  cartesian basis, and determine if v(*) is an eigenvector of the
c  total angular momentum operator.
c
c  input:
c  lmn = (l + m + n) where l, m, and n are the exponents of x, y, and z
c        respectively.  lmn=0 for cartesian "s" functions, 1 for
c        cartesian "p" functions, 2 for cartesian "d" functions, etc.
c  v(1:xyzdim) = input vector.  the elements are coefficients of the
c              cartesian basis functions (x**l * y**m * z**n).
c              where xyzdim = ((lmn + 1) * (lmn + 2))/2
c  lmnv(1:3,1:lmnvmx) = cartesian exponents of each basis function.
c  scr(1:lenscr) = scratch work vector.  This must be at least as large
c                  as the cartesian subspace dimension.
c                  lenscr >= xyzdim
c  output:
c  eval = vector expectation value = <v| L^2 |v> / <v|v>
c  leig = eigenvector return code.
c       = -4 for vector norm error.
c       = -3 for lmnv(*,*) inconsistency.
c       = -1 if v(*) is not an eigenvector of L^2
c       = principal quantum number if v(*) is an eigenvector.
c         0 for s-type vectors, 1 for p-type vectors, etc.
c         the eigenvalue is given by (leig*(leig+1)) to within numerical
c         precision.
c  28-oct-96 account for nonorthogonality of cartesian gaussians (rmp)
c  11-jun-91 written by ron shepard with suggestions from r. m. pitzer.
 
      USE RECP_IPT
      implicit logical(a-z)
 
      integer itl, itu, lmn, leig, lmnv(3,*)
!     real*8 eval, v(*)
      real*8 eval
      integer ics,igcs 
c     # local:
      integer l, m, n, i, j, it, jt, xyzdim, ll, mm, nn, p, q
      real*8 vi, rnorm, vnorm2
 
c     # small is used to determine if a vector is close enough to be
c     #       called an eigenvector.  this should be about 1000x the
c     #       machine epsilon to avoid sqrt() precision problems.
      real*8    zero,     one,     two,     fourth,        small
      parameter(zero=0d0, one=1d0, two=2d0, fourth=0.25d0, small=1d-10)
 
#include "argos_error.h"
 
      integer  iodfac
      external iodfac
 
c     # the [l,m,n] basis function is assigned the internal index:
c     #  (lx*(lx+1))/2 + n + 1   with lx=(lmn-l)=(m+n)
c     # this convention is inconsistent with the rest of this program,
c     # but it is used anyway since it is convenient and since the
c     # computed results are local.
      integer llp1
      llp1(l) = (l * (l + 1)) / 2 + 1
 
      leig = -1
      eval = -one
 
c     # check to make sure scr(*) is large enough.
      xyzdim = llp1( lmn) + lmn
      ALLOCATE(IPT_L2OPXV_SCR(xyzdim))

c     # initialize the scratch vector.
      do i = 1, xyzdim
         IPT_L2OPXV_SCR(i) = 0.0d0
      enddo 
 
c     # compute the matrix-vector product and the vector norm.
c     # matrix elements of lop(*,*) are computed from the operator
c     # definition in cartesian space.
      vnorm2 = zero
      do 30 it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
!        vi = v(i)
         vi = IPT_AO2SO(i,ics+1,igcs)
         do jt = itl, it-1
            j = jt - (itl - 1)
!           vnorm2 = vnorm2 + (two * vi) * v(j) *
            vnorm2 = vnorm2 + (two * vi) 
     &              * IPT_AO2SO(j,ics+1,igcs) *
     &       iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
!           print *,'v(',i,')',vi,'v(',j,')',v(j)
         enddo
         vnorm2 = vnorm2 + vi * vi * iodfac(l,l,m,m,n,n)
 
         if ( (l+m+n) .ne. lmn ) then
!           # inconsistent exponents in the basis function.
            call bummer('l2opxv: (l+m+n-lmn)=', (l+m+n-lmn), wrnerr )
            leig = -3
            return
         endif
 
!        # L^2 is sparse in this representation.  each vi contributes
!        # to, at most, only 7 elements in the matrix-vector product.

         if ( l .ge. 2 ) then
!           # include -l*(l-1) * ( [l-2,m,n+2] + [l-2,m+2,n] ) terms.
            ll = l * (l - 1)
            p = llp1( lmn - l + 2 ) + n
            IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (ll)
            p = p + 2
            IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (ll)
         endif
 
         if ( m .ge. 2 ) then
!           # include -m*(m-1) * ( [l,m-2,n+2] + [l+2,m-2,n] ) terms.
            mm = m * (m - 1)
            p = llp1( lmn - l ) + n + 2
            IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (mm)
            p = llp1( lmn - l - 2 ) + n
            IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (mm)
         endif
 
         if ( n .ge. 2 ) then
!           # include -n*(n-1) * ( [l,m+2,n-2] + [l+2,m,n-2] ) terms.
            nn = n * (n - 1)
            p = llp1( lmn - l ) + n - 2
            IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (nn)
            p = llp1( lmn - l - 2 ) + n - 2
            IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) - vi * (nn)
         endif
 
!        # include the 2*(l*m+l*n+m*n+l+m+n)*[l,m,n] diagonal term.
         p = llp1( lmn - l ) + n
         IPT_L2OPXV_SCR(p) = IPT_L2OPXV_SCR(p) 
     &       + vi * (2 * (l * (m + n) + m * n + l + m + n))
   30 continue
 
      if ( vnorm2 .le. small ) then
         call bummer('l2mxv: small vector norm, lnm=', lmn, wrnerr)
         leig = -4
         return
      endif
 
!     # compute the expectation value.
      eval  = zero
      do it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         p = llp1( lmn - l ) + n
         do jt = itl, itu
            j = jt - (itl - 1)
!           eval = eval + v(j) * IPT_L2OPXV_SCR(p) *
            eval = eval + IPT_AO2SO(j,ics+1,igcs) 
     &                   * IPT_L2OPXV_SCR(p) *
     &                  iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
         enddo
      enddo
      eval = eval / vnorm2
 
!     # compute the residual norm.
      rnorm = zero
 
      do it = itl, itu
         i = it - (itl - 1)
         l = lmnv(1,it)
         m = lmnv(2,it)
         n = lmnv(3,it)
         p = llp1( lmn - l ) + n
         do jt = itl, it-1
            j = jt - (itl - 1)
            q = llp1( lmn - lmnv(1,jt) ) + lmnv(3,jt)
!           rnorm = rnorm + two * ( IPT_L2OPXV_SCR(p) - eval * v(i) ) *
            rnorm = rnorm + two * ( IPT_L2OPXV_SCR(p) 
     &                             - eval * IPT_AO2SO(i,ics+1,igcs) ) *
!    &                            ( IPT_L2OPXV_SCR(q) - eval * v(j) ) *
     &                            ( IPT_L2OPXV_SCR(q) 
     &                             - eval * IPT_AO2SO(j,ics+1,igcs) ) *
     &                iodfac(l,lmnv(1,jt),m,lmnv(2,jt),n,lmnv(3,jt))
         enddo
!       rnorm = rnorm + (IPT_L2OPXV_SCR(p) - eval * v(i))**2 
        rnorm = rnorm + (IPT_L2OPXV_SCR(p) 
     &                    - eval * IPT_AO2SO(i,ics+1,igcs))**2 
     &                  * iodfac(l,l,m,m,n,n)
      enddo
 
!     # normalize w.r.t. |v|=1.
      rnorm = rnorm / vnorm2
      rnorm = sqrt( rnorm )

      if ( rnorm .gt. small ) then
!        # v(*) is not an eigenvector.
         leig = -1
      else
!        # v(*) is an eigenvector.
!        # determine leig such that eval = leig*(leig+1)
!        # the following assignment should truncate to the
!        # next smaller integer value...
         leig = ( sqrt( eval + fourth ) )
      endif
      DEALLOCATE(IPT_L2OPXV_SCR)
      return
      end



*deck hermit
      subroutine hermit(nn,x,a,eps)
c          calculates the zeros  x(i)  of the nn-th order
c        hermite polynomial.  the largest zero will be
c        stored in x(1).  also calculates the corresponding
c        coefficients  a(i)  of the nn-th order gauss-hermite
c        quadrature formula of degree 2*nn-1.  the factor of
c        sqrt(pi) has been removed from the a(i).
c
c  a. h. stroud & d. secrest, gaussian quadrature formulas,
c  prentice-hall, 1966
c
      implicit real*8 (a-h,o-z)
      parameter (sixth = 1.0d0/6.0d0)
      dimension x(*), a(*)

      IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : HERMIT ENTER * '
c
      fn = (nn)
      n1 = nn - 1
      n2 = (nn+1)/2
      cc = 1.0d0
      s  = 0.0d0
      do 15 i=1,n1
        s  = s + 0.5d0
        cc = s*cc
   15 continue
      s  = (2.0d0*fn+1.0d0)**sixth
      do 30 i=1,n2
        if( i.eq.1 ) then
c         # largest zero
          xt = s**3 - 1.85575d0/s
        elseif( i.eq.2 ) then
c         # second zero
          xt = xt - 1.14d0*fn**0.426d0/xt
        elseif( i.eq.3 ) then
c         # third zero
          xt = 1.86d0*xt - 0.86d0*x(1)
        elseif( i.eq.4 ) then
c         # fourth zero
          xt = 1.91d0*xt - 0.91d0*x(2)
        else
c         # all other zeros
          xt = 2.0d0*xt - x(i-2)
        endif
c
        call hroot(xt,nn,dpn,pn1,eps)
        x(i) = xt
        a(i) = cc/dpn/pn1
c       write (6,'(2i4,2d25.17)') nn, i, xt, a(i)
        ni = nn-i+1
        x(ni) = -xt
        a(ni) = a(i)
   30 continue

      IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : HERMIT EXIT  * '

      return
      end
*deck hroot
      subroutine hroot(x,nn,dpn,pn1,eps)
c     improves the approximate root  x
c     in addition we also obtain
c         dpn = derivative of h(n) at x
c         pn1 = value of h(n-1) at x

      implicit real*8 (a-h,o-z)

c     # iter = 5 sufficient for 8-byte accuracy up to nn = 7
      do 14 iter=1,10
        call hrecur_argos(p,dp,pn1,x,nn)
        d  = p/dp
        x  = x - d
        if( abs(d).le.eps ) go to 16
   14 continue
   16 dpn = dp
      return
      end

*deck hrecur_argos
      subroutine hrecur_argos(pn,dpn,pn1,x,nn)
c
      implicit real*8 (a-h,o-z)
c
      p1 = 1.0d0
      p  = x
      dp1 = 0.0d0
      dp  = 1.0d0
      do 20 j=2,nn
        fj = (j)
        fj2 = 0.5d0*(fj-1.0d0)
        q  = x*p - fj2*p1
        dq = x*dp + p - fj2*dp1
        p1 = p
        p  = q
        dp1 = dp
        dp  = dq
   20 continue
      pn  = p
      dpn = dp
      pn1 = p1
      return
      end

      subroutine cortab(eps,lmax,lmn1u,lproju,ndfac)
c     # tables for core potential and spin-orbit integrals.
      USE RECP_IPT
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0, a3=3.0d0)

c     # compute gauss-hermite points and weights for c, z integrals.
      igh =1
      nn = 5
      do 120 i = 1, 3
        call hermit(nn,IPT_SH_HPT(igh),IPT_SH_HWT(igh),eps)
        igh = igh + nn
        nn = 2*nn
  120 continue
c     # compute double factorials.
      IPT_SH_DFAC(1)=a1
      IPT_SH_DFAC(2)=a1
      fi=a0
      do 130 i=1,ndfac-2
        fi=fi+a1
        IPT_SH_DFAC(i+2)=fi * IPT_SH_DFAC(i)
  130 continue
c     # compute binomial coefficients.
      inew=1
      IPT_SH_BINOM(1)=a1
      do 150 j=1,lmn1u-1
        inew=inew+1
        IPT_SH_BINOM(inew)=a1
        do 140 i=1,j-1
          inew=inew+1
          IPT_SH_BINOM(inew)=((j-i+1)*IPT_SH_BINOM(inew-1))/(i)
  140   continue
        inew=inew+1
        IPT_SH_BINOM(inew)=a1
  150 continue
c     # compute tables by recursion for real spherical harmonics.  they
c     # are indexed by l, m and sigma.  the sequence number of the
c     # harmonic with quantum numbers l, m and sigma is given by
c     #            l**2+2*m+1-sigma
c     # IPT_SH_LMF(index) and IPT_SH_LML(index) hold the positions 
c     # of the first and last terms of the harmonic in the 
c     # IPT_SH_LMX, IPT_SH_LMY, IPT_SH_LMZ, and IPT_SH_ZLM arrays.
c     # the harmonics with angular momentum l are generated from those
c     # with angular momenta l-1 and l-2.
c     # for m = 0,1,2,...,l-1, the recursion relation
c     z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
c                  sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
c     # is used.
c     # for m = l, the recursion relation
c     x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) =
c                  sqrt((2*l))/((2*l+1)))*Z(l,l,s)
c     # is used.
c     # l=0
      IPT_SH_LMF(1) = 1
      IPT_SH_LML(1) = 1
      IPT_SH_LMX(1) = 0
      IPT_SH_LMY(1) = 0
      IPT_SH_LMZ(1) = 0
      IPT_SH_ZLM(1) = a1
c     # l=1
      IPT_SH_LMF(2) = 2
      IPT_SH_LML(2) = 2
      IPT_SH_LMX(2) = 0
      IPT_SH_LMY(2) = 0
      IPT_SH_LMZ(2) = 1
      IPT_SH_ZLM(2) = sqrt(a3)
      IPT_SH_LMF(3) = 3
      IPT_SH_LML(3) = 3
      IPT_SH_LMX(3) = 0
      IPT_SH_LMY(3) = 1
      IPT_SH_LMZ(3) = 0
      IPT_SH_ZLM(3) = IPT_SH_ZLM(2)
      IPT_SH_LMF(4) = 4
      IPT_SH_LML(4) = 4
      IPT_SH_LMX(4) = 1
      IPT_SH_LMY(4) = 0
      IPT_SH_LMZ(4) = 0
      IPT_SH_ZLM(4) = IPT_SH_ZLM(2)
      nterm=4
      do 270 lang=2,lmax
        do 240 mang=0,lang-1
          anum = ((2*lang-1)*(2*lang+1))
          aden = ((lang-mang)*(lang+mang))
          coef1 = sqrt(anum/aden)
          anum = ((lang+mang-1)*(lang-mang-1)*(2*lang+1))
          aden = (2*lang-3)*aden
          coef2 = sqrt(anum/aden)
          nsigma=min(1,mang)
          do 230 isigma=nsigma,0,-1
            indexh=lang**2+2*mang+1-isigma
            lone=lang-1
            ltwo=lang-2
            ione=lone**2+2*mang+1-isigma
            itwo=ltwo**2+2*mang+1-isigma
            IPT_SH_LMF(indexh)=IPT_SH_LML(indexh-1)+1
            IPT_SH_LML(indexh)=IPT_SH_LML(indexh-1)
            nxy=(mang-isigma+2)/2
            iu=IPT_SH_LMF(ione)+nxy-1
            do 200 i=IPT_SH_LMF(ione),iu
              IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
              j=IPT_SH_LML(indexh)
              IPT_SH_LMX(j)=IPT_SH_LMX(i)
              IPT_SH_LMY(j)=IPT_SH_LMY(i)
              IPT_SH_LMZ(j)=IPT_SH_LMZ(i)+1
              IPT_SH_ZLM(j)=IPT_SH_ZLM(i)*coef1
              nterm=nterm+1
  200       end do
            if(ltwo.ge.mang) then
              il=iu+1
              do 210 i=il,IPT_SH_LML(ione)
                IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
                j=IPT_SH_LML(indexh)
                k=IPT_SH_LMF(itwo)+i-il
                IPT_SH_LMX(j)=IPT_SH_LMX(k)
                IPT_SH_LMY(j)=IPT_SH_LMY(k)
                IPT_SH_LMZ(j)=IPT_SH_LMZ(k)
                IPT_SH_ZLM(j)=IPT_SH_ZLM(i)*coef1-IPT_SH_ZLM(k)*coef2
                nterm=nterm+1
  210         enddo
              il=IPT_SH_LML(itwo)-nxy+1
              if(mod(lang-mang,2).eq.0) then
                do 220 i=il,IPT_SH_LML(itwo)
                  IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
                  j=IPT_SH_LML(indexh)
                  IPT_SH_LMX(j)=IPT_SH_LMX(i)
                  IPT_SH_LMY(j)=IPT_SH_LMY(i)
                  IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
                  IPT_SH_ZLM(j)=-IPT_SH_ZLM(i)*coef2
                  nterm=nterm+1
  220           end do
              endif
            endif
  230     end do
  240   end do
        anum = (2*lang+1)
        aden = (2*lang)
        coef = sqrt(anum/aden)
        mang=lang
        isigma=1
        indexh=lang**2+2*mang+1-isigma
        IPT_SH_LMF(indexh)=IPT_SH_LML(indexh-1)+1
        IPT_SH_LML(indexh)=IPT_SH_LML(indexh-1)
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonic (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        k=IPT_SH_LMF(isigm)
        do 250 i=IPT_SH_LMF(isig),IPT_SH_LML(isig)
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          IPT_SH_LMX(j)=IPT_SH_LMX(i)+1
          IPT_SH_LMY(j)=IPT_SH_LMY(i)
          IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
          IPT_SH_ZLM(j)=(IPT_SH_ZLM(i)+IPT_SH_ZLM(k))*coef
          k=k+1
          nterm=nterm+1
  250   enddo
        if(mod(mang,2).eq.1) then
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          IPT_SH_LMX(j)=IPT_SH_LMX(k)
          IPT_SH_LMY(j)=IPT_SH_LMY(k)+1
          IPT_SH_LMZ(j)=IPT_SH_LMZ(k)
          IPT_SH_ZLM(j)=IPT_SH_ZLM(k)*coef
          nterm=nterm+1
        endif
        isigma=0
        indexh=lang**2+2*mang+1-isigma
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonuc (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        IPT_SH_LMF(indexh)=IPT_SH_LML(indexh-1)+1
        IPT_SH_LML(indexh)=IPT_SH_LMF(indexh)
        j=IPT_SH_LML(indexh)
        i=IPT_SH_LMF(isig)
        IPT_SH_LMX(j)=IPT_SH_LMX(i)+1
        IPT_SH_LMY(j)=IPT_SH_LMY(i)
        IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
        IPT_SH_ZLM(j)=IPT_SH_ZLM(i)*coef
        nterm=nterm+1
        k=IPT_SH_LMF(isigm)
        do 260 i=IPT_SH_LMF(isig)+1,IPT_SH_LML(isig)
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          IPT_SH_LMX(j)=IPT_SH_LMX(i)+1
          IPT_SH_LMY(j)=IPT_SH_LMY(i)
          IPT_SH_LMZ(j)=IPT_SH_LMZ(i)
          IPT_SH_ZLM(j)=(IPT_SH_ZLM(i)-IPT_SH_ZLM(k))*coef
          k=k+1
          nterm=nterm+1
  260   enddo
        if(mod(mang,2).eq.0) then
          IPT_SH_LML(indexh)=IPT_SH_LML(indexh)+1
          j=IPT_SH_LML(indexh)
          k=IPT_SH_LML(isigm)
          IPT_SH_LMX(j)=IPT_SH_LMX(k)
          IPT_SH_LMY(j)=IPT_SH_LMY(k)+1
          IPT_SH_LMZ(j)=IPT_SH_LMZ(k)
          IPT_SH_ZLM(j)=-IPT_SH_ZLM(k)*coef
          nterm=nterm+1
        endif
  270 end do
      ixy = 0
      iz = 0
      do 300 lang=1,lproju
        do 290 mang=0,lang-1
          nsigma=min(1,mang)
          ndelta=max(0,1-mang)
          anum = ((lang-mang)*(lang+mang+1))
          aden = (2*(2-ndelta))
          coef=sqrt(anum/aden)
          do 280 isigma=nsigma,0,-1
            isign=2*isigma-1
            ixy = ixy+1
            IPT_SH_FLMTX(1,ixy) = (isign)*coef
            IPT_SH_FLMTX(2,ixy) = coef
            if(mang.ne.0) then
              iz=iz+1
              IPT_SH_FLMTX(3,iz) = -(mang*isigma)
            endif
  280     enddo
  290   enddo
        iz=iz+1
        IPT_SH_FLMTX(3,iz) = -(lang)
  300 enddo
c     # column and row indices for angular momentum matrix elements.
      iadd = 1
      do 310 i=1,2*lproju-1
        IPT_SH_MC(1,i) = i
        IPT_SH_MC(2,i) = i
        IPT_SH_MC(3,i) = i+1
        IPT_SH_MR(1,i) = i+iadd
        IPT_SH_MR(2,i) = i+2
        IPT_SH_MR(3,i) = i+2
        iadd = 4-iadd
  310 continue
      return
      end

! iodfac -> FN_IODFAC 
      function iodfac(l1,l2,m1,m2,n1,n2)
c
c  compute the product of odd number factorials which gives the
c  one-center overlaps of cartesian gaussian aos.
c
      if(mod(l1+l2,2).eq.1 .or.
     &   mod(m1+m2,2).eq.1 .or.
     &   mod(n1+n2,2).eq.1) then
        iodfac = 0
      else
        iprd = 1
        do 20 i=3,l1+l2-1,2
          iprd = i*iprd
   20   continue
        do 30 i=3,m1+m2-1,2
          iprd = i*iprd
   30   continue
        do 40 i=3,n1+n2-1,2
          iprd = i*iprd
   40   continue
        iodfac = iprd
      endif
      return
      end
