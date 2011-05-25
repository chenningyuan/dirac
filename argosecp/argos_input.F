*deck syminp_argos
      subroutine argos_input( a,      c,          ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,       la,     lb,     lcr,
     &   lls,    lmnp1,  lmnv,   lmn1u,      lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,         mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,        nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,      nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,          y,      z,      zcr,
     &   zet,    lcore )
 
!  read and verify the user input: symmetry info, basis sets,
!  geometry, etc.
 
!  11-sep-91 parallel code added. -rls/rjh
!  12-jun-91 l2opxv() call added. -rls
!  01-mar-91 aoint2,fsplit,l1rec,l2rec input order changed. -rls/rmp
 
      implicit logical(a-z)

#include "argos_files.h" 
#include "argos_error.h" 
#include "argos_print.h" 

!     # /bufout/ holds some output integral file parameters.
      integer
     & info,     ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
      common /bufout/
     & info(10), ifmt1,   ifmt2,   itypea,  itypeb,
     & ibuf,     numout,  nrec,    ntape
 
      integer         only1e
      common /conly1/ only1e
 
      real*8         cutoff, tol
      common /parmr/ cutoff, tol
 
      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
 
      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)
 
      integer        ipq
      common /dim21/ ipq(256)
 
!     # dummy:
      integer lcore
      integer
     & ica(mcu,msu,*), idp(mstu,mstu,*), ipt(*),     la(mru,*),
     & lb(*),          lcr(*),           lls(*),     lmnp1(*),
     & lmnv(3,*),      lmn1u,            lproju,     maords(*),
     & mcons(*),       mcrs(*),          mgcs(*),    mnl(*),
     & ms(*),          nblpr(*),         nc(*),      ncon(*),
     & ncr(*),         ncru,             nct(*),     nd(*),
     & nf(*),          nir(*),           nkcrl(6,*), nkcru(6,*),
     & nklsl(4,*),     nklsu(4,*),       nrcr(*),    nso(*),
     & nsopr(*),       nt(*),            ntl(*),     ntu(*)
      real*8
     & a(*),               c(mctu,mcsu,*), ccr(*),   chg(*),
     & eta(mrcru,mconu,*), x(mcu,*),       y(mcu,*), z(mcu,*),
     & zcr(*),             zet(mconu,*)
      character*3 ityp(*), mtype(*)
c
c     # local:
      integer i, itol, icut, inrm, ncrs, fsplit, l1rec, l2rec,
     & ngen, ncons, naords, ist, ndpt, idpt, jst, kst, iprd,
     & iaords, iru, ir, igcs, icsu, ictu, ics, ict, ic, icons, iconu,
     & ircru, icon, ircr, kcru, icrs, lcru, llsu, lp1, nbfcr, kcrl,
     & kcr, l, ngenp1, isf, isfr, is, icu, ig, if, ilmnp1, lai,
     & nsf, itl, itu, jcts, it, jct, jt, ixyz, i2, i3, jconu,
     & jcon, js, jcu, jc, ngp1, iso, ibvtyp, l1recx,
     & l2recx, n1max, n2max, iau, ia, leig, lmn, lmnvmx
      integer ixyzir(3), lblso(16), iprdsm(3,7), ndptsm(0:9)
      real*8 twoc, prtint, sum, t, repnuc, chgprd, c1, ccc, eval

      character*1 ibl, ipc, isk, lblsh(21)

      real*8            a0,       a1s2,       a1,       a2,
     & aln10,           a4,       a10
      parameter(        a0=0.0d0, a1s2=0.5d0, a1=1.0d0, a2=2.0d0,
     & aln10=2.30258d0, a4=4.0d0, a10=10.0d0 )

c     # lrcmxp is the program default 1-e and 2-e record length.
c     # there is also an SIFS default, which may be different.
      integer    lrcmxp
      parameter( lrcmxp=4096 )

c     # for workspace allocation.
      integer  atebyt, forbyt
      external atebyt, forbyt

      integer  iodfac
      external iodfac

      data ibl, isk /' ','0'/, lblsh /'s','p','d','f','g','h','i','k',
     &  'l','m','n','o','q','r','t','u','v','w','x','y','z'/

c     # nst=4 and nst=8 product triples.
      data iprdsm / 2,3,4, 2,5,6, 2,7,8, 3,5,7, 3,6,8, 4,5,8, 4,6,7 /

c     # the number of product triples for each nst.
      data (ndptsm(i),i=0,9) /
     & -1, -1, 0, -1, 1, -1, -1, -1, 7, -1 /

      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,' * ECPINT : entering... ARGOS_INPUT'

c     # assign the default unit numbers and filenames.
      ninput = 5  ! RECP_INPUT :argos input file
      nlist  = 6  ! DIRAC output file 
      aoints = 4  ! aoints
      aoint2 = 8  ! aoints2
 
!    # open the input file.
      open(unit=ninput,status='old',file='RECP_INPUT')
 
      call ibummr(nlist)
 
      rewind ninput
 
!    # intitalize some input variables.
!    # [01-dec-90 aoint2, fsplit, l1rec, l2rec added. -rls]
      itol   =  0
      icut   =  0
      only1e =  0
      inrm   =  0
      ncrs   =  0
      l1rec  = -1
      l2rec  = -1
      fsplit =  2

      read (ninput,*)  ngen,   ns,     naords, ncons,  ngcs,   itol,
     & icut,   aoints, only1e, inrm,   ncrs,   l1rec,  l2rec,  aoint2,
     & fsplit
 
      if ( ns .gt. msu ) then
        call bummer('change msup (two places) to ',ns,faterr)
      elseif ( naords .gt. kaords ) then
        call bummer('change kaordp (one place) to ',naords,faterr)
      elseif ( ncons .gt. mconsu ) then
        call bummer('change mconsp (two places) to ',ncons,faterr)
      elseif ( ngcs .gt. mgcsu ) then
        call bummer('change mgcsup (one place) to ',ngcs,faterr)
      endif
 
!     # reset aoints, aoint2, only1e, and fsplit if necessary.
      if ( aoints .le. 0 ) aoints = 4
      if ( aoint2 .le. 0 ) aoint2 = 8
      if ( fsplit .ne. 1 ) fsplit = 2

      if ( only1e .ne. 0 ) then
c       # set only1e=0 for off, and only1e=1 for on.
        only1e = 1
        if ( fsplit .eq. 1 ) then
c         # this is inconsistent.  however, instead of treating
c         # this as an error, reset fsplit and continue.  the only
c         # harm is that the user may have to recompute the 1-e
c         # integrals later if this was not the right choice.
          write(nlist,*)'only1e is inconsistent with fsplit.'
     &     //' resetting fsplit=2 and continuing...'
          fsplit = 2
        endif
      endif
c
c     # read number of irreps, degeneracy, and character labels.
c
c     # *** the next time the input to this program is changed,
c     #     all of the items in this list should be read with
c     #     list-directed i/o.  ***
c
      read (ninput,'(i3,12(i3,a3))') nst,(nd(ist),ityp(ist),ist=1,nst)
c
      if ( nst .gt. mstu ) then
        call bummer('change mstup (two places) to ',nst,faterr)
      endif
c
      call izero_argos( nst,           nso, 1 )
      call izero_argos( mstu*mstu*nst, idp, 1 )
c
c     # assign the totally symmetric irrep binary products.
      do 104 ist = 1,nst
        idp(ist,ist,1) = 1
        idp(ist,1,ist) = 1
        idp(1,ist,ist) = 1
  104 continue
c
c     # read the number of product triples.
c
      ndpt = -1
      read(ninput,*) ndpt
c
      if ( ndpt .ge. 0 ) then
c
c       # read in the totally symmetric irrep product
c       # triples from the input file.
c
        do 108 idpt = 1, ndpt
c
          read (ninput,*) ist, jst, kst
c
          idp(ist,jst,kst) = 1
          idp(ist,kst,jst) = 1
          idp(jst,ist,kst) = 1
          idp(jst,kst,ist) = 1
          idp(kst,ist,jst) = 1
          idp(kst,jst,ist) = 1
  108   continue
      else
c
c       # assign default values based on nst.
c       # 04-dec-90 ndp<0 option added. -rls
c
        ndpt = ndptsm( min( max( 0, nst), 9 ) )
        if ( ndpt .lt. 0 ) then
          call bummer('syminp: unknown ndpt defaults for nst=',
     &     nst, faterr )
        else
          do 112 iprd = 1, ndpt
            ist = iprdsm(1,iprd)
            jst = iprdsm(2,iprd)
            kst = iprdsm(3,iprd)
c
            idp(ist,jst,kst) = 1
            idp(ist,kst,jst) = 1
            idp(jst,ist,kst) = 1
            idp(jst,kst,ist) = 1
            idp(kst,ist,jst) = 1
            idp(kst,jst,ist) = 1
  112     continue
        endif
      endif
c
      do 116 iaords = 1, naords
c
        read (ninput,*) iru, (la(ir,iaords), ir=1,iru)
c
        if ( iru .gt. mru ) then
          call bummer('change mrup (one place) to ',iru,faterr)
        endif
        nir(iaords) = iru
  116 continue
c
      do 128 igcs = 1,ngcs
c
        read (ninput,*) icsu, ictu, maords(igcs)
c
        if ( icsu .gt. mcsu ) then
          call bummer('change mcsup (one place) to ',icsu,faterr)
        elseif ( ictu .gt. mctu ) then
          call bummer('change mctup (one place) to ',ictu,faterr)
        endif
c
        do 124 ics = 1,icsu
          read (ninput,*) (c(ict,ics,igcs), ict=1,ictu)

          do ict = 1,ictu
            c1 = c(ict,ics,igcs)
            ccc = abs(c1)
            if ( ccc .ne. a0 .and. ccc .ne. a1 ) then
               c(ict,ics,igcs) = sign(sqrt(ccc),c1)
            endif
          enddo   
  124   continue
c
  128 continue
 
! ### read in exponents and contraction coefficients.
 
      do icons = 1,ncons

        read (ninput,*) iconu, lmnp1(icons), ircru
        if ( ircru .eq. 0 ) ircru = 1
        if ( ircru .gt. mrcru ) then
          call bummer('change mrcrup (one place) to ',ircru,faterr)
        elseif ( iconu .gt. mconu ) then
          call bummer('change mconup (one place) to ',iconu,faterr)
        endif
        nrcr(icons) = ircru
        ncon(icons) = iconu
 
        do icon = 1,iconu
          read (ninput,*) zet(icon,icons), (eta(ircr,icon,icons),
     &     ircr = 1,ircru)
        enddo
      enddo   
 
      if ( ncrs .eq. 0 ) go to 160
 
! ### read in expansions for core and spin-orbit potentials.
 
      ncru = -2
      lproju = 0
      kcru = 0
      do 156 icrs = 1,ncrs

!      # read number of core and spin-orbit
        read (ninput,*) lcru, llsu
        lproju = max(lproju, lcru-1, llsu)
 
!       if ( lcru .gt. 5 ) then
!         call bummer('lcru too large ',lcru,faterr)
!       elseif ( llsu .gt. 4 ) then
!         call bummer('llsu too large ',llsu,faterr)
!       endif
        lcr(icrs) = lcru
        lls(icrs) = llsu

!      # Read RECP (arep part)
        do lp1 = 1,(lcru+1)
          read (ninput,*) nbfcr
          kcrl = kcru+1
          kcru = kcru+nbfcr
          nkcrl(lp1,icrs) = kcrl
          nkcru(lp1,icrs) = kcru

          do kcr = kcrl,kcru
            read (ninput,*) ncr(kcr), zcr(kcr), ccr(kcr)
            if ( ncr(kcr) .lt. 0 ) then
              call bummer('ncr(kcr) too small ',ncr(kcr),faterr)
            endif
            ncru = max(ncru,ncr(kcr))
          enddo   
        enddo   

!      # Read RECP (spin-orbit part)
!      # Be careful with spin-orbit factor
        do l = 1,llsu 
          read (ninput,*) nbfcr
          kcrl = kcru+1
          kcru = kcru+nbfcr
          nklsl(l,icrs) = kcrl
          nklsu(l,icrs) = kcru

          do kcr = kcrl,kcru
            read (ninput,*) ncr(kcr), zcr(kcr), ccr(kcr)
            ccr(kcr)=ccr(kcr)*0.5  
            if ( ncr(kcr) .lt. -2 ) then
              call bummer('ncr(kcr) too small ',ncr(kcr),faterr)
            endif
            ncru = max(ncru,ncr(kcr))
          enddo
        enddo

  156 continue

      if ( kcru .gt. mcru ) then
        call bummer('change mcrup (one place) to ',kcru,faterr)
      endif
 
  160 continue
 
      ngenp1 = ngen + 1
      ng     = ngenp1
 
! ### read in atomic labels, charges, coordinates, basis sets,
! ### and core and spin-orbit potentials
 
      call izero_argos(ns,mcrs,1)
      lmn1u = 0
      isf = 0
      isfr = 0
      do 184 is = 1,ns
c
c       # *** the next time the input to this program is changed,
c       #     all of the items in this list should be read with
c       #     list-directed i/o.  ***
c
c       read (ninput,'(a3,2i3,f3.0)') mtype(is),nf(is),nc(is),chg(is)
        read (ninput,*) mtype(is),nf(is),nc(is),chg(is)
c
        icu = nc(is)
        if ( icu .gt. mcu ) then
          call bummer('change mcup (one place) to ',icu,faterr)
        endif

        do ic = 1,icu
          read (ninput,*) x(ic,is), y(ic,is), z(ic,is)
        enddo 
  
        if ( icu .gt. 1 ) then
          do ig = 2, ng
            read (ninput,*) (ica(ic,is,ig), ic=1,icu)
          enddo   
        endif
c
        do 180 if = 1,nf(is)
          isf = isf+1
c
          read (ninput,*) icons, igcs
c
          mcons(isf) = icons
          iaords = maords(igcs)
          iru =  nir(iaords)
          mgcs(isf) = igcs
          ilmnp1 = lmnp1(icons)
          lmn1u = max(lmn1u,ilmnp1)
          nt(isf) = ipq(ilmnp1+1)
          ntu(isf) = (nt(isf)*(ilmnp1+2))/3
          ntl(isf) = ntu(isf)-nt(isf)+1
          nct(isf) = icu*nt(isf)
          ircru = nrcr(icons)
          do 176 ircr = 1,ircru
            do 172 ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              nso(lai) = nso(lai)+1
              lb(isfr) = nso(lai)
  172       continue
  176     continue
  180   continue
        if ( ncrs .gt. 0 ) then
          read(ninput,*) mcrs(is)
        endif
  184 continue
      lmnvmx = (lmn1u*(lmn1u + 1)*(lmn1u + 2))/6
      ipt(2) = ipt(1) + forbyt( 3*lmnvmx )

c     # ipt(*)-->1:lmnv(1:3,1:lmnvmx)
      call mklmnv(lmn1u,lmnv)

c     # end of reading input file.
  192 continue
      close(unit=ninput)

!     ---------------------------

      nsf = isf
      if ( nsf .gt. msfu ) then
        call bummer('change msfup (two places) to ',nsf,faterr)
      elseif ( isfr .gt. msfru ) then
        call bummer('change msfrup (one place) to ',isfr,faterr)
      endif
c
c     # nbft = total number of symmetry orbital basis functions.
c
      nbft = isfr
c
c     # determine the number of types of one-electron integrals.
c
c     # initialize nu for s(*),t(*),v(*).
      nu = 3
c
      if ( ncrs .gt. 0 ) then
        do 196 is = 1,ns
          if ( mcrs(is) .ne. 0 ) then
            if ( lcr(mcrs(is)) .gt. 0 ) then
              nu = 4
              inam(4) = 4
              go to 200
            endif
          endif
  196   continue
c
  200   continue
        do 204 is = 1, ns
          if ( mcrs(is) .ne. 0 ) then
            if ( lls(mcrs(is)) .gt. 0 ) then
              nu = nu + 1
              inam(nu) = 5
              go to 208
            endif
          endif
  204   continue
c
  208   continue
      endif
c
      if ( inrm .ne. 0 ) then
c
c       # normalization of symmetry orbitals.
c
        do 244 igcs = 1,ngcs
          isf = 0
          do 240 is = 1,ns
            icu = nc(is)
            do 236 if = 1,nf(is)
              isf = isf+1
              if ( igcs .ne. mgcs(isf) ) go to 236
              itl = ntl(isf)
              itu = ntu(isf)
              iaords = maords(igcs)
              iru = nir(iaords)
              ilmnp1 = lmnp1(mcons(isf))
              ics = 0
              do 232 ir = 1,iru
                iau = nd(la(ir,iaords))
                do 228 ia = 1,iau
                  ics = ics+1
                  prtint = a0
                  ict = 0
                  do 220 ic = 1,icu
                    jcts = ict
                    do 216 it = itl,itu
                      ict = ict+1
                      if ( it .gt. itl ) then
                        twoc = a2*c(ict,ics,igcs)
                        jct = jcts
                        do 212 jt = itl,it-1
                          jct = jct+1
                          prtint = prtint + twoc*c(jct,ics,igcs)
     &                     * (iodfac(lmnv(1,it),lmnv(1,jt),lmnv(2,it),
     &                               lmnv(2,jt),lmnv(3,it),lmnv(3,jt)))
  212                   continue
                      endif
                      prtint = prtint + c(ict,ics,igcs)**2
     &                 * (iodfac(lmnv(1,it),lmnv(1,it),lmnv(2,it),
     &                           lmnv(2,it),lmnv(3,it),lmnv(3,it)))
  216                 continue
  220             continue
                  if ( prtint .ne. a1 ) then
                    prtint = sqrt(prtint)
                    do 224 jct = 1,ict
                      c(jct,ics,igcs) = c(jct,ics,igcs)/prtint
  224               continue
                  endif
  228           continue
  232         continue
              go to 244
  236       continue
  240     continue
  244   continue
      endif
c
      call izero_argos(3,ixyzir,1)
c
      if ( icut .le. 0 ) icut = 9
      cutoff = a1 / a10**icut
c
      if ( itol.le.0) itol = 20
      tol = aln10 * itol
c
c     # compute contracted orbital and lower-triangle-packed
c     # matrix symmetry offsets.
      nsopr(1) = 0
      nblpr(1) = 0
      do ist = 2,nst
         nsopr(ist) = nsopr(ist-1) + nso(ist-1)
         nblpr(ist) = nblpr(ist-1) + ipq( nso(ist-1) + 1 )
      enddo
c
c     # nnbft = total number of unique so function pairs.
c     # mpru = the number of elements in a diagonal,
c     #        symmetry-blocked, lower-triangle-packed matrix.  this
c     #        is also an upper bound to the size of any
c     #        symmetry-blocked array in which only unique
c     #        off-diagonal elements are stored and which is
c     #        associated with an operator that transforms as an
c     #        irreducible representation.
c
      nnbft = (nbft * (nbft + 1)) / 2
      mpru  = nblpr(nst) + ipq( nso(nst) + 1 )
c
c     # normalize the contraction coefficients.
c
      do 316 icons = 1, ncons
        iconu = ncon(icons)
        ircru = nrcr(icons)
        do 312 ircr =  1, ircru
          if ( iconu .eq. 1 ) then
            eta(1,1,icons) = a1
            go to 312
          endif
          prtint = eta(ircr,1,icons)**2
          do 304 icon = 2, iconu
            twoc  = a2 * eta(ircr,icon,icons)
            jconu = icon - 1
            do 300 jcon = 1, jconu
              t = (zet(icon,icons) + zet(jcon,icons)) * a1s2
              t = zet(icon,icons)/t * zet(jcon,icons)/t
              prtint = prtint + twoc * eta(ircr,jcon,icons)
     &         * sqrt(t**lmnp1(icons) * sqrt(t))
  300       continue
            prtint = prtint + eta(ircr,icon,icons)**2
  304     continue
          prtint = sqrt(prtint)
          do 308 icon = 1, iconu
            eta(ircr,icon,icons) = eta(ircr,icon,icons) / prtint
  308     continue
  312   continue
  316 continue
c
c     # calculate the nuclear repulsion energy.
c
      repnuc = a0
      do 332 is = 1, ns
        ic = nc(is)
        t = (ic) * chg(is)
c       # rearrangement theorem eliminates the ic loop in favor
c       # of the overall (ic) factor.
        do 328 js = 1, is
          chgprd = t * chg(js)
          if ( js .eq. is ) then
c           # repulsions within a set are double-counted, so
c           # compensate with a factor of 1/2.
            chgprd = a1s2 * chgprd
c           # don't include nc(is)=nc(js) self-repulsions.
            jcu = nc(js) - 1
          else
            jcu = nc(js)
          endif
          do 324 jc = 1, jcu
            repnuc = repnuc + chgprd / sqrt( (x(ic,is)-x(jc,js))**2
     1       + (y(ic,is)-y(jc,js))**2 + (z(ic,is)-z(jc,js))**2 )
  324     continue
  328   continue
  332 continue
c
      call gcentr( ica, nc )
c
      ngp1 = ng+1
      if ( ngp1 .gt. mgu ) then
        call bummer('change mgup (two places) to ',ngp1,faterr)
      endif
c

!     write (nlist,612) ngen,   ns,     naords, ncons,  ngcs,   itol,
!    &  icut,   aoints, only1e, inrm,   ncrs,   l1rec,  l2rec,  aoint2,
!    &  fsplit
! 612 format(/' argos input parameters :'/
!    & ' ngen   =',i3,' ns     =',i3,' naords =',i3,' ncons  =',i3,
!    & ' ngcs   =',i3,' itol   =',i3/
!    & ' icut   =',i3,' aoints =',i3,' only1e =',i3,' inrm   =',i3,
!    & ' ncrs   =',i3/
!    & ' l1rec  =',i10,           5x,' l2rec  =',i10,           5x,
!    & ' aoint2 =',i3,' fsplit =',i3/)

c
c     write (nlist,616) (ist, ist=1,nst)
c 616 format(//'irrep',9x,12(i4,4x))
c     write (nlist,620) (nd(ist), ist=1,nst)
c 620 format('degeneracy',12i8)
c     write (nlist,624) (ityp(ist), ist=1,nst)
c 624 format('label',11x,12(a3,5x))
c     write (nlist,628)
c 628 format(//'direct product table')
c     do 644 ist = 1,nst
c       do 640 jst = 1,nst
c         do 636 kst = 1,nst
c           if ( idp(ist,jst,kst) .eq. 1 )
c    1       write (nlist,632) ityp(ist), ityp(jst), ityp(kst)
c 632        format(3x,'(',a3,') (',a3,') = ',a3)
c 636     continue
c 640   continue
c 644 continue
c

      IF (PRINT_LEVEL .GE. 2) THEN
         write (nlist,'(/,a,f14.8)') 
     &         ' * nuclear repulsion energy', repnuc
         write (nlist,'(a,a,i2,a)') 
     &         ' * primitive ao integrals neglected', 
     &         'if exponential factor below 10**(-',itol,')'
         write (nlist,'(a,a,i2,a)') 
     &         ' * contracted ao and so integrals neglected ',
     &         'if value below 10**(-',icut,')'
      ENDIF

      isf = 0
      isfr = 0

      do 840 is = 1,ns
        icu = nc(is)

        IF (PRINT_LEVEL .GE. 2) THEN        
           write(nlist,656) mtype(is), chg(is),
     &                     (ic, x(ic,is), y(ic,is), z(ic,is), ic=1,icu)
  656      format(//36x,a3,' atoms'
     &            //30x,'nuclear charge',f7.2
     &            //11x,'center',12x,'x',15x,'y',15x,'z'
     &             /(i14,7x,3f16.8))
           if ( icu .gt. 1 ) then
              write (nlist,660) (ica(ic,is,1), ic=1,icu)
  660         format(' operator',6x,'center interchanges'/
     &               '  1(id.)',6x,24i3)

              do ig = 2,ngenp1
                 write (nlist,664) ig, (ica(ic,is,ig), ic=1,icu)
  664            format(i3,'(gen.)',5x,24i3)
              enddo   

              do ig = ngenp1+1,ng
                 write(nlist,'(i3,11x,24i3)') 
     &                ig,(ica(ic,is,ig), ic=1,icu)
              enddo   
           endif
        ENDIF

        do 796 if = 1,nf(is)
          isf = isf+1
          itl = ntl(isf)
          itu = ntu(isf)
          igcs = mgcs(isf)
          iaords = maords(igcs)
          iru = nir(iaords)
          ictu = nct(isf)
          icons = mcons(isf)
          ircru = nrcr(icons)
          iconu = ncon(icons)
          icsu = 0
          do 704 ir = 1,iru
            icsu = icsu + nd(la(ir,iaords))
  704     continue

          IF (PRINT_LEVEL .GE. 2) THEN        
              CALL ARGOS_INPUT_INFO2(
     &             LMNP1, ICONS, IRCRU, ICONU, ZET, ETA, ICSU, IGCS,
     &             ISFR, IRU, LA, ND, LB, ITYP, C, ICU, ITL, ITU, 
     &             LMNV, IAORDS)
          ENDIF

          do 792 ircr = 1,ircru
            ics = 0
            do 788 ir = 1,iru
              isfr = isfr+1
              lai = la(ir,iaords)
              iso = nsopr(lai)+lb(isfr)
              ms(iso) = is
c
c             # assign the mnl(*) values to the so functions.
c             # 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g,...
c             #
c             # compute the L^2-operator expectation value and
c             # determine if it is an eigenvalue.
c
              lmn = lmnp1(icons) - 1
c             numxyz = ((lmn + 1) * (lmn + 2)) / 2
c             # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:scr(1:xyzdim)
              call l2opxv( lmn, c(1,ics+1,igcs), itl, itu, lmnv,
     &         ipt(55)-ipt(2), a(ipt(2)), eval, leig )
              if ( leig .lt. 0 ) then

                IF (PRINT_LEVEL .GE. 2) THEN 
!                  # not an eigenvalue. print a warning.
                   write(nlist,768) iso, ityp(lai), eval
  768              format(' * warning: from l2opxv()',
     &               /,'   s.o. basis function',i3,a,' is not an ',
     &                 'L^2 eigenvector.',
     &               /,'   <v|L^2|v>/<v|v> =', f12.8 )
                ENDIF 
!               # just assign the maximum value.
                leig = lmn
              endif
              mnl(iso) = ( lmnp1(icons)**2 + 2 * leig + 5 ) / 4
              iau = nd(lai)
              do 784 ia = 1,iau
                ics = ics+1
 
                if ( inam(nu) .ge. 5 .and. leig .eq. 1 ) then
 
!                 # determine angular momentum symmetries.
                  do 780 ixyz = 1, 3
                    i2 = mod(ixyz,  3) + 1
                    i3 = mod(ixyz+1,3) + 1
                    sum = a0
                    ict = 0
                    do ic = 1, icu
                      do it = itl,itu
                        ict = ict + 1
                        if( mod(lmnv(ixyz,it),2) .eq. 1 .and.
     &                      mod(lmnv(i2,  it),2) .eq. 0 .and.
     &                      mod(lmnv(i3,  it),2) .eq. 0 ) then
                          sum = sum + c(ict,ics,igcs)
                        endif
                      enddo   
                    enddo   

                    if (  (sum .eq. a0)
     &               .or. (lai .eq. ixyzir(ixyz)) ) go to 780
                    if ( ixyzir(ixyz) .ne. 0 ) then
                      write (nlist,*) 'axis choice is unsuitable for'
     &                 //' spin-orbit integrals'
                      call bummer('syminp: ixyzir(ixyz)=',
     &                 ixyzir(ixyz),faterr)
                    endif
                    ixyzir(ixyz) = lai
  780             continue
                endif
  784         continue
  788       continue
  792     continue
  796   continue

c
        icrs = mcrs(is)
        if ( icrs .eq. 0 ) go to 840
        lcru = lcr(icrs)
        if ( lcru .lt. 0 ) go to 816
        kcrl = nkcrl(1,icrs)
        kcru = nkcru(1,icrs)

        IF (PRINT_LEVEL .GE. 2) THEN        
          WRITE (NLIST,1000)        
          write (nlist,804) mtype(is),lblsh(lcru+1),
     &          (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  804     format(//31x,a3,' core potential'//35x,a1,' potential'/
     &          19x,'powers',6x,'exponentials',4x,'coefficients'/
     &          (i22,7x,2(1pg16.7)))
          if ( lcru .eq. 0 ) go to 816
          do 812 lp1 = 1,lcru
             kcrl = nkcrl(lp1+1,icrs)
             kcru = nkcru(lp1+1,icrs)
             write (nlist,808) lblsh(lp1), lblsh(lcru+1),
     &             (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  808        format(/33x,a1,' - ',a1,' potential'/
     &             19x,'powers',6x,'exponentials',4x,'coefficients'/
     &             (i22,7x,2(1pg16.7)))
  812     continue
c
  816     continue
          llsu = lls(icrs)
          if ( llsu .le. 0 ) go to 840
          write (nlist,820) mtype(is)
  820     format(//28x,a3,' spin-orbit potential')
          do 836 l = 1,llsu
             kcrl = nklsl(l,icrs)
             kcru = nklsu(l,icrs)
             write (nlist,824) lblsh(l+1),
     &       (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  824        format(/35x,a1,' potential'/
     &             19x,'powers',6x,'exponentials',4x,'coefficients'/
     &             (i22,7x,2(1pg16.7)))
  836     continue
        ENDIF

  840 continue

!     -------------------------------
      WRITE(NLIST,1000)

      if ( inam(nu) .eq. 5 ) then
        do 860 ixyz = 1, 3
          i2 = mod(ixyz,  3) + 1
          i3 = mod(ixyz+1,3) + 1
          do 856 ist = 1, nst
            if ( idp(ist,ixyzir(i2),ixyzir(i3)) .ne. 0 ) then
              lxyzir(ixyz) = ist
              go to 860
            endif
  856     continue
  860   continue

        IF (PRINT_LEVEL.GE.5) THEN
           write (nlist,864) (ityp(lxyzir(ixyz)), ixyz=1,3)
  864      format(/'lx: ',a3,10x,'ly: ',a3,10x,'lz: ',a3)
        ENDIF
      endif
 
c     # set up the 1-e record parameters.
      ibvtyp = 0
      if ( l1rec .eq. 0 ) then
c       # use the default value for this program.
        l1recx = lrcmxp
      else
c       # use the input value. [(-1) invokes the SIFS default.]
        l1recx = l1rec
      endif
 
      call sifcfg( 1, l1recx, nbft, ibvtyp, ifmt1, l1rec, n1max, ierr )
      IF (PRINT_LEVEL.GE.5) THEN
         PRINT *, '   L1RECX :',L1RECX
         PRINT *, '   N1MAX  :',N1MAX
      ENDIF
 
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from 1-e sifcfg(), ierr=',ierr,faterr)
      endif
 
c     # set up the 2-e record parameters.
      ibvtyp = 1
      if ( l2rec .eq. 0 ) then
c       # use the default value for this program.
        l2recx = lrcmxp
      else
c       # use the input value. [(-1) envokes the SIFS default.]
        l2recx = l2rec
      endif
 
      call sifcfg( 2, l2recx, isfr, ibvtyp, ifmt2, l2rec, n2max, ierr )
      if ( ierr .ne. 0 ) then
        call bummer('syminp: from 2-e sifcfg(), ierr=',ierr,faterr)
      endif
 
!     # set up the info(*) array for the output file.
      info(1) = fsplit ! file-split parameter  fsplit=2
      info(2) = l1rec  ! maximum record length for 1-e integral records
      info(3) = n1max
      info(4) = l2rec  ! maximum record length for 1-e integral records
      info(5) = n2max
 
!     # write out the header information
!     # ipt(*)-->1:lmnv(1:3,1:lmnvmx), 2:map(1:2,1:nbft)
!     # push map(*) allocation onto the stack.
 
      call wthead( nst,    ns,     isfr,   repnuc, ityp,
     &   nso,      mtype,  ms,     mnl,    info,   aoints, nlist,
     &   a(ipt(2)) )
 
!     # multiply normalization constants into contraction coefficients.
      do icons = 1,ncons
         iconu = ncon(icons)
         ircru = nrcr(icons)
         do icon = 1,iconu
            do ircr = 1,ircru
               t = a4 * zet(icon,icons)
               eta(ircr,icon,icons) =
     &          sqrt( t**lmnp1(icons) / a2 * sqrt(t / a2) )
     &          * eta(ircr,icon,icons)
            enddo   
         enddo   
      enddo   

      IF (PRINT_LEVEL .GE. 2) THEN
         write (nlist,880) isfr, (ityp(ist), nso(ist), ist=1,nst)
  880    format(//i9,' symmetry orbitals,',3x,4(3x,a3,':',i4),
     &                                   (/31x,4(3x,a3,':',i4)))
      ENDIF
 
      IF (PRINT_LEVEL.GE.5) 
     &   PRINT *,' * ECPINT : exiting...  ARGOS_INPUT'

      return
 1000 FORMAT(/,10X,
     & '-----------------------------------------------------------')
      end

*deck mklmnv (lmnvgn)
      subroutine mklmnv(lmn1u,lmnv)
!     # generate cartesian gaussian exponent array.
!     # lmnv(*,*) = exponents of the cartesian gaussian basis functions.
!     #             s  p  d   f   g   h   i
!     #       lmn = 0  1  2   3   4   5   6
!     #    numxyz = 1, 3, 6, 10, 15  21  28      =((lmn+1)*(lmn+2))/2
      implicit logical (a-z)
      integer i, j, k, lmn, lmn1u, lmnv(3,*), ndx
 
      ndx = 0
      do 40 lmn = 0,lmn1u-1
        do 30 i = lmn,(lmn+2)/3,-1
          do 20 j = min(i,lmn-i),(lmn-i+1)/2,-1
            k = lmn - i -j
            if( i.eq.j .and. j.eq.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = i
              lmnv(3,ndx) = i
            elseif( i.eq.j .and.j.gt.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = i
              lmnv(3,ndx) = k
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = k
              lmnv(3,ndx) = i
              ndx = ndx + 1
              lmnv(1,ndx) = k
              lmnv(2,ndx) = i
              lmnv(3,ndx) = i
            elseif( i.gt.j .and. j.eq.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = j
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = i
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = j
              lmnv(3,ndx) = i
            elseif( i.gt.j .and. j.gt.k ) then
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = j
              lmnv(3,ndx) = k
              ndx = ndx + 1
              lmnv(1,ndx) = i
              lmnv(2,ndx) = k
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = i
              lmnv(3,ndx) = k
              ndx = ndx + 1
              lmnv(1,ndx) = k
              lmnv(2,ndx) = i
              lmnv(3,ndx) = j
              ndx = ndx + 1
              lmnv(1,ndx) = j
              lmnv(2,ndx) = k
              lmnv(3,ndx) = i
              ndx = ndx + 1
              lmnv(1,ndx) = k
              lmnv(2,ndx) = j
              lmnv(3,ndx) = i
            endif
   20     continue
   30   continue
   40 continue
civan-s
c  *********************
c     new ordering of 
c      d exponents
c  *********************
      lmnv(1,5) = 2     
      lmnv(2,5) = 0
      lmnv(3,5) = 0
      lmnv(1,6) = 1
      lmnv(2,6) = 1
      lmnv(3,6) = 0
      lmnv(1,7) = 1
      lmnv(2,7) = 0
      lmnv(3,7) = 1
      lmnv(1,8) = 0
      lmnv(2,8) = 2
      lmnv(3,8) = 0
      lmnv(1,9) = 0
      lmnv(2,9) = 1
      lmnv(3,9) = 1
      lmnv(1,10) = 0
      lmnv(2,10) = 0
      lmnv(3,10) = 2
c  *********************
c     new ordering of
c      f exponents  
c  *********************
      lmnv(1,11) = 3
      lmnv(2,11) = 0
      lmnv(3,11) = 0
      lmnv(1,12) = 2
      lmnv(2,12) = 1
      lmnv(3,12) = 0
      lmnv(1,13) = 2
      lmnv(2,13) = 0
      lmnv(3,13) = 1
      lmnv(1,14) = 1
      lmnv(2,14) = 2
      lmnv(3,14) = 0
      lmnv(1,15) = 1
      lmnv(2,15) = 1
      lmnv(3,15) = 1
      lmnv(1,16) = 1
      lmnv(2,16) = 0
      lmnv(3,16) = 2
      lmnv(1,17) = 0
      lmnv(2,17) = 3
      lmnv(3,17) = 0
      lmnv(1,18) = 0
      lmnv(2,18) = 2
      lmnv(3,18) = 1
      lmnv(1,19) = 0
      lmnv(2,19) = 1
      lmnv(3,19) = 2
      lmnv(1,20) = 0
      lmnv(2,20) = 0
      lmnv(3,20) = 3
c  *********************
c     new ordering of
c      g exponents  
c  *********************
      lmnv(1,21) = 4
      lmnv(2,21) = 0
      lmnv(3,21) = 0
      lmnv(1,22) = 3
      lmnv(2,22) = 1
      lmnv(3,22) = 0
      lmnv(1,23) = 3
      lmnv(2,23) = 0
      lmnv(3,23) = 1
      lmnv(1,24) = 2
      lmnv(2,24) = 2
      lmnv(3,24) = 0
      lmnv(1,25) = 2
      lmnv(2,25) = 1
      lmnv(3,25) = 1
      lmnv(1,26) = 2
      lmnv(2,26) = 0
      lmnv(3,26) = 2
      lmnv(1,27) = 1
      lmnv(2,27) = 3
      lmnv(3,27) = 0
      lmnv(1,28) = 1
      lmnv(2,28) = 2
      lmnv(3,28) = 1
      lmnv(1,29) = 1
      lmnv(2,29) = 1
      lmnv(3,29) = 2
      lmnv(1,30) = 1
      lmnv(2,30) = 0
      lmnv(3,30) = 3
      lmnv(1,31) = 0
      lmnv(2,31) = 4
      lmnv(3,31) = 0
      lmnv(1,32) = 0
      lmnv(2,32) = 3
      lmnv(3,32) = 1
      lmnv(1,33) = 0
      lmnv(2,33) = 2
      lmnv(3,33) = 2
      lmnv(1,34) = 0
      lmnv(2,34) = 1
      lmnv(3,34) = 3
      lmnv(1,35) = 0
      lmnv(2,35) = 0
      lmnv(3,35) = 4
civan-e
      return
      end

*deck wthead
      subroutine wthead( nsym,   ns,     nbft,
     &   repnuc, ityp,   nbpsy,  mtype,  ms,     mnl,    info,
     &   aoints, nlist,  map )
!  write the header info to the output aoints file.
!  input:
!   nsym = number of symmetry blocks.
!   ns   = number of symmetry inequivalent atom labels.
!   nbft = total number of basis functions.
!   repnuc = nuclear repulsion energy.
!   ityp(1:nsym) = input symmetry labels.
!   nbpsy(1:nsym) = number of basis functions in each symmetry block.
!   mtype(1:ns) = atom labels.
!   ms(1:nbft) = basis_function-to-center mapping vector.
!   mnl(1:nbft) = integer codes for basis function types.
!                 1:s, 2:p, 3:3s, 4:d, 5:4p, 6:f, 7:5s, 8:5d, 9:g
!   info(*) = info(*) array for the output file.
!   aoints = integral output file unit number.
!   nlist = listing file output unit.
!   map(1:2,1:nbft) = temporary scratch array for writing the output
!                     map(*) vectors.
!  01-dec-90 SIFS version written by ron shepard.

      implicit logical(a-z)
c
      integer    nbfmx,     nmap,   ninfo,   nenrgy,   nrtype
      parameter( nbfmx=255, nmap=2, ninfo=5, nenrgy=1, nrtype=-1 )

c     # dummy:
      integer nsym, ns, nbft, aoints, nlist
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
#include "argos_error.h"

      INTEGER  PRINT_LEVEL
 
      integer    bmin,   bmax
      parameter( bmin=0, bmax=17 )
c
      integer      clip, i, j
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
        WRITE (NLIST,1000) nbfmx,(nbft-nbfmx)
 1000   FORMAT(/
     &  1x,'In this version of ECP integral routine,'/
     &  1x,'the number of basis function are limited to',i4,/
     &  1x,i4,' basis functions are exceed')
        CALL QUIT('Number of basis function are exceed.')
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
c     # setup for a single nuclear repulsion energy.
c
      energy(1) = repnuc
      ietype(1) = nrtype


      IF (PRINT_LEVEL .GE. 2) THEN
c
c     # write the headers.
c
      write(nlist,6100)'RECP integral information:'
      write(nlist,6100)'output energy(*) values:'
      call sifpre( 6, nenrgy, energy, ietype )
c
      write(nlist,6200)'total core energy :',
     & sifsce( nenrgy, energy, ietype )

!
      write(nlist,6041) nsym
      write(nlist,6042) nbft
      write(nlist,6051) (i,i=1,nsym)
      write(nlist,6052) (slabel(i),i=1,nsym)
      write(nlist,6053) (nbpsy(i),i=1,nsym)

!    # the info(*) entries were passed in.
      write(nlist,'(/,a,i10)') 
     & ' file-split parameter                   (fsplit):', info(1)  
      write(nlist,'(a,i10)') 
     & ' maximum record length for 1-e integral (l1rec) :', info(2)
      write(nlist,'(a,i10)') 
     & '                                        (n1max) :', info(3)  
      write(nlist,'(a,i10)') 
     & ' maximum record length for 2-e integral (l2rec) :', info(4)
      write(nlist,'(a,i10)') 
     & '                                        (n2max) :', info(5)  

!    # orbital labels
      write(nlist,'(/,a)') ' output orbital labels (i:bfnlab(i))'
      write(nlist,'(5(i4,a1,a8))') (i,':',bfnlab(i),i=1,nbft)

!    # basis function to atomic center map
      write(nlist,'(/,a)') ' basis function to center (i:map(i))'
      write(nlist,'(5(i4,a1,i4,4x))') (i,':',map(i,1),i=1,nbft)

!    # basis function to orbital type map
      write(nlist,'(/,a)') ' basis function to orbital type (i:map(i))'
      write(nlist,'(5(i4,a1,i4,4x))') (i,':',map(i,2),i=1,nbft)

      ENDIF

6041  format(/1x,'number of symmetry operator = ',2x,i2)
6042  format( 1x,'number of basis functions   = ',i4)
6051  format(/1x,'symmetry number  = ',8i5)
6052  format( 1x,'symmetry label   = ',x,8(1x,a4))
6053  format( 1x,'basis function   = ',8i5)
6100  format(/1x,a,(1x,10i10))
6200  format(/1x,a,(1x,1p4e20.12))
c
      return
      end


      SUBROUTINE ARGOS_INPUT_INFO2(
     &           LMNP1, ICONS, IRCRU, ICONU, ZET, ETA, ICSU, IGCS,
     &           ISFR, IRU, LA, ND, LB, ITYP, C, ICU, ITL, ITU,
     &           LMNV, IAORDS) 
      IMPLICIT NONE

      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst

      integer
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & kaords, mconsu, mgcsu, mru,  mcsu, mctu, mcru, msfu, mgu,
     & msfru,  mnru,   ngcs,  nu,   mcxu2,
     & lxyzir(3),      inam(5),     nnam, mdum(32)

      INTEGER     NLIST

      INTEGER     LMNP1(*), I
      REAL*8      ZET(MCONU,*), ETA(MRCRU,MCONU,*), C(MCTU,MCSU,*) 
      INTEGER     ICONS, ICONU, ICON
      INTEGER     IRCRU, IRCRH, IRCR, IRCRL
      INTEGER     JSFR, ISFR, IR, IRU
      INTEGER     LAI, LA(MRU,*), IAORDS, IA, IAU, ND(*), LB(*)
      INTEGER     ICSR, ICOL, NCOL, LBLA(16), ICT, IC, IT
      INTEGER     ICU, ITL, ITU, LMNV(3,*), ICSU, IGCS 
      CHARACTER*3 ITYP(*), LTYP(16)
      CHARACTER*1 LBLSH(21)
      DATA        LBLSH
     &             /'s','p','d','f','g','h','i','k','l','m','n',
     &              'o','q','r','t','u','v','w','x','y','z'/

      NLIST = 6

      write (nlist,700)
!     # orbital exponents  contraction coefficients for each orbital
      write (nlist,708) LMNP1(ICONS), lblsh(lmnp1(icons))
  708 format(/i21,a1,' orbitals'//
     &       ' orbital exponents  contraction coefficients')
!     use larger value than 4 for output lines longer than 80
      ircrh = min(4,ircru)
      do 712 icon = 1,iconu
        write (nlist,'(8(1pg16.7))') zet(icon,icons),
     &    (eta(ircr,icon,icons), ircr=1,ircrh)
  712 continue

!     # symmetry orbital labels
      write (nlist,"(/21x,'symmetry orbital labels')")
      jsfr = isfr
      do 720 ir = 1,iru
        jsfr = jsfr + 1
        lai = la(ir,iaords)
        do 716 ia = 1,nd(lai)
          write (nlist,'(10x,7(i12,a3,i1))') (lb(jsfr+(ircr-1)*
     &      icsu),ityp(lai),ia, ircr=1,ircrh)
  716   continue
  720 continue

  724 if(ircrh.lt.ircru) then
        ircrl = ircrh + 1
        ircrh = min(ircrh+4,ircru)
        write (6,"(/20x,'contraction coefficients')")
        do 728 icon = 1,iconu
          write (nlist,'(16x,7(1pg16.7))')
     &      (eta(ircr,icon,icons), ircr=ircrl,ircrh)
  728   continue
        write (nlist,"(/21x,'symmetry orbital labels')")
        jsfr = isfr
        do 736 ir = 1,iru
          jsfr = jsfr + 1
          lai = la(ir,iaords)
          do 732 ia = 1,nd(lai)
            write (nlist,'(10x,7(i12,a3,i1))') (lb(jsfr+(ircr-1)*
     &        icsu),ityp(lai),ia, ircr=ircrl,ircrh)
  732     continue
  736   continue
        go to 724
      endif

      write (nlist,"(/11x,'symmetry orbitals')")
      icsr = 0
      ncol = 0
      do 760 ir = 1,iru
        lai = la(ir,iaords)
        iau = nd(lai)
        do 756 ia = 1,iau
          ncol = ncol+1
          ltyp(ncol) = ityp(lai)
          lbla(ncol) = ia
          if ( ncol .ne. 10 .and.
     &      (ir .ne. iru .or. ia .ne. iau) ) go to 756
          write (nlist,740) (ltyp(icol),lbla(icol), icol=1,ncol)
  740     format(' ctr, ao',16(3x,a3,i1))
          ict = 0
          do 752 ic = 1,icu
            do 748 it = itl,itu
              ict = ict+1
              write (nlist,744) ic, (lmnv(i,it), i=1,3),
     &          (c(ict,icsr+icol,igcs), icol=1,ncol)
  744         format(i3,', ',3i1,16f7.3)
  748       continue
  752     continue
          icsr = icsr+ncol
          ncol = 0
  756   continue
  760 continue
      
      WRITE(NLIST,*)
      RETURN

  700 FORMAT(/,10X,
     & '-----------------------------------------------------------')
      END
