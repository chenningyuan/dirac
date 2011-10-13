      subroutine recp_socfpd( mcxu,mcxu2,mblu,mgcsu,icxast,icxst,
     &   icxsv1, icxsv2, idp,  iprst,  la,
     &   lb,     maords, mau,  mcons,      mgcs,   mics,
     &   mjcs,   nblpr,  nc,   nd,         nf,     nir,
     &   npair,  nprir,  nrcr, nsopr,      nt,     ntl,
     &   ntu, CX_TEMP)
!    Compute symmetry orbital coefficient products.
!    on return from socfpd(*), cx(1:mcxu) was computed.

      USE RECP_IPT
      USE RECP_OUTPUT
      implicit logical(a-z)

!     il(*) packing factor. see also oneint() and wtint2().
      integer    ilfact
      parameter( ilfact=1024 )

#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"
#include "inc_print.h"

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
      integer mcxu,mcxu2,mblu,mgcsu,CX_TEMP, IPQ_TEMP

      integer
     & icxast(*),        icxst(2,*),   icxsv1(*),    icxsv2(mgcsu,*),
     & idp(mstu,mstu,*), iprst(*),     la(mru,*),
     & lb(*),            maords(*),    mau(*),       mcons(*),
     & mgcs(*),          mics(*),      mjcs(*),      nblpr(*),
     & nc(*),            nd(*),        nf(*),        nir(*),
     & npair(2,*),       nrcr(*),      nsopr(*),     nprir(2,mstu,*),
     & nt(*),            ntl(*),       ntu(*)

c     # local:
      integer mcxu1, icx, ipair, ijsf, isf, isfr, igcs, jgcs, is,
     & ifu, icu, isfis, isfris, jsf, jsfr, js, jsfjs, jsfrjs, jfu,
     & jcu, if, isfrif, iaords, iru, ircru, itl, itu, jf, jaords,
     & jtl, jsfrjf, jtu, jru, ist, jrcru, ircr, isfrib, jrcr,
     & jsfrjb, iesfb, npr, ipr, ics, ir, lai, jcs, jr, laj,
     & isv1, icxu, icxa, isv2, ict, ic, ictic, jct, jc, jctjc, it,
     & jt, ia, npru, I,J,K
      logical ea, easfb, es, esf, esfb, esfbc, msfbct

      CALL OUTPUT_LOC('ARGOS_SOCFPD','E')

      DO I = 1, 2
         DO J = 1, mstup 
            DO K = 1, mnsfup
               nprir(I,J,K) = 0
            ENDDO
         ENDDO
      ENDDO
!     # check the il(*) packing factor.
      if ( nbft .ge. ilfact ) then
        call bummer('socfpd: larger ilfact required, nbft=',
     &   nbft, faterr )
      endif
 
      ea = inam(nu).eq.5
 
!     # mcxu1 is the original mcxu offset.
      mcxu  = CX_TEMP 
      mcxu1 = mcxu
      icx   = 0
      ipair = 0
 
!     # ijsf = one dimensionalization of is, js, if, jf
!     # isfr = one dimensionalization of is, if, ir
      ijsf = 0
      isf  = 0
      isfr = 0
      do igcs = 1, ngcs
        icxsv1(igcs) = 0
        do jgcs = 1, ngcs
          icxsv2(jgcs,igcs) = 0
        enddo
      enddo

      do 340 is = 1, ns      ! icenter
        ifu = nf(is)
        if(ifu.eq.0) go to 340
        icu = nc(is)
        isfis = isf
        isfris = isfr
        jsf = 0
        jsfr = 0
        do 330 js = 1,is     ! jcenter
          jsfjs = jsf
          jsfrjs = jsfr
          es = is.eq.js
          if(.not.es) then
            jfu = nf(js)
            if(jfu.eq.0) go to 330
            jcu = nc(js)
          endif
          isf = isfis
          isfr = isfris
          do 324 if = 1,ifu   ! i-function
            isf = isf+1
            isfrif = isfr
            igcs = mgcs(isf)
            iaords = maords(igcs)
            iru = nir(iaords)
            ircru = nrcr(mcons(isf))
            itl = ntl(isf)
            itu = ntu(isf)
            jsf = jsfjs
            jsfr = jsfrjs
            if(es) jfu = if
            do 322 jf = 1,jfu ! j-function
              jsf = jsf+1
              jgcs = mgcs(jsf)
              jaords = maords(jgcs)
              jtl = ntl(jsf)
              jsfrjf = jsfr
              ijsf = ijsf+1
              esf = isf.eq.jsf
              if(esf) then
                if(ircru.ge.2) then
                  jtu = ntu(jsf)
                  jru = nir(jaords)
                else
                  npair(2,ijsf) = 0
                  icxst(2,ijsf) = 0
                  do 20 ist = 1,nst
                    nprir(2,ist,ijsf) = 0
20                continue
                endif
              else
                jrcru = nrcr(mcons(jsf))
                jtu   = ntu(jsf)
                jru   = nir(jaords)
              endif
              iprst(ijsf) = ipair
              isfr = isfrif
              do 320 ircr = 1,ircru
                isfrib = isfr
                jsfr = jsfrjf
                if(esf) jrcru = ircr
                do 310 jrcr = 1,jrcru
                  jsfrjb = jsfr
                  esfb = esf.and.ircr.eq.jrcr
                  easfb = ea.and.esfb
                  if(esfb) then
                    iesfb = 1
                    npru = (iru * (iru + 1)) / 2
                  else
                    iesfb = 2
                    npru = iru * jru
                  endif
 
                  if ( npru .gt. mnru ) then
                    call bummer('change mnrup (one place) to ',
     &               npru,faterr)
                  endif
 
                  npr = 0
                  do 75 ist = 1, nst
                    ipr = npr
                    ics = 0
                    isfr = isfrib
                    do 70 ir = 1, iru
                      isfr = isfr + 1
                      lai = la(ir,iaords)
                      if(esfb) jru = ir
                      jsfr = jsfrjb
                      jcs = 0
                      do 65 jr = 1, jru
                        jsfr = jsfr + 1
                        laj = la(jr,jaords)
!                print *,'socfpd: idp',lai,laj,ist,'=',idp(lai,laj,ist)
                        if ( idp(lai,laj,ist) .ne. 0 ) then
                          ipair = ipair + 1
                          IPT_IL(ipair) = (nsopr(lai) + lb(isfr))
     &                     * ilfact + (nsopr(laj) + lb(jsfr))
 
                          npr = npr + 1
                          mau(npr)  = nd(lai)
                          mics(npr) = ics
                          mjcs(npr) = jcs
                        endif
                        jcs = jcs + nd(laj)
65                    continue
                      ics = ics + nd(lai)
70                  continue
                    nprir(iesfb,ist,ijsf) = npr - ipr
75                continue
                  npair(iesfb,ijsf) = npr
 
                  if ( npr .eq. 0 ) go to 310
 
                  if ( esfb ) then
                    isv1 = icxsv1(igcs)
                    if ( isv1 .ne. 0 ) then
                      icxst(1,ijsf) = icxst(1,isv1)
                      if(ea) icxast(ijsf) = icxast(isv1)
                      go to 310
                    endif
                    icxsv1(igcs) = ijsf
                    icxst(1,ijsf) = icx
                    IPQ_TEMP = (icu*nt(isf)+1)*( (icu*nt(isf)+1)-1 )/2
                    icxu = icx+(npr*IPQ_TEMP)
!                   print *,'socfpd-icxu :',icxu
                    if ( ea ) then
                      IPQ_TEMP = (icu*nt(isf)+1)*((icu*nt(isf)+1)-1)/2
                      mcxu = mcxu-(npr*IPQ_TEMP)
!                     print *,'socfpd-mcxu :',mcxu
                      icxa = mcxu
                      icxast(ijsf) = icxa
                    endif
                  else
                    isv2 = icxsv2(jgcs,igcs)
                    if ( isv2 .ne. 0 ) then
                      icxst(2,ijsf) = icxst(2,isv2)
                      go to 310
                    endif
                    icxsv2(jgcs,igcs) = ijsf
                    icxst(2,ijsf) = icx
                    icxu = icx+npr*icu*nt(isf)*jcu*nt(jsf)
                  endif
                  if ( icxu .gt. mcxu ) then
                    call bummer('socfpd: (icxu-mcxu)=',
     &               (icxu-mcxu), faterr )
                  endif
 
!                 # compute new cx block
                  ict = 0
                  do 298 ic = 1, icu
                    ictic = ict
                    if(esfb) jcu = ic
                    jct = 0
                    do 296 jc = 1, jcu
                      jctjc = jct
                      esfbc = esfb .and. ic .eq. jc
                      ict = ictic
                      do 294 it = itl, itu
                        ict = ict + 1
                        if(esfbc) jtu = it
                        jct = jctjc
                        do 292 jt = jtl, jtu
                          jct = jct + 1
                          msfbct =  .not. esfb .or.
     &                     (ic.eq.jc .and. it.eq.jt)
                          do 290 ipr = 1, npr
                            icx = icx + 1
                            IPT_CX_TEMP(icx) = 0.0d0
                            if(easfb) then
                              icxa = icxa + 1
                              IPT_CX_TEMP(icxa) = 0.0d0
                            endif
                            ics = mics(ipr)
                            jcs = mjcs(ipr)
                            do 280 ia = 1, mau(ipr)
!                             print *, 'socfpd:',ic,jc,it,jt,ipr,ia
                              ics = ics + 1
                              jcs = jcs + 1
                              IPT_CX_TEMP(icx) = IPT_CX_TEMP(icx)
     &                         + (IPT_AO2SO(ict,ics,igcs) 
     &                           *IPT_AO2SO(jct,jcs,jgcs))
!                             print *,'cx1(',icx,')',IPT_CX_TEMP(icx)
                              if(easfb) IPT_CX_TEMP(icxa) 
     &                                  = IPT_CX_TEMP(icxa)
     &                                  + (IPT_AO2SO(ict,ics,igcs) 
     &                                    *IPT_AO2SO(jct,jcs,jgcs))
                              if(msfbct) go to 280
                              IPT_CX_TEMP(icx) = IPT_CX_TEMP(icx)
     &                         + (IPT_AO2SO(jct,ics,igcs) 
     &                           *IPT_AO2SO(ict,jcs,jgcs))
!                             print *,'cx2(',icx,')',IPT_CX_TEMP(icx)
                              if(easfb) IPT_CX_TEMP(icxa) 
     &                                  = IPT_CX_TEMP(icxa)
     &                                  - (IPT_AO2SO(jct,ics,igcs) 
     &                                    *IPT_AO2SO(ict,jcs,jgcs))
280                         continue
290                       continue
292                     continue
294                   continue
296                 continue
298               continue

310             continue
320           continue
322         continue
324       continue
330     continue
340   continue
 
      CALL RECP_SOCFPD_SHIFTCX(mblu,mcxu,mcxu1,mcxu2,icx,icxu,
     &     ijsf,icxast,ea) 

!     Allocate variables cx(*)
      ALLOCATE (IPT_CX(mcxu))
      DO I = 1, mcxu
         IPT_CX(I) = IPT_CX_TEMP(I)
      ENDDO
!     Deallocate variables cx_temp(*),c(*)
      DEALLOCATE(IPT_CX_TEMP)
      DEALLOCATE(IPT_AO2SO)

      CALL OUTPUT_LOC('ARGOS_SOCFPD','X')
      END

      SUBROUTINE RECP_SOCFPD_SHIFTCX(mblu,mcxu,mcxu1,mcxu2,icx,icxu,
     &           ijsf,icxast,ea)
!       # move antisymmetric coefficient products down.
 
!       #       offsets: 0               mcxu1  mblu
!       # before: a(*)-->................c(*)
!       #                0          mcxu mcxu1  mblu
!       # now:    a(*)-->cx(*)......x(*),c(*)
!       #                0    mcxu2 mcxu       (mblu+mcxu)
!       # after:  a(*)-->cx(*),x(*),........
        USE RECP_IPT
        IMPLICIT NONE
#include "inc_print.h"
!       global variable
        INTEGER mblu,mcxu,mcxu1,mcxu2,icx,icxu,ijsf,icxast(*)
        LOGICAL ea
!       local variable
        INTEGER I,icxl 
      if (ea) then
!       # on return, mcxu2 is the offset for the s.o. coefficients.
        mcxu2 = icx
        icxl  = icx + 1
        icxu  = mcxu - icx
        mcxu  = mcxu1 - icxu
!       print *,'socfpd :',icxl, mcxu
        do icx = icxl, mcxu
          IPT_CX_TEMP(icx) = IPT_CX_TEMP(icx+icxu)
!         print *,'socfpd:cx(',icx,')',IPT_CX_TEMP(icx)
        enddo

!       # make pointers relative to mcxu2 instead of 0.
        do I = 1, ijsf
          icxast(I) = icxast(I) - icxu
        enddo

!       # on return, mcxu and mblu divide the workspace as:
!       # a(*)-->[cx//x](1:mcxu),w(1:mblu)
        mblu = mblu - mcxu
        IF (PRINT_LEVEL.GE.3) write(RECP_OUT,6010) mcxu, mcxu2, mblu
6010    format(/' socfpd: mcxu=',i9,' mcxu2=',i9,' left=',i9)

      else
!       # on return, mcxu and mblu divide the workspace as:
!       # a(*)-->cx(1:mcxu),w(1:mblu)
        mcxu = icx
        mblu = mblu - mcxu
        IF (PRINT_LEVEL.GE.3) write(RECP_OUT,6020) mcxu, mblu
6020    format(/' socfpd: mcxu=',i9,' left=',i9)
      endif
      end

