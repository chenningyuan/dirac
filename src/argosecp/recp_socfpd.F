      subroutine recp_socfpd( mgcsu,icxast,icxst,
     &   icxsv1, icxsv2, idp,  iprst,  la,
     &   lb,     maords, mau,  mcons,      mgcs,   mics,
     &   mjcs,   nblpr,  nc,   nd,         nf,     nir,
     &   npair,  nprir,  nrcr, nsopr,      nt,     ntl,
     &   ntu, CX_TEMP)
!    Compute symmetry orbital coefficient products.
!    on return from socfpd(*), cx(1:mcxu) was computed.

      USE RECP_INP_ORBITAL
      USE RECP_IPT
      USE RECP_OUTPUT
      USE RECP_WRITE
      implicit logical(a-z)

!     il(*) packing factor. see also oneint() and wtint2().
      integer    ilfact
      parameter( ilfact=1024 )

#include "inc_mxvalue.h"
#include "inc_print.h"

      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      integer
     & mconsu, mru,   mcru, msfu, 
     & msfru,  ngcs,  nu,   
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & mconsu, mru,   mcru, msfu, 
     & msfru,  ngcs,  nu,   
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

!     # local:
      integer mcxu1, icx, ipair, ijsf, isf, isfr, igcs, jgcs, is,
     & ifu, icu, isfis, isfris, jsf, jsfr, js, jsfjs, jsfrjs, jfu,
     & jcu, if, isfrif, iaords, iru, ircru, itl, itu, jf, jaords,
     & jtl, jsfrjf, jtu, jru, ist, jrcru, ircr, isfrib, jrcr,
     & jsfrjb, iesfb, npr, ipr, ics, ir, lai, jcs, jr, laj,
     & isv1, icxu, icxa, isv2, ict, ic, ictic, jct, jc, jctjc, it,
     & jt,npru, I,J,K,FILE_NUM(3)
      LOGICAL ea, easfb, es, esf, esfb, esfbc, msfbct
      CHARACTER(10) FILE_NAME(3)
      DATA FILE_NUM  / 64, 65, 66 /
      DATA FILE_NAME / 'RECP_SOC_0','RECP_SOC_1','RECP_SOC_2' /

      CALL OUTPUT_LOC('ARGOS_SOCFPD','E')

!     write file in debug mode
      IF (RECP_DBG.GE.1) THEN
         CALL RECP_WRITE_FILEOPEN(FILE_NUM(2),FILE_NAME(2),'REPLACE')
         CALL RECP_WRITE_FILEOPEN(FILE_NUM(3),FILE_NAME(3),'REPLACE')
      ENDIF

      DO I = 1, 2
         DO J = 1, mstup 
            DO K = 1, mnsfup
               nprir(I,J,K) = 0
            ENDDO
         ENDDO
      ENDDO
!     # check the il(*) packing factor.
      if ( nbft .ge. ilfact ) then
        call bummer('socfpd: larger ilfact required, nbft=',nbft,2)
      endif
 
      ea = inam(nu).eq.5

!     mcxu,mblu 
      mblu   = 4500000  !temporary: total available workspace
      mcxu = mblu

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
 
                  if (npru.gt.mnrup) then
                    call bummer('change mnrup (one place) to ',npru,2)
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
     &               (icxu-mcxu),2)
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
                          msfbct =.not.esfb .or.
     &                     (ic.eq.jc .and. it.eq.jt)

                          CALL RECP_SOCFPD_CX0(npr,icx,icxa,mics,mjcs,
     &                     mau,ict,jct,igcs,jgcs,FILE_NUM,easfb,msfbct)
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
         IPT_CX(I) = IPT_CX0(I)
      ENDDO
!     Deallocate variables cx_temp(*),c(*)
      DEALLOCATE(IPT_CX0)
      DEALLOCATE(IPT_AO2SO)

      IF (RECP_DBG.GE.1) THEN
         CLOSE(FILE_NUM(2))
         CLOSE(FILE_NUM(3))
      ENDIF
      CALL OUTPUT_LOC('ARGOS_SOCFPD','X')
      END
