MODULE RECP_SOCFPD
CONTAINS
SUBROUTINE RECP_SOCFPD_MAIN( mgcsu,icxast,icxst,icxsv1,icxsv2,idp,iprst,la,   &
           lb,maords,mau,mcons,mgcs,mics,mjcs,nblpr,nc,nd,nf,nir,npair,nprir, &  
           nrcr,nsopr,nt,ntl,ntu,CX_TEMP)
! Compute symmetry orbital coefficient products.
! on return from socfpd(*), cx(1:mcxu) was computed.
  USE RECP_FUNCTION1
  USE RECP_FUNCTION2
  USE RECP_INP_ORBITAL
  USE RECP_IPT
  USE RECP_OUTPUT
  USE RECP_WRITE
  implicit logical(a-z)
! il(*) packing factor. see also oneint() and wtint2().
  integer    ilfact
  parameter( ilfact=1024 )
#include "inc_mxvalue.h"
#include "inc_print.h"
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  integer        mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,   inam,   nnam,mdum
  common /ntgr/  mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir(3),inam(5),nnam,mdum(32)
  integer mcxu,mgcsu,CX_TEMP, IPQ_TEMP

  integer icxast(*),icxst(2,*),icxsv1(*),icxsv2(mgcsu,*),idp(mstu,mstu,*),  &
          iprst(*),la(mru,*),lb(*),maords(*),mau(*),mcons(*),  &
          mgcs(*),mics(*),mjcs(*),nblpr(*),nc(*),nd(*),nf(*),nir(*),  &
          npair(2,*),nrcr(*),nsopr(*),nprir(2,mstu,*),nt(*),ntl(*),ntu(*)
! local variables
  integer mcxu1, icx, ipair, ijsf, isf, isfr, igcs, jgcs, is, &
          icu, isfis, isfris, jsf, jsfr, js, jsfjs, jsfrjs, jfu, &
          jcu, if, isfrif, iaords, iru, ircru, itl, itu, jf, jaords, &
          jtl, jsfrjf, jtu, jru, ist, jrcru, ircr, isfrib, jrcr, &
          jsfrjb,iesfb,npr,icxu,icxa,ict,ic,jct,jc,it, &
          jt,I,J,K,FILE_NUM(3),RETURNVALUE
  LOGICAL ea,es,esf,esfb,msfbct
  CHARACTER(10) FILE_NAME(3)
  DATA FILE_NUM  / 64, 65, 66 /
  DATA FILE_NAME / 'RECP_SOC_0','RECP_SOC_1','RECP_SOC_2' /

  CALL OUTPUT_LOC('RECP_SOCFPD_MAIN','E')

! write file in debug mode
  IF (RECP_DBG.GE.1) THEN
     CALL RECP_WRITE_FILEOPEN(FILE_NUM(2),FILE_NAME(2),'REPLACE')
     CALL RECP_WRITE_FILEOPEN(FILE_NUM(3),FILE_NAME(3),'REPLACE')
  ENDIF

  CALL RECP_SETZERO_I3(nprir,2,mstup,mnsfup)

! # check the il(*) packing factor.
  if (nbft.ge.ilfact) then
    call bummer('socfpd: larger ilfact required, nbft=',nbft,2)
  endif
 
  ea = inam(nu).eq.5
! mcxu
  mcxu  = CX_TEMP 
  mcxu1 = mcxu  !original mcxu offset
  icx   = 0
  ipair = 0
! # ijsf = one dimensionalization of is, js, if, jf
! # isfr = one dimensionalization of is, if, ir
  ijsf = 0
  isf  = 0
  isfr = 0

  CALL RECP_SETZERO_I1(icxsv1,ngcs)
  CALL RECP_SETZERO_I2(icxsv2,ngcs,ngcs)

  DO is = 1, ns      ! icenter
     if(nf(is).eq.0) goto 340
     icu = nc(is)
     isfis = isf
     isfris = isfr
     jsf = 0
     jsfr = 0
     DO js = 1,is     ! jcenter
        jsfjs = jsf
        jsfrjs = jsfr
        es = is.eq.js
        if(.not.es) then
          jfu = nf(js)
          if(nf(js).eq.0) goto 330
          jcu = nc(js)
        endif
        isf = isfis
        isfr = isfris
        DO if = 1,nf(is)   ! i-function
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
           DO jf = 1,jfu ! j-function
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
                  do ist = 1,nst
                    nprir(2,ist,ijsf) = 0
                  enddo
                endif
              else
                jrcru = nrcr(mcons(jsf))
                jtu   = ntu(jsf)
                jru   = nir(jaords)
              endif
              iprst(ijsf) = ipair
              isfr = isfrif
              DO ircr = 1,ircru
                 isfrib = isfr
                 jsfr = jsfrjf
                 if(esf) jrcru = ircr
                 DO jrcr = 1,jrcru
                    jsfrjb = jsfr
                    esfb = esf.and.ircr.eq.jrcr

                    CALL RECP_SOCFPD_SETNPR(nst,npr,isfr,jsfr,isfrib, &
                         jsfrjb,iru,jru,iaords,jaords,ipair,ijsf,mru,la, &
                         lb,mstu,idp,nsopr,mau,mics,mjcs,nd,nprir,npair,ilfact,esfb)
                  
                    IF (npr.NE.0) THEN 
                       CALL RECP_SOCFPD_SETICXU(RETURNVALUE,igcs,jgcs,ijsf,icx, &
                       IPQ_TEMP,icu,icxu,isf,jsf,npr,mcxu,icxa,mgcsu, &
                       icxsv1,icxsv2,icxst,icxast,nt,jcu,esfb,ea)
                       IF (RETURNVALUE.EQ.0) THEN
!                         compute new cx block
                          CALL RECP_SOCFPD_CALCX0(icu,jcu,itl,itu,jtl,jtu,npr,icx, &
                               icxa,mics,mjcs,mau,igcs,jgcs,FILE_NUM,ea,msfbct,esfb)
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO     
        330 CONTINUE
     ENDDO
     340 CONTINUE
  ENDDO 

  CALL RECP_SOCFPD_SHIFTCX(mcxu,mcxu1,icx,icxu,ijsf,icxast,ea) 

! Allocate variables cx(*)
  ALLOCATE (IPT_CX(mcxu))
  DO I = 1, mcxu
     IPT_CX(I) = IPT_CX0(I)
  ENDDO
! Deallocate variables cx_temp(*),c(*)
  DEALLOCATE(IPT_CX0)
  DEALLOCATE(IPT_AO2SO)

  IF (RECP_DBG.GE.1) THEN
     CLOSE(FILE_NUM(2))
     CLOSE(FILE_NUM(3))
  ENDIF
  CALL OUTPUT_LOC('RECP_SOCFPD_MAIN','X')
END SUBROUTINE RECP_SOCFPD_MAIN
END MODULE RECP_SOCFPD
