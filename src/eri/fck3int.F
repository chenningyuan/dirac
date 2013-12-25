!dalton_copyright_start
!      Copyright (c) 2013 by the authors of Dalton (see below).
!      All Rights Reserved.
!   
!      The source code in this file is part of
!      "Dalton, a molecular electronic structure program,
!       Release DALTON2013 (2013), see http://daltonprogram.org"
!   
!      This source code is provided under a written licence and may be
!      used, copied, transmitted, or stored only in accord with that
!      written licence.
!   
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the licence holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the licence holder may not take any part of Dalton,
!      including modified files, with him/her, unless that person has
!      obtained his/her own licence.
!   
!      For further information, including how to get a licence, see:
!         http://daltonprogram.org
!dalton_copyright_end

      SUBROUTINE GETINT3(RMAT,IODBCLS,IODBTCH,RODBTCH,
     &                   WORK,LWORK,IPRINT, AUXEL1, AUXEL2, 
     &                   I_LG, NDIMAUX, DOGAB)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3int.h"
CAREFUL
#include "erithr.h"
CAREFUL
      DIMENSION IODBCLS(5,*), IODBTCH(2,*),
     &     RODBTCH(7,*), WORK(*), RMAT(*), I_LG(*)
      LOGICAL GDER, AUXEL1, AUXEL2, DOGAB
      CALL QENTER('GETINT3')
CAREFUL
      THRSHSAVE  = THRSH
      THRSH      = 1.0D-25
CAREFUL

      CALL SETINT(AUXEL1, AUXEL2, DOGAB, NDIMAUX)

C     MAXDER for derivatives only, not implemented here
      MAXDER = 0 
      GDER   = .FALSE. 

      MAXKCK = MAXAONHK * (MAXAONHK + 1) / 2
      MAXNHK = MAXAONHK

      CALL PUTAO2ERICOM(MAXNHK,MAXKCK)

      JTOP  = 4*(MAXAONHK - 1)
      KC2MAX = MAXKCK**2
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
      NRDER  = 3*MAXDER

      CALL PUT2HERTOP(JTOP,NRTOP)

C     Dimensions for CSQ
C     ==================
C
      CALL SPHDIM
C
C     Allocate
C     ========
C
      KINDHR  = 1
      KINDSQ  = KINDHR  + ( JTOP3 - 1 ) / IRAT + 1
      KIODHR  = KINDSQ  + ( NRTOP - 1 ) / IRAT + 1
      KLMNVL  = KIODHR  + ( 8 * NRTOP - 1) / IRAT + 1
      KLMNSM  = KLMNVL  + ( 3 * MAXKCK * MAXNHK - 1) / IRAT + 1
      KCSQ    = KLMNSM  + (16 * MAXKCK * MAXNHK - 1) / IRAT + 1
      KIPNTUV = KCSQ    + NCSQ1 * NCSQ2
      KLAST   = KIPNTUV + (4 * KC2MAX * ( NRDER + 1 ) -1) / IRAT + 1 

C      write(LUPRI,*) 'Memory requirement in GETINT3 ', KLAST

      IF (KLAST .GT. LWORK) CALL STOPIT('GETINT3','HERPRP',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1

C     Spherical-harmonic transformation matrices
C     ==========================================
C

      CALL SPHMAT(WORK(KCSQ),NCSQ1,NCSQ2,MAXDER,GDER,WORK(KLAST),LWRK,
     &            IPRINT)

C
C     Arrays for Hermite integrals
C     ============================
C
      CALL HERPRP(WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR))
      CALL LMNPRP(WORK(KLMNVL),WORK(KLMNSM),IPRINT)
C
C     Main two-electron subroutine
C     ============================

      CALL GETINT31( RMAT,IODBCLS,IODBTCH,RODBTCH,
     &               WORK(KINDHR),WORK(KIODHR),WORK(KINDSQ),
     &               WORK(KLMNVL),WORK(KLMNSM),WORK(KCSQ),
     &               WORK(KIPNTUV),WORK(KLAST),LWRK,
     &               I_LG,NDIMAUX,IPRINT)
CAREFUL
      THRSH = THRSHSAVE
CAREFUL

      CALL QEXIT('GETINT3')
      RETURN
      END

      SUBROUTINE GETINT31(RMAT,IODBCLS,IODBTCH,RODBTCH,
     &                   INDHER,IODDHR,INDHSQ,LMNPWR,LMNSYM,
     &                   CSQ,IPNTUV,WORK,LWORK,I_LG,NDIMAUX,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "fck3mem.h"
#include "cbifck3.h"
#include "clsfmm.h"
#include "densfit.h"
#include "fck3int.h"
#include "dummy.h"
cSimen
#include "nuclei.h"
cSimen
      DIMENSION IAB0(3), ICD0(3), 
     &          IODBCLS(5,*), IODBTCH(2,*), RODBTCH(7,*),
     &          CSQ(NCSQ1,NCSQ2), RMAT(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHSQ(*),
     &          IODDHR(NRTOP,0:7),LMNPWR(MAXAOKCK,MAXAONHK,3),
     &          LMNSYM(MAXAOKCK,MAXAONHK,2,0:7),WORK(*), 
     &          KCCF(2,4), IODDKC(MAXAOKCK**2), IPNTUV(*), I_LG(*)
      LOGICAL GCNTA, GCNTB, GCNTC, GCNTD, GCNTAB, GCNTCD, 
     &        DODIST, GENCON, BOXACT
#include "ibtfun.h"
      CALL QENTER('GETINT31')

C     Monitor maximum memory requirement
      MAXLAST = 0

c     CALL INIT31(GENCON,NPQBCXFCK3,MAXNPQBCX,MAXDER,MDIMD,
c    &            NDER,NRDER,DODIST,MAXKCK,KC2MAX)

C     Set GENCON according to ERI      
      GENCON=.FALSE.

cSimen Number of density matrices
c     NDMAT      = 1
cSimen Because we have not implemented symmetry:
c     IHCBITFCK3 = 0
      MLTPZFCK3  = 1
      NPQBCSFCK3 = 1
      MAXNPQBCX  = 1
      NPQBCXFCK3 = MLTPZFCK3 * NPQBCSFCK3
cSimen Because we have no derivatives
      MAXDER     = 0
      MDIMD      = 1
      NDER       = 0
      NRDER      = 3 * MAXDER
cSimen Because we use batches directly, not distributions
      DODIST     = .FALSE.


      MAXKCK = MAXAONHK * (MAXAONHK + 1) / 2
      KC2MAX = MAXKCK**2

C     Initializes IABCLSSTART, IABCLSEND, ICDCLSSTART and ICDCLSEND
      CALL INILOOP

      IAOAUXSTART = 1
      IF (AUX1) THEN
        IAOAUXSTART = KAOBTCH(1,IODBTCH(2,IODBCLS(4,IABCLSSTART)))
      ELSE IF (AUX2) THEN
        IAOAUXSTART = KAOBTCH(1,IODBTCH(2,IODBCLS(4,ICDCLSSTART)))
      ENDIF

C     1st electron
      DO IABCLS = IABCLSSTART, IABCLSEND
        ICENTA    = IODBCLS(1,IABCLS)
        ICENTB    = IODBCLS(2,IABCLS)
        NHKTAB    = IODBCLS(3,IABCLS)
        KABSTBTCH = IODBCLS(4,IABCLS)
        NABBTCH   = IODBCLS(5,IABCLS)
        
        JMAXAB = NHKTAB - 1
        NTUVAB = NHKTAB*(NHKTAB+1)*(NHKTAB+2)/6
       
C       Distances between A and B
        CALL GETDIST( 1, ICENTA, ICENTB, NSYMOPR, NSYMOPS, 1, 
     &       IAB0, IABXYZ, AUX1)
        DO IABBTCH = KABSTBTCH, ( KABSTBTCH + NABBTCH - 1 ) 
         IAOBA = IODBTCH(1,IABBTCH)
         IAOBB = IODBTCH(2,IABBTCH)
C        Only enter when AB-batch is activated
         BOXACT = (.NOT.DOBOX).OR.GAB.OR.
     &            (ACTVAOBTCH(IAOBA).AND.ACTVAOBTCH(IAOBB))
         IF ((ACTVBTH(IAOBA,2).AND.ACTVBTH(IAOBB,1)).AND.BOXACT) THEN
          CALL GETAOB(NPRFA,NCTFA,NTYPA,NHKTA,JMAXA,KHKTA,KCKTA,
     &                KCSQA,KEXPA,KCCF(1,1),IAOBA)

          CALL GETAOB(NPRFB,NCTFB,NTYPB,NHKTB,JMAXB,KHKTB,KCKTB,
     &                KCSQB,KEXPB,KCCF(1,2),IAOBB)

          IAOA = KAOBTCH(1,IAOBA)
          IAOB = KAOBTCH(1,IAOBB)

          NPRFAB = NPRFA * NPRFB
          JMAXAB = JMAXA + JMAXB
          NTUVAB = ( (JMAXAB+1)*(JMAXAB+2)*(JMAXAB+3) ) / 6
          NCTFAB = NCTFA * NCTFB
          KCKTAB = KCKTA * KCKTB 
          KHKTAB = KHKTA * KHKTB
          TKMPAB = (IABXYZ.EQ.7).AND.
     &         (KHKTA.EQ.KHKTB).AND.(.NOT.DODIST)
          TCMPAB = TKMPAB .AND. .NOT. SPHERAB

          IF (TKMPAB) KHKTAB = KHKTA*(KHKTA + 1)/2
          IF (TCMPAB) KCKTAB = KCKTA*(KCKTA + 1)/2

          GCNTA  = GENCON .OR. NCTFA  .GT. 1
          GCNTB  = GENCON .OR. NCTFB  .GT. 1
          GCNTAB = GENCON .OR. NCTFAB .GT. 1


          SPHERA  = KHKTA .NE. KCKTA .OR. MAXDER .GT. 0
          SPHERB  = KHKTB .NE. KCKTB .OR. MAXDER .GT. 0
          SPHERAB = SPHERA .OR. SPHERB

C         Memory allocations for outer loop (1st electron)
          KPREFAB   = 1
          KIPRAB    = KPREFAB + NPRFAB
          KDATAP    = KIPRAB  + NPRFAB * 2 / IRAT + 1 
          KETIJAB   = KDATAP  + NPRFAB * 6 
          KLAST     = KETIJAB + NPRFAB * NHKTAB * NHKTA * NHKTB * 6
          
          IF (KLAST.GT.LWORK) 
     &         CALL STOPIT('GETINT3','GETINT31a',KLAST,LWORK)     
      
          CALL SETIPR(WORK(KIPRAB), NPRFA, NPRFB)
          CALL SETP(WORK(KDATAP),WORK(KIPRAB),
     &         COOR(1,0,1),COOR(1,0,2),EXPBTCH(KEXPA),
     &         EXPBTCH(KEXPB),NPRFA,NPRFB,NPRFAB,IPRINT)


C         Temporary allocations
          KCOORAB   = KLAST
          KEXPAB    = KCOORAB + NPRFAB * 3
          KHPIAB    = KEXPAB  + NPRFAB * 3
          KPA       = KHPIAB  + NPRFAB
          KPB       = KPA     + NPRFAB * 3
          KLAST1    = KPB     + NPRFAB * 3 
          IF (KLAST1.GT.LWORK) 
     &         CALL STOPIT('GETINT3','GETINT31b',KLAST1,LWORK)     
          MAXLAST  = MAX(MAXLAST,KLAST1)

          CALL SETETIJ(WORK(KETIJAB),WORK(KPREFAB),
     &             DIST12(1,0,0,1),DIST122(0,0,1),
     &             WORK(KIPRAB),EXPBTCH(KEXPA),EXPBTCH(KEXPB),
     &             WORK(KCOORAB),WORK(KEXPAB),WORK(KHPIAB),WORK(KPA),
     &             WORK(KPB),NPRFA,NHKTA,JMAXA,
     &             NPRFB,NHKTB,JMAXB,NPRFAB,NHKTAB,1,IPRINT)

          IF (GAB) THEN
            ICDCLSSTART = IABCLS
            ICDCLSEND   = IABCLS
          ELSEIF (ODTR12) THEN
            ICDCLSEND   = IABCLS
          END IF

C         2nd electron
          DO ICDCLS = ICDCLSSTART, ICDCLSEND
            ICENTC    = IODBCLS(1,ICDCLS)
            ICENTD    = IODBCLS(2,ICDCLS)
            NHKTCD    = IODBCLS(3,ICDCLS)
            KCDSTBTCH = IODBCLS(4,ICDCLS)
            NCDBTCH   = IODBCLS(5,ICDCLS)
            
            JMAXCD = NHKTCD - 1
            NTUVCD = NHKTCD*(NHKTCD+1)*(NHKTCD+2)/6
            
C           Distances between C and D
            CALL GETDIST( 2, ICENTC, ICENTD, NSYMOPR, NSYMOPS, 1, 
     &       ICD0, ICDXYZ, AUX2)
             
            IF (GAB) THEN
              ICDBTCHSTART = IABBTCH
              ICDBTCHEND   = IABBTCH
            ELSEIF (ODTR12 .AND. IABCLS .EQ. ICDCLS) THEN
              ICDBTCHSTART = KCDSTBTCH
              ICDBTCHEND   = IABBTCH
            ELSE
              ICDBTCHSTART = KCDSTBTCH
              ICDBTCHEND   = KCDSTBTCH + NCDBTCH - 1
            ENDIF

C           ************* INNER LOOP STARTS HERE ***************
            DO ICDBTCH = ICDBTCHSTART, ICDBTCHEND
             IAOBC = IODBTCH(1,ICDBTCH)
             IAOBD = IODBTCH(2,ICDBTCH)
C            Only enter when CD-batch is activated
             BOXACT = (.NOT.DOBOX).OR.GAB.OR.
     &                (ACTVAOBTCH(IAOBC).AND.ACTVAOBTCH(IAOBD))
             IF ((ACTVBTH(IAOBC,4).AND.ACTVBTH(IAOBD,3))
     &           .AND.BOXACT) THEN
              CALL GETAOB(NPRFC,NCTFC,NTYPC,NHKTC,JMAXC,KHKTC,KCKTC,
     &                    KCSQC,KEXPC,KCCF(1,3),IAOBC)
              CALL GETAOB(NPRFD,NCTFD,NTYPD,NHKTD,JMAXD,KHKTD,KCKTD,
     &                    KCSQD,KEXPD,KCCF(1,4),IAOBD)

              IAOC = KAOBTCH(1,IAOBC)
              IAOD = KAOBTCH(1,IAOBD)

              NPRFCD = NPRFC * NPRFD
              NPRF   = NPRFAB*NPRFCD
              NTYP   = NTYPA + NTYPB + NTYPC + NTYPD
              KCKTCD = KCKTC * KCKTD
              KHKTCD = KHKTC * KHKTD
              TKMPCD = (ICDXYZ.EQ.7).AND.
     &             (KHKTC.EQ.KHKTD).AND.(.NOT.DODIST)
              TCMPCD = TKMPCD .AND. .NOT. SPHERCD
              IF (TKMPCD) KHKTCD = KHKTC*(KHKTC + 1)/2
              IF (TCMPCD) KCKTCD = KCKTC*(KCKTC + 1)/2
              JMAX   = JMAXAB + JMAXCD
              NTUV   = ( (JMAX+1)*(JMAX+2)*(JMAX+3) ) / 6
              NCTFCD = NCTFC * NCTFD
              NCCPP   = NPQBCXFCK3 * NPRFAB * NCTFCD
              NCCCC   = NPQBCXFCK3 * NCTFAB * NCTFCD
              NCCT    = NPQBCXFCK3 * NCTFAB * NCTFCD
              NPPS    = NPQBCXFCK3 * NPRFAB * NPRFCD
              GCNTC  = GENCON .OR. NCTFC .GT. 1
              GCNTD  = GENCON .OR. NCTFD .GT. 1
              GCNTCD = GENCON .OR. NCTFCD .GT. 1
              SPHERC  = KHKTC .NE. KCKTC .OR. MAXDER .GT. 0
              SPHERD  = KHKTD .NE. KCKTD .OR. MAXDER .GT. 0
              SPHERCD = SPHERC .OR. SPHERD

C             Allocations  for the inner loop 
              KPREFCD = KLAST
              KIPRCD  = KPREFCD + NPRFCD 
              KETIJCD = KIPRCD  + NPRFCD * 2 / IRAT + 1 
              KFACINT = KETIJCD + NPRFCD * NHKTCD * NHKTC * NHKTD * 6
              KALPHA  = KFACINT + NPRFAB * NPRFCD
              KCOORPQ = KALPHA  + NPRFAB * NPRFCD 
              KLAST1  = KCOORPQ + NPRFAB * NPRFCD * 3
              IF (KLAST1.GT.LWORK) 
     &             CALL STOPIT('GETINT3','GETINT31c',KLAST1,LWORK)     

              CALL SETIPR(WORK(KIPRCD), NPRFC, NPRFD)
              CALL SETQR( WORK(KFACINT), WORK(KALPHA), WORK(KCOORPQ), 
     &                    WORK(KDATAP), WORK(KIPRAB), WORK(KIPRCD), 
     &                    COOR(1,0,3), COOR(1,0,4), EXPBTCH(KEXPC),
     &                    EXPBTCH(KEXPD),NPRFC,NPRFD,NPRFCD,
     &                    NPRF, IPQXYZFCK3, NPRFAB, 
     &                    JMAX, NTYP, NTUV, CCFBTCH(KCCF(1,1)), 
     &                    CCFBTCH(KCCF(1,2)), CCFBTCH(KCCF(1,3)), 
     &                    CCFBTCH(KCCF(1,4)), GCNTAB, GCNTCD)

              IHCBITFCK3 = IBTAND(IPQXYZFCK3,ICDXYZ)
              ICCXYZFCK3 = IBTAND(IPQXYZFCK3,IBTAND(IABXYZ,ICDXYZ))
              CALL PUT2ERI(NPRFA, NPRFB, NPRFC, NPRFD, 
     &                     NCTFA, NCTFB, NCTFC, NCTFD, 
     &                     KCKTA, KCKTB, KCKTC, KCKTD, 
     &                     KHKTA, KHKTB, KHKTC, KHKTD, 
     &                     NHKTA, NHKTB, NHKTC, NHKTD, 
     &                     NPQBCXFCK3, NPRFAB, NPRFCD, 
     &                     NCTFA*NCTFB, NCTFC*NCTFD, NPRFAB, NPRFCD, 
     &                     JMAXA, JMAXB, JMAXC, JMAXD,
     &                     NTUVAB, KCKTAB, KCKTCD, KHKTAB, KHKTCD,
     &                     KC2MAX, MAXAONHK, NRDER, SPHERA, SPHERB, 
     &                     SPHERC, SPHERD, TKMPAB, TKMPCD, 
     &                     GCNTA, GCNTB, GCNTC, GCNTD, 
     &                     ICD0, IAB0, ICDXYZ, IABXYZ)

              
              CALL PNTUV(IPNTUV, LMNPWR, INDHER, IPRINT)

C             Temporary allocations for SETETIJ
              KCOORCD   = KLAST1
              KEXPCD    = KCOORCD + NPRFCD * 3
              KHPICD    = KEXPCD  + NPRFCD * 3
              KPC       = KHPICD  + NPRFCD
              KPD       = KPC     + NPRFCD * 3
              KLAST2    = KPD     + NPRFCD * 3
              IF (KLAST2.GT.LWORK) 
     &             CALL STOPIT('GETINT3','GETINT31d',KLAST2,LWORK)     
              MAXLAST  = MAX(MAXLAST,KLAST2)

              CALL SETETIJ(WORK(KETIJCD),WORK(KPREFCD),
     &             DIST12(1,0,0,2),DIST122(0,0,2),
     &             WORK(KIPRCD),EXPBTCH(KEXPC),EXPBTCH(KEXPD),
     &             WORK(KCOORCD),WORK(KEXPCD),WORK(KHPICD),WORK(KPC),
     &             WORK(KPD),NPRFC,NHKTC,JMAXC,
     &             NPRFD,NHKTD,JMAXD,NPRFCD,NHKTCD,2,IPRINT)

C             Multiply FACINT with preexponetial factors to concur with ERI
              CALL FACEXP(WORK(KFACINT), WORK(KPREFAB), WORK(KPREFCD),
     &                    NPRFAB, NPRFCD)

C             Allocations for both GCDTUV and G_ABCD
              KHCINT     = KLAST1
              KHCCONT    = KHCINT + NCCPP * NTUVAB * KCKTCD * MAXNPQBCX
              KLAST2     = KHCCONT + NCCPP * NTUVAB * KCKTCD

C             Temporary allocations
              KHRINT1    = KLAST2
              KLAST3     = KHRINT1 + NPPS * NTUV
              KRJ000     = KLAST3
              KLAST4     = KRJ000 + NPPS * (JMAX + 1)
              IF (KLAST4.GT.LWORK) 
     &             CALL STOPIT('GETINT3','GETINT31e',KLAST4,LWORK)     
              MAXLAST = MAX(MAXLAST,KLAST4)

              LWRK    = LWORK - KLAST4
C             Code for Gaussian damped integrals; over the operator
C             exp(-a*r12*r12)/r12 instead of 1/r12. Turned off. \SR
              IF (DF_GAUDMP.AND.(AUX1.AND.AUX2.AND.(.NOT.GAB))) THEN
                NUCABCDT = NPRFAB * NPRFCD + 1000
                KHARGE = KLAST4
                KHARGF = KHARGE + NUCABCDT
                KHALPH = KHARGF + NUCABCDT
                KHBETA = KHALPH + NUCABCDT
                KHPRE1 = KHBETA + NUCABCDT
                KHPRE2 = KHPRE1 + NUCABCDT
                KHPRE3 = KHPRE2 + NUCABCDT
                KRO000 = KHPRE3 + NUCABCDT
                KRM000 = KRO000 + NUCABCDT * (JMAX + 1)
                KHRIN3 = KRM000 + NUCABCDT * (JMAX + 1)
                KHRIN4 = KHRIN3 + NUCABCDT * NTUV
                KHRIN5 = KHRIN4 + NUCABCDT * NTUV
                KHEXPP = KHRIN5 + NUCABCDT * NTUV
                KHEXPQ = KHEXPP + NUCABCDT
                KIPRCD = KHEXPQ + NUCABCDT
                KLASTT = KIPRCD  + NPRFCD * 2 / IRAT + 1 + 1000

                IF (KLASTT .GT. LWORK) CALL STOPIT('FCK3','FCK31'
     &                                             ,KLASTT,LWORK)     
                LWRKT   = LWORK - KLASTT + 1

                CALL SETEXPPQ(WORK(KIPRCD),NPRFC,NPRFD,NPRFAB,NPRFCD,
     &               EXPBTCH(KEXPC),EXPBTCH(KEXPD),WORK(KHEXPP),
     &               WORK(KHEXPQ),WORK(KDATAP))

                CALL ERIHER(WORK(KHRINT1),dummy,WORK(KRM000),
     &               WORK(KCOORPQ),INDHER,IODDHR(1,IPQXYZFCK3),
     &               WORK(KFACINT),WORK(KALPHA),WORK(KHARGE),
     &               WORK(KHARGF),WORK(KHALPH),WORK(KHBETA),
     &               WORK(KHPRE1),WORK(KHPRE2),WORK(KHPRE3),
     &               WORK(KRO000),WORK(KHRIN3),WORK(KHRIN4),
     &               WORK(KHRIN5),WORK(KHEXPP),WORK(KHEXPQ),
     &               WORK(KLASTT),LWRKT,IPRINT)
              ELSE
                CALL ERIGAM(WORK(KRJ000),WORK(KFACINT),WORK(KALPHA),
     &                      WORK(KCOORPQ),WORK(KLAST4),LWRK,IPRINT)

                CALL ERITUV(WORK(KHRINT1),WORK(KRJ000),
     &                      WORK(KCOORPQ),INDHER,IODDHR(1,IPQXYZFCK3),
     &                      WORK(KLAST4),LWRK,IPRINT)
              ENDIF
#if 0
              DO ICOMPC = 1, NTYPC
              IF (NTYPC.EQ.2) CALL SETJKKK(JC,KCSQC,KCKTC,KHKTC,ICOMPC)
              DO ICOMPD = 1, NTYPD
              IF (NTYPD.EQ.2) CALL SETJKKK(JD,KCSQD,KCKTD,KHKTD,ICOMPD)
C KHKTCD, KCKTCD, NTUVAB, NCCPP, NPRFCD, NCCCC, NCCT, KHKTAB
C KCKTAB, NCTFCD, NPRFAB,IAOA,IAOB,IAOC,IAOD,
#endif


              KCPINT  = KLAST3
              KFCINT  = KCPINT  + NCTFC * NPRFD * MAXNPQBCX
              KHCSINT = KFCINT  + NTUVAB * KHKTCD
              KFCSINT = KHCSINT + NCCPP * NTUVAB * KHKTD
              KHCPRIM = KFCSINT + NTUVAB * NHKTD
              KEUV    = KHCPRIM + NPRFCD * NPRFAB * NTUVAB
              KETUV   = KEUV    + NPRFCD
              KFRSTUV = KETUV   + NPRFCD
              KLAST4  = KFRSTUV + NTUVAB
              IF (KLAST4.GT.LWORK) 
     &             CALL STOPIT('GETINT3','GCDTUV',KLAST4,LWORK)     
              MAXLAST = MAX(MAXLAST,KLAST4)

C            Generates G^CD_tuv in spherical AO/contracted representation.
              CALL GCDTUV(WORK(KHCINT), KHKTCD, KCKTCD, NHKTC, NHKTD, 
     &             KCKTC, KCKTD, LMNPWR, IODDKC, 
     &             IODDHR, INDHER, WORK(KHRINT1), WORK(KHCPRIM), 
     &             WORK(KETIJCD), WORK(KEUV), WORK(KETUV), 
     &             WORK(KFRSTUV), WORK(KHCCONT), CCFBTCH(KCCF(1,3)), 
     &             CCFBTCH(KCCF(1,4)), WORK(KCPINT), INDHSQ,
     &             NCCPP, NTUVAB, WORK(KFCINT), WORK(KHCSINT), 
     &             WORK(KFCSINT), CSQ(KCSQC,1), CSQ(KCSQD,1), IPNTUV,
     &             GCNTCD, GCNTAB, NPRFCD, SPHRCD, IPRINT)

#if 0
              DO ICOMPA = 1, NTYPA
              IF (NTYPA.EQ.2) CALL SETJKKK(JA,KCSQA,KCKTA,KHKTA,ICOMPA)
              DO ICOMPB = 1, NTYPB
              IF (NTYPB.EQ.2) CALL SETJKKK(JB,KCSQB,KCKTB,KHKTB,ICOMPB)
#endif

              KAOINT  = KLAST2
              KFAOINT = KAOINT  + NCCT*KHKTAB*KHKTCD
              KCSVA   = KFAOINT + KHKTAB * KHKTCD
              KCSVB   = KCSVA   + NCCCC * KHKTA * KCKTA
              KHCCONT = KCSVB   + NCCCC * KHKTB * KCKTB 
              KCPINT  = KHCCONT + NCCCC * KHKTCD * KCKTAB * MDIMD
              KCSINT  = KCPINT  + NPRFB * NCTFA * NCTFCD * MAXNPQBCX
              KFCSINT = KCSINT  + NCCCC * KHKTCD * KCKTB
              KHCPRIM = KFCSINT + KHKTCD * KCKTB 
              KFRSTUV = KHCPRIM + MDIMD * NCCPP * KHKTCD
              KEUV    = KFRSTUV + KHKTCD
              KETUV   = KEUV    + NPRFAB
              KLAST3  = KETUV   + NPRFAB
              IF (KLAST3.GT.LWORK) CALL 
     &             STOPIT('GETINT31','G_ABCD',KLAST3,LWORK)     
              MAXLAST = MAX(MAXLAST,KLAST3)
              
              CALL DZERO(WORK(KAOINT),NCCT*KHKTAB*KHKTCD)

              CALL G_ABCD(WORK(KAOINT), WORK(KFAOINT), WORK(KHCINT),
     &             INDHER, IODDHR, IODDKC, WORK(KHCPRIM), WORK(KFRSTUV),
     &             WORK(KHCCONT), WORK(KCSINT), WORK(KFCSINT), 
     &             WORK(KETIJAB), WORK(KEUV), WORK(KETUV), LMNPWR, 
     &             IPNTUV, CSQ(KCSQA,1), CSQ(KCSQB,1), 
     &             WORK(KCPINT), MDIMD, KC2MAX, 
     &             NRDER, NCCCC,  NCCT, IABXYZ, KHKTAB, KHKTCD, KCKTAB, 
     &             NHKTA, NHKTB, KCKTA, KCKTB, WORK(KCSVA), WORK(KCSVB),
     &             CCFBTCH(KCCF(1,1)), CCFBTCH(KCCF(1,2)), IPRINT)
              
              KI_GL = KLAST3
              IF (DOBOX .AND. .NOT. GAB) THEN
               KLAST3 = KI_GL + NBASISAUX
               IF (KLAST3.GT.LWORK) CALL
     &             STOPIT('GETINT31','G_ABCD',KLAST3,LWORK)
               MAXLAST = MAX(MAXLAST,KLAST3)
               CALL INDEX_GLOBAL_LOC(WORK(KI_GL),NBASISAUX,I_LG,NDIMAUX)
              ENDIF

              CALL F_ABCD(RMAT,WORK(KAOINT),NCCT, 
     &             NCTFA,NCTFB,NCTFC,NCTFD,
     &             NHKTA,NHKTB,NHKTC,NHKTD, 
     &             KHKTA,KHKTB,KHKTC,KHKTD,
     &             IAOA,IAOB,IAOC,IAOD,
     &             KHKTAB,KHKTCD,IAOAUXSTART,
     &             WORK(KI_GL),IPRINT)

#if 0
               ENDDO
               ENDDO
               ENDDO
               ENDDO
#endif
C             Temporary allocations
c             KCPINT  = KLAST3
c             KFCINT  = KCPINT  + NCTFC * NPRFD * MAXNPQBCX
c             KHCSINT = KFCINT  + NTUVAB * KHKTCD
c             KFCSINT = KHCSINT + NCCPP * NTUVAB * KHKTD
c             KHCPRIM = KFCSINT + NTUVAB * NHKTD
c             KEUV    = KHCPRIM + NPRFCD * NPRFAB * NTUVAB
c             KETUV   = KEUV    + NPRFCD
c             KFRSTUV = KETUV   + NPRFCD
c             KLAST4  = KFRSTUV + NTUVAB
c             IF (KLAST4.GT.LWORK) 
c    &             CALL STOPIT('GETINT3','GCDTUV',KLAST4,LWORK)     
c             MAXLAST = MAX(MAXLAST,KLAST4)

c            Generates G^CD_tuv in spherical AO/contracted representation.
c             CALL GCDTUV(WORK(KHCINT), KHKTCD, KCKTCD, NHKTC, NHKTD, 
c    &             KCKTC, KCKTD, LMNPWR, IODDKC, 
c    &             IODDHR, INDHER, WORK(KHRINT1), WORK(KHCPRIM), 
c    &             WORK(KETIJCD), WORK(KEUV), WORK(KETUV), 
c    &             WORK(KFRSTUV), WORK(KHCCONT), CCFBTCH(KCCF(1,3)), 
c    &             CCFBTCH(KCCF(1,4)), WORK(KCPINT), INDHSQ,
c    &             NCCPP, NTUVAB, WORK(KFCINT), WORK(KHCSINT), 
c    &             WORK(KFCSINT), CSQ(KCSQC,1), CSQ(KCSQD,1), IPNTUV,
c    &             GCNTCD, GCNTAB, NPRFCD, SPHRCD, IPRINT)

c             KAOINT  = KLAST2
c             KFAOINT = KAOINT  + NCCT*KHKTAB*KHKTCD
c             KCSVA   = KFAOINT + KHKTAB * KHKTCD
c             KCSVB   = KCSVA   + NCCCC * KHKTA * KCKTA
c             KHCCONT = KCSVB   + NCCCC * KHKTB * KCKTB 
c             KCPINT  = KHCCONT + NCCCC * KHKTCD * KCKTAB * MDIMD
c             KCSINT  = KCPINT  + NPRFB * NCTFA * NCTFCD * MAXNPQBCX
c             KFCSINT = KCSINT  + NCCCC * KHKTCD * KCKTB
c             KHCPRIM = KFCSINT + KHKTCD * KCKTB 
c             KFRSTUV = KHCPRIM + MDIMD * NCCPP * KHKTCD
c             KEUV    = KFRSTUV + KHKTCD
c             KETUV   = KEUV    + NPRFAB
c             KLAST3  = KETUV   + NPRFAB
c             IF (KLAST3.GT.LWORK) CALL 
c    &             STOPIT('GETINT31','G_ABCD',KLAST3,LWORK)     
c             MAXLAST = MAX(MAXLAST,KLAST3)
c             
c             CALL DZERO(WORK(KAOINT),NCCT*KHKTAB*KHKTCD)

c             CALL G_ABCD(WORK(KAOINT), WORK(KFAOINT), WORK(KHCINT),
c    &             INDHER, IODDHR, IODDKC, WORK(KHCPRIM), WORK(KFRSTUV), 
c    &             WORK(KHCCONT), WORK(KCSINT), WORK(KFCSINT), 
c    &             WORK(KETIJAB), WORK(KEUV), WORK(KETUV), LMNPWR, 
c    &             IPNTUV, CSQ(KCSQA,1), CSQ(KCSQB,1), 
c    &             WORK(KCPINT), MDIMD, KC2MAX, 
c    &             NRDER, NCCCC,  NCCT, IABXYZ, KHKTAB, KHKTCD, KCKTAB, 
c    &             NHKTA, NHKTB, KCKTA, KCKTB, WORK(KCSVA), WORK(KCSVB),
c    &             CCFBTCH(KCCF(1,1)), CCFBTCH(KCCF(1,2)), IPRINT)
c             
c             CALL F_ABCD(RMAT,WORK(KAOINT),NCCT, 
c    &             NCTFA,NCTFB,NCTFC,NCTFD,
c    &             NHKTA,NHKTB,NHKTC,NHKTD, 
c    &             KHKTA,KHKTB,KHKTC,KHKTD,
c    &             IAOA,IAOB,IAOC,IAOD,
c    &             KHKTAB,KHKTCD,IAOAUXSTART,IPRINT)
             ENDIF ! ACTIVE CD-Batch
            END DO !ICDBTCH
          END DO !ICDCLS
         ENDIF ! ACTIVE AB-Batch
        END DO !IABBTCH
      END DO !IABCLS
      
C      write(LUPRI,*) 'Memory requirement in GETINT31 ', MAXLAST
      CALL QEXIT('GETINT31')
      RETURN 
      END

      SUBROUTINE SETEXPPQ(IPRCD,NPRFC,NPRFD,NPRFAB,NPRFCD,
     &                    EXPC,EXPD,EXPP,EXPQ,DATAP)
#include "implicit.h"
      DIMENSION IPRCD(2,NPRFCD),EXPC(NPRFC),EXPD(NPRFD),EXPP(NPRFAB),
     &          EXPQ(NPRFCD),DATAP(6,NPRFAB)

      CALL SETIPR(IPRCD, NPRFC, NPRFD)
      IJ = 0
      DO I=1,NPRFAB
        DO J=1,NPRFCD
          IJ = IJ + 1
          EXPP(IJ) = DATAP(1,I)
          IC = IPRCD(1,J)
          ID = IPRCD(2,J)
          C = EXPC(IC)
          D = EXPD(ID)
          Q = C + D
          EXPQ(IJ) = Q
        ENDDO
      ENDDO
      RETURN 
      END

C
      SUBROUTINE F_ABCD(RMAT,AOINT,NCCT, 
     &                  NCTFA,NCTFB,NCTFC,NCTFD,
     &                  NHKTA,NHKTB,NHKTC,NHKTD, 
     &                  KHKTA,KHKTB,KHKTC,KHKTD,
     &                  IAOA,IAOB,IAOC,IAOD,
     &                  KHKTAB,KHKTCD,NAUXSTART,
     &                  IND_GL,IPRINT)
C
C     Purpose: Copy AO-integrals (AOINTS) to the return-matrix (RMAT)
C              RMAT can be (alpha|beta), (alpha|cd), (ab|beta), (ab|cd)
C
#include "implicit.h"
#include "priunit.h"
#include "densfit.h"
#include "fck3int.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION AOINT(NCTFB,NCTFA,NCTFD,NCTFC,KHKTAB,KHKTCD),
     &          RMAT(NDIM(1),NDIM(2),NDIM(3),NDIM(4)), IND_GL(NBASISAUX)
C
      IF (IPRINT.GT.10) THEN
         CALL HEADER('AOINT in F_ABCD',-1)
         WRITE (LUPRI,'(1X,A,3I10)') ' NCCT, KHKTAB,KHKTCD ',
     &                                 NCCT, KHKTAB,KHKTCD 
         WRITE (LUPRI,'(1X,A,4I6)') ' NCTFC,NCTFD,NCTFA,NCTFB ',
     &                         NCTFC,NCTFD,NCTFA,NCTFB
         CALL OUTPUT(AOINT,1,NCCT,1,KHKTAB*KHKTCD,NCCT,
     &               KHKTAB*KHKTCD,1,LUPRI)
         CALL FLSHFO(LUPRI)
      END IF 
C
C     code
C
      IASTART = IAOA - 1
      IBSTART = IAOB - 1
      ICSTART = IAOC - 1
      IDSTART = IAOD - 1

      IF (AUX1) THEN
        IASTART = 0
        IBSTART = IBSTART - NAUXSTART + 1
        IF (DOBOX .AND. .NOT. GAB) THEN
          IBSTART = IND_GL(IBSTART+1) - 1
          IF (IBSTART.LT.0) CALL QUIT('IBSTART error in F_ABCD')
        ENDIF
      ENDIF

      IF (AUX2) THEN
        ICSTART = 0
        IDSTART = IDSTART - NAUXSTART + 1
        IF (DOBOX .AND. .NOT. GAB) THEN
          IDSTART = IND_GL(IDSTART+1) - 1
          IF (IDSTART.LT.0) CALL QUIT('IDSTART error in F_ABCD')
        ENDIF
      ENDIF


c      IF (DIATOM) THEN
c        IASTART = INDEX_ACT(IAOA)
c        IBSTART = INDEX_ACT(IAOB)
c        ICSTART = INDEX_ACT(IAOC)
c        IDSTART = INDEX_ACT(IAOD)
c        IF (AUX1) THEN
c          IASTART = 0
c        ENDIF
c        IF (AUX2) THEN
c          ICSTART = 0
c        ENDIF
c      ENDIF

      IF (GAB) THEN
        DO NCA=1,NCTFA
          IAOFF = (NCA - 1) * KHKTA
          DO NCB=1,NCTFB
            IBOFF = (NCB - 1) * KHKTB 
            IHKTAB = 0
            DO KA=1,KHKTA
              IA = IASTART + KA + IAOFF
              IBDIA = IBSTART + KA + IAOFF
              MAXB = KHKTB
              IF (TKMPAB) MAXB = KA
              DO KB=1,MAXB
                IHKTAB = IHKTAB + 1
                IB = IBSTART + KB + IBOFF
                IADIA = IASTART + KB + IBOFF
                SQRTAOINT = SQRT( AOINT(NCB,NCA,NCB,NCA,
     &                                  IHKTAB,IHKTAB)   )
                RMAT(IA,IB,1,1) = SQRTAOINT
                IF (.NOT.AUX1) RMAT(IB,IA,1,1) = SQRTAOINT
                IF (TKMPAB) THEN
                  IF ((IBDIA.EQ.0).OR.(IADIA.EQ.0)) 
     &              CALL QUIT('Programming error in F_ABCD')
                  RMAT(IBDIA,IADIA,1,1) = SQRTAOINT
                  RMAT(IADIA,IBDIA,1,1) = SQRTAOINT
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (AUX1.AND.AUX2) THEN
        DO NCA=1,NCTFA
          IAOFF = (NCA - 1) * KHKTA
          DO NCB=1,NCTFB
            IBOFF = (NCB - 1) * KHKTB 
            DO NCC=1, NCTFC
              ICOFF = (NCC - 1) * KHKTC
              DO NCD=1,NCTFD
                IDOFF = (NCD - 1) * KHKTD
                IHKTAB = 0
                DO KA=1,KHKTA
                  IA = IASTART + KA + IAOFF
                  DO KB=1,KHKTB
                    IHKTAB = IHKTAB + 1
                    IHKTCD = 0
                    IB = IBSTART + KB + IBOFF
                    DO KC=1,KHKTC
                      IC = ICSTART + KC + ICOFF
                      DO KD=1,KHKTD
                        IHKTCD = IHKTCD +  1
                        ID = IDSTART + KD + IDOFF
                        RMAT(IA,IB,IC,ID) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        RMAT(IC,ID,IA,IB) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
C     Would have to introduce TRIANGL=.FALSE. when alpha and beta do not
C     constitute the same set.
C                        IF (TRIANGL)
C     &                       RMAT(IC,ID,IA,IB) = AOINT(NCB,NCA,NCD,NCC,
C     &                         IHKTAB,IHKTCD)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (AUX1) THEN
        CALL QUIT('Given combination not tested. Compare with AUX2')
        DO NCA=1,NCTFA
          IAOFF = (NCA - 1) * KHKTA
          DO NCB=1,NCTFB
            IBOFF = (NCB - 1) * KHKTB 
            DO NCC=1, NCTFC
              ICOFF = (NCC - 1) * KHKTC
              DO NCD=1,NCTFD
                IDOFF = (NCD - 1) * KHKTD
                IHKTAB = 0
                DO KA=1,KHKTA
                  IA = IASTART + KA + IAOFF
                  DO KB=1,KHKTB
                    IHKTAB = IHKTAB + 1
                    IHKTCD = 0
                    IB = IBSTART + KB + IBOFF
                    DO KC=1,KHKTC
                      IC = ICSTART + KC + ICOFF
                      DO KD=1,KHKTD
                        IHKTCD = IHKTCD +  1
                        ID = IDSTART + KD + IDOFF
                        RMAT(IA,IB,IC,ID) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        RMAT(IA,IB,ID,IC) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (AUX2) THEN
        DO NCA=1,NCTFA
          IAOFF = (NCA - 1) * KHKTA
          DO NCB=1,NCTFB
            IBOFF = (NCB - 1) * KHKTB 
            DO NCC=1, NCTFC
              ICOFF = (NCC - 1) * KHKTC
              DO NCD=1,NCTFD
                IDOFF = (NCD - 1) * KHKTD
                IHKTAB = 0
                DO KA=1,KHKTA
                  IA = IASTART + KA + IAOFF
                  IBDIA = IBSTART + KA + IAOFF
                  MAXB=KHKTB
                  IF (TKMPAB) MAXB = KA
                  DO KB=1,MAXB
                    IHKTAB = IHKTAB + 1
                    IHKTCD = 0
                    IB = IBSTART + KB + IBOFF
                    IADIA = IASTART + KB + IBOFF
                    DO KC=1,KHKTC
                      IC = ICSTART + KC + ICOFF
                      DO KD=1,KHKTD
                        IHKTCD = IHKTCD +  1
                        ID = IDSTART + KD + IDOFF
                        RMAT(IA,IB,IC,ID) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        RMAT(IB,IA,IC,ID) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        IF (TKMPAB) THEN
                          RMAT(IBDIA,IADIA,IC,ID) = 
     &                       AOINT(NCB,NCA,NCD,NCC,IHKTAB,IHKTCD)
                          RMAT(IADIA,IBDIA,IC,ID) = 
     &                       AOINT(NCB,NCA,NCD,NCC,IHKTAB,IHKTCD)
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        WRITE(LUPRI,*) 'F_ABCD not tested for four center integrals',
     &                 'or auxilliary basis and '
        WRITE(LUPRI,*) 'F_ABCD implemented for GAB with normal',
     &                 'or auxilliary basis and '
        WRITE(LUPRI,*) 'for (alpha|beta) integrals.'
        WRITE(LUPRI,*) 'GAB, AUX1, AUX2 ', GAB, AUX1, AUX2
        CALL QUIT('Given combination not tested.')
        DO NCA=1,NCTFA
          IAOFF = (NCA - 1) * KHKTA
          DO NCB=1,NCTFB
            IBOFF = (NCB - 1) * KHKTB 
            DO NCC=1, NCTFC
              ICOFF = (NCC - 1) * KHKTC
              DO NCD=1,NCTFD
                IDOFF = (NCD - 1) * KHKTD
                IHKTAB = 0
                DO KA=1,KHKTA
                  IA = IASTART + KA + IAOFF
                  DO KB=1,KHKTB
                    IHKTAB = IHKTAB + 1
                    IHKTCD = 0
                    IB = IBSTART + KB + IBOFF
                    DO KC=1,KHKTC
                      IC = ICSTART + KC + ICOFF
                      DO KD=1,KHKTD
                        IHKTCD = IHKTCD +  1
                        ID = IDSTART + KD + IDOFF
                        RMAT(IA,IB,IC,ID) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        RMAT(IB,IA,IC,ID) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        RMAT(IA,IB,ID,IC) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                        RMAT(IB,IA,ID,IC) = AOINT(NCB,NCA,NCD,NCC,
     &                       IHKTAB,IHKTCD)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
        
      IF (IPRINT.GT.10) THEN
         CALL HEADER('RMAT in F_ABCD',-1)
         IF (GAB.AND..NOT.AUX1) THEN
           NDIMENSION1 = NDIM(1)
           NDIMENSION2 = NDIM(2)
         ELSE
           NDIMENSION1 = NDIM(1) * NDIM(2)
           NDIMENSION2 = NDIM(3) * NDIM(4)
         ENDIF
         WRITE (LUPRI,'(1X,A,2I5)') 'DIMENSION1 , DIMENSION2  ',
     &                               NDIMENSION1 , NDIMENSION2
         CALL OUTPUT(RMAT,1,NDIMENSION1,1,NDIMENSION2,
     &               NDIMENSION1,NDIMENSION2,1,LUPRI)
      END IF 

      RETURN
      END 


      SUBROUTINE GCDTUV(HCINT, KHKTCD, KCKTCD, NHKTC,
     &                  NHKTD, KCKTC, KCKTD, LMNPWR, IODDKC, IODDHR, 
     &                  INDHER, HERINT, HCPRIM, ETIJCD, EUV, ETUV,
     &                  FRSTUV, HCCONT, CMATC, CMATD, CPINT, INDHSQ, 
     &                  NCCPP, NTUVAB, FCINT, HCSINT, FCSINT,
     &                  CSQC, CSQD, IPNTUV, GCNTCD, GCNTAB, NPRFCD, 
     &                  SPHR12, IPRINT)
C
C     Generates G^CD_tuv in spherical AO representation. First we generate 
C             
C        G^cd_tuv = sum_{tau,nu,phi}F^cd_{tau,nu,phi} R_{t+tau,u+nu,v+phi}
C     
C     in primitive representation. We then contract to cartesian 
C     representation, and finally to spherical representation.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "fck3int.h"
#include "fck3odb.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
      INTEGER TUV
      LOGICAL FRSTUV, CONTRA, SPHRCD, GCNTCD, GCNTAB
      DIMENSION HCINT(NCCPP*NTUVAB,KHKTCD), LMNPWR(MAXAOKCK,MAXAONHK,3), 
     &          IODDKC(MAXAOKCK**2), IODDHR(NRTOP,0:7), 
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), HERINT(*), HCPRIM(*), 
     &          ETIJCD(*), EUV(*),  ETUV(*), FRSTUV(*), 
     &          HCCONT(NCCPP*NTUVAB, KCKTCD), 
     &          CMATC(*), CMATD(*), CPINT(*), INDHSQ(*), FCINT(*), 
     &          HCSINT(*), FCSINT(*), CSQC(*), CSQD(*), IPNTUV(*)
     &          
      CONTRA = (GCNTCD) .OR. (NPRFCD.NE.1)

      IF (KHKTCD .EQ. 1) THEN
        DO TUV = 1, NTUVAB
          FRSTUV(TUV) = IODDHR(TUV,IPQXYZFCK3) .NE. 0
        ENDDO
        CALL ERICT1(HERINT,HCINT,
     &              CMATC,CMATD,CPINT,IODDHR(1,IHCBITFCK3),
     &              0,FRSTUV,IPRINT)
      ELSE
        ICMPCD = 0
        DO ICOMPC = 1, KCKTC
          LC = LMNPWR(ICOMPC, NHKTC, 1)
          MC = LMNPWR(ICOMPC, NHKTC, 2)
          NC = LMNPWR(ICOMPC, NHKTC, 3)
          MAXD = KCKTD
          IF (TCMPCD) MAXD = ICOMPC
          DO ICOMPD = 1, MAXD
            ICMPCD = ICMPCD + 1
            LD = LMNPWR(ICOMPD, NHKTD, 1)
            MD = LMNPWR(ICOMPD, NHKTD, 2)
            ND = LMNPWR(ICOMPD, NHKTD, 3)
            IODDKC(ICMPCD) = IODDHR(INDHER(LC+LD,MC+MD,NC+ND),
     &                              IHCBITFCK3)
C           contracted spherical harmonics
C           ==============================
            IF (CONTRA .AND. SPHERCD) THEN
C             Returns GcdTUV in HCPRIM (primitive)
              CALL CR1UND(HERINT,DUMMY,HCPRIM,
     &             LC,MC,NC,LD,MD,ND,INDHER,INDHER,
     &             IODDHR(1,IPQXYZFCK3),INDHSQ,
     &             ETIJCD,EUV,ETUV,FRSTUV,ICOMPC,ICOMPC,
     &             DUMMY,IPRINT)
C             Returns contracted GcdTUV in HCCONT 
C             (contracted Cartesian)
              CALL ERICT1(HCPRIM,HCCONT(1,ICMPCD),
     &             CMATC,CMATD,CPINT,IODDHR(1,IHCBITFCK3),
     &             IODDKC(ICMPCD),FRSTUV,IPRINT)
C           contracted Cartesians
C           =====================
            ELSEIF (CONTRA) THEN
              CALL CR1UND(HERINT,DUMMY,HCPRIM,
     &             LC,MC,NC,LD,MD,ND,INDHER,INDHER,
     &             IODDHR(1,IPQXYZFCK3),INDHSQ,
     &             ETIJCD,EUV,ETUV,FRSTUV,ICOMPC,ICOMPC,
     &             DUMMY,IPRINT)
              CALL ERICT1(HCPRIM,HCINT(1,ICMPCD),
     &             CMATC,CMATD,CPINT,IODDHR(1,IHCBITFCK3),
     &             IODDKC(ICMPCD),FRSTUV,IPRINT)
C           primitive spherical harmonics
C           =============================
            ELSEIF (SPHERCD) THEN
              CALL CR1UND(HERINT,DUMMY,HCCONT(1,ICMPCD),
     &             LC,MC,NC,LD,MD,ND,INDHER,INDHER,
     &             IODDHR(1,IPQXYZFCK3),INDHSQ,
     &             ETIJCD,EUV,ETUV,FRSTUV,ICOMPC,ICOMPC,
     &             DUMMY,IPRINT)
              CALL ERICT1(HCCONT(1,ICMPCD),HCCONT(1,ICMPCD),
     &             CMATC,CMATD,CPINT,IODDHR(1,IHCBITFCK3),
     &             IODDKC(ICMPCD),FRSTUV,IPRINT)
C           primitive Cartesians
C           ====================
            ELSE
              CALL CR1UND(HERINT,DUMMY,HCINT(1,ICMPCD),
     &             LC,MC,NC,LD,MD,ND,INDHER,INDHER,
     &             IODDHR(1,IPQXYZFCK3),INDHSQ,
     &             ETIJCD,EUV,ETUV,FRSTUV,ICOMPC,ICOMPC,
     &             DUMMY,IPRINT)
              CALL ERICT1(HCINT(1,ICMPCD),HCINT(1,ICMPCD),
     &             CMATC,CMATD,CPINT,IODDHR(1,IHCBITFCK3),
     &             IODDKC(ICMPCD),FRSTUV,IPRINT)
            ENDIF
         ENDDO
        ENDDO

C       Returns GCDTUV in HCINT (contracted spherical)
        IF (SPHERCD) THEN
          CALL CR1SPH(HCCONT, HCINT, FCINT, HCSINT, FCSINT,
     &         CSQC, CSQD, IODDHR(1,IHCBITFCK3),
     &         IPNTUV, IODDKC, IPRINT)
        ENDIF  
      ENDIF
      
      IF (IPRINT .GT. 15) THEN
        CALL HEADER('HCINT in GCDTUV',-1)
        CALL HCPRINT(HCINT,IODDHR(1,IHCBITFCK3),IPNTUV,1)
      ENDIF
      RETURN
      END 

C     Generates g_ABCD in contracted spherical representation
      SUBROUTINE G_ABCD(AOINT, FAOINT, HCINT, INDHER, IODDHR, 
     &     IODDKC, CCPRIM, FRSTUV, CCONT,
     &     CSINT, FCSINT, ECOEF, EUV, ETUV, LMNPWR, IPNTUV, 
     &     CSQA, CSQB, CPINT, NDIMD, KC2MAX, NDER, NCCCC, 
     &     NCCT, IABXYZ, KHKTAB, KHKTCD, KCKTAB, NHKTA, NHKTB,
     &     KCKTA, KCKTB, CSVA, CSVB, CMATA, CMATB,
     &     IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "fck3int.h"
#include "fck3odb.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
      LOGICAL FRSTUV, FLEVEL, BDERA, BDERB
      DIMENSION AOINT(NCCT, KHKTAB, KHKTCD, NDIMD ),FAOINT(*),
     &          HCINT(*), INDHER(0:JTOP,0:JTOP,0:JTOP), 
     &          IODDHR(NRTOP,0:7), IODDKC(MAXAOKCK**2), CCPRIM(*), 
     &          FRSTUV(*), CCONT(NCCCC, KHKTCD, KCKTAB), CSINT(*),
     &          FCSINT(*), ECOEF(*), EUV(*), ETUV(*), 
     &          LMNPWR(MAXAOKCK,MAXAONHK,3), IPNTUV(KC2MAX, 0:NDER, 2), 
     &          CSQA(*), CSQB(*), CPINT(*),  CSVA(*), CSVB(*), 
     &          CMATA(*), CMATB(*)
cSimen CAREFUL, only when no derivatives:
      MAXDER = 0
      BDERA  = .FALSE.
      BDERB  = .FALSE.

      IELCT2 = 2

      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Output from G_ABCD',-1)
         WRITE (LUPRI,'(2X,A,3I5)') 'NCCCC, NHKTA, NHKTB ', 
     &                               NCCCC, NHKTA, NHKTB
         WRITE (LUPRI,'(2X,A,2I5)') 'KCKTA,NHKTB         ',KCKTA,NHKTB
         WRITE (LUPRI,'(2X,A,2I5)') 'MAXDER,KHKTAB       ',MAXDER,KHKTAB
         WRITE (LUPRI,'(2X,A,2I5)') 'KHKTCD,IELCT2       ',KHKTCD,IELCT2
         WRITE (LUPRI,'(2X,A,2I5)') 'KC2MAX, NDER        ',KC2MAX,NDER
      ENDIF

      IF ( (MAXDER.EQ.0) .AND. (KHKTAB .EQ. 1) ) THEN
        DO ICMPCD = 1, KHKTCD
          FRSTUV(ICMPCD) = 
     &         IODDHR(IPNTUV(ICMPCD, 0, IELCT2),IHCBITFCK3) .NE. 0
        CALL ERICT2(HCINT, AOINT, CMATA, CMATB, 
     &              CPINT, NCCT, 1, IODDHR(1,ICCXYZFCK3), IPNTUV, 
     &              0, FRSTUV, IPRINT)
        ENDDO
      ELSE
        FLEVEL = .TRUE.
        ICMPAB = 0
        DO ICOMPA = 1, KCKTA
          LA = LMNPWR(ICOMPA, NHKTA, 1)
          MA = LMNPWR(ICOMPA, NHKTA, 2)
          NA = LMNPWR(ICOMPA, NHKTA, 3)
          MAXB = KCKTB
          IF (TCMPAB) MAXB = ICOMPA
          DO ICOMPB = 1, MAXB
            ICMPAB = ICMPAB + 1
            LB = LMNPWR(ICOMPB, NHKTB, 1)
            MB = LMNPWR(ICOMPB, NHKTB, 2)
            NB = LMNPWR(ICOMPB, NHKTB, 3)

            IODDKC(ICMPAB) = 
     &           IODDHR(INDHER(LA+LB,MA+MB,NA+NB),IHCBITFCK3)

            CALL CR2UND(LA,MA,NA,LB,MB,NB,
     &                  HCINT,INDHER,IODDHR(1,IHCBITFCK3),CCPRIM,
     &                  FRSTUV,ECOEF,EUV,ETUV,IPNTUV,
     &                  ICOMPA,ICOMPB,IPRINT)

            IF (SPHERAB) THEN
              CALL ERICT2(CCPRIM,CCONT(1,1,ICMPAB),CMATA,
     &             CMATB,CPINT,NCCCC,1, IODDHR(1,ICCXYZFCK3),
     &             IPNTUV,IODDKC(ICMPAB),FRSTUV,IPRINT)
            ELSE
              CALL ERICT2(CCPRIM,AOINT(1,ICMPAB,1,1),CMATA,
     &             CMATB,CPINT,NCCT,KHKTAB,IODDHR(1,ICCXYZFCK3),
     &             IPNTUV,IODDKC(ICMPAB),FRSTUV,IPRINT)              
            ENDIF
          ENDDO
        ENDDO

        IODX = IODDHR(INDHER(0,0,0), ICCXYZFCK3)
        IF (SPHERAB) THEN
          CALL CR2SPH(CCONT,CSINT,FCSINT,
     &         AOINT,FAOINT,FLEVEL, CSQA, CSQB,
     &         CSAA,CSAB,CSVA,CSVB,BDERA,BDERB,NCENTR,
     &         IODDHR(1,ICCXYZFCK3),IPNTUV,IODDKC,0,IODX,IPRINT)
        ENDIF
        FLEVEL = .FALSE.
      ENDIF
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - G_ABCD',-1)
         DO 900 K = 1, NDIMD
            IF (NDER.GT.0) THEN
               INDX = MOD(K-1,3)+1
            ELSE
               INDX = 0 
            END IF 
            DO 910 I = 1, KHKTAB
            DO 910 J = 1, KHKTCD
cSimen CAREFUL not working properly
            IF (IODDHR(IPNTUV(I,0,1),ICCXYZFCK3) .EQ. 
     &          IODDHR(IPNTUV(J,INDX,2),ICCXYZFCK3)) THEN
               WRITE (LUPRI, '(/1X,A,3I3)') ' K, ICMPAB, ICMPCD ', K,I,J
               CALL OUTPUT(AOINT(1,I,J,K),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
            END IF
  910       CONTINUE
  900    CONTINUE
      END IF

      RETURN
      END 
C
C   --------------------------------------------------------------
C
      SUBROUTINE SETIPR(IPR12, NPRF1, NPRF2)
#include "implicit.h"
      DIMENSION IPR12(2,NPRF1*NPRF2)
      IPRF12 = 0
      DO IPRF1=1, NPRF1
        DO IPRF2=1, NPRF2
          IPRF12 = IPRF12 + 1
          IPR12(1,IPRF12) = IPRF1
          IPR12(2,IPRF12) = IPRF2
        ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE SETQR( FACINT, ALPHA, COORPQ, DATAP, IPRAB, IPRCD,
     &                  COORC, COORD, EXPBTC, EXPBTD, NPRFC, NPRFD,
     &                  NPRFCD, NPRF, IPQXYZ, NPRFAB, JMAX, NTYP, NTUV, 
     &                  CMATA, CMATB, CMATC, CMATD, GCNTAB, GCNTCD)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "cbifck3.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "fck3int.h"
#include "densfit.h"
      DIMENSION DATAP(6,*),COORPQ(NPRF,3),FACINT(*),ALPHA(*),IPQ0(3), 
     &          IPRAB(2,NPRFAB), IPRCD(2,NPRFCD), COORC(3), COORD(3), 
     &          EXPBTC(NPRFC), EXPBTD(NPRFD), 
     &          CMATA(*), CMATB(*), CMATC(*), CMATD(*)
      LOGICAL GCNTAB, GCNTCD

      I = 0 

      DO IPRFAB=1,NPRFAB

       EXPP  = DATAP(1,IPRFAB)
       EXPPI = DATAP(2,IPRFAB)
       FACPI = DATAP(3,IPRFAB)
       PX    = DATAP(4,IPRFAB)
       PY    = DATAP(5,IPRFAB)
       PZ    = DATAP(6,IPRFAB)


       IPRFA = IPRAB(1,IPRFAB)
       IPRFB = IPRAB(2,IPRFAB)

       DO IPRFCD=1,NPRFCD

         I = I + 1

         IPRFC = IPRCD(1,IPRFCD)
         IPRFD = IPRCD(2,IPRFCD)
         EXPC  = EXPBTC(IPRFC)
         EXPD  = EXPBTD(IPRFD)
           
         EXPQ  = EXPC + EXPD
         EXPQI = 1.0D0 / EXPQ
         QX = ( COORC(1)*EXPC + COORD(1)*EXPD ) * EXPQI
         QY = ( COORC(2)*EXPC + COORD(2)*EXPD ) * EXPQI
         QZ = ( COORC(3)*EXPC + COORD(3)*EXPD ) * EXPQI
         
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     
C     The definition of the alpha in the equations
C     is alpha_def = p*q/(p+q)
C     
C     However we have to put alpha = -2*alpha_def in order
C     to be compatible with way how with the erigam uses 
C     the alpha  
C     
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ALPHA(I)  = - 2.0D0 / (EXPPI + EXPQI) 
         FAC = FACPI*EXPQI/SQRT(EXPP+EXPQ)
         IF (.NOT.GCNTAB) FAC = FAC * CMATA(IPRFA) * CMATB(IPRFB)
         IF (.NOT.GCNTCD) FAC = FAC * CMATC(IPRFC) * CMATD(IPRFD)
         FACINT(I) = FAC

         COORPQ(I,1)= PX - QX
         COORPQ(I,2)= PY - QY
         COORPQ(I,3)= PZ - QZ
       END DO
      END DO

      CALL LOCSYM( COORPQ, NPRF, 1, IPQ0, IPQXYZ )

      CALL PUT2ERICOM2( NPRF, IPQ0, IPQXYZ, JMAX, NTUV )

      RETURN
      END


      SUBROUTINE FACEXP(FACINT, PREFAB, PREFCD, NPRFAB, NPRFCD)
#include "implicit.h"
      DIMENSION FACINT(*), PREFAB(*), PREFCD(*)
      I = 0 
      DO IPRFAB=1,NPRFAB
        DO IPRFCD=1,NPRFCD
          I = I + 1
          FACINT(I) = FACINT(I) * PREFAB(IPRFAB) * PREFCD(IPRFCD)
        ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE GETDIST( IELECTR, ICENT1, ICENT2, 
     &                    NSYMOP1, NSYMOP2, NACTCENT, I0, IXYZ, AUX )
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "fck3int.h"
      LOGICAL AUX
      DIMENSION I0(3)

      IF (IELECTR .EQ. 1) THEN
        IF (ICENT1.EQ.0) THEN
          COOR(1,0,1) = 0.0d0
          COOR(2,0,1) = 0.0d0
          COOR(3,0,1) = 0.0d0
        ELSE
          COOR(1,0,1) = CORD(1,ICENT1)
          COOR(2,0,1) = CORD(2,ICENT1)
          COOR(3,0,1) = CORD(3,ICENT1)
        ENDIF

        IF (ICENT2.EQ.0) THEN
          COOR(1,0,2) = 0.0d0
          COOR(2,0,2) = 0.0d0
          COOR(3,0,2) = 0.0d0
        ELSE
          COOR(1,0,2) = CORD(1,ICENT2)
          COOR(2,0,2) = CORD(2,ICENT2)
          COOR(3,0,2) = CORD(3,ICENT2)
        ENDIF
        DIST12(1,0,0,1) = COOR(1,0,1) - COOR(1,0,2)
        DIST12(2,0,0,1) = COOR(2,0,1) - COOR(2,0,2)
        DIST12(3,0,0,1) = COOR(3,0,1) - COOR(3,0,2)

        IF (AUX) THEN
          DO I=1,3
            I0(I) = 0
          ENDDO
          IXYZ = 0
        ELSE
          CALL LOCSYM(DIST12(1,0,0,1),1,1,I0,IXYZ)
        ENDIF

        DIST122(0,0,1)  = DIST12(1,0,0,1)**2 + DIST12(2,0,0,1)**2
     &                    + DIST12(3,0,0,1)**2 
      ELSE IF (IELECTR .EQ. 2) THEN 
        IF (ICENT1.EQ.0) THEN
          COOR(1,0,3) = 0.0d0
          COOR(2,0,3) = 0.0d0
          COOR(3,0,3) = 0.0d0
        ELSE
          COOR(1,0,3) = CORD(1,ICENT1)
          COOR(2,0,3) = CORD(2,ICENT1)
          COOR(3,0,3) = CORD(3,ICENT1)
        ENDIF
        IF (ICENT2.EQ.0) THEN
          COOR(1,0,4) = 0.0d0
          COOR(2,0,4) = 0.0d0
          COOR(3,0,4) = 0.0d0
        ELSE
          COOR(1,0,4) = CORD(1,ICENT2)
          COOR(2,0,4) = CORD(2,ICENT2)
          COOR(3,0,4) = CORD(3,ICENT2)
        ENDIF
        DIST12(1,0,0,2) = COOR(1,0,3) - COOR(1,0,4)
        DIST12(2,0,0,2) = COOR(2,0,3) - COOR(2,0,4)
        DIST12(3,0,0,2) = COOR(3,0,3) - COOR(3,0,4)

        IF (AUX) THEN
          DO I=1,3
            I0(I) = 0
          ENDDO
          IXYZ = 0
        ELSE
          CALL LOCSYM(DIST12(1,0,0,1),1,1,I0,IXYZ)
        ENDIF

        DIST122(0,0,2)  = DIST12(1,0,0,2)**2 + DIST12(2,0,0,2)**2
     &                    + DIST12(3,0,0,2)**2 
      ELSE
        CALL QUIT('Wrong value of IELECTR in GETDIST')
      ENDIF

      RETURN
      END 

      SUBROUTINE SETINT(AUXEL1, AUXEL2, DOGAB, NDIMAUX)
#include "implicit.h"
#include "priunit.h"
#include "fck3int.h"
#include "densfit.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "nuclei.h"
#include "shells.h"
#include "fck3aob.h"
      LOGICAL AUXEL1, AUXEL2, DOGAB
      INTEGER N_NORM(4), N_FIT(4), I_NORM(4), I_FIT(4) 
#include "ibtfun.h"
      AUX1 = AUXEL1
      AUX2 = AUXEL2
      GAB  = DOGAB

      ODTR12 = AUX1 .EQV. AUX2 .AND.(.NOT.GAB)

      NBASE  = 0
      DO I = 1, KMAX
         IF (IBTAND(K,ISTBAO(I)).EQ.0) NBASE = NBASE + KHKT(I)
      END DO
      
      DO J=1,4
        N_NORM(J) = NBASE
        N_FIT(J)  = NBASEAUX
        IF (DOBOX) N_FIT(J) = NDIMAUX
      ENDDO
      IF (DIATOM) THEN
        CALL QUIT('DIATOM not working.')
        DO J=1,4
          N_NORM(J) = 0
          N_FIT(J)  = 0
          I_NORM(J)   = 0
          DO ISH = 1, NAOBTCH
            IF (ACTVBTH(ISH,J)) THEN
              ICLASS = KCLSBTCH(ISH)
              IF (ICLASS .EQ. 1) THEN
                INDEX_ACT(KAOBTCH(1,ISH)) = I_NORM(J)
                NHK    = NHKTBTCH(ISH)
                KHK    = NHK * 2 - 1
                N_NORM(J) = N_NORM(J) + KHK
                I_NORM(J) = I_NORM(J) + KHK * NCTFBTCH(ISH)
              ELSEIF (ICLASS .EQ. 3) THEN
                INDEX_ACT(KAOBTCH(1,ISH)) = N_FIT(J)
                NHK    = NHKTBTCH(ISH)
                KHK    = NHK * 2 - 1
                N_FIT(J) = N_FIT(J) + KHK * NCTFBTCH(ISH)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      DO i=1,4
        NDIM(i) = N_NORM(i)
      END DO
      IF (AUX1) THEN
        NDIM(1) = 1
        NDIM(2) = N_FIT(1)
      END IF
      IF (AUX2) THEN
        NDIM(3) = 1        
        NDIM(4) = N_FIT(3)
      END IF

      IF (GAB) THEN
        NDIM(3) = 1
        NDIM(4) = 1
        IF (AUX1.NEQV.AUX2) 
     &       CALL QUIT('SETINT: Can use GAB=.TRUE. only when AUX1=AUX2')
      END IF

      RETURN
      END 


      SUBROUTINE INILOOP
#include "implicit.h"
#include "priunit.h"
#include "fck3int.h"
#include "densfit.h"

      IF (AUX1) THEN
        IABCLSSTART = NODCLSNORM + 1
        IABCLSEND   = NODCLSNORM + NODCLSAUX
      ELSE
        IABCLSSTART = 1
        IABCLSEND   = NODCLSNORM
      ENDIF

      IF (AUX2) THEN
        ICDCLSSTART = NODCLSNORM + 1
        ICDCLSEND   = NODCLSNORM + NODCLSAUX
      ELSE
        ICDCLSSTART = 1
        ICDCLSEND   = NODCLSNORM
      ENDIF

      RETURN
      END 

      
      SUBROUTINE PUT2ERI(NPRF3VAL, NPRF4VAL, NPRF1VAL, NPRF2VAL, 
     &                   NCTF3VAL, NCTF4VAL, NCTF1VAL, NCTF2VAL, 
     &                   KCKT3VAL, KCKT4VAL, KCKT1VAL, KCKT2VAL, 
     &                   KHKT3VAL, KHKT4VAL, KHKT1VAL, KHKT2VAL, 
     &                   NHKT3VAL, NHKT4VAL, NHKT1VAL, NHKT2VAL, 
     &                   NPQBCXVAL, NPP34VAL, NPP12VAL, 
     &                   NCTF34VAL, NCTF12VAL,
     &                   NPRF34VAL, NPRF12VAL, JMAX3VAL, JMAX4VAL, 
     &                   JMAX1VAL, JMAX2VAL, NTUV34VAL, 
     &                   KCKT34VAL, KCKT12VAL, KHKT34VAL, KHKT12VAL, 
     &                   KC2MAXVAL, NHKMAXVAL, NRDERVAL, 
     &                   SPHR3VAL, SPHR4VAL, SPHR1VAL, SPHR2VAL, 
     &                   TKMP34VAL, TKMP12VAL, 
     &                   GCON3VAL, GCON4VAL, GCON1VAL, GCON2VAL, 
     &                   I120VAL, I340VAL, I12XYZVAL, I34XYZVAL )
#include "implicit.h"
#include "priunit.h"
#include "ericom.h"
#include "eriao.h"
      DIMENSION I120VAL(3), I340VAL(3)
      LOGICAL SPHR1VAL, SPHR2VAL, SPHR3VAL, SPHR4VAL, 
     &        TKMP12VAL, TKMP34VAL, 
     &        GCON1VAL, GCON2VAL, GCON3VAL, GCON4VAL


      NPRF1   = NPRF1VAL
      NPRF2   = NPRF2VAL
      NPRF3   = NPRF3VAL
      NPRF4   = NPRF4VAL

      NCTF1   = NCTF1VAL
      NCTF2   = NCTF2VAL
      NCTF3   = NCTF3VAL
      NCTF4   = NCTF4VAL

      KCKT1   = KCKT1VAL
      KCKT2   = KCKT2VAL
      KCKT3   = KCKT3VAL
      KCKT4   = KCKT4VAL

      KCKTC   = KCKT1
      KCKTD   = KCKT2
      KCKTA   = KCKT3
      KCKTB   = KCKT4

      KHKT1   = KHKT1VAL
      KHKT2   = KHKT2VAL
      KHKT3   = KHKT3VAL
      KHKT4   = KHKT4VAL

      KHKTC   = KHKT1
      KHKTD   = KHKT2
      KHKTA   = KHKT3
      KHKTB   = KHKT4

      NHKT1   = NHKT1VAL
      NHKT2   = NHKT2VAL
      NHKT3   = NHKT3VAL
      NHKT4   = NHKT4VAL

      NHKTC   = NHKT1
      NHKTD   = NHKT2
      NHKTA   = NHKT3
      NHKTB   = NHKT4

      NPQBCX  = NPQBCXVAL
      NPP12   = NPP12VAL
      NPP34   = NPP34VAL

      NCTF12  = NCTF12VAL
      NCTF34  = NCTF34VAL

      NPRF12  = NPRF12VAL
      NPRF34  = NPRF34VAL
      JMAX1   = JMAX1VAL
      JMAX2   = JMAX2VAL
      JMAX3   = JMAX3VAL
      JMAX4   = JMAX4VAL

CAREFUL Because we have no derivatives BDER, GDER both = .FALSE.
      NDER3   = 0
      NDER4   = 0
      JMXD3   = JMAX3 + NDER3
      JMXD4   = JMAX4 + NDER4

      NTUV34  = NTUV34VAL
      KCKT12  = KCKT12VAL
      KHKT12  = KHKT12VAL

CAREFUL-filip believes it is fine
      KHKT34  = KHKT34VAL
      KCKT34  = KCKT34VAL

      KHKTAB  = KHKT34
      KHKTCD  = KHKT12
CAREFUL

      GCON1   = GCON1VAL
      GCON2   = GCON2VAL
      GCON3   = GCON3VAL
      GCON4   = GCON4VAL
      GCON12  = GCON1 .OR. GCON2
      GCON34  = GCON3 .OR. GCON4

      DO I=1,3
        I120(I) = I120VAL(I)
        I340(I) = I340VAL(I)
      ENDDO
      IABXYZ  = I12XYZVAL
      ICDXYZ  = I34XYZVAL
      KC2MAX  = KC2MAXVAL
      NHKMAX  = NHKMAXVAL
      NRDER   = NRDERVAL

      SPHR1   = SPHR1VAL
      SPHR2   = SPHR2VAL
      SPHR3   = SPHR3VAL
      SPHR4   = SPHR4VAL
      SPHR12  = SPHR1 .OR. SPHR2
      SPHR34  = SPHR3 .OR. SPHR4
      TKMP12  = TKMP12VAL
      TKMP34  = TKMP34VAL
      TCMP12  = TKMP12 .AND. .NOT. SPHR12
      TCMP34  = TKMP34 .AND. .NOT. SPHR34

C     CAREFUL, because we have IPATH = 2
      IPATH   = 2
      IELCT1  = 2
      IELCT2  = 1
      TKMPCD  = TKMP12
      TKMPAB  = TKMP34
      TCMPCD  = TCMP12
      TCMPAB  = TCMP34

C     Because of no symmetry
      MLTPZVAL= 1 

      NCCPP   = NPQBCXVAL * NCTF12 * NPRF34 * MLTPZVAL
      NCCCC   = NPQBCXVAL * NCTF12 * NCTF34 * MLTPZVAL
      NCCT    = MLTPZVAL * NPQBCXVAL * NCTF12 * NCTF34

      RETURN
      END 


C
C -----------------------------------------------------
C

      SUBROUTINE PUTAO2ERICOM(MAXAONHK,MAXAOKCK)
#include "implicit.h"
#include "ericom.h"
      NHKMAX = MAXAONHK
      KCKMAX = MAXAOKCK
      RETURN
      END


C     Calculates the charges of the auxilliary functions, and returns it in OAUX.
C     This means only the one-center integrals; contract with the fitted 
C     coefficients to get the charge. Only s-type integrals are non-zero in the
C     spherical-harmonic basis.
      SUBROUTINE AUXCHA(OAUX,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
      DIMENSION WORK(*), OAUX(*)
      CALL QENTER('AUXCHA')

C     Dimensions for CSQ
      CALL SPHDIM

C     Allocate
      KCSQ  = 1
      KLAST = KCSQ + NCSQ1 * NCSQ2
      IF (KLAST .GT. LWORK) CALL STOPIT('AUXCHA','KCSQ  ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1

      CALL AUXCH1( OAUX,WORK(KCSQ),WORK(KLAST),LWRK,IPRINT)
      CALL QEXIT('AUXCHA')
      RETURN
      END

      SUBROUTINE AUXCH1(OAUX,CSQ,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "nuclei.h"
#include "pi.h"
      DIMENSION CSQ(NCSQ1,NCSQ2), OAUX(*),WORK(*),KCCF(2)
      LOGICAL FIRST
      PARAMETER (PI32 = SQRTPI*PI)
      PR1CINT(A) = PI32/SQRT(A**3)
      CALL QENTER('AUXCH1')

      CALL DZERO(OAUX,NBASISAUX)
      FIRST = .TRUE.
      IAOIND = 0
      DO IAOBTCH = 1,NAOBTCH - 1
        ICLASS = KCLSBTCH(IAOBTCH)
        IAO    = KAOBTCH(1,IAOBTCH)
        IF (ICLASS .EQ. 3) THEN
          IF (FIRST) THEN 
C           The fitting basis starts after the normal basis (and Huckel);
C           add offset here (to start fitting basis index at 0).
            IOFF = KAOBTCH(1,IAOBTCH)
            FIRST = .FALSE.
          ENDIF
          CALL GETAOB(NPRF,NCTF,NTYP,NHKT,JMAX,KHKT,KCKT,
     &                KCSQ,KEXP,KCCF,IAOBTCH)
C         All total angular momentum greater than zero is zero by symmetry
          IF (NHKT.EQ.1) THEN
            DO ICTF=1,NCTF
              DO IPRF=1,NPRF
C               extract PGTO exponent
                A  = EXPBTCH(KEXP + IPRF -1)
C               extract PGTO contraction coefficient
                CC = CCFBTCH(KCCF(1) + (ICTF-1)*NCTF + IPRF - 1)
                OAUX(IAO-IOFF+ICTF)=OAUX(IAO-IOFF+ICTF) + CC*PR1CINT(A)
              ENDDO
            ENDDO
          ENDIF
        ENDIF    !LCLASS
      ENDDO    !IAOBTCH
      CALL QEXIT('AUXCH1')
      RETURN 
      END

