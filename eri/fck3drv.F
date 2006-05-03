C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C
C
C$Id$
C
C   FCK3: Program for fast evaluation of coulomb contribution to the Fock matrix
C

      SUBROUTINE FCK3DRV( FMAT, DMAT, NDMT, WORK, LWORK )
#include "implicit.h"
#include "priunit.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbihr2.h"
#include "cbieri.h"
Cdirac#include "cbione.h"
#include "cbifck3.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#include "fck3mem.h"
#include "erithr.h"
#include "clsfmm.h"

      DIMENSION TIMFCK3(5),FMAT(*), DMAT(*), WORK(*)
C
      WORD2MB = 8.0D-6
      STATISTICS(1) = 0.0D0
      STATISTICS(2) = 0.0D0
      STATISTICS(3) = 0.0D0
      STATISTICS(4) = 0.0D0
      STATISTICS(5) = 0.0D0

      IF ( FIXTHR ) THEN
        THRSH = THRFCK3_DEFAULT
      ELSE
CTROND...IFTHRS is not set, but I set FIXTHR to true
        THRSH = 10D0**(-IFTHRS)
      END IF

      THRPOD   = THRSH
      THRCOD   = THRSH
      THRNCL   = THRSH 
      IF ( NOSCRPRIM ) THRPOD = 0.0D0

#ifdef MOD_FMM
      IF (DOFMM) THEN
        NCLFCK3 = .TRUE.
        THRCLS  = THRNCL
      END IF

      CALL SETBFCFMM( BRFREE ) 
      FACCLS=ERFCIV(THRCLS)
#endif

      IF ( NDMT .NE. 1 ) THEN 
        CALL QUIT('In FCK3 internally is not implemented NDMT > 1') 
      END IF

      NDMAT = NDMT

      KFREE = 1 
      LFREE = LWORK

      CALL TIMER('START ',TIMSTR,TIMEND)

      CALL GETTIM(TIMFCK3(1),DUMTIM)
C
C     Set AO-batches
C
      CALL SETAOBTCH(SPEXP,CONTRACTED,WORK,LWORK,IPRFCK3)

      CALL GETTIM(TIMFCK3(2),DUMTIM)
C
C     Set overlap distribution batches
C
      CALL SETODBTCH(KIODBCLS,KIODBTCH,KRODBTCH,
     &               KFREE,LFREE,WORK,LWORK,IPRFCK3)
      MEMODB = LWORK - LFREE
#ifdef MOD_FMM
      IF ( DOFMM ) THEN 
        CALL SETMM(DMAT,NDMAT,WORK(KIODBCLS),WORK(KIODBTCH),
     &             WORK(KRODBTCH),WORK(KFREE),LFREE,IPRFCK3)
      END IF
#endif

CMAW **** SKIP EXPENSIVE 2-e NCL integrals via FTUV *****
C      WRITE(LUPRI,*) "*** maw: skipping all 2-e integrals in NCLFCK!"
C      RETURN

      CALL GETTIM(TIMFCK3(3),DUMTIM)

      NPARFTUV_NPAR = 1
     
      NPARFTUV_NSTEPODCLS = NPARFTUV_NPAR

      DO NPARFTUV_KSTARTODCLS = 1, NPARFTUV_NPAR

       IELECTR = 2

       CALL SETFTUV( KIFTUVCLS,KRFTUVCLS,KIFTUV,KRFTUV,KFTUV,
     &            KRFBTCH,KIFBTCH,WORK(KIODBCLS),WORK(KIODBTCH),
     &            DMAT,ISYMDM,NDMT,
     &            KFREE,LFREE,IELECTR,WORK,LWORK,
     &            NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS,IPRFCK3)

       CALL GETTIM(TIMFCK3(4),DUMTIM)

       CALL GETFCK3(FMAT,NDMT,WORK(KIODBCLS),WORK(KIODBTCH),
     &             WORK(KRODBTCH),WORK(KIFTUVCLS),WORK(KRFTUVCLS),
     &             WORK(KIFTUV),WORK(KRFTUV),WORK(KFTUV),WORK(KRFBTCH),
     &             WORK(KIFBTCH),WORK(KFREE),LFREE,IPRFCK3)

       CALL GETTIM(TIMFCK3(5),DUMTIM)

       IF ( IPRFCK3 .GE. 1 ) THEN 
        CALL HEADER('Settings in FCK3 :',-1)
        WRITE (LUPRI,'(/,6(5X,A,L3,/))') 
     &     'Screen prim. OD                 :', .NOT. NOSCRPRIM,
     &     'Skip classes of classical int.  :', .NOT. NOSKIPCL,
     &     'Skip classes by Cauchy-Schwarz  :', .NOT. NOSKIPCS,
     &     'Calculate only nonclassical part:', NCLFCK3,
     &     'Working in internally contracted basis :', CONTRACTED,
     &     'Branch-free code                :', BRFREE

        CALL HEADER('Thresholds in FCK3 :',-1)
        WRITE (LUPRI,'(/,2(5X,A,E16.2,/))') 
     &             'Primitive  OD overlap :', THRPOD,
     &             'Contracted OD overlap :', THRCOD
 
        IF ( NCLFCK3 ) THEN
         WRITE (LUPRI,'(/,2(5X,A,E16.2,/))') 
     &             'Nonclassical threshold:', THRNCL
        END IF

        IF ( SPEXP ) THEN
         WRITE (LUPRI,'(/,2(5X,A,E16.2,/))') 
     &             's-p exp.     threshold:', THRSPEXP 
        END IF

        CALL HEADER('Timing in FCK3 :',-1)
        WRITE (LUPRI,'(/,4(5X,A,F16.2,A,/))') 
     &             'SETAOBTCH :',TIMFCK3(2)-TIMFCK3(1),' seconds',
     &             'SETODBTCH :',TIMFCK3(3)-TIMFCK3(2),' seconds',
     &             'SETFTUV   :',TIMFCK3(4)-TIMFCK3(3),' seconds',
     &             'GETFCK3   :',TIMFCK3(5)-TIMFCK3(4),' seconds'
        MEMNDF = NDMT * 2 * ( NBASE * NBASE )
        CALL HEADER('Memory requirements in FCK3 :',-1)
        WRITE (LUPRI,'(A,F8.2,A)') 
     &  ' ( note - the conversion factor MW to MB used is ',
     &   WORD2MB * 1.0D6, ' ) '
        WRITE (LUPRI,'(/,7(5X,A,F16.2,A,/),5X,A,F16.2,A,/)') 
     &             'AOBTCH-com.b.:',MEMAOB   * WORD2MB,' MB',
     &             'ODBTCH       :',MEMODB   * WORD2MB,' MB',
     &             'FTUV - data  :',MEMFTUV1 * WORD2MB,' MB',
     &             'FTUV - temp  :',MEMFTUV2 * WORD2MB,' MB',
     &             'FCK3 - 1el   :',MEMFCK31 * WORD2MB,' MB',
     &             'FCK3 - 2el   :',MEMFCK32 * WORD2MB,' MB',
     &             'N*(DMAT+FMAT):',MEMNDF   * WORD2MB,' MB',
     &             'Available mem:',LWORK    * WORD2MB,' MB'

        CALL HEADER('Screening in FCK3 :',-1)

        A_ALL = STATISTICS(1) 
        A_CALCULATEDN = STATISTICS(2)/A_ALL
        A_CALCULATED0 = STATISTICS(3)/A_ALL
        A_SKIPPED1    = STATISTICS(4)/A_ALL
        A_SKIPPED2    = STATISTICS(5)/A_ALL

        WRITE (LUPRI,'(/,5(5X,A,F16.2,A,/))') 
     &         'Total number of primitive integrals in FCK3:',A_ALL,' ',
     &         ' - calculated nonzero     :', A_CALCULATEDN*100, ' % ',
     &         ' - calculated zero        :', A_CALCULATED0*100, ' % ',
     &         ' - skipped inside classes :', A_SKIPPED1*100, ' % ',
     &         ' - skipp by classes       :', A_SKIPPED2*100, ' % '

        CALL HEADER('Numbers of classes and batches in FCK3 :',-1)  
        WRITE (LUPRI,'(/,2(3(5X,A,I9,/),/))') 
     &         ' Number of AO   batches : ', NAOBTCH,
     &         ' Number of OD   batches : ', NODBTCH,
     &         ' Number of FTUV batches : ', NFTUVBTCH,
     &         ' Number of AO   classes : ', NAOBCLS,
     &         ' Number of OD   classes : ', NODBCLS,
     &         ' Number of FTUV classes : ', NFTUVCLS
       END IF 

      END DO

      CALL TIMER('FCK3  ',TIMSTR,TIMEND)
 
      CALL FLSHFO(LUPRI)

      CALL SETBFCFMM( .FALSE. ) 

      RETURN 
      END 

      SUBROUTINE GETFCK3(FMAT,NDMAT,IODBCLS,IODBTCH,RODBTCH,
     &                   IFTUVCLS,RFTUVCLS,IFTUV,RFTUV,FTUV,RFBTCH,
     &                   IFBTCH,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
      DIMENSION IODBCLS(5,*),IODBTCH(2,*),RODBTCH(5,*),
     &          IFTUVCLS(10,*),RFTUVCLS(4,*),IFTUV(*),RFTUV(*),
     &          FTUV(*),RFBTCH(5,*),IFBTCH(2,*),WORK(*)     
      LOGICAL GDER

      MAXDER = 0 
      GDER   = .FALSE. 

      JTOP  = 4*(MAXAONHK - 1)
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
     
      CALL PUT2HERTOP(JTOP,NRTOP)

C     Dimensions for CSQ
C     ==================
C
      CALL SPHDIM
C
C     Allocate
C     ========
C
      KINDHR = 1
      KINDSQ = KINDHR + ( JTOP3 - 1 ) / IRAT + 1
      KIODHR = KINDSQ + ( NRTOP - 1 ) / IRAT + 1
      KLMNVL = KIODHR + ( 8 * NRTOP - 1) / IRAT + 1
      KLMNSM = KLMNVL + ( 3 * MAXODKCK * MAXODNHK - 1) / IRAT + 1
      KCSQ   = KLMNSM + (16 * MAXODKCK * MAXODNHK - 1) / IRAT + 1
      KLAST  = KCSQ   + NCSQ1 * NCSQ2

      IF (KLAST .GT. LWORK) CALL STOPIT('FCK3','HERPRP',KLAST,LWORK)     
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

      CALL GETFCK31( FMAT,NDMAT,IODBCLS,IODBTCH,RODBTCH,
     &               IFTUVCLS,RFTUVCLS,IFTUV,RFTUV,FTUV,RFBTCH,
     &               IFBTCH,WORK(KINDHR),WORK(KIODHR),
     &               WORK(KLMNVL),WORK(KLMNSM),WORK(KCSQ),
     &               WORK(KLAST),LWRK,IPRINT)

      RETURN
      END

      SUBROUTINE GETFCK31(FMAT,NDMAT,IODBCLS,IODBTCH,RODBTCH,
     &                   IFTUVCLS,RFTUVCLS,IFTUV,RFTUV,FTUV,RFBTCH,
     &                   IFBTCH,INDHER,IODDHR,LMNPWR,LMNSYM,
     &                   CSQ,WORK,LWORK,IPRINT)
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
      DIMENSION IODBCLS(5,*), IODBTCH(2,*), RODBTCH(5,*),
     &          IFTUVCLS(10,*),RFTUVCLS(4,*),IFTUV(*),RFTUV(*),
     &          RFBTCH(5,*),IFBTCH(2,*),FTUV(*),CSQ(NCSQ1,NCSQ2),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHR(NRTOP,0:7),LMNPWR(MAXODKCK,MAXODNHK,3),
     &          LMNSYM(MAXODKCK,MAXODNHK,2,0:7),WORK(*), 
     &          KCCF(2,2)

      NOSKIPCL = ( ( .NOT. NCLFCK3 ) .OR. NOSKIPCL ) 
      IF ( BRFREE ) NCLFCK3 = .FALSE. 

      MAXNHK  = 4 * ( MAXAONHK - 1 ) + 1

      MAXPRF  = MAXODPRF * MAXFCLSPRF

      MAXNTUV = MAXNHK * ( MAXNHK + 1 ) * ( MAXNHK + 2 ) / 6

C
C   allocations  for the outer (one electron) loop 
C
      KETIJ    = 1
      KHRINT   = KETIJ    + MAXODPRF * MAXODNHK*(MAXAONHK*MAXAONHK)*3*2
      KPREFAB  = KHRINT   + MAXODPRF * MAXNTUV * NDMAT
      KIPRAB   = KPREFAB  + MAXODPRF
      KDATAP   = KIPRAB   + MAXODPRF * 2 / IRAT + 1 
      KCOORAB  = KDATAP   + MAXODPRF * 6 
      KEXPAB   = KCOORAB  + MAXODPRF * 3
      KHPI     = KEXPAB   + MAXODPRF * 3
      KPA      = KHPI     + MAXODPRF
      KPB      = KPA      + MAXODPRF * 3
      KFMATCP  = KPB      + MAXODPRF * 3
      KFMATSP  = KFMATCP  + MAXODPRF * MAXAOKCK * MAXAOKCK
      KFTMP    = KFMATSP  + MAXODPRF * MAXAOKHK * MAXAOKHK

C
C   allocations  for the inner (two-electron) loop 
C
      KRJ000   = KFTMP   + MAXAOKCK * MAXAOKHK
      KHRINT1  = KRJ000  + MAXPRF * MAXNHK
      KFACINT  = KHRINT1 + MAXPRF * MAXNTUV 
      KALPHA   = KFACINT + MAXPRF 
      KCOORPQ  = KALPHA  + MAXPRF 
      KFTUV2   = KCOORPQ + MAXPRF * 3
      KLAST    = KFTUV2  + MAXFCLSPRF * MAXODNHK**3

      MEMFCK31 = KRJ000 - 1 
      MEMFCK32 = KLAST  - KRJ000  

      IF (KLAST .GT. LWORK) CALL STOPIT('FCK3','FCK31',KLAST,LWORK)     
      LWRK   = LWORK - KLAST + 1

      IELECTR  = 1

      DO IABCLS = 1, NODBCLS

       ICENTA    = IODBCLS(1,IABCLS)
       ICENTB    = IODBCLS(2,IABCLS)
       NHKTAB    = IODBCLS(3,IABCLS)
       KABSTBTCH = IODBCLS(4,IABCLS)
       NABBTCH   = IODBCLS(5,IABCLS)

       JMAXAB = NHKTAB - 1
       NTUVAB = NHKTAB*(NHKTAB+1)*(NHKTAB+2)/6

       CALL GETR12( NSYMOPR, NSYMOPS, ICENTA, ICENTB, 1 )

       DO IABBTCH = KABSTBTCH, ( KABSTBTCH + NABBTCH - 1 ) 

        IAOBA = IODBTCH(1,IABBTCH)
        IAOBB = IODBTCH(2,IABBTCH)

        CALL GETAOB(NPRFA,NCTFA,NTYPA,NHKTA,JMAXA,KHKTA,KCKTA,
     &              KCSQA,KEXPA,KCCF(1,1),IAOBA)

        CALL GETAOB(NPRFB,NCTFB,NTYPB,NHKTB,JMAXB,KHKTB,KCKTB,
     &              KCSQB,KEXPB,KCCF(1,2),IAOBB)

        CALL SCREENPRIM12(NPRFAB,RODBTCH(1,IABBTCH),WORK(KIPRAB),.TRUE.,
     &                    NPRFA,NCTFA,NTYPA,NPRFB,NCTFB,NTYPB,
     &                    EXPBTCH(KEXPA),EXPBTCH(KEXPB),
     &                    CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &                    CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &                    R122(0,0),COOR1(1,0),COOR2(1,0),IPRINT)

        EXTPMAX = RODBTCH( 4, IABBTCH )
        SCRAB   = RODBTCH( 5, IABBTCH )

        CALL DZERO(WORK(KETIJ),NPRFAB*NHKTAB*NHKTA*NHKTB*3*2)

        CALL SETETIJ(WORK(KETIJ),WORK(KPREFAB),
     &               R12(1,0,0),R122(0,0),
     &               WORK(KIPRAB),EXPBTCH(KEXPA),EXPBTCH(KEXPB),
     &               WORK(KCOORAB),WORK(KEXPAB),WORK(KHPI),WORK(KPA),
     &               WORK(KPB),NPRFA,NHKTB,JMAXA,
     &               NPRFB,NHKTB,JMAXB,NPRFAB,NHKTAB,IELECTR,IPRINT)

        CALL SETP(WORK(KDATAP),WORK(KIPRAB),
     &            COOR1(1,0),COOR2(1,0),EXPBTCH(KEXPA),
     &            EXPBTCH(KEXPB),NPRFA,NPRFB,NPRFAB,IPRINT)

        CALL DZERO(WORK(KHRINT),NPRFAB*NTUVAB*NDMAT)

        DO IFCLS = 1, NFTUVCLS
 
         ICENTC      = IFTUVCLS(  9, IFCLS )
         ICENTD      = IFTUVCLS( 10, IFCLS ) 

         RCD         = RFTUVCLS(  1, IFCLS ) 
         EXTQMAX     = RFTUVCLS(  2, IFCLS ) 
         SCRCD       = RFTUVCLS(  4, IFCLS ) 

         IF ( NCLFCK3 ) 
     &    CALL CHECKRMIN(RMIN,RODBTCH(1,IABBTCH),ICENTC,ICENTD,RCD)

C
C   Skipping of the integrals for the whole FTUV class
C   ==================================================
C
C   NOSKIPCL - controls skipping of the classical integrals
C   NOSKIPCS - controls skipping based on Cauchy-Scharz parameters
C            

         NTFTUVPRIM  = IFTUVCLS(  2, IFCLS )

         NPRF0       = NPRFAB * NTFTUVPRIM
         NPRF        = NPRF0

         STATISTICS(1) = STATISTICS(1) + NPRF0

         IF ((( NOSKIPCL .OR. ( RMIN .LE. EXTPMAX + EXTQMAX ) ) .OR.
     &        ( BRFREE .AND. ( ( RMIN .LE. RPQMIN ) .OR.  
     &                       ( RMIN .LE. EXTPMAX + EXTQMAX ) ))) .AND.
     &        ( NOSKIPCS .OR. ( SCRAB * SCRCD .GT. THRCOD ) ) ) THEN

          NTFTUVBTCH  = IFTUVCLS(  1, IFCLS )
          KSIFTUV     = IFTUVCLS(  3, IFCLS )
          KSRFTUV     = IFTUVCLS(  4, IFCLS )
          KSFTUV      = IFTUVCLS(  5, IFCLS )
          JMAXF       = IFTUVCLS(  6, IFCLS )
          KSFTUVBTCH  = IFTUVCLS(  7, IFCLS )
          NTYPF       = IFTUVCLS(  8, IFCLS )

          JMAX = JMAXAB + JMAXF

          IF ( PACKFTUV ) THEN
            NTUVCD = ( JMAXF + 1 ) * ( JMAXF + 2 ) * ( JMAXF + 3 ) / 6
          ELSE 
            NTUVCD = ( JMAXF + 1 ) * ( JMAXF + 1 ) * ( JMAXF + 1 ) 
          END IF

          LSIZE = NDMAT * NTUVCD


          NTUV   = ( JMAX  + 1 ) * ( JMAX  + 2 ) * ( JMAX  + 3 ) / 6

          NTYP        = NTYPA + NTYPB + NTYPF 

          NTFTUVPRIM_TOT = NTFTUVPRIM

C          IF ( SCRFBTCH ) THEN
C            THRSCRF = THRCOD / SCRAB
C            CALL SCREENFBTCH( NTFTUVPRIM, IFBTCH(1,KSFBTCH),
C     &                   RFBTCH(1,KSFBTCH), THRSCRF, NTFTUVBTCH )
C          END IF
        
C
C         The main time consuming part !!!
C         ================================

          IF ( SKALLCL .AND.  NCLFCK3 ) THEN  
            CALL SETNPRF(NPRF,RODBTCH(1,IABBTCH),IFBTCH(1,KSFTUVBTCH),
     &                  RFBTCH(1,KSFTUVBTCH),NPRFAB,NTFTUVBTCH)
          END IF 

          IF ( NCLFCK3 .AND. ( NPRF .NE. NPRF0 ) ) THEN 
           CALL SETPQNCL1(WORK(KFACINT), WORK(KALPHA), WORK(KCOORPQ),
     &                   WORK(KFTUV2), NTFTUVPRIM, NPRF,IPQXYZ,
     &                   FTUV(KSFTUV), WORK(KDATAP), RODBTCH(1,IABBTCH),
     &                   IFTUV(KSIFTUV),RFTUV(KSRFTUV),
     &                   IFBTCH(1,KSFTUVBTCH), RFBTCH(1,KSFTUVBTCH), 
     &                   NPRFAB, NTFTUVBTCH, JMAX, NTYP, LSIZE, NTUV )    
           STATISTICS(2) = STATISTICS(2) + NPRF
          ELSE
           IF ( NCLFCK3 .AND. ( .NOT. SKALLCL ) ) THEN 
            CALL SETPQNCL2(WORK(KFACINT), WORK(KALPHA), WORK(KCOORPQ),
     &                   NTFTUVPRIM, NPRF,IPQXYZ,
     &                   WORK(KDATAP), RODBTCH(1,IABBTCH),
     &                   IFTUV(KSIFTUV),RFTUV(KSRFTUV),
     &                   IFBTCH(1,KSFTUVBTCH), RFBTCH(1,KSFTUVBTCH), 
     &                   NPRFAB, NTFTUVBTCH, JMAX, NTYP, LSIZE, NTUV )    
           ELSE 
            CALL SETPQ(WORK(KFACINT), WORK(KALPHA), WORK(KCOORPQ), NPRF,
     &               IPQXYZ, WORK(KDATAP), RODBTCH(1,IABBTCH),
     &               IFTUV(KSIFTUV),RFTUV(KSRFTUV),RFBTCH(1,KSFTUVBTCH), 
     &               NFTUVBTCH, NPRFAB, NTFTUVPRIM, JMAX, NTYP, NTUV )    
            STATISTICS(2) = STATISTICS(2) + NPRF
           END IF
          END IF

          IF ( NPRF .GT. 0 ) THEN 

          CALL ERIGAM(WORK(KRJ000),WORK(KFACINT),WORK(KALPHA),
     &                WORK(KCOORPQ),WORK(KLAST),LWRK,IPRINT)

          CALL ERITUV(WORK(KHRINT1),WORK(KRJ000),
     &                WORK(KCOORPQ),INDHER,IODDHR(1,IPQXYZ),
     &                WORK(KLAST),LWRK,IPRINT)

          IF ( NCLFCK3 .AND. ( NPRF .NE. NPRF0 ) ) THEN 
           IF ( PACKFTUV ) THEN
            CALL HR2FTUVNCL1(WORK(KHRINT),WORK(KHRINT1),WORK(KFTUV2),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMAXAB,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,
     &                  NPRFAB,NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP)
           ELSE
            CALL HR2FTUVNCL2(WORK(KHRINT),WORK(KHRINT1),WORK(KFTUV2),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMAXAB,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,
     &                  NPRFAB,NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP)
           END IF
          ELSE
           IF ( PACKFTUV ) THEN
            CALL HR2FTUV1(WORK(KHRINT),WORK(KHRINT1),FTUV(KSFTUV),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMAXAB,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,NPRFAB,
     &                  NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP )
           ELSE
            CALL HR2FTUV2(WORK(KHRINT),WORK(KHRINT1),FTUV(KSFTUV),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMAXAB,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,NPRFAB,
     &                  NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP )

           END IF
          END IF

          STATISTICS(4) = STATISTICS(4) + (NPRF0-NPRF)
C
C         End of the main time consuming part !!!
C         =======================================

         END IF

         ELSE 
          STATISTICS(5) = STATISTICS(5) + NPRF0
         END IF

        END DO 

        JA = JMAXA
        JB = JMAXB

        DO ICOMPA = 1, NTYPA
         IF ( NTYPA .EQ. 2 ) CALL SETJKKK(JA,KCSQA,KCKTA,KHKTA,ICOMPA)

         DO ICOMPB = 1, NTYPB 
          IF ( NTYPB .EQ. 2 ) CALL SETJKKK(JB,KCSQB,KCKTB,KHKTB,ICOMPB)           

           CALL ADD2FCK3(FMAT,WORK(KFMATCP),WORK(KFMATSP),
     &         WORK(KFTMP),NDMAT,
     &         WORK(KETIJ),WORK(KPREFAB),WORK(KIPRAB),
     &         WORK(KHRINT),INDHER,LMNPWR,IODDHR(1,IPQXYZ),
     &         CSQ(KCSQA,1),CSQ(KCSQB,1),
     &         CCFBTCH(KCCF(ICOMPA,1)),CCFBTCH(KCCF(ICOMPB,2)),
     &         ICOMPA,ICOMPB,NPRFAB,
     &         JA,JB,JMAXA, JMAXB, JMAXAB, NTUVAB,
     &         KHKTA, KHKTB, KCKTA, KCKTB,
     &         KCKMAX,NHKMAX,IAOBA,IAOBB,IPRINT )

         END DO
        END DO

       END DO 

      END DO

      RETURN 
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE SETP( DATAP, IPRAB, COORA, COORB,EXPBTA, EXPBTB, 
     &                 NPRFA, NPRFB, NPRFAB, IPRINT )
#include"implicit.h"
#include"priunit.h"
#include"pi.h"
      DIMENSION DATAP(6,NPRFAB), IPRAB(2,NPRFAB), COORA(3),COORB(3),
     &          EXPBTA(NPRFA), EXPBTB(NPRFB)

      FACTOR0 = 2.0D0*PI*PI*SQRT(PI) 

      DO I = 1, NPRFAB
       IPRFA = IPRAB(1,I)
       IPRFB = IPRAB(2,I)

       EXPA  = EXPBTA(IPRFA)
       EXPB  = EXPBTB(IPRFB)

       EXPP  = EXPA + EXPB
       EXPPI = 1.0D0 / ( EXPA + EXPB )

       DATAP(1,I) = EXPP
       DATAP(2,I) = EXPPI
       DATAP(3,I) = FACTOR0 * EXPPI 
       DATAP(4,I) = ( COORA(1)*EXPA + COORB(1)*EXPB ) * EXPPI
       DATAP(5,I) = ( COORA(2)*EXPA + COORB(2)*EXPB ) * EXPPI
       DATAP(6,I) = ( COORA(3)*EXPA + COORB(3)*EXPB ) * EXPPI
      END DO

      RETURN
      END 

C
C   --------------------------------------------------------------
C

      SUBROUTINE SETPQ( FACINT, ALPHA, COORPQ, NPRF, IPQXYZ,
     &                  DATAP, RODP, IFTUV, RFTUV, RODQ, 
     &                  NODQ, NPRFAB, NPRFCD, JMAX, NTYP, NTUV )
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
C      DIMENSION COORAB(3,NPRFAB),EXPABP(3,NPRFAB),
C     &          RFTUV(7,NPRFCD),IFTUV(2,NPRFCD),COORPQ(NPRF,3),
C     &          FACINT(*),ALPHA(*),RODP(4),RODQ(4,*),IPQ0(3)
      DIMENSION DATAP(6,*),RFTUV(7,*),IFTUV(2,*),COORPQ(NPRF,3),
     &          FACINT(*),ALPHA(*),RODP(4),RODQ(4,*),IPQ0(3)

      INTEGER X

      I=0 

      DO IPRFAB=1,NPRFAB

       EXPP  = DATAP(1,IPRFAB)
       EXPPI = DATAP(2,IPRFAB)
       FACPI = DATAP(3,IPRFAB)
       PX    = DATAP(4,IPRFAB)
       PY    = DATAP(5,IPRFAB)
       PZ    = DATAP(6,IPRFAB)

       DO IPRFCD=1,NPRFCD
        I = I + 1

        EXPQ  = RFTUV(1,IPRFCD)
        EXPQI = RFTUV(2,IPRFCD)
        QX    = RFTUV(3,IPRFCD)
        QY    = RFTUV(4,IPRFCD)
        QZ    = RFTUV(5,IPRFCD)
      
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  The definition of the alpha in the equations
C  is alpha_def = p*q/(p+q)
C  
C  However we have to put alpha = -2*alpha_def in order
C  to be compatible with way how with the erigam uses 
C  the alpha  
C
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ALPHA(I)  = - 2.0D0 / (EXPPI + EXPQI) 
        FACINT(I)= FACPI*EXPQI/SQRT(EXPP+EXPQ)

        COORPQ(I,1)= PX - QX
        COORPQ(I,2)= PY - QY
        COORPQ(I,3)= PZ - QZ

       END DO
      END DO


      IF ( NTYP .GT. 3 ) THEN 
       IPQXYZ  = 0
       IPQ0(1) = 0
       IPQ0(2) = 0
       IPQ0(3) = 0
      ELSE 
       CALL LOCSYM( COORPQ, NPRF, 1, IPQ0, IPQXYZ )
      END IF

      CALL PUT2ERICOM2( NPRF, IPQ0, IPQXYZ, JMAX, NTUV )

      RETURN
      END


C
C   --------------------------------------------------------------
C

C
C  This is not an efficient way - now only for first tests :
C

      SUBROUTINE SETNPRF(NPRF,RODBTCH,IFBTCH,RFBTCH,
     &                   NPRFAB,NFBTCH)
#include "implicit.h"
#include "priunit.h"
#include "fck3thr.h"
      DIMENSION RODBTCH(4),IFBTCH(2,*),RFBTCH(5,*)

      NPRFCDSUM = 0 

      PXBTCH = RODBTCH( 1 )
      PYBTCH = RODBTCH( 2 )
      PZBTCH = RODBTCH( 3 )
      EXTPBTCH = RODBTCH( 4 )

      DO IBTCH = 1, NFBTCH
       QXBTCH = RFBTCH( 1, IBTCH )
       QYBTCH = RFBTCH( 2, IBTCH )
       QZBTCH = RFBTCH( 3, IBTCH )
       EXTQBTCH = RFBTCH( 4, IBTCH )
       NPRFCD = IFBTCH( 1, IBTCH )

       RPQ = SQRT( ( QXBTCH - PXBTCH ) * ( QXBTCH - PXBTCH )
     &           + ( QYBTCH - PYBTCH ) * ( QYBTCH - PYBTCH )
     &           + ( QZBTCH - PZBTCH ) * ( QZBTCH - PZBTCH ) ) 

       IF ( RPQ .LE. EXTPBTCH + EXTQBTCH ) 
     &     NPRFCDSUM  = NPRFCDSUM + NPRFCD 

      END DO

      NPRF = NPRFAB * NPRFCDSUM

      RETURN
      END


C
C   --------------------------------------------------------------
C
      SUBROUTINE SETPQNCL1(FACINT,ALPHA,COORPQ,FTUV2,NPRFCD_TOT,NPRF,
     &             IPQXYZ,FTUV1,DATAP,RODBTCH,IFTUV,RFTUV,
     &             IFBTCH,RFBTCH,NPRFAB, NFBTCH, JMAX, NTYP, LSIZE,
     &             NTUV )
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
      DIMENSION COORAB(3,NPRFAB),EXPABP(3,NPRFAB),COORPQ(NPRF,3),
     &   RFTUV(7,NPRFCD_TOT),IFTUV(2,NPRFCD_TOT),FACINT(*),ALPHA(*),
     &   RODBTCH(4),IFBTCH(2,*),RFBTCH(5,*),IPQ0(3),
     &   FTUV1( NPRFCD_TOT, LSIZE ), FTUV2( NPRFCD_TOT, LSIZE ),
     &   DATAP(6,NPRFAB)

      INTEGER X

      I0   = 0
      I1   = 1
      I2   = 1
      IPRFCD= 1

      PXBTCH = RODBTCH( 1 )
      PYBTCH = RODBTCH( 2 )
      PZBTCH = RODBTCH( 3 )
      EXTPBTCH = RODBTCH( 4 )

      DO IBTCH = 1, NFBTCH
       QXBTCH = RFBTCH( 1, IBTCH )
       QYBTCH = RFBTCH( 2, IBTCH )
       QZBTCH = RFBTCH( 3, IBTCH )
       EXTQBTCH = RFBTCH( 4, IBTCH )
       NPRFCD = IFBTCH( 1, IBTCH )

       RPQ = SQRT( ( QXBTCH - PXBTCH ) * ( QXBTCH - PXBTCH )
     &           + ( QYBTCH - PYBTCH ) * ( QYBTCH - PYBTCH )
     &           + ( QZBTCH - PZBTCH ) * ( QZBTCH - PZBTCH ) )

       IF ( RPQ .LE. EXTPBTCH + EXTQBTCH ) THEN 

        DO IPRFAB=1,NPRFAB
         EXPP  = DATAP(1,IPRFAB)
         EXPPI = DATAP(2,IPRFAB)
         FACPI = DATAP(3,IPRFAB)
         PX    = DATAP(4,IPRFAB)
         PY    = DATAP(5,IPRFAB)
         PZ    = DATAP(6,IPRFAB)

         I = I0 + IPRFAB

         DO IPRFCDOFF=0,NPRFCD-1
          EXPQ  = RFTUV(1,IPRFCD+IPRFCDOFF)
          EXPQI = RFTUV(2,IPRFCD+IPRFCDOFF)
          QX    = RFTUV(3,IPRFCD+IPRFCDOFF)
          QY    = RFTUV(4,IPRFCD+IPRFCDOFF)
          QZ    = RFTUV(5,IPRFCD+IPRFCDOFF)

          ALPHA(I)  = - 2.0D0 / (EXPPI + EXPQI) 
          FACINT(I)= FACPI*EXPQI/SQRT(EXPP+EXPQ)

          COORPQ(I,1)= PX - QX
          COORPQ(I,2)= PY - QY
          COORPQ(I,3)= PZ - QZ

          I = I + NPRFAB     
         END DO
        END DO

        I0 = I0 + NPRFAB*NPRFCD

        DO IOFF = 0, NPRFCD-1
         DO J=1, LSIZE
          FTUV2( I2+IOFF, J ) = FTUV1( I1+IOFF, J )
         END DO 
        END DO

        I2 = I2 + NPRFCD
       END IF

       I1 = I1 + NPRFCD
       IPRFCD = IPRFCD + NPRFCD

      END DO

      NPRFCD_TOT = ( I2 - 1 )

      IF ( NPRF .NE. NPRFAB * ( I2 - 1 ) ) THEN 
       CALL QUIT('Error (FCK3:SETPQNCL): NPRF .NE. NPRFAB * (I2-1) !!!')
      END IF

      IF ( NTYP .GT. 3 ) THEN 
       IPQXYZ  = 0
       IPQ0(1) = 0
       IPQ0(2) = 0
       IPQ0(3) = 0
      ELSE 
       CALL LOCSYM( COORPQ, NPRF, 1, IPQ0, IPQXYZ )
      END IF

      CALL PUT2ERICOM2( NPRF, IPQ0, IPQXYZ, JMAX, NTUV )

      RETURN
      END

C
C   --------------------------------------------------------------
C
      SUBROUTINE SETPQNCL2(FACINT,ALPHA,COORPQ,NPRFCD_TOT,NPRF,
     &             IPQXYZ,DATAP,RODBTCH,IFTUV,RFTUV,
     &             IFBTCH,RFBTCH,NPRFAB, NFBTCH, JMAX, NTYP, LSIZE,
     &             NTUV )
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
      DIMENSION COORAB(3,NPRFAB),EXPABP(3,NPRFAB),COORPQ(NPRF,3),
     &   RFTUV(7,NPRFCD_TOT),IFTUV(2,NPRFCD_TOT),FACINT(*),ALPHA(*),
     &   RODBTCH(4),IFBTCH(2,*),RFBTCH(5,*),IPQ0(3),
     &   DATAP(6,NPRFAB)

      INTEGER X

      IPRFCD= 1

      PXBTCH = RODBTCH( 1 )
      PYBTCH = RODBTCH( 2 )
      PZBTCH = RODBTCH( 3 )
      EXTPBTCH = RODBTCH( 4 )

      DO IBTCH = 1, NFBTCH
       QXBTCH = RFBTCH( 1, IBTCH )
       QYBTCH = RFBTCH( 2, IBTCH )
       QZBTCH = RFBTCH( 3, IBTCH )
       EXTQBTCH = RFBTCH( 4, IBTCH )
       NPRFCD = IFBTCH( 1, IBTCH )

       RPQ = SQRT( ( QXBTCH - PXBTCH ) * ( QXBTCH - PXBTCH )
     &           + ( QYBTCH - PYBTCH ) * ( QYBTCH - PYBTCH )
     &           + ( QZBTCH - PZBTCH ) * ( QZBTCH - PZBTCH ) )

       IF ( RPQ .LE. EXTPBTCH + EXTQBTCH ) THEN 

        DO IPRFAB=1,NPRFAB
         EXPP  = DATAP(1,IPRFAB)
         EXPPI = DATAP(2,IPRFAB)
         FACPI = DATAP(3,IPRFAB)
         PX    = DATAP(4,IPRFAB)
         PY    = DATAP(5,IPRFAB)
         PZ    = DATAP(6,IPRFAB)

         I = IPRFCD + (IPRFAB-1) * NPRFCD_TOT

         DO IPRFCDOFF=0,NPRFCD-1
          EXPQ  = RFTUV(1,IPRFCD+IPRFCDOFF)
          EXPQI = RFTUV(2,IPRFCD+IPRFCDOFF)
          QX    = RFTUV(3,IPRFCD+IPRFCDOFF)
          QY    = RFTUV(4,IPRFCD+IPRFCDOFF)
          QZ    = RFTUV(5,IPRFCD+IPRFCDOFF)

          ALPHA(I)  = - 2.0D0 / (EXPPI + EXPQI) 
          FACINT(I)= FACPI*EXPQI/SQRT(EXPP+EXPQ)

          COORPQ(I,1)= PX - QX
          COORPQ(I,2)= PY - QY
          COORPQ(I,3)= PZ - QZ

          I = I + 1
         END DO
        END DO
        STATISTICS(2) = STATISTICS(2) + NPRFAB*NPRFCD
       ELSE
        DO IPRFAB=1,NPRFAB
         I0 = IPRFCD + (IPRFAB-1) * NPRFCD_TOT
         DO I = I0, I0 + NPRFCD-1
          FACINT(I) = 0.0D0
          ALPHA(I)  = 0.0D0
          COORPQ(I,1)= 0.0D0
          COORPQ(I,2)= 0.0D0
          COORPQ(I,3)= 0.0D0
         END DO
        END DO
        STATISTICS(3) = STATISTICS(3) + NPRFAB*NPRFCD
       END IF

       IPRFCD = IPRFCD + NPRFCD
      END DO

      IF ( NTYP .GT. 3 ) THEN 
       IPQXYZ  = 0
       IPQ0(1) = 0
       IPQ0(2) = 0
       IPQ0(3) = 0
      ELSE 
       CALL LOCSYM( COORPQ, NPRF, 1, IPQ0, IPQXYZ )
      END IF

      CALL PUT2ERICOM2( NPRF, IPQ0, IPQXYZ, JMAX, NTUV )

      RETURN
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE PUT2HERTOP(JTOPVAL,NRTOPVAL)
#include"implicit.h"
#include"hertop.h"

      JTOP  = JTOPVAL
      NRTOP = NRTOPVAL

      RETURN
      END


      SUBROUTINE PUT2ERICOM2( NPPXVAL, IPQ0VAL, IPQXYZVAL,
     &                        JMAXVAL, NTUVVAL )
#include"implicit.h"
#include"ericom.h"
      DIMENSION IPQ0VAL(3)

      NPPX    = NPPXVAL
      IPQ0(1) = IPQ0VAL(1)
      IPQ0(2) = IPQ0VAL(2)
      IPQ0(3) = IPQ0VAL(3)
      IPQXYZ  = IPQXYZVAL
      JMAX    = JMAXVAL
      NTUV    = NTUVVAL

      RETURN
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE HR2FTUV( HRINT, HRINT1, FTUV,  INDHER, IODDHR, NDMAT,
     &                    JMAXAB, JMAXCD, NTUVAB, NTUVCD, NTUV, 
     &                    NPRFAB, NPRFCD, NPRFCD_TOT,NPRF, JTOP )
#include "implicit.h"

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV
      DIMENSION HRINT(NPRFAB,NTUVAB,NDMAT),HRINT1(NPRF,NTUV),
     &       FTUV(NPRFCD_TOT,NTUVCD,NDMAT), 
     &       INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(*)

      INTEGER X

      DO IDMAT = 1, NDMAT
       DO JAB = 0, JMAXAB
        DO TAB = JAB, 0, -1
         DO UAB = JAB - TAB, 0, -1
          VAB = JAB - TAB - UAB
          TUVAB = INDHER(TAB,UAB,VAB)
          DO JCD=0,JMAXCD
           DO TCD = JCD, 0, -1 
            DO UCD = JCD - TCD, 0, -1
             VCD = JCD - TCD - UCD
             TUVCD = INDHER(TCD,UCD,VCD)
             TUV   = INDHER(TAB+TCD,UAB+UCD,VAB+VCD)

             IF (IODDHR(TUV) .EQ. 0) THEN
C               CALL DGEMV('n',NPRFAB,NPRFCD,1.0D0,
C     &                   HRINT1(1,TUV),NPRFCD,
C     &                   FTUV(1,TCD,UCD,VCD,IDMAT),1,
C     &                   HRINT(1,TUVAB,IDMAT),1)
              I=1
              DO IPRFAB = 1, NPRFAB 
               FSUM=0.0D0
               DO IPRFCD = 1, NPRFCD
                FSUM = FSUM + HRINT1(I,TUV)*FTUV(IPRFCD,TUVCD,IDMAT)
                I = I+1
               END DO
               HRINT(IPRFAB,TUVAB,IDMAT)=HRINT(IPRFAB,TUVAB,IDMAT)+FSUM
              END DO
             END IF

            END DO
           END DO
          END DO

         END DO
        END DO
       END DO
      END DO

      RETURN
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE HR2FTUV1( HRINT, HRINT1, FTUV,  INDHER, IODDHR, NDMAT,
     &                    JMAXAB, JMAXCD, NTUVAB, NTUVCD, NTUV, 
     &                    NPRFAB, NPRFCD, NPRFCD_TOT,NPRF, JTOP )
#include "implicit.h"

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV
      DIMENSION HRINT(NPRFAB,NTUVAB,NDMAT),HRINT1(NPRF,NTUV),
     &       FTUV(NPRFCD_TOT,NTUVCD,NDMAT), 
     &       INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(*)

      INTEGER X

      DO IDMAT = 1, NDMAT
       DO JAB = 0, JMAXAB
        DO TAB = JAB, 0, -1
         DO UAB = JAB - TAB, 0, -1
          VAB = JAB - TAB - UAB
          TUVAB = INDHER(TAB,UAB,VAB)
          DO JCD=0,JMAXCD
           DO TCD = JCD, 0, -1 
            DO UCD = JCD - TCD, 0, -1
             VCD = JCD - TCD - UCD
             TUVCD = INDHER(TCD,UCD,VCD)
             TUV   = INDHER(TAB+TCD,UAB+UCD,VAB+VCD)

             IF (IODDHR(TUV) .EQ. 0) THEN
C               CALL DGEMV('n',NPRFAB,NPRFCD,1.0D0,
C     &                   HRINT1(1,TUV),NPRFCD,
C     &                   FTUV(1,TCD,UCD,VCD,IDMAT),1,
C     &                   HRINT(1,TUVAB,IDMAT),1)
              I=1
              DO IPRFAB = 1, NPRFAB 
               FSUM=0.0D0
               DO IPRFCD = 1, NPRFCD
                FSUM = FSUM + HRINT1(I,TUV)*FTUV(IPRFCD,TUVCD,IDMAT)
                I = I+1
               END DO
               HRINT(IPRFAB,TUVAB,IDMAT)=HRINT(IPRFAB,TUVAB,IDMAT)+FSUM
              END DO
             END IF

            END DO
           END DO
          END DO

         END DO
        END DO
       END DO
      END DO

      RETURN
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE HR2FTUV2( HRINT, HRINT1, FTUV,  INDHER, IODDHR, NDMAT,
     &                    JMAXAB, JMAXCD, NTUVAB, NTUVCD, NTUV, 
     &                    NPRFAB, NPRFCD, NPRFCD_TOT,NPRF, JTOP )
#include "implicit.h"

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV
      DIMENSION HRINT(NPRFAB,NTUVAB,NDMAT),HRINT1(NPRF,NTUV),
     &       FTUV(NPRFCD_TOT,0:JMAXCD,0:JMAXCD,0:JMAXCD,NDMAT), 
     &       INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(*)

      INTEGER X

      DO IDMAT = 1, NDMAT
       DO JAB = 0, JMAXAB
        DO TAB = JAB, 0, -1
         DO UAB = JAB - TAB, 0, -1
          VAB = JAB - TAB - UAB
          TUVAB = INDHER(TAB,UAB,VAB)
          DO JCD=0,JMAXCD
           DO TCD = JCD, 0, -1 
            DO UCD = JCD - TCD, 0, -1
             VCD = JCD - TCD - UCD
             TUV   = INDHER(TAB+TCD,UAB+UCD,VAB+VCD)

             IF (IODDHR(TUV) .EQ. 0) THEN
C               CALL DGEMV('n',NPRFAB,NPRFCD,1.0D0,
C     &                   HRINT1(1,TUV),NPRFCD,
C     &                   FTUV(1,TCD,UCD,VCD,IDMAT),1,
C     &                   HRINT(1,TUVAB,IDMAT),1)
              I=1
              DO IPRFAB = 1, NPRFAB 
               FSUM=0.0D0
               DO IPRFCD = 1, NPRFCD
                FSUM = FSUM + HRINT1(I,TUV)*
     &             FTUV(IPRFCD,TCD,UCD,VCD,IDMAT)
                I = I+1
               END DO
               HRINT(IPRFAB,TUVAB,IDMAT)=HRINT(IPRFAB,TUVAB,IDMAT)+FSUM
              END DO
             END IF

            END DO
           END DO
          END DO

         END DO
        END DO
       END DO
      END DO

      RETURN
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE HR2FTUVNCL1(HRINT,HRINT1,FTUV,INDHER,IODDHR,NDMAT,
     &                    JMAXAB, JMAXCD, NTUVAB, NTUVCD, NTUV,
     &                    NPRFAB, NPRFCD, NPRFCD_TOT,NPRF, JTOP )
#include "implicit.h"

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV
      DIMENSION HRINT(NPRFAB,NTUVAB,NDMAT),HRINT1(NPRF,NTUV),
     &       FTUV(NPRFCD_TOT,NTUVCD,NDMAT), 
     &       INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(*)

      INTEGER X

      DO IDMAT = 1, NDMAT
       DO JAB = 0, JMAXAB
        DO TAB = JAB, 0, -1
         DO UAB = JAB - TAB, 0, -1
          VAB = JAB - TAB - UAB
          TUVAB = INDHER(TAB,UAB,VAB)
          DO JCD=0,JMAXCD
           DO TCD = JCD, 0, -1 
            DO UCD = JCD - TCD, 0, -1
             VCD = JCD - TCD - UCD
             TUVCD = INDHER(TCD,    UCD,    VCD)
             TUV   = INDHER(TAB+TCD,UAB+UCD,VAB+VCD)

             IF (IODDHR(TUV) .EQ. 0) THEN
              DO IPRFAB = 1, NPRFAB 
               FSUM=0.0D0
               I = IPRFAB
               DO IPRFCD = 1, NPRFCD
                FSUM = FSUM + 
     &                 HRINT1(I,TUV)*FTUV(IPRFCD,TUVCD,IDMAT)
                I = I+NPRFAB
               END DO
               HRINT(IPRFAB,TUVAB,IDMAT)=HRINT(IPRFAB,TUVAB,IDMAT)+FSUM
              END DO
             END IF

            END DO
           END DO
          END DO

         END DO
        END DO
       END DO
      END DO

      RETURN
      END


C
C   --------------------------------------------------------------
C

      SUBROUTINE HR2FTUVNCL2(HRINT,HRINT1,FTUV,INDHER,IODDHR,NDMAT,
     &                    JMAXAB, JMAXCD, NTUVAB, NTUVCD, NTUV,
     &                    NPRFAB, NPRFCD, NPRFCD_TOT,NPRF, JTOP )
#include "implicit.h"

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV
      DIMENSION HRINT(NPRFAB,NTUVAB,NDMAT),HRINT1(NPRF,NTUV),
     &       FTUV(NPRFCD_TOT,0:JMAXCD,0:JMAXCD,0:JMAXCD,NDMAT), 
     &       INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(*)

      INTEGER X

      DO IDMAT = 1, NDMAT
       DO JAB = 0, JMAXAB
        DO TAB = JAB, 0, -1
         DO UAB = JAB - TAB, 0, -1
          VAB = JAB - TAB - UAB
          TUVAB = INDHER(TAB,UAB,VAB)
          DO JCD=0,JMAXCD
           DO TCD = JCD, 0, -1 
            DO UCD = JCD - TCD, 0, -1
             VCD = JCD - TCD - UCD
             TUV   = INDHER(TAB+TCD,UAB+UCD,VAB+VCD)

             IF (IODDHR(TUV) .EQ. 0) THEN
              DO IPRFAB = 1, NPRFAB 
               FSUM=0.0D0
               I = IPRFAB
               DO IPRFCD = 1, NPRFCD
                FSUM = FSUM + 
     &                 HRINT1(I,TUV)*FTUV(IPRFCD,TCD,UCD,VCD,IDMAT)
                I = I+NPRFAB
               END DO
               HRINT(IPRFAB,TUVAB,IDMAT)=HRINT(IPRFAB,TUVAB,IDMAT)+FSUM
              END DO
             END IF

            END DO
           END DO
          END DO

         END DO
        END DO
       END DO
      END DO

      RETURN
      END

C
C   --------------------------------------------------------------
C

C
C  FMAT:    - AO is in Spherical Contracted (full)
C        CP - Cartesian Primitive (submatrix)
C        SP - Spherical Primitive (submatrix)
C
C  So far we do not take care about symmetry.
C  FMAT is assumed to be full dimension without any 
C  symmetry reduction
C
      SUBROUTINE ADD2FCK3(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQA,CSQB,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &        JA, JB, JMAXA, JMAXB, JMAXAB, NTUVAB,
     &        KHKTA, KHKTB, KCKTA, KCKTB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,IPRINT )
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
C
      INTEGER T,U,V,TUV
      DIMENSION FMAT(NBASE,NBASE,NDMAT),FMATCP(KCKTA,KCKTB,NPRFAB),
     &   FMATSP(KHKTA,KHKTB,NPRFAB),FTMP(KCKTA,KHKTB),
     &   ETIJ(NPRFAB,0:JMAXAB,0:JMAXA,0:JMAXB,3,*),
     &   PREFAB(NPRFAB),HRINT(NPRFAB,NTUVAB,NDMAT),IPRAB(2,NPRFAB),
     &   INDHER(0:JTOP,0:JTOP,0:JTOP),IODDHR(*),
     &   LMNPWR(MAXODKCK,MAXODNHK,3),
     &   CSQA(KHKTA,KCKTA), CSQB(KHKTB,KCKTB),
     &   CCFA(NPRFA,NCTFA), CCFB(NPRFB,NCTFB)

      IF ( IPRINT .GE. 10 ) THEN
        CALL HEADER('Output from ADDFCK3',-1)

        WRITE(LUPRI,*) 'JMAXA, KHKTA, KCKTA, NPRFA, NCTFA',
     &                  JMAXA, KHKTA, KCKTA, NPRFA, NCTFA
        WRITE(LUPRI,*) 'JMAXB, KHKTB, KCKTB, NPRFB, NCTFB',
     &       JMAXB, KHKTB, KCKTB, NPRFB, NCTFB
        WRITE(LUPRI,*) 'JMAXAB, NTUVAB, NDMAT, NBASE ',
     &                  JMAXAB, NTUVAB, NDMAT, NBASE 
      END IF 

      CALL DZERO(FMATCP, KCKTA*KCKTB*NPRFAB )
      CALL DZERO(FMATSP, KHKTA*KHKTB*NPRFAB )
      CALL DZERO(FTMP, KCKTA*KHKTB )

      IDMAT = 1 

      IREPA  = 0
      IREPB  = 0

      DO IPRFAB = 1, NPRFAB

       IPRFA = IPRAB(1,IPRFAB) 
       IPRFB = IPRAB(2,IPRFAB) 

       DO JAB = 0,JMAXAB
        DO T = 0, JAB
         DO U = 0, JAB - T 
          V = JAB - T - U
          TUV = INDHER(T,U,V)

           DO KCA = 1, KCKTA
            LA = LMNPWR(KCA,JA+1,1)
            MA = LMNPWR(KCA,JA+1,2)
            NA = LMNPWR(KCA,JA+1,3)
 
            DO KCB = 1, KCKTB
             LB = LMNPWR(KCB,JB+1,1)
             MB = LMNPWR(KCB,JB+1,2)
             NB = LMNPWR(KCB,JB+1,3)

             FAC = ETIJ(IPRFAB,T,LA,LB,1,1)*
     &         ETIJ(IPRFAB,U,MA,MB,2,1)*ETIJ(IPRFAB,V,NA,NB,3,1)
 
             FMATCP(KCA,KCB,IPRFAB) = 
     &            FMATCP(KCA,KCB,IPRFAB) + FAC*PREFAB(IPRFAB)
     &                      * HRINT(IPRFAB,TUV,IDMAT)
            END DO
           END DO

         END DO
        END DO
       END DO
     
       CALL DGEMM('N','T', KCKTA, KHKTB, KCKTB, 1.0D0,
     &             FMATCP(1,1,IPRFAB), KCKTA, CSQB, KHKTB, 
     &             0.0D0, FTMP, KCKTA  )

       CALL DGEMM('N','N', KHKTA, KHKTB, KCKTA, 1.0D0,
     &             CSQA, KHKTA, FTMP, KCKTA, 
     &             0.0D0, FMATSP(1,1,IPRFAB), KHKTA  )      

       IF ( IPRINT .GE. 15 ) THEN
        CALL HEADER('FMATCP',-1)
        WRITE( LUPRI, * ) 'IPRFAB, KCKTA, KCKTB: ',
     &                     IPRFAB, KCKTA, KCKTB   
        CALL OUTPUT(FMATCP(1,1,IPRFAB),1,KCKTA,1,KCKTB,
     &                      KCKTA,KCKTB,1,LUPRI)
        CALL HEADER('FTMP',-1)
        CALL OUTPUT(FTMP(1,1),1,KCKTA,1,KHKTB,
     &                     KCKTA,KHKTB,1,LUPRI)
        CALL HEADER('FMATSP',-1)
        WRITE( LUPRI, * ) 'IPRFAB, KHKTA, KHKTB: ',
     &                     IPRFAB, KHKTA, KHKTB   
        CALL OUTPUT(FMATSP(1,1,IPRFAB),1,KHKTA,1,KHKTB,
     &                      KHKTA,KHKTB,1,LUPRI)
       END IF

       DO ICTFA = 1, NCTFA
        DO ICTFB = 1, NCTFB

         IAOA = INDXBTCH(KNDXBTCH(ICOMPA,IAOBA) - 1 + ICTFA)
         IAOB = INDXBTCH(KNDXBTCH(ICOMPB,IAOBB) - 1 + ICTFB)

C
C Here it we should probably put summation over 
C  FMATSP(ISPA,ISPB,IPRFAB) separately outside 
C 

C         IF ( (IAOA .EQ. IAOB ) .OR. ( IAOBA .EQ. IAOBB )) THEN 
         IF ( ( IAOBA .EQ. IAOBB )) THEN 
          DO ISPA=1, KHKTA
           DO ISPB=1, KHKTB
             FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMAT) = 
     &            FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMAT) 
     &          + FMATSP(ISPA,ISPB,IPRFAB) 
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
           END DO
          END DO
         ELSE
          DO ISPA=1, KHKTA
           DO ISPB=1, KHKTB
              FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMAT) = 
     &            FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMAT) 
     &          + FMATSP(ISPA,ISPB,IPRFAB) 
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
              FMAT(IAOB+ISPB-1,IAOA+ISPA-1,IDMAT) = 
     &            FMAT(IAOB+ISPB-1,IAOA+ISPA-1,IDMAT) 
     &          + FMATSP(ISPA,ISPB,IPRFAB) 
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
           END DO
          END DO
         END IF


        END DO
       END DO

      END DO
 
      RETURN
      END 

      SUBROUTINE SCREENFBTCH( NACTFTUVPRIM, IFBTCH, RFBTCH, 
     &                        THRSCRF, NFBTCH )
#include"implicit.h"
#include"priunit.h"
      DIMENSION IFBTCH(2,NFBTCH), RFBTCH(5,NFBTCH)

      IBTCH         = 1
      SCRF          = RFBTCH(5,IBTCH)
      NACTFTUVPRIM  = IFBTCH(1,IBTCH)

      DO WHILE( ( SCRF .GT. THRSCRF ) .AND. (IBTCH .LT.NFBTCH))
       IBTCH = IBTCH + 1
       SCRF         = RFBTCH(5,IBTCH)        
       NACTFTUVPRIM = NACTFTUVPRIM + IFBTCH(1,IBTCH)
      END DO

      RETURN
      END

      SUBROUTINE FCK3INI
#include "implicit.h"
#include "cbifck3.h"
#include "fck3thr.h"

      NOSCRPRIM = .FALSE.
      NOSKIPCL  = .FALSE.
      NOSKIPCS  = .FALSE.
      SCRFBTCH  = .FALSE.
Cdirac      FIXTHR    = .FALSE.
Cdirac ..setting this true since IFTHRS is not set
      FIXTHR    = .TRUE.
      BRFREE    = .FALSE.
      FULL_J    = .FALSE.

      SPEXP     = .TRUE.
      PACKFTUV  = .TRUE.
      SKALLCL   = .TRUE.
      CONTRACTED= .TRUE.

      IPRFCK3   = 0 

      THRNCL    = THRFCK3_DEFAULT 
      THRPOD    = THRFCK3_DEFAULT 
      THRCOD    = THRFCK3_DEFAULT 
      THRSPEXP  = THRFCK3_DEFAULT 

      CALL SETBFCFMM( BRFREE ) 

      RETURN
      END

      SUBROUTINE FCK3INP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "cbifck3.h"
      PARAMETER ( NTABLE = 12 )
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DATA TABLE /'.FIXTHR','.PRINT ','.NOSPEX','.NOSCRP',
     &            '.NOSKCL','.NOSKCS','.SCRFBT','.NOPACK',
     &            '.NOSKAL','.FULL_J','.UNCONT','.BRFREE'/
C
C
CTROND
      WRITE(6,*) '..inside FCK3INP'      
      IPRINT=IPRFCK3
      NEWDEF = (WORD .EQ. '*JENGIN')
      IF(NEWDEF) THEN
        WORD1 = WORD
 100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)

            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
                  END IF
 200          CONTINUE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FCK3INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FCK3INP.')
 1    CONTINUE
C&&& FIXTHR: No dynamical threshold for screening
        FIXTHR     = .TRUE.
      GO TO 100
 2    CONTINUE
C &&& PRINT: Print level
        READ (LUCMD, '(I5)') IPRFCK3
      GO TO 100
 3    CONTINUE
C &&& NOSPEX: Do not exploit common s/p-exponents
        SPEXP      = .FALSE.
      GO TO 100
 4    CONTINUE 
C &&& NOSCRP: No screening of primitives in the overlap distribution
         NOSCRPRIM = .TRUE.
      GO TO 100
 5    CONTINUE
C &&& NOSKCL: Do not skip classes (ODB|FTUV) based on 
C             classical/non-classical contributions (for FMM)
         NOSKIPCL  = .TRUE.
      GO TO 100
 6    CONTINUE
C &&& NOSKCS: Do not skip classes (ODB|FTUV) based on 
C             Cauchy-Schwartz
         NOSKIPCS  = .TRUE.
      GO TO 100
 7    CONTINUE
C &&& SCRFBT: 
         SCRFBTCH  = .TRUE.
      GO TO 100
 8    CONTINUE
C &&& NOPACK: Do not pack Ftuv array
         PACKFTUV  = .FALSE.
      GO TO 100  
 9    CONTINUE
C &&& NOSKAL: Do not skip all classes (ODB|FTUV) based on 
C             Cauchy-Schwartz 
         SKALLCL   = .FALSE.
      GO TO 100  
 10   CONTINUE
C &&& FULL_J: combine 1-electron and 2-electron multipoles
C             (under construction)
         FULL_J    = .TRUE.
         FIXTHR    = .TRUE.
         SPEXP     = .FALSE.
      GO TO 100  
 11   CONTINUE
C &&& UNCONT: use internally decontracted basis sets
         CONTRACTED= .FALSE.
      GO TO 100  
 12   CONTINUE
C &&& BRFREE: Branch-free FMM
         BRFREE    = .TRUE.
         CONTRACTED= .FALSE.
      GO TO 100  
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FCK3INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 CALL QUIT('Illegal keyword in FCK3INP.')
            END IF
      ENDIF
 300  CONTINUE

      CALL HEADER('Settings for FCK3 calculation:',0)

      IF (IPRINT .NE. IPRFCK3) THEN
          WRITE (LUPRI,'(A,I5)') ' Print level in FCK3:',IPRFCK3
      END IF

      IF ( FIXTHR ) THEN
        WRITE (LUPRI,'(A,I5)') 'Threshold is fixed.'
      ELSE
        WRITE (LUPRI,'(A,I5)') 'Threshold set dynamically.'
      END IF

      IF ( NOSKIPCL ) THEN
        WRITE (LUPRI,'(A,I5)') 
     &  'No skipping of the nonclassical integrals (only set to 0).'
      END IF

      RETURN
      END

C
C  Sets the BFCFMM logical variable from cbieri common block 
C  for the use in ERIGAM subroutine 
C

      SUBROUTINE SETBFCFMM( BRFREE )
#include "cbieri.h" 
      LOGICAL BRFREE

      BFCFMM = BRFREE
  
      RETURN
      END  

      SUBROUTINE CHECKR2MIN( R2MIN, RODP, ICNTC,ICNTD, RCD )
#include"implicit.h"
#include"priunit.h"
#include"maxaqn.h"
#include"maxorb.h"
#include"mxcent.h"
#include"aobtch.h"
#include"nuclei.h"

C
C     This routine gives the minimal (square of) distance 
C    from point P to the the poing Q which we know is 
C    restricted to be on the line between points C and D. 
C

      DIMENSION RODP(4)

      IF ( ICNTC .EQ. ICNTD ) THEN
        R2MIN = 0.0D0 
        RETURN
      END IF

      PX = RODP(1)
      PY = RODP(2)
      PZ = RODP(3)

      CX = CORD(1,ICNTC)
      CY = CORD(2,ICNTC)
      CZ = CORD(3,ICNTC)

      DX = CORD(1,ICNTD)
      DY = CORD(2,ICNTD)
      DZ = CORD(3,ICNTD)

      RCD2 = RCD*RCD
      RPC2 = (PX-CX)*(PX-CX)+(PY-CY)*(PY-CY)+(PZ-CZ)*(PZ-CZ)
      RPD2 = (PX-DX)*(PX-DX)+(PY-DY)*(PY-DY)+(PZ-DZ)*(PZ-DZ)

C
C      First case is when the angle either angle PCD or PDC is 
C    bigger than 90 degrees ( in that case PQ > min(PC,PD) )
C      If both these angles are less than 90 degrees, than
C    we can actually calculate the minimal RPQ2 after
C    a few line derivation ...
C


      IF ( ABS( RPD2 - RPC2 ) .GT. RCD2 ) THEN
        R2MIN = MIN( RPC2, RPD2 )
      ELSE
        A     = RPC2 + RPD2 - RCD2
        R2MIN = ( 4*RPC2*RPD2 - A*A ) / (4*RCD2)
      END IF

      RETURN
      END 

      SUBROUTINE CHECKRMIN( RMIN, RODP, ICNTC,ICNTD, RCD )
#include"implicit.h"
#include"priunit.h"
#include"maxaqn.h"
#include"maxorb.h"
#include"mxcent.h"
#include"aobtch.h"
#include"nuclei.h"

      DIMENSION RODP(4)

C      AX = CORXBT(IA)
C      AY = CORYBT(IA)
C      AZ = CORZBT(IA)
C
C      BX = CORXBT(IB)
C      BY = CORYBT(IB)
C      BZ = CORZBT(IB)

      PX = RODP(1)
      PY = RODP(2)
      PZ = RODP(3)

      CX = CORD(1,ICNTC)
      CY = CORD(2,ICNTC)
      CZ = CORD(3,ICNTC)

      DX = CORD(1,ICNTD)
      DY = CORD(2,ICNTD)
      DZ = CORD(3,ICNTD)

      RPC = SQRT((PX-CX)**2 + (PY-CY)**2 + (PZ-CZ)**2)
      RPD = SQRT((PX-DX)**2 + (PY-DY)**2 + (PZ-DZ)**2)

      R = RPC + RPD - RCD
     
      IF ( R .GT. 0.0D0 ) THEN
         RMIN = 0.5D0*R
      ELSE
         RMIN = 0.0D0
      END IF
      RETURN
      END 

