C 
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#define MOD_DNF
      SUBROUTINE FCK3DRV( FMAT, DMAT, NDMT, GAB, GABFIT, CFIT, 
     &                    WORK, LWORK )
#include "implicit.h"
#include "priunit.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbihr2.h"
#include "cbieri.h"
#include "cbione.h"
#include "cbifck3.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3rel.h"
#include "fck3thr.h"
#include "fck3mem.h"
#include "erithr.h"
C
#include "dummy.h"
#ifdef MOD_PBC
#include "pbc.h"
#endif
C
#include "nuclei.h"
#if defined VAR_MPI
#include "infpar.h"
#endif
      DIMENSION TIMFCK3(5),FMAT(*), DMAT(*), WORK(*), GAB(*), 
     &          GABFIT(*), CFIT(*)

#if defined VAR_MPI
      CALL MPI3INI(DMAT, NDMT, GAB, GABFIT, CFIT, WORK, LWORK)
#endif
c     dofck3rel=.true.
      CALL INIFCK3REL()

      WORD2MB = 8.0D-6
      STATISTICS(1) = 0.0D0
      STATISTICS(2) = 0.0D0
      STATISTICS(3) = 0.0D0
      STATISTICS(4) = 0.0D0
      STATISTICS(5) = 0.0D0

      IF ( FIXTHR ) THEN
        THRSH = THRFCK3
      ELSE
c       The number of terms screened add up. For valinomycine (C54H90N6O18) 
c       this causes an error 2000 greater than the screening treshold (when
c       the treshold is 1E-7). For up to at least 10 000 atoms these 
c       tresholds should be fine if the molecules have the same (or smaller) 
c       geometrically distributed atomic density.
c        THRSH = 10D0**(-IFTHRS-6)
#if defined PRG_DIRAC
CCaspg, 2006-09-05
CC IFHTRS will dynamically change in dalton but not in dirac (it's not even used...)
CC so until dynamical thresholds are implemented in the scf, we replace the
CC dalton code, 
C          THRSH = SCRFCK*10.0d0**(-6)
CC with just the default/user-defined threshold... 
CC
          THRSH = SCRFCK
CCaspg
#else
          THRSH = 10D0**(-IFTHRS-6)
#endif
      END IF
Caspg
      write (LUPRI,*)
     & 'debug: thresholds in fck3drv: scrfck',SCRFCK,'thrsh',THRSH
Caspg

Careful / ErikT
C      if (dopbc .and. pbcactive) then
C        Allow only very strict  screening for PBC, at least
C        during the test phase. / ErikT
C         THRSH = 1.0D-18
C      end if
Careful
CAREFUL
c     THRSH = 1.0D-18
CAREFUL

      THRPOD   = THRSH
      THRCOD   = THRSH
      THRNCL   = THRSH
      THRCS    = THRSH

      IF ( NOSCRPRIM ) THRPOD = 0.0D0

        NCLFCK3 = .FALSE.
        THRCLS  = 0.0d0


      IF ( .NOT. NCLFCK3 ) THEN 
        NOSKIPCL = .TRUE.
      END IF

c     THRSH screens out induvidual integrals without any significant
c     effect on timings. Hardcoded to zero and reset afterwards.
      THRSHSAVE = THRSH
      THRSH  = 0.0D0

      FACCLS=ERFCIV(THRCLS)


      IF ( NDMT .NE. 1 ) THEN 
        CALL QUIT('In FCK3 internally is not implemented NDMT > 1') 
      END IF

      NDMAT = NDMT

      KFREE = 1 
      LFREE = LWORK

      CALL TIMER('START ',TIMSTR,TIMEND)

      CALL GETTIM(TIMFCK3(1),DUMTIM)

      CALL SETAOBTCH(SPEXP,CONTRACTED,WORK,LWORK,IPRFCK3)

C     Determines the two first dimensions of the density matrix depending 
C     on whether we are using denity-fitting or not.
      CALL DMATDIM

      CALL GETTIM(TIMFCK3(2),DUMTIM)

C     PBC comment: We ultimately want to compute contributions to
C
C      J_{ab}^{0l} = \sum_{mj} \sum_{cd} (ab|cd) D_{cd}^{0j}
C
C      for fixed vectors m and j (the sum is performed by calling this
C      routine several times with different m,j).
C      Here we call SETODBTCH to set up the cd overlap batches.
C      We need to add different translations (this info is passed
C      via variables in the common block include/pbc.h) to ab and
C      cd centers.
C      Translations for a and b: 0 and l, respectively.
C      Translations for c and d: m and m+j, respectively.
C      The last argument to SETODBTCH is a flag indicating
C      whether we have ab (flag=true) or cd (flag=false) batches. /ErikT
      CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &               KFREE,LFREE,WORK,LWORK,GAB,GABFIT,
     &               DMAT,CFIT,NDMT,IPRFCK3,.false.)

      MEMODB = LWORK - LFREE

c      CALL TIMER('FCKINI',TIMSTR,TIMEND)


c      CALL TIMER('FCKSMM',TIMSTR,TIMEND)

      CALL GETTIM(TIMFCK3(3),DUMTIM)

      NPARFTUV_NPAR = 1
     
      NPARFTUV_NSTEPODCLS = NPARFTUV_NPAR

      DO NPARFTUV_KSTARTODCLS = 1, NPARFTUV_NPAR
       IF ( CHCKNAI ) THEN
        IELECTR = 1
          CALL SETFTUV( KIFTUVCLS,KRFTUVCLS,KIFTUV,KRFTUV,KFTUV,
     &         KRFBTCH,KIFBTCH,WORK(KRODBCLS),WORK(KIODBCLS),
     &         WORK(KRODBTCH),WORK(KIODBTCH),DMAT,ISYMDM,NDMT,
     &         KFREE,LFREE,IELECTR,WORK,LWORK,
     &         NPARFTUV_KSTARTODCLS,NPARFTUV_NSTEPODCLS,IPRFCK3)
       END IF

       IELECTR = 2

       NPFTUVSTR = NPARFTUV_KSTARTODCLS
         CALL SETFTUV( KIFTUVCLS,KRFTUVCLS,KIFTUV,KRFTUV,KFTUV,
     &        KRFBTCH,KIFBTCH,WORK(KRODBCLS),WORK(KIODBCLS),
     &        WORK(KRODBTCH),WORK(KIODBTCH),DMAT,ISYMDM,NDMT,
     &        KFREE,LFREE,IELECTR,WORK,LWORK,
     &        NPFTUVSTR,NPARFTUV_NSTEPODCLS,IPRFCK3)

#ifdef MOD_PBC
       if (dopbc .and. pbcactive) then
C        PBC comment: See PBC comment above for context.
C         Here we call SETODBTCH with the last argument set to true,
C         indicating that we want ab batches.
C         Note that we overwrite the cd batches that were constructed
C         in the previous call to SETODBTCH, but this is fine since
C         this data has already been processed and the results stored
C         in variables with names that fit the regexp *FTUV*. /ErikT
          CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &      KFREE,LFREE,WORK,LWORK,GAB,GABFIT,
     &      DMAT,CFIT,NDMT,IPRFCK3,.true.)
       end if
#endif /* MOD_PBC */

       CALL GETTIM(TIMFCK3(4),DUMTIM)

c       CALL TIMER('FCKIN2',TIMSTR,TIMEND)

Careful / ErikT
C       write(LUPRI,*) 'PBC debug: calling GETFCK3',NPARFTUV_KSTARTODCLS
C        WRITE (LUPRI,'(/,2(3(5X,A,I9,/),/))') 
C     &         ' Number of AO   batches : ', NAOBTCH,
C     &         ' Number of OD   batches : ', NODBTCH,
C     &         ' Number of FTUV batches : ', NFTUVBTCH,
C     &         ' Number of AO   classes : ', NAOBCLS,
C     &         ' Number of OD   classes : ', NODBCLS,
C     &         ' Number of FTUV classes : ', NFTUVCLS
C        call flshfo(LUPRI)
Careful       

       CALL GETFCK3(FMAT,NDMT,WORK(KIODBCLS),WORK(KRODBCLS),
     &             WORK(KIODBTCH),
     &             WORK(KRODBTCH),WORK(KIFTUVCLS),WORK(KRFTUVCLS),
     &             WORK(KIFTUV),WORK(KRFTUV),WORK(KFTUV),WORK(KRFBTCH),
     &             WORK(KIFBTCH),WORK(KFREE),LFREE,IPRFCK3)

c       CALL TIMER('FCKGFC',TIMSTR,TIMEND)

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
     &         ' - skip by classes        :', A_SKIPPED2*100, ' % '

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

      CALL FLSHFO(LUPRI)


      THRSH = THRSHSAVE

      dofck3rel=.false.
      call inifck3rel()
      RETURN 
      END 

      SUBROUTINE GETFCK3(FMAT,NDMAT,IODBCLS,RODBCLS,IODBTCH,RODBTCH,
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
      DIMENSION IODBCLS(5,*),RODBCLS(2,*),IODBTCH(2,*),RODBTCH(7,*),
     &          IFTUVCLS(10,*),RFTUVCLS(6,*),IFTUV(*),RFTUV(*),
     &          FTUV(*),RFBTCH(6,*),IFBTCH(2,*),WORK(*), FMAT(*)
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

      CALL GETFCK31( FMAT,NDMAT,IODBCLS,RODBCLS,IODBTCH,RODBTCH,
     &               IFTUVCLS,RFTUVCLS,IFTUV,RFTUV,FTUV,RFBTCH,
     &               IFBTCH,WORK(KINDHR),WORK(KIODHR),
     &               WORK(KLMNVL),WORK(KLMNSM),WORK(KCSQ),
     &               WORK(KLAST),LWRK,IPRINT)


      RETURN
      END

      SUBROUTINE GETFCK31(FMAT,NDMAT,IODBCLS,RODBCLS,IODBTCH,RODBTCH,
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
#include "fck3rel.h"
#include "fck3thr.h"
#include "fck3mem.h"
#include "cbifck3.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C aspg, 2006-06-27 making ftnchek happy
C
#include "erithr.h"
C
#ifdef MOD_PBC
#include "pbc.h"
#endif
      LOGICAL AUX, SCREENCS, ACTIVE
      DIMENSION IODBCLS(5,*), RODBCLS(2,*), IODBTCH(2,*), RODBTCH(7,*),
     &          IFTUVCLS(10,*),RFTUVCLS(6,*),IFTUV(*),RFTUV(*),
     &          RFBTCH(6,*),IFBTCH(2,*),FTUV(*),CSQ(NCSQ1,NCSQ2),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          IODDHR(NRTOP,0:7),LMNPWR(MAXODKCK,MAXODNHK,3),
     &          LMNSYM(MAXODKCK,MAXODNHK,2,0:7),WORK(*), 
     &          KCCF(2,2), FMAT(*)
      dimension tvec1(3), tvec2(3), coord_nf_lat(3)
#if defined VAR_MPI
#include "infpar.h"
      NOFF  = MYNUM
      NSTEP = NODTOT + 1
      CALL GETTIM(TSFCK31,DUMTIM)
C     For the special case of parallel code run in serial NODTOT=-1
      IF (NSTEP.EQ.0) NSTEP = 1
#else
      PARAMETER(NSTEP = 1)
      NOFF  = 0
#endif

      DOSCREENCS = .FALSE.

CMAW  We need the next line if FCK3 is to do the fully classical
CMAW  interactions inside RPQMIN in the Branch-Free scheme:
CMAW      IF ( BRFREE ) NCLFCK3 = .FALSE. 

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

#ifdef MOD_PBC
      if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
C        R tensors referring to different replicas of dummy overlap
C        centers in the near-field should be added up first, and then
C        contracted. This means that we do the loop over near-field
C        translations in this routine instead of doing it in
C        the routine PBC_JOP_REALMAT in PBCDRV.F.
         call pbc_get_nfsize(nf_siz1,nf_siz2,nf_siz3)

C        Allocate memory for the summed R tensor
         khrint1_sum = KLAST
         KLAST = KLAST + MAXPRF * MAXNTUV
         if (KLAST .gt. LWORK) then
            write(LUPRI,*) 'Out memory in GETFCK31.'
            write(LUPRI,*) '  Needed    (Mwords): ',MAXPRF * MAXNTUV
            write(LUPRI,*) '  Available (Mwords):',LWORK-khrint1_sum
            call quit('Out memory in GETFCK31.')
         end if
         khrint1_aux = khrint1_sum
      else
#endif /* MOD_PBC */
C        R tensors will not be added up.
C        The loop over near-field translations occurs in PBC_JOP_REALMAT
C        in PBCDRV.F and not here.
         nf_siz1 = 0
         nf_siz2 = 0
         nf_siz3 = 0
         khrint1_aux = KHRINT1
#ifdef MOD_PBC
      end if
#endif /* MOD_PBC */

      MEMFCK31 = KRJ000 - 1 
      MEMFCK32 = KLAST  - KRJ000  

      IF (KLAST .GT. LWORK) CALL STOPIT('FCK3','FCK31',KLAST,LWORK)     
      LWRK   = LWORK - KLAST + 1

      IELECTR  = 1

      IABCLSSTART = 1
      IABCLSEND = NODCLSNORM
      AUX = .FALSE.

c     DO IABCLS = IABCLSSTART, IABCLSEND
      DO IABCLS = IABCLSSTART + NOFF, IABCLSEND, NSTEP
       ICENTA    = IODBCLS(1,IABCLS)
       ICENTB    = IODBCLS(2,IABCLS)
       NHKTAB    = IODBCLS(3,IABCLS)
       KABSTBTCH = IODBCLS(4,IABCLS)
       NABBTCH   = IODBCLS(5,IABCLS)

       CSABCLS   = RODBCLS(1,IABCLS)

       JMAXAB = NHKTAB - 1
       NHKABR = NHKTAB + NREL 
       NTUVAB = NHKABR*(NHKABR+1)*(NHKABR+2)/6

#ifdef MOD_PBC
       if (dopbc .and. pbcactive) then
C         We are ultimately interested in
C         J_{ab}^{0l} = \sum_{mj} \sum_{cd} (ab|cd) D_{cd}^{0j}.
C         Here we are dealing with ab overlaps and need to put in
C         the translation vector l to get the PBC case right.
          tvec1(1) = 0.0
          tvec1(2) = 0.0
          tvec1(3) = 0.0
          tvec2(1) = rspc_vec_std(1)
          tvec2(2) = rspc_vec_std(2)
          tvec2(3) = rspc_vec_std(3)
          CALL GETR12( NSYMOPR, NSYMOPS, ICENTA, ICENTB, 1,
     &                 tvec1, tvec2 )
       else
#endif /* MOD_PBC */
          tvec1(1) = 0.0
          tvec1(2) = 0.0
          tvec1(3) = 0.0
          tvec2(1) = 0.0
          tvec2(2) = 0.0
          tvec2(3) = 0.0          
          CALL GETR12( NSYMOPR, NSYMOPS, ICENTA, ICENTB, 1,
     &                 tvec1, tvec2 )
#ifdef MOD_PBC
       end if
#endif /* MOD_PBC */

       DO IABBTCH = KABSTBTCH, ( KABSTBTCH + NABBTCH - 1 ) 

        IAOBA = IODBTCH(1,IABBTCH)
        IAOBB = IODBTCH(2,IABBTCH)

c       IF (.TRUE.) THEN
       ACTIVE = ACTVAOBTCH(IAOBA).AND.ACTVAOBTCH(IAOBB)

        CSABBTCH = RODBTCH(6,IABBTCH)

        CALL GETAOB(NPRFA,NCTFA,NTYPA,NHKTA,JMAXA,KHKTA,KCKTA,
     &              KCSQA,KEXPA,KCCF(1,1),IAOBA)

        CALL GETAOB(NPRFB,NCTFB,NTYPB,NHKTB,JMAXB,KHKTB,KCKTB,
     &              KCSQB,KEXPB,KCCF(1,2),IAOBB)

C       PBC comment: SCREENPRIM12 uses coordinates and a distance set
C        by GETR12, so the necessary translations have already
C        been put in. Hence, SCREENPRIM12 works without modification for
C        PBC. One thing, though: In the molecular case, one screens out
C        just a few terms, but in the PBC case screening means neglecting
C        a converging infinite series. This means that the error
C        estimates may be a little different. / ErikT
        CALL SCREENPRIM12(NPRFAB,RODBTCH(1,IABBTCH),WORK(KIPRAB),.TRUE.,
     &                    NPRFA,NCTFA,NTYPA,NPRFB,NCTFB,NTYPB,
     &                    EXPBTCH(KEXPA),EXPBTCH(KEXPB),
     &                    CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &                    CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &                    R122(0,0),COOR1(1,0),COOR2(1,0), 1,
     &                    IPRINT)
        EXTPMAX = RODBTCH( 4, IABBTCH )
        SCRAB   = RODBTCH( 5, IABBTCH )

        CALL SETETIJ(WORK(KETIJ),WORK(KPREFAB),
     &               R12(1,0,0),R122(0,0),
     &               WORK(KIPRAB),EXPBTCH(KEXPA),EXPBTCH(KEXPB),
     &               WORK(KCOORAB),WORK(KEXPAB),WORK(KHPI),WORK(KPA),
     &               WORK(KPB),NPRFA,NHKTA,JMAXA,
     &               NPRFB,NHKTB,JMAXB,NPRFAB,NHKTAB,IELECTR,IPRINT)
        CALL MXETIJ(WORK(KETIJ),NPRFAB,JMAXA+NDER,JMAXB+NDER,ETIJAB)
        CALL SETP(WORK(KDATAP),WORK(KIPRAB),
     &            COOR1(1,0),COOR2(1,0),EXPBTCH(KEXPA),
     &            EXPBTCH(KEXPB),NPRFA,NPRFB,NPRFAB,IPRINT)

        CALL DZERO(WORK(KHRINT),NPRFAB*NTUVAB*NDMAT)

c
c       find largest exponent of the ab distribution
        CALL GETMAXEXPAB(EXPMAB,EXPBTCH(KEXPA),EXPBTCH(KEXPB),
     &                   NPRFA,NPRFB)

        IFCLSSTART = 1
        IFCLSEND   = NFTUVCLS


        DO IFCLS = IFCLSSTART, IFCLSEND
 
         ICENTC      = IFTUVCLS(  9, IFCLS )
         ICENTD      = IFTUVCLS( 10, IFCLS ) 

         RCD         = RFTUVCLS(  1, IFCLS ) 
         EXTQMAX     = RFTUVCLS(  2, IFCLS ) 
         FTUVMX      = RFTUVCLS(  3, IFCLS )
         SCRCD       = RFTUVCLS(  4, IFCLS ) 
         EXPMCD      = RFTUVCLS(  5, IFCLS ) 

         CSFTUVCLS   = RFTUVCLS(  6, IFCLS ) 
         RUBAB = R1UPPERB(JMAXA+JMAXB+NREL,EXPMAB)
         RUBCD = R1UPPERB(IFTUVCLS(6,IFCLS),EXPMCD)
         IF ( NCLFCK3 .OR. BRFREE ) 
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

CMAW  We need the full test below if FCK3 is to do the fully classical
CMAW  interactions inside RPQMIN in the Branch-Free scheme:
C         IF ( ( ( NOSKIPCL .OR. (  
C     &    ( .NOT. NCLFCK3 ).OR.( RMIN .LE. EXTPMAX + EXTQMAX ) ) .AND.
C     &    ( .NOT. BRFREE ).OR.( ( RMIN .LE. RPQMIN ) .OR.  
C     &                  ( RMIN .LE. EXTPMAX + EXTQMAX )) )  ) .AND.
C     &        ( NOSKIPCS .OR. ( SCRAB * SCRCD .GT. THRCOD ) ) ) THEN

         TOTFAC = SCRAB * ETIJAB * RUBAB * RUBCD * FTUVMX

         SCREENCS = DOSCREENCS .AND. CSFTUVCLS*CSABBTCH .LT. THRCS
     &                   .AND. (.NOT. NCSTUV)

C        Screening test
         IF ( ( NOSKIPCL .OR. ( RMIN .LE. ( EXTPMAX+EXTQMAX ) ) ) .AND.
     &        ( NOSKIPCS .OR. ( TOTFAC .GT. THRCOD ) ) .AND. 
     &        ( .NOT. SCREENCS ) ) THEN

          IF (IPRINT.GT.25) THEN
            write(LUPRI,*) 'NPRFAB, NTUVAB ',NPRFAB, NTUVAB
            CALL HEADER('HRINT before screen',-1)
            CALL OUTPUT(WORK(KHRINT),1,NPRFAB,1,NTUVAB,NPRFAB,NTUVAB,1,
     &                  LUPRI)
          ENDIF

          NTFTUVBTCH  = IFTUVCLS(  1, IFCLS )
          KSIFTUV     = IFTUVCLS(  3, IFCLS )
          KSRFTUV     = IFTUVCLS(  4, IFCLS )
          KSFTUV      = IFTUVCLS(  5, IFCLS )
          JMAXF       = IFTUVCLS(  6, IFCLS )
          KSFTUVBTCH  = IFTUVCLS(  7, IFCLS )
          NTYPF       = IFTUVCLS(  8, IFCLS )

          JMAX = JMAXAB + JMAXF + NREL

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

          IF ( NPRF .GT. 0 ) THEN 

#ifdef MOD_PBC
          if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
             call dzero(WORK(khrint1_sum), MAXPRF * NTUV)
          end if
#endif /* MOD_PBC */

          do mx = -nf_siz1, nf_siz1  ! PBC near-field loop (dim. 1)
          do my = -nf_siz2, nf_siz2  ! PBC near-field loop (dim. 2)
          do mz = -nf_siz3, nf_siz3  ! PBC near-field loop (dim. 3)

#ifdef MOD_PBC
          if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
C            See above for context.
C            Construct near-field translation vector, put it in
C            common block variables from pbc.h (these are checked
C            by the routine SETPQ below).
             nf_vec_index(1) = mx
             nf_vec_index(2) = my
             nf_vec_index(3) = mz
             coord_nf_lat(1) = dfloat(mx)
             coord_nf_lat(2) = dfloat(my)
             coord_nf_lat(3) = dfloat(mz)
             call pbc_lat2std_coord(coord_nf_lat, nf_vec_std)
          end if
#endif /* MOD_PBC */
       
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
c           IF(NPRF.GT.0)
            CALL SETPQ(WORK(KFACINT), WORK(KALPHA), WORK(KCOORPQ), NPRF,
     &               IPQXYZ, WORK(KDATAP),RFTUV(KSRFTUV), 
     &               NPRFAB, NTFTUVPRIM, JMAX, NTYP, NTUV )    
            STATISTICS(2) = STATISTICS(2) + NPRF
           END IF
          END IF          

            CALL ERIGAM(WORK(KRJ000),WORK(KFACINT),WORK(KALPHA),
     &           WORK(KCOORPQ),WORK(KLAST),LWRK,IPRINT)

            CALL ERITUV(WORK(KHRINT1),WORK(KRJ000),
     &           WORK(KCOORPQ),INDHER,IODDHR(1,IPQXYZ),
     &           WORK(KLAST),LWRK,IPRINT)

#ifdef MOD_PBC
          if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
C            Now we have the R tensor for the current near-field
C            translation. Add it to the total R tensor.
             do ir = 0, MAXPRF * NTUV - 1
                WORK(khrint1_sum+ir) = WORK(khrint1_sum+ir)
     &                               + WORK(KHRINT1+ir)
             end do             
          end if
#endif /* MOD_PBC */

          end do            ! PBC end near-field loop (dim. 3)
          end do            ! PBC end near-field loop (dim. 2)
          end do            ! PBC end near-field loop (dim. 1)


          JMXABR = JMAXAB + NREL
          IF ( NCLFCK3 .AND. ( NPRF .NE. NPRF0 ) ) THEN 
           IF ( PACKFTUV ) THEN
            CALL HR2FTUVNCL1(WORK(KHRINT),WORK(khrint1_aux),
     &                  WORK(KFTUV2),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMXABR,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,
     &                  NPRFAB,NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP)
           ELSE
            CALL HR2FTUVNCL2(WORK(KHRINT),WORK(khrint1_aux),
     &                  WORK(KFTUV2),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMXABR,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,
     &                  NPRFAB,NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP)
           END IF
          ELSE
           IF ( PACKFTUV ) THEN
            CALL HR2FTUV1(WORK(KHRINT),WORK(khrint1_aux),FTUV(KSFTUV),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMXABR,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,NPRFAB,
     &                  NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP )
           ELSE
            CALL HR2FTUV2(WORK(KHRINT),WORK(khrint1_aux),FTUV(KSFTUV),
     &                  INDHER,IODDHR(1,IPQXYZ),NDMAT,JMXABR,JMAXF,
     &                  NTUVAB,NTUVCD,NTUV,NPRFAB,
     &                  NTFTUVPRIM,NTFTUVPRIM_TOT,NPRF,JTOP )

           END IF
          END IF

          STATISTICS(4) = STATISTICS(4) + (NPRF0-NPRF)
C
C         End of the main time consuming part !!!
C         =======================================

         END IF ! NPRF

         IF (IPRINT.GT.25) THEN
           CALL HEADER('HRINT after screen',-1)
           CALL OUTPUT(WORK(KHRINT),1,NPRFAB,1,NTUVAB,NPRFAB,NTUVAB,1,
     &                 LUPRI)
         ENDIF

         ELSE 
          STATISTICS(5) = STATISTICS(5) + NPRF0
         END IF !  Screening test

        END DO  !IFCLS

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
     &         CSQ,KCSQA,KCSQB,
     &         CCFBTCH(KCCF(ICOMPA,1)),CCFBTCH(KCCF(ICOMPB,2)),
     &         ICOMPA,ICOMPB,NPRFAB,
     &         JA,JB,JMAXA, JMAXB, JMAXAB, NTUVAB,
     &         KHKTA, KHKTB, KCKTA, KCKTB,
     &         KCKMAX,NHKMAX,IAOBA,IAOBB,IPRINT )

         END DO
        END DO
       END DO !AB-batches
       
      END DO ! AB-calsses

      IF (IPRINT .GT. 20) THEN
        CALL HEADER('FMAT in FCK3DRV',-1)
        CALL OUTPUT(FMAT,1,NDIMA,1,NDIMB,NDIMA,NDIMB,1,LUPRI)
      ENDIF

#if defined VAR_MPI
      CALL GETTIM(TEFCK31,DUMTIM)
      CALL FCK3COLLECT(FMAT,NDMAT,TEFCK31-TSFCK31,WORK,LWORK,IPRINT)
#endif
      RETURN 
      END

C
C   --------------------------------------------------------------
C

      SUBROUTINE SETP( DATAP, IPRAB, COORA, COORB,EXPBTA, EXPBTB, 
     &                 NPRFA, NPRFB, NPRFAB, IPRINT )
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
C
#ifdef MOD_PBC
#include "pbc.h"
#endif
C
      DIMENSION DATAP(6,NPRFAB), IPRAB(2,NPRFAB), COORA(3),COORB(3),
     &          EXPBTA(NPRFA), EXPBTB(NPRFB)
      dimension pos_ab(3), pos_nf_lat(3), pos_nf_std(3)
      integer nf_pos_index(3)

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

#ifdef MOD_PBC
       if (dopbc .and. pbcactive) then
C        We now have the position P of the ab overlap and the
C        position Q of the cd overlap. However, we want the
C        electrostatic near-field to be centered in the cell nf_pos
C        containing P, so that the m (center of cell containing the
C        c orbital) and m' (center of cell containing the d orbital)
C        vectors refer to the neighbours of nf_pos. Hence, our Q is
C        really a vector relative to nf_pos, not relative to the
C        origin (center of the central cell). We therefore account
C        for nf_pos here. To save some CPU cycles and not mess up the
C        ab-cd screening, we also subtract nf_pos from P instead of
C        adding it to Q.
C
C        (Additional comment: Note that nf_pos depends on the a and b
C        orbitals. Since we treat many such orbitals at once it is
C        difficult to put in the vector nf_pos already in the m
C        vector. This is why we do it at such a late stage.) / ErikT
         pos_ab(1) = DATAP(4,I)
         pos_ab(2) = DATAP(5,I)
         pos_ab(3) = DATAP(6,I)
         call pbc_nf_center(pos_ab, nf_pos_index)
         pos_nf_lat(1) = dfloat(nf_pos_index(1))
         pos_nf_lat(2) = dfloat(nf_pos_index(2))
         pos_nf_lat(3) = dfloat(nf_pos_index(3))
         call pbc_lat2std_coord(pos_nf_lat, pos_nf_std)
         DATAP(4,I) = DATAP(4,I) - pos_nf_std(1)
         DATAP(5,I) = DATAP(5,I) - pos_nf_std(2)
         DATAP(6,I) = DATAP(6,I) - pos_nf_std(3)
       end if
#endif /* MOD_PBC */

      END DO

      RETURN
      END 

C
C   --------------------------------------------------------------
C

      SUBROUTINE SETPQ( FACINT, ALPHA, COORPQ, NPRF, IPQXYZ,
     &                  DATAP, RFTUV, NPRFAB, NPRFCD, JMAX, NTYP, NTUV )
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
#ifdef MOD_PBC
#include "pbc.h"
#endif
C
      DIMENSION DATAP(6,*),RFTUV(8,*),COORPQ(NPRF,3),
     &          FACINT(*),ALPHA(*),IPQ0(3)
      INTEGER X

      I=0 

      DO IPRFAB=1,NPRFAB

       EXPP  = DATAP(1,IPRFAB)
       EXPPI = DATAP(2,IPRFAB)
       FACPI = DATAP(3,IPRFAB)
       PX    = DATAP(4,IPRFAB)
       PY    = DATAP(5,IPRFAB)
       PZ    = DATAP(6,IPRFAB)

#ifdef MOD_PBC
       if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
C         See comments in GETFCK31 (file FCK3DRV.F) and
C         PBC_JOP_REALMAT (file PBCDRV.F) for context.
C         R tensors referring to different replicas should be
C         added up. The dummy overlaps (center Q) were constructed
C         for the central cell and need to have the near-field translation
C         added. As a small efficiency improvement, we subtract from P
C         instead of adding to Q. / ErikT
          PX = PX - nf_vec_std(1)
          PY = PY - nf_vec_std(2)
          PZ = PZ - nf_vec_std(3)
       end if
#endif /* MOD_PBC */


       DO IPRFCD=1,NPRFCD
        I = I + 1

        EXPQ    = RFTUV(1,IPRFCD)
        EXPQI   = RFTUV(2,IPRFCD)
        QX      = RFTUV(3,IPRFCD)
        QY      = RFTUV(4,IPRFCD)
        QZ      = RFTUV(5,IPRFCD)
        DWEIGHT = RFTUV(8,IPRFCD)
      
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
c       FACINT(I)= FACPI*EXPQI/SQRT(EXPP+EXPQ)
        FACINT(I)= FACPI*EXPQI/SQRT(EXPP+EXPQ)*DWEIGHT

        COORPQ(I,1)= PX - QX
        COORPQ(I,2)= PY - QY
        COORPQ(I,3)= PZ - QZ

       END DO
      END DO


      IF ( NTYP .GT. 3) THEN 
       IPQXYZ  = 0
       IPQ0(1) = 0
       IPQ0(2) = 0
       IPQ0(3) = 0
#ifdef MOD_PBC
      else if (dopbc .and. pbcactive .and. .not. pbc_sep_Rtensor) then
       IPQXYZ  = 0
       IPQ0(1) = 0
       IPQ0(2) = 0
       IPQ0(3) = 0
#endif /* MOD_PBC */
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
      DIMENSION RODBTCH(7),IFBTCH(2,*),RFBTCH(6,*)

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
#include "densfit.h"
      DIMENSION COORAB(3,NPRFAB),EXPABP(3,NPRFAB),COORPQ(NPRF,3),
     &   RFTUV(8,NPRFCD_TOT),IFTUV(2,NPRFCD_TOT),FACINT(*),ALPHA(*),
     &   RODBTCH(7),IFBTCH(2,*),RFBTCH(6,*),IPQ0(3),
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
#include "densfit.h"
      DIMENSION COORAB(3,NPRFAB),EXPABP(3,NPRFAB),COORPQ(NPRF,3),
     &   RFTUV(8,NPRFCD_TOT),IFTUV(2,NPRFCD_TOT),FACINT(*),ALPHA(*),
     &   RODBTCH(7),IFBTCH(2,*),RFBTCH(6,*),IPQ0(3),
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

      IF ( NTYP .GT. 3) THEN 
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
#include "implicit.h"
#include "hertop.h"

      JTOP  = JTOPVAL
      NRTOP = NRTOPVAL

      RETURN
      END


      SUBROUTINE PUT2ERICOM2( NPPXVAL, IPQ0VAL, IPQXYZVAL,
     &                        JMAXVAL, NTUVVAL )
#include "implicit.h"
#include "ericom.h"
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
C    INPUT: HRINT1, FTUV
C    OUTPUT: HRINT


      SUBROUTINE HR2FTUV1( HRINT, HRINT1, FTUV,  INDHER, IODDHR, NDMAT,
     &                    JMAXAB, JMAXCD, NTUVAB, NTUVCD, NTUV, 
     &                    NPRFAB, NPRFCD, NPRFCD_TOT,NPRF, JTOP )
#include "implicit.h"
#include "priunit.h"

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
C     Driver for adding contribution to the Fock-matrix
C
      SUBROUTINE ADD2FCK3(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQ,KCSQA,KCSQB,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &        JA, JB, JMAXA, JMAXB, JMAXAB, NTUVAB,
     &        KHKTA, KHKTB, KCKTA, KCKTB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,IPRINT)

C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3rel.h"
#include "densfit.h"
C
#include "pbc.h"
C
      DIMENSION FMAT(*),FMATCP(*),FMATSP(*),FTMP(*),ETIJ(*),
     &   PREFAB(*),HRINT(*),IPRAB(*),INDHER(*),IODDHR(*),LMNPWR(*),
     &   CSQ(NCSQ1,NCSQ2)

      IF (NDERLHS.EQ.0) THEN
C       Regular case
        CALL ADD2FCK31(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQ(KCSQA,1),CSQ(KCSQB,1),CCFA,CCFB,ICOMPA,ICOMPB,
     &        NPRFAB,JA,JB,JMAXA,JMAXB,JMAXAB,NTUVAB,
     &        KHKTA, KHKTB, KCKTA, KCKTB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,IPRINT )
      ELSEIF( NDERLHS.EQ.1) THEN
C       Relativistic case
        CALL ADD2FCK3REL(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQ,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &        JA, JB, JMAXA+NDERLHS, JMAXB+NDERLHS, JMAXAB, NTUVAB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,IPRINT)
      ELSE
        WRITE(LUPRI,'(1X,A,I2,A)') 'NDERLHS =',NDERLHS,
     &     ' not a case in ADD2FCK3!'
        CALL QUIT('ADD2FCK3 case not permitted')
      ENDIF

      RETURN
      END

C
C
C
      SUBROUTINE ADD2FCK3REL(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQ,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &        JA, JB, JMAXA, JMAXB, JMAXAB, NTUVAB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,IPRINT)
#include "implicit.h"
      DIMENSION FMAT(*),FMATCP(*),FMATSP(*),FTMP(*),ETIJ(*),
     &          PREFAB(*),HRINT(*),IPRAB(*),INDHER(*),IODDHR(*),
     &          LMNPWR(*),CSQ(*),CCFA(*),CCFB(*)
C
C     Large components
      CALL SETJKKK(JSA,IDUM,KCKTA,KHKTA,JA+1) 
      CALL SETJKKK(JSB,IDUM,KCKTB,KHKTB,JB+1) 
      CALL ADD2FCK3R1(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQ,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &        JSA, JSB, JMAXA, JMAXB, JSA+JSB, NTUVAB,
     &        KHKTA, KHKTB, KCKTA, KCKTB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,.FALSE.,IPRINT)
C     Small components regular
c     CALL SETJKKK(JSA,IDUM,KCKTA,KHKTA,JA+2) 
c     CALL SETJKKK(JSB,IDUM,KCKTB,KHKTB,JB+2) 
c     CALL ADD2FCK3R1(FMAT,FMATCP,FMATSP,FTMP,
c    &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
c    &        IODDHR,CSQ,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
c    &        JSA, JSB, JMAXA, JMAXB, JSA+JSB, NTUVAB,
c    &        KHKTA, KHKTB, KCKTA, KCKTB,
cAREFUL
c    &        KCKMAX, NHKMAX,IAOBA,IAOBB,.FALSE.,IPRINT)
c    &        KCKMAX, NHKMAX,IAOBA,IAOBB,.TRUE.,IPRINT)
CAREFUL
C     Small components modified
      IF ((JA.GT.0).AND.(JB.GT.0)) THEN
        CALL SETJKKK(JSA,IDUM,KCKTA,KHKTA,JA) 
        CALL SETJKKK(JSB,IDUM,KCKTB,KHKTB,JB) 
        CALL ADD2FCK3R1(FMAT,FMATCP,FMATSP,FTMP,
     &          NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &          IODDHR,CSQ,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &          JSA, JSB, JMAXA, JMAXB, JSA+JSB, NTUVAB,
     &          KHKTA, KHKTB, KCKTA, KCKTB,
     &          KCKMAX, NHKMAX,IAOBA,IAOBB,.TRUE.,IPRINT)
      ENDIF
C
      RETURN
      END

      SUBROUTINE ADD2FCK3R1(FMAT,FMATCP,FMATSP,FTMP,
     &        NDMAT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &        IODDHR,CSQ,CCFA,CCFB,ICOMPA,ICOMPB,NPRFAB,
     &        JA, JB, JMAXA, JMAXB, JMAXAB, NTUVAB,
     &        KHKTA, KHKTB, KCKTA, KCKTB,
     &        KCKMAX, NHKMAX,IAOBA,IAOBB,MODIFIED,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3rel.h"
#include "densfit.h"
C
#include "pbc.h"
C
      LOGICAL MODIFIED
      DIMENSION FMAT(NDIMA,NDIMB,NDMAT),FMATCP(KCKTA,KCKTB,NPRFAB),
     &   FMATSP(KHKTA,KHKTB,NPRFAB),FTMP(*),ETIJ(*),
     &   PREFAB(NPRFAB),HRINT(NPRFAB,NTUVAB,NDMAT),IPRAB(2,NPRFAB),
     &   INDHER(*),IODDHR(*),LMNPWR(*),CSQ(*),CCFA(*),CCFB(*)

      IF ( IPRINT .GE. 10 ) THEN
        CALL HEADER('Output from ADD2FCK3R1',-1)

        WRITE(LUPRI,*) 'JMAXA, KHKTA, KCKTA, NPRFA, NCTFA',
     &                  JMAXA, KHKTA, KCKTA, NPRFA, NCTFA
        WRITE(LUPRI,*) 'JMAXB, KHKTB, KCKTB, NPRFB, NCTFB',
     &       JMAXB, KHKTB, KCKTB, NPRFB, NCTFB
        WRITE(LUPRI,*) 'JMAXAB, NTUVAB, NDMAT, NBASE ',
     &                  JMAXAB, NTUVAB, NDMAT, NBASE 
      END IF 

      IDMAT = 1 

      IREPA  = 0
      IREPB  = 0

      DO IPRFAB = 1, NPRFAB
C       Contract E-coefficients with hermite-integrals (or extract
C       proper intermediates for integrals over the modified solid 
C       harmonic Gaussians)
        CALL A2F3HER(FMATCP(1,1,IPRFAB),HRINT(1,1,IDMAT),ETIJ,PREFAB,
     &               LMNPWR,INDHER,NPRFAB,IPRFAB,JMAXAB,JA,JB,
     &               KCKTA,KCKTB,IDMAT,NTUVAB,JMAXA,JMAXB,MODIFIED)

C      Shperical transformation
       CALL A2F3SPHR(FMATCP(1,1,IPRFAB),FMATSP(1,1,IPRFAB),CSQ,FTMP,
     &               JA,JB,KHKTA,KCKTA,KHKTB,KCKTB)

C      Print intermediates
       IF ( IPRINT .GE. 15 ) THEN
        CALL A2F3PRINT(FMATCP(1,1,IPRFAB),FMATSP(1,1,IPRFAB),FTMP,
     &                 IPRFAB,KHKTA,KCKTA,KHKTB,KCKTB)
       END IF

C      Contraction (with contraction coefficients) and distribute to
C      Fock matrix
       IF (.NOT.MODIFIED) THEN
       CALL A2F3FOCK(FMAT(1,1,IDMAT),FMATSP(1,1,IPRFAB),CCFA,CCFB,
     &               IPRAB(1,IPRFAB),NCTFA,NCTFB,NPRFA,NPRFB,
     &               KHKTA,KHKTB,IAOBA,IAOBB,ICOMPA,ICOMPB)
       ENDIF

      ENDDO
C     
      RETURN
      END

      SUBROUTINE A2F3FOCK(FMAT,FMATSP,CCFA,CCFB,IPRAB,NCTFA,NCTFB,
     &                    NPRFA,NPRFB,KHKTA,KHKTB,IAOBA,IAOBB,
     &                    ICOMPA,ICOMPB)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "densfit.h"
C
#include "pbc.h"
      DIMENSION FMAT(NDIMA,NDIMB),FMATSP(KHKTA,KHKTB),
     &          CCFA(NPRFA,NCTFA), CCFB(NPRFB,NCTFB), IPRAB(2)

      IPRFA = IPRAB(1)
      IPRFB = IPRAB(2)

      DO ICTFA = 1, NCTFA
       DO ICTFB = 1, NCTFB

         KNDX = KNDXBTCH(ICOMPA,IAOBA,MODA)
         IAOA = INDXBTCH(KNDX - 1 + ICTFA, MODA)
         
         IOFFSET = 0
         KNDX = KNDXBTCH(ICOMPB,IAOBB,MODB)
         IAOB = INDXBTCH(KNDX - 1 + ICTFB, MODB) - IOFFSET
C
C Here it we should probably put summation over 
C  FMATSP(ISPA,ISPB) separately outside 
C 

C         IF ( (IAOA .EQ. IAOB ) .OR. ( IAOBA .EQ. IAOBB )) THEN 
#ifdef MOD_PBC
        if (dopbc .and. pbcactive) then
         DO ISPA=1, KHKTA
          DO ISPB=1, KHKTB
            FMAT(IAOA+ISPA-1,IAOB+ISPB-1) = 
     &           FMAT(IAOA+ISPA-1,IAOB+ISPB-1) 
     &         + FMATSP(ISPA,ISPB) 
     &         * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
          END DO
         END DO
        else IF ( ( IAOBA .EQ. IAOBB )) THEN 
#else
        IF ( ( IAOBA .EQ. IAOBB )) THEN 
#endif /* MOD_PBC */
         DO ISPA=1, KHKTA
          DO ISPB=1, KHKTB
            FMAT(IAOA+ISPA-1,IAOB+ISPB-1) = 
     &           FMAT(IAOA+ISPA-1,IAOB+ISPB-1) 
     &         + FMATSP(ISPA,ISPB) 
     &         * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
          END DO
         END DO
        ELSE
         DO ISPA=1, KHKTA
          DO ISPB=1, KHKTB
             FMAT(IAOA+ISPA-1,IAOB+ISPB-1) = 
     &           FMAT(IAOA+ISPA-1,IAOB+ISPB-1) 
     &         + FMATSP(ISPA,ISPB) 
     &         * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
               FMAT(IAOB+ISPB-1,IAOA+ISPA-1) = 
     &              FMAT(IAOB+ISPB-1,IAOA+ISPA-1) 
     &              + FMATSP(ISPA,ISPB) 
     &              * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
           END DO
         END DO
       END IF

       END DO
      END DO

 
      RETURN
      END 
C
C
C
      SUBROUTINE A2F3PRINT(FMATCP,FMATSP,FTMP,IPRFAB,
     &                     KHKTA,KCKTA,KHKTB,KCKTB)
#include "implicit.h"
#include "priunit.h"
      DIMENSION FMATCP(*),FMATSP(*),FTMP(*)
C
      CALL HEADER('FMATCP',-1)
      WRITE( LUPRI, * ) 'IPRFAB, KCKTA, KCKTB: ',
     &                   IPRFAB, KCKTA, KCKTB   
      CALL OUTPUT(FMATCP,1,KCKTA,1,KCKTB,KCKTA,KCKTB,1,LUPRI)
      CALL HEADER('FTMP',-1)
      CALL OUTPUT(FTMP,1,KCKTA,1,KHKTB,KCKTA,KHKTB,1,LUPRI)
      CALL HEADER('FMATSP',-1)
      WRITE( LUPRI, * ) 'IPRFAB, KHKTA, KHKTB: ',
     &                   IPRFAB, KHKTA, KHKTB   
      CALL OUTPUT(FMATSP,1,KHKTA,1,KHKTB,KHKTA,KHKTB,1,LUPRI)
C
      RETURN
      END
C
C
C
      SUBROUTINE A2F3SPHR(FMATCP,FMATSP,CSQ,FTMP,JA,JB,
     &                    KHKTA,KCKTA,KHKTB,KCKTB)
#include "implicit.h"
#include "fck3odb.h"
      DIMENSION FMATCP(*),FMATSP(*),CSQ(NCSQ1,NCSQ2)

      CALL DGEMM('N','T', KCKTA, KHKTB, KCKTB, 1.0D0,
     &           FMATCP, KCKTA,CSQ(KSQADR(JB,0,0,0,0),1), KHKTB, 
     &           0.0D0, FTMP, KCKTA  )

      CALL DGEMM('N','N', KHKTA, KHKTB, KCKTA, 1.0D0,
     &           CSQ(KSQADR(JA,0,0,0,0),1), KHKTA, FTMP, KCKTA, 
     &           0.0D0, FMATSP, KHKTA  )      
      RETURN
      END

C
      SUBROUTINE A2F3HER(FMATCP,HRINT,ETIJ,PREFAB,LMNPWR,INDHER,
     &                   NPRFAB,IPRFAB,JMAXAB,JA,JB,KCKTA,KCKTB,IDMAT,
     &                   NTUVAB,JMAXA,JMAXB,MODIFIED)
#include "implicit.h"
#include "fck3rel.h"
#include "fck3odb.h"
      INTEGER T,U,V,TUV
      LOGICAL MODIFIED
      DIMENSION FMATCP(KCKTA,KCKTB),
     &   ETIJ(NPRFAB,0:JMAXA+JMAXB,0:JMAXA,0:JMAXB,3,*),
     &   PREFAB(NPRFAB),HRINT(NPRFAB,NTUVAB),
     &   INDHER(0:JTOP,0:JTOP,0:JTOP),
     &   LMNPWR(MAXODKCK,MAXODNHK,3)

C     Modified Hermite Gaussians (extraced from the regular Hermite 
C     Gaussians with angular momentum increased by two)
      IF (MODIFIED) THEN
        CALL A2F3HRX(FMATCP,LMNPWR,INDHER,JA,JB,KCKTA,KCKTB)
C     Regular Hermite Gaussians
      ELSE
        CALL DZERO(FMATCP, KCKTA*KCKTB)
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
     &          ETIJ(IPRFAB,U,MA,MB,2,1)*ETIJ(IPRFAB,V,NA,NB,3,1)

              FMATCP(KCA,KCB) = 
     &             FMATCP(KCA,KCB) + FAC*PREFAB(IPRFAB)
     &                       * HRINT(IPRFAB,TUV)
             END DO
            END DO

          END DO
         END DO
        END DO
      ENDIF
      RETURN
      END
C
      SUBROUTINE A2F3HRX(FMATCP,LMNPWR,INDHER,JA,JB,KCKTA,KCKTB)
#include "implicit.h"
#include "fck3rel.h"
#include "fck3odb.h"
      DIMENSION FMATCP(KCKTA,KCKTB),INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          LMNPWR(MAXODKCK,MAXODNHK,3)
      REAL*8, ALLOCATABLE :: FMATC2P(:,:)
      DIMENSION IA(3),IB(3)
C
      KCKTAP = (JA + 3) * (JA + 4)
      KCKTBP = (JB + 3) * (JB + 4)
      ALLOCATE(FMATC2P(KCKTAP,KCKTBP))
      CALL DCOPY(KCKTAP*KCKTBP,FMATCP,1,FMATCP2P,1)
      CALL DZERO(FMATCP,KCKTA*KCKTB)
      DO IDIRA = 1, 3
        DO IDIRB = 1, 3
          DO KCA = 1, KCKTA
            IA(1) = LMNPWR(KCA,JA+1,1)
            IA(2) = LMNPWR(KCA,JA+1,2)
            IA(3) = LMNPWR(KCA,JA+1,3)
            IA(IDIRA) = IA(IDIRA) + 2
            KCPA = INDHER(IA(1),IA(2),IA(3))
            DO KCB = 1, KCKTB
              IB(1) = LMNPWR(KCB,JB+1,1)
              IB(2) = LMNPWR(KCB,JB+1,2)
              IB(3) = LMNPWR(KCB,JB+1,3)
              IB(IDIRB) = IB(IDIRB) + 2
              KCPB = INDHER(IB(1),IB(2),IB(3))
              FMATCP(KCA,KCB) = FMATCP(KCA,KCB) + FMATC2P(KCPA,KCPB)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FMATC2P)
      RETURN
      END
C
C  FMAT:    - AO is in Spherical Contracted (full)
C        CP - Cartesian Primitive (submatrix)
C        SP - Spherical Primitive (submatrix)
C
C  So far we do not take care about symmetry.
C  FMAT is assumed to be full dimension without any 
C  symmetry reduction
C
      SUBROUTINE ADD2FCK31(FMAT,FMATCP,FMATSP,FTMP,
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
#include "densfit.h"
C
#include "pbc.h"
C
      INTEGER T,U,V,TUV
      DIMENSION FMAT(NDIMA,NDIMB,NDMAT),FMATCP(KCKTA,KCKTB,NPRFAB),
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

          KNDX = KNDXBTCH(ICOMPA,IAOBA,1)
          IAOA = INDXBTCH(KNDX - 1 + ICTFA,1)
          
          IOFFSET = 0
          KNDX = KNDXBTCH(ICOMPB,IAOBB,1)
          IAOB = INDXBTCH(KNDX - 1 + ICTFB,1) - IOFFSET
C
C Here it we should probably put summation over 
C  FMATSP(ISPA,ISPB,IPRFAB) separately outside 
C 

C         IF ( (IAOA .EQ. IAOB ) .OR. ( IAOBA .EQ. IAOBB )) THEN 
#ifdef MOD_PBC
         if (dopbc .and. pbcactive) then
          DO ISPA=1, KHKTA
           DO ISPB=1, KHKTB
             FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMAT) = 
     &            FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMAT) 
     &          + FMATSP(ISPA,ISPB,IPRFAB) 
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
           END DO
          END DO
         else IF ( ( IAOBA .EQ. IAOBB )) THEN 
#else
         IF ( ( IAOBA .EQ. IAOBB )) THEN 
#endif /* MOD_PBC */
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
     &               FMAT(IAOB+ISPB-1,IAOA+ISPA-1,IDMAT) 
     &               + FMATSP(ISPA,ISPB,IPRFAB) 
     &               * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
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
#include "implicit.h"
#include "priunit.h"
      DIMENSION IFBTCH(2,NFBTCH), RFBTCH(6,NFBTCH)

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
#include "fck3odb.h"

      CHCKNAI   = .FALSE.
      NOSCRPRIM = .FALSE.
      NOSKIPCL  = .FALSE.
      NOSKIPCS  = .FALSE.
      SCRFBTCH  = .FALSE.
      OLDMOM    = .FALSE.
      FIXTHR    = .FALSE.
      BRFREE    = .FALSE.
      FULL_J    = .FALSE.
      HERECF    = .FALSE.

      SPEXP     = .TRUE.
      PACKFTUV  = .TRUE.
      SKALLCL   = .TRUE.
      CONTRACTED= .TRUE.

      NCSODB    = .FALSE.
      NCSTUV    = .FALSE.

C
C aspg, 2006-06-27 making ftnchek happy
C
      IPRFDEF   = 0
C
      IPRFCK3   = 0 

      THRFCK3   = THRFCK3_DEFAULT
      THRNCL    = THRFCK3_DEFAULT 
      THRPOD    = THRFCK3_DEFAULT 
      THRCOD    = THRFCK3_DEFAULT 
      THRSPEXP  = THRFCK3_DEFAULT 


      SYMDENS = .TRUE.

      RETURN
      END

      SUBROUTINE FCK3INP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "cbifck3.h"
#include "fck3thr.h"
      PARAMETER ( NTABLE = 18 )
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DATA TABLE /'.FIXTHR','.THRFCK','.PRINT ','.NOSPEX','.NOSCRP',
     &            '.NOSKCL','.NOSKCS','.SCRFBT','.NOPACK','.NOSKAL',
     &            '.FULL_J','.UNCONT','.BRFREE','.CHCKNA','.NCSODB',
     &            '.NCSTUV','XXXXXXX','.HERMIT'/
C
      IPRINT=IPRFCK3
C
C aspg, 2006-06-27 making ftnchek happy
C
      IPRFDEF=0
C
      NEWDEF = (WORD .EQ. '*JENGIN')
      IF(NEWDEF) THEN
        WORD1 = WORD
C
 100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)

            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     &                      15,16,17,18),i
                  END IF
 200          CONTINUE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FCK3INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FCK3INP.')
 1    CONTINUE
        FIXTHR     = .TRUE.
      GO TO 100
 2    CONTINUE
        READ (LUCMD, * ) THRFCK3
      GO TO 100
 3    CONTINUE
        READ (LUCMD, '(I5)') IPRFCK3
      GO TO 100
 4    CONTINUE
        SPEXP      = .FALSE.
      GO TO 100
 5    CONTINUE 
         NOSCRPRIM = .TRUE.
      GO TO 100
 6    CONTINUE
         NOSKIPCL  = .TRUE.
      GO TO 100
 7    CONTINUE
         NOSKIPCS  = .TRUE.
      GO TO 100
 8    CONTINUE
         SCRFBTCH  = .TRUE.
      GO TO 100
 9    CONTINUE
         PACKFTUV  = .FALSE.
      GO TO 100  
 10   CONTINUE
         SKALLCL   = .FALSE.
      GO TO 100  
 11   CONTINUE
         FULL_J    = .TRUE.
         FIXTHR    = .TRUE.
         SPEXP     = .FALSE.
      GO TO 100  
 12   CONTINUE
         CONTRACTED= .FALSE.
      GO TO 100  
 13   CONTINUE
         BRFREE    = .TRUE.
         CONTRACTED= .FALSE.
      GO TO 100  
 14   CONTINUE
         CHCKNAI = .TRUE.
      GO TO 100  
 15   CONTINUE
         NCSODB  = .TRUE.
      GO TO 100  
 16   CONTINUE
         NCSTUV  = .TRUE.
      GO TO 100  
 17   CONTINUE
      GO TO 100
 18   CONTINUE
         HERECF  = .TRUE.
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

      IF (IPRINT .NE. IPRFDEF) THEN
          WRITE (LUPRI,'(A,I5)') ' Print level in FCK3:',IPRFCK3
      END IF

      IF ( FIXTHR ) THEN
        WRITE (LUPRI,'(A,I5)') 'Threshold is fixed.'
      ELSE
#if defined (PRG_DIRAC)
        WRITE (LUPRI,'(A,I5)') 
     &    'Dynamical threshold is not (yet) implemented in DIRAC'
#else
        WRITE (LUPRI,'(A,I5)') 'Threshold set dynamically.'
#endif
      END IF

      IF ( NOSKIPCL ) THEN
        WRITE (LUPRI,'(A,I5)') 
     &  'No skipping of the nonclassical integrals (only set to 0).'
      END IF

      RETURN
      END


      SUBROUTINE CHECKR2MIN( R2MIN, RODP, ICNTC,ICNTD, RCD )
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aobtch.h"
#include "nuclei.h"

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
c     
c     MXETIJ finds largest factor of x*y*z components of ETIJ matrix
c     and returns them in ETIJMX. We collapse
c     first four dimenstions of ETIJ into one for performance.
      SUBROUTINE MXETIJ(ETIJ,NPR12,JMAX1,JMAX2,ETIJMX)
#include "implicit.h"
      DIMENSION ETIJ(NPR12*(JMAX1+JMAX2+1)*(JMAX1+1)*
     &               (JMAX2+1),3)
c
      ETIJMX = 0D0
      NLEN = NPR12*(JMAX1+JMAX2+1)*(JMAX1+1)*(JMAX2+1)
      DO J= 1, NLEN
         ETIJMX = MAX(ETIJMX,
     &        ABS(ETIJ(J,1)*ETIJ(J,2)*ETIJ(J,3)))
      END DO
      END
c
c     R1UPPERB is a routine for finding an upper bound for R-integrals,
c     for given orbital exponent and 
      FUNCTION R1UPPERB(N,ALPHA)
#include "implicit.h"
#include "sdpre.h"
      UL  = 0D0
      ALPHA2  = ALPHA
      ALPHA2I = 1D0
      FAC = -1D0
      DO I = 0, N
         FAC     = FAC*SDPRE(2*I-1)
         UL      = MAX(UL,ABS(ALPHA2I*FAC/SDPRE(2*I+1)))
         ALPHA2I = ALPHA2I*ALPHA2
      END DO
      R1UPPERB = ul
      END
c     /* deck dmatdim */
      SUBROUTINE DMATDIM
#include "implicit.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"

        NDIMA = NBASE
        NDIMB = NBASE
        NDIMC = NBASE
        NDIMD = NBASE

      RETURN
      END

c     Picks the largest exponent of the ab-distribution
      SUBROUTINE GETMAXEXPAB(EXPMAB,EXPBT1,EXPBT2,NPRF1,NPRF2)
#include "implicit.h"
      DIMENSION EXPBT1(NPRF1),EXPBT2(NPRF2)
      EXPMAB=0.0D0
      DO I=1,NPRF1
        EXP1   = EXPBT1(I)
        DO J=1,NPRF2
          EXP2   = EXPBT2(J)
          EXP12  = EXP1 + EXP2 
          EXPMAB = MAX(EXPMAB,EXP12)
        ENDDO
      END DO
      RETURN
      END

      SUBROUTINE CHECKRMIN( RMIN, RODP, ICNTC,ICNTD, RCD )
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aobtch.h"
#include "nuclei.h"
#include "densfit.h"
C
#ifdef MOD_PBC
#include "pbc.h"
#endif /* MOD_PBC */
C
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

#ifdef MOD_PBC
      if (dopbc .and. pbcactive) then
C        Here PX,PY,PZ refer to the ab overlap distributions, while
C        CX,CY,CZ,DX,DY,DZ refer to the c and d orbitals.
C
C        The coordinates RODP already have the appropriate translations
C        due to previous calls to GETR12, so PX,PY,PZ are correct.
C        But CX,CY,CZ,DX,DY,DZ are taken from the global array CORD
C        and need to have translations added to them.
         CX = CX + nf_vec_std(1)
         CY = CY + nf_vec_std(2)
         CZ = CZ + nf_vec_std(3)
         DX = DX + nf_vec_std(1) + nf2_vec_std(1)
         DY = DY + nf_vec_std(2) + nf2_vec_std(2)
         DZ = DZ + nf_vec_std(3) + nf2_vec_std(3)
      end if
#endif /* MOD_PBC */

      RPC = SQRT((PX-CX)**2 + (PY-CY)**2 + (PZ-CZ)**2)
      RPD = SQRT((PX-DX)**2 + (PY-DY)**2 + (PZ-DZ)**2)

      R = RPC + RPD - RCD
      IF (RHSFIT) THEN
        R = RPD
      ENDIF
     
      IF ( R .GT. 0.0D0 ) THEN
         IF (RHSFIT) RMIN = R
         ELSE RMIN = 0.5D0*R
      ELSE
         RMIN = 0.0D0
      END IF
      RETURN
      END 

      SUBROUTINE INIFCK3REL()
#include "implicit.h"
#include "priunit.h"
#include "cbifck3.h"
#include "fck3rel.h"
C     FCK3PATH:  1 only on first (ab) electron
C                2 only on second (ftuv) electron
C                3 on both electrons electron
      NRELLHS = 0
      NRELRHS = 0
      NDERLHS = 0
      NDERRHS = 0
      IF (DOFCK3REL) THEN
        FCK3PATH = 1
        IF (FCK3PATH.EQ.1) THEN
          NRELLHS = 2
          NRELRHS = 0
          NDERLHS = 1
          NDERRHS = 0
        ELSE
          WRITE(LUPRI,'(1X,A,I2,A)') 
     &      'Wrong input in INIFCK3REL. FCK3PATH =',FCK3PATH,
     &      ' not a case for relativistic calculations.'
          CALL QUIT('Relativistic (LL|LL) and (SS|LL) only')
        ENDIF
        IF (.NOT.HERECF) 
     &     CALL QUIT('INIFCK3REL error, HERECF must be true!')
      ENDIF
      NREL = MAX(NRELLHS,NRELRHS)
      NDER = MAX(NDERLHS,NDERRHS)
      RETURN
      END
      
