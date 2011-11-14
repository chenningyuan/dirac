!dalton_copyright_start
!      Copyright (c) 2011 by the authors of Dalton (see below).
!      All Rights Reserved.
!   
!      The source code in this file is part of
!      "Dalton, a molecular electronic structure program,
!       Release DALTON2011 (2011), see http://daltonprogram.org"
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

#if defined (OLFDCK3DRV)
      SUBROUTINE OLDFCK3DRV(FMAT,DMAT,NDMT,WORK,LWORK)
C
C   Calculates the nonclassical coulomb contribution to the
C   Fock matrix
C
#include "implicit.h"
      DIMENSION FMAT(*),DMAT(*), WORK(*)
C
C
#include "aovec.h"
#include "maxash.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#if ! defined (PRG_DIRAC)
#include "infind.h"
#endif
#include "inftap.h"
#include "infdim.h"
#if ! defined (PRG_DIRAC)
#include "infpri.h"
#endif
#include "priunit.h"
#include "cbihr2.h"
#include "cbieri.h"
#include "cbifck3.h"
#include "erithr.h"
#include "clsfmm.h"
#if defined (VAR_MPI)
      INCLUDE "mpif.h"
#include "infpar.h"
#endif

C
      DIMENSION ISYMDM(2), IFCTYP(2)
      SAVE ISYMDM, IFCTYP
      DATA ISYMDM/1,1/, IFCTYP/13,13/
#include "memint.h"

      NDMAT = NDMT

      CALL DZERO(FMAT,NDMAT*N2BASX)
C
      IF ( FIXTHR ) THEN
        THRSH = 1.00D-15
      ELSE
        THRSH = 10D0**(-IFTHRS)
      END IF

      CALL MEMGET('REAL',KSAO,  NNBAST,         WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCCFBT,MXPRIM*MXCONT,  WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDXB,8*MXSHEL*MXCONT,WORK,KFREE,LFREE)

#ifdef VAR_MPI
c wake up slaves.
       IPRTYP = 5
       IPRINT = 0
       CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &      MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &      MPI_COMM_WORLD,IERR)
       CALL fck3_dims(NDMAT)
       CALL FCK3SN1(DMAT,NDMAT,ISYMDM,IFCTYP,WORK(KSAO))
#endif
       CALL HRFCK3(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,IPRFCK,
     &             WORK(KSAO),WORK(KCCFBT),WORK(KINDXB),
     &             WORK(KFREE),LFREE)

      IF ( IPRFCK .GE. 4 ) THEN
         CALL HEADER('Density and Fock matrices from NCLFCK',-1)
         KSTR = 1
         DO I = 1, NDMAT
          WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
          CALL OUTPUT(DMAT(KSTR),1,NBAST,1,NBAST,NBAST,
     &                  NBAST,1,LUPRI)
#ifdef VAR_MPI
          WRITE (LUPRI,'(//,1X,I2,A,I3)') MYNUM,' Fock matrix No.',I
#else
          WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
#endif
          CALL OUTPUT(FMAT(KSTR),1,NBAST,1,NBAST,NBAST,
     &                   NBASE,1,LUPRI)
          KSTR = KSTR + NBAST*NBAST
         END DO
      END IF

      CALL FLSHFO(LUPRI)

      RETURN
      END
#endif /* OLDFCK3DRV */

      SUBROUTINE HRFCK3(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     *                  SAO,CCFBT,INDXBT,WORK,LWORK)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
#include "ccom.h"
#include "cbieri.h"
#include "cbione.h"
#include "cbihr2.h"
#include "cbifck3.h"
#include "ericom.h"
#include "erithr.h"
#include "erimem.h"
#include "aobtch.h"
#include "odbtch.h"
#include "symmet.h"
#include "inforb.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "clsfmm.h"
#include "shells.h"

      LOGICAL FOUND
      DIMENSION FMAT(*), DMAT(*), IFCTYP(NDMT),ISYMDM(NDMT),
     &          SAO(NBAST,NBAST),CCFBT(*),INDXBT(*),WORK(LWORK)
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C   This should be checked to be consistent with the fmm part ...
C
      THRCLS = THRSH
#ifdef MOD_FMM
      IF (DOFMM) THEN
         NCLONE = .TRUE.
         NCLERI = .TRUE.
      END IF
#endif
C
C        Initialization in ER2INI
C

      CALL ER2INI

      ANCLS_CALC = 0.0D0
      ANCLS_SKIP = 0.0D0
      ANCLS_CALC0 = 0.0D0

      IPRINT = IPRFCK3

      NDMAT  = NDMT
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C
      CALL ODBCHS(KODCL1,KODCL2,
     &            KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C
C
C     Select integrals to be calculated
C     =================================
C
CC#ifndef PRG_DIRAC
      CALL PICKAO(IPRINT)
CC#else
CC      CALL PICKAO(NCENTS,ICENTS,ICLASSES,IPRINT)
CC#endif
C
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(INDXBT,IPRINT)

C
C  ( The value IELECTR=1 is only for testing of the nuclear
C    attraction value )

C      IELECTR = 1
C
C      CALL GETFTUV(KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
C     &             KRODQ,DMAT,ISYMDM,NDMT,
C     &             WORK(KODCL2),WORK(KODBC2),WORK(KRDBC2),
C     &             WORK(KODPP2),WORK(KRDPP2),
C     &             KFREE,LFREE,CCFBT,INDXBT,
C     &             IELECTR,WORK,LWORK,IPRINT)


C
C  ( The value IELECTR = 2 is for the real calculation
C    of the Fock matrix; FTUV than includes also
C    the factor (-1)^(T+U+V) which we do not include
C    in the calse IELECTR = 1. All the other changes
C    between electron 1 and 2 are done in expansion
C    coefficients ETIJ.
C

      IELECTR = 2

      CALL GETTIM(FTUVTIM_START,DUMTIM)

      CALL GETFTUV(KIFTUVCL,KRFTUVCL,KIFTUV,KRFTUV,KFTUV,
     &             KRODQ,DMAT,ISYMDM,NDMT,
     &             WORK(KODCL2),WORK(KODBC2),WORK(KRDBC2),
     &             WORK(KODPP2),WORK(KRDPP2),
     &             KFREE,LFREE,CCFBT,INDXBT,
     &             IELECTR,WORK,LWORK,IPRINT)

      CALL GETTIM(FTUVTIM_FTUV,DUMTIM)

      CALL GETFCK_ERI(FMAT,NDMT,IFCTYP,WORK(KIFTUVCL),WORK(KRFTUVCL),
     &            WORK(KIFTUV),WORK(KRFTUV),WORK(KFTUV),
     &            WORK(KRODQ),
     &            WORK(KODCL1),WORK(KODCL2),
     &            WORK(KODBC1),WORK(KODBC2),
     &            WORK(KRDBC1),WORK(KRDBC2),
     &            WORK(KODPP1),WORK(KODPP2),
     &            WORK(KRDPP1),WORK(KRDPP2),
     &            CCFBT,INDXBT,
     &            WORK(KFREE),LFREE,IPRINT)

      CALL FLSHFO(LUPRI)

      CALL GETTIM(FTUVTIM_END,DUMTIM)

      IF ( IPRINT .GE. -1 ) THEN
       WRITE (LUPRI,'(/,2(5X,A,E8.3,/))')
     &      'Screenig threshold in FCK3 : ', THRSH,
     &      'Nonclassical threshold     : ', THRCLS


       WRITE (LUPRI,'(/,2(5X,A,F16.2,A,/))')
     &             'Time in generating FTUV :',
     &              FTUVTIM_FTUV - FTUVTIM_START,' seconds',
     &             'Time in FCK3          :',
     &              FTUVTIM_END - FTUVTIM_FTUV,' seconds'
       ANCLS_TOTAL = ANCLS_SKIP + ANCLS_CALC
       ANCLS_CALC  = 1.0D2 * ANCLS_CALC  / ANCLS_TOTAL
       ANCLS_SKIP  = 1.0D2 * ANCLS_SKIP  / ANCLS_TOTAL
       ANCLS_CALC0 = 1.0D2 * ANCLS_CALC0 / ANCLS_TOTAL
       ANCLS_CALCN0= ANCLS_CALC - ANCLS_CALC0
       WRITE (LUPRI,'(/,A,/,3(5X,A,F10.4,A,/))')
     &          ' Percentage of nonclassical integrals: ',
     &          ' - calculated nonzero  : ', ANCLS_CALCN0, ' %',
     &          ' - calculated zero     : ', ANCLS_CALC0, ' %',
     &          ' - skipped             : ', ANCLS_SKIP, ' %'

       WRITE (LUPRI,'(/,A,E10.3,/)')
     &     'Total number of primitives OD processed ', ANCLS_TOTAL

      END IF

      IF ( IPRINT .GE. 2 ) THEN
        CALL HEADER('Density and Fock matrices in FCK3',-1)
        KSTR = 1
        DO I = 1, NDMT
          WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
          CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                  NBASE,1,LUPRI)
          WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
          CALL OUTPUT(FMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                   NBASE,1,LUPRI)
          KSTR = KSTR + NBASE*NBASE
        END DO
      END IF

      CALL TIMER('FCK3  ',TIMSTR,TIMEND)

      CALL FLSHFO(LUPRI)

      RETURN
      END

      SUBROUTINE GETFCK_ERI(FMAT,NDMT,IFCTYP,
     &             IFTUVCL,RFTUVCL,IFTUV,RFTUV,FTUV,
     &             RODQ,IODCL1, IODCL2,
     &             IODBC1,IODBC2,RODBC1,RODBC2,
     &             IODPP1,IODPP2,RODPP1,RODPP2,
     &             CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "ericom.h"
#include "cbieri.h"
#include "aobtch.h"
#include "symmet.h"
#include "odbtch.h"
#include "odclss.h"
#include "eriao.h"
C
      DIMENSION FMAT(*), IFCTYP(NDMT),
     &          IFTUVCL(10,*),RFTUVCL(2,*),
     &          IFTUV(*),RFTUV(*),FTUV(*),RODQ(4,*),
     &          IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          CCFBT(*), INDXBT(*), WORK(LWORK)

      IF ( IPRINT .GE. 2 ) THEN
        CALL HEADER('Entering GETFCK in eri2newfck ',-1)
      END IF

C      WRITE(LUPRI,*) 'RODQ',(RODQ(4,I),I=1,4)

      MAXDER = 0
      IF (GDER) MAXDER = 1
      IF (BDER) MAXDER = 1
C
      NHTYP = NHKMAX

      JTOP  = 4*(NHTYP - 1) + MAXDER
      IF (BDER) JTOP = JTOP + 1
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
C
      NHKMAX = NHTYP + MAXDER
      KCKMAX = NHKMAX*(NHKMAX + 1)/2
      KC2MAX = KCKMAX**2
      NRDER  = 3*MAXDER

C
C     Dimensions for CSQ
C     ==================
C
      CALL SPHDIM
C
C     Allocate
C     ========
C
      KINDHR = 1
      KINDSQ = KINDHR + (JTOP3 - 1)/IRAT + 1
      KIODHR = KINDSQ + (NRTOP - 1)/IRAT + 1
      KLMNVL = KIODHR + (8*NRTOP - 1)/IRAT + 1
      KLMNSM = KLMNVL + (3*NHKMAX*KCKMAX - 1)/IRAT + 1
      KCSQ   = KLMNSM + (16*NHKMAX*KCKMAX - 1)/IRAT + 1
      KLAST  = KCSQ   + NCSQ1*NCSQ2

C
      IF (KLAST .GT. LWORK) CALL STOPIT('ODCDRV','HERPRP',KLAST,LWORK)
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

      CALL GETFCK1(IODCL1,IODCL2,
     &             IODBC1,IODBC2,RODBC1,RODBC2,
     &             IODPP1,IODPP2,RODPP1,RODPP2,
     &             FMAT, DMAT,NDMT,IFCTYP,
     &             IFTUVCL,RFTUVCL,IFTUV,RFTUV,FTUV,RODQ,
     &             WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR),
     &             WORK(KLMNVL),WORK(KLMNSM),
     &             CCFBT,INDXBT,
     &             WORK(KCSQ),  WORK(KLAST),LWRK,IPRINT)

      RETURN
      END

C  /* Deck odcdr1 */
      SUBROUTINE GETFCK1(IODCL1,IODCL2,
     &                  IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  FMAT,DMAT,NDMT,IFCTYP,
     &                  IFTUVCL,RFTUVCL,IFTUV,RFTUV,FTUV,RODQ,
     &                  INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,
     &                  CCFBT,INDXBT,CSQ,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "eribuf.h"
#include "dummy.h"
#include "cbieri.h"
#include "erithr.h"
#include "odclss.h"
#include "ccom.h"
#include "aobtch.h"
#include "doxyz.h"
#include "eriao.h"
#include "ericom.h"
#include "erimem.h"
#include "erista.h"
#include "hertop.h"
#include "symmet.h"
#include "iratdef.h"
#include "pi.h"
#include "clsfmm.h"
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          FMAT(*),DMAT(*),IFCTYP(NDMT),
     &          IFTUVCL(10,*),RFTUVCL(4,*),RODQ(4,*),
     &          IFTUV(*),RFTUV(*),FTUV(*),
     &          CCFBT(*), INDXBT(*),
     &          CSQ(NCSQ1,NCSQ2), WORK(LWORK)
      LOGICAL CLASSIC
#if defined(VAR_MPI)
      INCLUDE 'mpif.h'
#include "infpar.h"
#include "ibtfun.h"
      NOFF  = MYNUM
      NSTEP = NODTOT+1
#else
      PARAMETER (NSTEP=1)
#include "ibtfun.h"
      NOFF=0
#endif

C      DO I=1,NODBC2
C        WRITE(LUPRI,'(A,I,4F12.8)') 'RODQ: ', I,
C     &               RODQ(1,I),RODQ(2,I),RODQ(3,I),RODQ(4,I)
C      END DO

      FAC0 = 2.0D0*PI*PI*SQRT(PI)
      FACCLS = ERFCIV(THRCLS)

      CALL GETTIM(tsGETFCK1,DUMTIM)
C      WRITE(LUPRI,'(A,6I6)') 'IFTUV-beg.: ', (IFTUV(I),I=1,6)

      KRODP = 1
      KWRK0 = KRODP  + 4
      DO 500 KODCL1=1,NODCL1

       NODTAB = IODCL1(KODCL1,1)
       KODSAB = IODCL1(KODCL1,2)

       NHKTA  = IODCL1(KODCL1,3)
       NHKTB  = IODCL1(KODCL1,4)
       NHKTAB = NHKTA + NHKTB - 1

       KHKTA  = IODCL1(KODCL1,5)
       KHKTB  = IODCL1(KODCL1,6)
       KCKTA  = IODCL1(KODCL1,7)
       KCKTB  = IODCL1(KODCL1,8)

       KCSQA = KSQADR(NHKTA-1,0,0,0,0)
       KCSQB = KSQADR(NHKTB-1,0,0,0,0)

       ICMATA = IODCL1(KODCL1,19)
       ICMATB = IODCL1(KODCL1,20)

       JMAXA  = NHKTA - 1
       JMAXB  = NHKTB - 1
       JMAXAB = JMAXA + JMAXB

       JMAX1  = JMAXA
       JMAX2  = JMAXB
       JMAX12 = JMAXAB

       NTUVAB = (JMAXAB + 1)*(JMAXAB + 2)*(JMAXAB + 3)/6

       NPRFA  = IODCL1(KODCL1,15)
       NPRFB  = IODCL1(KODCL1,16)

       NPRFAB  = IODCL1(KODCL1,11)
       NPRFABCLS = IODCL1(KODCL1,11)

       NCTFA  = IODCL1(KODCL1,17)
       NCTFB  = IODCL1(KODCL1,18)

       NPRFAB = NPRFA*NPRFB

       NPPAB  = NPRFAB
       NPP12  = NPRFAB

       KECOEFAB = KWRK0
       KHRINT   = KECOEFAB+ NPPAB * NHKTAB*NHKTA*NHKTB*3*2
       KCOORAB  = KHRINT  + NPPAB * NTUVAB*NDMT
       KCOORP   = KCOORAB + NPPAB * 3
       KEXPAB   = KCOORP  + NPPAB * 3
       KEXPABP  = KEXPAB  + NPPAB * 3
       KPREFAB  = KEXPABP + NPPAB * 3
       KHPIAB   = KPREFAB + NPPAB
       KIPRAB   = KHPIAB  + NPPAB
       KPA      = KIPRAB  + NPPAB * 2 / IRAT + 1
       KPB      = KPA     + NPPAB * 3
       KFMATCP  = KPB     + NPPAB * 3
       KFMATSP  = KFMATCP + NPPAB * KCKTA * KCKTB
       KFTMP    = KFMATSP + NPPAB * KHKTA * KHKTB
       KWRK1    = KFTMP   + KCKTA * KHKTB

C        WRITE(LUPRI,*) 'KECOEFAB, KCOORAB, KEXPAB, KPREFAB, ',
C     &                 'KHPIAB, KPA, KPB, KWRK1 ',
C     &                  KECOEFAB, KCOORAB, KEXPAB, KPREFAB,
C     &                  KHPIAB, KPA, KPB, KWRK1


       DO 400 KODBC1=KODSAB+NOFF,KODSAB+NODTAB-1,NSTEP

        IA  = IODBC1(KODBC1,2)
        IB  = IODBC1(KODBC1,3)

        CALL DZERO(WORK(KECOEFAB),NPPAB*NHKTAB*NHKTA*NHKTB*3*2)
        CALL DZERO(WORK(KHRINT),NPPAB*NTUVAB*NDMT)

        IELECTR1 = 1

        CALL SCREENAB(IA,IB,CCFBT(ICMATA),CCFBT(ICMATB),
     &                RAB,WORK(KIPRAB),WORK(KRODP),FACCLS,KODBC1)

C        WRITE(LUPRI,*) 'Before GETP:',NPRFAB,NPPAB,NPP12


        CALL GETP(WORK(KCOORP),WORK(KEXPAB),WORK(KEXPABP),
     &            WORK(KPREFAB),WORK(KIPRAB),
     &            WORK(KECOEFAB),WORK(KCOORAB),
     &            WORK(KHPIAB),WORK(KPA),WORK(KPB),
     &            FAC0,IODBC1,RODBC1,IODPP1,RODPP1,
     &            IELECTR1,KODBC1,IPRINT)

        EXTPMAX = WORK(KRODP+3)
        SCRAB = RODBC1(KODBC1,1)

        DO 300  IFCL=1,NFTUVCL

          NTFTUV = IFTUVCL(  1, IFCL )
          KRFTUV = IFTUVCL(  2, IFCL )
          KSFTUV = IFTUVCL(  3, IFCL )
          JMAXF  = IFTUVCL(  4, IFCL )
          NHKTF  = IFTUVCL(  5, IFCL )
          NCNTF  = IFTUVCL(  6, IFCL )
          ICNTC  = IFTUVCL(  7, IFCL )
          ICNTD  = IFTUVCL(  8, IFCL )
          KIFTUV = IFTUVCL(  9, IFCL )
          NODF   = IFTUVCL( 10, IFCL )

          RCD     = RFTUVCL( 1, IFCL )
          EXTQMAX = RFTUVCL( 2, IFCL )
          SCRCD   = RFTUVCL( 4, IFCL )

          CALL CHECKRMIN( RMIN, WORK(KRODP), ICNTC,ICNTD, RCD )

C           WRITE(6,*)   SCRAB,SCRCD, THRSH, 'SKIP'

          IF (.NOT.( NCLERI .AND. RMIN .GT. EXTPMAX+EXTQMAX )
     &         .AND. ( SCRAB*SCRCD .GT. THRSH ) )  THEN
C          IF ( SCRAB*SCRCD .GT. THRSH ) THEN

           NHKTCD = NHKTF
           JMAXCD = JMAXF

           NPPCD  = NTFTUV
           NPPX   = NPPAB*NPPCD

           ANCLS_CALC = ANCLS_CALC + NPPX

           JMAX = JMAXAB + JMAXCD

           NTUV   = (JMAX + 1)*(JMAX + 2)*(JMAX + 3)/6

           KRJ000   = KWRK1
           KHRINT1  = KRJ000  + NPPX * (JMAX+1)
           KFACINT  = KHRINT1 + NPPX * NTUV
           KALPHA   = KFACINT + NPPX
           KCOORPQ  = KALPHA  + NPPX
           KWRK2    = KCOORPQ + NPPX * 3

           LWRK2    = LWORK - (KWRK2-KWRK0)

           IF (LWRK2 .LE. 0) CALL QUIT('NOT ENOUGH MEM!')

           CALL GETPQ(KODBC1,WORK(KCOORP), WORK(KEXPABP),
     &               RFTUV(KRFTUV),IFTUV(KIFTUV),
     &               WORK(KCOORPQ),WORK(KFACINT),WORK(KALPHA),
     &               WORK(KRODP),RODQ,NODF,IPRINT)

C        Incomplete gamma function
C        =========================

          CALL ERIGAM(WORK(KRJ000),WORK(KFACINT),WORK(KALPHA),
     &               WORK(KCOORPQ),WORK(KWRK2),LWRK2,IPRINT)

C
C         Special cases
C         =============

C          IF ( JMAX .LE. 1 ) THEN
C            CALL HRFSPEC(WORK(KHRINT),WORK(KRJ000),WORK(KCOORPQ),
C     &             FTUV(KSFTUV),INDHER,IODDHR(1,IPQXYZ),NDMT,IPRINT)
C          ELSE

C        Hermite integrals
C        =================

          CALL ERITUV(WORK(KHRINT1),WORK(KRJ000),
     &               WORK(KCOORPQ),INDHER,IODDHR(1,IPQXYZ),
     &               WORK(KWRK2),LWRK2,IPRINT)


C        Contraction of Hermite integrals with FTUV
C        ==========================================


          CALL HRFTUV(WORK(KHRINT),WORK(KHRINT1),
     &               FTUV(KSFTUV),NDMT,
     &               INDHER,IODDHR(1,IPQXYZ),IPRINT)

C          END IF

         ELSE
          ANCLS_SKIP = ANCLS_SKIP + NPPAB*NPPCD
         END IF

 300    CONTINUE

C
C        Contribution to the Fock matrix ...
C

         CALL ADDFCK3(FMAT,WORK(KFMATCP),WORK(KFMATSP),
     &         WORK(KFTMP),NDMT,
     &         WORK(KECOEFAB),WORK(KPREFAB),WORK(KIPRAB),
     &         WORK(KHRINT),INDHER,LMNPWR,IODDHR(1,IPQXYZ),
     &         CSQ(KCSQA,1),CSQ(KCSQB,1),
     &         CCFBT(ICMATA),CCFBT(ICMATB),INDXBT,NPPAB,
     &         JMAXA, JMAXB, JMAXAB, NTUV,
     &         KHKTA, KHKTB, KCKTA, KCKTB,
     &         NPRFA, NPRFB, NCTFA, NCTFB,
     &         KCKMAX,NHKMAX,IA,IB,IPRINT )

 400  CONTINUE
        NOFF = MOD(NOFF+NODTAB,NSTEP)
 500  CONTINUE
      CALL GETTIM(teGETFCK1,DUMTIM)

#if defined(VAR_MPI)
      CALL fck3collect(FMAT,NDMAT,teGETFCK1-tsGETFCK1,WORK,LWORK,IPRINT)
#endif
      RETURN
      END

      SUBROUTINE GETPQ(KODBC1,COORP,EXPABP,RFTUV,IFTUV,
     &                 COORPQ,FACINT,ALPHA,RODP,RODQ,NODF,IPRINT)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "ericom.h"
#include "aobtch.h"
#include "odclss.h"
#include "cbieri.h"
#include "erithr.h"
      DIMENSION COORP(3,NPPAB),EXPABP(3,NPPAB),
     &          RFTUV(7,NPPCD),IFTUV(2,NPPCD),COORPQ(NPPX,3),
     &          FACINT(NPPX),ALPHA(NPPX),RODP(4),RODQ(4,*)

      INTEGER X

      IF (IPRINT .GE. 10) THEN
        CALL HEADER('OUTPUT FROM GETPQ',-1)
        WRITE (LUPRI,'(5X,6(A,I6,/),2(2(A,I6),/)/)')
     &         'KODCL2:',KODCL2,
     &         'NODTCD:',NODTCD,'KODSCD:',KODSCD,
     &         'NPPAB :',NPPAB,'NPPCD:',NPPCD,'NPPX:',NPPX,
     &         'NPRFA :',NPRFA,'NPRFB :',NPRFB,
     &         'NPRFC :',NPRFC,'NPRFD :',NPRFD

      END IF

      IPPX=1

      DO 300 IPPAB=1,NPPAB

       EXPP  = EXPABP(1,IPPAB)
       EXPPI = EXPABP(2,IPPAB)
       FACPI = EXPABP(3,IPPAB)

       DO IPPCD=1,NPPCD

        EXPQ  = RFTUV(1,IPPCD)
        EXPQI = RFTUV(2,IPPCD)

        QX = RFTUV(3,IPPCD)
        QY = RFTUV(4,IPPCD)
        QZ = RFTUV(5,IPPCD)

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

C        WRITE(LUPRI,'(A,5F10.6)') 'Q-coord',
C     &       QX,QY,QZ,EXPQ,EXPQI


         ALPHA(IPPX) = - 2.0D0 / (EXPPI + EXPQI)

         FACINT(IPPX)= FACPI*EXPQI/SQRT(EXPP+EXPQ)

         COORPQ(IPPX,1)=COORP(1,IPPAB)-QX
         COORPQ(IPPX,2)=COORP(2,IPPAB)-QY
         COORPQ(IPPX,3)=COORP(3,IPPAB)-QZ

         IPPX = IPPX + 1

        END DO

 300   CONTINUE

      IF ( NCLERI ) THEN


        IPPCD = 1

        DO IODF=1,NODF

         KQ   = IFTUV( 1, IPPCD )
         NPRF = IFTUV( 2, IPPCD )
C         WRITE(LUPRI,*) 'NODF: ', NODF, KQ, NPRF, IPPX,IPPCD

         RPQ2 = (RODP(1)-RODQ(1,KQ))**2 + (RODP(2)-RODQ(2,KQ))**2
     &        + (RODP(3)-RODQ(3,KQ))**2
         RPQ  = SQRT(RPQ2)
         EPQ  = RODP(4) + RODQ(4,KQ)

C         WRITE(LUPRI,*) 'RPQ,EPQ',RPQ,EPQ,KQ
C          WRITE(LUPRI,'(A,I5,10F12.8)') 'RPQ,EPQ...',KQ,
C     &                  RPQ,EPQ,RODP(4),RODQ(4,KQ),
C     &     RODP(1),RODP(2),RODP(3),RODQ(1,KQ),RODQ(2,KQ),RODQ(3,KQ)


         IF (RPQ.GT.EPQ) THEN
C          WRITE(LUPRI,'(A,10F12.8)') 'CLASSIC1',
C     &                  RPQ,EPQ,RODP(4),RODQ(4,KQ),
C     &     RODP(1),RODP(2),RODP(3),RODQ(1,KQ),RODQ(2,KQ),RODQ(3,KQ)

          IPPX=IPPCD
          DO IPPAB=1,NPPAB
           DO I=IPPX,IPPX+NPRF-1
            FACINT(I)=0.0D0
           END DO
           IPPX=IPPX+NPPCD
          END DO
          ANCLS_CALC0 = ANCLS_CALC0 + NPRF*NPPAB

         END IF

         IPPCD = IPPCD + NPRF
        END DO
C
C
C           WRITE(LUPRI,'(2(A,4F12.6))')
C     &      'NCLS:P:',RODP(1),RODP(2),RODP(3),RODP(4),
C     &      'NCLS:Q:',RODQ(1,KQ),RODQ(2,KQ),RODQ(3,KQ),RODQ(4,KQ)
C         WRITE(LUPRI,*) 'RPQ,EPQ',RPQ,EPQ
C          WRITE(LUPRI,'(A,4F12.8)') 'CLASSIC1',
C     &                  RPQ,EPQ,RODP(4),RODQ(4,KQ)

      END IF

C
C     Local symmetries
C     ================

      CALL LOCSYM(COORPQ,NPPX,1,IPQ0,IPQXYZ)

      RETURN
      END


      SUBROUTINE SCREENAB(IA,IB,CCFA,CCFB,RAB,IPRAB,RODP,FACCLS,KODBC1)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "eriao.h"
#include "ericom.h"
#include "aobtch.h"
#include "odclss.h"
#include "pi.h"
#include "erithr.h"

      DIMENSION CCFA(NPRFA,NCTFA), CCFB(NPRFB,NCTFB),
     &          IPRAB(2,NPPAB),RODP(4)

      IREPA = 0
      IREPB = 0

      AX = CORXBT(IA)
      AY = CORYBT(IA)
      AZ = CORZBT(IA)

      BX = CORXBT(IB)
      BY = CORYBT(IB)
      BZ = CORZBT(IB)

      RAB2 = (AX-BX)**2 + (AY-BY)**2 + (AZ-BZ)**2
      RAB = SQRT(RAB2)

      SM=0.0D0
      AM=0.0D0
      BM=0.0D0
      EXPPIMAX=0.0D0

      IPRFAB = 0

      DO IPRFA=1, NPRFA
       DO IPRFB=1, NPRFB

        KEXPA = KEXPBT(IA) - 1
        KEXPB = KEXPBT(IB) - 1

        EXPA  = EXPBT(KEXPA + IPRFA)
        EXPB  = EXPBT(KEXPB + IPRFB)
        EXPP  = EXPA + EXPB

        EXPPI = 1.0D0 / ( EXPA + EXPB )
        EXPMU = EXPA*EXPB/EXPP

        FAC=0.0D0
        DO ICTFA = 1, NCTFA
         DO ICTFB = 1, NCTFB
          FAC = FAC+ABS(CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB))
         END DO
        END DO

        FAC = R2PI52*EXP(-EXPMU*RAB2)*EXPPI * FAC

        IF ( ABS(FAC) .GT. THRSH ) THEN

         WT=0.0D0
         DO ICTFA = 1, NCTFA
          DO ICTFB = 1, NCTFB
           WT = WT + ABS(CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB))
          END DO
         END DO

         IPRFAB = IPRFAB + 1

         IPRAB(1,IPRFAB) = IPRFA
         IPRAB(2,IPRFAB) = IPRFB

         EXPPIMAX = MAX(EXPPIMAX,EXPPI)

         WT = WT * EXP(-EXPMU*RAB2)
         SM = SM + WT
         AM = AM + WT * EXPA * EXPPI
         BM = BM + WT * EXPB * EXPPI
C          WRITE(LUPRI,'(4F12.8)') WT,SM,AM,BM
C         WRITE(LUPRI,'(A,4F12.8)') 'EXPA,EXPB',EXPA,EXPB,EXPP,EXPPI
C         WRITE(LUPRI,'(A,7F12.8)') 'ScreenAB:',
C     &           FAC,EXP(-EXPMU*RAB2),WT,SM,AM,BM
        END IF

       END DO
      END DO

      AM = AM/SM
      BM = BM/SM

      RODP(1) = AM * AX + BM * BX
      RODP(2) = AM * AY + BM * BY
      RODP(3) = AM * AZ + BM * BZ
      RODP(4) = FACCLS * SQRT(EXPPIMAX)

C      WRITE(LUPRI,'(A,7F12.8)') 'ScreenAB:',AM,BM,SM,(RODP(I),I=1,4)
C
C      WRITE(LUPRI,'(A,4F12.8)') 'RODP: ',
C     &               RODP(1),RODP(2),RODP(3),RODP(4)

      NPRFAB = IPRFAB
      NPPAB  = IPRFAB
      NPP12  = IPRFAB

      RETURN
      END


      SUBROUTINE GETP(COORP,EXPAB,EXPABP,PREFAB,IPRAB,
     &                ETIJ,COORAB,HPIAB,PA,PB,FAC0,
     &                IODBC1,RODBC1,IODPP1,RODPP1,
     &                IELECTR,KODBC1,IPRINT)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "eriao.h"
#include "ericom.h"
#include "aobtch.h"
#include "odclss.h"
#include "pi.h"
#include "erithr.h"
      DIMENSION COORP(3,NPPAB),EXPABP(3,NPPAB),EXPAB(NPPAB,3),
     &          IODBC1(NODBC1,NITBC),RODBC1(NODBC1,NRTBC),
     &          IODPP1(NODPP1,NITPP),RODPP1(NODPP1,NRTPP),
     &          PREFAB(NPPAB),IPRAB(2,NPPAB),
     &          ETIJ(NPPAB,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COORAB(NPPAB,3),
     &          HPIAB(NPPAB), PA(NPPAB,3), PB(NPPAB,3)

      INTEGER X

      EXPPIMAX=0.0D0

      IREPA = 0
      IREPB = 0

      IA = IODBC1(KODBC1,2)
      IB = IODBC1(KODBC1,3)

      AX = CORXBT(IA)
      AY = CORYBT(IA)
      AZ = CORZBT(IA)

      BX = CORXBT(IB)
      BY = CORYBT(IB)
      BZ = CORZBT(IB)

      R2 = (AX-BX)**2 + (AY-BY)**2 + (AZ-BZ)**2

      COORAB(1,1) = AX - BX
      COORAB(1,2) = AY - BY
      COORAB(1,3) = AZ - BZ

      DO 100 I=1, NPPAB
        DO X=1,3
          COORAB(I,X) = COORAB(1,X)
        END DO
 100  CONTINUE

      CALL LOCSYM(COORAB,1,NPPAB,I120,IXYZ)

      IF (IPRINT .GE. 6) THEN

        CALL HEADER('OUTPUT FROM GETP',-1)

        KEXPA = KEXPBT(IA) - 1
C        WRITE (LUPRI,'(5X,A,3I6,F12.8/)')
C     &         'KODBC1, IA, IB: ',KODBC1, IA, IB, EXPBT(KEXPA + 1)

C        WRITE (LUPRI,'(2(5X,A,3F7.3,/),/)')
C     &        'COORA :', AX, AY, AZ,
C     &        'COORB :', BX, BY, BZ

      END IF

      DO IPRFAB=1,NPRFAB

        IPRFA = IPRAB(1,IPRFAB)
        IPRFB = IPRAB(2,IPRFAB)

        KEXPA = KEXPBT(IA) - 1
        KEXPB = KEXPBT(IB) - 1

        EXPA  = EXPBT(KEXPA + IPRFA)
        EXPB  = EXPBT(KEXPB + IPRFB)
        EXPP  = EXPA + EXPB

        EXPPI = 1.0D0 / ( EXPA + EXPB )

        EXPPIMAX = MAX(EXPPIMAX,EXPP)

        EXPABP(1,IPRFAB) = EXPP
        EXPABP(2,IPRFAB) = EXPPI
        EXPABP(3,IPRFAB) = FAC0*EXPPI

        EXPAB(IPRFAB,1) = EXPA * EXPPI
        EXPAB(IPRFAB,2) = EXPB * EXPPI
        EXPAB(IPRFAB,3) = EXPPI

        DO 150 X = 1,3
          PA(IPRFAB,X) = - EXPAB(IPRFAB,2) * COORAB(IPRFAB,X)
          PB(IPRFAB,X) =   EXPAB(IPRFAB,1) * COORAB(IPRFAB,X)
 150    CONTINUE

        PREFAB(IPRFAB) = exp( -EXPA*EXPB*EXPPI*R2 )

C        WRITE (LUPRI,'(5X,A,3F14.6,/)')
C     &         'EXPA, EXPB, EXPAB:',
C     &          EXPA, EXPB, EXPP

        COORP(1,IPRFAB) = (EXPA*AX+EXPB*BX)*EXPPI
        COORP(2,IPRFAB) = (EXPA*AY+EXPB*BY)*EXPPI
        COORP(3,IPRFAB) = (EXPA*AZ+EXPB*BZ)*EXPPI

C         WRITE(LUPRI,'(A,2I3,4F12.6)') '@KODBC1:',KODBC1,IPRFAB,EXPPI,
C     &           COORP(1,IPRFAB),COORP(2,IPRFAB),COORP(3,IPRFAB)

C        WRITE (LUPRI,'(5X,I4,3X,A,3F7.3,/)') IPRFAB,
C     &        'COORP :', ( COORP(X,IPRFAB),X=1,3 )

      END DO

C
C  IRUTIN - here just to print it
C
      IRUTIN = 1

      ISCAL1 = 0
      ISCAL2 = 0
C
C  level of derivatives calculated
C
      IDL    = 0
      NCNTAB = 0

      CALL EXPCFT(ETIJ,NPPAB,JMAX1,JMAX2,COORAB,EXPAB,I120,HPIAB,
     &                  PA,PB,IRUTIN,IELECTR,JMAX1,JMAX2,NCNTAB,IDL,
     &                  ISCAL1,ISCAL2,IPRINT)


C      DO IPRFAB=1,NPRFAB
C       WRITE(LUPRI,*) IPRAB(1,IPRFAB),IPRAB(2,IPRFAB)
C      END DO

      RETURN
      END

      SUBROUTINE HRFTUV( HRINT, HRINT1, FTUV, NDMT,
     &                  INDHER, IODDHR, IPRINT )
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "ericom.h"
#include "hertop.h"

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV
      DIMENSION HRINT(NPPAB,NTUVAB,NDMT),HRINT1(NPPX,NTUV),
     &       FTUV(NPPCD,0:JMAXCD,0:JMAXCD,0:JMAXCD,NDMT),
     &       INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(NRTOP)

      INTEGER X

C      CALL HEADER('HRINT in HRFTUV',-1)

      DO 600 IDMT = 1, NDMT
       DO 500 JAB = 0, JMAXAB
        DO 400 TAB = JAB, 0, -1
         DO 300 UAB = JAB - TAB, 0, -1
          VAB = JAB - TAB - UAB
          TUVAB = INDHER(TAB,UAB,VAB)
C          WRITE(LUPRI,*) 'TUVAB: ',TUVAB
C          CALL OUTPUT(HRINT(1,TUVAB,IDMT),1,NPPAB,1,1,
C     &                NPPAB,1,1,LUPRI)
           DO JCD=0,JMAXCD
            DO TCD = JCD, 0, -1
             DO UCD = JCD - TCD, 0, -1
              VCD = JCD - TCD - UCD
              TUV = INDHER(TAB+TCD,UAB+UCD,VAB+VCD)

              IF (IODDHR(TUV) .EQ. 0) THEN
C               CALL DGEMV('n',NPPAB,NPPCD,1.0D0,
C     &                   HRINT1(1,TUV),NPPCD,
C     &                   FTUV(1,TCD,UCD,VCD,IDMT),1,
C     &                   HRINT(1,TUVAB,IDMT),1)

C              IPPX = 1
C               DO 200 IPPCD = 1, NPPCD
C                DO 100 IPPAB = 1, NPPAB
C                 HRINT(IPPAB,TUVAB,IDMT) =
C     &             HRINT(IPPAB,TUVAB,IDMT) + HRINT1(IPPX,TUV)*
C     &                          FTUV(IPPCD,TCD,UCD,VCD,IDMT)
C                 IPPX = IPPX + 1
C 100            CONTINUE
C 200           CONTINUE

               IPPX=1
               DO 200 IPPAB = 1, NPPAB
                FSUM=0.0D0
                DO 100 IPPCD = 1, NPPCD
                 FSUM = FSUM +
     &                  HRINT1(IPPX,TUV)*FTUV(IPPCD,TCD,UCD,VCD,IDMT)
                 IPPX = IPPX+1
 100            CONTINUE
                 HRINT(IPPAB,TUVAB,IDMT)=HRINT(IPPAB,TUVAB,IDMT)+FSUM
 200           CONTINUE
              END IF
             END DO
            END DO
           END DO
 300     CONTINUE
 400    CONTINUE
 500   CONTINUE
 600  CONTINUE

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
      SUBROUTINE ADDFCK3(FMAT,FMATCP,FMATSP,FTMP,
     &             NDMT,ETIJ,PREFAB,IPRAB,HRINT,INDHER,LMNPWR,
     &             IODDHR,
     &             CSQA,CSQB,CCFA,CCFB,INDXBT,NPPAB,
     &             JMAXA, JMAXB, JMAXAB, NTUV,
     &             KHKTA, KHKTB, KCKTA, KCKTB,
     &             NPRFA, NPRFB, NCTFA, NCTFB,
     &             KCKMAX, NHKMAX,IA,IB,IPRINT )
C
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "aobtch.h"
#include "inforb.h"
#include "hertop.h"
#include "symmet.h"
C
      INTEGER T,U,V,TUV
      DIMENSION FMAT(NBAST,NBAST,NDMT),
     &          FMATCP(KCKTA,KCKTB,NPPAB),
     &          FMATSP(KHKTA,KHKTB,NPPAB),
     &          FTMP(KCKTA,KHKTB),
     &   ETIJ(NPPAB,0:JMAXAB,0:JMAXA,0:JMAXB,3,*),
     &   PREFAB(NPPAB),HRINT(NPPAB,NTUV,NDMT),
     &   IPRAB(2,NPPAB),
     &   INDHER(0:JTOP,0:JTOP,0:JTOP),
     &   IODDHR(NRTOP),
     &   LMNPWR(KCKMAX,NHKMAX,3),
     &   CSQA(KHKTA,KCKTA), CSQB(KHKTB,KCKTB),
     &   CCFA(NPRFA,NCTFA), CCFB(NPRFB,NCTFB),
     &   INDXBT(MXSHEL*MXCONT,0:7)

      IF ( IPRINT .GE. 10 ) THEN
        CALL HEADER('Output from ADDFCK3',-1)
      END IF

      CALL DZERO(FMATCP, KCKTA*KCKTB*NPPAB )
      CALL DZERO(FMATSP, KHKTA*KHKTB*NPPAB )
      CALL DZERO(FTMP, KCKTA*KHKTB )

      IDMT = 1

      IREPA  = 0
      IREPB  = 0

      IF ( IPRINT .GE. 15 ) THEN

       DO JAB = 0,JMAXAB
        DO T = 0, JAB
         DO U = 0, JAB - T
          V = JAB - T - U
          TUV = INDHER(T,U,V)
          DO IPPAB=1,NPPAB
           WRITE(LUPRI,*) 'JAB:',JAB, T,U,V,TUV,IPPAB,
     &                    HRINT(IPPAB,TUV,IDMT)
          END DO
         END DO
        END DO
       END DO
      END IF

      DO 500 IPRFAB = 1, NPPAB

       IPRFA = IPRAB(1,IPRFAB)
       IPRFB = IPRAB(2,IPRFAB)

        DO JAB = 0,JMAXAB
        DO T = 0, JAB
         DO U = 0, JAB - T
          V = JAB - T - U

           TUV = INDHER(T,U,V)

C           IF ( IODDHR(TUV) .EQ. 0) THEN
            DO KCA = 1, KCKTA

            LA = LMNPWR(KCA,JMAXA+1,1)
            MA = LMNPWR(KCA,JMAXA+1,2)
            NA = LMNPWR(KCA,JMAXA+1,3)

C            WRITE(LUPRI,*) 'A:', KCA, LA,MA,NA

            DO KCB = 1, KCKTB

             LB = LMNPWR(KCB,JMAXB+1,1)
             MB = LMNPWR(KCB,JMAXB+1,2)
             NB = LMNPWR(KCB,JMAXB+1,3)

C             WRITE(LUPRI,*) 'B:', KCB, LB,MB,NB
C             WRITE(LUPRI,*) 'Etij:',
C     &                 ETIJ(IPRFAB,T,LA,LB,1,1),
C     &                 ETIJ(IPRFAB,U,MA,MB,2,1),
C     &                 ETIJ(IPRFAB,V,NA,NB,3,1)


             FAC = ETIJ(IPRFAB,T,LA,LB,1,1)*
     &         ETIJ(IPRFAB,U,MA,MB,2,1)*ETIJ(IPRFAB,V,NA,NB,3,1)

             FMATCP(KCA,KCB,IPRFAB) =
     &            FMATCP(KCA,KCB,IPRFAB) + FAC*PREFAB(IPRFAB)
     &                      * HRINT(IPRFAB,TUV,IDMT)
C             WRITE(LUPRI,*) 'FMATCP:',FMATCP(KCA,KCB,IPRFAB),
C     &                      'T,U,V,TUV', T,U,V,TUV

           END DO
          END DO

C          END IF
          END DO
         END DO
        END DO

         CALL DGEMM('N','T', KCKTA, KHKTB, KCKTB, 1.0D0,
     &              FMATCP(1,1,IPRFAB), KCKTA, CSQB, KHKTB,
     &              0.0D0, FTMP, KCKTA  )
        CALL DGEMM('N','N', KHKTA, KHKTB, KCKTA, 1.0D0,
     &              CSQA, KHKTA, FTMP, KCKTA,
     &              0.0D0, FMATSP(1,1,IPRFAB), KHKTA  )

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


        DO 300 ICTFA = 1, NCTFA
         DO 200 ICTFB = 1, NCTFB

          IAOA = INDXBT(KNDXBT(IA) - 1 + ICTFA,IREPA)
          IAOB = INDXBT(KNDXBT(IB) - 1 + ICTFB,IREPB)

C           WRITE( LUPRI, * ) 'IA, IB, IAOA, IAOB: ',
C     &                       IA, IB, IAOA, IAOB

C
C Here it we should probably put summation over
C  FMATSP(ISPA,ISPB,IPRFAB) separately outside
C

          IF ( (IAOA .EQ. IAOB ) .OR. ( IA .EQ. IB )) THEN
           DO ISPA=1, KHKTA
            DO ISPB=1, KHKTB
              FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMT) =
     &            FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMT)
     &          + FMATSP(ISPA,ISPB,IPRFAB)
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
            END DO
           END DO
          ELSE
           DO ISPA=1, KHKTA
            DO ISPB=1, KHKTB
              FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMT) =
     &            FMAT(IAOA+ISPA-1,IAOB+ISPB-1,IDMT)
     &          + FMATSP(ISPA,ISPB,IPRFAB)
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
              FMAT(IAOB+ISPB-1,IAOA+ISPA-1,IDMT) =
     &            FMAT(IAOB+ISPB-1,IAOA+ISPA-1,IDMT)
     &          + FMATSP(ISPA,ISPB,IPRFAB)
     &          * CCFA(IPRFA,ICTFA)*CCFB(IPRFB,ICTFB)
            END DO
           END DO
          END IF


 200     CONTINUE
 300    CONTINUE

 500  CONTINUE

      RETURN
      END


      SUBROUTINE PRITMP(HRINT1,COORPQ,EXPAB,EXPCD,ALPHA,
     &                  INDHER,IODDHR,KODBC1, KODBC2)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "ericom.h"
#include "aobtch.h"
#include "odclss.h"
#include "hertop.h"

      INTEGER T,U,V,TUV
      DIMENSION HRINT1(NPPX,NTUV),COORPQ(NPPX,3),
     &          EXPAB(NPRFAB,3),EXPCD(NPRFCD,3),ALPHA(NPPX),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),IODDHR(NRTOP)

      IPPX=0

      JMAX = JMAXA + JMAXB + JMAXC + JMAXD

      CALL HEADER('PRITMP',-1)

      WRITE(LUPRI,*) ' KODBC1: ' ,KODBC1,
     &               ' KODBC2: ' ,KODBC2

      WRITE(LUPRI,*) ' NPRFA, NPRFB, JMAXA, JMAXB ',
     &                 NPRFA, NPRFB, JMAXA, JMAXB

      WRITE(LUPRI,*) ' NPRFC, NPRFD, JMAXC, JMAXD ',
     &                 NPRFC, NPRFD, JMAXC, JMAXD
      WRITE(LUPRI,*)

      DO IPRFCD=1,NPPCD
       DO IPRFAB=1,NPPAB

        IPPX=IPPX+1

        WRITE (LUPRI,'(5X,A,3I5,/,5X,A,3F12.8)')
     &     'IPRFAB, IPRFCD, IPPX :', IPRFAB, IPRFCD, IPPX,
     &     'COORPQ: ', (COORPQ(IPPX,IX),IX=1,3)

        EXPP = 1.0D0/ABS(EXPAB(IPRFAB,3))
        EXPQ = 1.0D0/ABS(EXPCD(IPRFCD,3))

        WRITE (LUPRI,'(5X,A,4F12.8,/)')
     &         'EXPP, EXPQ, ALPHA, ALPHA_DEF : ',
     &         EXPP,EXPQ,ALPHA(IPPX),EXPP*EXPQ/(EXPP+EXPQ)

        DO J=0,JMAX
         DO T = J, 0, -1
          DO U = J - T, 0, -1
           V = J - T - U
           TUV = INDHER(T,U,V)
C           IF (IODDHR(TUV) .EQ. 0) THEN
             WRITE (LUPRI,'(4X,3(A,I1),A,I4,F12.6)')
     &         'R(',T,',',U,',',V,')',TUV,HRINT1(IPPX,TUV)
C           END IF
          END DO
         END DO
        END DO

        WRITE(LUPRI,*)

       END DO
      END DO

      RETURN
      END

      SUBROUTINE HRFSPEC00(HRINT,RJ000,FTUV,
     &                     NDMT,NPPAB,NPPCD,NPPX)
#include "implicit.h"
#include "priunit.h"
      DIMENSION HRINT(NPPAB,NDMT),RJ000(NPPX),FTUV(NPPCD,NDMT)


      DO 300 IDMT=1,NDMT
      IPPX=1
       DO 200 IPPAB = 1, NPPAB
        FSUM=0.0D0
        DO IPPCD = 1, NPPCD
          FSUM=FSUM+RJ000(IPPX)*FTUV(IPPCD,IDMT)
          IPPX = IPPX + 1
        END DO
        HRINT(IPPAB,IDMT) =  HRINT(IPPAB,IDMT) + FSUM
 200   CONTINUE
 300  CONTINUE

      RETURN
      END


      SUBROUTINE HRFSPEC(HRINT,RJ000,RPQ,FTUV,
     &                   INDHER,IODDHR,NDMT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "hertop.h"
#include "ericom.h"
      DIMENSION HRINT(NPPAB,NTUVAB,NDMT),HRINT1(NPPX,NTUV),
     &          RJ000(NPPX,0:JMAX), RPQ(NPPX,3),
     &          FTUV(NPPCD,0:JMAXCD,0:JMAXCD,0:JMAXCD,NDMT),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), IODDHR(NRTOP)

      INTEGER TAB,UAB,VAB,TUVAB,TCD,UCD,VCD,TUVCD,TUV

      DO 1000 IDMT=1,NDMT

      IF (JMAX .EQ.0) THEN
      IPPX=1
       DO 200 IPPAB = 1, NPPAB
        FSUM=0.0D0
        DO IPPCD = 1, NPPCD
          FSUM=FSUM+RJ000(IPPX,0)*FTUV(IPPCD,0,0,0,IDMT)
          IPPX = IPPX + 1
        END DO
        HRINT(IPPAB,1,IDMT) =  HRINT(IPPAB,1,IDMT) + FSUM
 200   CONTINUE

      ELSE

      IF ( (JMAXAB .EQ. 0) .AND. (JMAXCD.EQ.1) ) THEN
      CALL HEADER('Special case JMAXAB=0, JMAXF=1',-1)
      IX=1
       DO IP = 1, NPPAB
        F000=0.0D0
        F100=0.0D0
        F010=0.0D0
        F001=0.0D0
        DO IQ = 1, NPPCD
          F000=F000+          RJ000(IX,0)*FTUV(IQ,0,0,0,IDMT)
          F100=F100+RPQ(IX,1)*RJ000(IX,1)*FTUV(IQ,1,0,0,IDMT)
          F010=F010+RPQ(IX,2)*RJ000(IX,1)*FTUV(IQ,0,1,0,IDMT)
          F001=F001+RPQ(IX,3)*RJ000(IX,1)*FTUV(IQ,0,0,1,IDMT)
          WRITE(LUPRI,'(3X,3I5,6F12.8)') IX,IP,IQ,
     &       RJ000(IX,0),RJ000(IX,1),F000,F100,F010,F001
          IX = IX + 1
        END DO
        HRINT(IP,1,IDMT) =  HRINT(IP,1,IDMT) + F000
        HRINT(IP,2,IDMT) =  HRINT(IP,2,IDMT) - F100
        HRINT(IP,3,IDMT) =  HRINT(IP,3,IDMT) - F010
        HRINT(IP,4,IDMT) =  HRINT(IP,4,IDMT) - F001

       END DO

      ELSE

      IF ( (JMAXAB .EQ. 1) .AND. (JMAXCD.EQ.0) ) THEN
      IX=1
       DO IP = 1, NPPAB
        F000=0.0D0
        F100=0.0D0
        F010=0.0D0
        F001=0.0D0
        DO IQ = 1, NPPCD
          F000=F000+          RJ000(IX,0)*FTUV(IQ,0,0,0,IDMT)
          F100=F100+RPQ(IX,1)*RJ000(IX,1)*FTUV(IQ,0,0,0,IDMT)
          F010=F010+RPQ(IX,2)*RJ000(IX,1)*FTUV(IQ,0,0,0,IDMT)
          F001=F001+RPQ(IX,3)*RJ000(IX,1)*FTUV(IQ,0,0,0,IDMT)
          IX = IX + 1
        END DO
        HRINT(IP,1,IDMT) =  HRINT(IP,1,IDMT) + F000
        HRINT(IP,2,IDMT) =  HRINT(IP,2,IDMT) - F100
        HRINT(IP,3,IDMT) =  HRINT(IP,3,IDMT) - F010
        HRINT(IP,4,IDMT) =  HRINT(IP,4,IDMT) - F001

       END DO

      END IF

      END IF

      END IF

 1000 CONTINUE

      RETURN
      END

c ===================================================================
c
c Parallelization part of the fck3 code.
c The variables to be broadcast are split into three groups:
c 1. variables to be broadcast on each run.
c 2. basis-set dependent variables (synced when NEWBAS).
c 3. variables to be synced only once (FIRST), as far as we can tell.
c
c fck3col collects the data back to master.
c FCK3_SLAVE_FOCK is the code executed by slaves.
c ===================================================================
#if defined (VAR_MPI)
      SUBROUTINE FCK3_SLAVE_FOCK(WRK,LWORK,IPRINT)
#include "implicit.h"
      DIMENSION WRK(LWORK)
#include "aovec.h"
#include "inforb.h"
#include "maxorb.h"
#include "primit.h"
#include "priunit.h"
#include "shells.h"
      INCLUDE "mpif.h"
#include "infpar.h"
c
#include "memint.h"
c
      CALL fck3_dims(NDMAT)
c     initialize some common blocks that are forgotten
      CALL carpow
      NSYMOP = 0
      CALL baspar(NSYMOP)
      Call MEMGET('REAL',KFMAT, NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT, NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT, NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSAO,  N2BASX,      WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCCFBT,MXPRIM*MXCONT,  WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDXBT,8*MXSHEL*MXCONT,WRK,KFREE,LFREE)
      CALL FCK3SN1(WRK(KDMAT),NDMAT,ISYMDM,IFCTYP,WRK(KSAO))
      CALL HRFCK3(WRK(KFMAT),WRK(KDMAT),NDMAT,ISYMDM,IFCTYP,IPRINT,
     &            WRK(KSAO),WRK(KCCFBT),WRK(KINDXBT),WRK(KFREE),LFREE)
c      CALL MEMREL('FCK3_SLAVE',WRK,KWORK,KWORK,KFREE,LFREE)
#if 0
      CALL HEADER('Density and Fock matrices from SLAVE',-1)
      KSTR = 0
      DO I = 1, NDMAT
         WRITE (LUPRI,'(//,1X,I2,A,I3)') MYNUM,' Density matrix No.',I
         CALL OUTPUT(WRK(KDMAT+KSTR),1,NBAST,1,NBAST,NBAST,
     &               NBAST,1,LUPRI)
         WRITE (LUPRI,'(//,1X,I2,A,I3)') MYNUM,' Fock matrix No.',I
         CALL OUTPUT(WRK(KFMAT+KSTR),1,NBAST,1,NBAST,NBAST,
     &               NBASE,1,LUPRI)
         KSTR = KSTR + N2BASX
      END DO
#endif
      END
c     fck2sn1 - broadcast variables each time.
      SUBROUTINE fck3sn1(DMAT,NDMT,ISYMDM,IFCTYP,SAO)
#include "implicit.h"
#include "inforb.h"
#include "siripc.h"
#if defined (PRG_DIRAC) /* PRG_DIRAC */
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif /* PRG_DIRAC */
c note: NDMT is always INTENT(IN)
      DIMENSION DMAT(N2BASX,NDMT), ISYMDM(NDMT),IFCTYP(NDMT),
     &          SAO(N2BASX)
      INCLUDE 'mpif.h'
#include "maxorb.h"
#include "infpar.h"
c
      CALL MPI_BCAST(DMAT,N2BASX*NDMT,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(SAO,N2BASX,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISYMDM,NDMT,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IFCTYP,NDMT,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOFMM,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      IF(NEWBAS) CALL fck3sn2
      IF(NEWGEO) CALL fck3sn3
      END
      SUBROUTINE fck3sn2
      END
      SUBROUTINE fck3sn3
      END
#if defined (FCK3COLLECT_NO_CONFLICT)
      SUBROUTINE fck3collect(FMAT,NDMAT,tm,WORK,LFREE,IPRINT)
#include "implicit.h"
#include "inforb.h"
      DIMENSION FMAT(N2BASX,NDMAT), WORK(*)
      INCLUDE 'mpif.h'
#include "maxorb.h"
#include "infpar.h"
#include "priunit.h"
c
      DIMENSION tms(NODTOT+1)
      NSIZE = N2BASX*NDMAT
      IF (LFREE .LT. NSIZE) THEN
         WRITE(LUPRI,"(I9,'*',I2,'=',I9,' needed, ',I9,' available.')")
     &        N2BASX,NDMAT,N2BASX*NDMAT,LFREE
        Call QUIT('Insufficient memory in fck3collect')
      END IF
      CALL DCOPY(NSIZE,FMAT,1,WORK,1)
      CALL MPI_REDUCE(WORK,FMAT,NSIZE,MPI_DOUBLE_PRECISION,MPI_SUM,
     &     MASTER,MPI_COMM_WORLD,IERR)
      iF(iprint.ge.0) THEN
c     for a good measure, collect times to to provide more info about
c     load balancing.
         CALL MPI_GATHER(tm,1,MPI_DOUBLE_PRECISION,
     &        tms,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,IERR)
         if(IERR.NE.0) WRITE(LUPRI,*) 'error collecting times:',IERR
         WRITE(LUPRI,"('Node times:',/,20(10F7.2,/))")
     &        (tms(i),i=1,NODTOT+1)
      END IF
      END
#endif /* FCK3COLLECT_NO_CONFLICT */
c
      SUBROUTINE FCK3_DIMS(NDMAT)
#include "implicit.h"
#include "erithr.h"
      INCLUDE 'mpif.h'
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "odclss.h"
#include "aovec.h"
#include "primit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "shells.h"
#include "sphtrm.h"
#include "nuclei.h"
#include "hertop.h"
#include "densfit.h"
c FIXME: ccom for NHTYP and SPHNRM which should really go away.
#include "ccom.h"
c     for ODBCH
      CALL MPI_BCAST(KMAX,1,my_MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NDMAT,1,my_MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBAST,1,my_MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(N2BASX,1,my_MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NODCL1,1,my_MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NODCL2,1,my_MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
c     for DENS_FIT
      CALL MPI_BCAST(KMAXTOT,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
c     for SETAOS
      KMAXX = KMAX
      IF (DOFIT) KMAXX = KMAXTOT
c     CALL MPI_BCAST(ISTBAO,KMAXX,my_MPI_INTEGER,MASTER,
c    &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MULT,8,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(JSTRT,KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NCENT,KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBCH, KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NUCO, KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NRCO, KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(KCKT, KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(KHKT, KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NHKT, KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LCLASS,KMAXX,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NHTYP,1 ,  my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(SPHNRM,1,  MPI_LOGICAL,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CENT, MXSHEL*6,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CORD, MXCENT*3,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PRIEXP,MXPRIM,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PRICCF, MXPRIM*MXCONT,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(THRSH,1,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PT,8,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
c     shptrm: FIXME: call sphinp instead
      CALL MPI_BCAST(ISPADR,MXQN,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CSP,NCSP,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
c     nuclei
      CALL MPI_BCAST(NATOMS,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
c     needed for HERPRP
      CALL MPI_BCAST(JTOP,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NRTOP,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      END
#endif
