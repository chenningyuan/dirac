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

C  /* Deck eriexc */
CC#ifndef PRG_DIRAC
      SUBROUTINE ERIEXC(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     *                  CCFBT,INDXBT,FCKTMP,WORK,LWORK)
CC#else
CC      SUBROUTINE ERIEXC(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
CC     *                  CCFBT,INDXBT,FCKTMP,
CC     *                  NCENTS,ICENTS,ICLASSES,WORK,LWORK)
CC#endif
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
#include "ccom.h"
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "erimem.h"
#include "aobtch.h"
#include "odbtch.h"
#include "symmet.h"
#include "infpar.h"
#include "inforb.h"
#include "shells.h"
#include "cbihr2.h"

#ifdef PRG_DIRAC
! HFXFAC in dcbham.h
#include "dcbham.h"
#endif

      DIMENSION FMAT(*), DMAT(*), WORK(LWORK),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*), 
     &	        FCKTMP(*), IFCSAV(NDMT)
C       
C and now comes the dimension statement for the class selection
C used in dirac
C
CC      DIMENSION NCENTS(4), ICENTS(2,4),ICLASSES(4), KOD(10)
      DIMENSION KOD(10)
      LOGICAL DOEXC
#include "ibtfun.h"
#include "memint.h"
C
      IF (SLAVE) THEN
         IPRINT = IPRFCK
      ELSE
         CALL TIMER('START ',TIMSTR,TIMEND)
C
C        Initialization in ER2INI
C
         CALL ER2INI
C
         IPRINT = MAX(IPRERI,IPRFCK)
      END IF
C     IF (NDMT.GT.1) CALL QUIT('ERIEXI works only with NDMT.eq.1')
C        
      DO I=1,NDMT
        IF(MOD(IFCTYP(I),10).LT.2) THEN
            CALL QUIT('Routine entered but no exchange requested')
        ENDIF
c	whatever was requested, we will compute exchange only
        IFCSAV(I) = IFCTYP(I)-MOD(IFCTYP(I),10)+2
      END DO
C 
      IF (MAXREP.GT.0) CALL QUIT('ERIEXC does not work with symmetry')
C
cSimen Multiplication by empirical factor 0.1 to get approximately the same
cSimen resuls as for hermit
      THRXCH = SCRFCK * 1.00D-1
      THRGAM = THRXCH * 1.00D-3
      THRSH  = THRXCH * 1.00D-5
      IF (THRXC1.NE.D0) THRXCH = THRXC1
      IF (THRXC2.NE.D0) THRGAM = THRXC2
      IF (THRXC3.NE.D0) THRSH  = THRXC3
      WRITE (LUPRI,'(2X,A,3(2X,E10.2))') 'Thresholds for exchange: ',
     &      THRXCH, THRGAM, THRSH
      NDMAT  = NDMT
C
      IF (HFXMU.GT.D0) ATTEXC = .TRUE. 
C
      WRTINT = .FALSE.
      FCKINT = .TRUE.
      EXCINT = .FALSE.
C
#ifndef PRG_DIRAC
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
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
      EXCINT = .TRUE.
      CALL ODBCHS(KODCL1,KODCL2,
     &            KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
      EXCINT = .FALSE.
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C
#else
C aspg, 2006-04-21
C
C eriprp_ is a function created to set-up the ao batches for dirac
C so wonder if a general case could be written and incorporated into
C the original eri code seen above in the ifndef clause
C
      DOEXC = .TRUE.
C
C aspg, 2006-06-27 making ftnchek happy
C
      LFREE = LWORK
C
      CALL ERIPRP(DOEXC,CCFBT,INDXBT,KOD,KMAX,IPRINT,WORK,KFREE,LFREE)
C
#endif /* ! PRG_DIRAC */
      ICALL = 0
      CALL GETDST(ICALL,ICALL,IPRINT)
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
      KGAB    = KFREE
      KDAB    = KGAB + NAOBCH**2
      KDNSCNT = KDAB + NAOBCH**2
      KDNSIDX = KDNSCNT + (NAOBCH      + 1)/IRAT + 1
      KSHLCNT = KDNSIDX + (NAOBCH**2   + 1)/IRAT + 1
      KSHLIDX = KSHLCNT + (NAOBCH      + 1)/IRAT + 1
      KCLSS   = KSHLIDX + (NAOBCH**2   + 1)/IRAT + 1
      KLAST   = KCLSS   + (2*NAOBCH**2 + 1)/IRAT + 1
      LWRK    = LFREE   - KLAST + 1
      IF (KLAST .GT. LFREE) CALL STOPIT('ERIEXC','RED',KLAST,LWRK)
#ifndef PRG_DIRAC
      CALL EXIRED(DMAT,NDMT,WORK(KDAB),WORK(KGAB),
     &            WORK(KDNSCNT),WORK(KDNSIDX),
     &            WORK(KSHLCNT),WORK(KSHLIDX),WORK(KCLSS),
     &            WORK(KODCL1),WORK(KODCL2),
     &            WORK(KODBC1),WORK(KODBC2),
     &            WORK(KRDBC1),WORK(KRDBC2),
     &            WORK(KODPP1),WORK(KODPP2),
     &            WORK(KRDPP1),WORK(KRDPP2),
     &            CCFBT,INDXBT,IPRINT,
     &            WORK(KLAST),LWRK)
#else
      CALL EXIRED(DMAT,NDMT,WORK(KDAB),WORK(KGAB),
     &            WORK(KDNSCNT),WORK(KDNSIDX),
     &            WORK(KSHLCNT),WORK(KSHLIDX),WORK(KCLSS),
     &            WORK(KOD(1)),WORK(KOD(2)),
     &            WORK(KOD(3)),WORK(KOD(4)),
     &            WORK(KOD(5)),WORK(KOD(6)),
     &            WORK(KOD(7)),WORK(KOD(8)),
     &            WORK(KOD(9)),WORK(KOD(10)),
     &            CCFBT,INDXBT,IPRINT,
     &            WORK(KLAST),LWRK)
#endif /* ! PRG_DIRAC */
C
C     Calculate integrals
C     ===================
C
      EXCINT = .TRUE.
      IF (SLAVE) THEN
#ifndef PRG_DIRAC
         CALL EXIDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               WORK(KDAB),  WORK(KGAB),
     &               WORK(KDNSCNT),WORK(KDNSIDX),
     &               WORK(KSHLCNT),WORK(KSHLIDX),WORK(KCLSS),
     &               FMAT,DMAT,NDMT,IFCSAV,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK(KLAST),LWRK,IPRINT)
#else /* ! PRG_DIRAC */
         CALL EXIDRV(WORK(KOD(1)),WORK(KOD(2)),
     &               WORK(KOD(3)),WORK(KOD(4)),
     &               WORK(KOD(5)),WORK(KOD(6)),
     &               WORK(KOD(7)),WORK(KOD(8)),
     &               WORK(KOD(9)),WORK(KOD(10)),
     &               WORK(KDAB),  WORK(KGAB),
     &               WORK(KDNSCNT),WORK(KDNSIDX),
     &               WORK(KSHLCNT),WORK(KSHLIDX),WORK(KCLSS),
     &               FMAT,DMAT,NDMT,IFCSAV,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK(KLAST),LWRK,IPRINT)
#endif /* ! PRG_DIRAC */
      ELSE
         IF (.NOT.INTSKP) THEN
#ifndef PRG_DIRAC
            CALL EXIDRV(WORK(KODCL1),WORK(KODCL2),
     &                  WORK(KODBC1),WORK(KODBC2),
     &                  WORK(KRDBC1),WORK(KRDBC2),
     &                  WORK(KODPP1),WORK(KODPP2),
     &                  WORK(KRDPP1),WORK(KRDPP2),
     &                  WORK(KDAB),  WORK(KGAB),
     &                  WORK(KDNSCNT),WORK(KDNSIDX),
     &                  WORK(KSHLCNT),WORK(KSHLIDX),WORK(KCLSS),
     &                  FMAT,DMAT,NDMT,IFCSAV,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK(KLAST),LWRK,IPRINT)
#else /* ! PRG_DIRAC */
            CALL EXIDRV(WORK(KOD(1)),WORK(KOD(2)),
     &                  WORK(KOD(3)),WORK(KOD(4)),
     &                  WORK(KOD(5)),WORK(KOD(6)),
     &                  WORK(KOD(7)),WORK(KOD(8)),
     &                  WORK(KOD(9)),WORK(KOD(10)),
     &                  WORK(KDAB),  WORK(KGAB),
     &                  WORK(KDNSCNT),WORK(KDNSIDX),
     &                  WORK(KSHLCNT),WORK(KSHLIDX),WORK(KCLSS),
     &                  FMAT,DMAT,NDMT,IFCSAV,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK(KLAST),LWRK,IPRINT)
#endif /* ! PRG_DIRAC */
         END IF
         CALL EXISYM(FMAT,IFCTYP,NDMT,NBASE)
C
C        Print densities and Fock matrix
C        ===============================
C
         IF (IPRINT.GT.4) THEN
            CALL HEADER('Density and Fock matrices in ERIFCK',-1)
            KSTR = 1 
            DO I = 1, NDMT
               WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
               CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
               CALL OUTPUT(FMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               KSTR = KSTR + NBASE*NBASE
            END DO
         END IF
C
         CALL TIMER('XBUILD',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
      END IF
      EXCINT = .FALSE.
      ATTEXC = .FALSE. 
C
      RETURN
      END
C  /* Deck exisym */
      SUBROUTINE EXISYM(FMAT,IFCTYP,NDMT,NBASE)
#include "implicit.h"
      PARAMETER (DP5 = 0.5D0)
      DIMENSION FMAT(NBASE,NBASE,NDMT), IFCTYP(NDMT)
C
      DO K = 1, NDMT
#if defined (PRG_DIRAC)
Caspg, modifications to allow anti-symmetric matrices
         IH = IFCTYP(K)/10
         IF (IH.EQ.1) THEN
Csymmetrix matrices 
#else
         IF (IFCTYP(K).GT.10) THEN
#endif
            DO J = 1, NBASE
            DO I = 1, J - 1 
               FIJ = DP5*(FMAT(I,J,K) + FMAT(J,I,K))
               FMAT(I,J,K) = FIJ 
               FMAT(J,I,K) = FIJ 
            END DO
            END DO
#if defined (PRG_DIRAC)
Canti-symmetric matrices
         ELSEIF (IH.EQ.2) THEN
              DO J = 1, NBASE
               DO I = 1, J - 1 
                  FIJ = DP5*(FMAT(I,J,K) - FMAT(J,I,K))
                  FMAT(I,J,K) =  FIJ 
                  FMAT(J,I,K) = -FIJ 
               END DO
               END DO
#endif
          END IF
      END DO
      RETURN
      END
C  /* Deck exidrv */
      SUBROUTINE EXIDRV(IODCL1,IODCL2,IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  DRED,GRED,IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,
     &                  FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                  CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER(MXNODT=30)
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          FMAT(*), DMAT(*), IFCTYP(*), D2MAT(*), ID2MAT(*),
     &          CCFBT(*), INDXBT(*), WORK(LWORK)
      DIMENSION DRED(*), GRED(*),
     &          IDNSCNT(*), IDNSIDX(*), ISHLCNT(*), ISHLIDX(*), ICLSS(*)
#include "cbieri.h"
#include "odclss.h"
#include "ccom.h"
#include "ericom.h"
#include "eridst.h"
#include "r12int.h"
#include "hertop.h"
#include "eribuf.h"
#include "nuclei.h"
#include "aobtch.h"
#include "iratdef.h"
#include "infpar.h"
C
      MAXDER = 0
      JTOP  = 4*(NHTYP - 1) + MAXDER
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
      NHKMAX = NHTYP + MAXDER
      KCKMAX = NHKMAX*(NHKMAX + 1)/2
      KC2MAX = KCKMAX**2
      NRDER  = 3*MAXDER
C
C     Buffers and bits
C     ================
C
      LBUF  = LBFINP
      NBUFS = NDISTR
      NPERTS = 1
      NCORS = 1
#if defined (SYS_CRAY) || defined (SYS_T90) || defined (INT_STAR8) || defined (SYS_NEC) 
      NIBUF = 1
      IF (NBASE .LE. 255) THEN
         NBITS = 8
      ELSE
         NBITS = 16
      END IF
#else
      IF (NBASE .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IF (DGABAB) THEN
         NIBUF = 1
         NBITS = 16
      END IF
C
      MAXBCH = 50
C
C     Dimensions for CSQ
C     ==================
C
      CALL SPHDIM
C
C     Allocate
C     ========
C
      KACTAB = 1
      KACTCD = KACTAB + (NODCL1 - 1)/IRAT + 1
      KSCRAB = KACTCD + (NODCL2 - 1)/IRAT + 1
      KSCRCD = KSCRAB + NODBC1
      KJDBC2 = KSCRCD + NODBC2
      KJDCL2 = KJDBC2 + (MXNODT*NAOBCH**2 - 1)/IRAT + 1
      KJDBC1 = KJDCL2 + (2*NODCL2*MXNODT - 1)/IRAT + 1
      KLAST  = KJDBC1 + (MXNODT - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('EXIDRV','HERPRP',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
C     Allocate ODCPAR
C
      CALL EXIDR1(IODCL1,IODCL2,IODBC1,IODBC2,RODBC1,RODBC2,
     &            IODPP1,IODPP2,RODPP1,RODPP2,
     &            DRED,GRED,IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,
     &            FMAT, DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &            WORK(KACTAB),WORK(KACTCD),WORK(KSCRAB),WORK(KSCRCD),
     &            CCFBT,INDXBT,
     &            WORK(KJDBC2),WORK(KJDCL2),WORK(KJDBC1),
     &            WORK(KLAST),LWRK,MXNODT,IPRINT)
C
      RETURN
      END
C
#if defined (VAR_SINGLE)
#include "single.h"
#endif
C
C  /* Deck exidr1 */
      SUBROUTINE EXIDR1(IODCL1,IODCL2,IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  DRED,GRED,IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,
     &                  FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                  NACTAB,NACTCD,SCRNAB,SCRNCD,
     &                  CCFBT,INDXBT,
     &                  JODBC2,JODCL2,JODBC1,
     &                  WORK,LWORK,MXNODT,IPRINT)
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
#include "ericom.h"
#include "erimem.h"
#include "erista.h"
#include "hertop.h"
#include "symmet.h"
#include "r12int.h"
#include "denfit.h"
#include "iratdef.h" 
C
      INTEGER A, B, C, D
C
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          NACTAB(NODCL1), NACTCD(NODCL2),
     &          SCRNAB(NODBC1), SCRNCD(NODBC2),
     &          FMAT(*), DMAT(*), IFCTYP(NDMT), D2MAT(*), ID2MAT(*),
     &          CCFBT(*), INDXBT(*), 
     &          WORK(LWORK)
      DIMENSION DRED(*), GRED(*),
     &          IDNSCNT(*), IDNSIDX(*), ISHLCNT(*), ISHLIDX(*) 
      DIMENSION IPNTCR(MAXBCH,4), IPNTPP(MAXBCH,2), CSQ(NCSQ1,NCSQ2),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2,2),
     &          BUF(1), IBUF(1), NCOUNT(1)
      DIMENSION ICLSS(NAOBCH,NAOBCH), JODBC2(NAOBCH**2,MXNODT),
     &          JODCL2(NODCL2,2,MXNODT), JODBC1(MXNODT)
#include "ibtfun.h"
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from subroutine ODCDR1',-1)
         WRITE (LUPRI,'(2X,A,2I4)') ' NODCL1, NODCL2 ', NODCL1, NODCL2
         WRITE (LUPRI,'(2X,A,2I4)') ' NODPP1, NODPP2 ', NODPP1, NODPP2
         WRITE (LUPRI,'(2X,A,2L4)') ' ODTRI1, ODTRI2 ', ODTRI1, ODTRI2
         IF (IPRINT .GT. 8) THEN
            CALL HEADER('IODCL1 in ODCDR1',1)
            DO I = 1, NODCL1
               WRITE (LUPRI,'(5X,I4,5X,16I4)') I,(IODCL1(I,J),J=1,NITCL)
            END DO
            CALL HEADER('IODCL2 in ODCDR1',1)
            DO I = 1, NODCL2
               WRITE (LUPRI,'(5X,I4,5X,16I4)') I,(IODCL2(I,J),J=1,NITCL)
            END DO
            CALL HEADER('IODBC1 in ODCDR1',1)
            DO I = 1, NODBC1
               WRITE (LUPRI,'(5X,I4,5X,7I4)') I,(IODBC1(I,J),J=1,NITBC)
            END DO
            CALL HEADER('IODBC2 in ODCDR1',1)
            DO I = 1, NODBC2
               WRITE (LUPRI,'(5X,I4,5X,7I4)') I,(IODBC2(I,J),J=1,NITBC)
            END DO
            CALL HEADER('RODBC1 in ODCDR1',-1)
            CALL OUTPUT(RODBC1,1,NODBC1,1,NRTBC,NODBC1,NRTBC,1,LUPRI)
            CALL HEADER('RODBC2 in ODCDR1',-1)
            CALL OUTPUT(RODBC2,1,NODBC2,1,NRTBC,NODBC2,NRTBC,1,LUPRI)
            CALL HEADER('RODPP1 in ODCDR1',-1)
            CALL OUTPUT(RODPP1,1,NODPP1,1,NRTPP,NODPP1,NRTPP,1,LUPRI)
            CALL HEADER('RODPP2 in ODCDR1',-1)
            CALL OUTPUT(RODPP2,1,NODPP2,1,NRTPP,NODPP2,NRTPP,1,LUPRI)
         END IF
      END IF
C
C
C     COMMON /DOXYZ/
C
      DOX = .TRUE.
      DOY = .TRUE.
      DOZ = .TRUE.
C
      FIRST = .TRUE.
      LAST  = .FALSE.
C
C     Spherical-harmonic transformation matrices
C     ==========================================
C
      CALL SPHMAT(CSQ,NCSQ1,NCSQ2,MAXDER,GDER,WORK,LWORK,IPRINT)
C
C     Arrays for Hermite integrals
C     ============================
C
C
C aspg, 2006-06-27 making ftnchek happy
C
C      CALL HERPRP(INDHER,INDHSQ,IODDHR,IPRINT)
      CALL HERPRP(INDHER,INDHSQ,IODDHR)
C
      CALL LMNPRP(LMNPWR,LMNSYM,IPRINT)
C
C     Active/inactive batches
C     =======================
C
      IELCTR = 1
      CALL ERIACT(NACTAB,SCRNAB,RODBC1,IODBC1,NODBC1,IODCL1,
     &            NODCL1,ODTRI1,IELCTR,IPRINT)
C
      IELCTR = 2
      CALL ERIACT(NACTCD,SCRNCD,RODBC2,IODBC2,NODBC2,IODCL2,
     &            NODCL2,ODTRI2,IELCTR,IPRINT)
C
      IPRKEP = IPRINT
C
      KVEC  = 1
      KTAK  = KVEC + (NAOBCH**2 - 1)/IRAT + 1
      KLAST = KTAK + (NAOBCH**2 - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('EXIDR1','EXISCD',KLAST,LWORK)
C
      DO IODAB = 1, NODCL1
      IF (NACTAB(IODAB) .GT. 0) THEN
         NODSAB = NACTAB(IODAB)
C
         NODTAB = IODCL1(IODAB,1)
         KODSAB = IODCL1(IODAB,2)
         NRUNS = (NODTAB - 1)/MXNODT + 1
         NBTCH = MIN(MXNODT,NODTAB)
         DO IRUN = 1, NRUNS
            ISTR = KODSAB + NBTCH*(IRUN - 1)
            IEND = MIN(ISTR + NBTCH,KODSAB + NODTAB) - 1
            JAB = 1
            DO I = ISTR, IEND
            IF (IODBC1(I,4).NE.0) THEN
               A = IODBC1(I,2)
               B = IODBC1(I,3)
               CALL EXISCD(A,B,DRED,GRED,NAOBCH,
     &                     IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,IODBC2,
     &                     LENGTH,JODBC2(1,JAB),JODCL2(1,1,JAB),
     &                     WORK(KVEC),WORK(KTAK))
               IF (LENGTH.GT.0) THEN
                  JODBC1(JAB) = I
                  JAB = JAB + 1
               END IF
            END IF 
            END DO
            NAB = JAB - 1
C
            IF (NAB.GT.0) THEN
               DO IODCD = 1, IODAB
               IF (NACTCD(IODCD) .GT. 0) THEN
                  NODSCD = NACTCD(IODCD)
                  CALL EXICOM(IODCL1,IODCL2)
                  IJ = 0
C
C aspg, 2006-06-27 making ftnchek happy
C
                  II = 0
C question: should II be incremented in the loops below??
C
                  IPASS  = 0
                  DO JAB = 1, NAB
                  DO JCD = JODCL2(IODCD,1,JAB), JODCL2(IODCD,2,JAB) 
                     I = JODBC1(JAB)
                     J = JODBC2(JCD,JAB)
                     IF (.NOT.(IODAB.EQ.IODCD.AND.J.GT.I)) THEN 
                     IF (IODBC2(J,4).NE.0) THEN
                        IJ = IJ + 1
                        IPNTCR(IJ,1) = IODBC1(I,2)
                        IPNTCR(IJ,2) = IODBC1(I,3)
                        IPNTCR(IJ,3) = IODBC1(J,2)
                        IPNTCR(IJ,4) = IODBC1(J,3)
                        IPNTPP(IJ,1) = IODBC1(I,4)
                        IPNTPP(IJ,2) = IODBC2(J,4)
                        IF (IJ.EQ.MAXBCH) THEN
                           CALL ODCPA2(IJ,II,IODPP1,IODPP2,
     &                                 RODPP1,RODPP2,
     &                                 IPNTCR,IPNTPP,INDHER,INDHSQ,
     &                                 IODDHR,LMNPWR,LMNSYM,IPNTUV,WORK,
     &                                 LWORK,FMAT,DMAT,NDMT,IFCTYP,
     &                                 D2MAT,ID2MAT,BUF,IBUF,CSQ,CCFBT,
     &                                 INDXBT,NCOUNT,
     &                                 IPASS,IPRINT)
                        END IF
                     END IF
                     END IF
                  END DO
                  END DO
                  IF (IJ.GT.0) THEN
                     CALL ODCPA2(IJ,II,IODPP1,IODPP2,RODPP1,RODPP2,
     &                           IPNTCR,IPNTPP,INDHER,INDHSQ,IODDHR,
     &                           LMNPWR,LMNSYM,IPNTUV,WORK,LWORK,
     &                           FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                           BUF,IBUF,CSQ,CCFBT,INDXBT,NCOUNT,IPASS,
     &                           IPRINT)
                  END IF
               END IF
               END DO
            END IF
         END DO 
      END IF
      END DO
C
      IPRINT = IPRKEP
C
      RETURN
      END
C  /* Deck exicom */
      SUBROUTINE EXICOM(IODCL1,IODCL2)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL)
#include "cbieri.h"
#include "symmet.h"
#include "ericom.h"
#include "odclss.h"
#include "ibtfun.h"
C
      IF (.TRUE.) THEN
         ISTBLA = 0
         ISTBLB = 0
         ISTBLC = 0
         ISTBLD = 0
         ISTBLR = 0 
         ISTBLS = 0 
         ISTBLT = 0 
         MLTPA = 1 
         MLTPB = 1 
         MLTPC = 1 
         MLTPD = 1 
         MLTPR = 1 
         MLTPS = 1 
         MLTPT = 1  
         MLTPX = 1 
C
         MLTPZ = MIN(MLTPX,NSPMAX)
         NREDZ = MLTPX/MLTPZ
C
         TKMPAB = .FALSE. 
         TKMPCD = .FALSE. 
         TCMPAB = .FALSE. 
         TCMPCD = .FALSE. 
      END IF
      NODTAB = IODCL1(IODAB, 1)
      KODSAB = IODCL1(IODAB, 2)
      NHKTA  = IODCL1(IODAB, 3)
      NHKTB  = IODCL1(IODAB, 4)
      KHKTA  = IODCL1(IODAB, 5)
      KHKTB  = IODCL1(IODAB, 6)
      KCKTA  = IODCL1(IODAB, 7)
      KCKTB  = IODCL1(IODAB, 8)
      NPRFAB = IODCL1(IODAB,11)
      NCTFAB = IODCL1(IODAB,12)
      NCNTAB = IODCL1(IODAB,13)
      NGTOAB = IODCL1(IODAB,14)
      NPRFA  = IODCL1(IODAB,15)
      NPRFB  = IODCL1(IODAB,16)
      NCTFA  = IODCL1(IODAB,17)
      NCTFB  = IODCL1(IODAB,18)
      ICMATA = IODCL1(IODAB,19)
      ICMATB = IODCL1(IODAB,20)
C
      NODTCD = IODCL2(IODCD, 1)
      KODSCD = IODCL2(IODCD, 2)
      NHKTC  = IODCL2(IODCD, 3)
      NHKTD  = IODCL2(IODCD, 4)
      KHKTC  = IODCL2(IODCD, 5)
      KHKTD  = IODCL2(IODCD, 6)
      KCKTC  = IODCL2(IODCD, 7)
      KCKTD  = IODCL2(IODCD, 8)
      NPRFCD = IODCL2(IODCD,11)
      NCTFCD = IODCL2(IODCD,12)
      NCNTCD = IODCL2(IODCD,13)
      NGTOCD = IODCL2(IODCD,14)
      NPRFC  = IODCL2(IODCD,15)
      NPRFD  = IODCL2(IODCD,16)
      NCTFC  = IODCL2(IODCD,17)
      NCTFD  = IODCL2(IODCD,18)
      ICMATC = IODCL2(IODCD,19)
      ICMATD = IODCL2(IODCD,20)
C
      NODSPQ = NODSAB*NODSCD
      NPRFPQ = NPRFAB*NPRFCD
      NCTFPQ = NCTFAB*NCTFCD
C
      JMAXA  = NHKTA - 1
      JMAXB  = NHKTB - 1
      JMAXC  = NHKTC - 1
      JMAXD  = NHKTD - 1
      JMAXAB = JMAXA + JMAXB + MAXDER
      JMAXCD = JMAXC + JMAXD + MAXDER
      JMAX   = JMAXA + JMAXB + JMAXC + JMAXD + MAXDER
C
      NTUVAB = (JMAXAB + 1)*(JMAXAB + 2)*(JMAXAB + 3)/6
      NTUVCD = (JMAXCD + 1)*(JMAXCD + 2)*(JMAXCD + 3)/6
      NTUV   = (JMAX   + 1)*(JMAX   + 2)*(JMAX   + 3)/6
C
      SPHRA  = KHKTA .NE. KCKTA .OR. MAXDER.GT.0
      SPHRB  = KHKTB .NE. KCKTB .OR. MAXDER.GT.0
      SPHRC  = KHKTC .NE. KCKTC .OR. MAXDER.GT.0
      SPHRD  = KHKTD .NE. KCKTD .OR. MAXDER.GT.0
C
      SPHRAB = SPHRA .OR. SPHRB
      SPHRCD = SPHRC .OR. SPHRD
C
      GCONA  = GENCON .OR. NCTFA .GT. 1 
      GCONB  = GENCON .OR. NCTFB .GT. 1 
      GCONC  = GENCON .OR. NCTFC .GT. 1 
      GCOND  = GENCON .OR. NCTFD .GT. 1 
C
      GCONAB = GCONA .OR. GCONB
      GCONCD = GCONC .OR. GCOND
C
      DIAGAB = NGTOAB .EQ. 1 .AND. ODTRI1
      DIAGCD = NGTOCD .EQ. 1 .AND. ODTRI2
C
      KHKTAB = KHKTA*KHKTB
      KHKTCD = KHKTC*KHKTD
      KCKTAB = KCKTA*KCKTB
      KCKTCD = KCKTC*KCKTD
C
      DOPTH1 = .TRUE.
      DOPTH2 = .FALSE.
C
      RETURN
      END
C  /* Deck exired */
      SUBROUTINE EXIRED(DMAT,NDMT,DRED,GRED,IDNSCNT,IDNSIDX,ISHLCNT,
     &                  ISHLIDX,ICLSS,IODCL1,IODCL2,
     &                  IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  CCFBT,INDXBT,IPRINT,WORK,LWORK)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "ccom.h"
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "erimem.h"
#include "aobtch.h"
#include "odbtch.h"
#include "symmet.h"
#include "infpar.h"
#include "inforb.h"
#include "shells.h"
#include "odclss.h"

#ifdef PRG_DIRAC
! HFXFAC in dcbham.h
#include "dcbham.h"
#endif
      LOGICAL WRTINT_SAVE, FCKINT_SAVE
C
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          DMAT(NBAST,NBAST,NDMT)
C
      DIMENSION ISHLIDX(NAOBCH,NAOBCH), ISHLCNT(NAOBCH),
     &          IDNSIDX(NAOBCH,NAOBCH), IDNSCNT(NAOBCH),
     &          DRED(NAOBCH,NAOBCH), GRED(NAOBCH,NAOBCH), 
     &          ICLSS(NAOBCH,NAOBCH)
C
      DIMENSION WORK(LWORK) 
C
C     Allocations
C
      KCSMAT = 1
      KLAST1 = KCSMAT + NBAST**2
      KGMAX  = 1
      KDMAX  = KGMAX + NAOBCH
      KLAST2 = KDMAX + NAOBCH
      KLAST  = MAX(KLAST1,KLAST2)
      IF (KLAST .GT. LWORK) CALL STOPIT('EXIRED',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      JPRINT = IPRINT
      THRXCH = THRXCH/HFXFAC
C
C     Cauchy-Schwarz factors
C
      DGABAB = .TRUE.
      WRTINT_SAVE = WRTINT
      FCKINT_SAVE = FCKINT
      WRTINT = .FALSE.
      FCKINT = .FALSE.
      CALL DZERO(WORK(KCSMAT),NBAST**2)
      CALL ODCDRV(IODCL1,IODCL2,IODBC1,IODBC2,
     &            RODBC1,RODBC2,IODPP1,IODPP2,
     &            RODPP1,RODPP2,
     &            DUMMY,DUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     &            WORK(KCSMAT),DUMMY,1,IDUMMY,DUMMY,IDUMMY,
     &            CCFBT,INDXBT,WORK(KLAST),LWRK,IPRINT)
      DGABAB = .FALSE.
      WRTINT = WRTINT_SAVE
      FCKINT = FCKINT_SAVE
C
      IF (JPRINT.GT.5) WRITE (LUPRI,'(1X,A,2I5)') 
     &   ' Full and reduced dimensions in EXIRED ', NBAST, NAOBCH
C
C     Reduced Cauchy-Schwarz and density matrices
C
      IF (JPRINT.GT.5)CALL HEADER('Cauchy-Schwarz matrix in EXIRED ',-1)
      CALL CSMEXT(WORK(KCSMAT),NBAST,1,GRED,NAOBCH,NORBBT,JPRINT)
C
      IF (JPRINT.GT.5)CALL HEADER('Density matrix in EXIRED ',-1)
      CALL CSMEXT(DMAT,NBAST,NDMT,DRED,NAOBCH,NORBBT,JPRINT)
C
C     Significant shell pairs
C
      CALL EXISSP(DRED,GRED,NAOBCH,NAOBCH,ISHLIDX,ISHLCNT,WORK(KGMAX),
     &            WORK(KDMAX))
      CALL EXICLS(ISHLCNT,ISHLIDX,ICLSS,IODCL1,IODCL2,IODBC1,IODBC2,
     &            ICLSS,NAOBCH)
C
C     Significant density elements
C
      CALL EXIOSD(NAOBCH,DRED,NAOBCH,WORK(KGMAX),IDNSCNT,IDNSIDX)
C
      IF (JPRINT.GT. 5) THEN
        CALL HEADER('Sorted list of significant shell pairs.',-1)
        CALL HEADER('I, ISHLCNT(I), GMAX(I), ISHLIDX(I,*) in EXIRED',-1)
        DO I = 1, NAOBCH
        IF (ISHLCNT(I).GT.0) THEN
           WRITE (LUPRI,'(2I5,F20.12,(5X,10I5),(/25X,10I5))') 
     &      I,ISHLCNT(I), WORK(KGMAX+I-1), (ISHLIDX(J,I),J=1,ISHLCNT(I))
        END IF
        END DO
        CALL HEADER('Sorted list of significant density elements.',-1)
        CALL HEADER('I, IDNSCNT(I), DMAX(I), IDNSIDX(I,*) in EXIRED',-1)
        DO I = 1, NAOBCH
        IF (IDNSCNT(I).GT.0) THEN
          WRITE (LUPRI,'(2I5,F20.12,(5X,10I5),(/25X,10I5))') 
     &     I,IDNSCNT(I), WORK(KDMAX+I-1), (IDNSIDX(J,I),J=1,IDNSCNT(I))
        END IF
        END DO
      END IF
      RETURN
      END
C  /* Deck exicls */
      SUBROUTINE EXICLS(ISHLCNT,ISHLIDX,ICLSS,IODCL1,IODCL2,
     &                  IODBC1,IODBC2,SIGLOG,NAOBCH)
#include "implicit.h"
#include "priunit.h"
C
#include "odclss.h"
C
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC)
      INTEGER   ISHLIDX(NAOBCH,NAOBCH), ISHLCNT(NAOBCH)
      INTEGER   ICLSS(NAOBCH,NAOBCH)
      LOGICAL   SIGLOG(NAOBCH,NAOBCH) 
C
      DO IMU = 1, NAOBCH
         DO INU = 1, NAOBCH
            SIGLOG(INU,IMU) = .TRUE.
         END DO 
         DO I = 1, ISHLCNT(IMU)
            SIGLOG(ISHLIDX(I,IMU),IMU) = .FALSE.
         END DO
      END DO
      DO I = 1, NODBC1
         IF (SIGLOG(IODBC1(I,2),IODBC1(I,3))) IODBC1(I,4) = 0
      END DO 
      DO I = 1, NODBC2
         IF (SIGLOG(IODBC2(I,2),IODBC2(I,3))) IODBC2(I,4) = 0
      END DO 
C
C     ICLSS
C
      DO IODCD = 1, NODCL2
      DO I = IODCL2(IODCD,2), IODCL2(IODCD,1) + IODCL2(IODCD,2) - 1
         ICLSS(IODBC2(I,2),IODBC2(I,3)) = I 
         ICLSS(IODBC2(I,3),IODBC2(I,2)) = I 
      END DO
      END DO
C
      RETURN
      END
C
C     select_significant_shell_pairs: given density matrix D and
C     Schwartz-Cauchy integrals G, select significant orbital pairs.
C     This can readily be reduced to shell pairs.
C
C  /* Deck exissp */
      SUBROUTINE EXISSP(D,G,LDA,NBATCH,ISHLIDX,ISHLCNT,GMAX,DMAX)
#include "implicit.h"
#include "priunit.h"
#include "erithr.h"
C
C     D - a (LDA, NBATCH) block of reduced density matrix (LDA,LDA)
c     the reduced density matrix contains 1-norms of respective
c     blocks of the density matrix in orbital basis.
C     G - ditto dimension-wise, but of Cauchy-Schwartz matrix.
C     ISHLIDX is an output parameter - sufficiently large! containing
c     GMAX is a upper limit on all values in a row.
C     ISHLCNT - will contain the list length.
C
      DIMENSION D(LDA, NBATCH), G(LDA, NBATCH), 
     &          GMAX(NBATCH), DMAX(NBATCH)
      INTEGER ISHLIDX(NBATCH,NBATCH), ISHLCNT(NBATCH)
      DIMENSION BVAL(NBATCH) 
C
      GMAXX = 0D0
      DO IB=1, NBATCH
         DMAX(IB) = 0D0
         GMAX(IB) = 0D0
         DO IA=1, NBATCH
            IF(D(IA,IB).GT.DMAX(IB)) DMAX(IB) = D(IA,IB)
            IF(G(IA,IB).GT.GMAX(IB)) GMAX(IB) = G(IA,IB)
         END DO
         IF(GMAX(IB).GT.GMAXX) GMAXX = GMAX(IB)
      END DO
C
      NTOT = 0
      DO IB=1, NBATCH
         ISHLCNT(IB) = 0
         DO IA=1, NBATCH
            DMM = MAX(DMAX(IB),DMAX(IA))
            GGD = G(IA,IB)*GMAXX*DMM
            IF(GGD.GT.THRXCH) THEN
               NTOT = NTOT + 1
               ISHLCNT(IB) = ISHLCNT(IB) +1
               ISHLIDX(ISHLCNT(IB),IB) = IA
               BVAL(ISHLCNT(IB)) = GGD
            END IF
         END DO
      CALL EXISRT(ISHLCNT(IB),ISHLIDX(1,IB),BVAL)
      END DO
      RETURN
      END
C
C  /* Deck exiosd */
      SUBROUTINE EXIOSD(NBATCH,DMAT,LDA,GMX,ISHLCNT, ISHLIDX)
#include "implicit.h"
#include "priunit.h"
#include "erithr.h"
C
      DIMENSION DMAT(LDA,NBATCH), GMX(NBATCH)
      DIMENSION ISHLCNT(NBATCH), ISHLIDX(NBATCH,NBATCH), DVAL(NBATCH)
C
      NTOT = 0
      DO IMU=1, NBATCH
         ISHLCNT(IMU) = 0
         DO INU=1, NBATCH
            IF(DMAT(INU,IMU)*GMX(INU)*GMX(IMU) .GT. THRXCH) THEN
               NTOT = NTOT + 1
               ISHLCNT(IMU)  = ISHLCNT(IMU) + 1
               ISHLIDX(ISHLCNT(IMU),IMU) = INU
               DVAL(ISHLCNT(IMU)) = DMAT(INU,IMU)*GMX(INU)
            END IF
         END DO
         CALL EXISRT(ISHLCNT(IMU),ISHLIDX(1,IMU),DVAL)
      END DO
      RETURN
      END
C
C  /* Deck exiscd */
      SUBROUTINE EXISCD(IA,IB,DRED,GRED,NAOBCH,
     &                  IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,
     &                  IODBC2,LENGTH,JODBC2,JODCL2,IVEC,ITAKEN)
#include "implicit.h"
#include "priunit.h"
#include "odclss.h"
C     
C     Determine significant elements (AB|CD) for given AB
C     
      DIMENSION DRED(NAOBCH,NAOBCH), GRED(NAOBCH,NAOBCH)
      DIMENSION IODBC2(NODBC2,NITBC), JODCL2(NODCL2,2)
      INTEGER IDNSIDX(NAOBCH,NAOBCH), IDNSCNT(NAOBCH)
      INTEGER ISHLIDX(NAOBCH,NAOBCH), ISHLCNT(NAOBCH)
      INTEGER ICLSS(NAOBCH,NAOBCH)
      INTEGER JODBC2(NAOBCH**2)
C
      INTEGER IVEC(NAOBCH**2), ITAKEN(NAOBCH**2)
C
      LENGTH = 0
      CALL IZERO(ITAKEN,NAOBCH**2)
      CALL EXIGCD(IA,IB,DRED,GRED,NAOBCH,
     &            IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,
     &            IVEC,LENGTH,ITAKEN)
      CALL EXIGCD(IB,IA,DRED,GRED,NAOBCH,
     &            IDNSCNT,IDNSIDX,ISHLCNT,ISHLIDX,ICLSS,
     &            IVEC,LENGTH,ITAKEN)
      IF (LENGTH.GT.0) THEN
         CALL IZERO(JODCL2(1,2),NODCL2)
         DO I = 1, LENGTH
            JODCL2(IODBC2(IVEC(I),1),2) = JODCL2( IODBC2(IVEC(I),1),2)+1
         END DO
         JODCL2(1,1) = 1 
         DO I = 2, NODCL2
            JODCL2(I,1) = JODCL2(I-1,1) + JODCL2(I-1,2)
         END DO
         DO I = 1, NODCL2
            JODCL2(I,2) = JODCL2(I,1) - 1
         END DO
         DO I = 1, LENGTH
            J = IVEC(I)
            K = IODBC2(J,1)
            L = JODCL2(K,2) + 1
            JODBC2(L) = J 
            JODCL2(K,2) = L 
         END DO
      END IF
C     IF (LENGTH.GT.0) THEN
C        CALL QSORT(JODBC2,LENGTH,4,NCLSRT,IODBC2)
C        DO K = 1, NODCL2
C           JODCL2(K,1) =  0
C           JODCL2(K,2) = -1 
C        END DO
C        ICL0 = 0
C        DO K = 1, LENGTH
C           ICLS = IODBC2(JODBC2(K,1),1)
C           IF (ICLS.NE.ICL0) THEN
C              JODCL2(ICLS,1) = K
C              IF (ICL0.GT.0) JODCL2(ICL0,2) = K - 1
C              ICL0 = ICLS
C           END IF
C        END DO
C        JODCL2(ICLS,2) = LENGTH 
C     END IF
      RETURN
      END
C  /* Deck exigcd */
      SUBROUTINE EXIGCD(A,B,DRED,GRED,NAOBCH,IDNSCNT,IDNSIDX,ISHLCNT,
     &                  ISHLIDX,ICLSS,IVEC,LENGTH,ITAKEN)
#include "implicit.h"
#include "priunit.h"
#include "erithr.h"
C
C     Determine significant elements (AB|CD) for given AB
C     when multiplied with D_AC
C     
      INTEGER A, B, C, D
C
      DIMENSION DRED(NAOBCH,NAOBCH), GRED(NAOBCH,NAOBCH)
      INTEGER IDNSIDX(NAOBCH,NAOBCH), IDNSCNT(NAOBCH)
      INTEGER ISHLIDX(NAOBCH,NAOBCH), ISHLCNT(NAOBCH)
      INTEGER ICLSS(NAOBCH,NAOBCH)
C
      INTEGER IVEC(NAOBCH**2)
      LOGICAL LTAKEN
      INTEGER ITAKEN(NAOBCH**2)
c
      THRAB = THRXCH/GRED(A,B)
      DO IC = 1, IDNSCNT(A)
         C = IDNSIDX(IC,A)
         LTAKEN = .FALSE.
         DO ID = 1, ISHLCNT(C)
            D = ISHLIDX(ID,C)
C        IF (DRED(A,C)*GRED(A,B)*GRED(C,D).LE.THRXCH) GO TO 199
         IF (DRED(A,C)*GRED(C,D).LE.THRAB) GO TO 199
            J = ICLSS(C,D)
            IF (ITAKEN(J).EQ.0) THEN
               LENGTH = LENGTH + 1
               ITAKEN(J) = 1
               IVEC(LENGTH) = J
            END IF
            LTAKEN = .TRUE.
         END DO
 199     CONTINUE
         IF(.NOT.LTAKEN) GO TO 201
      END DO
 201  CONTINUE
      RETURN
      END
C  /* Deck exisrt */
      SUBROUTINE EXISRT(ICNT,IVAL,DVAL)
#include "implicit.h"
c     sorts IVAL of length ICNT according to keys in DVAL
      DIMENSION IVAL(ICNT), DVAL(ICNT)
      DO I = 1, ICNT-1
         IMAX = I
         DO J=I+1, ICNT
            IF(DVAL(J).GT.DVAL(IMAX)) IMAX = J
         END DO
         TMP = DVAL(I)
         DVAL(I) = DVAL(IMAX)
         DVAL(IMAX) = TMP
         ITMP = IVAL(I)
         IVAL(I) = IVAL(IMAX)
         IVAL(IMAX) = ITMP
      END DO
      RETURN
      END
C  /* Deck exindf */
      SUBROUTINE EXINDF(FMAT,DMAT,NDMT,SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &                  INDXBT,IFCTYP,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER (DP5=0.50D0)
      LOGICAL AEQB, CEQD, PEQQ
      INTEGER A, B, C, D
      DIMENSION SO(NCCS,KHKTAB,KHKTCD),
     &          INDEX(NCCS,4), 
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          FMAT(NBASIS,NBASIS,NDMT), DMAT(NBASIS,NBASIS,NDMT),
     &          IFCTYP(NDMT)
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "nuclei.h"

#ifdef PRG_DIRAC
! HFXFAC in dcbham.h
#include "dcbham.h"
#endif

#include "aobtch.h"
#include "hertop.h"
#include "symmet.h"
C
#include "ibtfun.h"
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine EXINDF',-1)
#if ! defined (PRG_DIRAC)
      DO K = 1, NDMT
      IF (IFCTYP(K).NE.2.AND.IFCTYP(K).NE.12) THEN
         WRITE (LUPRI,'(2X,A,I5)') ' IFCTYP(K):', IFCTYP(K)
         CALL QUIT('Wrong value of IFCTYP in EXINDF')
      END IF
      END DO
#endif
C
      NCCS1 = (NPQBCS - NPPBCS)*NCTFAB*NCTFCD
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      CALL ERIPNT(INDEX,0,0,0,0,IPNTCR,INDXBT,1)
C
      FAC = -HFXFAC
      IF (DIAGAB) FAC = DP5*FAC
      IF (DIAGCD) FAC = DP5*FAC
      IF (IFCTYP(1).EQ.2) FAC = DP5*FAC
      ICMPAB = 0
      DO ICMPA = 0, KHKTA - 1
      DO ICMPB = 0, KHKTB - 1
         ICMPAB = ICMPAB + 1
         IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
         ICMPCD = 0
         DO ICMPC = 0, KHKTC - 1
         DO ICMPD = 0, KHKTD - 1
            ICMPCD = ICMPCD + 1
            IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
            IF (IODDAB .EQ. IODDCD) THEN
            IF ((IFCTYP(1).EQ.12).OR.(IFCTYP(1).EQ.22)) THEN
               DO I = 1, NCCS
                  EINT = FAC*SO(I,ICMPAB,ICMPCD)
                  IF (ABS(EINT) .GT. THRXCH) THEN
                     A = INDEX(I,1) + ICMPA
                     B = INDEX(I,2) + ICMPB
                     C = INDEX(I,3) + ICMPC
                     D = INDEX(I,4) + ICMPD
             DO K = 1, NDMT
                     FMAT(C,A,K) = FMAT(C,A,K) + EINT*DMAT(D,B,K)
                     FMAT(D,A,K) = FMAT(D,A,K) + EINT*DMAT(C,B,K)
                     FMAT(C,B,K) = FMAT(C,B,K) + EINT*DMAT(D,A,K)
                     FMAT(D,B,K) = FMAT(D,B,K) + EINT*DMAT(C,A,K)
             END DO
                  END IF
               END DO
            ELSE
               DO I = 1, NCCS
                  EINT = FAC*SO(I,ICMPAB,ICMPCD)
                  IF (ABS(EINT) .GT. THRXCH) THEN
                     A = INDEX(I,1) + ICMPA
                     B = INDEX(I,2) + ICMPB
                     C = INDEX(I,3) + ICMPC
                     D = INDEX(I,4) + ICMPD
             DO K = 1, NDMT
                     FMAT(C,A,K) = FMAT(C,A,K) + EINT*DMAT(D,B,K)
                     FMAT(D,A,K) = FMAT(D,A,K) + EINT*DMAT(C,B,K)
                     FMAT(C,B,K) = FMAT(C,B,K) + EINT*DMAT(D,A,K)
                     FMAT(D,B,K) = FMAT(D,B,K) + EINT*DMAT(C,A,K)
                     FMAT(A,C,K) = FMAT(A,C,K) + EINT*DMAT(B,D,K)
                     FMAT(A,D,K) = FMAT(A,D,K) + EINT*DMAT(B,C,K)
                     FMAT(B,C,K) = FMAT(B,C,K) + EINT*DMAT(A,D,K)
                     FMAT(B,D,K) = FMAT(B,D,K) + EINT*DMAT(A,C,K)
             END DO
                  END IF
               END DO
            END IF
            END IF
         END DO
         END DO
      END DO
      END DO
      RETURN
      END
C  /* Deck nclsrt */
      INTEGER FUNCTION NCLSRT(IA,IB,ICLPNT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION ICLPNT(*)
      IF (ICLPNT(IA).LE.ICLPNT(IB)) NCLSRT = -1
      IF (ICLPNT(IA).EQ.ICLPNT(IB)) NCLSRT =  0
      IF (ICLPNT(IA).GT.ICLPNT(IB)) NCLSRT =  1
      RETURN
      END
C  /* Deck exiexp */
      SUBROUTINE EXIEXP(EXPAB,EXPCD,ALPHA,FACPQ,HEXPP,HEXPQ,
     &                  RODPE1,RODPE2,IPNTPP,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "r12int.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0)
      INTEGER AB, CD, PQ, RST
      DIMENSION EXPAB(NPQBCS,NPRFAB,3),
     &          EXPCD(NPQBCS,NPRFCD,3),
     &          ALPHA(NPQBCS,NPRFPQ),
     &          HEXPP(NPQBCS,NPRFPQ),
     &          HEXPQ(NPQBCS,NPRFPQ),
     &          FACPQ(NPQBCS,NPRFPQ),
     &          RODPE1(NODPP1,  2), RODPE2(NODPP2,  2),
     &          IPNTPP(MAXBCH,2)
#include "cbieri.h"
#include "odclss.h"
#include "ericom.h"
#include "hertop.h"
#include "aobtch.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from EXIEXP','*',103)
      END IF
C
      DO I = 1, NPRFAB
      DO PQ = 1, NPQBCS
         AB = IPNTPP(PQ,1) + I - 1
         EXPA = RODPE1(AB,1)
         EXPB = RODPE1(AB,2)
         PINV = D1/(EXPA + EXPB)
         EXPAP = EXPA*PINV
         EXPBP = EXPB*PINV
         EXPAB(PQ,I,1) = EXPAP
         EXPAB(PQ,I,2) = EXPBP
         EXPAB(PQ,I,3) = PINV
      END DO
      END DO
C
      DO I = 1, NPRFCD
      DO PQ = 1, NPQBCS
         CD = IPNTPP(PQ,2) + I - 1
         EXPC = RODPE2(CD,1)
         EXPD = RODPE2(CD,2)
         QINV = D1/(EXPC + EXPD)
         EXPCQ = EXPC*QINV
         EXPDQ = EXPD*QINV
         EXPCD(PQ,I,1) = EXPCQ
         EXPCD(PQ,I,2) = EXPDQ
         EXPCD(PQ,I,3) = QINV
      END DO
      END DO
C
      IJ = 0
      DO J = 1, NPRFCD
      DO I = 1, NPRFAB
         IJ = IJ + 1
         DO PQ = 1, NPQBCS
            PINV = EXPAB(PQ,I,3) 
            QINV = EXPCD(PQ,J,3) 
            ALFA = D1/(PINV+QINV)
            ALPHA(PQ,IJ) = - ALFA - ALFA
            IF (IPNTPP(PQ,1).EQ.IPNTPP(PQ,2)) THEN
               FACPQ(PQ,IJ) = DP5*SQRT(PINV*QINV*ALFA)
            ELSE 
               FACPQ(PQ,IJ) = SQRT(PINV*QINV*ALFA)
            END IF
         END DO
      END DO
      END DO
      RETURN
      END
