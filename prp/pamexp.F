C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dirac, a relativistic ab initio electronic structure program,
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : pamexp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpexp */
      SUBROUTINE PRPEXP(EXPVAL,WORK,LWORK)
C*****************************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*),WORK(LWORK)
C
#include "dcbexp.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('PRPEXP')
#include "memint.h"
C
c     CALL MEMCHK('begin prpexp',WORK,1)
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KEXPST ,2*NEXPP   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT  ,NZ*N2BBASX*(NOPEN+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT  ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTRI  ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN ,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST ,NZ        ,WORK,KFREE,LFREE)
C
      CALL PRPEX1(EXPVAL,WORK(KEXPST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),WORK(KPTRI),
     &            WORK(KONEIN),WORK(KFIRST),WORK(KFREE),
     &            LFREE)
C
C
C     Memory deallocation
C
      CALL MEMREL('PRPEXP',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PRPEXP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex1 */
      SUBROUTINE PRPEX1(EXPVAL,LEXPST,CMO,EIG,DMAT,PMAT,PTRI,OP1INT,
     &                  FIRST,WORK,LWORK)
C*****************************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0)
C
      LOGICAL FIRST(NZ),LEXPST(2,*),FILEX
      DIMENSION EXPVAL(*),CMO(*),DMAT(*),PMAT(*),PTRI(*),
     &          OP1INT(*),EIG(*),WORK(LWORK)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   GETIOPEN
C
#include "dgroup.h"
ctec dcborb included
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
C
C
C     Initialize
C
#include "memint.h"
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Initialize
C     ==========
C
c     CALL MEMCHK('begin prpex1',WORK,1)
      CALL DZERO(EXPVAL,NEXPP)
      N2EXP = 2*NEXPP
      CALL LSET(N2EXP,.FALSE.,LEXPST)
C
C     Generate density matrix
C     =======================
C
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,1,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Expectation values normal procedure
C    ************************************
C
      IF(.NOT.ORBANA) THEN
        CALL GENDEN(DMAT,CMO,1,IPREXP)
        DO IOPEN = 1,NOPEN
           CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPREXP)
        END DO
        CALL TITLER('Expectation values','*',127)
        WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
        CALL PRSYMB(LUPRI,'-',76,0)
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          CALL PRPEX2(LU1INT,INDXPR,EXPVAL(I),LEXPST(1,I),
     &        DMAT,PMAT,PTRI,OP1INT,FIRST,ORBANA,WORK,LWORK,IPREXP)
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
     &       PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.',
     &       LEXPST(1,I),LEXPST(2,I)
          CALL FLSHFO(LUPRI)
        ENDDO
        CALL PRSYMB(LUPRI,'-',76,0)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
C
C     Expectation values by orbital
C    ************************************
C
      ELSE
        CALL TITLER('Expectation values','*',127)
        WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
        CALL PRSYMB(LUPRI,'-',76,0)
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          CALL PRSYMB(LUPRI,'-',60,4)
          WRITE(LUPRI,'(4X,A16,A44)')
     &          PRPNAM(INDXPR),
     &          ' : Expectation value for individual orbitals'
          EXPVAL(I) = D0
          DO 10 IFRP = 1,NFSYM
            IF(NORB(IFRP).EQ.0) GOTO 10
            DO J = 1,NOCC(IFRP)
              CALL DENORB(DMAT,J,1,IFRP,CMO,IPREXP)
              CALL PRPEX2(LU1INT,INDXPR,EXPBUF,LEXPST(1,I),
     &          DMAT,PMAT,PTRI,OP1INT,FIRST,ORBANA,WORK,LWORK,IPREXP)
              IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 20
              IOPEN = GETIOPEN(J,IFRP)
              EXPBUF = EXPBUF*DF(IOPEN)
              EXPVAL(I) = EXPVAL(I) + EXPBUF
              FMT = MXFORM(EXPBUF,15)
              WRITE(LUPRI,'(4X,A3,I5,3X,'//FMT//',1X,A4)')
     &          FREP(IFRP),J,EXPBUF,'a.u.'
              CALL FLSHFO(LUPRI)
            ENDDO
 10       CONTINUE
 20       CONTINUE
          CALL PRSYMB(LUPRI,'-',60,4)
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A8,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
     &     'Total : ',EXPVAL(I),'a.u.',
     &     LEXPST(1,I),LEXPST(2,I)
          CALL PRSYMB(LUPRI,'-',60,4)
        ENDDO
C        CALL PRSYMB(LUPRI,'-',76,0)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
      ENDIF
C
      CLOSE(LU1INT,STATUS = 'KEEP')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPEXP: Coefficient file not found !'
      CALL QUIT('PRPEXP: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex2 */
      SUBROUTINE PRPEX2(LU1INT,INDXPR,EXPVAL,LEXPST,
     &           DMAT,PMAT,PTRI,OP1INT,FIRST,ORBFLG,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
C
      LOGICAL FIRST(NZ),LEXPST(2),ORBFLG
      DIMENSION DMAT(N2BBASX,NZ,*),PMAT(N2BBASX,NZ),
     &          PTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
C     Check symmetry
C
      EXPVAL = D0
c       CALL MEMCHK('begin prpex2',WORK,1)
C
C     1. Molecular point group symmetry
C
      IREP = IPRPSYM(INDXPR)-1
      IF(IREP.NE.0) LEXPST(1) = .TRUE.
C
C     2. Time reversal symmetry
C
      ITIM = IPRPTIM(INDXPR)
      IF(ITIM.NE.1) LEXPST(2) = .TRUE.
C
C     If not totally symmetric, return with zero value
C
      IF (LEXPST(1).OR.LEXPST(2)) RETURN
C
C     Totally symmetric operator ==> proceed
C
      CALL LSET(NZ,.TRUE.,FIRST)
      CALL PRPMAO(LU1INT,INDXPR,PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C
C     Get expectation value, using only the matrices containing
C     integrals
C
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        IQ = IPQTOQ(IZ,IREP)
        IH = IHQMAT(IQ,ITIM)
        IF(IH.EQ.1) THEN
          CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
        ELSEIF(IH.EQ.2) THEN
          CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
        ENDIF
c       CALL MEMCHK('voor butobs',WORK,1)
        CALL BUTOBS(PMAT(1,IZ),1,WORK,LWORK)
        IF (.NOT. ORBFLG) THEN
           DO IOPEN = 0,NOPEN
              EXPVAL = EXPVAL + DF(IOPEN)
     &           *DDOT(N2BBASX,DMAT(1,IZ,IOPEN+1),1,PMAT(1,IZ),1)
            END DO
         ELSE
            EXPVAL = EXPVAL + DDOT(N2BBASX,DMAT(1,IZ,1),1,PMAT(1,IZ),1)
         END IF
      ENDIF
      ENDDO
      EXPVAL = EXPVAL + EXPVAL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_efg */
      SUBROUTINE NUC_EFG(IPRINT,DOATOM,AXIS,TMP)
C*****************************************************************************
C
C     Calculates nuclear contribution to electric field gradients
C
C     Written by T.Enevoldsen  - Sep 1996
C     Last revision by W.A. de Jong - April 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "mxorb.h"
      PARAMETER (D3 = 3.0D0)
      DIMENSION  TMP(3), ITMP(3), PVAL(9),
     &     AXIS(3,3,MXCENT),GEOM(3*8*MXCENT)
      LOGICAL    DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
      CALL DZERO(EFGNUC,9*MXCENT)
      KK = 1
      NATOM = 1
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               LL = 1
               DO 100 I2 = 1, NUCIND
                  IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
                  DO 110 ISYM = 0, MAXOPR
                     IF (IBTAND(ISYM,ISTBNU(I2)) .EQ. 0) THEN
                        IF (LL .NE. KK) THEN
                           XCOOR = GEOM(LL) - GEOM(KK)
                           YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                           ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                           R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                           R5 = R2*R2*SQRT(R2)
                           EFGNUC(1,1,NATOM) = EFGNUC(1,1,NATOM)
     &                             + CHARGE(I2)*(D3*XCOOR*XCOOR - R2)/R5
                           EFGNUC(2,2,NATOM) = EFGNUC(2,2,NATOM)
     &                             + CHARGE(I2)*(D3*YCOOR*YCOOR - R2)/R5
                           EFGNUC(3,3,NATOM) = EFGNUC(3,3,NATOM)
     &                             + CHARGE(I2)*(D3*ZCOOR*ZCOOR - R2)/R5
                           EFGNUC(1,2,NATOM) = EFGNUC(1,2,NATOM)
     &                             + CHARGE(I2)*D3*XCOOR*YCOOR/R5
                           EFGNUC(1,3,NATOM) = EFGNUC(1,3,NATOM)
     &                             + CHARGE(I2)*D3*XCOOR*ZCOOR/R5
                           EFGNUC(2,3,NATOM) = EFGNUC(2,3,NATOM)
     &                             + CHARGE(I2)*D3*YCOOR*ZCOOR/R5
                        END IF
                        LL = LL + 3
                     END IF
 110              CONTINUE
 100           CONTINUE
               EFGNUC(2,1,NATOM) = EFGNUC(1,2,NATOM)
               EFGNUC(3,1,NATOM) = EFGNUC(1,3,NATOM)
               EFGNUC(3,2,NATOM) = EFGNUC(2,3,NATOM)
               IF (IPRINT .GE. 2) THEN
                  WRITE (LUPRI ,'(/2X,2A,/)')
     &               'Nuclear EFG Tensor of nucleus: ',NAMDEP(NATOM)
                  CALL OUTPUT(EFGNUC(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
               END IF
               NATOM = NATOM + 1
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE
C
C     This part of the code is incorrect: one cannot split the
C     electronic and nuclear contributions and diagonalize them
C
C     diagonalize efg tensor
C     DO IATOM = 1, NUCIND
C        IF (DOATOM(I)) THEN
C           CALL DUNIT(AXIS(1,1,IATOM),3)
C           IJ = 1
C           DO I = 1, 3
C              DO J = 1, I
C                 PVAL(IJ) = EFGNUC(I,J,IATOM)
C                 IJ = IJ + 1
C              ENDDO
C           ENDDO
C           CALL JACO(PVAL,AXIS(1,1,IATOM),3,3,3,TMP,ITMP)
C           PVAL(2) = PVAL(3)
C           PVAL(3) = PVAL(6)
C           IZ = IDAMAX(3,PVAL(1),1)
C           CALL DSWAP(1,PVAL(3),1,PVAL(IZ),1)
C           IX = IDAMIN(3,PVAL(1),1)
C           CALL DSWAP(1,PVAL(1),1,PVAL(IX),1)
C           EFGNUC(1,1,IATOM) = PVAL(1)
C           EFGNUC(2,2,IATOM) = PVAL(2)
C           EFGNUC(3,3,IATOM) = PVAL(3)
C         END IF
C      ENDDO
C     print section
C     IF (IPRINT .GE. 4) THEN
C        WRITE(LUPRI,'(A/)')
C    &        '  Principal axis to nuclear contribution'
C        DO I = 1, NUCIND
C           IF (DOATOM(I)) THEN
C              WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(I)
C              CALL OUTPUT(AXIS(1,1,I),1,3,1,3,3,3,1,LUPRI)
C           END IF
C        ENDDO
C     END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck getiopen */
      INTEGER FUNCTION GETIOPEN(IORBITAL,IFRP)
C*****************************************************************************
C
C     Which shell does electron orbital IORBITAL belong to?
C
C     Written by J. Thyssen - Feb 11 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
      IF (IORBITAL .LT. 1)
     &   CALL QUIT('illegal value in GETIOPEN(neg)')
C
      IF (IORBITAL .LE. NISH(IFRP)) THEN
         GETIOPEN = 0
         RETURN
      END IF
C
      ILOW = NISH(IFRP)+1
      DO IOPEN = 1,NOPEN
         IHIGH = ILOW + NACSH(IFRP,IOPEN)
         IF ( (IORBITAL .GE. ILOW) .AND. (IORBITAL .LE. IHIGH) ) THEN
            GETIOPEN = IOPEN
            RETURN
         END IF
         ILOW = IHIGH+1
      END DO
      CALL QUIT('illegal value in GETIOPEN(pos)')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpexp_mp2 */
      SUBROUTINE PRPEXP_MP2(DMO,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate correlation contribution to static first order properties
C     defined in /CBIEXP/
C
C     Based on PRPEXP
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Joost van Stralen, april 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION DMO(*)
      DIMENSION WORK(LFREE)
C
#include "dcbexp.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('PRPEXP_MP2')
C
C     We come from RELCCSD, The total memory available for the WORK array,
C     LWORK is equal to LFREE.
C
      KFRSAV=KFREE
C
C     Memory allocation
C     -----------------
      CALL MEMGET('REAL',KCMO ,N2BBASXQ   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT      ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KEXPST,2*NEXPP   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,NNBBASX   ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
      CALL MEMGET('REAL',KEXP  ,NEXPP     ,WORK,KFREE,LFREE)
C
      CALL MEMCHK('voor prpex1_mp2',WORK,1)
      CALL PRPEX1_MP2(WORK(KEXP),WORK(KEXPST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),WORK(KPTRI),
     &            WORK(KONEIN),DMO,WORK(KFIRST),WORK,
     &            KFREE,LFREE)
C
C
C     Memory deallocation
C     -------------------
      CALL MEMREL('PRPEXP_MP2',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPEXP_MP2')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex1_mp2 */
      SUBROUTINE PRPEX1_MP2(EXPVAL,LEXPST,CMO,EIG,DMAT,PMAT,PTRI,OP1INT,
     &                      DMO,FIRST,WORK,KFREE,LWORK)
C*****************************************************************************
C
C     Calculate correlation contribution to static first order properties
C     defined in /CBIEXP/
C
C     Based on PRPEXP1
C
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Joost van Stralen, april 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0)
C
      LOGICAL FIRST(NZ),LEXPST(2,*),FILEX
      DIMENSION EXPVAL(*),CMO(*),DMAT(*),PMAT(*),PTRI(*),
     &          OP1INT(*),EIG(*),DMO(*),WORK(*)
      CHARACTER MXFORM*6,FMT*6
C
#include "dgroup.h"
ctec dcborb included
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
C
      CALL QENTER('PRPEX1_MP2')
C
C     Initialize
C     ==========
C
#include "memint.h"
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
      CALL MEMCHK('voor dzero',WORK,1)
      CALL DZERO(EXPVAL,NEXPP)
      CALL MEMCHK('na dzero',WORK,1)
      N2EXP = 2*NEXPP
      CALL LSET(N2EXP,.FALSE.,LEXPST)
C
C     Open coefficient file
C     =====================
C
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP_MP2')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,1,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C
C     Transform DMO to 'orthonormal SO' basis
C     =======================================
C
      CALL MEMCHK('voor genden_mp2',WORK,1)
      CALL GENDEN_MP2(DMO,DMAT,CMO,WORK,KFREE,LFREE,IPREXP)
      CALL MEMCHK('na genden_mp2',WORK,1)
C  
C     Calculate correlation contribution to property
C     ==============================================
C
      CALL TITLER('MP2 correlation contribution to property','*',115)
      WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
      CALL PRSYMB(LUPRI,'-',76,0)
      DO I = 1,NEXPP
        INDXPR = LEXPP(I)
        CALL PRPEX2(LU1INT,INDXPR,EXPVAL(I),LEXPST(1,I),
     &      DMAT,PMAT,PTRI,OP1INT,FIRST,ORBANA,WORK,LWORK,IPREXP)
        FMT = MXFORM(EXPVAL(I),15)
        WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')   
     &     PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.',
     &     LEXPST(1,I),LEXPST(2,I)
        CALL FLSHFO(LUPRI)
      ENDDO
      CALL MEMCHK('na verwachtingswaarde',WORK,1)
C
C     Write expectation values to file
C     ================================
C
      CALL WEXPVTOF(EXPVAL)
c     CALL READEXPV(EXPVAL)
C
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(4X,A)')
     & 's0 = T : Expectation value zero by point group symmetry.'
      WRITE(LUPRI,'(4X,A)')
     & 't0 = T : Expectation value zero by time reversal symmetry.'
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
      CALL QEXIT('PRPEX1_MP2')
C
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPEXP_MP2: Coefficient file not found !'
      CALL QUIT('PRPEXP_MP2: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wexpvtof */
      SUBROUTINE WEXPVTOF(EXPVAL)
C*****************************************************************************
C
C     Write expectation values to file
C
C     Written by Joost van Stralen, June 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*)
#include "dcbexp.h"
#include "dcbprp.h"
C
      CALL QENTER('WEXPVTOF')
C
      OPEN (LUCORP,FILE='CORPROP',FORM='UNFORMATTED',STATUS='NEW')
      WRITE (LUCORP) (EXPVAL(I),I=1,NEXPP)
      CLOSE (LUCORP,STATUS='KEEP')
C
      CALL QEXIT('WEXPVTOF')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readexpv */
      SUBROUTINE READEXPV(EXPVAL)
C*****************************************************************************
C
C     read expectation values from file
C
C     Written by Joost van Stralen, June 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*)
#include "dcbexp.h"
#include "dcbprp.h"
C
      CALL QENTER('READEXPV')
C
      CALL DZERO(EXPVAL,NEXPP)
C
      OPEN (LUCORP,FILE='CORPROP',FORM='UNFORMATTED',STATUS='OLD')
      REWIND (LUCORP)
      READ (LUCORP) (EXPVAL(I),I=1,NEXPP)
      CLOSE (LUCORP)
C
      CALL QEXIT('READEXPV')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck  */
      SUBROUTINE EXPTOT(EXPVAL,EXP_CORR,EXP_TOT)
C*****************************************************************************
C
C     put HF and correlation contribution to expectation value together
C
C     Written by Joost van Stralen, June 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*), EXP_CORR(*), EXP_TOT(*)
#include "dcbexp.h"
#include "dcbprp.h"
C
      CALL QENTER('EXPTOT')
C
      CALL DZERO(EXP_TOT,NEXPP)
C
      DO I = 1,NEXPP
         EXP_TOT(I) = EXPVAL(I) + EXP_CORR(I)
      ENDDO
C
      CALL QEXIT('EXPTOT')
      RETURN
      END
