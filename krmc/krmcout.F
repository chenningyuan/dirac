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
C FILE    : dirout.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE KRMCOUT(WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Driver to KRMCOU1
C
C     Written by J. Thyssen - Nov 30 1998
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WORK(*)
C
#include "dcbkrmc_itinfo.h"
#include "dcbopt.h"
C
      CALL QENTER('KRMCOUT')
      KFRSAV = KFREE
C
      CALL MEMGET('REAL',KDINFO,LDINFO*MOPT_MXMACRO,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIINFO,LIINFO*MOPT_MXMACRO,WORK,KFREE,LFREE)
C
      CALL KRMCOU1(WORK(KDINFO),LDINFO,WORK(KIINFO),LIINFO,
     &             WORK,KFREE,LFREE)
C
      CALL QEXIT('KRMCOUT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE KRMCOU1(DINFO,LDINFO,IINFO,LIINFO,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Output from wave function optimization
C
C     Written by J. Thyssen - Nov 30 1998
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "frame.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), DINFO(LDINFO,*), IINFO(LIINFO,*)
C
C     local variables
      LOGICAL MOISNO
      CHARACTER*72 VECPRI(2)
C
      CALL QENTER('KRMCOU1')
      KFRSAV = KFREE
C
C
C
      CALL HEADER('Optimization statistics for Dirac',-1)
C
      CALL TSTAMP(' ',LUPRI)
C
      REWIND LUITFO
      ITM = 0
  100 CONTINUE
         ITM = ITM + 1
         IF (ITM .GE. MOPT_MXMACRO) GOTO 110
         READ(LUITFO,ERR=110,END=120)
     &      (DINFO(I,ITM),I=1,LDINFO), (IINFO(I,ITM),I=1,LIINFO)
      GOTO 100
  110 CONTINUE
      WRITE(LUPRI,9000)
  120 ITM = ITM - 1
      IF (ITM .EQ. 0) THEN
         WRITE(LUPRI,9001)
         GOTO 200
      END IF
 9000 FORMAT(/'KRMCOU1: I/O error reading file KRMC ITINFO')
 9001 FORMAT(/'KRMCOU1: No statistics found on KRMC ITINFO')
C
      EMCSCF = DINFO(3,ITM)
      GRDNRM = DINFO(15,ITM)
C
      WRITE(LUPRI,9100) (IINFO(1,I),IINFO(2,I),DINFO(3,I),DINFO(15,I),
     &                    DINFO(6,I),DINFO(8,I),I=1,ITM)
      WRITE(LUPRI,9110) (IINFO(1,I),
     &                    DINFO(11,I),DINFO(12,I),
     &                    DINFO(13,I),DINFO(14,I),
     &                    DINFO(15,I),I=1,ITM)
C
      WRITE(LUPRI,9120) 
      DO I = 1,ITM
         WRITE(LUPRI,9130) IINFO(1,I),IINFO(2,I),
     &        IINFO(3,I),IINFO(4,I),IINFO(5,I),IINFO(6,I),
     &        DINFO(16,I),DINFO(18,I),DINFO(17,I)
      END DO
C
      WRITE(LUPRI,9140) (IINFO(1,I),DINFO(1,I),DINFO(2,I),
     &                   DINFO(3,I),I=1,ITM)
      WRITE(LUPRI,9150) (IINFO(1,I),DINFO(4,I),DINFO(5,I),
     &                   DINFO(6,I),I=1,ITM)
      WRITE(LUPRI,9160) (IINFO(1,I),DINFO(7,I),DINFO(10,I),
     &                   DINFO(9,I),I=1,ITM)
      WRITE(LUPRI,9170) (IINFO(1,I),DINFO(8,I),DINFO(19,I),DINFO(20,I),
     &                   DINFO(21,I),I=1,ITM)
C
 9100 FORMAT(//'  ITER ITMIC',T18,'EMCSCF',T35,'GRDNRM',
     &       T49,'RATIO',T60,'STPLNG'/1X,69('-'),
     &       /,(2I5,F20.12,F15.10,F10.6,F15.10))
 9110 FORMAT(1P//'  ITER',T11,'GCINRM',T26,'GOEENRM',T41,'GOEPNRM',
     &       T56,'GOBNRM',T71,'GRDNRM'
     &       /1X,78('-')/(I5,5(D15.5)))
 9120 FORMAT(//'  ITER ITMIC NZRED',T41,'TIMMAC',
     &       T51,'TIMITR',T61,'TIMMIC',/,69('-'))
 9130 FORMAT(I5,I6,I6,' (',I3,' C;',I3,' E;',I3,' P)',3F10.2)
 9140 FORMAT(//' ITER',T15,'EMY',T35,'EACTIV',T55,'EMCSCF',
     &     /69('-'),
     &     /(I5,3F20.12))
 9150 FORMAT(//' ITER',T15,'DEPRED',T35,'DEACT',T55,'RATIO',
     &     /69('-'),
     &     /(I5,3F20.12))
 9160 FORMAT(//' ITER',T10,'BETA',T25,'GAMMA',T43,'RTRUST',
     &     /69('-'),
     &     /(I5,F16.8,F12.8,F20.12))
 9170 FORMAT(//' ITER',T10,'STPLEN',T26,'STPLCI',T42,'STPLEE',
     &     T58,'STPLEP',
     &     /69('-'),
     &     /(I5,4F16.12))

C
C
C
  200 CONTINUE
C
C
C     ******************************
C     *** Print detailed timings ***
C     ******************************
C
C
      CALL HEADER('Detailed timings for 2nd order optimization',-1)
C
      WRITE(LUPRI,9600)
      WRITE(LUPRI,9610) 'ROPTST','Start guess',CPUOPT(1),WLLOPT(1)
      WRITE(LUPRI,9610) '  RCIST','CI start guess',CPUOPT(2),WLLOPT(2)
C
      WRITE(LUPRI,9610) 'RTRACTL','MO 4-index transformation',
     &     CPUOPT(7),WLLOPT(7)
C
      WRITE(LUPRI,9610) 'RGRAD','Gradient',CPUOPT(9),WLLOPT(9)
      WRITE(LUPRI,9610) '  RMAKDM','Density matrices',
     &     CPUOPT(8),WLLOPT(8)
      WRITE(LUPRI,9610) '  RFCKMAT','FV and FC',CPUOPT(12),WLLOPT(12)
      WRITE(LUPRI,9610) '  RGETH2','FQ and H2AC',CPUOPT(3),WLLOPT(3)
      WRITE(LUPRI,9610) '  RCIGRAD','CI Gradient',CPUOPT(10),WLLOPT(10)
C
      WRITE(LUPRI,9610) 'XRSSVC','Sigma vectors',
     &     CPUOPT(18),WLLOPT(18)
      WRITE(LUPRI,9610) '  XRSSEP','OC, CO, and OO svc',
     &     CPUOPT(13),WLLOPT(13)
      WRITE(LUPRI,9610) '    XRSTDM','DVT and PVT',
     &     CPUOPT(14),WLLOPT(14)
      WRITE(LUPRI,9610) '    GMOLITX','FCX, FVX, and FVT',
     &     CPUOPT(5),WLLOPT(5)
      WRITE(LUPRI,9610) '    FMOLI','One index transform',
     &     CPUOPT(6),WLLOPT(6)
      WRITE(LUPRI,9610) '    RGETH2TX','FQX, FQT, and H2ACX',
     &     CPUOPT(4),WLLOPT(4)
      WRITE(LUPRI,9610) '      RTR1H1','one index transforms',
     &     CPUOPT(15),WLLOPT(15)
      WRITE(LUPRI,9610) '    RSIGOC, RSIGOO','OC and OO sigma vectors',
     &     CPUOPT(16),WLLOPT(16)
      WRITE(LUPRI,9610) '    RSIGCO','CO sigma vectors',
     &     CPUOPT(17),WLLOPT(17)
      WRITE(LUPRI,9610) '  XRSSCI','CC sigma vectors',
     &     CPUOPT(11),WLLOPT(11)
C
 9600 FORMAT(' Subroutine',T23,'Description',
     &     T47,'CPU time',T59,'Wall time',/)
 9610 FORMAT(T2,A,T23,A,T44,F10.2,T56,F10.2)
C
C     ************************
C     *** Print energy etc ***
C     ************************
C
C
      CALL HEADER('Final results from Dirac',-1)
C
      CALL TSTAMP(' ',LUPRI)
C
      WRITE(LUPRI,'(/)')
      IF (MCTYPE .GE. JMCMIN .OR. IOPT_ISTATE .GT. 1) THEN
C
C        Write out KR-KMCSCF energy
C
         WRITE(LUPRI,9200) IOPT_ISTATE,IOPT_SYMMETRY,
     &        EMCSCF,POTNUC,EMCSCF-POTNUC
         WRITE(LUPRI,9220) GRDNRM
      ELSE
C
C        Write out KR-KMCSCF energy
C
         WRITE(LUPRI,9210) EMCSCF,POTNUC,EMCSCF-POTNUC
         WRITE(LUPRI,9220) GRDNRM
      END IF
 9200 FORMAT(  T6,'State number         :',T31,I5,
     &       /,T6,'Symmetry             :',T31,I5,
     &      //,T6,'Final KR-MCSCF energy:',T31,F20.12,
     &       /,T6,'Nuclear repulsion    :',T31,F20.12,
     &       /,T6,'Electronic energy    :',T31,F20.12)
 9210 FORMAT(/,T6,'Final DHF energy     :',T31,F20.12,
     &       /,T6,'Nuclear repulsion    :',T31,F20.12,
     &       /,T6,'Electronic energy    :',T31,F20.12)
 9220 FORMAT(/,T6,'Final gradient norm  :',T31,F20.12)
      WRITE(LUPRI,'(/)')
C
C
C     *********************************************
C     *** Print orbitals and occupation numbers ***
C     *********************************************
C
C
C     Get coefficients and CI vector
C     ------------------------------
C
      CALL MEMGET('REAL',KUNO,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KOCCNO,NORBT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KOCCNO),NORBT)
      CALL MEMGET('REAL',KCREF,NZCONFQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C
      CALL REAKRMC(LUKRMC,'NEWORB  ',WORK(KCMO),NCMOTQ)
      IF (NZCONF .GT. 0) 
     &   CALL REAKRMC(LUKRMC,'CREF    ',WORK(KCREF),NZCONFQ)
C
C     Calculate natural orbital occ.
C     ------------------------------
C
C     Also, transform orbitals to NO.
C
      CALL RGETNO(WORK(KCREF),WORK(KOCCNO),WORK(KUNO),WORK(KCMO),
     &     .TRUE.,MOISNO,IPROPT,WORK,KFREE,LFREE)
C
C     Write orbital back to KRMCSCF
C
C     FIXME: -"-
C
      IF ( IPROPT .GE. 2 ) THEN
C
C        Print orbitals.
C
         WRITE(LUPRI,9400)
C
         CALL RPRORB(WORK(KCMO),.TRUE.,.FALSE.)
C        CALL RPRORB(CMO,PROCC,PRSML)
C
 9400    FORMAT(//' (KRMCOUT) Final orbitals'/1X,25('-')/)
      END IF
C
      CALL MEMREL('KRMCOU1.NO',WORK,1,KFRSAV,KFREE,LFREE)
C
C
C     ****************************
C     *** Output of CI vectors ***
C     ****************************
C
C
      IF ( NZCONF .GT. 0 .AND. IPROPT .GE. 1) THEN
         CALL MEMGET('REAL',KCREF,NZCONFQ,WORK,KFREE,LFREE)
         CALL REAKRMC(LUKRMC,'CREF    ',WORK(KCREF),NZCONFQ)
C
         WRITE(LUPRI,9500)
         CALL RPRVEC(NZCONF,WORK(KCREF),1,NZ,IPQTOQ(1,0),
     &        THRPCI,200,LUPRI)
C
         IF (OPT_CIPROGRAM .EQ. 'GOSCI') THEN
            CALL GOSCI_ANACI(WORK(KCREF),WORK(KZCONF),WORK(KFREE),LFREE)
         ELSE IF (OPT_CIPROGRAM .EQ. 'LUCIAREL') THEN
C
C           Analyze CI vectors
C
            CALL LUCI_ANAL(WORK(KFREE),LFREE)
C
         END IF
C
         CALL MEMREL('KRMCOU1.CI',WORK,1,KFRSAV,KFREE,LFREE)
 9500    FORMAT(//' (KRMCOUT) Final CI vector'/1X,25('-')/)
      END IF
C
C
      CALL QEXIT('KRMCOU1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rprorb */
      SUBROUTINE RPRORB(CMO,PROCC,PRSML)
C*****************************************************************************
C
C     Print molecular orbital coefficients.
C
C     Written by J. Thyssen - Mon Mar 26 11:58:53 MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxorb.h"
C
#include "dcblab.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION CMO(*)
      LOGICAL PROCC, PRSML
C 
      DIMENSION CMOB(7)
C
C
C *** FORMAT statements
C
 2400 FORMAT(/5X,'Molecular orbitals for symmetry species',I4)
 2600 FORMAT(/' Orbital  ',4X,7I9)
 2800 FORMAT(I4,2X,A12,7F9.4)
C
C *** Print molecular orbitals
C
      ISTBAS = 0
      DO 400 ISYM = 1,NFSYM
         NORBI = NORB(ISYM)
         IF (PROCC) THEN
            NENDI = NOCC(ISYM)
         ELSE
            NENDI = NORBI
         END IF
      IF (NENDI.EQ.0) GO TO 300
         NBASI = NFBAS(ISYM,0)
         IF (PRSML) THEN
            IBAEND = NBASI
         ELSE
            IBAEND = NFBAS(ISYM,1)
         END IF
         WRITE(LUPRI,2400) ISYM
C
         ICMOI = ICMOQ(ISYM)
         ISTORB = IORB(ISYM)
         IEND=0
  100       IST =IEND+1
            ISTMO=IEND*NBASI+ICMOI
            IEND=IEND+7
            IF(IEND.GT.NENDI) IEND=NENDI
            IEMO=NBASI*(IEND-1)+ICMOI
            WRITE(LUPRI,2600) (I,I=IST,IEND)
            DO 200 I=1,IBAEND
               JSMO=ISTMO+I
               JEMO=IEMO+I
               JJ = 0
               DO J = JSMO, JEMO, NBASI
                  JJ = JJ + 1
                  JJJ = J
                  CMOB(JJ) = CMO(JJJ)**2
                  DO IZ = 2, NZ
                     JJJ = JJJ + NBASI * NORBI
                     CMOB(JJ) = CMOB(JJ) + CMO(JJJ)**2
                  END DO
                  CMOB(JJ) = SQRT(CMOB(JJ))
               END DO
               WRITE(LUPRI,2800) I,PLABEL(IPLAB(I,1),1),
     *            (CMOB(J),J=1,IEND-IST+1)
  200       CONTINUE
         IF (IEND.NE.NENDI) GO TO 100
C
  300 CONTINUE
        ISTBAS = ISTBAS + NBASI
  400 CONTINUE
C
C *** End of subroutine PRORB
C
      RETURN
      END
