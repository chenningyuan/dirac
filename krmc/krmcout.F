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
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE KRMCOUT(WORK,LWORK)
C*****************************************************************************
C
C     Driver to KRMCOU1
C
C     Written by J. Thyssen - Nov 30 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
C
#include <dcbkrmc_itinfo.h>
#include <dcbopt.h>
C
#include <memint.h>
C
      CALL QENTER('KRMCOUT')
C
      CALL MEMGET('REAL',KDINFO,LDINFO*MOPT_MXMACRO,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIINFO,LIINFO*MOPT_MXMACRO,WORK,KFREE,LFREE)
C
      CALL KRMCOU1(WORK(KDINFO),LDINFO,WORK(KIINFO),LIINFO,
     &             WORK(KFREE),LFREE)
C
      CALL QEXIT('KRMCOUT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiopt */
      SUBROUTINE KRMCOU1(DINFO,LDINFO,IINFO,LIINFO,WORK,LWORK)
C*****************************************************************************
C
C     Output from wave function optimization
C
C     Written by J. Thyssen - Nov 30 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbopt.h>
#include <dcborb.h>
#include <frame.h>
C
      DIMENSION WORK(LWORK), DINFO(LDINFO,*), IINFO(LIINFO,*)
C
#include <memint.h>
C
      CALL QENTER('KRMCOU1')
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
         IF (IINFO(2,I) .NE. 0) THEN
            WRITE(LUPRI,9130) IINFO(1,I),IINFO(2,I),
     &           IINFO(3,I),
     &           DINFO(16,I),DINFO(18,I),DINFO(17,I),
     &           (DINFO(17,I)/IINFO(2,I))
         ELSE
            WRITE(LUPRI,9130) IINFO(1,I),IINFO(2,I),
     &           IINFO(3,I),
     &           DINFO(16,I),DINFO(18,I),DINFO(17,I) 
         END IF
      END DO
C
      WRITE(LUPRI,9140) (IINFO(1,I),DINFO(1,I),DINFO(2,I),
     &                   DINFO(3,I),I=1,ITM)
      WRITE(LUPRI,9150) (IINFO(1,I),DINFO(4,I),DINFO(5,I),
     &                   DINFO(6,I),I=1,ITM)
      WRITE(LUPRI,9160) (IINFO(1,I),DINFO(7,I),DINFO(10,I),
     &                   DINFO(8,I),DINFO(9,I),I=1,ITM)

 9100 FORMAT(//,'  ITER ITMIC',T18,'EMCSCF',T35,'GRDNRM',
     &       T49,'RATIO',T60,'STPLNG',/,1X,69('-'),
     &       /,(2I5,F20.12,F15.10,F10.6,F15.10))
 9110 FORMAT(1P,//,'  ITER',T11,'GCINRM',T26,'GOBNRM',T41,'GOEENRM',
     &       T56,'GOEPNRM',T71,'GRDNRM',
     &       /,1X,78('-'),/,(I5,5(D15.5)))
 9120 FORMAT(//,'  ITER ITMIC NZRED      ',T28,'TIMMAC',
     &       T38,'TIMITR',T48,'TIMMIC',T58,'      ',
     &       T68,'TIMMIC/ITMIC',/,1X,78('-'),/)
 9130 FORMAT(I5,2I6,5X,3F10.2,5X,F10.2)
 9140 FORMAT(//,' ITER',T15,'EMY',T35,'EACTIV',T55,'EMCSCF',
     &       //,(I5,3F20.12))
 9150 FORMAT(//,' ITER',T15,'DEPRED',T35,'DEACT',T55,'RATIO',
     &       //,(I5,3F20.12))
 9160 FORMAT(//,' ITER',T10,'BETA',T25,'GAMMA',T43,'STPLNG',
     &       T63,'RTRUST',//,(I5,F16.8,F12.8,2F20.12))

C
C
C
  200 CONTINUE
C
C
      CALL HEADER('Final results from Dirac',-1)
C
      CALL TSTAMP(' ',LUPRI)
C
      WRITE(LUPRI,'(/)')
      IF (NASHT .GT. 1 .OR. IOPT_ISTATE .GT. 1) THEN
C
C        Write out KR-KMCSCF energy
C
         WRITE(LUPRI,9200) IOPT_ISTATE,EMCSCF,POTNUC,EMCSCF-POTNUC
         WRITE(LUPRI,9220) GRDNRM
      ELSE
C
C        Write out KR-KMCSCF energy
C
         WRITE(LUPRI,9210) EMCSCF,POTNUC,EMCSCF-POTNUC
         WRITE(LUPRI,9220) GRDNRM
      END IF
 9200 FORMAT(  T6,'State number         :',T31,I5,
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
      CALL QEXIT('KRMCOU1')
      RETURN
      END
