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
#if defined (VAR_SINGLE)
#include "single.h"
#endif
C
C  /* Deck ericsm */
      SUBROUTINE ERICSM(SO,CSMAT,IPNTCR,INDXBT,WORK,LWORK,IPRINT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "eridst.h"
#include "eribuf.h"
      DIMENSION SO(NAOINT,*), CSMAT(*), IPNTCR(MAXBCH,4), INDXBT(*),
     &          WORK(LWORK) 
#include "cbieri.h"
#include "ericom.h"
#include "symmet.h"
#include "hertop.h"
#ifndef PRG_DIRAC
#include "inftap.h"
#else
#include "hrunit.h"
#endif
C
C     Allocation for ERIOUT
C
      LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
      KNDORB = 1 
      KLAST  = KNDORB + (2*LBIN - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIOUT',' ',KLAST,LWORK)
      CALL ERICS1(SO,CSMAT,WORK(KNDORB),IPNTCR,INDXBT,IPRINT)
      RETURN
      END
C  /* Deck erics1 */
      SUBROUTINE ERICS1(SO,CSMAT,INDORB,IPNTCR,INDXBT,IPRINT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER (D0 = 0.0D0, TEN14 = 1.1D14)
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,*),
     &          INDORB(NCCS,2), IPNTCR(MAXBCH,4),
     &          IPNRST(0:7,3), INDXBT(MXSHEL*MXCONT,0:7),
     &          IADCMP(MXAQN,MXAQN,2), CSMAT(NBAST,NBAST)
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "aobtch.h"
#include "hertop.h"
#include "inforb.h"
#include "symmet.h"
C
#include "ibtfun.h"
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERICSM',-1)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      NBITX = 16
      IBITA = 2**NBITX
      IBITB = 1
      IBITC = 2**NBITX
      IBITD = 1
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
         WRITE (LUPRI,'(2X,A,4I5)') 
     &                    ' NCTFX         ',NCTFA,NCTFB,NCTFC,NCTFD 
         WRITE (LUPRI,'(2X,A,2I5)') ' NCCS2  ',NCCS2
      END IF
C
C
      IF (NCCS2 .GT. 0) THEN
C
         DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            C = A
            D = B
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 200
            IF (DIAGCD .AND. D.GT.C) GO TO 200
C
            IADR = 0
            DO N = 1, NPQBCS
               NA = KNDXBT(IPNTCR(N,1)) - 1
               NB = KNDXBT(IPNTCR(N,2)) - 1
               NC = KNDXBT(IPNTCR(N,3)) - 1
               ND = KNDXBT(IPNTCR(N,4)) - 1
               DO I = 1, NCTFA 
               DO J = 1, NCTFB 
               DO K = 1, NCTFC 
               DO L = 1, NCTFD 
                  IADR = IADR + 1
                  INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                           + (INDXBT(NB + J,B)-1)
                  INDORB(IADR,2) = (INDXBT(NC + K,C)-1)*IBITC
     &                           + (INDXBT(ND + L,D)-1)
               END DO 
               END DO 
               END DO 
               END DO
            END DO
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
C
            IA   = 0
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB = IAB + IBITB
                  IC = 0
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  DO ICMPC = 1, ICMPA
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC  = IC + IBITC
                     ICD = IC
                     IF (ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ICD  = ICD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IF (ICMPAB .EQ. ICMPCD) THEN
                           DO I = 1, NCCS
                              KAB = INDORB(I,1) + IAB
                              KCD = INDORB(I,2) + ICD
                              IF (KAB.EQ.KCD) THEN
                                 KA = IBTAND(IBTSHR(KAB,NBITX),IBITA-1)
                                 KB = IBTAND(       KAB,       IBITA-1)
                                 CSI = SQRT(SO(I,R,S,T,ICMPAB,ICMPCD,1))
                                 CSMAT(KA,KB) = CSI
                                 CSMAT(KB,KA) = CSI
                              END IF
                           END DO
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  200       CONTINUE
C
         END IF
         END DO
         END IF
         END DO
      END IF
C
      RETURN
      END
C  /* Deck ericso */
      SUBROUTINE ERICSO(SO,ITYPE,INDORB,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,INT,LBIN,NIBUF,INDXBT,NBITS,IPRINT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER (D0 = 0.0D0, TEN14 = 1.1D14)
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,*),
     &          INDORB(NCCS,NIBUF),
     &          IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          BIN(LBIN), IBIN(LBIN,NIBUF),
     &          IPNRST(0:7,3), INDXBT(MXSHEL*MXCONT,0:7),
     &          IADCMP(MXAQN,MXAQN,2)
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "aobtch.h"
#include "hertop.h"
#include "symmet.h"
C
#include "ibtfun.h"
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERICSO',-1)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      NBITX = 16
      IBITA = 2**NBITX
      IBITB = 1
      IBITC = 2**NBITX
      IBITD = 1
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
         WRITE (LUPRI,'(2X,A,4I5)') 
     &                    ' NCTFX         ',NCTFA,NCTFB,NCTFC,NCTFD 
         WRITE (LUPRI,'(2X,A,2I5)') ' NCCS2  ',NCCS2
         WRITE (LUPRI,'(2X,A,2I5)') ' NIBUF, NBITX  ',NIBUF, NBITX 
      END IF
C
      INT = 0
C
      IF (NCCS2 .GT. 0) THEN
C
         DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            C = A
            D = B
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 200
            IF (DIAGCD .AND. D.GT.C) GO TO 200
C
            IADR = 0
            DO N = 1, NPQBCS
               NA = KNDXBT(IPNTCR(N,1)) - 1
               NB = KNDXBT(IPNTCR(N,2)) - 1
               NC = KNDXBT(IPNTCR(N,3)) - 1
               ND = KNDXBT(IPNTCR(N,4)) - 1
               DO I = 1, NCTFA 
               DO J = 1, NCTFB 
               DO K = 1, NCTFC 
               DO L = 1, NCTFD 
                  IADR = IADR + 1
                  INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                           + (INDXBT(NB + J,B)-1)
                  INDORB(IADR,2) = (INDXBT(NC + K,C)-1)*IBITC
     &                           + (INDXBT(ND + L,D)-1)
               END DO 
               END DO 
               END DO 
               END DO
            END DO
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
C
            IA   = 0
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB = IAB + IBITB
                  IC = 0
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  DO ICMPC = 1, ICMPA
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC  = IC + IBITC
                     ICD = IC
                     IF (ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ICD  = ICD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IF (ICMPAB .EQ. ICMPCD) THEN
                           DO I = 1, NCCS
                              KAB = INDORB(I,1) + IAB
                              KCD = INDORB(I,2) + ICD
C                             KA = IBTAND(IBTSHR(KAB,NBITX),IBITA-1)
C                             KB = IBTAND(       KAB,       IBITA-1)
C                             KC = IBTAND(IBTSHR(KCD,NBITX),IBITA-1)
C                             KD = IBTAND(       KCD,       IBITA-1)
C                             IJ = MAX(KA,KB)*IBITA + KA + KB
C                             KL = MAX(KC,KD)*IBITA + KC + KD
C                             IF (IJ.EQ.KL) THEN
                              IF (KAB.EQ.KCD) THEN
                                 CSI = SQRT(SO(I,R,S,T,ICMPAB,ICMPCD,1))
                                 IF (ABS(CSI) .GT. THRSH) THEN
                                    INT = INT + 1
                                    BIN (INT  ) = CSI 
                                    IBIN(INT,1) = KAB 
                                 END IF
                              END IF
                           END DO
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  200       CONTINUE
C
         END IF
         END DO
         END IF
         END DO
      END IF
C
      RETURN
      END
C  /* Deck csmext */
      SUBROUTINE CSMEXT(FMAT,NDIMF,NMATF,RMAT,NDIMR,IPNT,IPRINT)
C
C     October 2004 tuh
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.D0)
      DIMENSION FMAT(NDIMF,NDIMF,NMATF), RMAT(NDIMR,NDIMR), IPNT(NDIMR)
C
      IF (IPRINT.GT.10) THEN
         CALL HEADER('Full matrix FMAT',-1)
         DO N = 1, NMATF
            CALL OUTPUT(FMAT(1,1,N),1,NDIMF,1,NDIMF,NDIMF,NDIMF,1,LUPRI)
         END DO
      END IF
C
      IOFF = 0
      DO I = 1, NDIMR
         JOFF = 0
         DO J = 1, NDIMR
            AMAX = D0
            DO N = 1, NMATF
            DO L = IOFF + 1, IOFF + IPNT(I)
            DO K = JOFF + 1, JOFF + IPNT(J)
               AMAX = MAX(ABS(FMAT(K,L,N)),AMAX)
            END DO
            END DO
            END DO
            RMAT(I,J) = AMAX
            JOFF = JOFF + IPNT(J)
         END DO
         IOFF = IOFF + IPNT(I)
      END DO
C
      IF (IPRINT.GT.5) THEN
         CALL HEADER('Reduced matrix',-1)
         CALL OUTPUT(RMAT,1,NDIMR,1,NDIMR,NDIMR,NDIMR,1,LUPRI)
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C aspg, 2006-04-19, merging latest eri code into dirac
C
C replacing csmext for its new incarnation, old code follows
C when i confirm its ok to delete it, it goes away...
C
CC  /* Deck csmext */
C      SUBROUTINE CSMEXT(CSMAT)
CC
CC     October 2004 tuh
CC
C#include "implicit.h"
C#include "priunit.h"
C      PARAMETER (D0 = 0.D0)
C#include "maxaqn.h"
C#include "mxcent.h"
C#include "maxorb.h"
C#include "aovec.h"
C      DIMENSION CSMAT(NBAST,NBAST), CSMSHL(KMAX,KMAX)
C#include "inforb.h"
C#include "shells.h"
C      IOFF = 0
C      DO ISHELL = 1, KMAX
C         JOFF = 0
C         DO JSHELL = 1, KMAX
C            CSMAX = D0
C            DO ICMP = IOFF + 1, IOFF + KHKT(ISHELL)
C               DO JCMP = JOFF + 1, JOFF + KHKT(JSHELL)
C                  CSMAX = MAX(CSMAT(ICMP,JCMP),CSMAX)
C               END DO
C            END DO
C            CSMSHL(ISHELL,JSHELL) = CSMAX
C            JOFF = JOFF + KHKT(JSHELL)
C         END DO
C         IOFF = IOFF + KHKT(ISHELL)
C      END DO
C      IF (.FALSE.) THEN
C         CALL HEADER('Subroutine CSMAT',-1) 
C         CALL OUTPUT(CSMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C         CALL HEADER('Subroutine CSMSHL',-1)
C         CALL OUTPUT(CSMSHL,1,KMAX,1,KMAX,KMAX,KMAX,1,LUPRI)
C      END IF
C      RETURN
C      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
