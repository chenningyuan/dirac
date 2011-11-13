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
C VERSION : $Revision: 1.21.2.1 $
C DATE    : $Date: 2004/09/01 08:08:50 $
C FILE    : diropen.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MKMOFK */
      SUBROUTINE MKMOFK(FMO,FAO,DMAT,TMAT,CMO,SAOMO,DTEMP1,DTEMP2,
     &                  WORK,LWORK)
C****************************************************************************
C     PURPOSE:
C        Construct MO Fock matrix from Fock matrices in AO-basis
C
C        Written by J. Thyssen
C        Last revision: Dec 30 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( DP5 = 0.50D00 , D0 = 0.00D00, D1 = 1.00D00,
     &            DM1 = -1.00D00, D2 = 2.00D00, DFTHR = 1.0D-12)
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION FMO(*),FAO(N2BBASXQ,NFMAT),TMAT(*),CMO(*),
     &          DTEMP1(N2BBASXQ,*),DTEMP2(N2BBASXQ),
     &          DMAT(N2BBASXQ,NFMAT),SAOMO(*),WORK(*)
      CHARACTER CTMP*3
C
#include "memint.h"
C
      CALL QENTER('MKMOFK')
C
C     One or more open shells that do NOT couple:
C     -------------------------------------------
C
C     Note: t stands for hermitian conjungation
C     Fmo = FDmo + DVmo * FVmo * (DCmo - DSmo)t + H.C.
C         = FDmo + Ct Sao DVao FVao (beta DCao + gamma DSao)t Sao C + H.C
C     where beta  = alpha*f
C     and   gamma = a-1
C
C     Eg. for a single open shell (NACSH = 1, NAELEC = 1)
C        beta  = 1 and 
C        gamma = -1
C
C     Open shells that do couple:
C     ---------------------------
C
C     See paper by J. Thyssen
C     (reference: manuscript)
C
C        Transform FD to MO-basis
C
      DO I = 1,NFSYM
      IF(NORB(I).GT.0) THEN
        IF(SUB_BL) THEN
          CALL TSUBBL(I,FAO,FMO,TMAT,TMAT,WORK(KFREE),LFREE,IPRSCF)
        ELSE           
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NTMO(I),NTMO(I),
     &            FAO(I2BASX(I,I)+1,1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2TMOT(I)+1),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            TMAT(1+I2TMT(I)),NFBAS(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRSCF)
        ENDIF
        IF(IPRSCF.GE.10) THEN
           WRITE(CTMP,'(I1)') I
           CALL HEADER('MKMOFK: D Fock matrix in MO-basis'//
     &                 ', corep '//CTMP,-1)
           CALL PRQMAT(FMO(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &                 NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
      END IF
      END DO
C
C     If not .trivec. then we don't have coefficients to form a density 
C     matrix, so return.
C
      IF (.NOT. TRIVEC) THEN
C
C        In next iteration we'll have some coefficients
C
         TRIVEC = .TRUE.
         GOTO 30
      END IF
C
C
C     Construct coupling operator L
C
C
      DO 10 IOPEN = 1,NOPEN
C
         BETA  = DALPHA(IOPEN)*DF(IOPEN)
         GAMMA = DA(IOPEN)-D1
         IF (IPRSCF.GE.5) THEN
            WRITE(LUPRI,'(A,I1,A,F10.6)') 
     &         'MKMOFK: BETA  factor, open shell ',IOPEN,' : ',BETA
            WRITE(LUPRI,'(A,I1,A,F10.6)') 
     &         'MKMOFK: GAMMA factor, open shell ',IOPEN,' : ',GAMMA
         END IF
C
C
C        Get secondary density matrices as projection operator
C           in AO-basis
C
         CALL GENDEN(DTEMP1,CMO,12,IPRSCF)
C
C        Construct total secondary density matrix DSao
C
         CALL DAXPY(N2BBASXQ,D1,DTEMP1(1,2),1,DTEMP1(1,1),1)
C
C        Construct the projection operator:
C         BETA DCao + GAMMA DSao
C           = alpha*f DCao + (a_o-1)*DSao
C
         CALL DSCAL(N2BBASXQ,GAMMA,DTEMP1(1,1),1)
         CALL DAXPY(N2BBASXQ,BETA,DMAT(1,1),1,DTEMP1(1,1),1)
C
C
         IOFFSAOMO = 1
         DO 20 I = 1,NFSYM
C
            IF (NACSH(I,IOPEN).GT.0) THEN
C
C              QVao (beta DCao + gamma DSao) (store in DTEMP1(1,2))
C
C              QVao = 2*f/2 * Fock(TWOINT) = f * Fock(TWOINT)
C
               CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),DF(IOPEN),
     &            'N','N',IPQTOQ(1,0),FAO(1+I2BASX(I,I),1+IOPEN),
     &            NTBAS(0),NTBAS(0),NZ,
     &            'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),1),
     &            NTBAS(0),NTBAS(0),NZ,
     &            D0,IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &            NTBAS(0),NTBAS(0),NZ)
C
C              Active-Active elements
C
               DO 25 JOPEN = IOPEN+1,NOPEN
                  IF (NACSH(I,JOPEN).EQ.0) GOTO 25
C
C                 We use DTEMP1(xxx,3) as temporary storage
C
C                 fi = fj?
C
                  CALL DZERO(DTEMP1(1,3),N2BBASXQ)
C
                  IF (ABS(DF(IOPEN)-DF(JOPEN)).LT.DFTHR) THEN
C
C                    This is a tough case...
C                    The gradient is correct, but there might
C                    be serious trouble converging the wave-function,
C                    because differences of diagonal elements are not
C                    propertional to the orbital hessian.
C
C                    If fi = fj = 1 then it corresponds to inactive shells
C                    so don't make the correction below.
C
                     IF (ABS(DF(IOPEN)-D1).GE.DFTHR) THEN
C
                        CALL DAXPY(N2BBASXQ,-D1,FAO,1,DTEMP1(1,3),1)
C
C                       add two-electron Fock matrices...
C
C                       ( + (ai-1) fi QViao
C                         - (aj-1) fj QVjao ) DVjao
C
                        DSCALE = DF(IOPEN)*
     &                           (DA(IOPEN)-D1)
C
C                       QViao = 2*fi * 1/2 * Fock(TWOINT)
C
                        DSCALE = DSCALE * DF(IOPEN)
C
                        CALL DAXPY(N2BBASXQ,DSCALE,
     &                     FAO(1,1+IOPEN),1,DTEMP1(1,3),1)
                        DSCALE = -DF(JOPEN)*
     &                           (DA(JOPEN)-D1)
C
C                       QVjao = 2*fj * 1/2 * Fock(TWOINT)
C
                        DSCALE = DSCALE * DF(JOPEN)
C
                        CALL DAXPY(N2BBASXQ,DSCALE,
     &                     FAO(1,1+JOPEN),1,DTEMP1(1,3),1)
                        CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &                     'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),3),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     'N','N',IPQTOQ(1,0),
     &                        DMAT(1+I2BASX(I,I),1+JOPEN),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     D1,IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &                     NTBAS(0),NTBAS(0),NZ)
C
                     END IF
C
C           
                  ELSE
C
C                    ( + (ai-1)/(fi-fj) fi QViao
C                      - (aj-1)/(fi-fj) fj QVjao ) DVjao
C
                     DSCALE = DF(IOPEN)*
     &                        (DA(IOPEN)-D1)/(DF(IOPEN)-DF(JOPEN))
C
C                    QViao = 2*fi * 1/2 * Fock(TWOINT)
C
                     DSCALE = DSCALE * DF(IOPEN)
C
                     CALL DAXPY(N2BBASXQ,DSCALE,
     &                  FAO(1,1+IOPEN),1,DTEMP1(1,3),1)
                     DSCALE = -DF(JOPEN)*
     &                        (DA(JOPEN)-D1)/(DF(IOPEN)-DF(JOPEN))
C
C                    QVjao = 2*fj * 1/2 * Fock(TWOINT)
C
                     DSCALE = DSCALE * DF(JOPEN)
C
                     CALL DAXPY(N2BBASXQ,DSCALE,
     &                  FAO(1,1+JOPEN),1,DTEMP1(1,3),1)
                     CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &                  'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),3),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  'N','N',IPQTOQ(1,0),DMAT(1+I2BASX(I,I),1+JOPEN),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  D1,IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &                  NTBAS(0),NTBAS(0),NZ)
C
                  END IF
C
 25            CONTINUE
C
C              Multiply the lot with DViao (store in DTEMP2)
C
C
               CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &            'N','N',IPQTOQ(1,0),DMAT(1+I2BASX(I,I),1+IOPEN),
     &            NTBAS(0),NTBAS(0),NZ,
     &            'N','N',IPQTOQ(1,0),DTEMP1(1+I2BASX(I,I),2),
     &            NTBAS(0),NTBAS(0),NZ,
     &            D0,IPQTOQ(1,0),DTEMP2(1+I2BASX(I,I)),
     &            NTBAS(0),NTBAS(0),NZ)
C              
               IF ( IPRSCF .GE. 15 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: coupling operator in AO-basis'//
     &               ', corep '//CTMP,-1)
                  CALL PRQMAT(DTEMP2(1+I2BASX(I,I)),
     &               NFBAS(I,0),NFBAS(I,0),
     &               NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Transform with halftransformed Saomo from left and right
C              Store the result in DTEMP1(1,1)
C
               CALL QTRANS('AOMO','S',D0,
     &           NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &           DTEMP2(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DTEMP1(1+I2ORBT(I),1),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &           SAOMO(IOFFSAOMO),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &           SAOMO(IOFFSAOMO),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRSCF)
               IF ( IPRSCF .GE. 15 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: Correction matrix in MO-basis'//
     &               ' step 1, corep '//CTMP,-1)
                  CALL PRQMAT(DTEMP1(1+I2ORBT(I),1),NORB(I),NORB(I),
     &               NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Calculate A + At
C
               CALL SUMMMH(DTEMP1(1+I2ORBT(I),1),NORB(I),NZ,
     &                  NORB(I),NORB(I))
               IF ( IPRSCF .GE. 10 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: Correction matrix in MO-basis'//
     &               ' step 2, corep '//CTMP,-1)
                  CALL PRQMAT(DTEMP1(1+I2ORBT(I),1),NORB(I),NORB(I),
     &               NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
C
C              Sum up FD + A
C
               CALL QAXPY(N2ORB(I),'N','N',D1,IPQTOQ(1,0),1,
     &            DTEMP1(1+I2ORBT(I),1),N2ORB(I),1,
     &            IPQTOQ(1,0),NZ,
     &            FMO(1+I2ORBT(I)),N2ORB(I),1,
     &            IPQTOQ(1,0),NZ)
               IF ( IPRSCF .GE. 10 ) THEN
                  WRITE(CTMP,'(I1)') I
                  CALL HEADER('MKMOFK: Total Fock matrix in MO-basis'//
     &               ', corep '//CTMP,-1)
                  CALL PRQMAT(FMO(1+I2ORBT(I)),NORB(I),NORB(I),
     &               NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
               END IF
            END IF
            IOFFSAOMO = IOFFSAOMO + NZ*NORB(I)*NFBAS(I,0)
 20      CONTINUE
 10   CONTINUE
C
 30   CONTINUE
      CALL QEXIT('MKMOFK')
      CALL FLSHFO(LUPRI)
      RETURN
      END
