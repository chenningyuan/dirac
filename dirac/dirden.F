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
C VERSION : $Revision: 1.11.4.1 $
C DATE    : $Date: 2004/09/01 08:08:46 $
C FILE    : dirden.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckden */
      SUBROUTINE FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get backtransformed (contravariant) density matrices.
C     Based on DENSTY
C
C     Written by J. Thyssen
C     Last revision: 1997/09/22 - jth
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
      DIMENSION DCAO(*),DVAO(*),CMO(*),DV(*)
      DIMENSION WORK(*)
      LOGICAL GETDC,GETDV
C
      CALL QENTER('FCKDEN')
      IF ( IPRINT .GE. 9 ) CALL HEADER('Output from FCKDEN',-1)
C
C     Get active density matrix
C
      IF (GETDV) THEN
C         CALL DENSDV(DVAO,CMO,DV,IPRINT,WORK,LWORK)
         CALL GENDEN(DVAO,CMO,2,IPRINT)
         CALL DSCAL(N2BBASXQ,0.50D00,DVAO,1)
      END IF
C
C     Get inactive density matrix
C   
      CALL DZERO(DCAO,N2BBASXQ)
      IF (GETDC) THEN
         CALL DENSTY(DCAO,CMO,IPRINT)
      END IF
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('FCKDEN')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densdv */
      SUBROUTINE DENSDV(DVAO,CMO,DV,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get backtransformed (contravariant) active density.
C     Based on DALTON/FCKDEN
C
C     Written by J. Thyssen
C     Last revision: 1997/09/22 - jth
C                    2002/09/26 - tf : minor cleanup
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION DVAO(*),CMO(*),DV(NASHT,NASHT,NZ)
      DIMENSION WORK(LWORK)
      PARAMETER ( D1 = 1.0 D00, D0 = 0.0 D00 )
C
      CALL QENTER('DENSDV')
#include "memint.h"
C
      IF ( IPRINT .GE. 10 ) CALL HEADER('Output from DENSDV',-1) 
C
C     Get active density matrix
C
      IF (NASHT .EQ. 0) THEN
         WRITE(LUPRI,*) 'WARNING!!! Active density matrix'//
     &                  ' requested but there are no active electrons.'
         WRITE(LUPRI,*) 'So DVAO = 0'
         GOTO 1000
      END IF
      CALL DZERO(DVAO,N2BBASXQ)
C 
C     Allocate memory for temp. DVMO matrix
C
      CALL MEMGET('REAL',KTEMP,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTEMP),N2ORBXQ)
C
C
      CALL DVSCT(DV,WORK(KTEMP))
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('DENSDV "Expanded" DVMO',-1)
         CALL PRQMAT(WORK(KTEMP),NORBT,NORBT,NORBT,NORBT,
     &        NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      DO I = 1,NFSYM
         IF ( NASH(I) .GT. 0 ) THEN
            CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NORB(I),NORB(I),
     &               DVAO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),
     &                   NZ,IPQTOQ(1,0),
     &               WORK(KTEMP+I2ORBX(I,I)),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &               CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &               CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('DENSDV: Active density matrix',-1)
         DO I = 1,NFSYM
            WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
            CALL PRQMAT(DVAO(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
      CALL MEMREL('DENSDV',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
 1000 CONTINUE
      CALL QEXIT('DENSDV')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dgetdv */
      SUBROUTINE DGETDV(DV,NZ)
C*****************************************************************************
C
C     Get active density matrix 
C     This is a dummy routine; it just returns DV=[[1.0]]
C
C     Written by J. Thyssen
C     Last revision: 1997/10/22 - jth
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION DV(*)
C
      DV(1) = 0.50D00
      DO I = 2,NZ
         DV(I) = 0.00D00
      END DO
      RETURN
      END
C
C Dummy routine: writtes matrix to file for use by jth python program
C
      SUBROUTINE PYTHONMATRIX(A,NR,NC,NDIMR,NDIMC,FILENAME)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(NDIMR,NDIMC)
      CHARACTER FILENAME*(*)
      WRITE(6,*) FILENAME
      OPEN(100, FILE = FILENAME, STATUS = 'NEW', FORM = 'FORMATTED' )
      WRITE(100,'(I10)') NR
      WRITE(100,'(I10)') NC
      DO J = 1,NC
         DO I = 1,NR
            WRITE(100,'(D22.12)') A(I,J)
         END DO
      END DO
      CLOSE (100)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck genden */
      SUBROUTINE GENDEN(DMAT,CMO,IFLAG,IPRINT)
C*****************************************************************************
C
C     Get (projection) projection density matrices (normed to unity)
C     IFLAG is a 6 bit pattern:
C       000000: No density matrices
C       000001: Inactive density matrix
C       000010: Active density matrix (see bit 6-?)
C       000100: Secondary electron density matrix
C       001000: Secondary positron density matrix
C       010000: C Ct (that is 1 in AO-basis)
C       100000: QED inactive density matrix: 
C               inactive density matrix + Secondary positron density matrix
C
C     Important: Act. density matrix works only for one open shell.
C
C     Written by J. Thyssen
C     Last revision: 1997/11/23 - jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL LBIT
      DIMENSION DMAT(*), CMO(*)
      CHARACTER CHR*80
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
#include "ibtfun.h"
C
      CALL QENTER('GENDEN')
C
      IOFF = 1
      DO 10 I = 1,6
      IF (LBIT(IFLAG,I)) THEN
        CALL DZERO(DMAT(IOFF),N2BBASXQ)
        DO 20 ISYM = 1,NFSYM
          IF (I.EQ.1) THEN
            ITEST  = NISH(ISYM)
            ISTART = NPSH(ISYM)+1
            ILOOP  = NISH(ISYM)
            CHR    = 'GENDEN: Inactive density matrix, corep '
          ELSEIF (I.EQ.2) THEN
            ITEST  = NASH(ISYM)
            ISTART = NPSH(ISYM)+NISH(ISYM)+1
            ILOOP  = NASH(ISYM)
            CHR    = 'GENDEN: Active density matrix, corep '
          ELSEIF (I.EQ.3) THEN
            ITEST  = NSSH(ISYM)
            ISTART = NPSH(ISYM)+NISH(ISYM)+NASH(ISYM)+1
            ILOOP  = NSSH(ISYM)
            CHR    = 'GENDEN: Sec. (electron) density matrix, corep '
          ELSEIF (I.EQ.4) THEN
            ITEST  = NPSH(ISYM)
            ISTART = 1
            ILOOP  = NPSH(ISYM)
            CHR    = 'GENDEN: Sec. (positron) density matrix, corep '
          ELSEIF (I.EQ.5) THEN
            ITEST  = NORB(ISYM)
            ISTART = 1
            ILOOP  = NORB(ISYM)
            CHR    = 'GENDEN: C Ct (1 in AO-basis), corep '
          ELSEIF (I.EQ.6) THEN
            ITEST  = NPSH(ISYM)+NISH(ISYM)
            ISTART = 1
            ILOOP  = NPSH(ISYM)+NISH(ISYM)
            CHR    = 'GENDEN: QED inactive density matrix, corep '
          END IF
          IF (ITEST.EQ.0) GOTO 20
          CALL DENST1(DMAT(I2BASX(ISYM,ISYM)+IOFF),NTBAS(0),NTBAS(0),
     &              CMO(ICMOQ(ISYM)+1),NFBAS(ISYM,0),NORB(ISYM),
     &              ISTART,ILOOP,NFBAS(ISYM,0))
          IF (IPRINT.GE.10) THEN
            WRITE(LUPRI,'(/3X,A,I2/)') CHR,ISYM
            CALL PRQMAT(DMAT(I2BASX(ISYM,ISYM)+IOFF),
     &                  NFBAS(ISYM,0),NFBAS(ISYM,0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          END IF
 20     CONTINUE           
        IOFF = IOFF + N2BBASXQ
      ENDIF
 10   CONTINUE
C
C   
      CALL QEXIT('GENDEN')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck actden */
      SUBROUTINE ACTDEN(DMAT,CMO,IOPEN,IPRINT)
C*****************************************************************************
C
C     Get active density matrix number IOPEN
C
C     Written by J. Thyssen
C     Last revision: 1997/12/27 - jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DMAT(*), CMO(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      CALL QENTER('ACTDEN')
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 20 ISYM = 1,NFSYM
         ITEST  = NACSH(ISYM,IOPEN)
         ISTART = NPSH(ISYM)+NISH(ISYM)+1
         ILOOP  = NACSH(ISYM,IOPEN)
         DO JOPEN = 1,IOPEN-1
            ISTART = ISTART + NACSH(ISYM,JOPEN)
         END DO
         IF (ITEST.EQ.0) GOTO 20
         CALL DENST1(DMAT(I2BASX(ISYM,ISYM)+1),NTBAS(0),NTBAS(0),
     &              CMO(ICMOQ(ISYM)+1),NFBAS(ISYM,0),NORB(ISYM),
     &              ISTART,ILOOP,NFBAS(ISYM,0))
         IF (IPRINT.GE.10) THEN
            WRITE(LUPRI,'(/3X,A,I2,A,I2/)') 
     &         'ACTDEN: Active density matrix number ',IOPEN,
     &         ', corep ',ISYM
            CALL PRQMAT(DMAT(I2BASX(ISYM,ISYM)+1),
     &                NFBAS(ISYM,0),NFBAS(ISYM,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
         END IF
 20   CONTINUE
C
      CALL QEXIT('ACTDEN')
      RETURN
      END
