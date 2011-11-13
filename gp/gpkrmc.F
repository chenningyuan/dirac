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
C VERSION : $Revision: 1.4.2.1 $
C DATE    : $Date: 2004/09/01 08:10:55 $
C FILE    : gpkrmc.F
C
C Purpose: General Purpose utility routines originally made for krmc,
C          but some are now also used other places (e.g. prp/).
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prdnz3 */
      SUBROUTINE PRDNZ3(AMAT,NDIM1,NDIM2,NZ,IP,LUNIT)
C***********************************************************************
C
C     Print a matrix in the Dirac (NZ,3) format.
C
C     Input:
C       AMAT  - the matrix
C       NDIM1, NDIM2, NZ - dimensions.
C       IP    - index vector from packed quarternion to quaternion.
C       LUNIT - print unit
C
C     Output:
C       None
C
C     Written by J. Thyssen - Nov 24 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NDIM1,NDIM1,NDIM2,NZ,3)
      CHARACTER M(4)*2
      DIMENSION IP(4)
      DATA M /'AR', 'AI', 'BR', 'BI' /
C
      IF (NDIM1 * (NDIM1+1) / 2 .NE. NDIM2) THEN
         WRITE(LUPRI,9000) NDIM1, NDIM2
         CALL QUIT('*** ERROR in PRDNZ3 ***')
      END IF
C
 9000 FORMAT(/,'*** ERROR in PRDNZ3 ***',/,
     &     'Inconsistent dimensions: NDIM1 = ',I5,
     &     ', NDIM2 = ',I5)
C
C     Print matrix.
C
 9001 FORMAT(/3X,'Indices',5X,'Matrix',4X,'Elements',
     &       /1X, 65('-'))
 9010 FORMAT(1P,1X,4I3,4X,A2,6X,3(D14.7,1X))
 9011 FORMAT(1P,17X,A2,6X,3(D14.7,1X))
C
      WRITE(LUNIT,9001)
C
      IJ = 0
      DO I = 1,NDIM1
         DO J = 1,I
            IJ = IJ + 1
            DO K = 1,NDIM1
               DO L = 1,NDIM1
                  WRITE(LUNIT,9010)
     &                 K,L,I,J,
     &                 M(IP(1)),
     &                 (AMAT(K,L,IJ,1,ICL),ICL = 1,3)
                  DO IZ = 2, NZ
                     WRITE(LUNIT,9011) M(IP(IZ)),
     &                    (AMAT(K,L,IJ,IZ,ICL),ICL = 1,3)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rprkap */
      SUBROUTINE RPRKAP(XKAPEE,XKAPEP,JXOPE,NZXOPE,JXOPP,NZXOPP,
     &                  NZ,IP,LUPRI)
C***********************************************************************
C
C     Print Kappa vector
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
C
      DIMENSION XKAPEE(NZXOPE,NZ),XKAPEP(NZXOPP,NZ),
     &          JXOPE(2,NZXOPE),JXOPP(2,NZXOPP),IP(4)
      CHARACTER M(4)*1
      DATA M /'1','i','j','k'/
C
      IF (NZXOPE .EQ. 0) GOTO 100
      WRITE(LUPRI,9000) 'e-e'
      WRITE(LUPRI,9005) (M(IP(I)),I=1,NZ)
      WRITE(LUPRI,9006) ('---------------  ',I=1,NZ)
      DO IG = 1,NZXOPE
         WRITE(LUPRI,9010) IG,JXOPE(1,IG),JXOPE(2,IG),
     &                     (XKAPEE(IG,I),I=1,NZ)
      END DO
      DN = DNORM2(NZXOPE*NZ,XKAPEE,1)
      WRITE(LUPRI,9015) DN
  100 CONTINUE
      IF (NZXOPP .EQ. 0) GOTO 200
      WRITE(LUPRI,9000) 'e-p'
      WRITE(LUPRI,9005) (M(I),I=1,NZ)
      WRITE(LUPRI,9006) ('---------------  ',I=1,NZ)
      DO IG = 1,NZXOPP
         WRITE(LUPRI,9010) IG,JXOPP(1,IG),JXOPP(2,IG),
     &                     (XKAPEP(IG,I),I=1,NZ)
      END DO
      DN = DNORM2(NZXOPP*NZ,XKAPEP,1)
      WRITE(LUPRI,9015) DN
 9000 FORMAT(/,1X,A,' part of matrix:')
 9005 FORMAT(1P,/,1X,'Index(r,s)   r    s   ',4(A1,'-part           '))
 9006 FORMAT(1X, '----------  ---  ---  ',4A17)
 9010 FORMAT(1P,5X,I4,4X,I3,2X,I3,2X,4(D15.8,'  '))
 9015 FORMAT(/,1X,'Norm of vector: ',1P,D15.8)
C
  200 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtkrmc */
      SUBROUTINE WRTKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Write VECTOR to LUKRMC under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "iratdef.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
C
C     Write label at end or - if already existing - after SODLABEL
C
      CALL LABKRMC(LUNIT,LABEL)
C
C     Write vector
C
      LEN = MAX(4,ISIZE)
      CALL WRITSQ(LUNIT,IRAT*LEN,VECTOR)
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck iwrtkrmc */
      SUBROUTINE IWRTKRMC(LUNIT,LABEL,IVECTOR,ISIZE)
C***********************************************************************
C
C     Write VECTOR to LUKRMC under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "iratdef.h"
C
      DIMENSION IVECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
C
C     Write label at end or - if already existing - after SODLABEL
C
      CALL LABKRMC(LUNIT,LABEL)
C
C     Write vector
C
      LEN = MAX(4*IRAT,ISIZE)
      CALL WRITI(LUNIT,LEN,IVECTOR)
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck labkrmc */
      SUBROUTINE LABKRMC(LUNIT,LABEL)
C***********************************************************************
C
C     Write label LABEL to LUNIT;
C       If label LABEL exists:
C         Seek to label SODLABEL and write from there.
C       Else:
C         Seek to label EOFLABEL, backspace and write from there
C       End If
C
C     Input:
C        LABEL   - the label
C
C     Output:
C
C     Written by J. Thyssen - Dec 28 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
      CHARACTER LABEL*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C     Label exists.
C     Seek to label SODLABEL and write from there.
C
         REWIND(LUNIT)
         IF (.NOT. FNDLAB('SODLABEL',LUNIT)) THEN
            CALL QUIT('LABKRMC: wrong structure in file - no SODLABEL')
         END IF
      ELSE
C
C     Label does not exist
C     Seek to label EOFLABEL and write from there
C
         REWIND(LUNIT)
         IF (FNDLAB('EOFLABEL',LUNIT))
     &      BACKSPACE LUNIT
C
      END IF
C
C     Write label
C
      CALL NEWLAB(LABEL,LUNIT,LUPRI)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reakrmc */
      SUBROUTINE REAKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUKRMC under label LABEL.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "iratdef.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         CALL READT(LUNIT,ISIZE,VECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     &        '*** ERROR in REAKRMC ***',
     &        'Label <',LABEL,'> does not exist on unit ',LUNIT
         CALL QUIT('*** ERROR in REAKRMC ***')
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ireakrmc */
      SUBROUTINE IREAKRMC(LUNIT,LABEL,IVECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUKRMC under label LABEL.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "iratdef.h"
C
      DIMENSION IVECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         CALL READI(LUNIT,ISIZE,IVECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     $        '*** ERROR in IREAKRMC ***',
     $        'Label <',LABEL,'> does not exists on unit ',LUNIT
         CALL QUIT('*** ERROR in IREAKRMC ***')
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reakrmc */
      SUBROUTINE REAKRMC2(LUNIT,FILENAME,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUNIT under label LABEL.
C
C     Open & close files if necesary.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "iratdef.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
      CHARACTER FILENAME*(*)
      LOGICAL   FNDLAB, FEXIST, FOPEN
C
C
C     Check if file exists
C
      INQUIRE(FILE=FILENAME,EXIST=FEXIST)
      IF (FEXIST) THEN
         INQUIRE(FILE=FILENAME,OPENED=FOPEN)
         IF (.NOT. FOPEN) THEN
            OPEN(LUNIT,FILE=FILENAME,STATUS='UNKNOWN',
     &           FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         END IF
      ELSE
         FOPEN = .FALSE.
         OPEN(LUNIT,FILE=FILENAME,STATUS='NEW',
     &        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         REWIND(LUNIT)
         CALL NEWLAB('SODLABEL',LUNIT,LUPRI)
      END IF
C
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         LEN = MAX(4,ISIZE)
         CALL READT(LUNIT,ISIZE,VECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     $        '*** ERROR in REAKRMC ***',
     $        'Label <',LABEL,'> does not exists on unit ',LUNIT
         CALL QUIT('*** ERROR in REAKRMC ***')
      END IF
C
      IF (.NOT. FOPEN) CLOSE(LUNIT)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rupkwop */
      SUBROUTINE RUPKWOP(WOPEE,WOPEP,UNWOP,NORBT,
     &                   JXOPE,NZXOPE,JXOPP,NZXOPP,NZ)
C***********************************************************************
C
C     Unpack WOP into UNWOP
C     WOP is a rotation type vector and UNWOP is a square matrix.
C
C     Input:
C        WOPEE   - e-e orbital rotation type vector
C        WOPEP   - e-p orbital rotation type vector
C        JXOPE   - e-e rotations
C        JXOPP   - e-p rotations
C
C     Output:
C        UNWOP   - unpacked matrix
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WOPEE(NZXOPE,NZ),WOPEP(NZXOPP,NZ),
     &          UNWOP(NORBT,NORBT,NZ),JXOPE(2,*),JXOPP(2,*)
C
C     e-e rotations
C
      DO IG = 1,NZXOPE
         K = JXOPE(1,IG)
         L = JXOPE(2,IG)
         UNWOP(K,L,1) =  WOPEE(IG,1)
         UNWOP(L,K,1) = -WOPEE(IG,1)
         DO IZ = 2,NZ
            UNWOP(K,L,IZ) =  WOPEE(IG,IZ)
            UNWOP(L,K,IZ) =  WOPEE(IG,IZ)
         END DO
      END DO
C
C     e-p rotations
C
      DO IG = 1,NZXOPP
         K = JXOPP(1,IG)
         L = JXOPP(2,IG)
         UNWOP(K,L,1) =  WOPEP(IG,1)
         UNWOP(L,K,1) = -WOPEP(IG,1)
         DO IZ = 2,NZ
            UNWOP(K,L,IZ) =  WOPEP(IG,IZ)
            UNWOP(L,K,IZ) =  WOPEP(IG,IZ)
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotcmo */
      SUBROUTINE ROTCMO(CMO,NBASI,NORBI,NZ,XKAPPA,NORBT,IP,IPRINT,
     &                  WORK,LWORK)
C***********************************************************************
C
C     Rotate coefficients
C
C     CMO_{new} = CMO_{old} exp( -\kappa )
C     We use a trick a la Sostrup notes, chapter 3.1.5
C
C     Input:
C        CMO     - old coefficients
C        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
C        IPRINT  - print level
C
C     Output:
C        CMO     - new coefficients
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
C
      PARAMETER (D1 = 1.00 D00, D0 = 0.00 D00, DM1 = -D1)
      PARAMETER (DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00)
      PARAMETER (D1P5 = 1.5D00, D2 = 2.0D00, THREQL = 1.0D-12)
C
      DIMENSION CMO(NBASI,NORBI,NZ), XKAPPA(NORBT,NORBT,NZ),
     &          IP(4), WORK(*)
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12
C
#include "memint.h"
C
      CALL QENTER('ROTCMO')
      CALL GETTIM(CPU1,WALL1)
      IF (IPRINT .GE. 5) CALL HEADER('Output from ROTCMO',-1)
C
#ifndef OLDCODE
C
      NORBXQ = NORBI*NORBI*NZ
      CALL MEMGET('REAL',KEIG,NORBI,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KKSQ,NORBI*NBASI*NZ,WORK,KFREE,LFREE)
      KCOS = KKSQ
      CALL MEMGET('REAL',KSIN,NORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KV  ,NORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP,NORBXQ,WORK,KFREE,LFREE)
C
C     Calculate -\kappa^2 = \kappa^{\dagger} \kappa
C     ---------------------------------------------
C
      CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &           'H','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           'N','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           D0,IP,WORK(KKSQ),NORBI,NORBI,NZ)
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: -kappa^2',-1)
         CALL PRQMAT(WORK(KKSQ),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Diagonalize -\kappa^2
C     ---------------------------------------------
C
      CALL QDIAG(NZ,NORBI,WORK(KKSQ),NORBI,NORBI,
     &           WORK(KEIG),1,WORK(KV),NORBI,NORBI,
     &           WORK(KFREE),LFREE,IERR)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: V matrix (eigenvector)',-1)
         CALL PRQMAT(WORK(KV),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO: eigenvalues',-1)
         CALL PRQMAT(WORK(keig),1,NORBI,1,NORBI,1,IP,LUPRI)
      END IF
C
C     Make cos(d) and sin(d)/d ; d = sqrt(eigenvalues)
C     ------------------------------------------------
C
      CALL GETTIM(CPU2,WALL2)
      CALL DZERO(WORK(KCOS),NORBXQ)
      CALL DZERO(WORK(KSIN),NORBXQ)
      DO I = 1,NORBI
         IADDR = (I-1) + (I-1)*NORBI
         EIG = WORK(KEIG+I-1)
         IF (EIG .LT. D0) THEN
            IF (IPRINT .GE. 3 .AND. ABS(EIG).GT. THRZER)
     &         WRITE(LUPRI,9010) EIG
            EIG = D0
            WORK(KEIG+I-1)   = D0
         END IF
         IF (EIG .LT. DTHRS) THEN
C
C           Use Taylor expansion to avoid division by zero and/or
C           to speed up evaluation (avoid call to sin/cos/sqrt).
C           For EIG < 1.0D-8 the 3. order Taylor expansion is correct
C             to machine precision.
C
            WORK(KCOS+IADDR) = D1-DP5*EIG
            WORK(KSIN+IADDR) = D1-DSIXTH*EIG
         ELSE
            SQEIG = SQRT(EIG)
            WORK(KCOS+IADDR) = COS(SQEIG)
            WORK(KSIN+IADDR) = SIN(SQEIG)/SQEIG
         END IF
      END DO
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: sin(d)/d',-1)
         CALL PRQMAT(WORK(KSIN),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO: cos(d)',-1)
         CALL PRQMAT(WORK(KCOS),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU3,WALL3)
 9010 FORMAT('WARNING from ROTCMO: Negative eigenvalue ',1P,D15.6)
C
C     Make V cos(d) V^{\dagger}
C     (use KTMP as temporatory storage)
C     ------------------------------------
C
      CALL QTRANS('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI,
     &            WORK(KTMP),NORBI,NORBI,NZ,IP,
     &            WORK(KCOS),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KFREE),LFREE,IPROPT)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: V cos(d) V^{dagger}',-1)
         CALL PRQMAT(WORK(KTMP),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Make V sin(d)/d V^{\dagger}
C     (use KCOS as temporatory storage)
C     ------------------------------------
C
      CALL QTRANS('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI,
     &            WORK(KCOS),NORBI,NORBI,NZ,IP,
     &            WORK(KSIN),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KFREE),LFREE,IPROPT)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(WORK(KCOS),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
C     ---------------------------------------------------------
C
      CALL QGEMM(NORBI,NORBI,NORBI,DM1,
     &           'N','N',IP,WORK(KCOS),NORBI,NORBI,NZ,
     &           'N','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           D1,IP,WORK(KTMP),NORBI,NORBI,NZ)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: exp(-kappa)',-1)
         CALL PRQMAT(WORK(KTMP),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Check if unitary
C
c#ifdef KRMC_DEBUG
      CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &           'N','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           'H','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KCOS),NORBI,NORBI,NZ)
      ioff = 0
      i_max = 0
      j_max = 0
      d_max = D0
      do i = 1,norbi
      do j = 1,norbi
         if (i.eq.j) then
            if (work(kcos+ioff)-d1.gt.d_max) then
               i_max = i
               j_max = j
               d_max = work(kcos+ioff)-d1
            end if
         else
            if (work(kcos+ioff).gt.d_max) then
               i_max = i
               j_max = j
               d_max = work(kcos+ioff)
            end if
         end if
         ioff = ioff + 1
      end do
      end do
      if (d_max .gt. d0) then
         write(6,*) 'ROTCMO: maximum deviation from unitary matrix'
         write(6,*) 'Element ',i_max,j_max,' value :',d_max
      end if
C#endif
C
C
C     Calculate new orbitals: CMO_{old} exp(-kappa)
C     ---------------------------------------------
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &           'N','N',IP,CMO,NBASI,NORBI,NZ,
     &           'N','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KCOS),NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,WORK(KCOS),1,CMO,1)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('ROTCMO: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (IPRINT .GE. 2) THEN
         CPUTID = SECTID(CPU4-CPU1)
         WALLTID = SECTID(WALL4-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID
      END IF
 9020 FORMAT(/'ROTCMO: total CPU (WALL) time       : ',A12,' (',A12,')')
 9021 FORMAT( 'ROTCMO: sin/sqrt/cos CPU (WALL) time: ',A12,' (',A12,')')
C
#else
      call quit('this code does not work')
C     fixme: dette ser ikke ud til at virke for complex groups.
C
C     Memory alloc
C
      CALL MEMGET('REAL',KXKOLD,NORBI*NORBI*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXKNEW,NORBI*NORBI*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP  ,NBASI*NBASI*NZ,WORK,KFREE,LFREE)
C
C
C     Calculate exp(-K) = 1 - K
C
      K = 0
      DO IZ = 1,NZ
         DO J = 1,NORBI
            DO I = 1,NORBI
               K = K + 1
               WORK(KXKNEW + K - 1) = -XKAPPA(I,J,IZ)
               IF ( I .EQ. J .AND. IZ .EQ. 1)
     $              WORK(KXKNEW + K - 1 ) = WORK(KXKNEW + K - 1 ) +  D1
            END DO
         END DO
      END DO
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('ROTCMO: 1 - Kappa',-1)
         CALL PRQMAT(WORK(KXKNEW),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      DMX = D0
      ITS = 0
C
C     Symmetric Orthogonalization with N-R iterative algorithm:
C
C     C(p+1) = 1/2( 3C(p) - C(p) S C^{\dagger}(p) C(p))
C
C            = 3/2 C(p) + ( -1/2 C(p) C^{\dagger}(p) ) C(p)
C
C     (here S = unit matrix)
C
 100  CONTINUE
         ITS = ITS + 1
         CALL DCOPY(NORBI*NORBI*NZ,WORK(KXKNEW),1,WORK(KXKOLD),1)
         CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &        'N','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        'H','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        D0,     IP,WORK(KTMP),  NORBI,NORBI,NZ)
         CALL QGEMM(NORBI,NORBI,NORBI,-DP5,
     &        'N','N',IP,WORK(KTMP),  NORBI,NORBI,NZ,
     &        'N','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        D1P5,   IP,WORK(KXKNEW),NORBI,NORBI,NZ)
C
         DMX1 = DMX
         DMX = DABSMAXDIFF(WORK(KXKNEW),WORK(KXKOLD),NORBI*NORBI*NZ)
         IF (DMX .GT. D2 * DMX1 .AND. DMX1 .GT. D0) THEN
            WRITE(LUPRI,'(/1X,2A)')
     &           '*** ERROR in ROTCMO ***: ',
     &           'Divergence in symmetric ortho.'
            CALL QUIT('*** ERROR in ROTCMO: Divergence ***')
         END IF
         IF (DMX .GT. THREQL) GOTO 100
C
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: symm. ortho. rot. mat.',-1)
         CALL PRQMAT(WORK(KXKNEW),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C
C     Calculate new orbitals: CMO_{old} * (symm.ortho.rot.mat)
C     --------------------------------------------------------
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &           'N','N',IP,CMO,NBASI,NORBI,NZ,
     &           'N','N',IP,WORK(KXKNEW),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KTMP),NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,WORK(KTMP),1,CMO,1)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('ROTCMO: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
#endif
      CALL MEMREL('ROTCMO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('ROTCMO')
      RETURN
      END
      DOUBLE PRECISION FUNCTION DABSMAXDIFF(A,B,N)
#include "implicit.h"
      PARAMETER ( D0 = 0.0D00 )
      DOUBLE PRECISION A(N),B(N)
      DMX = D0
      DO I = 1,N
         DEL = ABS(A(I) - B(I))
         DMX = MAX(DMX,DEL)
      END DO
      DABSMAXDIFF = DMX
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdbcidx */
      SUBROUTINE SETDCBIDX()
C***********************************************************************
C
C     Set common block DCBIDX (index vectors, orbital symmetries etc.)
C
C     Written by J. Thyssen - Nov 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbopt.h"
C
C
      IORBI = 0
C
      IESHI = 0
      IPSHI = 0
C
      IISHI = 0
      IASHI = 0
      ISSHI = 0
C
      DO IFSYM = 1, NFSYM
         DO I = 1, NPSH(IFSYM)
            IORBI = IORBI + 1
            IPSHI = IPSHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTPOSI
C
            IDXG2E(IORBI) = -1
            IDXG2P(IORBI) = IPSHI
            IDXG2I(IORBI) = -1
            IDXG2U(IORBI) = -1
            IDXG2V(IORBI) = -1
C
         END DO
         DO I = 1, NISH(IFSYM)
            IORBI = IORBI + 1
            IESHI = IESHI + 1
            IISHI = IISHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTINAC
C
            IDXG2E(IORBI) = IESHI
            IDXG2P(IORBI) = -1
            IDXG2I(IORBI) = IISHI
            IDXG2U(IORBI) = -1
            IDXG2V(IORBI) = -1
C
            IDXE2G(IESHI) = IORBI
            IDXI2G(IISHI) = IORBI
C
         END DO
         DO I = 1, NASH(IFSYM)
            IORBI = IORBI + 1
            IESHI = IESHI + 1
            IASHI = IASHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTACT
C
            IDXG2E(IORBI) = IESHI
            IDXG2P(IORBI) = -1
            IDXG2I(IORBI) = -1
            IDXG2U(IORBI) = IASHI
            IDXG2V(IORBI) = -1
C
            IDXU2G(IASHI) = IORBI
            IDXE2G(IESHI) = IORBI
C
         END DO
         DO I = 1, NSSH(IFSYM)
            IORBI = IORBI + 1
            IESHI = IESHI + 1
            ISSHI = ISSHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTSEC
C
            IDXG2E(IORBI) = IESHI
            IDXG2P(IORBI) = -1
            IDXG2I(IORBI) = -1
            IDXG2U(IORBI) = -1
            IDXG2V(IORBI) = ISSHI
C
            IDXE2G(IESHI) = IORBI
            IDXV2G(ISSHI) = IORBI
         END DO
      END DO
C
C
      IF (IPROPT .GE. 10) THEN
         WRITE(LUPRI,'(1X,A)')
     &        '(SETDCBIDX) I, IFSMO, IOBTYP, IDXG2U, IDXG2E'
         DO I = 1, NORBT
            WRITE(LUPRI,'(8X,I6,I7,I8,2I12)')
     &           I, IFSMO(I), IOBTYP(I), IDXG2U(I), IDXG2E(I)
         END DO
         WRITE(LUPRI,'(1X,A)')
     &        '(SETDCBIDX) E and IDXE2G'
         DO I = 1, NESHT
            WRITE(LUPRI,'(8X,I6,I11)') I, IDXE2G(I)
         END DO
         WRITE(LUPRI,'(1X,A)')
     &        '(SETDCBIDX) U and IDXU2G'
         DO I = 1, NASHT
            WRITE(LUPRI,'(8X,I6,I11)') I, IDXU2G(I)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck rprvec */
      SUBROUTINE RPRVEC(NDIM,VEC,INCVEC,NZ,IPQ,THRESH,MAXLIN,LUOUT)
C
C 19-Aug-1989 hjaaj
C   print VEC(1:NDIM*INCVEC:INCVEC)
C
C   NDIM      : Number of elements in vector VEC
C   VEC(:)    : Vector to be printed
C   INCVEC    : Increment between each element in vector VEC
C   NZ        : real, complex, or quaternion
C   IPQ       : index packed quaternion to quaternion
C   THRESH    : Print threshold for vector with unit norm
C               (if THRESH .lt. 0 then -THRESH is used without
C                renormalization).
C   MAXLIN    : max. lines of output with vector elements
C   LUOUT     : output unit
C
C
#include "implicit.h"
      DIMENSION VEC(NDIM*INCVEC,NZ)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (D1LOW = D1 - 1.0D-10, D1HGH = D1 + 1.0D-10)
      LOGICAL   VSCALE
      DIMENSION IPQ(4)
#include "dcbqal.h"
C
C     Test input
C
      NERR = 0
      IF (NDIM   .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &     'No print from PRVEC because NDIM .le. 0'
         NERR = NERR + 1
      END IF
      IF (INCVEC .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &      'No print from PRVEC because INCVEC .le. 0'
         NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) RETURN
C
C
C
      C2NRM = DNORM2(NDIM*NZ,VEC,INCVEC)
      IF (THRESH .LE. D0) THEN
         THRES2 = -THRESH
         VSCALE = .FALSE.
      ELSE
         THRES2 =  THRESH * C2NRM
         VSCALE = (C2NRM .LT. D1LOW .OR. C2NRM .GT. D1HGH)
      END IF
      IF (VSCALE) THEN
         WRITE (LUOUT,'(/2A,1P,D12.4,/)')
     &      ' Print of vector elements (vector scaled to unit norm)',
     &      ' larger than',ABS(THRESH)
         SCALE = D1 / C2NRM
      ELSE
         WRITE (LUOUT,'(/A,1P,D12.4,/)')
     &      ' Print of vector elements larger than',ABS(THRESH)
         SCALE = D1
      END IF
C
      IPR   = 0
      IZER  = 0
      C2SUM = D0
C
      IF (VSCALE) THEN
         WRITE(LUOUT,9000)
      ELSE
         WRITE(LUOUT,9001)
      END IF
 9000 FORMAT(5X,'Element',3X,'Coefficient (scaled to unit norm)',
     &     /,1X,78('-'))
 9001 FORMAT(5X,'Element',3X,'Coefficient',/,1X,78('-'))
C
      DO 300 I = 1, NDIM
         NA = 1 + (I-1)*INCVEC
         DNNA = D0
         DO IZ = 1, NZ
            DNNA = DNNA + VEC(NA,IZ)**2
         END DO
         IF (SQRT(DNNA).LE.THRES2 .OR. IPR .GE. MAXLIN) THEN
            C2SUM = C2SUM + DNNA
            IF (DNNA .EQ. D0) IZER = IZER + 1
         ELSE
            IF (MOD(IPR,5) .EQ. 0 .AND. IPR .NE. 0) WRITE (LUOUT,'()')
            IPR = IPR + 1
            IF (VSCALE) THEN
               IF (NZ .EQ. 1) THEN
                  WRITE(LUOUT,50)
     &                 I,VEC(NA,1),
     &                 SCALE*VEC(NA,1)
               ELSE IF (NZ .EQ. 2) THEN
                  WRITE(LUOUT,51)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ),
     &                 SCALE*VEC(NA,1),
     &                 (SCALE*VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               ELSE
                  WRITE(LUOUT,52)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ),
     &                 SCALE*VEC(NA,1),
     &                 (SCALE*VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               END IF
            ELSE
               IF (NZ .EQ. 1) THEN
                  WRITE(LUOUT,60)I,VEC(NA,1)
               ELSE IF (NZ .EQ. 2) THEN
                  WRITE(LUOUT,61)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               ELSE
                  WRITE(LUOUT,62)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               END IF
            END IF
 50         FORMAT(1X,I10,3X,1P,D10.2,0P,'(',F10.6,')')
 51         FORMAT(1X,I10,3X,1P,D10.2,2X,D10.2,1X,A1,2X,
     &           0P,'(',F10.6,2X,F10.6,1X,A1,')')
 52         FORMAT(1X,I10,3X,1P,D10.2,2X,3(D10.2,1X,A1,2X),
     &           0P,'(',F10.6,2X,3(F10.6,1X,A1,2X),')')
 60         FORMAT(1X,I10,3X,F20.10)
 61         FORMAT(1X,I10,3X,F20.10,2X,F20.10,1X,A1,2X)
 62         FORMAT(1X,I10,3X,F20.10,2X,3(F20.10,1X,A1,2X))
         END IF
  300 CONTINUE
      IF (IPR .GE. MAXLIN) THEN
C
C     *** We have reached the print limit
C
         WRITE (LUOUT,910) IPR
      END IF
  910 FORMAT(/' Print limit of',I6,' elements has been reached.')
      C2SUM = SQRT(C2SUM)
      IF (IZER .EQ. NDIM) THEN
         WRITE (LUOUT,920) NDIM
      ELSE
         WRITE (LUOUT,930) NDIM,NDIM-IPR,IZER,C2SUM,C2NRM
      END IF
  920 FORMAT(/' Length of vector                      :',I10,
     *       /' All elements are zero.',/)
  930 FORMAT(/' Length of vector                      :',I10,
     *       /' Number of elements not printed        :',I10,
     *       /' Number of zero elements               :',I10,
     *       /' Total norm of coefficients not printed:',F10.6,
     *       /' (the coefficients are normalized to    ',F10.6,')',/)
      RETURN
C
C End of PRVEC
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnorm */
      SUBROUTINE QNORM(S,VC,NBASI,NORBI,NZ,IP,IRETUR,SCR1,SCR2)
C***********************************************************************
C
C     Gram-Schmidt orthogonalize coefficients.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
#include "consts.h"
      PARAMETER (THRRND = 0.9D0, THNORM = 1.0D-6)
C
      DIMENSION S(NBASI,NBASI)
      DIMENSION VC(NBASI,NORBI,NZ)
      DIMENSION IP(*)
      DIMENSION SCR1(NBASI,NZ)
      DIMENSION SCR2(NORBI)
C
      DO 20 I = 1, NORBI
         ITURN = 0
 1       ITURN = ITURN + 1
         IROUND = 0
C
C        Calculate S * VC(I)
C
         CALL QGEMM(NBASI,1,NBASI,D1,
     &        'N','N',IP,S,NBASI,NBASI,1,
     &        'N','N',IP,VC(1,I,1),NBASI,NORBI,NZ,
     &        D0,IP,SCR1,NBASI,1,NZ)
C
         TNORM = D0
         DO IZ = 1, NZ
            TNORM = TNORM + DDOT(NBASI,VC(1,I,IZ),1,SCR1(1,IZ),1)
         END DO
D        write(lupri,*) 'tnorm = ',tnorm
         IF (TNORM .LT. THRZER) THEN
            IRETUR = -1
            RETURN
         END IF
C
C        Normalize vector
C
         TNORM = D1 / SQRT(TNORM)
         DO IZ = 1, NZ
            CALL DSCAL(NBASI,TNORM,VC(1,I,IZ),1)
         END DO
C
C        Calculate S * VC(I)
C
         CALL QGEMM(NBASI,1,NBASI,D1,
     &        'N','N',IP,S,NBASI,NBASI,1,
     &        'N','N',IP,VC(1,I,1),NBASI,NORBI,NZ,
     &        D0,IP,SCR1,NBASI,1,NZ)
C
         TNORM = D0
         DO IZ = 1, NZ
            TNORM = TNORM + DDOT(NBASI,VC(1,I,IZ),1,SCR1(1,IZ),1)
         END DO
C
C        Begin normalization loop
C
         DO J = 1, I - 1
            T = D0
            DO IZ = 1, NZ
               T = T +  DDOT(NBASI,VC(1,J,IZ),1,SCR1(1,IZ),1)
            END DO
            TNORM = TNORM - T * T
            SCR2(J) = -T
         END DO
C
         IF (TNORM .LT. THRZER) THEN
C           ... zero vector after orthogonalization
            IRETUR = -1
            RETURN
         END IF
C
         IF (TNORM .LT. THNORM) IROUND = IROUND + 1
C
         TNORM = D1 / SQRT (TNORM)
         DO J = 1, I - 1
            SCR2(J) = SCR2(J) * TNORM
         END DO
         SCR2(I) = TNORM
C
C        Replace VC(*,I)
C
C        SCR1 = \sum_{k=1,i} SCR2(k) * VC(*,k)
C
         CALL QGEMM(NBASI,1,I,D1,
     &        'N','N',IP,VC,NBASI,NORBI,NZ,
     &        'N','N',IP,SCR2,NORBI,1,1,
     &        D0,IP,SCR1,NBASI,1,NZ)
C
         DO IZ = 1, NZ
            CALL DCOPY(NBASI,SCR1(1,IZ),1,VC(1,I,IZ),1)
         END DO
C
         IF (ITURN .EQ. 1 .AND. IROUND .GT. 0) THEN
D           write(LUPRI,*) 'INFO: second round in QNORM'
            GOTO 1
         END IF
         IF (IROUND .GT. 0) CALL QUIT('QNORM internal error')
 20   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck matgat */
      SUBROUTINE MATGAT(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,NZ)
C***********************************************************************
C
C     Construct BMAT from AMAT elements using IDX.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ), BMAT(NRB,NCB,NZ), IDX(*)
C
      DO IZ = 1, NZ
         DO J = 1, NCB
            JJ = IDX(J)
            DO K = 1, NRB
               KK = IDX(K)
               BMAT(K,J,IZ) = AMAT(KK,JJ,IZ)
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck matsct */
      SUBROUTINE MATSCT(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,NZ)
C***********************************************************************
C
C     Construct BMAT from AMAT elements using IDX.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ), BMAT(NRB,NCB,NZ), IDX(*)
C
      DO IZ = 1, NZ
         DO J = 1, NCB
            JJ = IDX(J)
            DO K = 1, NRB
               KK = IDX(K)
               AMAT(KK,JJ,IZ) = BMAT(K,J,IZ)
            END DO
         END DO
      END DO
      RETURN
      END
