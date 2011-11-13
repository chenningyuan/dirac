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
C FILE    : gpsaue.F
C
C
#if defined (PRG_DIRAC)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      BLOCK DATA BDQALGEBRA
C
C Initialize /QALGEBRA/
C
#include "implicit.h"
#include "dcbqal.h"
C
      DATA IQMULT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1,
     &              2, 1, 4, 3, 1, 2, 3, 4, 4, 3, 2, 1, 3, 4, 1, 2,
     &              3, 4, 1, 2, 4, 3, 2, 1, 1, 2, 3, 4, 2, 1, 4, 3,
     &              4, 3, 2, 1, 3, 4, 1, 2, 2, 1, 4, 3, 1, 2, 3, 4/
      DATA IQPHASE/ 1, 1, 1, 1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1, 1,-1,
     &              1,-1,-1, 1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1, 1, 1,
     &              1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1, 1,
     &              1,-1, 1,-1, 1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1/
      DATA IQSIGN / 1, 1, 1, 1, 1,-1,-1,-1,
     &              1, 1,-1,-1, 1,-1, 1, 1,
     &              1,-1, 1,-1, 1, 1,-1, 1,
     &              1,-1,-1, 1, 1, 1, 1,-1/
      DATA QUNIT  / '1','i','j','k'/
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prqmat */
      SUBROUTINE PRQMAT(QMAT,NROW,NCOL,LRQ,LCQ,NZ,IQP,IUNIT)
#include "implicit.h"
      DIMENSION QMAT(LRQ,LCQ,NZ),IQP(*)
      CHARACTER M(4)*6
      M(1) = 'A real'
      M(2) = 'A imag'
      M(3) = 'B real'
      M(4) = 'B imag'
      DO 10 IZ = 1,NZ
        IQ = IQP(IZ)
        WRITE(IUNIT,'(/A,A6,A,2I4/)')
     &  '*** ',M(IQ),' part ***;    NROW,NCOL:',NROW,NCOL
        CALL OUTPUT(QMAT(1,1,IZ),1,NROW,1,NCOL,LRQ,LCQ,1,IUNIT)
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck iqpack */
      SUBROUTINE IQPACK(IQPA,NZA,IQPB,NZB,IQPC,NZC)
C***********************************************************************
C
C     IQPA and IQPB are pointers to quaternion units 1,i,j and k
C     IQPACK considers what quaternion units span the quaternion product
C     A * B and makes packs these in pointer IQPC as generated.
C
C     Written by T.Saue - June 27 1996
C     Last revision: June 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "dcbqal.h"
C
      LOGICAL LBUF(4)
      DIMENSION IQPA(NZA),IQPB(NZB),IQPC(*)
      EXTERNAL BDQALGEBRA
C
C     Initialize LBUF
C     ===============
C
      LBUF(1) = .FALSE.
      LBUF(2) = .FALSE.
      LBUF(3) = .FALSE.
      LBUF(4) = .FALSE.
C
C     Determine quaternion units that contribute to A*B
C     =================================================
C
      NZC = 0
      DO IZA = 1,NZA
        IQA = IQPA(IZA)
        DO IZB = 1,NZB
          IQB = IQPB(IZB)
          IQC = IQMULT(IQA,IQB,1)
          IF(.NOT.LBUF(IQC)) THEN
            NZC = NZC + 1
            IQPC(NZC) = IQC
            LBUF(IQC) = .TRUE.
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qgemm */
      SUBROUTINE QGEMM(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA,
     &                             FB,TB,IQPB,B,LRB,LCB,NZB,
     &                              BETA,IQPC,C,LRC,LCC,NZC)
C*****************************************************************************
C
C     Performs matrix-matrix operations on general matrices:
C
C     C := alpha * op(A) * op(B) + beta * C
C
C     where alpha and beta are real constants. A,B and C
C     are matrices that may be real(NZ = 1), complex(NZ = 2)
C     or quaternionic (NZ = 4). Calls BLAS routine DGEMM.
C
C     LRX and LCX are leading rows and columns of matrix X.
C     op(A) is an M by K matrix, op(B) is a K by N matrix
C     and C an M by N matrix.
C
C     op(X) is determined by the corresponding character variable
C     FX:
C       FX = 'N'        Normal matrix
C       FX = 'T'        Transpose..
C       FX = 'C'        Complex conjugate
C       FX = 'H'        Hermitian conjugate
C     and TX
C       TX = 'N'        No transform    a+bj
C       TX = 'I'        i-transform   -i(a+bj)i = a - bj
C       TX = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TX = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     AFUL and BFUL indicates on input what components of
C     matrices A and B are non-zero, whereas CFUL on output
C     indicates what components of C are non-zero.
C
C     Written by T.Saue, University of Oslo, Nov. 1994
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D1=1.0D0)
C
C     Global variables
C
      DIMENSION A(LRA*LCA,NZA),B(LRB*LCB,NZB),C(LRC*LCC,NZC),
     &          IQPA(NZA),IQPB(NZB),IQPC(NZC)
      CHARACTER  FA*1,FB*1,TA*1,TB*1,OA*1,OB*1
C
C     Local variables
C
      DIMENSION IPB(4)
      LOGICAL   LBUF(4)      
C
      EXTERNAL BDQALGEBRA
C
C     Initialize LBUF
C     ===============
C
      CALL LSET(4,.FALSE.,LBUF)
C
C     Make pointer IPB
C     ================
C
      CALL IZERO(IPB,4)
      DO IZB = 1,NZB
        IPB(IQPB(IZB)) = IZB
      ENDDO
C
C     MATRIX A
C

C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FA.EQ.'N') THEN
        OA  = 'N'
        IFA = 1
C     Transpose...
      ELSEIF(FA.EQ.'T') THEN
        OA  = 'T'
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        OA  = 'N'
        IFA = 2
C     Hermitian conjugate...
      ELSEIF(FA.EQ.'H') THEN
        OA  = 'T'
        IFA = 2
      ELSE
        CALL QUIT('QGEMM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
C     MATRIX B
C
C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FB.EQ.'N') THEN
        OB  = 'N'
        IFB = 1
C     Transpose...
      ELSEIF(FB.EQ.'T') THEN
        OB  = 'T'
        IFB = 1
C     Complex conjugate....
      ELSEIF(FB.EQ.'C') THEN
        OB  = 'N'
        IFB = 2
C     Hermitian conjugate...
      ELSEIF(FB.EQ.'H') THEN
        OB  = 'T'
        IFB = 2
      ELSE
        CALL QUIT('QGEMM:Unknown format '//FA//' of matrix B!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TB.EQ.'N') THEN
        ITB = 1
C     i-transform
      ELSEIF(TB.EQ.'I') THEN
        ITB = 2
C     j-transform
      ELSEIF(TB.EQ.'J') THEN
        ITB = 3
C     k-transform
      ELSEIF(TB.EQ.'K') THEN
        ITB = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TB '//TB//' of matrix B!')
      ENDIF
C
      DO 10 IZC = 1,NZC
        IQC = IQPC(IZC)
        DO 20 IZA = 1,NZA
          IQA = IQPA(IZA)
          IQB = IQMULT(IQA,IQC,1)
          IZB = IPB(IQB)
          IF(IZB.EQ.0) GOTO 20
          FACA = ALPHA*(IQSIGN(IQA,IFA,ITA)*IQSIGN(IQB,IFB,ITB)
     &                *IQPHASE(IQA,IQB,1))
          IF(LBUF(IZC)) THEN
            FACB = D1
          ELSE
            FACB = BETA
          ENDIF
          CALL DGEMM(OA,OB,M,N,K,FACA,A(1,IZA),LRA,B(1,IZB),LRB,
     &               FACB,C(1,IZC),LRC)
          LBUF(IZC) = .TRUE.
   20   CONTINUE
   10 CONTINUE
C
C     Check that matrix C is filled
C
      DO IZC = 1,NZC
      IF(.NOT.LBUF(IZC)) 
     &    WRITE(LUPRI,'(A,I5)') '* WARNING(QGEMM): Not filled: ',IZC
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qgeold */
      SUBROUTINE QGEOLD(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA,
     &                             FB,TB,IQPB,B,LRB,LCB,NZB,
     &                              BETA,IQPC,C,LRC,LCC,NZC)
C*****************************************************************************
C
C     Performs matrix-matrix operations on general matrices:
C
C     C := alpha * op(A) * op(B) + beta * C
C
C     where alpha and beta are real constants. A,B and C
C     are matrices that may be real(NZ = 1), complex(NZ = 2)
C     or quaternionic (NZ = 4). Calls BLAS routine DGEMM.
C
C     LRX and LCX are leading rows and columns of matrix X.
C     op(A) is an M by K matrix, op(B) is a K by N matrix
C     and C an M by N matrix.
C
C     op(X) is determined by the corresponding character variable
C     FX:
C       FX = 'N'        Normal matrix
C       FX = 'T'        Transpose..
C       FX = 'C'        Complex conjugate
C       FX = 'H'        Hermitian conjugate
C     and TX
C       TX = 'N'        No transform    a+bj
C       TX = 'I'        i-transform   -i(a+bj)i = a - bj
C       TX = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TX = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     AFUL and BFUL indicates on input what components of
C     matrices A and B are non-zero, whereas CFUL on output
C     indicates what components of C are non-zero.
C
C     Written by T.Saue, University of Oslo, Nov. 1994
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D1=1.0D0)
C
C     Global variables
C
      DIMENSION A(LRA*LCA,NZA),B(LRB*LCB,NZB),C(LRC*LCC,NZC),
     &          IQPA(NZA),IQPB(NZB),IQPC(NZC)
      CHARACTER  FA*1,FB*1,TA*1,TB*1,OA*1,OB*1
C
C     Local variables
C
      LOGICAL   LBUF(4)      
      DIMENSION IPC(4)
C
      EXTERNAL BDQALGEBRA
C
C
C     Initialize LBUF
C     ===============
C
      CALL LSET(4,.FALSE.,LBUF)
C
C     Make pointer IPC
C     ================
C
      CALL IZERO(IPC,4)
      DO IZC = 1,NZC
        IPC(IQPC(IZC)) = IZC
      ENDDO
      IF(NZA.GT.NZC) CALL QUIT('QGEMM:NZA.GT.NZC')
      IF(NZB.GT.NZC) CALL QUIT('QGEMM:NZB.GT.NZC')
C
C     MATRIX A
C

C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FA.EQ.'N') THEN
        OA  = 'N'
        IFA = 1
C     Transpose...
      ELSEIF(FA.EQ.'T') THEN
        OA  = 'T'
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        OA  = 'N'
        IFA = 2
C     Hermitian conjugate...
      ELSEIF(FA.EQ.'H') THEN
        OA  = 'T'
        IFA = 2
      ELSE
        CALL QUIT('QGEMM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
C     MATRIX B
C
C     Determine FX:
C     =============
C     Normal matrix...
      IF    (FB.EQ.'N') THEN
        OB  = 'N'
        IFB = 1
C     Transpose...
      ELSEIF(FB.EQ.'T') THEN
        OB  = 'T'
        IFB = 1
C     Complex conjugate....
      ELSEIF(FB.EQ.'C') THEN
        OB  = 'N'
        IFB = 2
C     Hermitian conjugate...
      ELSEIF(FB.EQ.'H') THEN
        OB  = 'T'
        IFB = 2
      ELSE
        CALL QUIT('QGEMM:Unknown format '//FA//' of matrix B!')
      ENDIF
C     Determine TX:
C     =============
C     No transform
      IF    (TB.EQ.'N') THEN
        ITB = 1
C     i-transform
      ELSEIF(TB.EQ.'I') THEN
        ITB = 2
C     j-transform
      ELSEIF(TB.EQ.'J') THEN
        ITB = 3
C     k-transform
      ELSEIF(TB.EQ.'K') THEN
        ITB = 4
      ELSE
        CALL QUIT('QGEMM:Unknown TB '//TB//' of matrix B!')
      ENDIF
C
      DO 10 IZA = 1,NZA
        IQA = IQPA(IZA)
        DO 20 IZB = 1,NZB
          IQB = IQPB(IZB)
          IQC = IQMULT(IQA,IQB,1)
          IZC = IPC(IQC)
          IF(IZC.EQ.0) THEN
            WRITE(LUPRI,'(A)') 'QGEMM: Error in IPC !'
            WRITE(LUPRI,'(A,4I5)') 'IQPC: ',(IQPC(I),I=1,NZC)
            WRITE(LUPRI,'(A,4I5)') 'IPC : ',(IPC(I), I=1,4)
            CALL QUIT('QGEMM: Error in IPC !')
          ENDIF
          FACA = ALPHA*(IQSIGN(IQA,IFA,ITA)*IQSIGN(IQB,IFB,ITB)
     &                *IQPHASE(IQA,IQB,1))
          IF(LBUF(IZC)) THEN
            FACB = D1
          ELSE
            FACB = BETA
          ENDIF
          CALL DGEMM(OA,OB,M,N,K,FACA,A(1,IZA),LRA,B(1,IZB),LRB,
     &               FACB,C(1,IZC),LRC)
          LBUF(IZC) = .TRUE.
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qdiag */
      SUBROUTINE QDIAG(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,
     &                 WORK,LWORK,IERR)
C***********************************************************************
C
C       This is a driver routine for the diagonalization of
C         * real symmetric matrices                     NZ = 1
C         * complex Hermitian matrices                  NZ = 2
C         * quaternionic Hermitian matrices             NZ = 4
C
C       INPUT:
C
C       A       matrix to be diagonalized
C
C       Control parameter MATZ:
C               MATZ = 0        Only eigenvalues desired.
C               MATZ = 1        Eigenvalues and eigenvectors
C
C       OUTPUT:
C
C       EIG     eigenvalues in ascending order.
C       VEC     eigenvectors when MATZ=1.
C
C       Error parameter IERR:
C               IERR = 0        Normal completion
C               IERR.NE.0       Erroneous completion
C
C       TEMPORARY STORAGE:      FV1,FV2,TAU
C       Written by T.Saue November 1994 - Odense
C       Last revision : Nov 15 1994 - tsaue
C***********************************************************************
C
#include "implicit.h"
      DIMENSION A(LRA,LCA,NZ),EIG(N),VEC(LRV,LCV,NZ),WORK(LWORK)
      CALL QENTER('QDIAG')
#include "memint.h"
      IF((NZ.EQ.1.OR.NZ.EQ.2).AND.(LRA.NE.LRV)) THEN
        CALL QUIT('QDIAG: LRA.ne.LRV')
      ENDIF
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFV1,N   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFV2,N   ,WORK,KFREE,LFREE)
      IF(NZ.GT.1) THEN
        NTAU = N*NZ
        CALL MEMGET('REAL',KTAU,NTAU,WORK,KFREE,LFREE)
      ENDIF
C
C       Householder transform to real symmetric tridiagonal matrix
C       ==========================================================
C
      IF (NZ.EQ.4) THEN
        CALL  QHTRID(A,N,LRA,LCA,EIG,
     &               WORK(KFV1),WORK(KFV2),WORK(KTAU))
        IF(MATZ.NE.0) GO TO 10
      ELSEIF(NZ.EQ.2) THEN
        CALL  HTRIDI(LRA,N,A(1,1,1),A(1,1,2),EIG,
     &               WORK(KFV1),WORK(KFV2),WORK(KTAU))
        IF(MATZ.NE.0) GO TO 10
      ELSEIF(NZ.EQ.1) THEN
        IF(MATZ.EQ.0) THEN
          CALL  TRED1(LRA,N,A(1,1,1),EIG,WORK(KFV1),WORK(KFV2))
        ELSE
          CALL  TRED2(LRA,N,A(1,1,1),EIG,WORK(KFV1),VEC(1,1,1))
          GOTO 20
        ENDIF
      ELSE
        CALL QUIT('QDIAG: Illegal value of NZ !')
      ENDIF
C
C       FIND EIGENVALUES ONLY
C       =====================
C
      CALL  TQLRAT(N,EIG,WORK(KFV2),IERR)
      GO TO 30
C
C       FIND BOTH EIGENVALUES AND EIGENVECTORS
C       ======================================
C
   10 CONTINUE
      CALL DUNIT2(VEC,N,LRV,LCV,NZ)
C
C       Find eigenvalues/vectors of real symmetric
C       tridiagonal matrix
C
   20 CONTINUE
      CALL  TQL2(LRV,N,EIG,WORK(KFV1),VEC(1,1,1),IERR)
      IF (IERR.NE.0) GO TO 30
C
C       Backtransform to find eigenvectors of
C       quaternionic Hermitian matrix
C
      IF    (NZ.EQ.4) THEN
        CALL  QHTRBK(A,N,LRA,LCA,WORK(KTAU),VEC,N,LRV,LCV)
      ELSEIF(NZ.EQ.2) THEN
        CALL  HTRIBK(LRA,N,A(1,1,1),A(1,1,2),WORK(KTAU),N,
     &               VEC(1,1,1),VEC(1,1,2))
      ENDIF
   30 CONTINUE
C
C     Memory deallocation
C
      CALL MEMREL('QDIAG.QDIAG',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('QDIAG')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qmgath */
      SUBROUTINE QMGATH(AIN,NRW1,NCL1,AOUT,NRW2,NCL2,NZ,
     &                  NROW,NCOL,IROW,ICOL)
C***********************************************************************
C
C     Gather operations for generally quaternionic matrices
C
C     Written by T.Saue November 1994
C     Last revision: Nov 15 1994
C***********************************************************************
#include "implicit.h"
      DIMENSION AIN(NRW1,NCL1,NZ),AOUT(NRW2,NCL2,NZ)
      DIMENSION IROW(*),ICOL(*)
C
      DO 10 IZ = 1,NZ
        DO 20 J = 1, NCOL
          JJ = ICOL(J)
          DO 30 I = 1, NROW
           AOUT(I,J,IZ) = AIN(IROW(I),JJ,IZ)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qhtrid */
      SUBROUTINE QHTRID(A,N,LRA,LCA,D,E,E2,TAU)
C
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
      DIMENSION A(LRA,LCA,4),D(N),E(N),E2(N),TAU(4,N)
C***********************************************************************
C

C       Trond Saue, University of Tromsoe, Norway , December 1993
C
C       This subroutine is a quaternionic analogue of the
C       EISPACK routine HTRIDI. QHTRID reduces a quaternionic
C       Hermitian matrix to a realsymmetric tridiagonal matrix.
C       The unitary transformation is done in two steps:
C
C       1. Reduction to quaternionic Hermitian tridiagonal
C       matrix using (N-2) Householder matrices.
C
C       2.Reduction to real symmetric tridiagonal matrix
C       using a quaternionic diagonal matrix TAU.
C
C       Householder matrices has the general form:
C
C       P = 1 - (uu(+))/H       ;       H = u(+)u/2
C
C       INPUT:
C
C       [AR+iAJ+jAJ+kAK]        is a quaternionic Hermitian matrix of
C                               order N and with leading dimension LRA.
C
C       OUTPUT:
C
C       The lower triangle of A contains the row vectors u
C       from the (n-2) Householder matrices. The strict upper
C       triangle and the diagonal of A is unaltered.
C
C       D       diagonal elements of the tridiagonal matrix
C       E       off-diagonal elements of tridiagonal matrix
C       E2      squares of elements in E
C       TAU     the diagonal of TAU
C     ------------------------------------------------------------------
C
      TAU(1,N) = D1
      TAU(2,N) = D0
      TAU(3,N) = D0
      TAU(4,N) = D0
C
C       Store diagonal elements of A temporarily in D
C
      DO 100 I = 1,N
        D(I) = A(I,I,1)
  100 CONTINUE
C
C     For I = N step -1 until 1 DO
C
      DO 300 II = 1,N
         I = N + 1 - II
         L = I - 1
         H = D0
         SCALE = D0
         IF (L.LT.1) GO TO 130
C
C     Evaluate scale factor SCALE
C
         DO 120 K = 1, L
            SCALE = SCALE + ABS(A(I,K,1)) + ABS(A(I,K,2))
     +                    + ABS(A(I,K,3)) + ABS(A(I,K,4))
  120    CONTINUE
C
         IF (SCALE.NE.D0) GO TO 140
C
C     SCALE = 0: Skip transformation
C
         TAU(1,L) = D1
         TAU(2,L) = D0
         TAU(3,L) = D0
         TAU(4,L) = D0
  130    CONTINUE
         E(I)  = D0
         E2(I) = D0
         GO TO 290
C
C     SCALE.NE.0 : Do transformation
C
  140    CONTINUE
         SCLINV = D1/SCALE
C
C     Scale row and evaluate its norm
C
         DO 150 K = 1, L
            A(I,K,1) = SCLINV*A(I,K,1)
            A(I,K,2) = SCLINV*A(I,K,2)
            A(I,K,3) = SCLINV*A(I,K,3)
            A(I,K,4) = SCLINV*A(I,K,4)
            H = H+A(I,K,1)*A(I,K,1)+A(I,K,2)*A(I,K,2)
     +           +A(I,K,3)*A(I,K,3)+A(I,K,4)*A(I,K,4)
  150    CONTINUE
C
         GR    = SQRT(H)
         E(I)  = SCALE*GR
         E2(I) = SCALE*SCALE*H
C
C     Form next diagonal element of matrix T: TAU(L) = - TAU(I)*A(I,L)/F
C     where F is the absolute value of A(I,L)
C
         FR = PYTHAG(PYTHAG(A(I,L,1),A(I,L,2)),
     &               PYTHAG(A(I,L,3),A(I,L,4)))
         IF (FR.EQ.D0) GO TO 160
         FRINV = D1/FR
         TAU(1,L) = FRINV*(TAU(4,I)*A(I,L,4)+TAU(3,I)*A(I,L,3)
     +                    +TAU(2,I)*A(I,L,2)-TAU(1,I)*A(I,L,1))
         SI       = FRINV*(TAU(1,I)*A(I,L,2)+TAU(2,I)*A(I,L,1)
     +                    +TAU(3,I)*A(I,L,4)-TAU(4,I)*A(I,L,3))
         SJ       = FRINV*(TAU(1,I)*A(I,L,3)-TAU(2,I)*A(I,L,4)
     +                    +TAU(3,I)*A(I,L,1)+TAU(4,I)*A(I,L,2))
         SK       = FRINV*(TAU(1,I)*A(I,L,4)+TAU(2,I)*A(I,L,3)
     +                    -TAU(3,I)*A(I,L,2)+TAU(4,I)*A(I,L,1))
         IF (L .EQ. 1) THEN
           TAU(1,L) = -TAU(1,L)
           SI       = -SI
           SJ       = -SJ
           SK       = -SK
           GO TO 270
         ENDIF
C
C     Form element L in vector U
C        U(L) = A(I,L) + A(I,L)[SQRT(H)/QVAL(A(I,L))]
C

         H       = H  + FR   *GR
         GR      = D1 + FRINV*GR
         A(I,L,1) = GR * A(I,L,1)
         A(I,L,2) = GR * A(I,L,2)
         A(I,L,3) = GR * A(I,L,3)
         A(I,L,4) = GR * A(I,L,4)
C         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    CONTINUE
         TAU(1,L) = -TAU(1,I)
         SI       =  TAU(2,I)
         SJ       =  TAU(3,I)
         SK       =  TAU(4,I)
         A(I,L,1)  = GR
  170    CONTINUE
         FR = D0
C
         DO 240 J = 1, L
            GR = D0
            GI = D0
            GJ = D0
            GK = D0
C
C     Form element of A*U
C       (A*U)(j) is stored in (GR,GI,GJ,GK)
C
            DO 180 K = 1, J
               GR = GR+A(J,K,1)*A(I,K,1)+A(J,K,2)*A(I,K,2)
     +                +A(J,K,3)*A(I,K,3)+A(J,K,4)*A(I,K,4)
               GI = GI-A(J,K,1)*A(I,K,2)+A(J,K,2)*A(I,K,1)
     +                -A(J,K,3)*A(I,K,4)+A(J,K,4)*A(I,K,3)
               GJ = GJ-A(J,K,1)*A(I,K,3)+A(J,K,2)*A(I,K,4)
     +                +A(J,K,3)*A(I,K,1)-A(J,K,4)*A(I,K,2)
               GK = GK-A(J,K,1)*A(I,K,4)-A(J,K,2)*A(I,K,3)
     +                +A(J,K,3)*A(I,K,2)+A(J,K,4)*A(I,K,1)
  180       CONTINUE
C
            JP1 = J+1
            IF(L.LT.JP1) GO TO 220
C
C     Exploiting Hermitian symmetry to avoid using upper triangle
C
            DO 200 K = JP1, L
               GR = GR+A(K,J,1)*A(I,K,1)-A(K,J,2)*A(I,K,2)
     +                -A(K,J,3)*A(I,K,3)-A(K,J,4)*A(I,K,4)
               GI = GI-A(K,J,1)*A(I,K,2)-A(K,J,2)*A(I,K,1)
     +                +A(K,J,3)*A(I,K,4)-A(K,J,4)*A(I,K,3)
               GJ = GJ-A(K,J,1)*A(I,K,3)-A(K,J,2)*A(I,K,4)
     +                -A(K,J,3)*A(I,K,1)+A(K,J,4)*A(I,K,2)
               GK = GK-A(K,J,1)*A(I,K,4)+A(K,J,2)*A(I,K,3)
     +                -A(K,J,3)*A(I,K,2)-A(K,J,4)*A(I,K,1)
  200       CONTINUE
C
C     Form element of P
C        P(j) = (A*U)(j)/H
C             - stored in (E(j),TAU(2,j),TAU(3,j),TAU(4,j))
C
  220       CONTINUE
            E(J)     = GR/H
            TAU(2,J) = GI/H
            TAU(3,J) = GJ/H
            TAU(4,J) = GK/H
            FR       = FR+A(I,J,1)*E(J) -A(I,J,2)*TAU(2,J)
     +               - A(I,J,3)*TAU(3,J)-A(I,J,4)*TAU(4,J)
  240    CONTINUE
C
C     Evaluate HH = SUM[j](U*(j)P(j)/H)
C
         HH=FR/(H+H)
C
C     Form reduced A
C
C     Computational formula: A'(j,k) = A(j,k)-Q(j)U*(k)-U(j)Q*(k)
C        where Q(j) = P(j) - HH*U(j)
C
         DO 260 J = 1,L
C
C     Q(j) is stored in (GR,GI,GJ,GK)
C     Q*(k) is stored in (E(k),TAU(2,j),TAU(3,j),TAU(4,j))
C     U(j) is stored in (FR,FI,FJ,FK)
C     U*(k) is stored in (AR(i,k),AI(i,k),AJ(i,k),AK(i,k))
C
            FR       =  A(I,J,1)
            GR       =  E(J)     - HH*FR
            E(J)     =  GR
            FI       = -A(I,J,2)
            GI       =  TAU(2,J) - HH*FI
            TAU(2,J) = -GI
            FJ       = -A(I,J,3)
            GJ       =  TAU(3,J) - HH*FJ
            TAU(3,J) = -GJ
            FK       = -A(I,J,4)
            GK       =  TAU(4,J) - HH*FK
            TAU(4,J) = -GK
C
            DO 260 K = 1, J
              A(J,K,1) = A(J,K,1)
     +            -FR*E(K)+FI*TAU(2,K)+FJ*TAU(3,K)+FK*TAU(4,K)
     +            -GR*A(I,K,1)+GI*A(I,K,2)+GJ*A(I,K,3)+GK*A(I,K,4)
              A(J,K,2) = A(J,K,2)
     +            -FR*TAU(2,K)-FI*E(K)-FJ*TAU(4,K)+FK*TAU(3,K)
     +            -GR*A(I,K,2)-GI*A(I,K,1)-GJ*A(I,K,4)+GK*A(I,K,3)
              A(J,K,3) = A(J,K,3)
     +            -FR*TAU(3,K)+FI*TAU(4,K)-FJ*E(K)-FK*TAU(2,K)
     +            -GR*A(I,K,3)+GI*A(I,K,4)-GJ*A(I,K,1)-GK*A(I,K,2)
              A(J,K,4) = A(J,K,4)
     +            -FR*TAU(4,K)-FI*TAU(3,K)+FJ*TAU(2,K)-FK*E(K)
     +            -GR*A(I,K,4)-GI*A(I,K,3)+GJ*A(I,K,2)-GK*A(I,K,1)
  260    CONTINUE
C
  270    CONTINUE
         DO 280 K = 1,L
            A(I,K,1) = SCALE*A(I,K,1)
            A(I,K,2) = SCALE*A(I,K,2)
            A(I,K,3) = SCALE*A(I,K,3)
            A(I,K,4) = SCALE*A(I,K,4)
  280    CONTINUE
C
         TAU(2,L) = -SI
         TAU(3,L) = -SJ
         TAU(4,L) = -SK
  290    CONTINUE
         HH      = D(I)
         D(I)    = A(I,I,1)
         A(I,I,1) = HH
         A(I,I,2) = SCALE*SQRT(H)
  300 CONTINUE
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qhtrbk */
      SUBROUTINE QHTRBK(A,N,LRA,LCA,TAU,Z,M,LRZ,LCZ)
C
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION A(LRA,LCA,4),TAU(4,N),Z(LRZ,LCZ,4)
C***********************************************************************
C
C       Trond Saue, University of Tromsoe, Norway , December 1993
C
C       This subroutine is a quaternionic analogue of the
C       EISPLACK routine HTRIBK. It forms the eigenvectors
C       of a quaternionic Hermitian matrix by backtransforming
C       those of the corresponding real symmetric matrix
C       determined by QHTRID.
C
C       INPUT:
C
C       [AR+iAJ+jAJ+kAK]        is a quaternionic Hermitian matrix of
C                               order N and with leading dimension LRA.
C                               The lower triangle contains the
C                               u-vectors of the (N-2) Householder
C                               matrices used in the reduction of
C                               A to a quaternionic Hermitian tridiagonal
C                               matrix. The diagonal and upper triangle
C                               of the original matrix is preserved.
C       TAU                     TAU is the diagonal quaternionic elements
C                               of the diagonal matrix transforming
C                               the quaternionic Hermitian tridiagonal
C                               matrix to a real tridiagonal matrix.
C       OUTPUT:
C
C       [ZR+iZI+jZJ+kZK]        eigenvectors
C
C     Note that the last component of each returned vector
C     is real and that vector Euclidean norms are preserved.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C
C     Transform the eigenvector of the real symmetric tridiagonal matrix
C     to those of the quaternionic Hermitian tridiagonal matrix
C
      DO 50 K = 1, N
         DO 50 J = 1, M
            Z(K,J,4) = -TAU(4,K)*Z(K,J,1)
            Z(K,J,3) = -TAU(3,K)*Z(K,J,1)
            Z(K,J,2) = -TAU(2,K)*Z(K,J,1)
            Z(K,J,1) =  TAU(1,K)*Z(K,J,1)
   50 CONTINUE
C
      IF(N.EQ.1) GO TO 200
C
C     Recover and apply the Householder matrices
C
      DO 140 I = 3,N
         L = I-1
         H = A(I,I,2)
         IF (H.EQ.D0) GO TO 140
C
         DO 130 J = 1,M
            SR = D0
            SI = D0
            SJ = D0
            SK = D0
            DO 110 K = 1,L
               SR = SR+A(I,K,1)*Z(K,J,1)-A(I,K,2)*Z(K,J,2)
     +                -A(I,K,3)*Z(K,J,3)-A(I,K,4)*Z(K,J,4)
               SI = SI+A(I,K,1)*Z(K,J,2)+A(I,K,2)*Z(K,J,1)
     +                +A(I,K,3)*Z(K,J,4)-A(I,K,4)*Z(K,J,3)
               SJ = SJ+A(I,K,1)*Z(K,J,3)-A(I,K,2)*Z(K,J,4)
     +                +A(I,K,3)*Z(K,J,1)+A(I,K,4)*Z(K,J,2)
               SK = SK+A(I,K,1)*Z(K,J,4)+A(I,K,2)*Z(K,J,3)
     +                -A(I,K,3)*Z(K,J,2)+A(I,K,4)*Z(K,J,1)
  110       CONTINUE
C
C     Double divisions to avoid possible underflow
C
            SR = (SR/H)/H
            SI = (SI/H)/H
            SJ = (SJ/H)/H
            SK = (SK/H)/H
C
C     Do second part
C
            DO 120 K = 1, L
               Z(K,J,1) = Z(K,J,1)-SR*A(I,K,1)-SI*A(I,K,2)
     +                            -SJ*A(I,K,3)-SK*A(I,K,4)
               Z(K,J,2) = Z(K,J,2)+SR*A(I,K,2)-SI*A(I,K,1)
     +                            -SJ*A(I,K,4)+SK*A(I,K,3)
               Z(K,J,3) = Z(K,J,3)+SR*A(I,K,3)+SI*A(I,K,4)
     +                            -SJ*A(I,K,1)-SK*A(I,K,2)
               Z(K,J,4) = Z(K,J,4)+SR*A(I,K,4)-SI*A(I,K,3)
     +                            +SJ*A(I,K,2)-SK*A(I,K,1)
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
  200 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qopvec */
      SUBROUTINE QOPVEC(FA,TA,IQPA,A,N,NZA)
C***********************************************************************
C
C     Modify quaternion vector as follows:
C
C     FA:
C       FA = 'N'        Normal array
C       FA = 'C'        Complex conjugate
C     and TA
C       TA = 'N'        No transform    a+bj
C       TA = 'I'        i-transform   -i(a+bj)i = a - bj
C       TA = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TA = 'K'        k-transform   -k(a+bj)k = a*-b*j
C
C     Written by T.Saue July 1996
C     Last revision: July 14 196 - tsaue
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D1 = 1.0D0)
C
      CHARACTER FA*1, TA*1
      DIMENSION A(N,NZA),IQPA(*)
C
      EXTERNAL BDQALGEBRA
C
C     Determine FA:
C     =============
C     Normal ARRAY...
      IF    (FA.EQ.'N') THEN
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        IFA = 2
      ELSE
        CALL QUIT('QVECOP:Unknown FA '//FA//' of matrix A!')
      ENDIF
C
C     Determine TA:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QVECOP:Unknown TA '//TA//' of matrix A!')
      ENDIF
C
      DO IZA = 1,NZA
        IQA = IQPA(IZA)
        FACA = IQSIGN(IQA,IFA,ITA)
        IF(FACA.NE.D1) CALL DSCAL(N,FACA,A(1,IZA),1)
      ENDDO      
C
      RETURN
      END
#endif
C  /* Deck ran1 */
      FUNCTION RAN1(IDUM)
#include "implicit.h"
      PARAMETER (IA=16807,IM=2147483647,AM=1.0D0/IM,
     +           IQ=127773,IR=2386,NTAB=32,NDIV=1+(IM-1)/NTAB,
     +           EPS=1.0D-7,RNMX=1.0D0-EPS)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     "Minimal" random number generator of Park and Miller with
C     Bays-Durham shuffle and added safeguards. Return a uniform random
C     deviate between 0.0 and 1.0 (exclusive of the end point values).
C     Call with IDUM a negative integer to initialize; thereafter,
C     do not alter IDUM between succesive deviates in a sequence. RNMX
C     should approximate the largest floating value that is less than 1.
C
C     References:
C       Press,Teukolsky,Vetterling & Flannery: 
C         "Numerical Recipes in FORTRAN",
C         Cambridge University Press 2nd edition 1992
C       Park & Miller (1988), Communications of the ACM, vol.31,pp.1192-1201
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER IV(NTAB)
      SAVE IV,IY
      DATA IV/NTAB*0/, IY /0/
C
C     Initialization:
C        Be sure to prevent IDUM = 0
C        Load the shuffle table(after 8 warm-ups)
C
      IF(IDUM.LE.0.OR.IY.EQ.0) THEN
        IDUM = MAX(-IDUM,1)
        DO 10 J = NTAB+8,1,-1
          K = IDUM/IQ
          IDUM = IA*(IDUM-K*IQ)-IR*K
          IF(IDUM.LT.0) IDUM = IDUM+IM
          IF(J.LE.NTAB) IV(J) = IDUM
   10   CONTINUE
        IY = IV(1)
      ENDIF
C
C     Start here when not initializaing
C
      K = IDUM/IQ
C
C     Compute IDUM = MOD(IA*IDUM,IM) without overflows by
C     Schrage's method
C
      IDUM = IA*(IDUM-K*IQ)-IR*K
      IF(IDUM.LT.0) IDUM = IDUM+IM
      J = 1+IY/NDIV
      IY = IV(J)
      IV(J) = IDUM
      RAN1 = MIN(AM*IY,RNMX)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ranelm */
      FUNCTION RANELM(RMIN,RMAX,IDUM)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Pick a random number in the range <RMIN,RMAX>.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
      R = RAN1(IDUM)
      RANELM = RMIN + (RMAX-RMIN)*R
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ranmat */
      SUBROUTINE RANMAT(TYP,A,LDA,N,RMIN,RMAX,IDUM)
C**********************************************************************
C
C	Generation of test matrices with random elements with
C	values in the range <RMIN,RMAX>.
C	  TYP = S	Symmetric matrix	
C	  TYP = A	Antisymmetric matrix
C	  TYP = G	General matrix
C	On initialization of random number generation IDUM should
C	be a negative number
C**********************************************************************
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      CHARACTER TYP*1
      DIMENSION A(LDA,N)
      IF(TYP.EQ.'G') THEN
        DO 10 J = 1,N
          DO 20 I = 1,N
            A(I,J) = RANELM(RMIN,RMAX,IDUM)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(TYP.EQ.'S') THEN
        DO 30 J = 1,N
          DO 40 I = (J+1),N
            A(I,J) = RANELM(RMIN,RMAX,IDUM)
            A(J,I) = A(I,J)
   40     CONTINUE
          A(J,J) = RANELM(RMIN,RMAX,IDUM)
   30   CONTINUE
      ELSEIF(TYP.EQ.'A') THEN
        DO 50 J = 1,N
          DO 60 I = (J+1),N
            A(I,J) = RANELM(RMIN,RMAX,IDUM)
            A(J,I) = -A(I,J)
   60     CONTINUE
        A(J,J) = D0
   50   CONTINUE
      ELSE
        CALL QUIT('RANMAT: UNKNOWN MATRIX TYPE !')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fulmat */
      SUBROUTINE FULMAT(MATTYP,LDA,NDIM,A)
#include "implicit.h"
#include "priunit.h"
      CHARACTER MATTYP*1
      DIMENSION A(LDA,NDIM)
C*SYMM  : Real symmetric matrix
C*ASYM  : Real antisymmetric matrix
C
C   ON INPUT: A - lower triangular matrix
C   ON OUTPUT:A - full matrix
      IF(MATTYP.EQ.'S') THEN
      DO 10 J = 1,NDIM
        DO 20 I = J,NDIM
          A(I,J) = A(I,J) + A(J,I)
          A(J,I) = A(I,J)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(MATTYP.EQ.'A') THEN
      DO 30 J = 1,NDIM
        DO 40 I = J,NDIM
          A(I,J) = A(I,J) - A(J,I)
          A(J,I) = -A(I,J)
   40     CONTINUE
   30   CONTINUE
      ELSE
        WRITE(LUPRI,'(A,A6,A)') '** ERROR in FULMAT **  : Keyword ',
     +    MATTYP,' not recognized!'
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck offtst */
      FUNCTION OFFTST(A,N,IZ)
C
C	Checks the absolute value of off-diagonal elements
C	of a matrix
C
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
      DIMENSION A(N,N,IZ)
      DIA = D0
      ELM = D0
      DIATST = D0
      OFFTST = D0
      DO 10 K = 1,IZ
        DO 20 J = 1,N
          DO 30 I = 1,N
            ELM = ELM + D1
            OFFTST = ((ELM-D1)*OFFTST + ABS(A(I,J,K)))/ELM
   30     CONTINUE
   20   CONTINUE
        DO 40 J = 1,N
          DIA = DIA + D1
          DIATST = ((DIA-D1)*DIATST + ABS(A(J,J,K)))/DIA
   40   CONTINUE
   10 CONTINUE
      OFFTST = (ELM*OFFTST - DIA*DIATST)/(ELM-DIA)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mxform */
      FUNCTION MXFORM(X,NDIG)
C***********************************************************************
C
C	Find the format that gives optimum precision for real
C	variable X
C
C***********************************************************************
#include "implicit.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER MXFORM*6
C
      MXFORM = '      '
      IF(X.EQ.D0) THEN
        NDEC = NDIG
      ELSE
        NDEC = INT(LOG10(ABS(X)))
      ENDIF
      IF((NDEC.LT.0).OR.(NDEC.GT.(NDIG-3))) THEN
        MXFORM(1:1) = 'D'
        NDEC = NDIG - 7
      ELSE
        MXFORM(1:1) = 'F'
        NDEC = NDIG - NDEC - 3
      ENDIF
      IF(NDEC.LE.0) CALL QUIT('MXFORM: NDEC.LE.0')
      IF(NDIG.LT.10) THEN
        IDIG = ICHAR('0') + NDIG
        MXFORM(2:2) = CHAR(IDIG)
        IOFF= 0
      ELSE
        IDIG = ICHAR('0') + INT(NDIG/10)
        MXFORM(2:2) = CHAR(IDIG)
        IDIG = ICHAR('0') + MOD(NDIG,10)
        MXFORM(3:3) = CHAR(IDIG)
        IOFF = 1
      ENDIF
      IND = 3+IOFF
      MXFORM(IND:IND) = '.'
      IF(NDEC.LT.10) THEN
        IDIG = ICHAR('0') + NDEC
        IND  = 4+IOFF
        MXFORM(IND:IND) = CHAR(IDIG)
        IOFF = IOFF+ 0
      ELSE
        IDIG = ICHAR('0') + 1
        IND = 4 + IOFF
        MXFORM(IND:IND) = CHAR(IDIG)
        IND = IND + 1
        IDIG = ICHAR('0') + MOD(NDEC,10)
        MXFORM(IND:IND) = CHAR(IDIG)
        IOFF = IOFF + 1
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck opnfil */
      SUBROUTINE OPNFIL(IUNIT,NAME,STATUS,SUB)
C***********************************************************************
C
C	Open sequential, unformatted file.
C	Based on OPENDX
C	T.Saue, University of Tromsoe, Norway, Feb.-94
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) NAME,STATUS,SUB
      IF (STATUS .EQ. 'NEW') GO TO 10
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 100
C
C	Open old file
C
      OPEN(IUNIT,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',
     +     ACCESS='SEQUENTIAL',ERR=10)
C
C     Just in case....
      REWIND IUNIT
      GOTO 20
   10 CONTINUE
      IF(STATUS.EQ.'OLD') GOTO 110
C
C	Open new file
C
      OPEN(IUNIT,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     +     ACCESS='SEQUENTIAL')
   20 CONTINUE
      RETURN
C
C	Error messages
C
  100 CONTINUE
      WRITE(LUPRI,'(6A,I2,A)') SUB,': Invalid STATUS keyword ',STATUS,
     +     ' for file ',NAME,', (unit ',IUNIT,').'
      CALL QUIT('OPNFIL:Invalid STATUS keyword')
  110 CONTINUE
      WRITE(LUPRI,'(4A,I2,A)') SUB,': Old file ',NAME,' (unit',IUNIT,
     +    ') not found.'
      CALL QUIT('OPNFIL: Old file not found')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck putltr */
      SUBROUTINE PUTLTR(NM,N,ATRI,B)
C***********************************************************************
C
C	Copy elements from a rowwise lower triangularly packed matrix 
C       into a square matrix
C	
C	ATRI 	- rowwise lower triangularly packed matrix
C	B    	- square matrix
C	NM	- leading dimension of B
C	N	- order of matrix B
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION ATRI(N*(N+1)/2),B(NM,N)
C
      DO 10 I = 1,N
        IOFF = I*(I-1)/2
        DO 20 J = 1,I
          B(I,J) = ATRI(IOFF+J)
   20   CONTINUE
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck putsqr */
      SUBROUTINE PUTSQR(NM,NROW,NCOL,ASQR,B)
C***********************************************************************
C
C	Copy elements from columnwise array into rectangular matrix
C	
C	ASQR 	- rowwise lower triangularly packed matrix
C	B    	- rectangular matrix
C	NM	- leading dimension of B
C	NROW	- number of rows in B
C	NCOL    - number of columns in B
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION ASQR(NROW*NCOL),B(NM,NCOL)
C
      DO 10 J = 1,NCOL
        IOFF = (J-1)*NROW
        DO 20 I = 1,NROW
          B(I,J) = ASQR(IOFF+I)
   20   CONTINUE
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck unfld */
      SUBROUTINE UNFLD(NDIM,AR,AI)
C**********************************************************************
C	A Kramer symmetric matrix has the following structure:
C	
C  		-------------
C		|  A  |  B  |
C		-------------
C               | -B* |  A* |
C		-------------
C
C	in which A is a complex Hermitian matrix and B a complex 
C       antisymmetric matrix.
C
C	This subroutines assumes on input that only the lower triangle 
C       of blocks A and B are given, abd will on output give the full 
C       blocks A and B.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AR(NDIM,NDIM/2),AI(NDIM,NDIM/2)
      IF(MOD(NDIM,2).EQ.1) THEN
        WRITE(LUPRI,'(A,I10,A)') 'Error in UNFLD: NDIM = ',NDIM,
     +  ' is not an even number.'
        CALL QUIT('Uneven dimension in UNFLD')
      ENDIF
      NHALF = NDIM/2
      IOFF = NHALF + 1
C** Unfold real part:
      CALL UNFLDR('SYMM',NDIM,NHALF,AR)
      CALL UNFLDR('ASYM',NDIM,NHALF,AR(IOFF,1))
C** Unfold imaginary part:
      CALL UNFLDR('ASYM',NDIM,NHALF,AI)
      CALL UNFLDR('ASYM',NDIM,NHALF,AI(IOFF,1))
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck unfldr */
      SUBROUTINE UNFLDR(MATTYP,LDA,NDIM,A)
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
      CHARACTER*4 MATTYP
      DIMENSION A(LDA,NDIM)
C*SYMM  : Real symmetric matrix
C*ASYM  : Real antisymmetric matrix
C
C   ON INPUT: A - lower triangular matrix
C   ON OUTPUT:A - full matrix
      IF(MATTYP.EQ.'SYMM') THEN
      DO 10 J = 1,NDIM
        DO 20 I = (J+1),NDIM
          A(J,I) = A(I,J)
   20     CONTINUE
   10   CONTINUE
      ELSEIF(MATTYP.EQ.'ASYM') THEN
      DO 30 J = 1,NDIM
        DO 40 I = (J+1),NDIM
          A(J,I) = -A(I,J)
   40     CONTINUE
        A(J,J) = D0
   30   CONTINUE
      ELSE
        WRITE(LUPRI,'(A,A6,A)') '** ERROR in UNFLD **  : Keyword ',
     +    MATTYP,' not recognized!'
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sectid */
      FUNCTION SECTID(TID)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Given time in seconds this function return a character string
C     with time in suitable units
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
      CHARACTER*12 SECTID
C
      ITID = TID
      IF(ITID.GE.86400) GOTO 3
      IF(ITID.GE.3600)  GOTO 2
      IF(ITID.GE.60)    GOTO 1
C
C     Format : Seconds
C
      WRITE(SECTID,'(F11.8,A1)') TID,'s'
      RETURN
C
C     Format: Minutes and seconds
C
    1 CONTINUE
      IMIN = ITID/60
      TID  = TID - IMIN*60
      WRITE(SECTID,'(I2,A3,F6.3,A1)') IMIN,'min',TID,'s'
      RETURN
C
C     Format: Hours, minutes and seconds
C
    2 CONTINUE
      IHOUR = ITID/3600
      ITID  = MOD(ITID,3600)
      IMIN  = ITID/60
      ITID  = 3600*IHOUR + 60*IMIN
      ITID   = TID - ITID
      WRITE(SECTID,'(I3,A1,I2,A3,I2,A1)')
     &     IHOUR,'h',IMIN,'min',ITID,'s'
      RETURN
C
C     Format: Days, hours and minutes
C
    3 CONTINUE
      IDAY  = ITID/86400
      ITID  = MOD(ITID,86400)
      IHOUR = ITID/3600
      ITID  = MOD(ITID,3600)
      IMIN  = ITID/60
      WRITE(SECTID,'(I3,A1,I2,A1,I2,A3)')
     &      IDAY,'d',IHOUR,'h',IMIN,'min'
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck logset */
      SUBROUTINE LOGSET(N,LVAL,LVEC)
C***********************************************************************
C     This routine will initialize a logical array
C     Written by T.Saue , January 1995
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL LVAL,LVEC(*)
      DO 10 I = 1,N
        LVEC(I) = LVAL
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lbit */
C##############lbit##############################################
      LOGICAL FUNCTION LBIT(I,N)
C	Based on an analogous DISCO routine(Jan Almloef)
#include "implicit.h"
#include "ibtfun.h"
      LBIT = .FALSE.
      IF (IBTAND(IBTSHR(I,N-1),1).EQ.1) LBIT = .TRUE.
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lset */
      SUBROUTINE LSET(N,LVAL,LVEC)
#include "implicit.h"
      LOGICAL LVEC(N),LVAL
      DO 10 I = 1,N
        LVEC(I) = LVAL
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (SYS_CRAY) || defined (SYS_T3D)
      REAL FUNCTION GAMMP(A,X)
#else
      DOUBLE PRECISION FUNCTION GAMMP(A,X)
#endif
#include "implicit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C****************************************************************C
C     SUBFUNCTION TO EVALUATE THE NORMALIZED INCOMPLETE GAMMA    C
C     FUNCTION USING SERIES EXPANSION FOR SMALL X AND CONTINUED  C
C     FRACTION REPRESENTATION FOR 'LARGE' X (I.E.  X > A+1)      C
C                                                                C
C     FROM THE NUMERICAL RECIPES LIBRARY, BUT MODIFIED BY AN     C
C     ADDITIONAL FACTOR OF GAMMA(A) SO THAT GAMMAP RETURNS THE   C
C     UN-NORMALIZED GAMMA-FUNCTION (6.5.2), AND NOT THE          C
C     P-FUNCTION (6.5.1). EQUATION NUMBERS AND DEFINITIONS ARE   C
C     FROM ABRAMOWITZ AND STEGUN                                 C
C                                                                C
C****************************************************************C
      IF(X.LT.D0.OR.A.LE.D0) THEN
         GAMMP=D0
         CALL QUIT('GAMMP: Negative A or X..')
      ELSEIF(X.LT.A+D1) THEN
         CALL GSER(GAMMP,A,X,GLN)
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP = D1 - GAMMCF
      ENDIF
C
C     RETURN THE UNNORMALIZED INCOMPLETE GAMMA FUNCTION
C
      GAMMP = GAMMP * EXP(GLN)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GSER(GAMSER,A,X,GLN)
#include "implicit.h"
C******************************************************************C
C     SERIES EXPANSION REPRESENTATION                              C
C******************************************************************C
      PARAMETER (ITMAX=100,EPS=1.D-12,D1 = 1.0D0)
      GLN=GAMMLN(A)
      AP=A
      SUM=D1/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+D1
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      CALL QUIT('GSER: ITMAX Too Small')
C
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
#include "implicit.h"
C*******************************************************************C
C     CONTINUED FRACTION REPRESENTATION                             C
C*******************************************************************C
      PARAMETER (ITMAX=100,EPS=1.D-12,D0 = 0.0D0,D1 = 1.0D0)
#include "sdpre.h"
      GLN=GAMMLN(A)
      GOLD=D0
      A0=D1
      A1=X
      B0=D0
      B1=D1
      FAC=D1
      DO 11 N=1,ITMAX
        AN=SDPRE(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.D0)THEN
          FAC=D1/A1
          G=B1*FAC
          IF(DABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      CALL QUIT('GCF: ITMAX Too Small')
C
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION GAMMLN(XX)
#include "implicit.h"
C***********************************************************************C
C
C     GAMMALN RETURNS THE NATURAL LOGARITHM OF THE GAMMA FUNCTION       C
C     (FROM "NUMERICAL RECIPES")                                        C
C                                                                       C
C***********************************************************************C
      DIMENSION COF(6)
      DATA HALF,ONE,FPF/0.50D0,1.0D0,5.5D0/
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     1 -1.231739516D0,0.120858002D-2,-0.536382D-5,2.50662827465D0/
C
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
       X=X+ONE
       SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DFACTR
#include "implicit.h"
      PARAMETER(D1 = 1.0D0)
C***************************************************************C
C     DFACTR INITIALIZES A VECTOR OF DOUBLE-FACTORIAL FUNCTIONS C
C     SUBJECT TO THE CONVENTIONS                                C
C                                                               C
C     DFACT(I) = (I-2)!!                                        C
C     (-1)!! = 1                                                C
C     0!!    = 1                                                C
C     N!!    = N*(N-2)!!                                        C
C                                                               C
C***************************************************************C
      COMMON/FACTRL/DFACT(30)
C
#include "sdpre.h"
      DFACT(1)=D1
      DFACT(2)=D1
      DO 10 M=3,30
      DFACT(M)=SDPRE(M-2)*DFACT(M-2)
10    CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsjaco */
      SUBROUTINE RSJACO(LDM,NDIM,NROWE,AMAT,EIG,IJOB,IORDER,IPACK,EVEC,
     &                  WORK,LWORK)
C***********************************************************************
C
C     Diagonalize real symmetric matrix by the Jacobi method
C     This is basically a driver for the routine JACO.
C
C     Control parameter IORDER:
C       IORDER.EQ.-2 : Sort eigenvalues in absolute descending order
C       IORDER.EQ.-1 : Sort eigenvalues in descending order
C       IORDER.EQ. 1 : Sort eigenvalues in ascending order
C       IORDER.EQ. 2 : Sort eigenvalues in absolute ascending order
C       IORDER.EQ. 0 : No sorting of eigenvalues and eigenvectors
C     Control parameter IJOB:
C               IJOB = 0        Only eigenvalues desired.
C               IJOB = 1        Eigenvalues and eigenvectors in 
C                               orthonormal basis. EVEC set to identity 
C                               matrix.NROWE is dimension of basis 
C                               equal to NDIM.
C               IJOB = 2        Eigenvalues and eigenvectors in another 
C                               basis.EVEC equal to transformation to 
C                               other basis. NROWE is dimension of 
C                               basis and may differ from NDIM.
C     Control parameter IPACK:
C               IPACK = 0   : AMAT given as square matrix
C               IPACK = 1   : AMAT given as row-packed lower triangle
C
C     Written by T.Saue - Aug 15 1996
C     Last revision : Aug 15 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AMAT(*),EIG(NDIM),EVEC(LDM,NDIM),
     &          WORK(LWORK)
C
      CALL QENTER('RSDIAG')
#include "memint.h"
C
      CALL MEMGET('REAL',KBUF ,NDIM ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KIBUF,NDIM ,WORK,KFREE,LFREE)            
C
C     Determine job
C     =============
C
      IF    (IJOB.EQ.0) THEN
        NEVEC = 0
      ELSEIF(IJOB.EQ.1) THEN
        NEVEC = NDIM
        CALL DUNIT2(EVEC,NDIM,LDM,NDIM,1)
      ELSEIF(IJOB.EQ.2) THEN
        NEVEC = NDIM
      ELSE
        CALL QUIT('RSJACO: Illegal value of IJOB !')
      ENDIF
C
      IF(IPACK.EQ.0) THEN
C
C     AMAT is a full square matrix
C     ============================
C
C       Pack lower triangle of AMAT into ATRI
C
        NNDIM = (NDIM*(NDIM+1))/2
        CALL MEMGET('REAL',KATRI,NNDIM,WORK,KFREE,LFREE)
        CALL DAMATR(NDIM,AMAT,LDM,WORK(KATRI))
C
C       Diagonalize ATRI
C
        CALL JACO2(WORK(KATRI),EVEC,NDIM,NDIM,NROWE,LDM,
     &            WORK(KBUF),WORK(KIBUF))
C
C       Extract eigenvalues from ATRI
C
        CALL XTRCDI(WORK(KATRI),EIG,NDIM,1)
C
      ELSEIF(IPACK.EQ.1) THEN
C
C       AMAT is a row-packed lower triangular matrix
C       ============================================
C
C       Diagonalize AMAT
C
        CALL JACO(AMAT,EVEC,NDIM,NDIM,NEVEC,
     &            WORK(KBUF),WORK(KIBUF))
C
C       Extract eigenvalues from AMAT
C
        CALL XTRCDI(AMAT,EIG,NDIM,1)
      ELSE
        CALL QUIT('RSJACO: Illegal value of IPACK !')
      ENDIF
C
C     Order eigenvalues/vectors in ascending/descending order
C
      IF    (IORDER.NE.0) THEN
        CALL ORDER3(EVEC,EIG,LDM,NDIM,NEVEC,IORDER)
      ENDIF
C
      CALL MEMREL('RSJACO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RSDIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DAMATR(N,AMAT,LDM,ATRI)
C***********************************************************************
C
C     Row-pack lower triangle of matrix AMAT into ATRI
C
C     Written by T.Saue - Aug 16 1996
C     Last revision: Aug 16 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
      DIMENSION AMAT(LDM,N), ATRI(*)
C
#include "nosync.h"
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ATRI(JOFF+I) = AMAT(I,J)
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dunit2 */
      SUBROUTINE DUNIT2(A,N,LR,LC,NZ)
C***********************************************************************
C
C     Subroutine DUNIT2 sets a general square matrix A equal
C     to a unit matrix.
C     /VER 2/ 14-Sep-1983 hjaaj
C     Include leading dimension - T.Saue Aug 16 1996
C
C***********************************************************************
#include "implicit.h"
      DIMENSION A(LR,LC,NZ)
      PARAMETER (D1=1.0D00, D0=0.0D00)
C
      IF (N.GT.LR .OR. N.GT.LC) THEN
        CALL QUIT('Inconsisten N,LR,LC in DUNIT2 call')
      ENDIF
      DO IZ = 1,NZ
        DO J = 1,N
          DO I = 1,N
            A(I,J,IZ) = D0
          ENDDO
        ENDDO
      ENDDO
      DO J = 1,N
        A(J,J,1) = D1
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck jaco2 */
      SUBROUTINE JACO2(F,V,NB,NMAX,NROWV,LDM,BIG,JBIG)
C
C Revisions:
C
C  26-Aug-1984 tsaue  Introduced a new parameter LDM, which
C                     is the leading dimension of V and must
C                     be larger or equal to NROWV.
C   2-Nov-1984 hjaaj (new parameter NROWV such that
C                     dim(V) = (NROWV,NMAX). This makes
C                     it possible to solve eigenproblem
C                     in a reduced basis but get the
C                     eigenvectors in the original full
C                     basis, e.g. less mo's than ao's)
C  23-Feb-1989 hjaaj  Note that if NROWV = 0 then only
C                     eigenvalues will be calculated,
C                     V matrix will not be referenced.
C  27-Jul-1990 hjaaj  Changed -CX,+SX transformation to +CX,-SX
C                     transformation; probably the -CX,+SX
C                     transformation was responsible for that
C                     the eigenvectors easily changed sign.
C                     Changed initial test on NB. Changed SD.
C                     Optimized IR loop.
C     Jun-1992 ov     Parameters for 0.5, 1.5, ... (for Cray)
C  20-Jul-1992 hjaaj  Changed C1,C2 to THRZER
C  30-oct-1992 hjaaj  zero f(ab) to avoid round-off errors
C                     absolute conv.threshold SD=C1
C
#include "implicit.h"
      DIMENSION F(*),V(*)
      DIMENSION BIG(*) ,JBIG(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, ROOT2 = 0.707106781186548D0)
      PARAMETER(DP5 = 0.5D0, D1P5 = 1.5D0, D1P375 = 1.375D0,
     *          D3P875 = 3.875D0, DP25 = 0.25D0)
#include "thrzer.h"
C      DATA C1,C2,C3,C4,C5,C6/THRZER,THRZER,1.D-20,1.D-14,1.D-9,1.D-5/
COLD   DATA C1,C2,C3,C4,C5,C6/1.D-12,1.D-12,1.D-20,1.D-14,1.D-9,1.D-5/
      DATA C1,C2,C3,C4,C5,C6/1.D-15,1.D-15,1.D-20,1.D-14,1.D-9,1.D-5/
      IF (NB.LE.1 .OR. NMAX.LE.0) RETURN
Cwas: IF (NB.EQ.1) RETURN !900727-hjaaj
      IF(LDM.LT.NROWV) CALL QUIT('JACO2: LDM.LT.NROWV')
      DO 190 I=1,NB
         JBIGI=0
         J=MIN(I-1,NMAX)
         IF (J .GT. 0) THEN
            II = (I*I-I)/2
            ABIGI=D0
            DO 18 K=1,J
            IF (ABIGI .GE. ABS(F(II+K))) GO TO  18
               ABIGI=ABS(F(II+K))
               JBIGI=K
   18       CONTINUE
         END IF
         IF (JBIGI .GT. 0) THEN
            JBIG(I) = JBIGI
            BIG(I)  = F(II+JBIGI)
         ELSE
            JBIG(I) = 0
            BIG(I)  = D0
         END IF
  190 CONTINUE
C
#if defined (VAR_OLDCODE)
C 900727-hjaaj:
C SD calculation was done in every Jacobi iteration.
C Now the largest absolute element in F is found once and
C the SD based on that value is used in every iteration.
  410 SD=1.05D 00
      DO 220 J=1,NMAX
         DAB=ABS(F(J*(J+1)/2))
CHJ-861103: commented out next line, it seems to make the loop
C           meaningless (setting SD equal to J=NMAX value always!)
C        IF (SD .GT. DAB) SD=DAB
  220    SD=MAX(SD,DAB)
      SD=MAX(C1,C2*SD)
#else
C 921030-hjaaj: SD = C1 now
      NNB = (NB*NB+NB)/2
C     SD = 1.05D0
C     DO 220 J = 1,NNB
C        SD = MAX(SD, ABS(F(J)) )
C 220 CONTINUE
C     SD=MAX(C1,C2*SD)
      SD=C1
C
      MXITJA = 50*NNB
      ITJACO = 0
  410 ITJACO = ITJACO + 1
      IF (ITJACO .GT. MXITJA) THEN
         CALL QUIT('ERROR: JACO did not converge ...')
      END IF
#endif
      T = D0
      DO 230 I=2,NB
      IF (T .GE.  ABS(BIG(I))) GO TO 230
         T = ABS(BIG(I))
         IB= I
  230 CONTINUE
      IF(T.LT.SD) GO TO 420
         IA =JBIG(IB)
         IAA=IA*(IA-1)/2
         IBB=IB*(IB-1)/2
         DIF=F(IAA+IA)-F(IBB+IB)
         IF( ABS(DIF) .GT. C3) GO TO 271
            SX=ROOT2
            CX=ROOT2
         GO TO 270
  271       T2X2 =BIG(IB)/DIF
            T2X25=T2X2*T2X2
         IF(T2X25 .GT. C4) GO TO 240
            CX=1.D 00
            SX=T2X2
         GO TO 270
  240    IF(T2X25 .GT. C5) GO TO 250
            SX=T2X2*(D1 - D1P5*T2X25)
            CX=D1 - DP5*T2X25
         GO TO 270
  250    IF(T2X25 . GT . C6) GO TO 260
            CX=D1+T2X25*(T2X25*D1P375 - DP5 )
            SX= T2X2*(D1 + T2X25*(T2X25*D3P875 - D1P5))
         GO TO 270
  260       T = DP25  / SQRT(DP25   + T2X25)
            CX= SQRT(DP5   + T)
            SX= SIGN( SQRT(DP5 - T),T2X2)
  270    CONTINUE
#include "ivdep.h"
         DO 275 IR=1,IA
            T        = F(IAA+IR)*SX
            F(IAA+IR)= F(IAA+IR)*CX+F(IBB+IR)*SX
            F(IBB+IR)=-T           +F(IBB+IR)*CX
  275    CONTINUE
         IEAA=IAA+IA
         IEAB=IBB+IA
         TT  =F(IEAB)
         F(IEAB)=BIG(IB)
         IF (JBIG(IA) .EQ. 0) THEN
            IRST = IA   + 1
            IEAR = IEAA + IA
            IEBR = IEAB + 1
         ELSE
            IRST = IA
            IEAR = IEAA
            IEBR = IEAB
         END IF
         DO 390 IR = IRST,NB
#if !defined (VAR_OLDCODE)
            IF (IR .EQ. IA) GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
#else
            IF (IR .EQ. IA) THEN
               GO TO 360
C              ... we have checked above that JBIG(IA) .ne. 0
C              IF(JBIG(IR)) 360,380,360
            END IF
#endif
            T      = F(IEAR)*SX
            F(IEAR)= F(IEAR)*CX+F(IEBR)*SX
            F(IEBR)=-T         +F(IEBR)*CX
            T   =F(IEAR)
            IT  =IA
            IF(IR-IB) 340,310,320
  310          F(IEAA)=F(IEAA)*CX+F(IEAB)*SX
C              921030+hjaaj: zero f(ab) to avoid round-off errors
C              F(IEAB)=     TT*CX+F(IEBR)*SX
               F(IEAB)=     D0
               F(IEBR)=    -TT*SX+F(IEBR)*CX
            GO TO 360
  320       IF(ABS(T) .GE.  ABS(F(IEBR))) GO TO 340
               T   =F(IEBR)
               IT  =IB
  340       IF(ABS(T) .LT.  ABS(BIG(IR))) GO TO 350
               BIG(IR)  = T
               JBIG(IR) = IT
            GO TO 380
  350       IF(IA .NE. JBIG(IR) .AND. IB .NE. JBIG(IR))  GO TO 380
  360          K= IEAR - IA
               JBIGI = 0
               IR1=MIN (IR-1,NMAX)
               IF (IR1 .GT. 0) THEN
                  ABIGI = D0
                  DO 370 I=1,IR1
                  IF(ABIGI .GE. ABS(F(K+I)))  GO TO 370
                     ABIGI = ABS(F(K+I))
                     JBIGI =I
  370             CONTINUE
               END IF
               IF (JBIGI .GT. 0) THEN
                  JBIG(IR) = JBIGI
                  BIG(IR)  = F(K+JBIGI)
               ELSE
                  JBIG(IR) = 0
                  BIG(IR)  = D0
               END IF
  380          CONTINUE
               IEAR = IEAR + IR
               IF (IR .GE. IB) THEN
                  IEBR = IEBR + IR
               ELSE
                  IEBR = IEBR + 1
               END IF
  390       CONTINUE
Cwas         JAA=(IA-1)*NROWV
Cwas         JBB=(IB-1)*NROWV
         JAA=(IA-1)*LDM
         JBB=(IB-1)*LDM
#include "ivdep.h"
         DO 400 I=1,NROWV
            T=V(JBB+I)*SX
            V(JBB+I)=-V(JAA+I)*SX + V(JBB+I)*CX
  400       V(JAA+I)= V(JAA+I)*CX + T
      GO TO 410
  420 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck order3 */
      SUBROUTINE ORDER3(EVEC,EVAL,LDM,N,NEVEC,IJOB)
C***********************************************************************
C
C Purpose: order the N values in EVAL and their associated vectors 
C     in EVEC:
C     IJOB.EQ. 2  -  absolute ascending order  
C                    ABS(EVAL(i+1)) .ge. ABS(EVAL(i))
C     IJOB.EQ. 1  -  ascending order   EVAL(i+1) .ge. EVAL(i)
C     IJOB.EQ.-1  -  descending order   EVAL(i+1) .le. EVAL(i)
C     IJOB.EQ.-2  -  absolute descending order  
C                    ABS(EVAL(i+1)).le. ABS(EVAL(i))
C     IJOB.EQ.0  -  return
C          
C          Generalization of routines ORDER/ORDER2
C          Written by T.Saue - Aug 26 1996
C          Last revision: Aug 26 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
      DIMENSION EVEC(*),EVAL(*)
      IF (N.LE.1.OR.IJOB.EQ.0) RETURN
      II = IJOB + 3
      GOTO (1,2,3,4,5) II
1     CONTINUE
C***********************************************************************
C IJOB.EQ.-2  -  absolute descending order 
C                ABS(EVAL(i+1)).le. ABS(EVAL(i))
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (ABS(EVAL(J)) .GE. ABS(EMIN)) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)   =EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
2     CONTINUE
C***********************************************************************
C      IJOB.EQ.-1  -  descending order   EVAL(i+1) .le. EVAL(i)
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (EVAL(J) .GE. EMIN) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
3     CONTINUE
C***********************************************************************
C   IJOB.EQ. 0  -  return
C***********************************************************************
      RETURN
4     CONTINUE
C***********************************************************************
C   IJOB.EQ. 1  -  ascending order   
C                  EVAL(i+1) .ge. EVAL(i)
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
5     CONTINUE
C***********************************************************************
C   IJOB.EQ. 2  -  absolute ascending order  
C                  ABS(EVAL(i+1)) .ge. ABS(EVAL(i))
C***********************************************************************
      IN = 1
      DO I=1,N-1
         EMIN = EVAL(I)
         IMIN = I
         DO J=I+1,N
            IF (ABS(EVAL(J)).LT.ABS(EMIN)) THEN
               EMIN = EVAL(J)
               IMIN = J
            ENDIF
         ENDDO
         IF (IMIN.NE.I) THEN
            EVAL(IMIN)=EVAL(I)
            EVAL(I)=EMIN
            IF (NEVEC .GT. 0) THEN
              CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*LDM+1),1)
            ENDIF
         ENDIF
         IN = IN + LDM
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatpak */
      SUBROUTINE GATPAK(NA,A,NB,B,IND,IOFF)
C***********************************************************************
C
C     This routine will gather elements of the same "kin" in the
C     rowpacked lower triangular matrix A and place them in the
C     rowpacked lower triangular matrix B
C
C     "Kinship" is determined by the integer array IND pointing from
C     elements of A to elements of B
C
C     Written by T.Saue  -  July 23 1994
C     LAST VERSION: July 23 1994
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(*),B(*),IND(*)
C
      KA = 0
      DO 10 IA = 1,NA
        IB = IND(IA) - IOFF
        DO 20 JA = 1,IA
          JB = IND(JA) - IOFF
          KA = KA + 1
          KB = (IB*(IB-1)/2)+JB
          B(KB) = B(KB) + A(KA)
   20   CONTINUE
   10 CONTINUE
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igath */
      SUBROUTINE IGATH(IV1,N1,IV2,N2,N,IND)
C***********************************************************************
C
C     Reorder an integer array
C
C     Written by T..Saue November 1994
C     Last revision: Nov 16 1994 - tsaue
C***********************************************************************
#include "implicit.h"
      DIMENSION IV1(N1),IV2(N2),IND(*)
      DO 10 I = 1,N
        IV2(I) = IV1(IND(I))
   10 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck privec */
      SUBROUTINE PRIVEC(NAME,IVEC,NVEC)
C***********************************************************************
C
C     Print integer array with indices in compact form
C
C     Written by T.Saue, November 1994
C     Last revised: Nov 24 1994 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER NAME*6
      DIMENSION IVEC(NVEC)
      NROW = (NVEC+15)/16
      J1 = 1
      DO 10 J = 1,NROW
        J2 = MIN(J1+15,NVEC)
        WRITE(LUPRI,'(3X,A3,4X,20I6)') 'IND',(JJ, JJ = J1,J2)
        WRITE(LUPRI,'(3X,A6,1X,20I6)')
     &            NAME,(IVEC(JJ),JJ = J1,J2)
        J1 = J1 + 16
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mtcomp */
      SUBROUTINE MTCOMP(COMP,AMAT,NDIM,NZ,CMAT)
C***********************************************************************
C
C     This subroutine will pick out symmetric('S') or antisymmetric('A')
C     component of a quadratic matrix AMAT and place it in CMAT.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DP5 = 0.5D0)
      DIMENSION AMAT(NDIM,NDIM,NZ),CMAT(NDIM,NDIM,NZ)
      CHARACTER COMP*1
      N2DIMQ = NDIM*NDIM*NZ
C
C     Symmetric component
C
      IF    (COMP.EQ.'S') THEN
        DO 10 IZ = 1,NZ
          DO 20 J = 1,NDIM
            DO 30 I = 1,NDIM
              CMAT(I,J,IZ) = AMAT(I,J,IZ)+AMAT(J,I,IZ)
 30         CONTINUE
 20       CONTINUE
 10     CONTINUE
        CALL DSCAL(N2DIMQ,DP5,CMAT,1)
      ELSEIF(COMP.EQ.'A') THEN
        DO 40 IZ = 1,NZ
          DO 50 J = 1,NDIM
            DO 60 I = 1,NDIM
              CMAT(I,J,IZ) = AMAT(I,J,IZ)-AMAT(J,I,IZ)
 60         CONTINUE
 50       CONTINUE
 40     CONTINUE
        CALL DSCAL(N2DIMQ,DP5,CMAT,1)
      ELSE
        CALL QUIT('MTCOMP:Unknown component '//COMP)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck difmmh */
      SUBROUTINE DIFMMH(AMAT,NDIM,NZ,LRA,LCA)
C***********************************************************************
C
C     Take difference between a general quaternionic matrix AMAT and
C     its Hermitian conjugate; return in AMAT. Note that this gives
C     an anti-Hermitian result !
C
C     Written by T.Saue, August 23 1995
C     Last revision: tsaue - aug23-95
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION AMAT(LRA,LCA,NZ)
C
C     Take difference
C
      DO J = 1,NDIM
        DO I = 1,(J-1)
          AMAT(I,J,1) = AMAT(I,J,1)-AMAT(J,I,1)
          AMAT(J,I,1) = -AMAT(I,J,1)
        ENDDO
        AMAT(J,J,1) = D0
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(I,J,IZ) = AMAT(I,J,IZ)+AMAT(J,I,IZ)
            AMAT(J,I,IZ) = AMAT(I,J,IZ)
          ENDDO
          AMAT(J,J,IZ) = AMAT(J,J,IZ)+AMAT(J,J,IZ)
        ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE QGETHE(AMAT,NDIM,NZ,LRA,LCA)
C***********************************************************************
C
C     Quaternian GEneral matrix To HErmitian matrix:
C     generate the Hermitian component of the general matrix AMAT;
C     return in AMAT.
C
C     Written by H.J.Aa. Jensen, June 2002
C     (based on DIFMMH)
C
C***********************************************************************
#include <implicit.h>
      PARAMETER(D0=0.0D0, DP5=0.5D0)
      DIMENSION AMAT(LRA,LCA,NZ)
C
C     Generate A + A(dagger)
C
      DO J = 1,NDIM
        DO I = 1,(J-1)
          AMAT(I,J,1) = DP5*(AMAT(I,J,1)+AMAT(J,I,1))
          AMAT(J,I,1) = AMAT(I,J,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(I,J,IZ) = DP5*(AMAT(I,J,IZ)-AMAT(J,I,IZ))
            AMAT(J,I,IZ) = -AMAT(I,J,IZ)
          ENDDO
          AMAT(J,J,IZ) = D0
        ENDDO
      ENDDO
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtdac */
      SUBROUTINE WRTDAC (IT,N,X,IREC)
#include "implicit.h"
      DIMENSION X(N)
      WRITE (IT, REC=IREC) X
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readac */
      SUBROUTINE READAC (IT,N,X,IREC)
#include "implicit.h"
      DIMENSION X(N)
      READ (IT, REC=IREC, IOSTAT=IOS, ERR=10) X
      RETURN
 10   CONTINUE
      CALL QUIT('READAC: Error reading file')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prmutc */
      SUBROUTINE PRMUTC(AMAT,NCOL,NCTL,LUPRI)
C***********************************************************************
C
C  Based on OUTPAK by Hans Jorgen Aa. Jensen.
C
C PRMUTC prints a real symmetric matrix stored in column-packed upper
C triangular form (see diagram below) in formatted form with numbered
C rows and columns.  The input is as follows:
C
C        AMAT(*).............packed matrix
C        NCOL................number of columns to be output
C        NCTL................carriage control flag: 1 for single space,
C                                                   2 for double space,
C                                                   3 for triple space.
C
C The matrix elements are arranged in storage as follows:
C
C      1   2   4   7  11
C          3   5   8  12
C              6   9  13
C                 10  14
C                     15
C       and so on.
C
C PRMUTC is set up to handle 6 columns/page with a 6F20.14 format
C for the columns.  If a different number of columns is required, change
C formats 1000 and 2000, and initialize kcol with the new number of
C columns.
#include "implicit.h"
      DIMENSION AMAT(*)
      INTEGER BEGIN
      CHARACTER*1 ASA(3),BLANK,CTL,FX,FPRI
      CHARACTER   PFMT*27, COLUMN*8, LFMT*13,XFMT*6,RFMT*6
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
C     First check whether the matrix is all zero
C
      J = NCOL*(NCOL+1)/2
      AMAX = ZERO
      AMIN = ZERO
      DO 5 I=1,J
         AMAX = MAX(AMAX,ABS(AMAT(I)))
         AMIN = MIN(AMIN,ABS(AMAT(I)))
    5 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 200
      END IF
C
C     Determine output format
C
      IF (AMIN.GT.FFMIN.AND.AMAX.LE.FFMAX) THEN
C        use F output format
         LFMT = '(   A1,I7,2X,'
         XFMT = '(15X),'
         RFMT = 'F15.8)'
      ELSE
C        use 1PD output format
         LFMT = '(1P,A1,I7,2X,'
         XFMT = '(15X),'
         RFMT = 'D15.8)'
      END IF
C
C LAST is the last column number in the row currently being printed
C
      LAST = MIN(NCOL,KCOL)
C
C BEGIN is the first column number in the row currently being printed.
C
C.....BEGIN NON STANDARD DO LOOP.
      BEGIN= 1
 1050 CONTINUE
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 40 IROW = 1,LAST
            ICOL = MAX(IROW,BEGIN)
            IOFF = (ICOL*(ICOL-1))/2 + IROW
            NPRI = LAST - ICOL + 1
            DO 10 J = 0,(NPRI-1)
               IF (AMAT(IOFF+ICOL*J+(J*(J-1))/2).NE.ZERO) GO TO 20
   10       CONTINUE
            GO TO 40
   20       CONTINUE
            IDIG = ICHAR('0') + NPRI
            FPRI = CHAR(IDIG)
            NX   = ICOL - BEGIN
            IDIG = ICHAR('0') + NX
            FX   = CHAR(IDIG)
            IF(NX.GT.0) THEN
              PFMT = LFMT//FX//XFMT//FPRI//RFMT
            ELSE
              PFMT = LFMT//'      '//FPRI//RFMT
            ENDIF
            WRITE (LUPRI,PFMT) CTL,IROW,
     &           (AMAT(ICOL*J+(J*(J-1))/2+IOFF),J=0,(NPRI-1))
   40    CONTINUE
         LAST = MIN(LAST+KCOL,NCOL)
         BEGIN= BEGIN + KCOL
      IF (BEGIN.LE.NCOL) GO TO 1050
  200 CONTINUE
      RETURN
C
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtinfo */
      SUBROUTINE GTINFO(TEXT)
C
C Written by T.Saue Sept 1 1995
C Based on TSTAMP by Hans Joergen Aa. Jensen 9-Jul-1990
C
#include "implicit.h"
      CHARACTER TEXT*(*)
C
#if defined (SYS_IBM)
      CHARACTER DAOFWK*9
#endif
#if defined (SYS_UNIX) || defined (SYS_AIX) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_HPUX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_PARAGON) || defined (SYS_SX) || defined (SYS_LINUX)
      CHARACTER*(24) FDATE
#endif
C
      TEXT = '                        '
#if defined (SYS_VAX)
      CALL DATE(TEXT(1))
      CALL TIME(TEXT(12))
#endif
#if defined (SYS_IBM)
      CALL DATE(DAOFWK,TEXT(1))
      CALL TIME(TEXT(12))
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
      WRITE(TEXT(1:8),'(A8)') DATE()
      WRITE(TEXT(12:19),'(A8)') CLOCK()
#endif
#if defined (SYS_AIX) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_HPUX) || defined (SYS_IRIX) || defined (SYS_PARAGON) || defined (SYS_SX) || defined (SYS_LINUX)
      TEXT = FDATE()
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prsymb */
      SUBROUTINE PRSYMB(IUNIT,SYMB,N,IX)
C***********************************************************************
C
C     Print a symbol/character N times on the same line
C
C     Written by T.Saue Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER SYMB*1,FMT*9
      FMT = '         '
      IF(IX.EQ.0) THEN
        WRITE(FMT(1:5),'(I3,A1)') N,'A1'
      ELSE
        WRITE(FMT,'(I3,A2,I3,A1)') IX,'X,',N,'A1'
      ENDIF
      WRITE(IUNIT,'('//FMT//')') (SYMB,I=1,N)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readns */
      SUBROUTINE REDSYM(IRREP,IGROUP,KSOP,NELM)
C***********************************************************************
C
C     Determine point group of reduced symmetry for a non-symmetric
C     operator
C
C     Written by T.Saue - May 15 1996
C     Last revision: May 15 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
C
      DIMENSION KSOP(0:7),KBUF(0:7)
      DO I = 0,7
        KBUF(I) = 0
      ENDDO
      NELM = -1
C
C     Go through characters for the irrep to find the
C     totally symmetric operations
C     KSOP is a pointer to those....
C
      DO I = 0,MAXREP
      II = JSOP(I)
      IF(IXVAL(IRREP,II).EQ.1) THEN
        NELM = NELM + 1
        KSOP(NELM) = II
        KBUF(II)  = 1
      ENDIF
      ENDDO
      NELM = NELM+1
C
C     Sum up number of rotations, reflections and inversion to
C     determine the reduced group
C
      KROTS = KBUF(3)+KBUF(5)+KBUF(6)
      KREFL = KBUF(4)+KBUF(2)+KBUF(1)
      KINVC = KBUF(7)
      IGROUP = MIN(7,NINT((4*KROTS+8*KINVC+6*KREFL)/3.0))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatmat */
      SUBROUTINE GATMAT(IOPT,LDA,N,A,B,IND,NIND)
C***********************************************************************
C
C     This routine will collect columns elements of square matrix A
C     in groups based on pointer index IND.
C     IOPT discerns two options:
C	IOPT = 0	pick out the elements of largest absolute value
C	IOPT = 1        sum up elements
C
C	Written by T.Saue, June 1994
C	Last revision: October 9 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(LDA,N),B(NIND,NIND),IND(N)
C
      CALL QENTER('GATMAT')
C
CTROND No initialization !!!!      CALL DZERO(B,NIND*NIND)
C
C     Pick out largest elements for each group
C     ========================================
C
      IF    (IOPT.EQ.0) THEN
        DO 10 J = 1,N
          JJ = IND(J)
          DO 20 I = 1,N
            II = IND(I)
            IF(ABS(A(I,J)).GT.ABS(B(II,JJ))) B(II,JJ) = ABS(A(I,J))
   20     CONTINUE
   10   CONTINUE
C
C     Sum up elements within each group
C     =================================
C
      ELSEIF(IOPT.EQ.1) THEN
        DO 30 J = 1,N
          JJ = IND(J)
          DO 40 I = 1,N
            II = IND(I)
            B(II,JJ) = B(II,JJ) + A(I,J)
   40     CONTINUE
   30   CONTINUE
      ELSE
        CALL QUIT('GATMAT:Unknown option IOPT!')
      ENDIF
C
      CALL QEXIT('GATMAT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gatvcc */
      SUBROUTINE GATVCC(IOPT,NVEC,A,N,B,IND,NIND)
C***********************************************************************
C
C     This routine will collect columns elements of vector A
C     in groups based on pointer index IND.
C     IOPT discerns two options:
C	IOPT = 0	pick out the elements of largest absolute value
C	IOPT = 1        sum up elements
C
C	Written by T.Saue, June 1994
C	Last revision: October 9 1995 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(N,NVEC),B(NIND,NVEC),IND(N)
CTROND No initialization !!!!      CALL DZERO(B,NIND*NIND)
C
C     Pick out largest elements for each group
C     ========================================
C
      IF    (IOPT.EQ.0) THEN
        DO J = 1,NVEC
          DO I = 1,N
            II = IND(I)
            B(II,J) = MAX(B(II,J),ABS(A(I,J)))
          ENDDO
        ENDDO
C
C     Sum up elements within each group
C     =================================
C
      ELSEIF(IOPT.EQ.1) THEN
        DO J = 1,N
          DO I = 1,N
            II = IND(I)
            B(II,J) = B(II,J) + A(I,J)
          ENDDO
        ENDDO
      ELSE
        CALL QUIT('GATVCC:Unknown option IOPT!')
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qaxpy */
      SUBROUTINE QAXPY(N,FA,TA,A,IQPA,NZA,
     &                 X,NX,INCX,IQPX,NZX,Y,NY,INCY,IQPY,NZY)
C***********************************************************************
C     Computes a constant times a vector plus av vector
C
C        Y := A*X + Y
C
C     The vectors and constant may be real(NZ=1),complex(NZ=2),
C     or quaternionic(NZ = 4).
C     NX and NY are the total lengths of vectors X and Y,
C     whereas INCX and INCY are increments.
C
C     op(A) is determined by the corresponding character variable
C     FA:
C       FA = 'N'        Normal 
C       FA = 'C'        Complex conjugate
C     and TX
C       TA = 'N'        No transform    a+bj
C       TA = 'I'        i-transform   -i(a+bj)i = a - bj
C       TA = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TA = 'K'        k-transform   -k(a+bj)k = a*-b*j
C     Calls BLAS routine DAXPY
C
C     Written by T.Saue, University of Oslo, Apr 7 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0=0.0D0)
C
C     Global variables
C
      DIMENSION A(NZA),X(NX,NZX),Y(NY,NZY),
     &          IQPA(NZA),IQPX(NZX),IQPY(NZY)
      CHARACTER FA*1,TA*1
C
C     Local variables
C
      DIMENSION IPA(4)
C
      EXTERNAL BDQALGEBRA
      IF (N.LE.0)       RETURN
C
      IF(N*INCX.GT.NX) CALL QUIT('QAXPY:N*INCX.GT.NX')
      IF(N*INCY.GT.NY) CALL QUIT('QAXPY:N*INCY.GT.NY')
C
C     Make pointer IPA
C     ================
C
      CALL IZERO(IPA,4)
      DO IZA = 1,NZA
        IPA(IQPA(IZA)) = IZA
      ENDDO
C
C     Constant A
C

C     Determine FX:
C     =============
C     Normal number...
      IF    (FA.EQ.'N') THEN
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        IFA = 2
      ELSE
        CALL QUIT('QAXPY:Unknown FA '//FA//' of scalar A!')
      ENDIF
C
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QAXPY:Unknown TA '//TA//' of scalar A!')
      ENDIF
C
      DO 10 IZY = 1,NZY
        IQY = IQPY(IZY)
        DO 20 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQA = IQMULT(IQX,IQY,1)
          IZA = IPA(IQA)
          IF(IZA.EQ.0) GOTO 20
          IF(A(IZA).EQ.D0) GOTO 20
          AA = A(IZA)*(IQSIGN(IQA,IFA,ITA)*IQPHASE(IQA,IQX,1))
          CALL DAXPY(N,AA,X(1,IZX),INCX,Y(1,IZY),INCY)
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qxapy */
      SUBROUTINE QXAPY(N,FA,TA,A,IQPA,NZA,
     &                         X,NX,INCX,IQPX,NZX,
     &                         Y,NY,INCY,IQPY,NZY)
C***********************************************************************
C     Computes a vector times a constant plus av vector
C
C        Y := X*A + Y
C
C     The vectors and constant may be real(NZ=1),complex(NZ=2),
C     or quaternionic(NZ = 4).
C     NX and NY are the total lengths of vectors X and Y,
C     whereas INCX and INCY are increments.
C
C     op(A) is determined by the corresponding character variable
C     FA:
C       FA = 'N'        Normal 
C       FA = 'C'        Complex conjugate
C     and TX
C       TA = 'N'        No transform    a+bj
C       TA = 'I'        i-transform   -i(a+bj)i = a - bj
C       TA = 'J'        j-transform   -j(a+bj)j = a*+b*j
C       TA = 'K'        k-transform   -k(a+bj)k = a*-b*j
C     Calls BLAS routine DAXPY
C
C     Written by T.Saue, University of Oslo, Apr 7 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0=0.0D0)
C
C     Global variables
C
      DIMENSION A(NZA),X(NX,NZX),Y(NY,NZY),
     &          IQPA(NZA),IQPX(NZX),IQPY(NZY)
      CHARACTER FA*1,TA*1
C
C     Local variables
C
      DIMENSION IPA(4)
C
      EXTERNAL BDQALGEBRA
      IF (N.LE.0)       RETURN
C
      IF(N*INCX.GT.NX) CALL QUIT('QXAPY:N*INCX.GT.NX')
      IF(N*INCY.GT.NY) CALL QUIT('QXAPY:N*INCY.GT.NY')
C
C     Make pointer IPA
C     ================
C
      CALL IZERO(IPA,4)
      DO IZA = 1,NZA
        IPA(IQPA(IZA)) = IZA
      ENDDO
C
C     Constant A
C

C     Determine FX:
C     =============
C     Normal number...
      IF    (FA.EQ.'N') THEN
        IFA = 1
C     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        IFA = 2
      ELSE
        CALL QUIT('QXAPY:Unknown FA '//FA//' of scalar A!')
      ENDIF
C
C     Determine TX:
C     =============
C     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
C     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
C     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
C     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        CALL QUIT('QXAPY:Unknown TA '//TA//' of scalar A!')
      ENDIF
C
      DO 10 IZY = 1,NZY
        IQY = IQPY(IZY)
        DO 20 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQA = IQMULT(IQX,IQY,1)
          IZA = IPA(IQA)
          IF(IZA.EQ.0) GOTO 20
          IF(A(IZA).EQ.D0) GOTO 20
          AA = A(IZA)*(IQSIGN(IQA,IFA,ITA)*IQPHASE(IQX,IQA,1))
          CALL DAXPY(N,AA,X(1,IZX),INCX,Y(1,IZY),INCY)
   20   CONTINUE
   10 CONTINUE
      RETURN
C
      END
C  /* Deck fndamx */
      SUBROUTINE FNDAMX(VEC,NDIM,MXELMN,IPLACE,NELMN,IPRT,THRES)
C
C 890205 hjaaj. Based on FNDMN3 by Jeppe Olsen.
C Last revision 900405 hjaaj.
C
C FIND MXELMN/NELMN ELEMENTS IN VEC with highest absolute value.
C NELMN IS THE LARGEST NUMBER LOWER THAN MXELMN THAT DOES NOT
C SPLIT DEGENERATE PAIRS.
C ORIGINAL PLACES OF THE LOWEST ELEMENTS ARE STORED IN IPLACE.
C
#include "implicit.h"
      DIMENSION VEC(NDIM),IPLACE(*)
C
      PARAMETER ( D1 = 1.0D0 )
C
C. FIRST OCCURANCE OF LOWEST ELEMENT AND LARGEST ELEMENT
      XMIN = ABS(VEC(1))
      XMAX = ABS(VEC(1))
      IMIN = 1
      IMAX = 1
      DO 100 I = 2,NDIM
        AVECI = ABS(VEC(I))
        IF( AVECI .GT. XMAX) THEN
           XMAX = AVECI
           IMAX = I
        ELSE IF( AVECI .LT. XMIN) THEN
           XMIN = AVECI
           IMIN = I
        END IF
  100 CONTINUE
C
      IF (IPRT.GT.0) WRITE(6,*) ' >> output from FNDAMX <<'
      IF(IPRT .GE. 5 ) THEN
         WRITE(6,*) ' LOWEST VALUE AND PLACE ',XMIN,IMIN
         WRITE(6,*) ' HIGHST VALUE AND PLACE ',XMAX,IMAX
      END IF
C
      IPLACE(1) = IMAX
      NDEG = 1
      IF(IPRT .GE. 5 ) WRITE(6,'(A,I8,1P,D15.5,I8,I4)')
     &    'IELMNT value IMAX NDEG ', 1,VEC(IMAX),IMAX,NDEG
C
      ITOP = MIN(NDIM,MXELMN+1)
      DO 200 IELMNT = 2,ITOP
        XMAXPR = XMAX
        IMAXPR = IMAX
        XMAX   = -D1
        IMAX   = -1
        DO 150 I = 1,NDIM
          AVECI = ABS(VEC(I))
          IF (AVECI .GT. XMAX) THEN
            IF (AVECI .LT. XMAXPR) THEN
               IMAX = I
               XMAX = AVECI
            ELSE IF (AVECI .EQ. XMAXPR .AND. I .GT. IMAXPR) THEN
               IMAX = I
               XMAX = AVECI
               GO TO 151
            END IF
          END IF
  150   CONTINUE
  151   CONTINUE
        IF(XMAXPR-XMAX .LE. THRES )THEN
          NDEG = NDEG + 1
        ELSE
          NDEG = 1
        END IF
C
        IF( IELMNT .LE. MXELMN ) IPLACE(IELMNT) = IMAX
        IF( IPRT .GE. 5 ) WRITE(6,'(A,I8,1P,D15.5,I8,I4)')
     &    ' IELMNT value IMAX NDEG ',IELMNT,VEC(IMAX),IMAX,NDEG
  200 CONTINUE
C
C CHECK DEGENERACY ON LAST VALUE
      IF(MXELMN .LT. NDIM ) THEN
        NELMN = MXELMN+1-NDEG
      ELSE
        NELMN = NDIM
      END IF
      IF (IPRT .GT. 0)
     &   WRITE(6,*) ' NUMBER OF ELEMENTS OBTAINED IN FNDAMX ',NELMN
C
C
      IF( IPRT  .GE. 3 ) THEN
        WRITE(6,*) ' FROM FNDAMX : '
        WRITE(6,*) '   PLACES OF ELEMENTS of highest absolut value'
        CALL IWRTMA(IPLACE,1,NELMN ,1,NELMN )
      END IF
C
      IF( IPRT .GE. 1 ) THEN
       WRITE(6,*)
     & ' MAX AND MIN IN SELECTED SUPSPACE ',
     &   VEC(IPLACE(1)),VEC(IPLACE(NELMN))
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Extelm */
      SUBROUTINE EXTELM(FA,A,LDA,NA,B,LDB,NBR,NBC,IBR,IBC)
C***********************************************************************
C
C     Extract elements of matrix A(LDA,NA) and
C     insert into matrix B(LDB,NBC) using 
C     row array IBR and column array IBC.
C
C     Written by T.Saue Mar 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER FA*1
      DIMENSION A(LDA,NA),B(LDB,NBC),IBR(NBR),IBC(NBC)
C
C     Extract from normal matrix
C
      IF(FA.EQ.'N') THEN
        DO J = 1,NBC
          JA = IBC(J)
          DO I = 1,NBR
            B(I,J) = A(IBR(I),JA)
          ENDDO
        ENDDO
C
C     Extract from transposed matrix
C
      ELSEIF(FA.EQ.'T') THEN
        DO J = 1,NBC
          JA = IBC(J)
          DO I = 1,NBR
            B(I,J) = A(JA,IBR(I))
          ENDDO
        ENDDO
      ELSE
        CALL QUIT('EXTELM:Unknown FA '//FA//' of matrix A!')
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lmnrep */
      SUBROUTINE LMNREP(IREP,NHKTA,KHKTA,LVALUE,MVALUE,NVALUE)
C***********************************************************************
C
C     Get out xyz-powers for a given irrep IREP.
C     Based on LMNVAL by tuh March 87
C
C     Written by T.Saue - Mar 24 1998
C***********************************************************************
#include "implicit.h"
#include "maxmom.h"
#include "xyzpow.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION LVALUE(KHKTA), MVALUE(KHKTA), NVALUE(KHKTA)
#include "symmet.h"
C
      ICOMP = 0
      DO 100 I = 1, KHKTA
      IF(ISYMAO(NHKTA,I).EQ.IREP) THEN
        ICOMP = ICOMP + 1
        LVALUE(ICOMP) = NHKTA - ISTEP(I)
        MVALUE(ICOMP) = MVAL(I)
        NVALUE(ICOMP) = NVAL(I)
      ENDIF
  100 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck qchold */
      SUBROUTINE QCHOLD(A,N,NZ,LRA,LCA,D,STOL,NEFF,JOB,IPIVOT)
C***********************************************************************
C
C     This routine performs upper triangle Cholesky decomposition of
C         * real symmetric matrices                     NZ = 1
C         * complex Hermitian matrices                  NZ = 2
C         * quaternionic Hermitian matrices             NZ = 4
C     On input/output A contains the lower triangle and diagonal of
C     the Hermitian matrix to be decomposed.
C     On output the upper triangle of the Cholesky matrix is stored
C     in the upper triangle of A, and the diagonal in D.
C     With pivoting(JOB.ne.0) the Choleksy matrix has been scrambled
C     according to the ordering indicated by IPIVOT.
C     
C     Written by T.Saue March 25 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION A(LRA,LCA,NZ),D(N),IPIVOT(*),P(4),IP(4)
C
C     Initialization
C
      DO IZ = 1,NZ
        IP(IZ) = IZ
      ENDDO
      NA = LRA*LCA
      NEFF = 0
      J = 1
      KMAX = J-1
C
C     Decomposition without pivoting
C
      IF(JOB.EQ.0) THEN
        IF(A(J,J,1).LE.STOL) RETURN
        NEFF = J
        D(J) = SQRT(A(J,J,1))
        DO I = 2,N
          A(J,I,1) = A(I,J,1)/D(1)
          DO IZ = 2,NZ
            A(J,I,IZ) = -A(I,J,IZ)/D(1)
          ENDDO
        ENDDO
        DO J = 2,N
          KMAX = J-1
          CALL QDOT(D(J),IP,1,KMAX,
     &               'H',IP,A(1,J,1),NA,NZ,1,
     &               'N',IP,A(1,J,1),NA,NZ,1)
          D(J) = A(J,J,1) - D(J)
          IF(D(J).LE.STOL) RETURN
          NEFF = J
          D(J) = SQRT(D(J))
          DO I = (J+1),N
            CALL QDOT(P,IP,NZ,KMAX,
     &               'H',IP,A(1,I,1),NA,NZ,1,
     &               'N',IP,A(1,J,1),NA,NZ,1)
            A(J,I,1) = (A(I,J,1)-P(1))/D(J)
            DO IZ = 2,NZ
              A(J,I,IZ) = (P(IZ)-A(I,J,IZ))/D(J)
            ENDDO
          ENDDO
        ENDDO
      ELSE
C
C     Cholesky decomposition with full pivoting
C
C...intialize pivot array
        DO I = 1,N
          IPIVOT(I) = I
        ENDDO
C...scan diagonal  of matrix A for largest element
        D(J) = A(J,J,1)
        JJ   = J
        DO K = 2,N
          IF(A(K,K,1).GT.D(J)) THEN
            D(J) = A(K,K,1)
            JJ   = K
          ENDIF
        ENDDO
        IF(D(J).LE.STOL) RETURN
        IF(JJ.NE.J) THEN
          NEFF = IPIVOT(JJ)
          IPIVOT(JJ) = IPIVOT(J)
          IPIVOT(J)  = NEFF
          JJ = NEFF
        ENDIF
        NEFF = J
C...generate elements of Cholesky matrix
        D(J) = SQRT(D(J))
        DO I = 2,N
          II = IPIVOT(I)
          IF(II.GE.JJ) THEN
            A(J,I,1) = A(II,JJ,1)/D(J)
            DO IZ = 2,NZ
              A(J,I,IZ) = -A(II,JJ,IZ)/D(J)
            ENDDO
          ELSE
            A(J,I,1) = A(JJ,II,1)/D(J)
            DO IZ = 2,NZ
              A(J,I,IZ) = A(JJ,II,IZ)/D(J)
            ENDDO
          ENDIF
        ENDDO
C...iterate
        DO J = 2,N
          KMAX = J-1
C...      generate diagonal elements of Cholesky matrix
C...      and find largest element
          JJ = IPIVOT(J)
          CALL QDOT(D(J),IP,1,KMAX,'H',IP,A(1,J,1),NA,NZ,1,
     &                             'N',IP,A(1,J,1),NA,NZ,1)
          D(J) = A(JJ,JJ,1) - D(J)
          DO K = (J+1),N
            KK = IPIVOT(K)
            CALL QDOT(D(K),IP,1,KMAX,'H',IP,A(1,K,1),NA,NZ,1,
     &                               'N',IP,A(1,K,1),NA,NZ,1)
            D(K) = A(KK,KK,1) - D(K)
            IF(D(K).GT.D(J)) THEN
              D(J) = D(K)
              JJ   = K
            ENDIF
          ENDDO
          IF(D(J).LE.STOL) RETURN
C..       if pivoting, interchange columns of previous rows of Cholesky matrix
          IF(JJ.NE.IPIVOT(J)) THEN
            DO IZ = 1,NZ
              CALL DSWAP(KMAX,A(1,J,IZ),1,A(1,JJ,IZ),1)
            ENDDO
            NEFF = IPIVOT(JJ)
            IPIVOT(JJ) = IPIVOT(J)
            IPIVOT(J)  = NEFF
            JJ = NEFF
          ENDIF
          NEFF = J
C...      generate elements of Cholesky matrix
          D(J) = SQRT(D(J))
          DO I = (J+1),N
            II = IPIVOT(I)
            CALL QDOT(P,IP,NZ,KMAX,
     &               'H',IP,A(1,I,1),NA,NZ,1,
     &               'N',IP,A(1,J,1),NA,NZ,1)
            IF(II.GE.JJ) THEN
              A(J,I,1) = (A(II,JJ,1)-P(1))/D(J)
              DO IZ = 2,NZ
                A(J,I,IZ) = (P(IZ)-A(II,JJ,IZ))/D(J)
              ENDDO
            ELSE
              A(J,I,1) = (A(JJ,II,1)-P(1))/D(J)
              DO IZ = 2,NZ
                A(J,I,IZ) = (P(IZ)+A(JJ,II,IZ))/D(J)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qdot */
      SUBROUTINE QDOT(P,IQPP,NZP,N,FX,IQPX,X,LRX,NZX,INCX,
     &                             FY,IQPY,Y,LRY,NZY,INCY)
C***********************************************************************
C
C     Performs generally real,complex or quaternion dot product
C     (determined by NZ)
C     FX/FY = 'N' : Normal vectors
C     FX/FY = 'H' : Hermitian conjugate
C
C     Written by T.Saue March 25 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0=0.0D0)
      DIMENSION X(LRX,NZX),Y(LRY,NZY),P(NZP),
     &          IQPX(NZX),IQPY(NZY),IQPP(NZP)
      CHARACTER  FX*1,FY*1
      DIMENSION IPY(4)
C
      EXTERNAL BDQALGEBRA
C
C     Initialization
C
      IFX = 1
      IF(FX.EQ.'H') IFX = 2
      IFY = 1
      IF(FY.EQ.'H') IFY = 2
C
C     Make pointer IPY
C     ================
C
      CALL IZERO(IPY,4)
      DO IZY = 1,NZY
        IPY(IQPY(IZY)) = IZY
      ENDDO
C
      DO IZP = 1,NZP
        IQP = IQPP(IZP)
        P(IZP) = D0
        DO 10 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQY = IQMULT(IQX,IQP,1)
          IZY = IPY(IQY)
          IF(IZY.EQ.0) GOTO 10
          FAC = IQSIGN(IQX,IFX,1)*IQSIGN(IQY,IFY,1)*IQPHASE(IQX,IQY,1)
          P(IZP) = P(IZP) + FAC*DDOT(N,X(1,IZX),INCX,Y(1,IZY),INCY)
 10     CONTINUE
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qdotm */
      SUBROUTINE QDOTM(P,IQPP,NZP,FX,IQPX,X,NX,LRX,LCX,NZX,INCX,
     &                            FY,IQPY,Y,NY,LRY,LCY,NZY,INCY,S,LDS)
C***********************************************************************
C
C     Performs generally real,complex or quaternion dot product
C     (determined by NZ)
C     FX/FY = 'N' : Normal vectors
C     FX/FY = 'H' : Hermitian conjugate
C     with REAL metric S
C
C     Written by T.Saue Nov 3 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0=0.0D0)
      DIMENSION X(LRX,LCX,NZX),Y(LRY,LCY,NZY),S(LDS,NY),P(NZP),
     &          IQPP(NZP),IQPX(NZX),IQPY(NZY)
      CHARACTER  FX*1,FY*1
      DIMENSION IPY(4)
C
      EXTERNAL BDQALGEBRA
C
C     Found out INC[XY] is not used! /hjaaj 4-Nov-2001
C
      IF (INCX .NE. 1 .OR. INCY .NE. 1) THEN
         CALL QUIT('QDOTM not implemented with INCX/Y .ne. 1')
      END IF
C
C     Initialization
C
      IFX = 1
      IF(FX.EQ.'H') IFX = 2
      IFY = 1
      IF(FY.EQ.'H') IFY = 2
C
C
C     Make pointer IPY
C     ================
C
      CALL IZERO(IPY,4)
      DO IZY = 1,NZY
        IPY(IQPY(IZY)) = IZY
      ENDDO
      DO IZP = 1,NZP
        IQP = IQPP(IZP)
        P(IZP) = D0
        DO 10 IZX = 1,NZX
          IQX = IQPX(IZX)
          IQY = IQMULT(IQX,IQP,1)
          IZY = IPY(IQY)
          IF(IZY.EQ.0) GOTO 10
          FAC = IQSIGN(IQX,IFX,1)*IQSIGN(IQY,IFY,1)*IQPHASE(IQX,IQY,1)
          P(IZP) = P(IZP)+FAC*DDOTM(NX,NY,LDS,X(1,1,IZX),S,Y(1,1,IZY))
 10     CONTINUE
      ENDDO
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnrm2 */
      FUNCTION QNRM2(N,NZ,X,NX,INCX)
C***********************************************************************
C
C     Gives squared norm of real,complex or quaternion vector
C     (determined by NZ)
C
C     Written by T.Saue March 25 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
      DIMENSION X(NX,NZ)
C
      QNRM2 = D0
      DO IZ = 1,NZ
        QNRM2 = QNRM2 + DDOT(N,X(1,IZ),INCX,X(1,IZ),INCX)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck qchols */
      SUBROUTINE QCHOLS(A,N,NEFF,NSIM,NZ,LRA,LCA,D,
     &                  B,LRB,LCB,X,LRX,LCX,JOB,IPIVOT,W)
C***********************************************************************
C
C     This routine solves the NSIM simultaneous equations
C         A(n,n,nz)*x(n,nsim,nz) = b(n,nsim,nz)
C     where A and D is output from the Cholesky decomposition
C     routine QCHOLD. The upper triagonal of A contains the 
C     upper triangular Cholesky matrix and the the diagonal
C     and lower triangle is the original matrix A
C         * a real symmetric matrix                     NZ = 1
C         * a complex Hermitian matrix                  NZ = 2
C         * a quaternion Hermitian matrix               NZ = 4
C     The vector D contains the diagonal elements of the
C     Cholesky matrix.
C
C     The system of linear equations is written as
C
C        Ax = U(dagger)Ux = U(dagger)y = b ; Ux = y
C
C     JOB.NE.0 means that full pivoting should be used....
C     W(N) is a work array needed for sorting in case of pivot.
C     
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0, D0=0.0D0)
      DIMENSION A(LRA,LCA,NZ),D(N),B(LRB,LCB,NZ),X(LRX,LCX,NZ),
     &          IPIVOT(*),W(*)
      DIMENSION P(4),IQA(4),IQB(4)
C
C     Initialization
C
      DO IZ = 1,NZ
        IQA(IZ) = IZ
        IQB(IZ) = IZ
      ENDDO
      NA = LRA*LCA
c     NB = LRB*LCB
      NX = LRX*LCX
C
C     First solve U(dagger)y = b by forward substitution
C
      DO J = 1,NSIM
C
C       In case of pivot, reorder coefficient array b
C
        IF(JOB.NE.0) THEN
          DO IZ = 1,NZ
            CALL DCOPY(N,B(1,J,IZ),1,W,1)
            DO I = 1,NEFF
              B(I,J,IZ) = W(IPIVOT(I))
            ENDDO
          ENDDO
        ENDIF
        DO IZ = 1,NZ
          X(1,J,IZ) = B(1,J,IZ)/D(1)
        ENDDO
        DO I = 2,NEFF
          KMAX = I-1
            CALL QDOT(P,IQA,NZ,KMAX,
     &               'H',IQA,A(1,I,1),NA,NZ,1,
     &               'N',IQA,X(1,J,1),NX,NZ,1)
          DO IZ = 1,NZ
            X(I,J,IZ) = (B(I,J,IZ)-P(IZ))/D(I)
          ENDDO
        ENDDO
      ENDDO
C
C     Solve Ux = y by backsubstitution
C
      DO J = 1,NSIM
        DO IZ = 1,NZ
          X(NEFF,J,IZ) = X(NEFF,J,IZ)/D(NEFF)
        ENDDO
      ENDDO
      DO I = (NEFF-1),1,-1
        P(2) = D1/D(I)
        P(1) = -P(2)
        NELM = NEFF-I
        KINT = I+1
        CALL QGEMM(1,NSIM,NELM,P(1),
     &             'N','N',IQA,A(I,KINT,1),LRA,LCA,NZ,
     &             'N','N',IQA,X(KINT,1,1),LRX,LCX,NZ,
     &                P(2),IQB,X(I,1,1),LRX,LCX,NZ)
      ENDDO
C
C     In case of pivot, reorder solution vectors
C
      IF(JOB.NE.0) THEN
        DO J = 1,NSIM
          DO IZ = 1,NZ
            CALL DCOPY(NEFF,X(1,J,IZ),1,W,1)
            DO I = 1,NEFF
              X(IPIVOT(I),J,IZ) = W(I)
            ENDDO
            DO I = (NEFF+1),N
              X(IPIVOT(I),J,IZ) = D0
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* qhmrst */
      SUBROUTINE QHMRST(A,N,NZ,LRA,LCA)
C***********************************************************************
C
C     In input is given a (real/complex/quaternion) Hermitian
C     matrix where the upper diagonal has been destroyed. Restore.
C
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LRA,LCA,NZ)
C
      DO J = 1,N
        DO I = 1,(J-1)
          A(I,J,1) = A(J,I,1)
        ENDDO
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,N
          DO I = 1,(J-1)
            A(I,J,IZ) = -A(J,I,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QSWAP(N,DX,ISYMX,FAC,DAGGER)
C***********************************************************************
C
C     Swap upper and lower part of quaternionic vector and scale
C     with a factor FAC.
C     If DAGGER then also complex conjugate vector.
C
C     Calls BLAS routine DSWAP.
C
C     Written by panor 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      PARAMETER ( DM1 = -1.0D0 )
C
      LOGICAL DAGGER
      DIMENSION DX(N,NZ)
C
      IF (N.LE.0) RETURN
C
      DO IZ=1,NZ
         CALL DSWAP(N/2,DX(1,IZ),1,DX(N/2+1,IZ),1)
         IF (DAGGER .AND. IPQTOQ(IZ,ISYMX-1).GT.1) 
     &        CALL DSCAL(N,DM1,DX(1,IZ),1)
      END DO
      CALL DSCAL(N*NZ,FAC,DX,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck curve */
      SUBROUTINE CURVE(NORDER,NPOINTS,A,B,X,Y,CHISQ,KVPT)
C***********************************************************************
C
C     Least square polynomial fit solving normal equations
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
      PARAMETER (D0=0.0D0)
      DIMENSION A(*),B(*),KVPT(*)
      DIMENSION X(NPOINTS),Y(NPOINTS)
      NDIM = NORDER + 1
C* Stabilize by subtracting mean from Y-vector:
      YMEAN = D0
      DO I = 1,NPOINTS
        YMEAN = YMEAN + Y(I)
      ENDDO
      YMEAN = YMEAN/NPOINTS
      DO I = 1,NPOINTS
        Y(I) = Y(I) - YMEAN
      ENDDO
C* Calculate upper triangle column packed matrix A :
      IOFF = 0
      DO 10 J = 1,NDIM
        DO 20 I = 1,J
          IOFF = IOFF + 1
          IORDER = I + J - 2
          A(IOFF) = XN(IORDER,NPOINTS,X)
   20 CONTINUE
   10 CONTINUE  
C* Calculate column vector b:
      DO 30 I = 1,NDIM
        IORDER = I - 1
        B(I) = XNY(IORDER,NPOINTS,X,Y)
   30 CONTINUE  
C* Solve linear equations
      CALL DSPSOL(NDIM,1,A,B,KVPT,INFO)      
      DO I = 1,NPOINTS
        Y(I) = Y(I) + YMEAN
      ENDDO
      B(1) = B(1) + YMEAN
C*    Calculate chisquare
      CHISQ = D0
      DO I = 1,NPOINTS
        DIFF  = Y(I) - POLVAL(NORDER,B,X(I))
        CHISQ = CHISQ + DIFF*DIFF
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Polsvd */
      SUBROUTINE POLSVD(ND,NP,A,B,X,Y,W,RV,CHISQ,ISKIP)
C***********************************************************************
C
C     Least square polynomial fit using singular value decomposition
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER(TOL=1.0D-16)
      DIMENSION A(NP,ND),B(NP),X(NP),Y(NP),W(ND),RV(ND)
C
      NORDER = ND-1
C* Stabilize by subtracting mean from Y-vector:
      YMEAN = D0
      DO I = 1,NP
        YMEAN = YMEAN + Y(I)
      ENDDO
      YMEAN = YMEAN/NP
      DO I = 1,NP
        Y(I) = Y(I) - YMEAN
      ENDDO
C     Set up design matrix A 
      CALL DCOPY(NP,D1,0,A(1,1),1)
      CALL DCOPY(NP,X,1,A(1,2),1)
      DO J = 3,ND
        DO I = 1,NP
          A(I,J)=A(I,J-1)*A(I,2)
        ENDDO
      ENDDO
C     Set up response vector B
      CALL DCOPY(NP,Y,1,B,1)
C     Perform singular value decomposition
      CALL MINFIT(NP,NP,ND,A,W,1,B,IERR,RV)
C     Edit the singular values: Note that this is done to B !!
      WA = D0
      DO I = 1,ND
        IF(W(I).GT.WA) WA=W(I)
      ENDDO
      WA = TOL*WA
      ISKIP = 0
      DO I = 1,ND
        IF(W(I).LT.WA) THEN
          B(I)=D0
          ISKIP = ISKIP + 1
        ELSE
          B(I)=B(I)/W(I)
        ENDIF
      ENDDO
C     Backsubsitution to get parameters      
      CALL DGEMV('N',ND,ND,D1,A,NP,B,1,D0,W,1)
      DO I = 1,NP
        Y(I) = Y(I) + YMEAN
      ENDDO
      W(1) = W(1) + YMEAN
C     Calculate chisquare
      CHISQ = D0
      DO I = 1,NP
        DIFF  = Y(I) - POLVAL(NORDER,W,X(I))
        CHISQ = CHISQ + DIFF*DIFF
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE NEWRAP(NORDER,DTOL,B,XVAL,NITER,IERR)
#include "implicit.h"
#include "priunit.h"
      DIMENSION B(NORDER + 1)
      IERR = 0
      NSTEP = 0
   10 CONTINUE
      NSTEP = NSTEP + 1
      XOLD = XVAL  
      XFAC = POLDER(NORDER,B,XOLD)/POL2DER(NORDER,B,XOLD)
      XVAL = XOLD - XFAC
      IF((ABS((XVAL-XOLD)).GT.DTOL).AND.(NSTEP.LE.NITER)) GOTO 10
      IF(NSTEP.GT.NITER) THEN
        IERR = 1
      ENDIF
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION XN(NORDER,NPOINTS,X)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION X(NPOINTS)
      XN = D0
      DO 10 I = 1,NPOINTS
        XN = XN + X(I)**NORDER
   10 CONTINUE    
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION XNY(NORDER,NPOINTS,X,Y)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION X(NPOINTS),Y(NPOINTS)
      XNY = D0
      FAC = NORDER
      DO 10 I = 1,NPOINTS
        XNYBUF = Y(I)
        XBUF = X(I)
        DO 20 J = 1,NORDER
          XNYBUF = XNYBUF*XBUF
   20   CONTINUE    
C      XNY = XNY + XNYBUF
      XNY = XNY + Y(I)*(X(I)**FAC)
   10 CONTINUE    
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLVAL(NORDER,B,XVAL)
#include "implicit.h"
      DIMENSION B(NORDER+1)
      POLVAL = B(1)
      XBUF = 1
      DO 10 I = 2,(NORDER+1)
        XBUF = XBUF*XVAL
        POLVAL = POLVAL + XBUF*B(I)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLDER(NORDER,B,XVAL)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.1) THEN
        POLDER = D0
        RETURN
      ENDIF
      POLDER = B(2)
      DO 10 I = 2,NORDER
        POLDER = POLDER + (XVAL**(I-1))*I*B(I+1)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POL2DER(NORDER,B,XVAL)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D2 = 2.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.2) THEN
        POL2DER = D0
        RETURN
      ENDIF
      POL2DER = D2*B(3)
      DO 10 I = 3,NORDER
        POL2DER = POL2DER + (XVAL**(I-2))*I*(I-1)*B(I+1)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POL3DER(NORDER,B,XVAL)
#include "implicit.h"
      PARAMETER(D0 = 0.0D0,D6 = 6.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.3) THEN
        POL3DER = D0
        RETURN
      ENDIF
      POL3DER = D6*B(4)
      DO 10 I = 4,NORDER
        POL3DER = POL3DER + (XVAL**(I-3))*(I-2)*(I-1)*I*B(I+1)
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION POLNDER(NORDER,B,XVAL,N)
#include "implicit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION B(NORDER+1)
      IF(NORDER.LT.N) THEN
        POLNDER = D0
        RETURN
      ENDIF
      FAC = D1
      DO J = 2,N
        FAC = FAC*J
      ENDDO
      POLNDER = FAC*B(N+1)
      DO 10 I = N+1,NORDER
        FAC = I
        DO J = 1,N-1
          FAC = FAC*(I-J)
        ENDDO
        POLNDER = POLNDER + (XVAL**(I-N))*B(I+1)*FAC
   10 CONTINUE     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION MIND(NPOINTS,Y)
#include "implicit.h"
      DIMENSION Y(NPOINTS)
      YVAL = Y(1)
      DO 10 I = 2,NPOINTS
        IF(Y(I).LE.YVAL) THEN
          YVAL = Y(I)
          MIND = I
        ENDIF
   10 CONTINUE
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qtrans */
      SUBROUTINE QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  BUF,LBUF,IPRINT)
C***********************************************************************
C
C     This routine performs unitary transformations indicated by TYP:
C
C       TYP = AOMO : AO-to-MO-transformation  FMO = (C+)FAOC
C       TYP = MOAO : MO-to-AO-transformation  FAO = CFMO(C+)
C
C     FADD gives the possibility of adding the results
C
C     TREV indicates symmetry of F-matrix under timereversal
C
C       TREV = 'S' - symmetric     (t =  1)
C       TREV = 'A' - anti-symmtric (t = -1)
C
C     The unitary C-matrix has a time-symmetric structure.
C
C     The AO-to-MO-transformation can be set up quaternionically as
C
C       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
C
C       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
C
C     and the MO-to-Ao-transformation can be set up as:
C
C       FAO = [(Ca)+t(Cb)j][(FMOa)+(FMOb)j][(Ca+)-(Cb)j]
C
C       F(NORB,NORB,NZ) --->  F(NBAS,NBAS,NZ)
C
C     Written by T.Saue,January 1995
C     Last revision: 29 January 1997 : Luuk Visscher
C**********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      CHARACTER TREV*1,TRA*1,TYP*4
      DIMENSION FAO(LRAO,LCAO,*),IQAO(*),
     &          FMO(LRMO,LCMO,*),IQMO(*),
     &          TM1(LR1,LC1,*),IQTM1(*),
     &          TM2(LR2,LC2,*),IQTM2(*),
     &          BUF(LBUF)
      DIMENSION IQW(4)
C
      CALL QENTER('QTRANS')
C
C     Add the results on to the contents of the MO matrix ?
C
#ifdef DEBUG
      print*, 'QTRANS: The AO matrix'
      CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IQAO,LUPRI)
      print*, 'QTRANS: The transformation matrix TM1'
      CALL PRQMAT(TM1,NRAO,NCMO,LR1,LC1,NZTM1,IQTM1,LUPRI)
      print*, 'QTRANS: The transformation matrix TM2'
      CALL PRQMAT(TM2,NRAO,NCMO,LR2,LC2,NZTM2,IQTM2,LUPRI)
#endif
C
C     Time reversal symmetry of F-matrix
C
      IF    (TREV.EQ.'S') THEN
        TRA = 'N'
      ELSEIF(TREV.EQ.'A') THEN
        TRA = 'I'
      ELSE
        CALL QUIT('QTRANS: Unknown keyword TREV '//TREV)
      ENDIF
C
C     *******************************************************
C     *****  AO - to - MO -transformation : F' = (C+)FC *****
C     *******************************************************
C
      IF(TYP.EQ.'AOMO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF (NCAO.GE.NRAO) THEN
C
C         Perform the two steps:
C           1. W  = FC
C           2. F' = (C+)W
C
          CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = FC
C         ------------------------------------
C
          CALL QGEMM(NRAO,NCMO,NCAO,D1,
     &              'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &              'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &               D0,IQW,BUF,NRAO,NCMO,NZW)
C
C         Second part of transformation: H = (C+)W
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCMO,NRAO,D1,
     &               'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &               'N','N',IQW,BUF,NRAO,NCMO,NZW,
     &               FADD,IQMO,FMO,LRMO,LCMO,NZMO)
        ELSE
C
C         Perform the two steps:
C           1. W  = (C+)F
C           2. F' = WC
C
          CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF

C
C         First part of transformation: W = (C+)F
C         ----------------------------------------
C
          CALL QGEMM(NRMO,NCAO,NRAO,D1,
     &       'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,
     &       'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &        D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: H = WC
C         -------------------------------------
C
          CALL QGEMM(NRMO,NCMO,NCAO,D1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
C
        ENDIF
C        IF(IPRINT.GE.5) THEN
C          CALL HEADER('QTRANS: AO-to-MO-transformed Fock matrix',-1)
C          CALL PRQMAT(FMO,NRMO,NCMO,LRMO,LCMO,NZMO,IPQTOQ(1,0,1),LUPRI)
C        ENDIF
C
C     *******************************************************
C     *****  MO - to - AO -transformation : F' = CF(C+) *****
C     *******************************************************
C
      ELSEIF(TYP.EQ.'MOAO') THEN
C
C       Look for the most efficient transformation
C       ------------------------------------------
C
        IF(NCMO.GE.NRMO) THEN
          CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = F(C+)
C
          CALL QGEMM(NRMO,NCAO,NCMO,D1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       D0,IQW,BUF,NRMO,NCAO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL QGEMM(NRAO,NCAO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ELSE
          CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (LUPRI,'(/A,I10/A,I10)')
     &   ' >>> QTRANS error, need work space    ',NBUF,
     &   '                   current work space ',LBUF
          ENDIF
C
C         First part of transformation: W = CF
C
          CALL QGEMM(NRAO,NCMO,NRMO,D1,
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,
     &       D0,IQW,BUF,NRAO,NCMO,NZW)
C
C         Second part of transformation: H = CW
C
          CALL QGEMM(NRAO,NCAO,NCMO,D1,
     &       'N','N',IQW ,BUF,NRAO,NCMO,NZW,
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ENDIF
C       IF(IPRINT.GE.5) THEN
C          CALL HEADER('QTRANS: MO-to-AO-transformed Fock matrix',-1)
C          CALL PRQMAT(FAO,NRAO,NCAO,LRAO,LCAO,NZAO,IPQTOQ(1,0,1),LUPRI)
C        ENDIF
      ELSE
        CALL QUIT('FTRANS: Unknown keyword TYP '//TYP)
      ENDIF
C
      CALL QEXIT('QTRANS')
      RETURN
C
      END
