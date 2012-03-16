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
      use quaternion_algebra
#include "implicit.h"
C
      LOGICAL LBUF(4)
      DIMENSION IQPA(NZA),IQPB(NZB),IQPC(*)
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
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
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
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
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

#ifdef LAPACK_QDIAG
C     Use the fast DSYEVR for real case
      IF (NZ.EQ.1.AND.MATZ.EQ.1) THEN
            DDUM = 0.0D0
            IDUM = 0
            WORKDUM = 0.0D0
            IWORKDUM = 0
C     eigval and eigvec
c     first let lapack estimate memory use
CMI ...  I have to initialize it for intel mkl-lapack
            INFO = 0
            IERR = 0
            CALL DSYEVR('V','A','U',N,A,LRA,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &           M,EIG,VEC,LRV,IDUM,WORKDUM,-1,IWORKDUM,-1,INFO)
            LLWORK = NINT(WORKDUM)
            LILWORK = IWORKDUM
            CALL MEMGET('REAL',KLWORK, LLWORK, WORK, KFREE, LFREE)
            CALL MEMGET('INTE',KLIWORK, LILWORK, WORK, KFREE, LFREE)
            CALL MEMGET('INTE',KISUPPZ, 2*N, WORK, KFREE, LFREE)
            CALL DSYEVR('V','A','U',N,A,LRA,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &           M,EIG,VEC,LRV,WORK(KISUPPZ),WORK(KLWORK),LLWORK,
     &           WORK(KLIWORK),LILWORK,IERR)
C     DSYEVR does not guarantee a sorting order, but I have never seen it
C     do anything other than ascending, so the sorting should be a noop.
            CALL ORDER3(VEC,EIG,LRV,N,N,1)
            CALL QEXIT('QDIAG')
            RETURN
      ENDIF
#endif
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

      SUBROUTINE QDIAG90(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,IERR)
C     Like QDIAG but uses f90 memory allocation internally.
      use memory_allocator
#include "implicit.h"
      real(8), allocatable :: work(:)

      lenid = 3 ! lenid used in MEMGET in QDIAG
      lwork = 4*lenid + 2*N
      if(nz.gt.1) lwork = lwork + N*NZ

#ifdef LAPACK_QDIAG
      if(matz.eq.1)then
        DDUM = 0.0D0
        IDUM = 0
        WORKDUM  = 0.0D0
        IWORKDUM = 0
c       first let lapack estimate memory use
        INFO = 0
        CALL DSYEVR('V','A','U',N,A,LRA,DDUM,DDUM,IDUM,IDUM,0.0D0,
     &              M,EIG,VEC,LRV,IDUM,WORKDUM,-1,IWORKDUM,-1,INFO)
        LLWORK  = NINT(WORKDUM)
        LILWORK = IWORKDUM
        lwork = LLWORK + LILWORK + 2*N + 6*lenid
      end if
#endif
      if(lwork.le.0)then
        call quit(' *** error in QDIAG90: buffer space allocation not
     & possible for length <= 0.***')
      end if
      call alloc(work,lwork, id='QDIAG90-buffer-space')
      work = 0
      call qdiag(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,work,lwork,IERR)
      call dealloc(work, id='QDIAG90-buffer-space')
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
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
      CHARACTER FA*1, TA*1
      DIMENSION A(N,NZA),IQPA(*)
C
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
