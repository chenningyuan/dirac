      SUBROUTINE COPDIA(A,VEC,NDIM,IPACK)
*
* Copy diagonal of matrix A into vector VEC
*
*   IPACK = 0 : Full matrix
*   IPACK = 1 : Lower triangular matrix
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION A(*),VEC(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          VEC(I) = A((I-1)*NDIM+I)
  100   CONTINUE
      ELSE
        DO 200 I = 1, NDIM
          VEC(I) = A(I*(I+1)/2)
  200   CONTINUE
      END IF
*
      RETURN
      END
      SUBROUTINE SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
*
* Calculate square root of positive definite symmetric matrix A
* if(ITASK .EQ. 2 ) Inverted square root matrix is also calculated
* In case of singularities in A A -1/2 is defined to have the same
* singularity
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
*
      DIMENSION A(NDIM,NDIM)
      DIMENSION ASQRT(NDIM,NDIM),AMSQRT(NDIM,NDIM)
      DIMENSION SCR(*)
* Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2
      KLFREE = 1
*
      KLASYM = KLFREE
      KLAVAL = KLASYM
      KLFREE = KLASYM + NDIM*(NDIM+1)/2
*
      KLAVEC = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      NTEST = 0
*
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(A,SCR(KLASYM),1,NDIM,NDIM)
      CALL EIGEN(SCR(KLASYM),SCR(KLAVEC),NDIM,0,1)
      CALL COPDIA(SCR(KLASYM),SCR(KLAVAL),NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT(SCR(KLAVAL),NDIM,1,NDIM,1)
      END IF
*. Check for negative eigenvalues
*
      DO 100 I = 1,NDIM
        SCR(KLAVAL-1+I) = SQRT(SCR(KLAVAL-1+I))
  100 CONTINUE
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(ASQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
*
      IF(ITASK .EQ. 2 ) THEN
        DO 200 I = 1,NDIM
          IF(SCR(KLAVAL-1+I) .GT. 1.0D-13 ) then
            SCR(KLAVAL-1+I) = 1.0D0/SCR(KLAVAL-1+I)
          ELSE
            SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
          END IF
  200   CONTINUE
        CALL XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
      END IF
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Info from SQRTMT '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' Input matrix to SQRTMT '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Square root of matrix '
        CALL WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM)
        IF(ITASK .EQ. 2 ) THEN
          WRITE(6,*) ' Inverse square root of matrix '
          CALL WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM)
        END IF
      END IF
*
      RETURN
      END

      SUBROUTINE XDIAXT(XDX,X,DIA,NDIM,SCR)
*
* Obtain XDX = X * DIA * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XDX(NDIM,NDIM)
      DIMENSION X(NDIM,NDIM),DIA(NDIM)
      DIMENSION SCR(NDIM,NDIM)
*
* DIA * X(transposed)
      DO 100 I=1,NDIM
        CALL COPVEC(X(1,I),SCR(1,I),NDIM)
        CALL SCALVE(SCR(1,I),DIA(I),NDIM)
  100 CONTINUE
* X * DIA * X(Transposed)
      CALL MATML4(XDX,X,SCR,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,2)
*
      RETURN
      END

      SUBROUTINE EXTKOP
*
* Use Extended Koopmaans' theorem to obtain ionization energies
*
*
* Jeppe Olsen , LUCIA version adapeted from LUCAS, Nov 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "lucinp.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2


*
      WRITE(6,*) '        *****************************************'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *  Welcome to the magic(?) world of the *'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *      Extended Koopmans Theorem        *'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *       Version of Dec 12 1996          *'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *****************************************'
C     KFREEO = KFREE
      IDUM = 0
      WRITE(6,*) ' NTOOB = ', NTOOB
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXTKOP')
      CALL MEMMAN(KDEN1SM,2*NINT1,'ADDL ',2,'DEN1SM')
      LENGTH = NTOOB ** 2
      CALL MEMMAN(KLSSQ,LENGTH,'ADDL  ',2,'SSQ   ')
      CALL MEMMAN(KLSSQI,LENGTH,'ADDL  ',2,'SSQI  ')
      CALL MEMMAN(KLFP  ,LENGTH,'ADDL  ',2,'FP    ')
      CALL MEMMAN(KLFEXP,LENGTH,'ADDL  ',2,'KLFEXP')
      CALL MEMMAN(KLSCR ,4*LENGTH,'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KFOCK ,LENGTH ,'ADDL  ',2,'KFOCK ')

*. Density matrices and integrals assumed in place
*
*. Construct Fock matrix over all orbitals
      CALL FOCK_MAT(WORK(KFOCK))
*. Multiply Fock matrix with a factor of 2 to improve results
*. - a deeper argument does exist
COLD  TWO = 2.0D0
*. will give problems if we reduve size of FOCL matrix
*   to symmetry blocked form
COLD  CALL SCALVE(WORK(KFOCK),TWO,NTOOB**2)
      WRITE(6,*) ' Fock matrix '
      WRITE(6,*) ' ============'
      ISM = 0
      CALL APRBLM2(WORK(KFOCK),NTOOBS,NTOOBS,NSMOB,ISM)
*. Obtain density matrix in symmetry blocks
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFFM=1
          IOFFO=1
        ELSE
          IOFFM = IOFFM + NTOOBS(ISM-1)** 2
          IOFFO = IOFFO + NTOOBS(ISM-1)
        END IF
        DO IOB = IOFFO,IOFFO + NTOOBS(ISM)-1
           LOB = NTOOBS(ISM)
           DO JOB = IOFFO,IOFFO + NTOOBS(ISM)-1
*. Corresponding type indeces
             IOBP = IREOST(IOB)
             JOBP = IREOST(JOB)
             WORK(KDEN1SM-1+IOFFM-1+(JOB-IOFFO)*LOB+IOB-IOFFO+1)
     &     = WORK(KRHO1-1+(IOBP-1)*NTOOB+JOBP)
           END DO
        END DO
        WRITE(6,*) ' Symmetry packed one body density matrix '
        IMATSM = 0
        CALL APRBLM2(WORK(KDEN1SM),NTOOBS,NTOOBS,NSMOB,IMATSM)
      END DO
*. Loop over Symmetries and  and solve F C = E S C for each symmetry
      DO 200 ISM = 1, NSMOB
        N = NTOOBS(ISM)
        IF(ISM.EQ.1) THEN
           IBASE = 1
           IBASE2 = 1
           IBASE2P = 1
        ELSE
           IBASE = IBASE + NTOOBS(ISM-1)
           IBASE2P = IBASE2P + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
           IBASE2 = IBASE2 + NTOOBS(ISM-1)**2
        END IF
*. Pack block of Fock matrix to complete form
C       CALL TRIPAK(WORK(KLFEXP),WORK(KFOCK-1+IBASE2P),2,N,N)
*. Already in complete form so
        CALL COPVEC(WORK(KFOCK-1+IBASE2),WORK(KLFEXP),N**2)
        WRITE(6,*) ' Outpacked block of FOCK '
        CALL WRTMAT(WORK(KLFEXP),n,n,n,n)

*. Solve FC = eSC with Ct S C = 1 as F'C' = e S'C'( S is one body density)
* with F' = S-(1/2)F S(-1/2), C' = S(1/2)C
*. S(-1/2),S(1/2)
        CALL SQRTMT(WORK(KRHO1+IBASE2-1),N,2,WORK(KLSSQ),WORK(KLSSQI),
     &              WORK(KLSCR))
        write(6,*) ' S-1/2 matrix '
        call wrtmat(work(klssqi),n,n,n,n)
*. S(-1/2)F S(-1/2)
        CALL TRNMAT(WORK(KLFEXP),WORK(KLSSQI),WORK(KLSCR),N,N)
        write(6,*) ' Transformed matrix'
        call wrtmat(work(KLFEXP),n,n,n,n)
*. Diagonalize
        CALL TRIPAK(WORK(KLFEXP),WORK(KLSCR),1,N,N)
        CALL EIGEN(WORK(KLSCR),WORK(KLFP),N,0,1)
        WRITE(6,*)
     &  ' Eigenvalues of Extended Koopmans Equation for ISM = ',ISM
        DO 80 I = 1, N
          WRITE(6,'(F20.15)') WORK(KLSCR-1+I*(I+1)/2)
   80   CONTINUE
  200 CONTINUE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXTKOP')
      RETURN
      END
      SUBROUTINE PTEKT(LU0,LUN,N,ISM,ISPC)
*
* Perturbation analysis of the Extended Koopmans' theorem
* through general order.
*
* It is assumed that this calculation is preceded  by
* a call to the perturbation routine to obtain the
* wave function corrections to the neutral state.
*
* Input
*       LUN : File containing wave function corrections
*       LU0 : File containing reference wave funcrtion
*         N : Max order of expansion
*      ISM : Symmetry of reference state
*      ISPC : Space of referencestate
*
*. Specification of zero order states to be analyzed is
* supplied through crun :
*    NPTEKT : Number of EKT per sym to be analyzed
*    LPTEKT   EKT's to be analyzed : oorbital and symmetry of zero order state
* Jeppe and Dage, Nov11, 1995
* Completed       Jan31  1997 ( A long time for a one-week project!)
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
*. Should not be called with ICISTR = 1
*. Nor with inactive orbitals
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "lucinp.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      NTEST = 100
*
      WRITE(6,*)
      WRITE(6,*) ' ======================== '
      WRITE(6,*) ' PTEKT is now in CONTROL '
      WRITE(6,*) ' ======================== '
      WRITE(6,*)

* a bit on files :
* LUSC36 is LUN.
* Two additional scratch files to be used are  LUSC1 and LUSC2
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PTEKT ')
*
*
*     ========================
* 1 : Local memory allocation
*     ========================
*
*. Allocate space for two vector chunks
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(IPRCIX.GE.1)
     &WRITE(6,*) ' PTEKT  : LCSBLK = ', LCSBLK
      LSCR1 = MAX(LCSBLK,LSCR1)
      CALL MEMMAN(KLVEC1,LSCR1,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,LSCR1,'ADDL  ',2,'VEC2  ')
* one-body Density matrices through order n
      NMAT = N+1
      LENGTH = NMAT * NTOOB ** 2
      CALL MEMMAN(KLDEN1,LENGTH,'ADDL  ',2,'DENN1 ')
*. and Two-body density matrices
       LENGTH = NMAT * NTOOB ** 2*(NTOOB**2+1)/2
      CALL MEMMAN(KLDEN2,LENGTH,'ADDL  ',2,'DENN2 ')
*. Fock matrices through order n
      LENGTH = NMAT * 2*NINT1
      CALL MEMMAN(KLFN,LENGTH,'ADDL  ',2,'F(N)  ')
*. A scratch matrix ( not a nice thing to say about a matrix )
      LENGTH =  2*NINT1
      CALL MEMMAN(KLFSCR,LENGTH,'ADDL  ',2,'FSCR  ')
*. Space for Fock matrices, density matrix beloinging to a
*  given symmetry
      LENGTH = NMAT * 2*NINT1
      CALL MEMMAN(KLFNS,LENGTH,'ADDL  ',2,'F(N)S ')
      LENGTH = NMAT * 2*NINT1
      CALL MEMMAN(KLDNS,LENGTH,'ADDL  ',2,'D(N)S ')
*. Wave function corrections and energy corrections
      LENGTH = NMAT * NTOOB
      CALL MEMMAN(KLCN,LENGTH,'ADDL  ',2,'C(N)  ')
      CALL MEMMAN(KLEN,N     ,'ADDL  ',2,'E(N)  ')
*. and three scratch vectors for solver
      CALL MEMMAN(KLSCR1,NTOOB,'ADDL  ',2,'SCR1 ')
      CALL MEMMAN(KLSCR2,NTOOB,'ADDL  ',2,'SCR2 ')
      CALL MEMMAN(KLSCR3,NTOOB,'ADDL  ',2,'SCR3 ')

*
* ===============================================
* 2 : Construct density matrices through order N
* ===============================================
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      DO K = 0, N
        CALL PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
     &              WORK(KLDEN1+(K-0)*LRHO1),
     &              WORK(KLDEN2+(K-0)*LRHO2),LUSC1,LUSC2)
      END DO
*
* =====================================================
* 3 : Construct n'th order contribution to Fock matrix
* =====================================================
*
*. Contains two terms : Zero order Hamiltonian with N'th order densities
*                       Perturbation with N-1'th order densities
*. Here is where things would be messed up with inactive orbitals
*
*. The Zero'th order operator is assumed in place in KFI
      WRITE(6,*) ' Operator in KFI '
      CALL APRBLM2(WORK(KFI),NTOOBS,NTOOBS,NSMOB,1)
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB **2*(NTOOB**2+1)/2
      LFOCK = 2*NINT1
*
      ONE = 1.0D0
      ONEM = -1.0D0
      DO K = 0, N
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLFN+(K-0)*LFOCK),ZERO,LFOCK)
        IF(K.GT.0) THEN
*. Full Hamiltonian with K-1 order density
          CALL COPVEC(WORK(KLDEN1+(K-1-0)*LRHO1),WORK(KRHO1),LRHO1)
          CALL COPVEC(WORK(KLDEN2+(K-1-0)*LRHO2),WORK(KRHO2),LRHO2)
          CALL FOCK_MAT(WORK(KLFSCR),2)
          CALL COPVEC(WORK(KLFSCR),WORK(KLFN+(K-0)*LFOCK),LFOCK)
*. Subtract Zero order operator with K-1 order densities
          CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
          CALL FOCK_MAT(WORK(KLFSCR),1)
          CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
          CALL VECSUM(WORK(KLFN+(K-0)*LFOCK),
     &                WORK(KLFN+(K-0)*LFOCK),
     &                WORK(KLFSCR),ONE,ONEM,LFOCK)
        END IF
*. Zero order hamiltonian with k'th order densities
        CALL COPVEC(WORK(KLDEN1+(K-0)*LRHO1),WORK(KRHO1),LRHO1)
        WRITE(6,*) ' Density in RHO1 '
        CALL WRTMAT(WORK(KRHO1),NTOOB,NTOOB,NTOOB,NTOOB)
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        CALL FOCK_MAT(WORK(KLFSCR),1)
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        CALL VECSUM(WORK(KLFN+(K-0)*LFOCK),
     &              WORK(KLFN+(K-0)*LFOCK),
     &              WORK(KLFSCR),ONE,ONE,LFOCK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 'Correction to Fock matrix of order =',K
          CALL APRBLM2(WORK(KLFN+(K-0)*LFOCK),
     &                 NOCOBS,NOCOBS,NSMOB,0)
        END IF
*
      END DO
*
* =================================
* 4 Perturbation expansions of EKT
* =================================
*
      DO IEKT = 1, NPTEKT
        WRITE(6,*) ' Zero order w.f. for EKT '
        WRITE(6,*) ' ========================'
        IORB = LPTEKT(1,IEKT)
        IOSYM = LPTEKT(2,IEKT)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry = ', IOSYM
        WRITE(6,*) ' Orbital  = ', IORB
        WRITE(6,*)
*. Obtain density and fock matrices for just this sym
        DO JSYM = 1, IOSYM
          IF(JSYM .EQ. 1 ) THEN
            IOFF = 1
            IOFFS = 1
          ELSE
             IOFF = IOFF + NOCOBS(JSYM-1)**2
             IOFFS= IOFFS+ NOCOBS(JSYM-1)*(NOCOBS(JSYM-1)+1)/2
          END IF
        END DO
        LORB = NOCOBS(IOSYM)
        WRITE(6,*) ' IOFF = ', IOFF
        LSYMBLK = LORB**2
        DO I =0, N
          CALL COPVEC(WORK(KLFN-1+(I-0)*LFOCK+IOFF),
     &                WORK(KLFNS+(I-0)*LSYMBLK),LSYMBLK)
C  TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
C. In general the one-body density is not symmetry packed
C  I'll assume that we just a single symmetry active !!
          WRITE(6,*) ' Density copied from '
          CALL WRTMAT(WORK(KLDEN1-1+(I-0)*LRHO1+IOFF),
     &                LORB,LORB,LORB,LORB)
          CALL COPVEC(WORK(KLDEN1-1+(I-0)*LRHO1+IOFF),
     &                WORK(KLDNS+(I-0)*LSYMBLK),LSYMBLK)
        END DO
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Symmetry blocks of F and S '
          DO I = 0, N
            WRITE(6,*) ' Fock matrix of order ',I
            CALL WRTMAT(WORK(KLFNS+(I-0)*LSYMBLK),LORB,LORB,LORB,LORB)
            WRITE(6,*) ' density matrix of order ',I
            CALL WRTMAT(WORK(KLDNS+(I-0)*LSYMBLK),LORB,LORB,LORB,LORB)
          END DO
        END IF
*. Initialize zero order vector
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLCN),ZERO,LORB)
        WORK(KLCN-1+IORB) = 1.0D0
*. And transfer control and responsibility
        CALL GENPERT(WORK(KLFNS),WORK(KLDNS),LORB,N,WORK(KLEN),
     &               WORK(KLCN),WORK(KLSCR1),WORK(KLSCR2),
     &               WORK(KLSCR3))
      END DO
*.Finito
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'PTEKT ')
*
      RETURN
      END
C       PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
C    &              WORK(KLDEN1+(N-0)*LRHO1),
C    &              WORK(KLDEN2+(N-0)*LRHO2),LUSC1,LUSC2)
      SUBROUTINE PERTDN
     &(N,LU0,LUN,ISM,ISPC,VEC1,VEC2,RHO1N,RHO2N,LUSC1,LUSC2)
*
* Construct one body density matrix of order N
*
*      Jeppe + Dage, Nov. 11 1995
*                    Debugged Jan 31 '97
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
*. Should not be called with ICISTR = 1
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
C     INCLUDE 'clunit.inc'
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Output
      DIMENSION RHO1N(*),RHO2N(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'PERTDN')
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      CALL MEMMAN(KLDEN1,LRHO1,'ADDL  ',2,'KLDEN1')
      CALL MEMMAN(KLDEN2,LRHO2,'ADDL  ',2,'KLDEN2')
*
      LBLK = -1
      ZERO = 0.0D0
      CALL SETVEC(RHO1N,ZERO,LRHO1)
      CALL SETVEC(RHO2N,ZERO,LRHO2)
*
      DO L = 0, N
C?      write(6,*) ' Will load next pair of vectors '
        NMINL = N - L
CTOBE   IF(L.LE.NMINL) THEN
*. put correction vector L and NMINL on LUSC1 and LUSC2, respectively
          IF(L.EQ.0) THEN
             CALL COPVCD(LU0,LUSC1,VEC1,1,LBLK)
          ELSE
             CALL SKPVCD(LUN,L-1,VEC1,1,LBLK)
             CALL REWINO(LUSC1)
             CALL COPVCD(LUN,LUSC1,VEC1,0,LBLK)
          END IF
*
          IF(NMINL.EQ.0) THEN
             CALL COPVCD(LU0,LUSC2,VEC1,1,LBLK)
          ELSE
             CALL SKPVCD(LUN,NMINL-1,VEC1,1,LBLK)
             CALL REWINO(LUSC2)
             CALL COPVCD(LUN,LUSC2,VEC1,0,LBLK)
          END IF
C?      write(6,*) ' next pair of vectors loaded '
* Do the densi
          LEQR = 0
C DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR)
          I12 = 2
C?        WRITE(6,*) ' Calling DENSI2 '
          CALL DENSI2(I12,WORK(KLDEN1),WORK(KLDEN2),VEC1,VEC2,
     &                LUSC1,LUSC2,0)
          WRITE(6,*) ' Home from DENSI2 '
          WRITE(6,*) ' Densities fresh from DENSI2 '
C?        WRITE(6,*) ' NTOOB = ', NTOOB
          CALL WRTMAT(WORK(KLDEN1),NTOOB,NTOOB,NTOOB,NTOOB)
          CALL PRSYM(WORK(KLDEN2),NTOOB**2)
*
CTOBE     IF(L.NE.NMINL) THEN
*. The matrix <L! E !NMINL> was calculated, add <NMINL! E ! L>
*. as simple transposition
CTOBE        CALL TRPAD(WORK(KLDEN),ONE,NTOOB)
CTOBE     END IF
          ONE = 1.0D0
          CALL VECSUM(RHO1N,RHO1N,WORK(KLDEN1),ONE,ONE,LRHO1)
          CALL VECSUM(RHO2N,RHO2N,WORK(KLDEN2),ONE,ONE,LRHO2)
CTOBE   END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Density matrix of order in perturbation ', N
        WRITE(6,*) ' ==========================================='
        WRITE(6,*)
        WRITE(6,*) ' One-body density '
        WRITE(6,*) ' ================ '
        CALL WRTMAT(RHO1N,NTOOB,NTOOB,NTOOB,NTOOB)
        WRITE(6,*) ' Two-body density '
        WRITE(6,*) ' ================ '
        CALL PRSYM(RHO2N,NTOOB**2)
      END IF
*
      RETURN
      END
      SUBROUTINE GENPERT(F,S,NDIM,NORD,EN,C)
*
* The matrices of a  generalized eigenvalue problem FC = E SC
* are separated into orders
*
* F = Sum(k=0,NORD) F(K)
* S = Sum(K=0,NORD) S(K)
*
*. Obtain corrections to energy and wawe functions
*
*. The normalization condition used is C(K)T S(0) C(0) = 0
*
* The energy corrections become
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
* and the wave function corrections are
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(I=0,N-K)E(N-K-L)S(L)C(K))
*
* The zero order matrices F(0),S(0) are assumed diagonal
*
* Jeppe and Dage, Jan 31 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL INPROD
*. Input
      DIMENSION F(NDIM**2,*),S(NDIM**2,*)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(*)
*. Scratch
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
*
*. Zero order energy
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(F,C,VEC1,NDIM,NDIM,0)
      C0FC0 = INPROD(VEC1,C,NDIM)
*
      CALL MATVCB(S,C,VEC1,NDIM,NDIM,0)
      C0SC0 = INPROD(VEC1,C,NDIM)
*
      WRITE(6,*) 'C(0)T F(0) C(0) = ', C0FC0
      WRITE(6,*) 'C(0)T S(0) C(0) = ', C0SC0
      E0 = C0FC0/C0SC0
*. Save diagonal of F(0) - E(0)C(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = F((I-1)*NDIM+I,1)-E0*S((I-1)*NDIM+I,1)
      END DO
      WRITE(6,*) ' Zero order diagonal '
      CALL WRTMAT(VEC3,1,NDIM,1,NDIM)
*. And then start the iterations
      DO IORD = 1, NORD
*
*  =================
*. Energy correction
*  =================
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      - Sum(I=0,N-1)Sum(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
        EN(IORD) = 0.0D0
        DO I = 1, IORD
          CALL MATVCB(F(1,I+1),C(1,IORD-I+1),VEC1,NDIM,NDIM,0)
          EN(IORD) = EN(IORD)+INPROD(C,VEC1,NDIM)
        END DO
        DO I = 0,IORD-1
          DO J = 1, IORD-I
            CALL MATVCB(S(1,J+1),C(1,I+1),VEC1,NDIM,NDIM,0)
            EN(IORD) = EN(IORD)+EN(IORD-I-J)*INPROD(C,VEC1,NDIM)
          END DO
        END DO
        EN(IORD) = EN(IORD)/C0SC0
        WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ==========================
*. Wave function corrections
*  ==========================
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(L=0,N-K)E(N-K-L)S(L)C(K))
        ZERO = 0.0D0
        CALL SETVEC(VEC2,ZERO,NDIM)
        DO K = 0, IORD -1
          CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
        END DO
*
        DO K = 0, IORD -1
          DO L = 0, IORD -K
            CALL MATVCB(S(1,L+1),C(1,K),VEC1,NDIM,NDIM,0)
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
          END DO
        END DO
*. Check overlap with S(0) times zero order state ( should be zero )
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
        CALL MATVCB(S(1,1),C(1,1),VEC1,NDIM,NDIM,0)
        OVLAP = INPROD(VEC1,VEC2,NDIM)
        FACTOR = -OVLAP/C0SC0
        WRITE(6,*) ' Vector before ortho'
        CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
        CALL VECSUM(VEC2,VEC2,C(1,1),ONE,FACTOR,NDIM)
        WRITE(6,*) ' OVLAP = ', OVLAP
        WRITE(6,*) ' Vector before DIAVC2'
        CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
*. Multiply with (F(0)-E(0)S(0))-1
C  DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        ZERO = 0.0D0
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
        CALL COPVEC(VEC1,C(1,IORD+1),NDIM)
        WRITE(6,*) ' Eigenfunction correction ', IORD
        CALL WRTMAT(C(1,IORD+1),1,NDIM,1,NDIM)
      END DO
*
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order             Correction '
      WRITE(6,*) ' ===================================='
      DO IORD = 1, NORD
        WRITE(6,'(1H ,3X,I3,E20.8)')IORD,EN(IORD)
      END DO
*
      RETURN
      END









*













