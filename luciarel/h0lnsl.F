      SUBROUTINE H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS

*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVER
*
#include "clunit.inc"
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL QENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          CALL TRIPAK(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMAT(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
C               INVMAT(A,B,MATDIM,NDIM)
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================

         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUSC99,-1)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LUSC99)
*. Initial Guess
        CALL REWINE(LUSC96,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LUSC96)
*. Right hand side
        CALL REWINE(LUSC97,-1)
        CALL TODSC(RHS,NPQDM,-1,LUSC97)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM,LUSC96,LUSC97,LUSC98,LUSC99,
     &              SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LUSC96,-1)
        CALL FRMDSC(SCR(KLV1),NPQDM,-1,LUSC96,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      CALL QEXIT('H0LNS')
      RETURN
      END
