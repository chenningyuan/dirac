      SUBROUTINE H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal Diag
* and a block matrix of the form
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
*
* Note : The diagonal elements in DIAG corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , May 1990

*
*
* =====
* Input
* =====
* DIAG : Diagonal of matrix
* VECIN : Input vector
* NVAR : Dimension of full matrix
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
* ======
* Output
* ======
* VECUT : Output vector (you guessed ?? ), can occupy same space
*         as VECIN or DIAG
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      REAL * 8  INPROD
*
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
C?    write(6,*) ' H0M1TV , NPQDM = ', NPQDM
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
        CALL GATVEC(WORK(KLV1),VECIN,IPNTR,NPQDM)
* X(T)(DIAG+SHIFT)-1 X in subspace, for later subtraction
        CALL GATVEC(WORK(KLV2),DIAG,IPNTR,NPQDM)
        CALL DIAVC3(WORK(KLV2),WORK(KLV1),
     &       WORK(KLV2),SHIFT,NPQDM,X1)
       ELSE
         X1 = 0.0D0
       END IF
*
      CALL DIAVC3(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
*
      IF(NPQDM .NE. 0 ) THEN
C                H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
C    &           X,RHS,S,SCR)
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
C?     write(6,*) ' KLPHP KLPHQ KLQHQ ',KLPHP,KLPHQ,KLQHQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
         X3 = INPROD(WORK(KLV1),WORK(KLV2),NPQDM)
         CALL SCAVEC(VECUT,WORK(KLV2),IPNTR,NPQDM)
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(6,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3


*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TV '
        WRITE(6,*) ' ========================= '
        CALL WRTMAT(VECUT,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
