      SUBROUTINE H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
* Disc version
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  LUOUT       LUDIA        LUIN
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal on file LUDIA
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           approximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements on LUDIA corresponding to
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
* Jeppe Olsen , September 1993
*
*
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0
* LUIN  : File Containing input vector
* LBLK : Defines format of files
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
COLD  REAL * 8  INPROD
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10)
     & write(6,*) ' H0M1TV , NPQDM = ', NPQDM
*
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
*. Obtain subspace components of input vector
C            GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,NTESTG)
        IZERO = 0
        CALL GATVCD(LUIN,LBLK,NPQDM,IPNTR,WORK(KLV1),VEC1,
     &              IZERO)
*. Solve linear equations in subspace
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
      END IF
*
*. Calculate inverse diagonal and scatter results from subspace,
*. Write to file LUOUT
C     DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
C    &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
      CALL DMTVDS(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &            IPNTR,WORK(KLV2),NPQDM,LBLK,XH0PSX)
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TV '
        WRITE(6,*) ' ========================= '
*. Note : works only if result vector is first file on LUOUT
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
        WRITE(6,*) ' Overlap between input and output vector',
     &               XH0PSX
      END IF
*
      RETURN
      END
