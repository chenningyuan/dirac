      SUBROUTINE FETCH_CSPEXP_BLK(CSPEXP,IB,NAOFSM,NMOFSM,
     &           IRI,IBU,ISPSM,IAB,CBLK,LROW,LCOL)
Ctimo  We need IAB instead of IAOSM, so that IAOSM can be determined
Ctimo  from IMOSF_SP matrix.
*
* Fetch block of Spinor expansion matrix defined by
* IRI,IBU,ISPSM,IAOSM,IAB,
* IRI = 1 => Real part of expansion
*     = 2 => Imag part of expansion
*
* IBU = 1 => Barred spinor
*     = 2 => unbarred spinor
*
* ISPSM  : Spinor irrep (in dbg sym)
* IAOSM  : Sym    irrep (in pnt sym)
*
* IAB  : Alpha-or beta part of spinor expansion
*
* ============
*. Other input
* ============
*
* CSPEXP : expansion matrix
* IB     : Pointer to start of real/imaginary,barred/unbarred part
* NAOFSM : Number of AO's per sym
* NMOFSM : Number of MO's per sym (in pntgrp sym)
*. Symmetry info from dbgsym
*
* ========
*. Output
* ========
*
* CBLK : BLock of interest ordered as
* CBLK(IASO,ISP), where ISAO is atomic spin orbital and ISP is spinor.
* LROW : Number of rows of output routine
* LCOL : Number of cols of output routine

* A given column of CBLK is organized as  1 : alpha spin atomic orbitals
*                                             of correct symmetry
*                                         2 : beta spin atomic orbitals
*                                             of correct symmetry
* Jeppe + Timo November 1997
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. Input
      DIMENSION CSPEXP(*)
      INTEGER IB(2,2), NMOFSM(*),NAOFSM(*)
*. General input, symmetry info
#include "symm.inc"
#include "input.inc"
*. Output
      DIMENSION CBLK(*)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
*. Offset to this block
      IOFF = IB(IRI,IBU)
Ctimo  determine AO symmetry from given spinor symmetry and spin
      IAOSM = ISPSF_MO(ISPSM,IAB)
      if (IAB.eq.1) IUBDIM = 2
      if (IAB.eq.2) IUBDIM = 1
      IUBSM = ISPSF_MO(ISPSM,IUBDIM)
C?    write(6,*) 'IAOSM is   ',IAOSM
*
      LROW = 0
      LCOL = 0
      DO JSPSM = 1, ISPSM
C?       write(6,*) 'JSPSM is   ',JSPSM
*
         IAOSM_A = ISPSF_MO(JSPSM,1)
         IAOSM_B = ISPSF_MO(JSPSM,2)
C?       write(6,*) 'IAOSM_A is   ',IAOSM_A
C?       write(6,*) 'IAOSM_B is   ',IAOSM_B
         IF(IAOSM_A.NE.0) THEN
           LAO_A = NAOFSM(IAOSM_A)
           LMO_A = NMOFSM(IAOSM_A)
         ELSE
           LAO_A = 0
           LMO_A = 0
         END IF
         IF(IAOSM_B.NE.0) THEN
           LAO_B = NAOFSM(IAOSM_B)
           LMO_B = NMOFSM(IAOSM_B)
         ELSE
           LAO_B = 0
           LMO_B = 0
         END  IF
         LAO_AB = LAO_A+LAO_B
C?       write(6,*) 'LAO_AB is   ',LAO_AB
*
         IF(IBU.EQ.1) THEN
            LSP = LMO_B
         ELSE
            LSP = LMO_A
         END IF
*
         LAO = NAOFSM(IAOSM)
c        LSP = NMOFSM(IUBSM)
*
         IF(JSPSM.NE.ISPSM) THEN
           IOFF = IOFF + LSP*LAO_AB
         ELSE
           IF(IAB.EQ.1.AND.IAOSM_A.EQ.IAOSM) THEN
             DO IAO = 1, LSP
               IOFFE = IOFF+(IAO-1)*LAO_AB
               IOFFP = 1   +(IAO-1)*LAO_A
               CALL COPVEC(CSPEXP(IOFFE),CBLK(IOFFP),LAO_A)
             END DO
           ELSE IF (IAB.EQ.2.AND.IAOSM_B.EQ.IAOSM) THEN
             DO IAO = 1, LSP
               IOFFE = IOFF + (IAO-1)*LAO_AB + LAO_A
               IOFFP = 1    + (IAO-1)*LAO_B
               CALL COPVEC(CSPEXP(IOFFE),CBLK(IOFFP),LAO_B)
             END DO
           END IF
*          ^ End if IAB =1,2
         END IF
*          ^ End of test whether ISPSM = JSPSM
         LROW = LAO
         LCOL = LSP
      END DO
*     ^ End of loop over spinor symmetries
*
      IF(NTEST.gt.2) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' FETCH_CSPEXP_BLK'
        WRITE(6,*) ' ================'
*
        WRITE(6,*) ' Specified output block : '
        WRITE(6,'(A,I4,I4,I4,I4)')
     &        ' IBU IRI IAOSM ISPSM',IBU,IRI,IAOSM,ISPSM
        WRITE(6,*)
        IF(LCOL.NE.0.AND.LROW.NE.0) THEN
          WRITE(6,*) ' output block '
          CALL WRTMAT(CBLK,LROW,LCOL,LROW,LCOL)
        ELSE
          WRITE(6,'(A,I4,I4)') ' Trivial vanishing block, NROW,NCOL= ',
     &    NROW,NCOL
        END IF
      END IF
*
      RETURN
      END
