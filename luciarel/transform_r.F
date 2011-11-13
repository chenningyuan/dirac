!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE CSPAOSF_FROM_MOAO(CSPAOSF,NBSPFSM,NUSPFSM,IBSPEXP,
     &                             CMOAO,NAOFSM,NMOFSM)
*
* Obtain spinor-ao*spinfunction transformation from
* mo-ao transformation matrix
*
* Jeppe + Timo, November 1997
*
*. The spinor-ao_spinfunction matrices are ordered as
*
* Loop over spin-down/spin-up spinors
*  Loop over real/imaginary part
*   Loop over symmetries of spinors
*    Loop over AO's * alpha functions with correct spinor sym
*    Loop over AO's * beta  functions with correct spinor sym
*   End of loop over symmetries of spinors
*  End of loop over real/imaginary parts
* End of loop over Spin-down/Spin up
*
* Spin down spinors are spinor set1, spin up spinors are spinor set2
* (B as in Barred comes before U as in Unbarred)
*
* There will be 4 spinor-aosf expansions
* 1 : real part, barred spinors
* 2 : Imag part, barred spinors
* 3 : real part, unbarred spinors
* 4 : Imag part, unbarred spinors
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* ======
*. Input
* ======
*
#include "mxpdim.inc"
*. Info on double group symmetry
#include "symm.inc"
#include "input.inc"
      DIMENSION CMOAO(*)
      INTEGER NMOFSM(*),NAOFSM(*)
*. Number of barred and unbarred spinors per symmetry 
*  - actually recomputed
      INTEGER NBSPFSM(*), NUSPFSM(*)
*
* =======
*. Output
* =======
*. Offset to the four expansions
      INTEGER IBSPEXP(2,2)
*. And the spinor expansions
      DIMENSION CSPAOSF(*)
*
*. Internal switch : We can initialize with either real or imaginary
*  spinors, IRI_ACT = 1 => Real      initialization
*           IRI_ACT = 2 => Imaginary initialization
*
* barred spinors initialized  with spin down-spin orbitals,
* unbarred spinors initialized with spin-up spin orbitals
*
      NTESTL = 0
      NTEST = max(IPRINT,NTESTL)
*
      if (NTEST.ge.3) then
         WRITE(6,*) ' NBSPFSM '
         CALL IWRTMA(NBSPFSM,1,NIRR_DG,1,NIRR_DG)
         WRITE(6,*) ' NUSPFSM '
         CALL IWRTMA(NUSPFSM,1,NIRR_DG,1,NIRR_DG)
         WRITE(6,*) ' NAOFSM '
         CALL IWRTMA(NAOFSM,1,NIRR_DG/2,1,NIRR_DG/2)
      end if
*

      IRI_ACT = 1
      IOFFSM = 1
      DO IBU=1,2
       DO IRI = 1, 2
         IBSPEXP(IRI,IBU) = IOFFSM
         DO ISPSM = 1, NIRR_DG
*. Sym of AO's that will combine with alpha to give double group sym ISPSM
           IAOSM_A = ISPSF_MO(ISPSM,1)
C?         write(6,*) 'ISPSF_MO has value   ',ISPSF_MO(ISPSM,1)
*. Sym of AO's that will combine with beta  to give double group sym ISPSM
           IAOSM_B = ISPSF_MO(ISPSM,2)
C?         write(6,*) 'ISPSF_MO has value   ',ISPSF_MO(ISPSM,2)
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
           END IF
C          NUSPFSM(ISPSM) = LMO_A
C          NBSPFSM(ISPSM) = LMO_B
*
           LAO_AB = LAO_A + LAO_B
           if (NTEST.ge.3) then
              write(6,*) ' ISPSM = ', ISPSM
              write(6,*) ' IAOSM_A, IAOSM_B', IAOSM_A, IAOSM_B
              write(6,*) ' LMO_A, LMO_B, LAO_A,LAO_B',
     &                     LMO_A, LMO_B, LAO_A,LAO_B
           end if
*
*. Alpha-part
*
C?      write(6,*) ' alpha part '
           IF(LAO_A.NE.0) THEN
             IF(IBU.EQ.1.OR.IRI.NE.IRI_ACT) THEN
*. No alpha contributions to initial barred spinors,
               DO ISP = 1, LMO_A
                 ZERO = 0.0D0
                 CALL SETVEC(CSPAOSF(IOFFSM+(ISP-1)*LAO_AB),
     &                       ZERO,LAO_A)
               END DO
             ELSE IF (IBU .EQ. 2.AND.IRI.EQ.IRI_ACT) THEN
*. Offset for this MO-AO block in the CMOAO expansion matrix
               IBCMOAO = IOFF_BLOCK_MAT(NAOFSM,NMOFSM,IAOSM_A)
               DO ISP = 1, LMO_A
                 CALL COPVEC(CMOAO(IBCMOAO+(ISP-1)*LAO_A),
     &                CSPAOSF(IOFFSM+(ISP-1)*LAO_AB),LAO_A)
C?              write(6,*) ' spinor, copying from, to '
C?              write(6,*) ISP,IBCMOAO+(ISP-1)*LAO,
C?   &                     IOFFSM+(ISP-1)*LAO_AB
               END DO
             END IF
           END IF
*          ^ End of alpha part
*
*. beta  part
*
C?      write(6,*) ' beta part '
           IF(LAO_B.NE.0) THEN
             IF(IBU.EQ.2.OR.IRI.NE.IRI_ACT) THEN
*. No beta contribution to unbarred spinors
               DO ISP = 1, LMO_B
                 CALL SETVEC(CSPAOSF(IOFFSM+(ISP-1)*LAO_AB+LAO_A),
     &                       ZERO,LAO_B)
               END DO
             ELSE IF (IBU .EQ. 1.AND.IRI.EQ.IRI_ACT) THEN
*. Offset for this MO-AO block in the CMOAO expansion matrix
               IBCMOAO = IOFF_BLOCK_MAT(NAOFSM,NMOFSM,IAOSM_B)
               DO ISP = 1, LMO_B
                 CALL COPVEC(CMOAO(IBCMOAO+(ISP-1)*LAO_B),
     &                CSPAOSF(IOFFSM+(ISP-1)*LAO_AB+LAO_A),LAO_B)
C?              write(6,*) ' spinor, copying from, to '
C?              write(6,*) ISP,IBCMOAO+(ISP-1)*LAO_B,
C?   &                     IOFFSM+(ISP-1)*LAO_AB+LAO_A
               END DO
             END IF
           END IF
*          ^ End of beta part
*. Update pointer to start of symmetry block
           IF(IBU.EQ.1) THEN
             IOFFSM = IOFFSM + LMO_B*LAO_AB
           ELSE
             IOFFSM = IOFFSM + LMO_A*LAO_AB
           END IF
         END DO
*.       ^ End of loop over spinor symmetries
       END DO
*      ^ End of real/imaginary switch
      END DO
*     ^ End of barred/unbarred switch
*
      if (NTEST.ge.2) then
        WRITE(6,*)
        WRITE(6,*) ' Initial spinor expansion '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        CALL WRT_SPEXP(CSPAOSF,IBSPEXP,NAOFSM,NBSPFSM,NUSPFSM)
C            WRT_SPEXP(CSPSAO, ICSPSAO,NAOFSM,NBSPFSM,NUSPFSM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
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
#include "ipoist8.inc"
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
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTERM_BLK(NBLOCKO,IBLOCKO,INDEXO,LISTO,IRIO,
     &           ISPITRA,ISPITRA_RC,LSPINOR,NRORC,
     &           NLISTO,N_SPINR_EXP,
     &           N1BL,I1BL,N2BL,I2BL,N3BL,I3BL,N4BL,I4BL,MXBLK)
*
* The NBLOCKO blocks IBLOCKO are to  be obtained from the
* transformation of a single input symmetry block .
*
* By using the same intermediate blocks for several
* different output blocks, significant savings in the
* total transformation can be realized. The price is
* - additional routines like this one.
*
*
* We use a very simple algorithm to decide whether
* a given block has 1 or 2 components from
* real/imag components : NRORC = 1 => 1 component
*                        NRORC = 2 => 2 components
*
* In this routine the real and complex parts of
* a block are considered as a single block. This
*
* Jeppe Olsen + Timo Fleig, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "input.inc"
*
*. General input
*. List of spinor transformations in each output list
      INTEGER ISPITRA(4,NLISTO)
*. spinor tranformations, real or complex
      INTEGER ISPITRA_RC(NLISTO)
*. Size of each spinor basis
      INTEGER LSPINOR(N_SPINR_EXP)
*. Specific input
      INTEGER IBLOCKO(4,NBLOCKO),INDEXO(4,NBLOCKO)
      INTEGER LISTO(NBLOCKO),IRIO(NBLOCKO)
*. Output ( see explanation below)
      INTEGER I1BL(6,MXBLK), I2BL(6,MXBLK)
      INTEGER I3BL(6,MXBLK), I4BL(6,MXBLK)
*  Should be added to argument list
      integer IXBL(4,1000)
*
*
*..1 : one index transformed blocks
*
*      N1BL : Number of blocks with one index transformed
*      I1BL : Specification of one-index transformed blocks as :
*        I1BL(1,*) : Spinor transformation
*        I1BL(2,*) : Symmetry
*        I1BL(3,*) : real or imaginary part
*        I1BL(4,*) : Transformation of index 3 or 4
*        I1BL(5,*) : Number of two index transformed blocks
*                    this one-index transformed block will
*                    contribute to
*        I1BL(6,*) : First two-index transformed block this
*                    block will contribute to
*.. 2 : Two index transformed blocks
*      N2BL : Number of blocks with two indeces transformed
*      I2BL : Specification of two-index transformed blocks as :
*        I2BL(1,*) : Spinor transformation
*        I2BL(2,*) : Symmetry of spinor transformation
*        I2BL(3,*) : real or imaginary part
*        I2BL(4,*) : Transformation of index 3 or 4
*        I2BL(5,*) : Number of three index transformed blocks
*                    this two-index transformed block will
*                    contribute to
*        I2BL(6,*) : First three-index transformed block this
*                    block will contribute to
*.. 3 : Three index transformed blocks
*      N3BL : Number of blocks with three indeces transformed
*      I3BL : Specification of three-index transformed blocks as :
*        I3BL(1,*) : Spinor transformation
*        I3BL(2,*) : Symmetry of spinor transformation
*        I3BL(3,*) : real or imaginary part
*        I3BL(4,*) : Transformation of index 1 or 2
*        I3BL(5,*) : Number of four index transformed blocks
*                    this three -index transformed block will
*                    contribute to
*        I3BL(6,*) : First three-index transformed block this
*                    block will contribute to
*.. 4 : four  index transformed blocks
*      N3BL : Number of blocks with three indeces transformed
*      I4BL : Specification of three-index transformed blocks as :
*        I4BL(1,*) : Spinor transformation
*        I4BL(2,*) : Symmetry of spinor transformation
*        I4BL(3,*) : real or imaginary part
*        I4BL(4,*) : Transformation of index 1 or 2
*        I4BL(5,*) : Number of four index transformed blocks
*                    this three -index transformed block will
*                    contribute to
*        I4BL(6,*) : ORIGINAL number of this competely transformed
*                    integral
*
* The idea of our transformation is to reuse the intermediate terms.
* The reason for producing a second list of completely transformed
* integrals as in I4BL is that this list groups transformed
* blocks that are obtained from the same three-index transformed
* list together
*
*
* we transform first the last indeces so the
* transformation can be considered as
*
* ( p q ! r s ) => (p q ! K L )  => ( I J ! K L )
*
* Each half transformation can formally be realized in
* two ways, dependent of what index is transformed first
*
* First half transformation :
* ============================
*
* Route 1 : ( p q ! r s ) => ( p q ! K s ) => ( p q ! K L )
* Route 2 : ( p q ! r s ) => ( p q ! r L ) => ( p q ! K L )
* (alpha beta ! gamma delta )
*
* Second half transformation :
* ===================================
*
* Route 1 : ( p q ! K L ) => ( p J ! K L ) => ( I J ! K L )
* Route 2 : ( p q ! K L ) => ( I q ! K L ) => ( I P ! K L )
*
* since the program holds complete lists of ( p q ! * * )
* and ( * * ! K L ) in core at any given time, both ways are
* feasible, and the actual choice is to minimize the number
* of constructions. In general it is advantageous to transform
* to the smallest set first. For this selection we have included
* an array LSPINOR, containing an integer for each spinor
* transformation. The transformation will be made so the
* index with the smallest transformation expansion is transformed first.
* If the expansions are identical, the last index is transformed
* first.
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Welcome to INTERM_BLK '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
      END IF
*
* =======================================================
*. Part one : Blocks with one index transformed
* =======================================================
*
      N1BL = 0
      DO JBLOCKO = 1, NBLOCKO
*. Index three or index four to be transformed first
        JLISTO = LISTO(JBLOCKO)
        J4 = ISPITRA(4,JLISTO)
        J3 = ISPITRA(3,JLISTO)
        IF(LSPINOR(J4).GE.LSPINOR(J3)) THEN
          JFIRST = 4
          JSPITRA = J4
        ELSE
          JFIRST = 3
          JSPITRA = J3
        END IF
        JSYM = IBLOCKO(JFIRST,JBLOCKO)
*. Has first index of this block been obtained before
        I_AM_NEW = 1
        DO J1BL = 1, N1BL
*         I1BL : Specification of one-index transformed blocks as :
*         I1BL(1,*) : Spinor transformation
*         I1BL(2,*) : Symmetry
*         I1BL(3,*) : real or imaginary part
*         I1BL(4,*) : Transformation of index 3 or 4
*         I1BL(5,*) : Number of two index transformed blocks
*                     this one-index transformed block will
*                     contribute to
*         I1BL(6,*) : First two-index transformed block this
*                     block will contribute to
          IF(I1BL(1,J1BL).EQ.JSPITRA.AND.
     &       I1BL(2,J1BL).EQ.JSYM   .AND.
     &       I1BL(4,J1BL).EQ.JFIRST       ) THEN
             I_AM_NEW = 0
             K1BL = J1BL
          END IF
        END DO
*
        IF(I_AM_NEW.EQ.0) THEN
* record the parentage of this block
          IXBL(1,JBLOCKO) = K1BL
        ELSE
*. Enroll a new block
          N1BL = N1BL + 1
          I1BL(1,N1BL) = JSPITRA
          I1BL(2,N1BL) = JSYM
Ctimo          I1BL(3,N1BL) = NRORC
          I1BL(3,N1BL) = IRIO(JBLOCKO)
          I1BL(4,N1BL) = JFIRST
          IXBL(1,JBLOCKO) = N1BL
        END IF
*       /\ End if block is new
      END DO
*     /\ End of loop over output blocks
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*) '****************************************'
        WRITE(6,*) ' Information about trans of first index '
        WRITE(6,*) '****************************************'
        WRITE(6,*)
        WRITE(6,*) ' Number of unique blocks', N1BL
        WRITE(6,*)
        WRITE(6,*) ' Spinor transformation      Symmetry  Index '
        WRITE(6,*) ' ==========================================='
C2011        DO J1BL = 1, N1BL
C2011          WRITE(6,'(1H ,8X,I3,18X,I3,7X,I2)')
C2011     &    I1BL(1,J1BL),I1BL(2,J1BL),I1BL(4,J1BL)
C2011        END DO
      END IF
*
* =======================================================
*. Part two  : Blocks with two indeces transformed
* =======================================================
*
*. Find blocks with given first transformation
      N2BL= 0
      IB2BL = 1
      DO J1BL = 1, N1BL
        N2BLF1 = 0
        I1BL(6,J1BL) = N2BL + 1
C?     write(6,*) ' I1BL(6,J1BL) = ',I1BL(6,J1BL)
*. What index remains to be transformed
        IF(I1BL(4,J1BL).EQ.4) THEN
         JINDEX = 3
        ELSE
         JINDEX = 4
        END IF
        DO JBLOCKO = 1, NBLOCKO
        IF(IXBL(1,JBLOCKO).EQ.J1BL) THEN
C2011          write(6,*) ' Inside for J1BL JBLOCKO',J1BL,JBLOCKO
*. This output block belongs to the J1BL for
*. the first transformed index. Has this type of
* transformation for the second index been observed?
          I_AM_NEW = 1
*
          JLISTO=LISTO(JBLOCKO)
          JSPITRA=ISPITRA(JINDEX,JLISTO)
          JSYM = IBLOCKO(JINDEX,JBLOCKO)
C2011          WRITE(6,*) ' JLISTO, JSPITRA JSYM',
C2011     &    JLISTO,JSPITRA,JSYM
            do J2BL = I1BL(6,J1BL), I1BL(6,J1BL) + N2BLF1 - 1
            IF(I2BL(1,J2BL).EQ.JSPITRA.AND.
     &         I2BL(2,J2BL).EQ.JSYM   .AND.
     &         I2BL(4,J2BL).EQ.JINDEX       ) THEN
               I_AM_NEW = 0
               K2BL = J2BL
            END IF
          END DO
Ctimo     WRITE(6,*) ' I_AM_NEW = ', I_AM_NEW
*
          IF(I_AM_NEW.EQ.0) THEN
* record the parentage of this block
            IXBL(2,JBLOCKO) = K2BL
          ELSE
*. Enroll a new block
            N2BL = N2BL + 1
            N2BLF1 =  N2BLF1 + 1
            I2BL(1,N2BL) = JSPITRA
            I2BL(2,N2BL) = JSYM
Ctimo            I2BL(3,N2BL) = NRORC
            I2BL(3,N2BL) = IRIO(JBLOCKO)
            I2BL(4,N2BL) = JINDEX
            IXBL(2,JBLOCKO) = N2BL
          END IF
*         /\ End if block is new
        END IF
*       /\ End if first index belongs to J1BL
        END DO
*       /\ End of loop over output blocks
        I1BL(5,J1BL) = N2BLF1
      END DO
*     /\ End of loop over blocks with index 1 transformed
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*) '****************************************'
        WRITE(6,*) ' Information about trans of second index '
        WRITE(6,*) '****************************************'
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks ', N2BL
        WRITE(6,*)
C2011        DO J1BL = 1, N1BL
C2011          WRITE(6,*) ' First group =', J1BL
C2011          WRITE(6,*)  ' Number of second index groups ',
C2011     &    I1BL(5,J1BL)
C2011          WRITE(6,*)
C2011          WRITE(6,*) ' Spinor transformation      Symmetry  Index '
C2011          WRITE(6,*) ' ==========================================='
C2011          DO J2BL = I1BL(6,J1BL),I1BL(6,J1BL)+I1BL(5,J1BL)-1
C?          write(6,*) ' J2BL = ',J2BL
C2011            WRITE(6,'(1H ,8X,I3,18X,I3,7X,I2)')
C2011     &      I2BL(1,J2BL),I2BL(2,J2BL),I2BL(4,J2BL)
C2011          END DO
C2011        END DO
      END IF
*
* =======================================================
*. Part three  : Blocks with three indeces transformed
* =======================================================
*
*. Blocks with given first transformation of dirst two indeces
      N3BL= 0
      IB3BL = 1
*
      DO J1BL = 1, N1BL
      DO J2BL = I1BL(6,J1BL),I1BL(6,J1BL)+I1BL(5,J1BL)-1
        N3BLF12 = 0
        I2BL(6,J2BL) = N3BL + 1
C2011        write(6,*) ' Inside for J1BL J2BL',J1BL,J2BL
        DO JBLOCKO = 1, NBLOCKO
C2011          write(6,*) ' JBLOCKO '
C2011          write(6,*) 'IXBL(1,JBLOCKO),IXBL(2,JBLOCKO)',
C2011     &    IXBL(1,JBLOCKO),IXBL(2,JBLOCKO)
        IF(IXBL(1,JBLOCKO).EQ.J1BL.AND.
     &     IXBL(2,JBLOCKO).EQ.J2BL     ) THEN
*. Should this block be transformed first in index 3 or 4 ?
          JLISTO = LISTO(JBLOCKO)
          J2 = ISPITRA(2,JLISTO)
          J1 = ISPITRA(1,JLISTO)
          IF(LSPINOR(J2).GE.LSPINOR(J1)) THEN
            JFIRST = 2
            JSPITRA = J2
          ELSE
            JFIRST = 1
            JSPITRA = J1
          END IF
          JSYM = IBLOCKO(JFIRST,JBLOCKO)
* has transformation for the third index been included
          I_AM_NEW = 1
C JOJULY  DO J3BL = 1, N3BL
          do J3BL = I2BL(6,J2BL), I2BL(6,J2BL) + N3BLF12 - 1
*
            IF(I3BL(1,J3BL).EQ.JSPITRA.AND.
     &         I3BL(2,J3BL).EQ.JSYM   .AND.
     &         I3BL(4,J3BL).EQ.JFIRST      ) THEN
               I_AM_NEW = 0
               K3BL = J3BL
            END IF
          END DO
*
C?      write(6,*) ' I_am_new',I_AM_NEW
          IF(I_AM_NEW.EQ.0) THEN
* record the parentage of this block
            IXBL(3,JBLOCKO) = K3BL
          ELSE
*. Enroll a new block
            N3BL = N3BL + 1
            N3BLF12 =  N3BLF12 + 1
            I3BL(1,N3BL) = JSPITRA
            I3BL(2,N3BL) = JSYM
Ctimo            I3BL(3,N3BL) = NRORC
            I3BL(3,N3BL) = IRIO(JBLOCKO)
            I3BL(4,N3BL) = JFIRST
            IXBL(3,JBLOCKO) = N3BL
          END IF
*
        END IF
        END DO
        I2BL(5,J2BL) = N3BLF12
      END DO
      END DO
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*) '****************************************'
        WRITE(6,*) ' Information about trans of third index '
        WRITE(6,*) '****************************************'
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks ',N3BL
        WRITE(6,*)
C2011        DO J1BL = 1, N1BL
C2011        DO J2BL = I1BL(6,J1BL),I1BL(6,J1BL)+I1BL(5,J1BL)-1
C2011          WRITE(6,*) ' group one and two ', J1BL, J2BL
C2011          WRITE(6,*) ' Number of third index groups ',
C2011     &    I2BL(5,J2BL)
C2011          WRITE(6,*)
C2011          WRITE(6,*) ' Spinor transformation      Symmetry  Index '
C2011          WRITE(6,*) ' ==========================================='
C2011          DO J3BL = I2BL(6,J2BL),I2BL(6,J2BL)+I2BL(5,J2BL)-1
C2011            WRITE(6,'(1H ,8X,I3,18X,I3,7X,I2)')
C2011     &      I3BL(1,J3BL),I3BL(2,J3BL),I3BL(4,J3BL)
C2011          END DO
C2011        END DO
C2011        END DO
      END IF
*
* =======================================================
*. Part four   : Blocks with four indeces transformed
* =======================================================
*
*. Blocks with given transformations  of first three indeces
      N4BL= 0
      IB4BL = 1
*
      DO J1BL = 1, N1BL
      DO J2BL = I1BL(6,J1BL),I1BL(6,J1BL)+I1BL(5,J1BL)-1
      DO J3BL = I2BL(6,J2BL),I2BL(6,J2BL)+I2BL(5,J2BL)-1
        N4BLF123 = 0
        I3BL(6,J3BL) = N4BL + 1
        DO JBLOCKO = 1, NBLOCKO
        IF(IXBL(1,JBLOCKO).EQ.J1BL.AND.
     &     IXBL(2,JBLOCKO).EQ.J2BL.AND.
     &     IXBL(3,JBLOCKO).EQ.J3BL     ) THEN
C2011           write(6,*) 'JBLOCKO      ', JBLOCKO
C2011           write(6,*) 'IXBL(1,JBLOCKO),IXBL(2,JBLOCKO),IXBL(3,JBLOCKO)',
C2011     &                 IXBL(1,JBLOCKO),IXBL(2,JBLOCKO),IXBL(3,JBLOCKO)
*. Should this block be transformed in index 1 or 2 ?
          IF(I3BL(4,J1BL).EQ.2) THEN
           JINDEX = 1
          ELSE
           JINDEX = 2
          END IF
          JSYM = IBLOCKO(JINDEX,JBLOCKO)
          JLISTO = LISTO(JBLOCKO)
          JSPITRA = ISPITRA(JINDEX,JLISTO)
*. This block has to be new, since it is a
* completely transformed block.
          N4BL = N4BL + 1
          N4BLF123 =  N4BLF123 + 1
          I4BL(1,N4BL) = JSPITRA
          I4BL(2,N4BL) = JSYM
Ctimo          I4BL(3,N4BL) = NRORC
          I4BL(3,N4BL) = IRIO(JBLOCKO)
          I4BL(4,N4BL) = JINDEX
          I4BL(6,N4BL) = JBLOCKO
        END IF
        END DO
        I3BL(5,J3BL) = N4BLF123
      END DO
      END DO
      END DO
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*) '****************************************'
        WRITE(6,*) ' Information about trans of fourth index '
        WRITE(6,*) '****************************************'
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks ', N4BL
        WRITE(6,*)
C2011        DO J1BL = 1, N1BL
C2011    DO J2BL = I1BL(6,J1BL),I1BL(6,J1BL)+I1BL(5,J1BL)-1
C2011        DO J3BL = I2BL(6,J2BL),I2BL(6,J2BL)+I2BL(5,J2BL)-1
C2011          WRITE(6,*) ' group 1,2,3 ', J1BL, J2BL, J3BL
C2011          WRITE(6,*) ' Number of fourth  index groups ',
C2011     &    I3BL(5,J3BL)
C2011          WRITE(6,*)
C2011          WRITE(6,*) ' Spinor transformation      Symmetry  Index '
C2011          WRITE(6,*) ' ==========================================='
C2011          DO J4BL = I3BL(6,J3BL),I3BL(6,J3BL)+I3BL(5,J3BL)-1
C2011            WRITE(6,'(1H ,8X,I3,18X,I3,7X,I2)')
C2011     &      I4BL(1,J4BL),I4BL(2,J4BL),I4BL(4,J4BL)
C2011          END DO
C2011        END DO
C2011        END DO
C2011        END DO
      END IF
*
* That's all she writes, folks
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Distribute integrals to partial lists according to
*  chosen GAS spaces
*
*  Order is symmetry - type at this stage
*  Within one block, the order is UU,UB,BU,BB
*  and within these blocks as in memory.
*
      subroutine intgas3(GASINT,
     &                   OINT,NGSOB,NGAS,NIRR_PN,INTLIST,
     &                   ITOTDIM,IDIM,ISMIND,IGSOFF,IGSIND,
     &                   IACSYM,IOUTLST,MXGSDIM,
     &                   ISPSF_MO,IDBG)
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "onel.inc"
#include "input.inc"
*
      dimension GASINT(4*MXGSDIM**2,NGAS,NGAS)
      dimension OINT(ITOTDIM,ITOTDIM),
     &          NGSOB(MXPIRR,MXPNGAS)
      dimension ISMIND(ITOTDIM),IGSIND(ITOTDIM,IACSYM:IACSYM+1)
      dimension IDIM(2),IGSOFF(2*NGAS),ISPSF_MO(MXNDGIRR,2)
      dimension IUBOFF(4,MXPNGAS,MXPNGAS),NFCS(4,MXPNGAS,MXPNGAS)
*
      NTESTL = 0
      NTEST = MAX(IPRINT,NTESTL)
*
      if (NTEST.ge.2) then
         write(6,*) 'This is intgas3.'
         write(6,*) 'Checking NGSOB matrix'
         do I=1,NIRR_PN,1
            write(6,*) (NGSOB(I,J),J=1,NGAS,1)
         end do
      end if
*
*  Calculate input offsets for start of GAS blocks
      JCT = 1
      IDGGASOFF(1,1) = 1
      ICOUNTER = 1
      IUB = 1
      do JGAS=2,NGAS,1
         if (IDBG.ge.4.and.IDBG.le.6) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
           JCT = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
         else if (IDBG.ge.7.and.IDBG.le.8) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IOUTLST,JGAS-1)
           JCT = JCT + NGSOB(IOUTLST,JGAS-1)
         end if
         ICOUNTER = JGAS
      end do
      if (IDBG.ge.4.and.IDBG.le.6) then
        JCT = JCT + NGSOB(IACSYM+IUB-1,ICOUNTER)
      else if (IDBG.ge.7.and.IDBG.le.8) then
        JCT = JCT + NGSOB(IOUTLST,ICOUNTER)
      end if
      IDGGASOFF(1,2) = JCT
      IUB = 2
      do JGAS=2,NGAS,1
         if (IDBG.ge.4.and.IDBG.le.6) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
           JCT = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
         else if (IDBG.ge.7.and.IDBG.le.8) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IOUTLST,JGAS-1)
           JCT = JCT + NGSOB(IOUTLST,JGAS-1)
         end if
      end do
* Check
      if (NTEST.ge.2) then
        write(6,*) 'Checking offset matrix IDGGASOFF(GAS,IUB) ...'
        call iwrtma(IDGGASOFF,NGAS,2,MXPNGAS,2)
      end if
*
*  Determine output offsets with respect to number of functions
*  in UB/GAS block
      IDBGSM = IOUTLST + NIRR_PN
      do IGAS=1,NGAS,1
        do JGAS = 1,NGAS,1
          IUBOFF(1,IGAS,JGAS) = 1
          do IIJUB = 1,4,1
            if (IIJUB.eq.1) then
             IUB = 1
             JUB = 1
            else if (IIJUB.eq.2) then
             IUB = 1
             JUB = 2
            else if (IIJUB.eq.3) then
             IUB = 2
             JUB = 1
            else
             IUB = 2
             JUB = 2
            end if
*
            IPGSM = ISPSF_MO(IDBGSM,IUB)
            JPGSM = ISPSF_MO(IDBGSM,JUB)
            NFCS(IIJUB,IGAS,JGAS) =
     &           NGSOB(IPGSM,IGAS) * NGSOB(JPGSM,JGAS)
*
            if (IIJUB.le.3) then
              IUBOFF(IIJUB+1,IGAS,JGAS) = IUBOFF(IIJUB,IGAS,JGAS)
     &                                    + NFCS(IIJUB,IGAS,JGAS)
            end if
          end do
        end do
      end do
*
      if (NTEST.ge.3) then
        do IGAS=1,NGAS,1
          do JGAS=1,NGAS,1
            do IIJUB=1,4,1
              write(6,*) 'IIJUB is  ',IIJUB
              write(6,'(A,1X,I6,I3,I3,A,1X,I6)')
     &                   'NFCS(',IIJUB,IGAS,JGAS,') = ',
     &                    NFCS(IIJUB,IGAS,JGAS)
              write(6,'(A,1X,I6,I3,I3,A,1X,I6)')
     &                   'IUBOFF(',IIJUB,IGAS,JGAS,') = ',
     &                    IUBOFF(IIJUB,IGAS,JGAS)
            end do
          end do
        end do
      end if
*
*  Distribute integrals
*
      call evenodd(IEOSYM,IOUTLST)
      IGSCT = 0
      IUBCT = 0
      if (IEOSYM.eq.1) then
*  Unbarred up left
        do JGAS=1,NGAS,1
           do IGAS=1,NGAS,1
              do IUB=1,2,1
                 do JUB=1,2,1
                    ISTART = IDGGASOFF(IGAS,IUB)
                    IEND = IDGGASOFF(IGAS,IUB) +
     &                      NGSOB(IOUTLST+IUB-1,IGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      IEND = IDGGASOFF(IGAS,IUB) +
     &                       NGSOB(IOUTLST,IGAS) - 1
                    end if
                    JSTART = IDGGASOFF(JGAS,JUB)
                    JEND = IDGGASOFF(JGAS,JUB) +
     &                      NGSOB(IOUTLST+JUB-1,JGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      JEND = IDGGASOFF(JGAS,JUB) +
     &                       NGSOB(IOUTLST,JGAS) - 1
                    end if
                    if (NTEST.ge.5) then
                      write(6,*) 'ISTART =  ',ISTART
                      write(6,*) 'IEND =  ',IEND
                      write(6,*) 'JSTART =  ',JSTART
                      write(6,*) 'JEND =  ',JEND
                    end if
*
                    if (IUB.eq.1.and.JUB.eq.1) NIJUB = 1
                    if (IUB.eq.1.and.JUB.eq.2) NIJUB = 2
                    if (IUB.eq.2.and.JUB.eq.1) NIJUB = 3
                    if (IUB.eq.2.and.JUB.eq.2) NIJUB = 4
                    IOUTOFF = IUBOFF(NIJUB,IGAS,JGAS)
                    if (NTEST.ge.5) then
                      write(6,*) 'Output offset is ',IOUTOFF
                    end if
*
                    do J=JSTART,JEND,1
                      do I=ISTART,IEND,1
                        GASINT(IOUTOFF+IUBCT,IGAS,JGAS)
     &                       = OINT(I,J)
                        if (NTEST.ge.5) then
                          write(6,*) 'Input integral ',OINT(I,J)
                          write(6,*) 'Output position ',IOUTOFF+IUBCT
                        end if
                        IUBCT = IUBCT + 1
                        IGSCT = IGSCT + 1
                      end do
                    end do
                    IPUBPGS(NIJUB,IGAS,JGAS,INTLIST)
     &                     = IUBCT
                    IUBCT = 0
                 end do
              end do
              if (NTEST.ge.5) then
                 do I=1,IGSCT,1
                    write(6,*) I,GASINT(I,IGAS,JGAS)
                 end do
                 write(6,*) '... end of block ...'
              end if
              IGSCT = 0
           end do
        end do
      else if (IEOSYM.eq.2) then
*  Barred up left: reorder!
        do JGAS=1,NGAS,1
           do IGAS=1,NGAS,1
              do IBU=2,1,-1
                 do JBU=2,1,-1
                    ISTART = IDGGASOFF(IGAS,IBU)
                    IEND = IDGGASOFF(IGAS,IBU) +
     &                      NGSOB(IACSYM+IBU-1,IGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      IEND = IDGGASOFF(IGAS,IBU) +
     &                       NGSOB(IOUTLST,IGAS) - 1
                    end if
                    JSTART = IDGGASOFF(JGAS,JBU)
                    JEND = IDGGASOFF(JGAS,JBU) +
     &                      NGSOB(IACSYM+JBU-1,JGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      JEND = IDGGASOFF(JGAS,JBU) +
     &                       NGSOB(IOUTLST,JGAS) - 1
                    end if
                    if (NTEST.ge.3) then
                      write(6,*) 'ISTART =  ',ISTART
                      write(6,*) 'IEND =  ',IEND
                      write(6,*) 'JSTART =  ',JSTART
                      write(6,*) 'JEND =  ',JEND
                    end if
*
                    if (IBU.eq.1.and.JBU.eq.1) NIJUB = 4
                    if (IBU.eq.1.and.JBU.eq.2) NIJUB = 3
                    if (IBU.eq.2.and.JBU.eq.1) NIJUB = 2
                    if (IBU.eq.2.and.JBU.eq.2) NIJUB = 1
                    IOUTOFF = IUBOFF(NIJUB,IGAS,JGAS)
                    if (NTEST.ge.3) then
                      write(6,*) 'Output offset is ',IOUTOFF
                    end if
*
                    do J=JSTART,JEND,1
                       do I=ISTART,IEND,1
                          GASINT(IOUTOFF+IUBCT,IGAS,JGAS)
     &                         = OINT(I,J)
                          if (NTEST.ge.5) then
                            write(6,*) 'Input integral ',OINT(I,J)
                            write(6,*) 'Output position ',IOUTOFF+IUBCT
                          end if
                          IUBCT = IUBCT + 1
                          IGSCT = IGSCT + 1
                       end do
                    end do
                    IPUBPGS(NIJUB,IGAS,JGAS,INTLIST)
     &                      = IUBCT
                    IUBCT = 0
                 end do
              end do
              if (NTEST.ge.2) then
                 do I=1,IGSCT,1
                    write(6,*) I,GASINT(I,IGAS,JGAS)
                 end do
                 write(6,*) '... end of block ...'
              end if
              IGSCT = 0
           end do
        end do
      end if
*
      if (NTEST.ge.2) then
         do IGAS=1,NGAS,1
            do JGAS=1,NGAS,1
               write(6,*) 'IPUBPGS for IGAS',IGAS,'and JGAS',JGAS
               call iwrtma(IPUBPGS(1,IGAS,JGAS,INTLIST),4,1,4,1)
            end do
         end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Arrange AO integral matrix for transformation.
*  Read integrals are on AOTRISM in packed triangle form,
*  all symmetries subsequently. This routine picks the
*  desired blocks and unpacks them to AOINT matrix.
*  This matrix finally is transformed by oneltra.
*
*!!!! CAUTION !!!
*!!!! CAUTION !!!
*
*  If spiorb = true, then HELPMAT2 contains the input
*  integrals. Otherwise, i.e. for spin-free integrals,
*  the input integrals are on AOTRISM
*  (This is terrible, but pragmatic ...)
*
*!!!! CAUTION !!!
*!!!! CAUTION !!!
*
      subroutine intsel(AOINT,AOTRISM,HELPMAT2,
     &                  PROD,IPRODDIM,
     &                  IPGIDIM,IDGIDIM,ICURRDIM,IDIMGRP,NIRR_PN,
     &                  ITRIA,ITRIASUM,ITRIADB,ITRIADBSUM,
     &                  SPIORB,IDBG,NTRIAEL,ENVIRO,COORD)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "input.inc"
*
      logical SPIORB
      character*72 ENVIRO
      character*1 COORD
*
      dimension IPGIDIM(8),IDGIDIM(2),ITRIA(8),ITRIADB(4)
      dimension AOINT(ICURRDIM,ICURRDIM)
      dimension AOTRISM(ITRIASUM),
     &          HELPMAT2(NTRIAEL),PROD(IPRODDIM)
*
*  Spin-free or spin-orbit integrals:
*  The whole stuff only has to be done if we have symmetry
*  blocks of integrals which have to be picked and expanded.
*  In the SO case, this is obsolete.
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      ZERO = 0.0D0
*  Erase AOINT (just to make sure)
      do ICOL=1,IDGIDIM(IDIMGRP),1
        call setvec(AOINT(1,ICOL),ZERO,IDGIDIM(IDIMGRP))
      end do
*
*---------------------------------------------------------------
      if (SPIORB) goto 500
*
*  Processing spin-orbit integrals (HELPMAT2 contains input ints)
*---------------------------------------------------------------
*
      if (NTEST.gt.2) then
      write(6,*) 'Checking matrix AOTRISM ...'
         do I=1,ITRIASUM,1
            write(6,'(1G18.12)') AOTRISM(I)
         end do
      end if
*
*  Erase matrix HELPMAT2 (it has been used in rdcoef)
      do I=1,NTRIAEL,1
         HELPMAT2(I) = ZERO
      end do
*
*  Determine offset for relevant symmetry block(s)
      ISTART = 1
      if (IDIMGRP.eq.2) then
         do J=1,2,1
            ISTART = ISTART + ((IPGIDIM(J)**2 + IPGIDIM(J)) /2)
         end do
         if (IDBG.eq.7)
     &       ISTART = 1 + ((IPGIDIM(1)**2 + IPGIDIM(1)) /2)
      end if
C?    write(6,*) 'ISTART is    ',ISTART
*
*  Expand to full triangle matrix
      ISYM = 2*IDIMGRP - 1
      if (IDBG.eq.7) ISYM = IDIMGRP
      IEND = ISTART + ((IPGIDIM(ISYM)**2 + IPGIDIM(ISYM)) /2) - 1
C?    write(6,*) 'IEND is    ',IEND
      if (IDIMGRP.eq.2) then
         IENDEND = ITRIASUM
      else if (IDIMGRP.eq.1) then
         IENDEND = 0
         do J=1,2,1
            IENDEND = IENDEND + ((IPGIDIM(J)**2 + IPGIDIM(J)) /2)
         end do
         if (IDBG.eq.7)
     &       IENDEND = ((IPGIDIM(1)**2 + IPGIDIM(1)) /2)
      end if
      NSKIP = IPGIDIM(ISYM)
*
      if (IDBG.ge.4.and.IDBG.le.6) then
        J = 1
        do I=ISTART,IEND,1
          HELPMAT2(J) = AOTRISM(I)
          J = J + 1
        end do
*
        INSYM = IPGIDIM(ISYM+1)
        KOUT = J
        INP = IEND + 1
        do ILINE=1,INSYM,1
          KOUT = KOUT + NSKIP
          do I=INP,INP+ILINE-1,1
            HELPMAT2(KOUT) = AOTRISM(I)
            if (NTEST.gt.2) then
              write(6,*) 'HELPMAT2(KOUT) is    ',
     &                    HELPMAT2(KOUT),KOUT
              write(6,*) 'AOTRISM(I) is    ',AOTRISM(I),I
            end if
            KOUT = KOUT + 1
          end do
          INP = INP + ILINE
        end do
      else if (IDBG.ge.7.and.IDBG.le.8) then
        if (IDIMGRP.eq.1) then
          J = 1
          do I=ISTART,IEND,1
            HELPMAT2(J) = AOTRISM(I)
            J = J + 1
          end do
        else if (IDIMGRP.eq.2) then
          INSYM = IPGIDIM(2)
          KOUT = 1
          do I=ISTART,IEND,1
            HELPMAT2(KOUT) = AOTRISM(I)
            if (NTEST.gt.2) then
              write(6,*) 'HELPMAT2(KOUT) is    ',
     &                    HELPMAT2(KOUT),KOUT
              write(6,*) 'AOTRISM(I) is    ',AOTRISM(I),I
            end if
            KOUT = KOUT + 1
          end do
        end if
      end if
*
 100  continue
*
*  Check matrix HELPMAT2:
      if (NTEST.gt.2) then
      write(6,*) 'Checking matrix HELPMAT2:'
         do I=1,NTRIAEL,1
            write(6,'(1G18.12)') HELPMAT2(I)
         end do
      end if
*
*------------------------------------------
 500  continue
*------------------------------------------
*
*  Unfold to full matrix for transformation
      if (SPIORB) then
*
        if (NTEST.ge.5) then
          write(6,*) 'Checking input SO integrals in intsel:'
          call wrtmat(HELPMAT2,NTRIAEL,1,NTRIAEL,1)
        end if
*
*  Erase PROD (just to make sure)
        call setvec(PROD,ZERO,IPRODDIM)
*
        if (IDIMGRP.eq.2) then
           if (ENVIRO(1:6).eq.'MOLCAS') then
*  Pick second dimension group:
*  Check for pathological cases ...
*  (See LUCIAREL notes for derivation of these formulae!)
             ILENREQ = IPGIDIM(3)*IPGIDIM(4)
             if (ILENREQ.gt.IPRODDIM) then
                write(6,*) 'Problem in intsel (SO)'
                write(6,*) 'Required length:  ',ILENREQ
                write(6,*) 'Length usable: ',IPRODDIM
                Call Abend2( 'Quitting.' )
             end if
*
             if (IDBG.ge.4.and.IDBG.le.6) then
               J = 1
               INP = ITRIADB(1) + IDGIDIM(1) + 1
               do ISYM=1,IDGIDIM(2),1
                do I=INP,INP+ISYM-1,1
                     PROD(J) = HELPMAT2(I)
                     J = J + 1
                  end do
                  INP = INP + ISYM + IPGIDIM(1) + IPGIDIM(2)
               end do
               IGLOFAC = +1
               call unfold(AOINT(1,1),PROD(1),IDGIDIM(IDIMGRP),
     &                     -1,ITRIADB(IDIMGRP),IGLOFAC)
             else if (IDBG.eq.7) then
               J = 1
               INP = ITRIA(1) + IPGIDIM(1) + 1
               do ISYM=1,IPGIDIM(2),1
                 do I=INP,INP+ISYM-1,1
                   PROD(J) = HELPMAT2(I)
                   J = J + 1
                 end do
                 INP = INP + ISYM + IPGIDIM(1)
               end do
               IGLOFAC = +1
               call unfold(AOINT(1,1),PROD(1),IDGIDIM(2),
     &                     -1,ITRIA(2),IGLOFAC)
               if (COORD.eq.'Z') then
                 do ICOL=1,IPGIDIM(2),1
                   call copvec(AOINT(1,ICOL),
     &               AOINT(IPGIDIM(2)+1,IPGIDIM(2)+ICOL),IPGIDIM(2))
                   call setvec(AOINT(1,IPGIDIM(2)+ICOL),
     &                         ZERO,IPGIDIM(2))
                   call setvec(AOINT(IPGIDIM(2)+1,ICOL),
     &                         ZERO,IPGIDIM(2))
                 end do
               else if (COORD.eq.'X'.or.COORD.eq.'Y') then
                 do ICOL=1,IPGIDIM(2),1
                   call copvec(AOINT(1,ICOL),
     &               AOINT(1,IPGIDIM(2)+ICOL),IPGIDIM(2))
                   call copvec(AOINT(1,ICOL),
     &               AOINT(IPGIDIM(2)+1,ICOL),IPGIDIM(2))
                   call setvec(AOINT(1,ICOL),ZERO,IPGIDIM(2))
                   call setvec(AOINT(IPGIDIM(2)+1,IPGIDIM(2)+ICOL),
     &                         ZERO,IPGIDIM(2))
                 end do
               end if
             end if
*
           else if (ENVIRO(1:6).eq.'MOLECU') then
*  Pick blocks from MOLECULE ordering (Ag,Au,Bu,Bg)
             IAUBUDIM = IPGIDIM(3) + IPGIDIM(4)
             IAUSTA = ITRIA(1) + IPGIDIM(1) + 1
             do I=1,IAUBUDIM,1
               do II=1,I,1
                 PROD(J) = HELPMAT2(IAUSTA+II-1)
                 J = J + 1
               end do
               IAUSTA = IAUSTA + IPGIDIM(1) + I
             end do
           end if
*          ^ end ENVIRO
        else if (IDIMGRP.eq.1) then
           if (ENVIRO(1:6).eq.'MOLCAS') then
             if (IDBG.ge.4.and.IDBG.le.6) then
               IGLOFAC = +1
               call unfold(AOINT(1,1),HELPMAT2(1),IDGIDIM(IDIMGRP),
     &                     -1,ITRIADB(IDIMGRP),IGLOFAC)
             else if (IDBG.ge.7.and.IDBG.le.8) then
               IGLOFAC = +1
               if (COORD.eq.'Z') then
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),1,1)
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),
     &                          IPGIDIM(1)+1,IPGIDIM(1)+1)
               else if (COORD.eq.'X'.or.COORD.eq.'Y') then
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),
     &                          IPGIDIM(1)+1,1)
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),
     &                          1,IPGIDIM(1)+1)
               end if
             end if
*
           else if (ENVIRO(1:6).eq.'MOLECU') then
             I3TRIADIM = IPGIDIM(1) + IPGIDIM(3) + IPGIDIM(4)
             I3TRIA = (I3TRIADIM**2 + I3TRIADIM) / 2
             IBGSTA = I3TRIA + 1
             IBGDIM = IPGIDIM(2)
             IAGDIM = IPGIDIM(1)
             IAGEND = ITRIA(1)
             J = 1
             do I=1,IBGDIM,1
               do K=1,IAGDIM,1
                 PROD(J) = HELPMAT2(IBGSTA+K-1)
                 J = J + 1
               end do
               do II=1,I,1
                 PROD(J) = HELPMAT2(IBGSTA+I3TRIADIM+II-1)
                 J = J + 1
               end do
             end do
             K = 1
             do JJ=1,J-1,1
               HELPMAT2(IAGEND+K) = PROD(JJ)
               K = K + 1
             end do
             IGLOFAC = +1
             call unfold(AOINT(1,1),HELPMAT2(1),IDGIDIM(IDIMGRP),
     &                   -1,ITRIADB(IDIMGRP),IGLOFAC)
           end if
*          ^ end ENVIRO
        else
           write(6,*) 'Problem in intsel.'
           write(6,*) 'IDIMGRP is  ',IDIMGRP
           Call Abend2( 'Quitting.' )
        end if
*       ^ end IDIMGRP
*
*  IGLOFAC is the global factor for l integrals, as pure r x nabla
*  is calculated by genmfint/amfi. The l integrals should be
*  h/i = -ih, so assume imaginary integrals and multiply all
*  SO ints with -1: no, we're not doing that right now ...
*
      else
         IGLOFAC = 1
         call unfold(AOINT(1,1),HELPMAT2(1),IDGIDIM(IDIMGRP),1,
     &               ITRIA(IDIMGRP),IGLOFAC)
         if (IDBG.ge.7.and.IDBG.le.8) then
           do ICOL=1,IPGIDIM(IDIMGRP),1
             call copvec(AOINT(1,ICOL),
     &              AOINT(IPGIDIM(IDIMGRP)+1,IPGIDIM(IDIMGRP)+ICOL),
     &              IPGIDIM(IDIMGRP))
           end do
         end if
      end if
*     ^ end SPIORB
*
      if (NTEST.ge.3) then
         write(6,*) 'Checking AOINT in intsel:'
         if (SPIORB) then
            write(6,*) 'These are spin-orbit integrals'
         else
            write(6,*) 'These are spin-free integrals'
         end if
         call wrtmat(AOINT,IDGIDIM(IDIMGRP),IDGIDIM(IDIMGRP),
     &               IDGIDIM(IDIMGRP),IDGIDIM(IDIMGRP))
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTTRA(I1234,I12,I34,ICCSM,IOPSMPN,IOPSMDG,NAOFSM,
     &                  NSMSP,NSMAO,MULTSP,MULTAO,
     &                  NLISTO,I12O,I34O,I1234O,ICCSMO,
     &                  ISPITRA,ISPITRA_RC,IAOSFSP,ISPSFAO,LSPINOR )
*
* General integral transformation for CI with general
* spinors
*
* On input we must have nformation about
*     1  The list of input ao. integrals
*     2  sp-ao transformation matrix( matrices)
*     3  The required list of transformed integrals
*     4  symmetry information,
*
* +++
*  1
* +++
*
* The list a.o. integrals are specified by
* a) The number of a.o per sym  NAOFSM
* b) The symmetry of the operator in the pointgroup ISMPN
* c) The permutational symmetry of the 4 indeces as specified by
*    I1234,I12,I34
*
* I1234 : Symmetry between particle 1 and 2
*
*     1 =>   (ij!kl) = (kl!ij)
*    -1 =>   (ij!kl) =-(kl!ij)
*     0 =>  No symmetry between particle 1 and 2
*
* I12   : Symmetry between first two indeces (indeces of particle 1 )
*
*    1 => (ij!kl) = (ji!kl)
*   -1 => (ij!kl) =-(ji!kl)
*    0 => No symmetry between the first two indeces
*
* I34   : Symmetry between last two indeces (indeces of particle 2)
*
*    1 => (ij!kl) = (ij!lk)
*   -1 => (ij!kl) =-(ij!lk)
*    0 => No symmetry between the last two indeces
*
*    We will not assume anything about the specific order of
*    the a.o integral list or the block sizes ...
*    We will instead invoke the a.o integrals through subroutine
*    calls where these details for specific implementations are
*    hidden
*
* +++
*  2
* +++
*
* There can be more than one sp-ao transformation matrix
* ( barred and unbarred spinors for example), and each
* transformation matrix is in general complex
*
*
* So we have as input
*  NSPAO : Number of mo-ao transformations
*  NSPFSMG(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, general expansion. The symmetry refers to
*  general spinor symmetry, say double groups.
*  NSPFSMR(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, restricted expansion,The symmetry refers to
*  general spinor symmetry, say double groups.
*  CSPAO : List of all SP-AO transformation matrices
*  ICSPAO : ICSPAO(IRI,IEXP)  Offset for Real (IRI = 1)
*  or imaginary part (IRI = 2 ) of SPAO expansion IEXP
*
* The  real or imaginary part of a given sp-ao expansion is
* ordered as
*
* Loop over symmetries
*  Loop over spinors of given symmetry
*   Loop over alpha functions * spatial functions of correct symmetry
*   End alpha functions * spatial functions of correct symmetry
*   Loop over beta  functions * spatial functions of correct symmetry
*   End beta  functions * spatial functions of correct symmetry
*  End of loop over spinors of given symmetry
* End of loop over symmetries
*
* The real/imaginary part of a given spao expansion  contains thus
* in each column the expansion of a given spinor. In the given column
* the are first alpha spin orbitals and the beta spin orbitals.
* The symmetries of the ao orbitals that combine with a given
* spin function to give a certain double group symmetry is
* given by the symmetry tables ( see these )
*
*
*
* +++
*  3
* +++
*
* I guess the following is obsolete, since the transformation
* to restricted sets is realized in a different way.
* One defines restricted spinor transformations,
* and includes in the input setup, the relevant types
* of reduced inegrals lsists
*
* N_GEN_IND
*
*    In general one is not interested in a complete list of integrals
*    rather, one or several indeces are restricted to a subset
*    This is specified by the N_GEN_IND parameter. Using the
*    notation that g stands for a general index, whereas
*    r stands for a restricted index, we have the following
*
*    0 => (rr!rr) integrals ( for a CI in restricted orbital space)
*    1 => (rr!rg) integral list ( for fisrst order casscf)
*    2 => (rr!gg)+(rg!rg) integral list ( for second order mcscf
*          and linear response, CASPT2)
*    3 => for quadratic response (not extremely relevant)
*    4 => for general CI in all orbitals
*
* Well, to make life (atleast programming) really complicated
* we are in general interested in more than one output integral list
* so another input parameter is
*
* NLISTO : Number of output integral lists
*
* The  output list are defined by
*
* I1234O(ILIST), I12O(ILIST), I34O(ILIST),and ICCSMO
* giving permutational
* symmetries of output list ILIST. The newcomer is ICCSMO
* that indicated whether a list ishould be restricted
* by complex conjugation restrictions, (12!34) = (21!43)*.
* This will define the index restrictions on the output integrals
*
*
* ISPITRA(INDEX,ILIST) : Index number INDEX(1,2,3,4) in output
* list ILIST refers to spao expansion ISPITRA (INDEX,ILIST)
*
* ISMDGG : Symmetry of the operator in the double group
* representation
*
* +++
*  4
* +++
*
* The information  about symmetry is
*
* NSMSP : Number of symmetries of spinors ( number of iireps in
*         double group )
* NSMAO : Number of symmetries of atomic orbitals ( Number of
*         irreps of pointgroup)
* ISPSFAO(ISPSM,IAB) : Spinor symmetry X spinfunction => ao symmetry
*                      The point group symmetry to be multiplied
*                      with alpha (IAB=1) or beta (IAB = 2)  to
*                      give spinor symmetry ISPSM
* IAOSFSP(IAOSM,IAB) : ao symmetry X spin function => spinor symmetry
*                      The spinor  symmetry to be multiplied
*                      with alpha (IAB=1) or beta (IAB = 2)  to
*                      give point group symmetry ISPSM
*
* Timo Fleig + Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. input
*.. Symmetry information
      INTEGER ISPSFAO(NSMSP,2), IAOSFSP(NSMAO,2)
*.. Information about output integral lists
      INTEGER ISPITRA(4,NLISTO)
      INTEGER I12O(NLISTO),I34O(NLISTO),I1234O(NLISTO)
*  NSPAO : Number of mo-ao transformations
*  NSPFSMG(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, general expansion. The symmetry refers to
*  general spinor symmetry, say double groups.
*  NSPFSMR(ISYM,IEXP) : Number of spinors of symmetry ISYM in
*  expansion IEXP, restricted expansion,The symmetry refers to
*  general spinor symmetry, say double groups.
*  CSPAO : List of all SP-AO transformation matrices
*  ICSPAO : ICSPAO(IRI,IEXP)  Offset for Real (IRI = 1)
*  or imaginary part (IRI = 2 ) of SPAO expansion IEXP
*
C      INTTRA(I1234,I12,I34,ICCSM,IOPSMPN,NAOFSM,
C     &                  NSMSP,NSMAO,MULTSP,MULTAO,
C     &                  NLISTO,I12O,I34O,I1234O,ICCSMO,
C     &                  ISPITRA,ISPITRA_RC,IAOSFSP,ISPSFAO         )
*
      WRITE(6,*) ' ISPITRA'
      CALL IWRTMA(ISPITRA,4,NLISTO,4,NLISTO)
*
      CALL SPITRA2(I1234,I12,I34,NSMAO,MULTAO,NSMSP,MULTSP,
     &             IOPSMPN,IOPSMDG,NLISTO,I1234O,I12O,I34O,ICCSMO,
C    &             ISPITRA,ISPITRA_RC,IAOSFSP,ISPSFAO,LSPINOR)
     &             ISPITRA)
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Control section for integral transformation
*
      subroutine itrctl(INTSPACE,ITRA_TYPE,INTITER,
     &                  ENVIRO,ISPINFREE,IFOCKEX,NORB_FULL,
     &                  NL2D,IUB2D,ISM2D,ISYM_T,SIGDEN_ROUTE,
     &                  F1,F2,T_BUFF,T_BUFF_D,IT_TTPL,IBOSYM,
     &                  MCSCF,CIRUN)
*
*---------------------------------------------------
*  Written by Timo Fleig (including all subroutines)
*---------------------------------------------------
*
*  Implementation of MPI "shared" memory version 
*  (use of T_BUFF) by Stefan Knecht  
*
*
      use memory_allocator
#include "implicit.inc"
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "krciprop.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "symm.inc"
#include "glbbas.inc"
#include "input.inc"
#include "clunit.inc"
#include "oper.inc"
#include "ctcc.inc"
#include "integrals_off.inc"
#if defined (VAR_MPI2)
#include "krmc_shmem.h"
#include "mxdim_mpi2.h"
      POINTER (MY_CHECK_F_PTR, ICHECK_FILE(IDUMMY_BUFF_SZ))
#endif
*
#include "dgroup.h"
#include "dcborb.h"
*
      character*72 ENVIRO
      character*6 SIGDEN_ROUTE
      character*6 CIRUN
      character*1 COORD
      LOGICAL MCSCF, EX, ALLOC_INT, ALLOC_INT_SCR, fileopen, ex_mc
      dimension INTEG(1000), T_BUFF(*), T_BUFF_D(*), IT_TTPL(*)
      dimension NORB_FULL(MXPIRR)
      integer*8 IDUMMY_LENGTH
      integer, allocatable :: ibeig(:)
*
* ITRA_TYPE       Defines the type of transformation
*---------------------------------------------------
* ITRA_TYPE = 1 : Normal nonrelativistic transformation
* ITRA_TYPE = 2 : 2/4-component relativistic transformation for CI,
*                 all nonredundant integral lists
* ITRA_TYPE = 3 : 2/4-component relativistic transformation for CI,
*                 Excluding integrals with an odd number of bars
* ITRA_TYPE = 4 : 2/4-component relativistic transformation for MCSCF,
*                 or density matrix spinor optimization
*                 all nonredundant integral lists
* ITRA_TYPE = 5 : 2/4-component relativistic transformation for MCSCF,
*                 or density matrix spinor optimization
*                 Excluding integrals with an odd number of bars
*
C
      if (ENVIRO(1:6).eq.'MOLCAS') then
        FNONEAO = 'ONEINT'
        FNONEMO = 'TRAONE'
        INTIMP = 5
      else if (ENVIRO(1:6).eq.'DIRAC ') then
        DIRACONE = 'MRCONEE'
        DIRACTWO = 'MDCINT'
        INTIMP = 6
        if (MCSCF) INTIMP = 7
      end if
*
      CALL QENTER('ITRCTRL')
*
      NTESTL = 00
      NTEST = max(IPRINT,NTESTL)
      NTEST = 000
*
      IF( NTEST .ge. 1 ) THEN
      WRITE(LUWRT,*) ' ************************************************'
      WRITE(LUWRT,*) ' *                                              *'
      WRITE(LUWRT,*) ' *  Transformer man                             *'
      WRITE(LUWRT,*) ' *  You run the show                            *'
      WRITE(LUWRT,*) ' *  Remote control                              *'
      WRITE(LUWRT,*) ' *  Direct the action with the push of a button *'
      WRITE(LUWRT,*) ' *                                              *'
      WRITE(LUWRT,*) ' *                               Neil Young     *'
      WRITE(LUWRT,*) ' *                                              *'
      WRITE(LUWRT,*) ' ************************************************'
      END IF
*
*  CI or MCSCF transformation?
*  Note: We will do the one-particle transformation always with
*        the full orbital space, it doesn't kill us. Afterwards,
*        a routine will select a matrix of integrals over the
*        active spaces for the GAS routines.
*
      NTRIAEL = ((SYMDIM(1)+SYMDIM(2))**2
     &          +(SYMDIM(1)+SYMDIM(2))) /2
*
      NTRSM_FULL = 0
      NDIM_TOT = 0
      do ISYM=1,NIRR_PN,1
        NDIM_TOT = NDIM_TOT + SYDI(ISYM)
      end do
*
      call isetvc(SYMDIM_F,0,NIRR_DG)
      if (DOUGRP.lt.7.and.DOUGRP.ge.4) then
        I = 1
        do IRR=1,NIRR_PN,1
          if (IRR.eq.3.or.IRR.eq.5.or.IRR.eq.7) I = I + 1
          SYMDIM_F(I) = SYMDIM_F(I) + NORB_FULL(IRR)
        end do
      else if (DOUGRP.eq.7) then
        SYMDIM_F(1) = 2*NORB_FULL(1)
        SYMDIM_F(2) = 2*NORB_FULL(2)
      else if (DOUGRP.eq.8) then
        SYMDIM_F(1) = 2*NORB_FULL(1)
      end if 
*
      if (IFOCKEX.eq.1.and.ENVIRO(1:6).eq.'DIRAC ') then
        write(LUWRT,*) 'Problem in itrctl.'
        write(LUWRT,*) 'IFOCKEX is  ',IFOCKEX
        write(LUWRT,*) 'and environment is ',ENVIRO(1:6)
        write(LUWRT,*) 'This is not allowed. Omit FOCKEX variable.'
        Call Abend2( 'Quitting.' )
      end if
*
*
*  1) Static memory
*
*  Transformed integrals, 32 lists
      call findhigh(NGSOBPI,MXPIRR*MXPNGAS,MXGSDIM)
      if (NTEST.ge.1) write(LUWRT,*) 'Highest GAS dimension is ',MXGSDIM
      call findhigh(SYMDIM(1),NIRR_DG,IMXDDIM)
      call findhigh(NORB_FULL(1),NIRR_PN,IMXDIM_F)
      if (IMXDDIM.gt.MXNSPIPDMGRP) then
        write(LUWRT,*) 'Spinor expansion too long for 
     & double dim. arrays.'
        write(LUWRT,*) 'MXNSPIPDMGRP is  ',MXNSPIPDMGRP
        write(LUWRT,*) 'IMXDDIM is ',IMXDDIM
        write(LUWRT,*) 'Adjust MXNSPIPDMGRP in mxpdim.inc or 
     & reconsider.'
        Call Abend2( 'Quitting.' )
      end if
      call findhigh(SYDI(1),NIRR_DG,IMXDIM)
      INTSPACE = (2*MXGSDIM)**2*NGASPI**2*32
      ITRALINE = (2*MXGSDIM)**2*NGASPI**2
*
      if (ENVIRO(1:6).eq.'DIRAC ') then
        call memmar(KBOSYM,2*NDIM_TOT,'ADDS  ',1,'IBOSUM')
        call memmar(KISREO,2*NDIM_TOT,'ADDS  ',1,'IISREO')
      end if
*
      if (SIGDEN_ROUTE.eq.'SIGDEN') then
!
!       calculate # of diagonal integrals
        IRIMLT = MIN(NZ,2)
        NDIAMINT = (2*NDIM_TOT + 3*((2*NDIM_TOT)**2)) * IRIMLT/2
!
        IF( NDIAMINT .gt. MXINT_DIA ) THEN
          WRITE(LUWRT,*) ' STOP in subroutine itrctl !!!'
          WRITE(LUWRT,*) ' Number of diagonal integrals larger'
          WRITE(LUWRT,*) ' than MXINT_DIA: ', NDIAMINT, MXINT_DIA
          WRITE(LUWRT,*) ' increase MXINT_DIA (mxpdim.inc): to', 
     &                     NDIAMINT
          CALL QUIT(
     &    ' *** Error in itrctl *** too low value of MXINT_DIA')
        END IF
C
C       ... property run. This subroutine is called twice. 
C       1st time: do reordering - no slaves enter this routine
C       2nd time: all CPUs enter this routine - calculate 
C                 necessary information and allocate required arrays. 
C                 Reordered property matrix is communicated outside 
C                 this routine.
        IF( CIRUN .eq. 'PROP1 ')THEN
           ALLOC_INT = .TRUE.
           IF( .NOT. DOPROPREOD ) ALLOC_INT = .FALSE.
           GOTO 100
        END IF
C       ... other run types
        ALLOC_INT = .TRUE.
C
C       allocate integral/density array (if necessary)
C
C       IT_SHL == - 2    :  if MPI-shared memory AND CIRUN == .IJKLRO
C       IT_SHL == - 1    :  if not MPI-shared memory
C       IT_SHL == 0 or 1 :  if MPI-shared memory
C
        IF( IT_SHL .eq. - 2 )THEN
C
C         what may happen if we do a little DIRTY trick?
C         simulating local disks but having a global disk!!!
C         will it crash??? SK - Feb08
C
          IF( NFLGRPS_REL .eq. 1 ) THEN
            IF( MYPROC .ne. MASTER ) ALLOC_INT = .FALSE.
          ELSE
            IF( MYNEW_ID .ne. N_MASTER ) ALLOC_INT = .FALSE.
          END IF
C
        ELSE IF( IT_SHL .eq. - 1 )THEN
C
          IF( CIRUN .eq. 'IJKLRO' )THEN
            IF( NFLGRPS_REL .eq. 1 ) THEN
              IF( MYPROC .ne. MASTER ) ALLOC_INT = .FALSE.
            ELSE
              IF( MYNEW_ID .ne. N_MASTER ) ALLOC_INT = .FALSE.
            END IF
          ELSE
            IF( SPLIT_IJKL ) ALLOC_INT = .FALSE.
          END IF
C
        ELSE
          ALLOC_INT = .FALSE.
        END IF
C
        EX = .FALSE.
C
        IF( CIRUN .eq. 'IJKLRO' )THEN
          IF( ALLOC_INT )THEN
            INQUIRE(FILE='IJKL_REOD',EXIST=EX)
C
            IF( EX ) THEN
              WRITE(LUWRT,*)' '
              WRITE(LUWRT,*)' using provided integral file '
              WRITE(LUWRT,*)' IJKL_REOD! '
              WRITE(LUWRT,*)' '
            END IF
            IXXX_CLR = 1
          ELSE  
            IXXX_CLR = 2
          END IF
CSK       WRITE(6,*) ' my EX result: ALLOC_INT',myproc,ex,ALLOC_INT
#if defined(VAR_MPI2)
          IXXX_KEY  = MYPROC
          IXCOMM_SZ = 0
          IXCOMM_RK = 0
          IXCOMM    = 0
          CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,IXXX_CLR,
     &                         IXXX_KEY,IXCOMM,IERR)
          CALL MPI_COMM_SIZE(IXCOMM,IXCOMM_SZ,IERR)
          CALL MPI_COMM_RANK(IXCOMM,IXCOMM_RK,IERR)
          IF( ALLOC_INT )THEN
            MY_CHECK_F_LEN = IXCOMM_SZ
            CALL MPIXMEM_ALLOC( MY_CHECK_F_PTR, 2, MY_CHECK_F_LEN,
     &                          MPI_INFO_NULL, .FALSE. )
C
            CALL IZERO(ICHECK_FILE,IXCOMM_SZ)
C
            IEX_GAT   = 0 
C
            IF( EX ) THEN
              IEX_GAT = 1
            ELSE
              IEX_GAT = 0
            END IF
C
CSK            CALL MPI_BARRIER(IXCOMM,IERR)
CSK            WRITE(6,*) ' MYPROC: IEX_GAT is',MYPROC, IEX_GAT
C
            CALL MPI_GATHER(IEX_GAT,1,my_MPI_INTEGER,ICHECK_FILE,1,
     &                      my_MPI_INTEGER,0,IXCOMM,IERR)
CSK            CALL MPI_BARRIER(IXCOMM,IERR)
            NPTESTVAR = 0
CSK         NPTESTVAR = 20
            IF( NPTESTVAR .ge. 20 )THEN
              CALL MPI_BARRIER(IXCOMM,IERR)
              DO IJK = 1, IXCOMM_SZ
                CALL MPI_BARRIER(IXCOMM,IERR)
                IF( IXCOMM_RK .eq. IJK - 1 )THEN
                  WRITE(6,*) ' my ICHECK_FILE array', 
     &                         MYPROC, MYNEW_ID, IXCOMM_RK
                  DO JJJ = 1, IXCOMM_SZ
                     WRITE(6,*) ' ICHECK_FILE(',jjj,') = ',
     &                            ICHECK_FILE(jjj)
                  END DO
                END IF
              END DO
            END IF
            IF( IXCOMM_RK .eq. 0 )
     &      CALL CHECK_CONSIST( ICHECK_FILE, IXCOMM_SZ )
            CALL MPIXMEM_FREE( ICHECK_FILE )
          END IF
          CALL MPI_COMM_FREE(IXCOMM,IERR)
C         consistency is checked; update MPI_COMM_WORLD
          CALL MPI_BCAST(EX,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
#endif
        END IF
C       CIRUN == IJKLRO
        IF( EX ) THEN
          GOTO 1011
        END IF
C
 100    CONTINUE
C
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
        IF( CIRUN .eq. 'IJKLRO' ) THEN
          ALLOC_INT_SCR = ALLOC_INT
          ALLOC_INT     = .FALSE.
        END IF
#endif
        IF( CIRUN .eq. 'IJKLRO' )THEN
!         using low-memory sorting scheme (disk-I/O driven)
          IF(LOWSRT_IJKL) THEN 
            CALL MEMMAR(KT_CC,           0,'ADDS  ',2,'KT_CCX')
            CALL MEMMAR(KMDIASV,4*NDIAMINT,'ADDS  ',1,'MDIAXX')
            CALL IZERO(WORK(KMDIASV),4*NDIAMINT)
            CALL IZERO(IPOSDIA,NDIAMINT)
            GOTO 200
          END IF
        END IF
!
!       allocate array for diagonal integrals
        IF( ALLOC_INT )THEN
          call memmar(KMDIASV,4*NDIAMINT,'ADDS  ',1,'MDIASV')
          CALL IZERO(IPOSDIA,NDIAMINT)
          CALL IZERO(WORK(KMDIASV),4*NDIAMINT)
          IF(CIRUN .ne. 'PROP1 ' ) WRITE(LUWRT,'(A,I16,A/)') 
     &         '  Allocating for ',NDIAMINT,' diagonal integrals.'
        ELSE
          call memmar(KMDIASV, 0,'ADDS  ',1,'MDIASV')
          CALL IZERO(IPOSDIA,NDIAMINT)
        END IF
!
!       allocate full integral array
        IF( ALLOC_INT ) THEN
          IF(CIRUN .ne. 'PROP1 ' ) WRITE(LUWRT,'(/A,I16,A)') 
     &         '  Allocating for ',NALLINT,' integrals.'
          CALL MEMMAR_I8(KT_CC,NALLINT,'ADDS  ',2,'KT_CC ')
          CALL DZERO8(WORK(KT_CC),NALLINT)
        ELSE
          IDUMMY_LENGTH = 0
          CALL MEMMAR_I8(KT_CC,IDUMMY_LENGTH,'ADDS  ',2,'KT_CC ')
        END IF
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
        IF( CIRUN .eq. 'IJKLRO' ) ALLOC_INT = ALLOC_INT_SCR
#endif
 200    CONTINUE
      end if
C     ^ SIGDEN route
*
*
************************************************
* Set mark for flushing local memory
************************************************
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'ITRCTL')
*
*  2) Local memory
      if (ENVIRO(1:6).eq.'DIRAC ') then
        NUM_ONEI = (2*NDIM_TOT)**2
        call memmar(KDUMAR,NDIM_TOT*2,'ADDL  ',2,'DUMARR')
        call memmar(KIDUMAR,NDIM_TOT*2,'ADDL  ',1,'IDUMAR')
        call memmar(KDIONEE,2*NUM_ONEI,'ADDL  ',2,'IDIONE')
        NUM_TWEI = (2*NDIM_TOT)**4
        if (INTIMP.eq.6) then
          call memmar(KDITWOE,2*NUM_TWEI,'ADDL  ',2,'IDITWO')
          call memmar(KINDK,NUM_TWEI,'ADDL  ',2,'INDK  ')
          call memmar(KINDL,NUM_TWEI,'ADDL  ',2,'INDL  ')
          call memmar(KIOFF12,4*(2*NDIM_TOT)**2,'ADDL  ',1,'IOFF12')
          call memmar(KNUMINT,4*(2*NDIM_TOT)**2,'ADDL  ',1,'NUMINT')
          call memmar(KTPLOC,(2*NDIM_TOT)**2,'ADDL  ',1,'ITPLOC')
        end if
        call memmar(KKR,(2*NDIM_TOT+1),'ADDL  ',1,'IKR   ')
        call memmar(KINVREO,2*NDIM_TOT+1,'ADDL  ',1,'INVREO')
      end if
*
*
      if (IFOCKEX.eq.0) then
        if (ENVIRO(1:6).eq.'DIRAC ') then
          if (NTEST.ge.1) then
           write(LUWRT,*)
           write(LUWRT,*) '********************************************'
           write(LUWRT,*) '* Distributing MOLTRA transformed          *'
           write(LUWRT,*) '* dirac one- (and two-) electron integrals *'
           write(LUWRT,*) '* using input integral reordering          *'
           write(LUWRT,*) '********************************************'
          end if
*
          IRIOFF = (2*NDIM_TOT)**2
          NDIM_TOT2 = NDIM_TOT*2
          if (SIGDEN_ROUTE.eq.'TRADIT') then
            call dist_ints_dirac(WORK(KPICHE2),WORK(KDIONEE),
     &                           WORK(KDIONEE+IRIOFF),
     &                           WORK(KPICI1),WORK(KBOSYM),
     &                           WORK(KISREO),WORK(KINVREO),
     &                           NL2D,IUB2D,ISM2D,
     &                           ITRALINE,NDIM_TOT,NDIM_TOT2,
     &                           IHAM12,ISPINFREE,IMXDIM,DOUGRP)
C
C         distribute 1-/2- e- integrals in SIGDEN formalism
C
          else if ( SIGDEN_ROUTE .eq. 'SIGDEN' ) then
C
            call memmar(KISTR_CA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_CA')
            call memmar(KISTR_CB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_CB')
            call memmar(KISTR_AA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_AA')
            call memmar(KISTR_AB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_AB')
C
C           orbital symmetry info
            call alloc(ibeig,norbt)
            call izero(ibeig,norbt)
#if defined (VAR_MPI)
            if(CIRUN.ne.'PROP1 ')then
              call mpi_bcast(linear,1,MPI_LOGICAL,master,
     &                       mpi_comm_world,ierr)
            end if
#endif
            if(linear)then
              if(myproc.eq.master)then
                luni = -1
                inquire(file='KRMCSCF',exist=ex_mc,opened=fileopen,
     &                  number=luni)
                if(ex_mc.and.fileopen)then
                  rewind(luni)
                else if(ex_mc.and.(.not.fileopen))then
                  luni = 99
                  call opnfil(luni,'KRMCSCF','OLD    ','itrctl')
                  rewind(luni)
                end if
                call ireakrmc(luni,'MJVEC   ',ibeig,norbt)
                if(.not.fileopen) close(luni,status='keep')
              end if
#if defined (VAR_MPI)
              if(CIRUN.ne.'PROP1 ')then
                call mpi_bcast(ibeig,norbt,my_MPI_INTEGER,master,
     &                         mpi_comm_world,ierr)
              end if
#endif
            end if
            
 
C
            IF( REORD_IJKL ) THEN
              IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' ) THEN
#if defined (VAR_MPI2)
C
                  WRITE(LUWRT,'(/A)') 
     &            '  ===================================='
                  WRITE(LUWRT,'(A)') 
     &            '  Reading from integral file IJKL_REOD'
                  WRITE(LUWRT,'(A/)') 
     &            '  ===================================='
#endif
C
                CALL READ_IJKL_GAS(WORK(KT_CC),T_BUFF,T_BUFF_D, 
     &                             IT_TTPL,WORK(KMDIASV),NDIAMINT,
     &                             NSPOBEX_TP,N1ELINT)
C
                call dealloc(ibeig)
                GOTO 1001
C
              ELSE IF( CIRUN .eq. 'IJKLRO' ) THEN
C
#if defined (VAR_MPI2)
C
                IF( ALLOC_INT )THEN
                  IXXX_CLR = 1
                ELSE  
                  IXXX_CLR = 2
                END IF
                IXXX_KEY = MYPROC
                CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,IXXX_CLR,
     &                              IXXX_KEY,IXCOMM,IERR)
C
                IF( IXXX_CLR .eq. 2 ) GOTO 797
#endif
C
                WRITE(LUWRT,*)' '
                WRITE(LUWRT,*)' Resorting integrals...'
                WRITE(LUWRT,*)' '
C
                CALL WRITE_IJKL_GAS( WORK(KT_CC), T_BUFF,
     &                               INTIMP,ISPINFREE,ISYM_T,
     &                               WORK(KLSOBEX),WORK(KISTR_CA),
     &                               WORK(KISTR_CB),WORK(KISTR_AA),
     &                               WORK(KISTR_AB),NDIM_TOT2,
     &                               F1,WORK(KMDIASV),IBOSYM,
     &                               WORK(KLABEXTP),CIRUN,NTEST,
     &                               NDIAMINT,ENVIRO,ibeig,SIGDEN_ROUTE,
     &                               INTSPACE)
C
 797            CONTINUE
C
#if defined (VAR_MPI2)
                CALL MPI_COMM_FREE(IXCOMM,IERR)
#endif
C
                call dealloc(ibeig)
                GOTO 1001
              END IF
C             ^ CIRUN ???
            END IF
C           ^ REORD_IJKL
C
            IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' )THEN
C
              IF( IT_SHL .ge. 0 .or. SPLIT_IJKL )THEN
                WRITE(LUWRT,*) ' Running in MPI-shared memory mode     '
                WRITE(LUWRT,*) ' or making use of revised integral     '
                WRITE(LUWRT,*) ' allocation (.IJKLSP) requires         '
                WRITE(LUWRT,*) ' the integral file IJKL_REOD.          '
                WRITE(LUWRT,*) ' Programm will abort now.              '
                WRITE(LUWRT,*) ' Restart job with option .IJKLRO !     '
                CALL QUIT(
     &          '*** Error in itrctl *** Option .IJKLRO required. ')
              END IF
C
              CALL MEMMAR_I8(KINTI12,LEN_ALL_INT,'ADDL  ',2,'XTIT12')
C
              if(iiomod_rel.eq.1 .or. myproc.eq.master)then

!               get dimensions for scratch matrices
                call rgeth2_dim(ndim_ij,ndim_kl,ndim_h12_tmp)

                call memmar(kdumx12,idum,'MARK  ',idum,'geth12')
!               allocate scratch matrices
                call memmar(khmat,ndim_h12_tmp,'ADDL  ',2,'hmat12')
                call memmar(kh2xy,ndim_h12_tmp,'ADDL  ',2,'h12xy ')
                call memmar(kh2yx,ndim_h12_tmp,'ADDL  ',2,'h12yx ')

                call memmar(kdxb12, 2*ndim_ij,'ADDL  ', 1,'idxb12')
                call memmar(kdxb34, 2*ndim_kl,'ADDL  ', 1,'idxb34')
                call memmar(kdxkr12,2*ndim_ij,'ADDL  ', 1,'idxk12')

!               get integrals on work(kinti12)
                call rgeth21(.false.,.true.,dummy,dummy,work(kinti12),
     &                       work(khmat),work(kh2xy),work(kh2yx),
     &                       work(kdxkr12),work(kdxb12),work(kdxb34),
     &                       ibeig,ntest,luwrt)

!               release scratch space
                call memmar(kdumx12,idum,'FLUSM ',idum,'geth12')
              end if
!             ^ iiomod_rel == 1 or myproc == master
#if defined (VAR_MPI2)
              IF( IIOMOD_REL .eq. 0 )THEN
*               ... co-workers get the integrals
                LEN_INT_BCAST = 0
#if defined (INT_STAR8)
                LEN_INT_BCAST = LEN_ALL_INT
#else
                XLENALL = REAL(LEN_ALL_INT)
                LEN_INT_BCAST = IGIVE_I_B(XLENALL)
#endif
                CALL MPI_BCAST(WORK(KINTI12),LEN_INT_BCAST,MPI_REAL8,
     &                              MASTER,MPI_COMM_WORLD,IERR)
              END IF
#endif
C
C             integrals on WORK(KINTI12) / F1-F2 --> 
C             distribute to WORK(KT_CC) resp. T_BUFF!
C
              call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                              WORK(KLSOBEX),WORK(KT_CC),
     &                              T_BUFF,
     &                              WORK(KISTR_CA),WORK(KISTR_CB),
     &                              WORK(KISTR_AA),WORK(KISTR_AB),
     &                              NDIM_TOT2,F1,WORK(KINTI12),
     &                              WORK(KMDIASV),IBOSYM,
     &                              WORK(KLABEXTP),CIRUN,
     &                              NTEST)
            ELSE IF( CIRUN .eq. 'PROP1 ') THEN
              IF( DOPROPREOD )THEN
CSK           NTEST = 20
              call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                              WORK(KLSOBEX),WORK(KT_CC),
     &                              T_BUFF,
     &                              WORK(KISTR_CA),WORK(KISTR_CB),
     &                              WORK(KISTR_AA),WORK(KISTR_AB),
     &                              NDIM_TOT2,F1,DUMMY,
     &                              WORK(KMDIASV),IBOSYM,
     &                              WORK(KLABEXTP),CIRUN,
     &                              NTEST)
CSK           NTEST = 00
              END IF
            ELSE
              call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                              WORK(KLSOBEX),WORK(KT_CC),
     &                              T_BUFF,
     &                              WORK(KISTR_CA),WORK(KISTR_CB),
     &                              WORK(KISTR_AA),WORK(KISTR_AB),
     &                              NDIM_TOT2,F1,F2,
     &                              WORK(KMDIASV),IBOSYM,
     &                              WORK(KLABEXTP),CIRUN,
     &                              NTEST)
C
            END IF
C           ^ CIRUN eq CIINII .or. CIRUN .eq. 'KR-CI '?
            call dealloc(ibeig)
          end if
C         ^ SIGDEN_ROUTE ?
        end if ! (ENVIRO(1:6).eq.'DIRAC ')
      else
        WRITE(LUWRT,'(/A,I3)') ' Problem in itrctl: IFOCKEX is ',IFOCKEX
        CALL QUIT( ' *** ERROR in ITRCTL: unknown IFOCKEX value ***' )
      end if
C     ^ end if external Fock matrix active or not
C
      if (SIGDEN_ROUTE.eq.'SIGDEN' .and. CIRUN .ne. 'PROP1 ' ) then
        write(LUWRT,'(/A,I16)') ' Number of integrals/CC amplitudes: ',
     &                            NALLINT
        write(LUWRT,'(A/)') ' (including imaginary parts) '
      end if
C
1001  CONTINUE
C
C     flush local memory to mark set at beginning of do loop
C
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'ITRCTL')
C
1011  CALL QEXIT('ITRCTRL')
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  This piece of code is supposed to transform any type
*  of one-electron integrals from a basis of atomic functions
*  blocked by spatial symmetry to a basis of molecular
*  spinors.
*
      subroutine oneltra(SPININTS,
     &                   AOTRISM,SOINXLINE,SOINYLINE,SOINZLINE,
     &                   AOINT,HELPMAT2,TRAINT,COEF,
     &                   ADJCO,PROD,
     &                   SOINTZ,SOINTX,SOINTY,
     &                   ISMIND,IGSOFF,IGSIND,
     &                   ITRALINE,IMXDIM,MXGSDIM,
     &                   IOUTLST,ICURRDIM,NTRIAEL,
     &                   INTITER,ITRA_TYPE,
     &                   INTSPACE,ENVIRO,IFOCKEX)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
*  for lblas, COORD for x,y,z spin-orbit
      character*1 transa,transb,COORD
      character*72 ENVIRO
*
      logical SPIORB
*
      dimension SPININTS(ITRALINE,32)
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      dimension AOINT(ICURRDIM,ICURRDIM),
     &     COEF(ICURRDIM,ICURRDIM,2),
     &     ADJCO(ICURRDIM,ICURRDIM,2),
     &     PROD(ICURRDIM,ICURRDIM,2),
     &     AOTRISM(ITRIASUM),
     &     HELPMAT2((2*IMXDIM)**2),
     &     TRAINT(ICURRDIM,ICURRDIM,4),
     &     SOINTZ(ICURRDIM,ICURRDIM),
     &     SOINTX(ICURRDIM,ICURRDIM),
     &     SOINTY(ICURRDIM,ICURRDIM),
     &     SOINZLINE(NTRIAEL),
     &     SOINXLINE(NTRIAEL),
     &     SOINYLINE(NTRIAEL)
*
      dimension ISMIND(ICURRDIM),
     &     IGSOFF(2*NGASPI),
     &     IGSIND(ICURRDIM,NIRR_PN),
     &     IAD2(3)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      write(*,*) '****************************'
      write(*,*) '                            '
      write(*,*) ' This is oneltra in action  '
      write(*,*) '                            '
      write(*,*) '****************************'
*
*     We are transforming integrals with two general indices
*
 300  format(4G18.12)
*
*  transpose coefficients ...
      IMRECO = 1
      call adjoint(ADJCO(1,1,1),COEF(1,1,1),
     &     SYMDIM(IDIMGRP),
     &     SYMDIM(IDIMGRP),IMRECO)
*
*  ... and imaginary parts of coefficients
      IMRECO = -1
      call adjoint(ADJCO(1,1,2),COEF(1,1,2),
     &     SYMDIM(IDIMGRP),
     &     SYMDIM(IDIMGRP),IMRECO)
*
*  What type of operator?
      if (IOPER.eq.1.or.IOPER.eq.3) then
*
*  Select and arrange integrals
         SPIORB = .FALSE.
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(AOINT(1,1),AOTRISM(1),HELPMAT2(1),
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
 10      continue
*
         if (NTEST.ge.1) then
            write(*,*) 'Full AO integral matrix for this'
            write(*,*) 'transformation run reads:'
            write(*,*) 'AOINT'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (AOINT(I,J),J=1,SYMDIM(IDIMGRP))
            end do
         end if
*
*     HEAVY DUTY SECTION ...
*     Transform
*
*     Do integrals x real
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,AOINT,SYMDIM(IDIMGRP),
     &              COEF(1,1,1),SYMDIM(IDIMGRP),beta,PROD(1,1,1),
     &              SYMDIM(IDIMGRP))
*
*     Do integrals x imag.
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,AOINT,SYMDIM(IDIMGRP),
     &              COEF(1,1,2),SYMDIM(IDIMGRP),beta,PROD(1,1,2),
     &              SYMDIM(IDIMGRP))
*
         if (NTEST.gt.2) then
            write(*,*) 'The halftransformed integrals read (real):'
            write(*,*) 'PROD'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (PROD(I,J,1),J=1,SYMDIM(IDIMGRP))
            end do
            write(*,*) 'The halftransformed integrals read (imag):'
            write(*,*) 'PROD * i'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (PROD(I,J,2),J=1,SYMDIM(IDIMGRP),1)
            end do
         end if
*
*     Do second half transformation
*
*     First block: real x real
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,1),
     &              SYMDIM(IDIMGRP),PROD(1,1,1),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,1),SYMDIM(IDIMGRP))
*
*     Second block: imag. x imag.
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,2),
     &              SYMDIM(IDIMGRP),PROD(1,1,2),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,2),SYMDIM(IDIMGRP))
*
*     Include factor i*i in second real block:
*
         do I=1,SYMDIM(IDIMGRP),1
            do J=1,SYMDIM(IDIMGRP),1
               TRAINT(J,I,2) = - TRAINT(J,I,2)
            end do
         end do
*
*     Add real parts
*
         do I=1,SYMDIM(IDIMGRP),1
            do J=1,SYMDIM(IDIMGRP),1
               TRAINT(J,I,1) = TRAINT(J,I,1) + TRAINT(J,I,2)
            end do
         end do
*
*     Third block: real x imag.
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,1),
     &              SYMDIM(IDIMGRP),PROD(1,1,2),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,3),SYMDIM(IDIMGRP))
*
*     Fourth block: imag. x real
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,2),
     &              SYMDIM(IDIMGRP),PROD(1,1,1),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,4),SYMDIM(IDIMGRP))
*     add
         do I=1,SYMDIM(IDIMGRP),1
            do J=1,SYMDIM(IDIMGRP),1
               TRAINT(J,I,2) = TRAINT(J,I,3) + TRAINT(J,I,4)
            end do
         end do
*
         if (NTEST.ge.1) then
            write(*,*) 'The fully transformed integrals read:'
            write(*,*) 'TRAINT (real)'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (TRAINT(I,J,1),J=1,SYMDIM(IDIMGRP),1)
            end do
            write(*,*) 'TRAINT * i (imag.)'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (TRAINT(I,J,2),J=1,SYMDIM(IDIMGRP),1)
            end do
         end if
*
*     Calculate appropriate integral output list and
*     write transformed integrals to list
*
         do K=1,2,1
            INTLIST = 2*IOUTLST - 2 + K
*
*  Fold and write to appropriate integral output list!
*
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(*,*) 'are being passed to intgas'
             end if
            end if
*
            if (ITRA_TYPE.le.3) then
*  Distribute to GAS matrices
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Provide space for transformed integrals for MCSCF run
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
      end if
*
* =====================================
*  Now the spin-orbit integrals ...
* =====================================
*
      if (IOPER.eq.2.or.IOPER.eq.3) then
         SPIORB = .TRUE.
         COORD = 'Z'
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(SOINTZ,AOTRISM,SOINZLINE,
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
*  Call transformation routine
*
         call sotra(SOINTZ,TRAINT,ADJCO,COEF,PROD,
     &              IOUTLST,SYMDIM(IDIMGRP),
     &              ITOTDM(1),COORD,DOUGRP)
*
*     Write to appropriate output list
*
         do K=1,2,1
            INTLIST = 2*IOUTLST + K + 22
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             end if
            end if
            if (ITRA_TYPE.le.3) then
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
*
         SPIORB = .TRUE.
         COORD = 'X'
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(SOINTX,AOTRISM,SOINXLINE,
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
         call sotra(SOINTX,TRAINT,ADJCO,COEF,PROD,
     &              IOUTLST,SYMDIM(IDIMGRP),
     &              ITOTDM(1),COORD,DOUGRP)
*
         do K=1,2,1
            INTLIST = 2*IOUTLST + K + 6
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             end if
            end if
            if (ITRA_TYPE.le.3) then
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
*
         SPIORB = .TRUE.
         COORD = 'Y'
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(SOINTY,AOTRISM,SOINYLINE,
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
         call sotra(SOINTY,TRAINT,ADJCO,COEF,PROD,
     &              IOUTLST,SYMDIM(IDIMGRP),
     &              ITOTDM(1),COORD,DOUGRP)
*
         do K=1,2,1
            INTLIST = 2*IOUTLST + K + 14
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             end if
            end if
            if (ITRA_TYPE.le.3) then
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
      end if
*
      if (IOPER.le.0.or.IOPER.gt.3) then
         write(*,*) 'IOPER out of range.'
         write(*,*) 'IOPER =  ',IOPER
         Call Abend()
      end if
*
 1000 return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
****************************************************************
***                     PICASSO                              ***
***                                                          ***
***   P rogram						     ***
***     for 						     ***
***   I ntelligent                                           ***
***   C omplete                                              ***
***   A ctive                                                ***
***   S pace                                                 ***
***   S pin                                                  ***
***   O rbit                                                 ***
***     calculations                                         ***
***         (at least what it was once planned to become...) ***
***                                                          ***
****************************************************************
***
*** Written by Timo Fleig, Bonn
***        and Jeppe Olsen, Lund
***
*** October 1996 - February 1998
***
*** (poor picasso is merely one of LUCIAREL's subroutines ... :-(
***
      subroutine picasso(INTSPACE,ISOOP,IPRTRA,MCSCF,
     &                   NL2D,IUB2D,ISM2D,F1,F2,T_BUFF,
     &                   T_BUFF_D,IT_TTPL,IBOSYM,CIRUN,RUN_WO_T_BUFF)
*
*******************************************************************
*  !!!  INATIT states iteration point for natural spinor method !!!
*******************************************************************
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      integer POI(30)
      logical MCSCF
      logical RUN_WO_T_BUFF
      common /WORKPIC/IPOI
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
#include "onel.inc"
#include "crun.inc"
C
      IF( RUN_WO_T_BUFF ) GOTO 9999
C
* for integral summation
      ISUMNEW = 1
* for global purposes
      INTITER = 1
      N_SPINR_EXP = 2
*
C     write(6,*) '     *******************************************'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *  This is PICASSO,                       *'
C     write(6,*) '     *  (Transformation section for            *'
C     write(6,*) '     *   one-electron sf and SO integrals      *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *******************************************'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *        0.1   stoneage version           *'
C     write(6,*) '     *        0.2   GAS version                *'
C     write(6,*) '     *        0.2.1  AO ints from MOLECULE     *'
C     write(6,*) '     *        0.2.2  SOMF ints from AMFI       *'
C     write(6,*) '     *        0.2.2.1 GAS-DG dependent         *'
C     write(6,*) '     *        0.2.3  general trafo CI/MCSCF    *'
C     write(6,*) '     *        0.2.3.1 SGI portable version     *'
C     write(6,*) '     *        0.3   MEMMAN allocation included *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *******************************************'
*
C     write(6,*) '     *******************************************'
C     write(6,*) '     *                                         *'
C     write(6,*) '     * written by                              *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *      TIMO FLEIG, University of Bonn     *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *           and                           *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *      JEPPE OLSEN, University of Lund    *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *******************************************'
*
*     determine double group information
      call sspir(ISOI)
*
*     Move to integral transformation section:
*     Define type of transformation
*     Here for now: complete two-component trans for CI
*
      if (MCSCF) then
        ITRA_TYPE = 4
      else
        if (INATIT.eq.1) then
          ITRA_TYPE = 2
        else if (INATIT.ge.2) then
          ITRA_TYPE = 4
          ORBFILE = .true.
        end if
      end if
*
      call itrctl(INTSPACE,ITRA_TYPE,INTITER,ENVIRO,
     &            ISPINFREE,IFOCKEX,NORB_F,NL2D,IUB2D,ISM2D,ISYM_T,
     &            SIGDEN_ROUTE,F1,F2,T_BUFF,T_BUFF_D,IT_TTPL,IBOSYM,
     &            MCSCF,CIRUN)
*
9999  return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Reorder SO integrals according to matrices INDREO
*
      subroutine reo_soint(SOIN,HELPMAT,HELPMAT2,IDIM,
     &                     INDREO,NTRIAIN,NMAXBAS,IPRINT)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      dimension SOIN(NTRIAIN),HELPMAT(IDIM,IDIM),
     &          HELPMAT2(IDIM,IDIM)
      dimension INDREO(NMAXBAS)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      call unfold(HELPMAT,SOIN,IDIM,1,NTRIAIN,1)
      call matsort(HELPMAT,HELPMAT2,INDREO,INDREO,IDIM,IDIM)
*
      if (NTEST.ge.3) then
        write(6,*) 'Checking reordered SO integral matrix: '
        call wrtmat(HELPMAT,IDIM,IDIM,IDIM,IDIM)
      end if
*
      call fold(SOIN,HELPMAT,IDIM,1,NTRIAIN)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Transform spin-orbit integrals
*  l_z, l_x, l_y
*
      subroutine sotra(SOINT,TRAINT,ADJCO,COEF,PROD,
     &                 IOUTLST,IDIM,ITOTDM,
     &                 COORD,IDBG)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
      character*2 transa,transb
      character*1 COORD
*
#include "mxpdim.inc"
#include "input.inc"
*
      dimension SOINT(IDIM,IDIM),
     &          TRAINT(IDIM,IDIM,4),
     &          ADJCO(IDIM,IDIM,2),
     &          COEF(IDIM,IDIM,2),
     &          PROD(IDIM,IDIM,2)
*
      dimension ITOTDM(MXNSPIR)
*
      NTESTL = 00
      NTEST = MAX(NTESTL,IPRINT)
*
      if (NTEST.ge.10) then
        write(6,*) 'Checking input block in sotra:'
        call wrtmat(SOINT,IDIM,IDIM,IDIM,IDIM)
      end if
*
 300  format(4G18.12)
      ZERO = 0.0D0
*
*  Erase TRAINT
      do INUM=1,4,1
        do ICOL=1,IDIM,1
          call setvec(TRAINT(1,ICOL,INUM),ZERO,IDIM)
        end do
      end do
*
*     Write to symmetry matrix (see notes for structure)
*
*     include spin integration factor according to
*       <alpha| s(z) |alpha> = 1
*       <beta| s(z) |beta> = -1
*
*       <alpha| s(x) |beta> = 1
*       <beta| s(x) |alpha> = 1
*
*       <alpha| s(y) |beta> = -i
*       <beta| s(y) |alpha> = i
*
* This scheme also applies for quaternionic double groups.
* It is not relevant which spatial symmetries are coupled by
* the spin-orbit operator. In this transformation, the structure
* of the trafo coefficient matrices is always the same, i.e.
* it has the same U/B alpha/beta structure. So the l integral
* matrices are set up in the same way, no matter which group
* is considered. sf and SO(z) SAO integrals over the l operator
* remain in the diagonal blocks, SO(x,y) integrals off-diagonal.
* The difference in coupled spatial symmetries is already taken
* care of in the atomic l integrals, and the corresponding
* expansion coefficients are copied into the transformation
* matrix accordingly.
*
* Caution: When the spinor symmetries G4 (or G6, G8 with respect
*          to C_2v^*) are being treated, the coefficient matrices
*          are rearranged, so the sign has to be applied in the
*          upper instead of the lower blocks!
*
      IACDIM = 1
      if (IOUTLST.gt.2) IACDIM = 3
      if (IDBG.eq.7.and.IOUTLST.eq.2) IACDIM = 2
*
      call evenodd(IEOD,IOUTLST)
      if (COORD.eq.'Y'.or.COORD.eq.'Z') then
        if (IEOD.eq.1) then
          do K=1,2,1
            do I=1,IDIM,1
              do J=ITOTDM(IACDIM)+1,IDIM,1
                TRAINT(J,I,K+2) = - COEF(J,I,K)
              end do
            end do
            do I=1,IDIM,1
              do J=1,ITOTDM(IACDIM),1
                TRAINT(J,I,K+2) = COEF(J,I,K)
              end do
            end do
          end do
        else if (IEOD.eq.2) then
          do K=1,2,1
            do I=1,IDIM,1
              do J=1,ITOTDM(IACDIM),1
                TRAINT(J,I,K+2) = - COEF(J,I,K)
              end do
            end do
            do I=1,IDIM,1
              do J=ITOTDM(IACDIM)+1,IDIM,1
                TRAINT(J,I,K+2) = COEF(J,I,K)
              end do
            end do
          end do
        end if
      else if (COORD.eq.'X') then
        do K=1,2,1
          do I=1,IDIM,1
            do J=1,IDIM,1
              TRAINT(J,I,K+2) = COEF(J,I,K)
            end do
          end do
        end do
      end if
*
*     Transform
*
*     Do integrals x real part
*     Result is imaginary
*
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,SOINT(1,1),
     &           IDIM,TRAINT(1,1,3),IDIM,
     &           beta,PROD(1,1,1),IDIM)
*     Do integrals x imag.
*     Result is (-1) real!!
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,SOINT(1,1),
     &           IDIM,TRAINT(1,1,4),IDIM,
     &           beta,PROD(1,1,2),IDIM)
*     include factor (-1)
      do I=1,IDIM,1
         do J=1,IDIM,1
            PROD(J,I,2) = - PROD(J,I,2)
         end do
      end do
*
*     Transform second index
*
*     Real x real
*
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,1),
     &           IDIM,PROD(1,1,2),IDIM,
     &           beta,TRAINT(1,1,1),IDIM)
*     imag. x imag.
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,2),
     &           IDIM,PROD(1,1,1),IDIM,
     &           beta,TRAINT(1,1,2),IDIM)
*
      do I=1,IDIM,1
         do J=1,IDIM,1
            TRAINT(J,I,2) = - TRAINT(J,I,2)
         end do
      end do
*
      do I=1,IDIM,1
         do J=1,IDIM,1
            TRAINT(J,I,1) = TRAINT(J,I,1) + TRAINT(J,I,2)
         end do
      end do
*     Real x imag.
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,1),
     &           IDIM,PROD(1,1,1),IDIM,
     &           beta,TRAINT(1,1,3),IDIM)
*     Imag. x real
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,2),
     &           IDIM,PROD(1,1,2),IDIM,
     &           beta,TRAINT(1,1,4),IDIM)
*
      do I=1,IDIM,1
         do J=1,IDIM,1
            TRAINT(J,I,2) = TRAINT(J,I,3) + TRAINT(J,I,4)
         end do
      end do
*
      if (NTEST.ge.3) then
        if (COORD.eq.'Y') then
          write(6,*) 'SO-Y before factor:'
          call wrtmat(TRAINT(1,1,2),IDIM,IDIM,
     &                IDIM,IDIM)
        end if
      end if
*
*  Include factor of i for the SO(y) operator from spin integration
*  These are the imaginary parts, so include factor (-1)
      if (COORD.eq.'Y') then
         do I=1,IDIM,1
            do J=1,IDIM,1
               TRAINT(J,I,3) = - TRAINT(J,I,2)
            end do
         end do
*
*  now interchange real and imaginary parts
         do I=1,IDIM,1
            do J=1,IDIM,1
               TRAINT(J,I,2) = TRAINT(J,I,1)
            end do
         end do
         do I=1,IDIM,1
            do J=1,IDIM,1
               TRAINT(J,I,1) = TRAINT(J,I,3)
            end do
         end do
      end if
*
      if (NTEST.ge.1) then
         write(*,*) 'The fully transformed SO ints. read:'
         write(*,*) 'TRAINT (real)  ',COORD
         call wrtmat(TRAINT(1,1,1),IDIM,IDIM,
     &               IDIM,IDIM)
         write(*,*) 'TRAINT (imag.)  ',COORD
         call wrtmat(TRAINT(1,1,2),IDIM,IDIM,
     &               IDIM,IDIM)
      end if
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine spitra2(I12O,I34O,I1234O,
     &                   ICCSMO,ISPITRA,ISPITRA_RC,
     &                   LSPINOR,IRIO,ICCSM,
     &                   NLISTO,I1234I,I12I,I34I,IOPSMIN,IOPSMOUT)
*
* Input
* ======
*
* LBATCHI : Length of batch of input integrals
* LBATCHO : Length of batch of output integrals
*
* I1234I : Permutational symmetry between particle and and two
* I12I  : Permutational symmetry between first two orbital indeces
* I34I  : Permutational symmetry between last two orbital indeces
* NIRR_PN : Number of irreps in input symmmetry group
*
* IOPSMIN : Symmetry of operator in input group
* IOPSMOUT : Symmetry of operator in output group
* NLISTO : Number of output lists
* I1234O : Permutational symmetry between particle one and two
*          for each output list
* I12O  : Permutational symmetry between the first two indeces
*         for each output list
* I34O : Permutational symmetry between the last two indeces
*        for each output list
*
*
* Main routine for the transformartion of
* two-electron integrals to from a.o. basis
* to spinors in double group symmetry
*
* Timo Fleig + Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
#include "input.inc"
*. Input
      INTEGER I1234O(NLISTO),I12O(NLISTO),I34O(NLISTO)
      INTEGER ICCSMO(NLISTO)
      INTEGER ISPITRA(4,NLISTO)
*
      INTEGER IBLOCKO(4,MXBLKO),LISTO(MXBLKO)
      INTEGER INDEXO(4,MXBLKO),ISIGNO(MXBLKO),IRIO(MXBLKO)
      INTEGER IBL1(6,MXBLKO),IBL2(6,MXBLKO),IBL3(6,MXBLKO)
      INTEGER IBL4(6,MXBLKO)
*
      dimension ISPITRA_RC(NLISTO)
      dimension LSPINOR(N_SPINR_EXP)

C     CALL QENTER('DBGT2')
      IDUM = 0
C     CALL MEMMAR(KDUM,IDUM,'MARK',IDUM,'SPITR2')
*
*
      WRITE(6,*) ' Control has been transferred to SPITRA2'
      WRITE(6,*)
      WRITE(6,*) ' I12I, I34I,I1234I',I12I,I34I,I1234I
      WRITE(6,*) 'ISPITRA : '
      CALL IWRTMA(ISPITRA,4,NLISTO,4,NLISTO)
*
*. Loop over symmetry-unique blocks of input ao integrals
*
      DO IPSM = 1, NIRR_PN
*. Symmetry Range of second index
       IF(I12I.EQ.1.OR.I12I.EQ.-1) THEN
         IMXQSM = IPSM
       ELSE
         IMXQSM = NIRR_PN
       END IF
       DO IQSM = 1, IMXQSM
*. Symmetry range of third index
        IF(I1234I.EQ.1 .OR. I1234I.EQ.-1) THEN
          IMXRSM = IPSM
        ELSE
          IMXRSM = NIRR_PN
        END IF
        DO IRSM = 1, IMXRSM
*. Symmetry range of fourth index
         IMXSSM = NIRR_PN
         IF (IPSM.EQ.IRSM.AND.(I1234I.EQ.1.OR.I1234I.EQ.-1))
     &   IMXSSM = IQSM
         IF(I34I.EQ.1.OR.I34I.EQ.-1) IMXSSM = MIN(IMXSSM,IRSM)
*. Required symmetry of index 4 to produce total operator of
*. correct symmetry
         IPQSM =  IPGMULT(IPSM,IQSM)
         IPQRSM = IPGMULT(IPQSM,IRSM)
         ISSM = IPGMULT(IPQRSM,IOPSMIN)
         WRITE(6,*) ' IPQSM IPQRSM ISSM',IPQSM,IPQRSM,ISSM
         IF(1.LE.ISSM.AND.ISSM.LE.IMXSSM) THEN
*
            IF(IPRINT.GE.2) THEN
              WRITE(6,*) ' Next unique input symmetry block'
              WRITE(6,*) ' IPSM,IQSM,IRSM,ISSM',
     &        IPSM,IQSM,IRSM,ISSM
            END IF
*. Which output  blocks corresponds to this input block
            CALL INTBLK_IN_TO_OUT(IPSM,IQSM,IRSM,ISSM,
     &           I1234I,I12I,I34I,NLISTO,I1234O,I12O,I34O,
     &           ICCSMO,IMOSF_SP,ISPSF_MO,NIRR_PN,NIRR_DG,
     &           ISPITRA,ISPITRA_RC,
     &           NBLOCKO,IBLOCKO(1,1),LISTO(1),INDEXO(1,1),
     &           ISIGNO,IRIO,NCOMP)
*
*     INPUT
*    -------
*           IPSM,IQSM,IRSM,ISSM	: spatial index symmetries
*           I1234I,I12I,I34I	: inlist particle and complexity symmetries
*           NLISTO		: number of output lists
*           ISPITRA		: spinor expansion of index and output list
*           ISPITRA_RC		: r/i/c type of transformation
*
*     OUTPUT
*    --------
*           I1234O,I12O,I34O	: outlist particle and complexity symmetries
*           NBLOCKO		: number of output blocks
*           IBLOCKO		: dg symmetry for index in part. block
*           LISTO		: output list of part. block
*           INDEXO		: permutation of index in this block
*           IRIO		: real or imaginary part of this block
*
*. Intermediate blocks
*. Should not be here but...
*
*
            NRORC = NCOMP
            CALL INTERM_BLK(NBLOCKO,IBLOCKO,INDEXO,LISTO,
     &           IRIO,ISPITRA,ISPITRA_RC,LSPINOR,NRORC,
     &           NLISTO,N_SPINR_EXP,
     &           NBL1,IBL1,NBL2,IBL2,NBL3,IBL3,NBL4,IBL4,
     &           MXBLKO)
*
*     INPUT
*    -------
*           NBLOCKO,IBLOCKO,INDEXO,
*           LISTO,IRIO			: output from intblk_in_to_out
*
*           LSPINOR			: size of spinor expansion sets
*           NRORC			: number of r/i components
*
*     OUTPUT
*    --------
*           NBL1,NBL2,NBL3,NBL4	: number of blocks with X indeces transformed
*           IBL1,IBL2,IBL3,IBL4	: specification of one-index transformed
*				  block as
*				   (1,*) spinor transformation
*				   (2,*) dg symmetry
*				   (3,*) real or imaginary part
*				   (4,*) transformation of index
*					   3 or 4  if X=1,2
*					   1 or 2  if X=3,4
*				   (5,*) # of X+1 index transformed blocks
*				         this X index transf. bl. will
*				         contribute to
*				   (6,*) first X+1 index transf. bl. this
*				         bl. will contribute to
*
         END IF
        END DO
       END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SPOBEX_TO_ABOBEX(ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &                            IFLAG,NAOBEX_TP,NBOBEX_TP,IAOBEX_TP,
     &                            IBOBEX_TP)
*
* Split spin-orbital excitations into alpha and beta-orbital excitations
*
* IFLAG = 1 : Find only number of alpha- and beta- orbital excitations
*
* Jeppe Olsen, July 2000
*
#include "implicit.inc"
#include "ipoist8.inc"
*
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
*. Output
      INTEGER IAOBEX_TP(2*NGAS,*),IBOBEX_TP(2*NGAS,*)
*
C?    WRITE(6,*) ' SPOBEX_TO_ABOBEX : NSPOBEX_TP = ', NSPOBEX_TP
      DO IAB = 1, 2
        LEX_TP = 0
        DO JSPOBEX_TP = 1, NSPOBEX_TP
*. Has this a or b excitation been observed before
          I_AM_OLD_HAT = 0
          DO KSPOBEX_TP = 1, LEX_TP
            I_DENTICAL = 1
            DO ICA = 1, 2
              IEXP_OFF = (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
              DO IGAS = 1, NGAS
                IF(ISPOBEX_TP(IGAS+IEXP_OFF-1,JSPOBEX_TP).NE.
     &             ISPOBEX_TP(IGAS+IEXP_OFF-1,KSPOBEX_TP)    )
     &             I_DENTICAL = 0
              END DO
            END DO
            IF(I_DENTICAL.EQ.1) I_AM_OLD_HAT = 1
          END DO
          IF(I_AM_OLD_HAT.EQ.0) THEN
            LEX_TP = LEX_TP + 1
            IF(IFLAG.EQ.0) THEN
              DO ICA = 1, 2
               IEXP_OFF =  (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
               IAB_OFF  =  (ICA-1)*NGAS + 1
               IF(IAB.EQ.1) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IAOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               ELSE IF (IAB.EQ.2) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IBOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               END IF
              END DO
            END IF
          END IF
*         ^ End of I_AM_OLD_HAT = 0
        END DO
*       ^ End of loop over JSPOBEX_TP
      IF(IAB.EQ.1) THEN
        NAOBEX_TP = LEX_TP
      ELSE
        NBOBEX_TP = LEX_TP
      END IF
      END DO
*     ^ End of loop over IAB
*
      NTEST = 000
      IF(NTEST.GE.3) THEN
       WRITE(6,*) ' Number of alpha-excitation operators ', NAOBEX_TP
       WRITE(6,*) ' Number of beta-excitation operators ',  NBOBEX_TP
      END IF
      IF(NTEST.GE.5.AND.IFLAG.EQ.0) THEN
        WRITE(6,*) ' Alpha-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NAOBEX_TP  
         WRITE(6,*)
         WRITE(6,*) ' alphaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)') 
     &   ' Creation      :',  (IAOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)') 
     &   ' Annihilation  :',  (IAOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' beta-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NBOBEX_TP  
         WRITE(6,*)
         WRITE(6,*) ' betaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)') 
     &   ' Creation      :',  (IBOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)') 
     &   ' Annihilation  :',  (IBOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
      END IF
*     ^ End of print is active
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Transform spin-orbit integrals only
*
      subroutine tra_so(SPININTS,HE_F,
     &                  AOTRISM_F,SOINPUT_F,
     &                  HELPMAT2,TRAINT,TRAINT_F,
     &                  COEF_F,ADJCO_F,PROD,
     &                  PROD_F,SOINT_F,
     &                  ISMIND,IGSOFF,IGSIND,
     &                  ITRALINE,IMXDIM,IMXDIM_F,MXGSDIM,
     &                  IOUTLST,ICURRDIM,ICURRDIM_F,NTRIAEL,
     &                  INTITER,ITRA_TYPE,
     &                  INTSPACE,ENVIRO,COORD,
     &                  NORB_F,
     &                  NTRSM_FULL,IFOCKEX)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
*  for lblas, COORD for x,y,z spin-orbit
      character*1 transa,transb,COORD
      character*72 ENVIRO
*
      logical SPIORB
*
      dimension SPININTS(ITRALINE,32)
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      dimension ISMIND(ICURRDIM),
     &     IGSOFF(2*NGASPI),
     &     IGSIND(ICURRDIM,NIRR_PN),
     &     IAD2(3),NORB_F(NIRR_PN),
     &     NOSYM(2),NOSYM_F(2)
*
      dimension COEF_F(ICURRDIM_F,ICURRDIM_F,2),
     &     HE_F((4*IMXDIM_F)**2),
     &     ADJCO_F(ICURRDIM_F,ICURRDIM_F,2),
     &     PROD(ICURRDIM,ICURRDIM,2),
     &     PROD_F(ICURRDIM_F,ICURRDIM_F,2),
     &     AOTRISM_F(ITRIASUM_F),
     &     HELPMAT2((2*IMXDIM)**2),
     &     TRAINT(ICURRDIM,ICURRDIM,4),
     &     TRAINT_F(ICURRDIM_F,ICURRDIM_F,4),
     &     SOINT_F(ICURRDIM_F,ICURRDIM_F),
     &     SOINPUT_F(NTRSM_FULL)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      write(*,*) '****************************'
      write(*,*) '                            '
      write(*,*) ' This is tra_so in action   '
      write(*,*) ' Using external Fock matrix '
      write(*,*) '                            '
      write(*,*) '****************************'
*
*     We are transforming integrals with two general indices
*
      if (NTEST.ge.10) then
        write(6,*) 'Coefficients in tra_so after rdcoef_m'
        do K=1,2,1
          if (K.eq.1) write(6,*) 'REAL'
          if (K.eq.2) write(6,*) 'IMAGINARY'
          call wrtmat(COEF_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                SYMDIM_F(IDIMGRP), SYMDIM_F(IDIMGRP),
     &                SYMDIM_F(IDIMGRP))
        end do
      end if
*
*  transpose coefficients ...
      IMRECO = 1
      call adjoint(ADJCO_F(1,1,1),COEF_F(1,1,1),
     &     SYMDIM_F(IDIMGRP),
     &     SYMDIM_F(IDIMGRP),IMRECO)
*
*  ... and imaginary parts of coefficients
      IMRECO = -1
      call adjoint(ADJCO_F(1,1,2),COEF_F(1,1,2),
     &     SYMDIM_F(IDIMGRP),
     &     SYMDIM_F(IDIMGRP),IMRECO)
*
* =====================================
*  Do spin-orbit integrals ...
* =====================================
*
      SPIORB = .TRUE.
      IHELEN = (4*IMXDIM_F)**2
      call intsel(SOINT_F,AOTRISM_F,SOINPUT_F,
     &            HE_F,IHELEN,
     &            NORB_F(1),SYMDIM_F(1),SYMDIM_F(IDIMGRP),
     &            IDIMGRP,NIRR_PN,
     &            ITRIA_F(1),ITRIASUM_F,
     &            ITRIADB_F(1),ITRIADBSUM_F,
     &            SPIORB,DOUGRP,NTRSM_FULL,ENVIRO,
     &            COORD)
*
*  Call transformation routine
*
      call sotra(SOINT_F,TRAINT_F,ADJCO_F,
     &           COEF_F,PROD_F,
     &           IOUTLST,
     &           SYMDIM_F(IDIMGRP),NORB_F(1),
     &           COORD,DOUGRP)
*
*     Write to appropriate output list
*
      if (COORD.eq.'X') then
        LISTIND = 6
        if (NTEST.ge.5) then
          do K=1,2,1
            write(6,*)
            write(6,'(A,A,A)') 'Writing SO-',COORD,' integrals:'
            if (K.eq.1) write(6,*) 'Real transformed ints'
            if (K.eq.2) write(6,*) 'Imaginary transformed ints'
            call wrtmat(TRAINT_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),SYMDIM_F(IDIMGRP))
          end do
        end if
      else if (COORD.eq.'Y') then
        LISTIND = 14
        if (NTEST.ge.5) then
          do K=1,2,1
            write(6,*)
            write(6,'(A,A,A)') 'Writing SO-',COORD,' integrals:'
            if (K.eq.1) write(6,*) 'Real transformed ints'
            if (K.eq.2) write(6,*) 'Imaginary transformed ints'
            call wrtmat(TRAINT_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),SYMDIM_F(IDIMGRP))
          end do
        end if
      else if (COORD.eq.'Z') then
        LISTIND = 22
        if (NTEST.ge.5) then
          do K=1,2,1
            write(6,*)
            write(6,'(A,A,A)') 'Writing SO-',COORD,' integrals:'
            if (K.eq.1) write(6,*) 'Real transformed ints'
            if (K.eq.2) write(6,*) 'Imaginary transformed ints'
            call wrtmat(TRAINT_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),SYMDIM_F(IDIMGRP))
          end do
        end if
      end if
*
      do K=1,2,1
*
*  Now pick integrals over active/virtual orbitals
        if (DOUGRP.eq.7.or.DOUGRP.eq.8) then
          NOSYM_F(1) = NORB_F(IOUTLST)
          NOSYM_F(2) = NORB_F(IOUTLST)
          NOSYM(1) = ITOTDM(IOUTLST)
          NOSYM(2) = ITOTDM(IOUTLST)
          call pick_active(2,NOSYM_F,NOSYM,ADJCO_F(1,1,1),
     &                     ADJCO_F(1,1,2),TRAINT_F(1,1,K),
     &                     TRAINT(1,1,K),1,1,1,1,
     &                     SYMDIM_F(IDIMGRP),SYMDIM(IDIMGRP),
     &                     2,2,NTEST)
        else
          ISYMOFF = 1
          if (DOUGRP.eq.4.and.IOUTLST.ge.3) ISYMOFF = 3
          call pick_active(2,NORB_F(ISYMOFF),ITOTDM(ISYMOFF),
     &                     ADJCO_F(1,1,1),
     &                     ADJCO_F(1,1,2),TRAINT_F(1,1,K),
     &                     TRAINT(1,1,K),1,1,1,1,
     &                     SYMDIM_F(IDIMGRP),SYMDIM(IDIMGRP),
     &                     2,2,NTEST)
        end if
*
        if (NTEST.ge.1) then
          write(6,*) 'Spin-orbit integrals after pick_active:'
          if (K.eq.1) then
           write(6,'(A5,I3,A11,A2,A5)')
     &                'List ',IOUTLST,',coordinate ',COORD,', real:'
           call wrtmat(TRAINT(1,1,1),SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &                               SYMDIM(IDIMGRP),SYMDIM(IDIMGRP))
          else if (K.eq.2) then
           write(6,'(A5,I3,A11,A2,A5)')
     &                'List ',IOUTLST,',coordinate ',COORD,',imag.:'
           call wrtmat(TRAINT(1,1,2),SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &                               SYMDIM(IDIMGRP),SYMDIM(IDIMGRP))
          end if
        end if
*
        INTLIST = 2*IOUTLST + K + LISTIND
        if (NTEST.ge.3) then
          if (K.eq.1) then
            write(*,*) 'Real part of Integrals of list  '
            write(*,*) IOUTLST,'INTLIST is   ',INTLIST
            write(6,*) 'coordinate ',COORD
            write(*,*) 'are being passed to intgas'
          else
            write(*,*) 'Im. part of Integrals of list  '
            write(*,*) IOUTLST,'INTLIST is   ',INTLIST
            write(6,*) 'coordinate ',COORD
            write(*,*) 'are being passed to intgas'
          end if
        end if
        if (ITRA_TYPE.le.3) then
          IACSYM = 1
          if (IOUTLST.gt.2) IACSYM = 3
          call intgas3(SPININTS(1,INTLIST),
     &                 TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                 INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                 ISMIND,IGSOFF,IGSIND,IACSYM,
     &                 IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
        else
*  Store for MCSCF run ...
          if (K.eq.1) LSM = 1
          if (K.eq.2) LSM = -1
          call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &              SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
        end if
      end do
*
* Check transformed active integral matrix for
* hermiticity/antihermiticity:
      do K=1,2,1
        if (K.eq.1) KF = 1
        if (K.eq.2) KF = -1
        XMXDIFF = 1.0D-11
        DO I = 1,SYMDIM(IDIMGRP)
          DO J = 1, I
            if (NTEST.ge.20) then
              write(6,'(A,1F18.10)')
     &          'DIFF = ',ABS(TRAINT(I,J,K) - KF * TRAINT(J,I,K))
              write(6,'(5X,2F18.10)') TRAINT(I,J,K),TRAINT(J,I,K)
            end if
            IF(ABS(TRAINT(I,J,K) - KF * TRAINT(J,I,K)).
     &             GT.XMXDIFF) THEN
              write(6,*) 'Transformed SO ints not (anti)hermitian.'
              if (K.eq.1) then
                write(6,*) 'Real integral I,J,TRAINT(I,J) '
                write(6,*) I,J,TRAINT(I,J,1)
                write(6,*) 'does not match'
                write(6,*) 'Real integral J,I,TRAINT(J,I) '
                write(6,*) J,I,TRAINT(J,I,1)
                stop 'Problem in tra_so.'
              else if (K.eq.2) then
                write(6,*) 'Imag integral I,J,TRAINT(I,J) '
                write(6,*) I,J,TRAINT(I,J,1)
                write(6,*) 'does not match'
                write(6,*) 'Imag integral J,I,TRAINT(J,I) '
                write(6,*) J,I,-TRAINT(J,I,1)
                stop 'Problem in tra_so.'
              end if
            END IF
          END DO
        END DO
      end do
*
C     if (IOUTLST.eq.2) stop 'Timo forced me to stop in tra_so.'
*
 1000 return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine tweltra(I12O,I34O,I1234O,ICCSMO,
     &                   IRIO,ISPITRA_RC,ISPITRA,NAOFSM,
     &                   LSPINOR,ITRA_TYPE,ITER)
*
* Test general spinor transformation program
*
* Jeppe + Timo, April-July 1997
*
* ITRA_TYPE Defines the type of transformation
* ITRA_TYPE = 1 : Normal nonrelativistic transformation
* ITRA_TYPE = 2 : 2-component relativistic transformation for CI,
*                 all nonredundant integral lists
* ITRA_TYPE = 3 : 2-component relativistic transformation for CI,
*                 Excluding integrals with an odd number of bars
* ITRA_TYPE = 4 : 2-component relativistic transformation for MCSCF,
*                 all nonredundant integral lists
* ITRA_TYPE = 5 : 2-component relativistic transformation for MCSCF,
*                 Excluding integrals with an odd number of bars
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
*
*. output lists of integrals
      dimension I12O(*),I34O(*),I1234O(*)
      dimension ICCSMO(*),IRIO(MXBLKO)
      dimension ISPITRA(4,*),ISPITRA_RC(*)
*
*. size of spinor expansions
      dimension LSPINOR(N_SPINR_EXP)
*
#include "clunit.inc"
*
777   format(1X,4I3)
888   format(1X,2I3)
999   format(1X,8I3)
700   format(6I3)
*
*
      WRITE(6,*) ' ***************************** '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' *     Welcome to TWELTRA    * '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' ***************************** '
*
*  Test
*
      if (DOUGRP.eq.4) then
         write(*,*) 'Point group mult. matrix'
         do I=1,NIRR_PN,1
            write(6,'(10X,4I3)') (IPGMULT(I,J),J=1,NIRR_PN,1)
         end do
         write(*,*) 'Double group mult. matrix'
         do I=1,NIRR_DG,1
            write(6,'(10X,8I3)') (IDBGMULT(I,J),J=1,NIRR_DG,1)
         end do
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
         do I=1,NIRR_PN,1
            write(6,'(10X,2I3)') (IPGMULT(I,J),J=1,NIRR_PN,1)
         end do
         write(*,*) 'Double group mult. matrix'
         do I=1,NIRR_DG,1
            write(6,'(10X,4I3)') (IDBGMULT(I,J),J=1,NIRR_DG,1)
         end do
      end if
*
*. Normal list of input integrals
      I12 = 1
      I34 = 1
      I1234 = 1
      ICCSM = 1
      IOPSMPN = 1
      IOPSMDG = 1
*
*  Read input from file 40:
      read(LU40,'(A)')
      read(LU40,*) NLISTO
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,700) I12O(I)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) I34O(I)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) I1234O(I)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) ICCSMO(I)
      end do
*
      read(LU40,'(A)')
      do I=1,4,1
         read(LU40,700) (ISPITRA(I,J),J=1,NLISTO,1)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) ISPITRA_RC(I)
      end do
*
      write(6,*) ' tweltra: ITRA_TYPE = ' ,ITRA_TYPE
      if (ITRA_TYPE.eq.1) then
*  Good old non-relativistic integral transformation
*. A normal integral transformation, spinors have just
* pn symmetry, and symmetry effects of spin functions
* are neglected
         NSPINOR = 1
C         NLISTO = 1
C         I12O(1) = 1
C         I34O(1) = 1
C         I1234O(1) = 1
C         ICCSMO(1) = 1
C         ISPITRA(1,1) = 1
C         ISPITRA(2,1) = 1
C         ISPITRA(3,1) = 1
C         ISPITRA(4,1) = 1
C         ISPITRA_RC(1) = 1
      else if (ITRA_TYPE.eq.2) then
         NSPINOR = 2
*  Relativistic transformation to double group spinors,
*  CI with all classes of integrals included.
*  Unbarred trans is 1, bar trans is 2.
*
* First list : (K L_ | M N_) where _ = bar
* Second list : (K L_ | M N)
* Third list : (K L | M N)
* Fourth list : (K_ L | M N_)
* Fifth list : (K_ L | M N)
* Sixth list : (K_ L | M_ N)
*
      else if(ITRA_TYPE.EQ.3) then
*. relativistic transformation to double group spinors,
*. CI with integrals including an even number of bars
         NSPINOR = 2
*
* First list : (K L_ | M N_) Where _ = bar
* second list : (K L | M N)
* Third list : (K_ L | M N_)
* Fourth list : (K_ L | M_ N)
      else if (ITRA_TYPE.eq.4) then
         write(6,*) ' Add types for MCSCF transformation'
      ELSE IF (ITRA_TYPE.EQ.5) THEN
         WRITE(6,*) ' Add types for MCSCF transformation'
      END IF
*
*. Relative size of spinor expansions.
*  As an example, in the case of C2(*) there are two different spinor
*  expansions, one for alpha and one for beta spinors. The size of the
*  expansions of two particular spatial symmetries might be different, so that
*  there is an overall difference of the expansion sizes of two indices in
*  question. LSPINOR assumes the value 1 for a given expansion and another
*  1 if the expansion compared to it has the same size, a value of 2 if the
*  second expansion is smaller than the first. The expansion compared to
*  each other are characterized by the value of ISPITRA(index,list), which
*  in this case is either 1 or two, i.e. there are 2 expansions, and in this
*  case they are of equal size.
*
      DO JSPINOR = 1,NSPINOR,1
         LSPINOR(JSPINOR) = 1
      END DO

      CALL SPITRA2(I12O,I34O,I1234O,
     &             ICCSMO,ISPITRA,ISPITRA_RC,
     &             LSPINOR,IRIO,ICCSM,
     &             NLISTO,I1234,I12,I34,IOPSMPN,IOPSMDG)
*
      WRITE(6,*) ' ***************************** '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' *     TWELTRA says goodbye  * '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' ***************************** '
*
      return
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WRT_SPEXP(CSPSAO,ICSPSAO,NAOFSM,NBSPFSM,NUSPFSM)
*
* Write Spinor -SAO expansion CSPSAO
*
* Jeppe + Timo, December 1997
*
      IMPLICIT REAL*8  (A-H,O-Z)
#include "ipoist8.inc"
*.General input :  Double group - point group translation :
#include "mxpdim.inc"
#include "symm.inc"
*. Specific Input
      DIMENSION CSPSAO(*)
      INTEGER   ICSPSAO(2,2),NAOFSM(*),NBSPFSM(*),NUSPFSM(*)
*
      write(6,*) ' In WRT_SPEXP: '
      WRITE(6,*) ' NAOFSM : ', NAOFSM(1), NAOFSM(2)
      write(6,*) ' NBSPFSM: ',NBSPFSM(1),NBSPFSM(2)
      write(6,*) ' NUSPFSM: ',NUSPFSM(1),NUSPFSM(2)
      IOFF = 1
      DO IBU = 1,2
      DO IRI = 1,2
        IF(IBU.EQ.1) THEN
          WRITE(6,*)
          IF(IRI.EQ.1) THEN
            WRITE(6,*) ' ==================================='
            WRITE(6,*) ' Barred spinor expansion, real part '
            WRITE(6,*) ' ==================================='
          ELSE
            WRITE(6,*) ' ========================================='
            WRITE(6,*) ' Barred spinor expansion, imaginary  part '
            WRITE(6,*) ' ========================================='
          END IF
        ELSE IF(IBU.EQ.2) THEN
          IF(IRI.EQ.1) THEN
            WRITE(6,*) ' ====================================='
            WRITE(6,*) ' Unbarred spinor expansion, real part '
            WRITE(6,*) ' ====================================='
          ELSE
            WRITE(6,*) ' ===================================='
            WRITE(6,*) ' Unbarred spinor expansion, imag part '
            WRITE(6,*) ' ===================================='
          END IF
          WRITE(6,*)
        END IF
        DO ISPSM = NIRR_DG/2+1, NIRR_DG
            WRITE(6,'(1H ,A,I3)') '   Spinors of symmetry ', ISPSM
            WRITE(6,'(1H , A)'  ) '   **************************'
            WRITE(6,*)
*. ao * alpha
          IAOSM_A = ISPSF_MO(ISPSM,1)
*. ao * beta
          IAOSM_B = ISPSF_MO(ISPSM,2)
          IF(IAOSM_A.EQ.0) THEN
            Write(6,*) ' (No alpha component ) '
            LAO_A = 0
          ELSE
            WRITE(6,*) ' Alpha * AO''s with symmetry ', IAOSM_A
            LAO_A = NAOFSM(IAOSM_A)
          END IF
          IF(IAOSM_B.EQ.0) THEN
            Write(6,*) ' (No beta  component ) '
            LAO_B = 0
          ELSE
            WRITE(6,*) ' Beta  * AO''s with symmetry ', IAOSM_B
            LAO_B = NAOFSM(IAOSM_B)
          END IF
C         WRITE(6,*) ' IAOSM_A, IAOSM_B ', IAOSM_A, IAOSM_B
*
          LAO_AB = LAO_A + LAO_B
          IF(IBU.EQ.1) THEN
            LSPINOR = NBSPFSM(ISPSM)
          ELSE
            LSPINOR = NUSPFSM(ISPSM)
          END IF
*
          CALL WRTMAT(CSPSAO(IOFF),LAO_AB,LSPINOR,LAO_AB,LSPINOR)
          IOFF = IOFF + LSPINOR*LAO_AB
        END DO
*       ^ End of loop over spinor symmetries
      END DO
      END DO
*     ^ End of loop over barred/unbarred, real/imaginary parts
      RETURN
      END
