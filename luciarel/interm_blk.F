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
