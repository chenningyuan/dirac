      SUBROUTINE INUM_FOR_OCC2_REL(IOCC,INUM,NDIFF,
     &                             WORK,KFREE,LFREE)
*
* An operator is defined by OCC(NGAS,4). 
* Obtain the type number of this operator in list of CC operators
* A -1 indicates a nontrivial excitation not included in the list
* NDIFF is the number of smallest number of operators needed 
* to bring IOCC into the space 
*
* Difference to INUM_FOR_OCC2 : NDIFF added
*
* Jeppe Olsen, November 2000
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "ctcc.inc"
#include "ctccp.inc"
*
      DIMENSION WORK(*)
*
      CALL INUM_FOR_OCC21(IOCC,INUM,NGAS,WORK(KLSOBEX_CC),
     &                   NSPOBEX_TP_CC,NDIFF)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Coupled cluster operator '
        CALL WRT_SPOX_TP_CC_KRCC(IOCC,1)
        WRITE(6,*) ' Corresponding type number ', INUM 
        WRITE(6,*) ' NDIFF = ', NDIFF
      END IF
*
      RETURN
      END
*
      SUBROUTINE INUM_FOR_OCC21(IOCC,INUM,NGAS,ISPOBEX_TP,
     &                         NSPOBEX_TP,NDIFF)
*
* Type number of CAAB operator
*. ISPOBEX_TP is assumed to be extended list also containing 
*. zero-particle unit operator
*
* Jeppe Olsen, May 2000
*
#include "implicit.inc"
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP+1)
      INTEGER IOCC(4*NGAS)
*
      INUM = -1
      DO ITP = 1, NSPOBEX_TP+1
        IDIFF = 0
        DO I = 1, 4*NGAS
          IDIFF = IDIFF + ABS(IOCC(I)-ISPOBEX_TP(I,ITP))
        END DO
        IF(ITP.EQ.1) THEN
          NDIFF = IDIFF
        ELSE
          NDIFF = MIN(NDIFF,IDIFF)
        END IF
        IF(IDIFF.EQ.0) INUM = ITP
      END DO
*
      RETURN
      END

