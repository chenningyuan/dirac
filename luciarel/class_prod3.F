      SUBROUTINE CLASS_PROD3(VEC1,VEC2,IBLOCK_OFF,NBLOCK,IBLOCK,
     &                       IBLTOCLS,NCLS,CLSVEC)
*
* Two vectors in blocked form are given.
* Find contributions of product to each occupation class
*
* Current version (PROD3 !!) uses IBLTOCLS to give relation
* between blocks and classes
      IMPLICIT REAL*8(A-H,O-Z)
*
      REAL*8 INPROD
*. Input
      DIMENSION VEC1(*),VEC2(*)
      INTEGER IBLOCK(8,*), IBLTOCLS(*)
*. Input/output
      DIMENSION CLSVEC(*)

*
      IOFF = 1
      DO JBLOCK = IBLOCK_OFF,IBLOCK_OFF-1+NBLOCK
      IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        JATP = IBLOCK(1,JBLOCK)
        JBTP = IBLOCK(2,JBLOCK)
        NELMNT = IBLOCK(8,JBLOCK)
        JCLS = IBLTOCLS(JBLOCK)
        XTERM = INPROD(VEC1(IOFF),VEC2(IOFF),NELMNT)
C?      WRITE(6,*)
C?   &  ' CLASS_PROD : CLASS and XTERM = ', JCLS,XTERM
          CLSVEC(JCLS) = CLSVEC(JCLS) + XTERM
COLD    END IF
        IOFF = IOFF + NELMNT
      END IF
      END DO
*
      NTEST = 0
      IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Updated CLSVEC '
         CALL WRTMAT(CLSVEC,1,NCLS,1,NCLS)
      END IF
*
      RETURN
      END
