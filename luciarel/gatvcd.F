
      SUBROUTINE GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,IPRT)
*
* Gather elements from a file LU
*
* XGAT(I) = Vector(IGAT(I))
*
* Jeppe Olsen, September 1993
*
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input
      INTEGER IGAT(NGAT)
*. Output
      DIMENSION XGAT(NGAT)
*. Scratch
      DIMENSION SEGMNT(*)
*
      CALL REWINE(LU,-1)
*
      IBASE = 1
      IBLOCK = 0
*
*. Loop over blocks of file
*
 1000 CONTINUE
        IBLOCK = IBLOCK + 1
        CALL NEXREC(LU,LBLK,SEGMNT,IEND,LENGTH)
        IF(IPRT.GE.10)
     &  WRITE(6,*) LENGTH, ' elements in block ',IBLOCK
        IF(IEND.EQ.0) THEN
          IFIRST = IBASE
          ILAST = IBASE + LENGTH - 1
          DO 100 I = 1, NGAT
            IF(IFIRST .LE. IGAT(I) .AND. IGAT(I) .LE. ILAST )
     &      XGAT(I) = SEGMNT(IGAT(I)-IFIRST+1)
C?          IF(IFIRST .LE. IGAT(I) .AND. IGAT(I) .LE. ILAST )
C?   &      write(6,*) ' Catch I IGAT(I) XGAT(I) ',
C?   &                         I,IGAT(I),XGAT(I)
  100     CONTINUE
          IBASE = IBASE + LENGTH
      IF(LBLK.LT.0) GOTO 1000
        END IF
*
      NTEST = 0
      NTEST = MAX(IPRT,NTEST)
      IF(NTEST.GE.5) THEN
       WRITE(6,*) ' Gathered vector from GATVCD '
       CALL WRTMAT(XGAT,1,NGAT,1,NGAT)
      END IF
*
      RETURN
      END
