
      SUBROUTINE INTDIM(IPRNT)
*
* Number of integrals and storage mode
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "crun.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
* =======
*. Output
* =======
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*.1 : Number of one-electron integrals
      NINT1 =  NSXFSM(NSMOB,MXPOBS,NTOOBS,NTOOBS,ITSSX,ADSXA,1,IPRNT)
*.2 : Number of two-electron integrals
C     IF(INCORE.EQ.1.OR.EXTSPC.EQ.0) THEN
        IF(PNTGRP.EQ.1) THEN
*. Full eightfold symmetry can be used
          I12S = 1
          I34S = 1
          I1234S = 1
        ELSE
*. Only symmetry between 12 and 34
          I12S = 0
          I34S = 0
          I1234S = 1
        END IF
        NINT2 = NDXFSM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,
     &                  NTOOBS,ITSDX,ADSXA,SXDXSX,I12S,I34S,I1234S,
     &                  IPRNT )
C     END IF
*. Number of symmetry blocks of one- and two-electron integrals
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
Ctimo
*  Write number of integrals
      write(6,*) '# 1-particle integrals is   ',NINT1
      write(6,*) '# 2-particle integrals is   ',NINT2
      RETURN
      END
