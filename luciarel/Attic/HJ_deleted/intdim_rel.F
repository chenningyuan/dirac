      SUBROUTINE INTDIM_REL(IPRNT)
*
* Number of integrals and storage mode
*
*. Last modifications : Jan 98, relativistic integrals added
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
#include "cgas.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),
     &              SXDXSX(2*MXPOBS,4*MXPOBS)
*. Local arrays
      INTEGER NNORB(MXNDGIRR,4)
*
* =======
*. Output
* =======
*. Atmost 10 double group integral lists
C     COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2,
C    &              NL1D,IUB1D(2,10),L1D(10),
C    &              NL2D,IUB2D(4,10),ISM2D(3,10),L2D(10)
#include "cintfo.inc"
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.2)
     &    WRITE(6,'(A25,I3)') ' INTDIM speaking, IDBG = ', IDBG
*
* =================================
* integrals in point group symmetry
* =================================
*
      IPNTGRP_SAVE = PNTGRP
*. Force the use of point group symmetry
      IF(IDBG.EQ.1) PNTGRP = 1
*.1 : Number of one-electron integrals
C     NINT1 =  NSXFSM(NSMOB,MXPOBS,NTOOBS,NTOOBS,ITSSX,ADSXA,1,IPRNT)
      NINT1 =  NSXFSM2(NSMSH,NGSSHT,NGSSHT,ITSSX,1,IPRNT)
C              NSXFSM2(NSMOB,NO1PS,NO2PS,ISXSM,ISYM,IPRNT)
*.2 : Number of two-electron integrals
*. Full eightfold symmetry can be used
      I12S = 1
      I34S = 1
      I1234S = 1
      NINT2 = NDXFSM2(NSMSH,NGSSHT,NGSSHT,NGSSHT,NGSSHT,
     &                I12S,I34S,I1234S,IPRNT )
C     NDXFSM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
C    &         IS12,IS34,IS1234,IPRNT)
*
* =======================
* Double group integrals
* =======================
*
      IF(IDBG.EQ.0) THEN
        NL1D = 0
        NL2D = 0
      ELSE
        PNTGRP = IPNTGRP_SAVE
*. Number of double group integral lists
        CALL DBG_INT_LISTS(NL1D,IUB1D,NL2D,IUB2D,ISM2D,IDENSI,
     &                     ISPINFREE)
C           DBG_INT_LISTS(NL1,IUB1,NL2,IUB2,ISM2,IDENSI,
C    &                    ISO1,ISO2)
C    &               NL1D,IUB1D(2,10),L1D(10),
C    &               NL2D,IUB2D(4,10),ISM2D(3,10),L2D(10)
*
        DO IL1 = 1, NL1D
          DO INDEX = 1, 2
            IF(IUB1D(INDEX,IL1).EQ.1) THEN
              CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
            ELSE
              CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
            END IF
          END DO
          L1D(IL1) = NSXFSM2(NSMOB,NNORB(1,1),NNORB(1,2),ITSSX,1,IPRNT)
        END DO
*
        DO IL2 = 1, NL2D
          DO INDEX = 1, 4
            IF(IUB2D(INDEX,IL2).EQ.1) THEN
              CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
            ELSE
              CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
            END IF
          END DO
          L2D(IL2) = NDXFSM2(NSMOB,
     &               NNORB(1,1),NNORB(1,2),NNORB(1,3),NNORB(1,4),
     &               ISM2D(1,IL2),ISM2D(2,IL2),ISM2D(3,IL2),
     &               IPRNT )
        END DO
      END IF
*. Largest block of two-electron integrals
      LMAXE2LST = NINT2
      DO IL2 = 1, NL2D
        LMAXE2LST = MAX(LMAXE2LST,L2D(IL2))
      END DO
*
      write(6,*)
      write(6,'(A46,I10)')
     &      'The largest block of two-electron integrals is',LMAXE2LST
*
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
*
      RETURN
      END
