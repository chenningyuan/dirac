      SUBROUTINE INTPNT_DBG
*
* Pointers to symmetry blocks of integrals for DOUBLE GROUP
*
* Jeppe Olsen, Jan 98

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
#include "cintfo.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
C     COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2,
C    &              NL1D,IUB1D(2,10),L1D(10),
C    &              NL2D,IUB2D(4,10),ISM2D(3,10),L2D(10),
C    &              LMAXE2LST
*
* =======
*. Output
* =======
*
*. Local scratch
      DIMENSION NNORB(MXNDGIRR,4)
*
*.1 : Offset for one-electron integrals
      DO IL1 = 1, NL1D
        DO INDEX = 1,2,1
         IF(IUB1D(INDEX,IL1).EQ.1) THEN
           CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
         ELSE
           CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
         END IF
        END DO
*
        CALL PNT2DM2(1,NSMOB,NNORB(1,1),NNORB(1,2),ITSSX,
     &       WORK(KPINT1D(IL1)) )
C            PNT2DM2(I12SM,NSMOB,IPSM,JPSM,IJSM,IPNTR)
      END DO
*
*.2 : two-electron integrals
      DO IL2 = 1, NL2D
        DO INDEX = 1,4,1
         IF(IUB2D(INDEX,IL2).EQ.1) THEN
           CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
         ELSE
           CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
         END IF
        END DO
*
        CALL PNT4DM2(NSMOB,
     &       NNORB(1,1),NNORB(1,2),NNORB(1,3),NNORB(1,4),
     &       ITSDX,ISM2D(1,IL2),ISM2D(2,IL2),ISM2D(3,IL2),
     &       WORK(KPINT2D(IL2)))
      END DO
*
      CALL LMEMCHK
      RETURN
      END
