      SUBROUTINE GET_BLOCKS_FROM_DISC_REL
     &           (LU,NBLOCK,IOFF,IBLOCK,NBLOCKT,C,IREW)
*
* Obtain blocks IOFF - IOFF + NBLOCK from file LU
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(*)
      DIMENSION IBLOCK(8,*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "csm.inc"
#include "stinf.inc"
*
*. Local Scratch
*. Should be uncommented if NEW_OR_OLD is set to 2
*
C     PARAMETER (MXPNBLK = 33000)
C     DIMENSION ISCR(3*MXPNBLK)
*
C?    WRITE(6,*) ' ENTERING GET_BLOCKS_FROM_DISC'
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
CM    IREW = 1
      ISCAL = 0
      CALL GET_TTS_BATCHN(C,NBLOCK,IBLOCK(1,IOFF),NBLOCKT,IBLOCK,
     &                     NOCTPA,NOCTPB,NSMST,
     &                     WORK(KNSTSO(IATP)), WORK(KNSTSO2(IBTP)),
     &                     IDC,LU,IREW,ISCAL)

*
      RETURN
      END
