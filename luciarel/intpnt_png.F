      SUBROUTINE INTPNT_PNG(IPNT1,IPNT2)
*
* Pointers to symmetry blocks of integrals for POINT GROUP
* IPNT1 : Pointer to given one-electron block
* IPNT2 : Pointer to given two-electron block
*
* Note : Generation of symmetry of last index (ISL1, ISL2) eliminated

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
#include "cgas.inc"
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
      INTEGER IPNT1(NSMOB)
      INTEGER IPNT2(NSMOB,NSMOB,NSMOB)
*
* if the actual group is a double group, hide this for the moment
* and use the corresponding pointgroup
      IPNTGRP_SAVE = PNTGRP
      IF(IDBG.EQ.1) THEN
        PNTGRP = 1
      END IF
*
*.1 : Number of one-electron integrals
      CALL PNT2DM2(1,NSMSH,NGSSHT,NGSSHT,ITSSX,IPNT1)
C          PNT2DM2(I12SM,NSMOB,IPSM,JPSM,IJSM,IPNTR)

*.2 : two-electron integrals
      I12S = 1
      I34S = 1
      I1234S = 1
      CALL PNT4DM2(NSMSH,NGSSHT,NGSSHT,NGSSHT,NGSSHT,
     &            ITSDX,I12S,I34S,I1234S,IPNT2)
C          PNT4DM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
C    &           IDXSM,IS12,IS34,IS1234,IPNTR)

      PNTGRP = IPNTGRP_SAVE
*
      CALL LMEMCHK
      RETURN
      END
