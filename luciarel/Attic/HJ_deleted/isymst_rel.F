      FUNCTION ISYMST_REL(STRING,NEL)
*
* Master routine for symmetry of string
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input ( PNTGRP is used )
C     INTEGER PNTGRP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C
#include "mxpdim.inc"
#include "lucinp.inc"
*. Specific input
      INTEGER STRING(*)
      IF(PNTGRP.EQ.1) THEN
*.D2h
        ISYMST_REL = ISYMS1_REL(STRING,NEL)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
*.Cinfv Dinfh O3
        ISYMST_REL = ISYMS2_REL(STRING,NEL)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        ISYMST_REL = ISYMSDG(STRING,NEL)
      ELSE
        WRITE(6,*) ' Sorry PNTGRP option not programmed ', PNTGRP
        WRITE(6,*) ' Enforced stop in ISYMST '
        Call Abend1( 5 )
      END IF
*
      RETURN
      END
