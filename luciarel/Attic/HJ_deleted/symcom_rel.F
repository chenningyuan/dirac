      SUBROUTINE SYMCOM_REL(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1 I2 = I12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I2
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* IOBJ = 1 : I1,I2 are strings I12 determinant
* ( Other things can follow )
* IOBJ = 2 : I1,I2,I3 are externals
* IOBJ = 3 : I1 is an external, I2,I3 are dets
* IOBJ = 4 : I1 is orbital, I2 is string,l, I12 is string
* IOBJ = 5 : I1 is single excitation, I2 is string,l, I12 is string
* IOBJ = 6 : I1 is orbital, I2 is Orbital I12 is single excitation
*
* If obtained symmetry I1 or I2 is outside bounds,
* zero is returned.
*
* Jeppe Olsen , Spring of 1991
*
* ================
*. Driver routine
* ================
*.LUCINP (PNTGRP is used )
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
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "mxpdim.inc"
#include "lucinp.inc"
*
      NTEST = 0000
*
      if (NTEST.ge.2000) then
        write(6,*)
        write(6,*) 'symcom test:'
        write(6,*) 'ITASK : ',ITASK
        write(6,*) 'IOBJ  : ',IOBJ
        write(6,*) 'I1    : ',I1
        write(6,*) 'I2    : ',I2
        write(6,*) 'I12   : ',I12
      end if
*
      IF(PNTGRP.EQ.1) THEN
        CALL SYMCM1(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
        CALL SYMCM2(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        CALL SYMCMDG(ITASK,IOBJ,I1,I2,I12)
      ELSE
        WRITE(6,*) ' PNTGRP parameter out of bounds ', PNTGRP
        WRITE(6,*) ' Enforced stop in SYMCOM '
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
