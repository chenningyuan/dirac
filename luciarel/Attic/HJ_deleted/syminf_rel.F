      SUBROUTINE SYMINF_REL(IPRNT)
*
* Information about number of symmetries
*
* Input : /LUCINP/,/ORBINP
* Output : /CSM/,/CSMPRO/
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     INTEGER PNTGRP,EXTSPC
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

*. Output
* NSMSX : number of symmetries of single excitations
* NSMDX : Number of symmetries of double excitations
* NSMST : Number of symmetries of strings
* NSMCI : NUmber of symmetries of CI spaces
* ITSSX : Total symmetrix single excitation
* ITSDX : Total symmetrix double excitation
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
* ADASX : symmetry of orbs i and i => symmetry of a+iaj
* ASXAD : symmetry of orb j and excit a+iaj => symmetry of i
* ADSXA : symmetry of orb i and excit a+iaj => symmetry of j
*
* SXSXDX : Symmetry of two single excitations
*          => symmetry of double  excitation
* SXDXSX : Symmetry of single excitation and double excitation
*          => symmetry of single  excitation

*.
      IF(PNTGRP.EQ.1) THEN
* =====
* D 2 h
* =====
        CALL ZSYM1_REL(NIRREP,IPRNT)
      ELSE IF(PNTGRP.EQ.2) THEN
* ========
* C inf V
* ========
C            ZNONAB_REL(INVCEN,MAXMLO,NSMOB)
        CALL ZNONAB_REL(0,MAXML,NSMOB,IPRNT)
        CALL ZSYM2_REL(IPRNT)
      ELSE IF(PNTGRP.EQ.3.OR.PNTGRP.EQ.4) THEN
* ===========
* D inf H O3
* ===========
C            ZNONAB_REL(INVCEN,MAXMLO,NSMOB)
        CALL ZNONAB_REL(1,MAXML,NSMOB,IPRNT)
        CALL ZSYM2_REL(IPRNT)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
* ===========
* Double group
* ============
*
        CALL ZSYMDG(IPRNT)
      ELSE
        WRITE(6,*) ' You are too early , sorry '
        WRITE(6,*) ' Illegal PNTGRP in SYMINF ',PNTGRP
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
