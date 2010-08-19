      SUBROUTINE WRT_CC_VEC2(CC,LU,CCTYPE)
*
* Print vector of CC amplitudes, type is governed by CCTYPE
*     
#include "implicit.inc"
*. General input 
#include "mxpdim.inc"
#include "crun.inc"
#include "ctcc.inc"
#include "wrkspc.inc"
C     COMMON/CTCC/KLSOBEX,NSPOBEX_TP,KLLSOBEX,KLIBSOBEX,LEN_T_VEC,
C    &             MX_ST_TSOSO,MX_ST_TSOSO_BLK

      CHARACTER*6 CCTYPE
      DIMENSION CC(*) 
*     
      IF(CCTYPE(1:2).EQ.'CC') THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Single excitation part '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        CALL WRT_SX(CC,1)
*
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Double excitation part '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        CALL WRT_DX1(CC(1+NSXE),1)
      ELSE IF (CCTYPE(1:6).EQ.'GEN_CC') THEN
        CALL WRTBLKN_CC(CC,NSPOBEX_TP,WORK(KLLSOBEX))
      END IF
*
      RETURN
      END
*
      SUBROUTINE WRT_CC_VEC_REL(CC,LU)
*
* Print vector of CC amplitudes 
* Timo Fleig, Winter 2009/2010
*
#include "implicit.inc"
*. Generel input 
#include "mxpdim.inc"
#include "crun.inc"
      DIMENSION CC(*)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Single excitation part '
      WRITE(6,*) ' CAUTION : This is not yet adapted to the rel. case !'
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_SX_REL(CC,1)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Double excitation part '
      WRITE(6,*) ' CAUTION : This is not yet adapted to the rel. case !'
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_DX_REL(CC(1+NSXE),1)
*
      RETURN
      END
*
