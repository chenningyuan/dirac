      SUBROUTINE DIM_T1T2_TO_T12_MAP_REL(LEN_T1T2_STRING,LEN_T1T2_TCC,
     &                                   WORK,KFREE,LFREE)
*
* Dimension of T1*T2 => T12 mappings of strings
*
* LEN_T1T2_STRING : Max dimension of MAT(T1,T2) where T1 and 
*                   T2 are spin strings, and T1*T2 is allowed 
*                   spin strings
*
* Jeppe Olsen, May 2000 
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
!#include "wrkspc.inc"
#include "cgas.inc"
*
      DIMENSION WORK(*)
*
      CALL DIM_T1T2_TO_T12_S_REL(LEN_T1T2_STRING,LEN_T1T2_TCC,
     &                        NSPOBEX_TP_CC,WORK(KLSOBEX_CC),
     &                        NGAS)
      RETURN
      END
