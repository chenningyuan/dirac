      subroutine sigden_ctrl(C,HC,LUC,LUHC,T,ISIGDEN)
*
* Outer routine for sigma vector generation
* in SIGDEN formalism
*
* Modified mv7 routine 
*  Timo Fleig, Halloween 2001  =:-()
*
#include "implicit.inc"
#include "mxpdim.inc"
*
* ======
*.Output
* ======
*
#include "wrkspc.inc"
#include "cstate.inc"
*
      CALL QENTER('SIGDC')
      IDUM = 0
*
*. Arrays for partitioning of sigma
* (*_MS2 arrays are in CSTATE )
      JCMBSPC = 1
      call z_blkfo(IDC,NMS2VAL,JCMBSPC,IREFSM,
     &             KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,KLSBLTP,
     &             NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2)
*
      call sigden_rel(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                WORK(KLSI1BT),WORK(KLSIBT),LUC,LUHC,T,ISIGDEN)
*
      CALL QEXIT('SIGDC')
*
      RETURN
      END
