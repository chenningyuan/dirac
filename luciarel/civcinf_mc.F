*  Fill common block with communication information
*  KRMC - LUCIAREL.
*
*  Timo Fleig        Mar 17, 2002
*
      subroutine civcinf_mc(IBATS,LBATS,I1BATS,LEBATS,NBAT_MS2,
     &                      NMS2VAL,NBLOCK)
*
#include "implicit.inc"
#include "mxpdim.inc"
*
#include <civecinf.h>
*
      dimension IBATS(8,NBLOCK),LBATS(NBLOCK),I1BATS(NBLOCK),
     &          LEBATS(NBLOCK),NBAT_MS2(NMS2VAL)
*
      do ITP = 1,NBLOCK,1
        LBATS_MC(ITP) = LBATS(ITP)
        I1BATS_MC(ITP) = I1BATS(ITP)
        LEBATS_MC(ITP) = LEBATS(ITP)
        do JTP = 1,8,1
          IBATS_MC(JTP,ITP) = IBATS(JTP,ITP)
        end do
      end do
*
      NMK2VAL = NMS2VAL
      do IM = 1,NMS2VAL,1
        NBAT_MS2_MC(IM) = NBAT_MS2(IM)
      end do
*
      return
      end
