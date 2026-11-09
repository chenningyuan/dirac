!-----------------------------------------
! Module with variables previously defined
! at blocks.h include file (COMMON /BLOCKS/)
!-----------------------------------------
   module blocks_mod
   use maxorb_mod, only: MXCORB, MXSHEL
   use aovec_mod, only: MXAOVC
   implicit none
   real*8 :: CENTSH(MXSHEL,3)
   integer :: MAXSHL, NLRGBL, NSMLBL, NHKTSH(MXSHEL), KHKTSH(MXSHEL), KCKTSH(MXSHEL), &
              ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL), NSTRSH(MXSHEL), &
              NCNTSH(MXSHEL), NSETSH(MXSHEL,2), JSTRSH(MXSHEL,MXAOVC,2), &
              NPRIMS(MXSHEL,MXAOVC,2), NCONTS(MXSHEL,MXAOVC,2), IORBSH(MXSHEL,MXAOVC), &
              IORBSB(0:MXCORB-1), NRCSH(MXSHEL), LCLASH(MXSHEL), NLRBL, ISYMBL(MXSHEL,8), &
              NSYMBL
   logical :: BIGVEC, SEGMEN, SEGMSH(MXSHEL), SPHRSH(MXSHEL)
   end module blocks_mod
