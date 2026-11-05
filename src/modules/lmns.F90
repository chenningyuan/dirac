!-----------------------------------------
! Module with variables previously defined
! at lmns.h include file (COMMON /LMNS/)
!-----------------------------------------
   module lmns_mod
   use maxaqn_mod, only: MXAQN
   implicit none
   integer :: LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN), LVALUB(MXAQN), MVALUB(MXAQN), &
              NVALUB(MXAQN)
   end module lmns_mod
