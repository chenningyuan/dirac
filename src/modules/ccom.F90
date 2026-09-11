!-----------------------------------------
! Module with variables previously defined
! at ccom.h include file
!-----------------------------------------
   module ccom_mod
   use maxaqn_mod, only: MXQN
   use kinds_mod, only: dp
   implicit none
!
   LOGICAL DOCART, SPH(MXQN), SPHNRM, DOSPHE
   REAL(dp) THRS
   INTEGER NHTYP, IUTYP, KHK(MXQN), KCK(MXQN), NHKOFF(MXQN)
   CHARACTER*4 GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6)
!
   end module ccom_mod