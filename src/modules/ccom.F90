!-----------------------------------------
! Module with variables previously defined
! at ccom.h include file (COMMON /CCOM/, /CCOMC/)
!-----------------------------------------
   module ccom_mod
   use maxaqn_mod, only: MXQN
   implicit none
   real*8 :: THRS
   integer :: NHTYP, IUTYP, KHK(MXQN), KCK(MXQN), NHKOFF(MXQN)
   logical :: DOCART, SPHNRM, SPH(MXQN), DOSPHE
   character(len=4) :: GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6)
   end module ccom_mod
