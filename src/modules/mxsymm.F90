!-----------------------------------------
! Module with variables previously defined
! at mxsymm.h include file (COMMON /MXSYMM/)
!-----------------------------------------
!  Symmetry detection (values set and used in hersym.F)
   module mxsymm_mod
   implicit none
   real*8 :: TOLLRN, TOLLRN_2, ZERTOL
   integer :: MAXAXS, MAXMIR, IPRSYM
   end module mxsymm_mod
