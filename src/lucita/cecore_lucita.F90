!-----------------------------------------
! Module with variables previously defined
! at src/lucita/cecore.inc include file (COMMON /CECORE/)
!-----------------------------------------
!comdeck cecore.inc $Revision$ 
   module cecore_lucita_mod
   implicit none
! ECORE_HEX: only in the inline COMMON/CECORE/ of some lucita units (same block)
   real*8 :: ECORE, ECORE_ORIG, ECORE_H, ECORE_HEX
   end module cecore_lucita_mod
