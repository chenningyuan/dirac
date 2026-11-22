!-----------------------------------------
! Module with variables previously defined
! at src/luciarel/cands.inc include file (COMMON /CANDSR/)
!-----------------------------------------
!  Definition of c and sigma spaces
   module cands_luciarel_mod
   implicit none
   integer :: ICSM, ISSM, ICSPC, ISSPC, IASTP, IBSTP, IACTP, IBCTP, IELCI
!  IA/BS/CTP : Type (number of elecs in A/B,S/C)
!  IELCI     : One- or many-electron CI (2 means many!)
   end module cands_luciarel_mod
