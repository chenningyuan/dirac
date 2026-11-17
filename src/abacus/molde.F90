!-----------------------------------------
! Module with variables previously defined
! at molde.h include file (COMMON /MOLDEN/)
!-----------------------------------------
   module molde_mod
   use maxorb_mod, only: MXCORB
   implicit none
   real*8 :: OREN(MXCORB), OROC(MXCORB), EMCEN(MXCORB)
   integer :: IXYZPRINTS
   logical :: MOLDEN, DONEIT, DONEIU, DONEIV, DONEIW, DONEIX
   end module molde_mod
