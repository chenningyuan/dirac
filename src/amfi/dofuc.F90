!-----------------------------------------
! Module with variables previously defined
! at dofuc.h include file (COMMON /DOFUC/)
!-----------------------------------------
!Comdeck dofuc $Revision$
   module dofuc_mod
   use para_mod, only: NDFMX
   implicit none
   real*8 :: DF(0:ndfmx), DFFRAC(0:ndfmx,0:ndfmx)
   end module dofuc_mod
