!-----------------------------------------
! Module with variables previously defined
! at infvar.h include file (COMMON /INFVAR/)
!-----------------------------------------
   module infvar_mod
   use maxorb_mod, only: MXCORB
   implicit none
#if defined (VAR_SIRBIG)
   integer, parameter :: MAXWOP = 30000, MAXOCC = 120
#else
   integer, parameter :: MAXWOP = 5000, MAXOCC = 120
#endif
   integer :: NCONF, NWOPT, NVAR, JWOPSY, NWOP(8), NWOPH, NVARH
#if defined (VAR_STAR2)
   integer*2 :: JWOP(2,MAXWOP), KLWOP(MAXOCC,MXCORB)
#else
   integer :: JWOP(2,MAXWOP), KLWOP(MAXOCC,MXCORB)
#endif
   end module infvar_mod
