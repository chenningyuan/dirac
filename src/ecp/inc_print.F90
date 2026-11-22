!-----------------------------------------
! Module with variables previously defined
! at src/ecp/inc_print.h include file (COMMON /RECP_VAR_PRINT/)
!-----------------------------------------
!  RECP_FNUM(10) and RECP_FNAM(10) (I/O file numbers and names) were dropped:
!  RECP_SETPARM set them with DATA and nothing ever read them.
   module inc_print_mod
   implicit none
   integer :: PRINT_LEVEL
   integer :: RECP_OUT      ! RECP output
   integer :: RECP_DBG      ! For debug program
   end module inc_print_mod
