!-----------------------------------------
! Module with variables previously defined
! at cbilnr.h include file (COMMON /LNLBL/, /CBILNR/)
!-----------------------------------------
!  Used to control linear response solver in ABACUS
   module cbilnr_mod
   implicit none
   integer, parameter :: MXFR = 101
   integer, parameter :: MAXLN = 80
   real*8 :: THCLNR, FRVAL(MXFR)
   integer :: NFRVAL, IPRLNR, NABALN
   logical :: ALFA, ROAA, ROAG, STATIC
   character(len=8) :: LABALN(MAXLN)
   end module cbilnr_mod
