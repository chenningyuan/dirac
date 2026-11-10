!-----------------------------------------
! Module with variables previously defined
! at inftap.h include file (COMMON /INFTAP/, /CHRTAP/)
!-----------------------------------------
! include/inftap.h (from Dalton, needed for ABACUS routines)
   module inftap_mod
   implicit none
   integer :: LUINP, LUINTA, LUINTM, LUINTR, LUSUPM, LUONEL, LUIT1, LUTEMP, LUDA, LUINDX, &
              LUPROP, LUINFO, LUSOL, LUSIFC, LBINTM
   character(len=8) :: FNONEL, FNSOL, FNSUPM, FNSIFC, LBSIFC
   end module inftap_mod
