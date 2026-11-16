!-----------------------------------------
! Module with variables previously defined
! at dcbhoc.h include file (COMMON /DCCHOC/, /DCRHOC/)
!-----------------------------------------
!** Atomic Huckel start ***
   module dcbhoc_mod
   implicit none
   integer, parameter :: MAXHOC = 15
   real*8 :: HUCPAR
   character(len=6) :: HOCFIL(MAXHOC)
   character(len=72) :: VECHOC(MAXHOC)
   end module dcbhoc_mod
