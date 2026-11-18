!-----------------------------------------
! Module with variables previously defined
! at cbivib.h include file (COMMON /CBIVIB/, /CBIVBC/)
!-----------------------------------------
   module cbivib_mod
   use mxcent_mod, only: MXCENT, MXCOOR
   implicit none
   integer, parameter :: MAXSUB = 20
   real*8 :: SCALE(MXCOOR), COEF(MXCOOR)
   integer :: IATOMS(4,MXCOOR), NCARD, NINTCM, IPRINT, MAXDIF, NISOTP, ISOTP(MAXSUB,MXCENT)
   logical :: DIPOL, SKIP, DOVCD, HESFIL, HESPUN, ROACID
   character(len=1) :: KWORD(MXCOOR)
   character(len=4) :: ITYPCD(MXCOOR), ITYPCM(MXCOOR)
   end module cbivib_mod
