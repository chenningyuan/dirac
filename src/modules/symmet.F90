!-----------------------------------------
! Module with variables previously defined
! at symmet.h include file
!-----------------------------------------
   module symmet_mod
   use maxaqn_mod, only: MXQN, MXAQN
   use mxcent_mod, only: MXCENT
   use primitives_mod, only: MXCORB
   implicit none
!
   integer, parameter :: NSYMMTR = 16
!
   ! ---- Variables from old COMMON /SYMMTI/ ----
   integer :: I1_SYMMTI, MAXREP, MAXOPR
   integer :: MULT(0:7)
   integer :: ISYMAX(3, 2)
   integer :: ISYMAO(MXQN, MXAQN)
   integer :: NPARSU(8)
   integer :: NAOS(8)
   integer :: NPARNU(8, 8)
   integer :: IPTSYM(MXCORB, 0:7)
   integer :: IPTCNT(3*MXCENT, 0:7, 2)
   integer :: NCRREP(0:7, 2)
   integer :: IPTCOR(3*MXCENT, 2)
   integer :: NAXREP(0:7, 2)
   integer :: IPTAX(3, 2)
   integer :: IPTXYZ(3, 0:7, 2)
   integer :: IPTNUC(MXCENT, 0:7)
   integer :: ISOP(0:7)
   integer :: NROTS, NINVC, NREFL
   integer :: IXVAL(0:7, 0:7)
   integer :: NCOS(8, 2)
   integer :: ICLASS(MXCORB)
   integer :: ICNTAO(MXCORB)
   integer :: IPAR(0:7)
   integer :: JSOP(0:7)
   integer :: ICOS(8, 2)
   integer :: I2COSX(8, 8, 2, 2)
   integer :: I2_SYMMTI
!
   ! ---- Variables from old COMMON /SYMMTR/ ----
   real*8 :: FMULT(0:7)
   real*8 :: PT(0:7)
!
   end module symmet_mod
