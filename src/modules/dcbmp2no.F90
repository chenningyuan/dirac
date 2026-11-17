!-----------------------------------------
! Module with variables previously defined
! at dcbmp2no.h include file (COMMON /CBLMP2NO/, /CBIMP2NO/, /CBCMP2NO/)
!-----------------------------------------
! dcbmp2no.h -- info for MP2 natural orbitals
   module dcbmp2no_mod
   implicit none
   integer :: IPRMP2NO, INTMP2NO, IF2MP2NO, IAOMOSCHEME, LMP2_DEN_INP, MP2NOFN_L, IMO_NAT, &
              IOCC_NAT, LUNATORB, NACTMP2NO(2)
   logical :: RUN_CCMOD, MP2_DENRUN, SKIP_MP2DEN, SKIP_ITRAFO, MP2NATPOP, MP2NATPOP2
   character(len=14) :: MP2NOFN
   character(len=6) :: SELPOPNAT
   character(len=62) :: INFO_STRING(3)
   character(len=72) :: VECMP2NO(2), VECNOSEL(2)
   end module dcbmp2no_mod
