!-----------------------------------------
! Module with variables previously defined
! at comdis.h include file (COMMON /COMDIS/, /DISBUF/, /DIS2AB/)
!-----------------------------------------
!*** Information for generating two-electron distributions
!    with two indices fixes
!     Length of common block : keep up to date !!
   module comdis_mod
   implicit none
   integer :: NSPCK(0:7,0:3), ISPCK(0:7,0:7,3), NPASS, NGBFSZ, LGFIL, IRECG, KLGREC, &
              KLGBUF, KRGBUF, KIGBUF, ISHLA, ISHLB, ISETA, ISETB, NINSHA, NINSHB
   end module comdis_mod
