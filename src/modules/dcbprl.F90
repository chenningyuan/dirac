!-----------------------------------------
! Module with variables previously defined
! at dcbprl.h include file (COMMON /XCBPRL/, /XCBLBL/)
!-----------------------------------------
!  PROPERTY LABELS:
!     PRPLBL  - property label
!     IPRLREP - boson irrep of scalar operator (0 - 7)
!     IPRLTYP - symmetry of scalar operator under matrix transpose (-1,0,1)
   module dcbprl_mod
   implicit none
   integer, parameter :: MAXPRPLBL = 3200
   integer :: NPRPLBL, IPRLREP(MAXPRPLBL), IPRLTYP(MAXPRPLBL)
   character(len=4) :: PDOINT(MAXPRPLBL)
   character(len=8) :: PRPLBL(MAXPRPLBL)
   end module dcbprl_mod
