!-----------------------------------------
! Module with variables previously defined
! at dcblab.h include file (COMMON /LABNAM/, /LABELS/)
!-----------------------------------------
!*** Labels ***
!  PLABEL(*,1), IATTR(*,1), IPLAB(*,1) for AO basis
!  PLABEL(*,2), IATTR(*,2), IPLAB(*,2) for SO basis
!  NPLAB(0:2): number of labels total/large comp/small comp
   module dcblab_mod
   use maxorb_mod, only: MXCORB
   implicit none
   integer, parameter :: MAXLAB = 12000
   integer :: IATTR(MAXLAB,2), IPLAB(MXCORB,2), NPLAB(0:2)
   character(len=12) :: PLABEL(MAXLAB,2)
   end module dcblab_mod
