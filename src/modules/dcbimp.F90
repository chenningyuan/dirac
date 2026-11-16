!-----------------------------------------
! Module with variables previously defined
! at dcbimp.h include file (COMMON /DCBIMP/)
!-----------------------------------------
!----------------------------------------------------------------------!
!     variables needed for the MO import
!----------------------------------------------------------------------!
   module dcbimp_mod
   implicit none
   logical :: DOIMPMOS
   character(len=64) :: QCPACK, FILIMP(2)
   end module dcbimp_mod
