!-----------------------------------------
! Module with variables previously defined
! at mp1stpr.h include file (COMMON /MP2PRP/, /MP2CPR/)
!-----------------------------------------
!  MP1STPR:  Parameters for calculation of 1st order properties at
!            the MP2 level
   module mp1stpr_mod
   implicit none
   integer, parameter :: LUMPPR = 92
   integer :: NTMATR1, NTMATR2, IOFFTM1(2,0:7,2), IOFFTM2(2,0:7,2)
   character(len=72) :: MPR_INDSTR(4,2)
   end module mp1stpr_mod
