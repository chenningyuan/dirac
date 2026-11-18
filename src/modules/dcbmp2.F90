!-----------------------------------------
! Module with variables previously defined
! at dcbmp2.h include file (COMMON /DCLMP2/, /DCCMP2/, /DCBMP2/, /DCRMP2/)
!-----------------------------------------
!  DCBMP2:  Parameters for MP2 calculation
!     MP2_NSTR(OCC/VIR, Fermion irrep)   - number of occ/virt pairs
!     MP2_INDSTR(OCC/VIR, Fermion irrep) - string of active pairs
   module dcbmp2_mod
   implicit none
   real*8 :: DMP2_VIRTHR, SCRMP2, EMP2, TPRI34, TPRI44
   integer :: IPRMP2, ISTRMP2, MP2_INTFLG, MP2_NSTR(2,2)
   logical :: MP2_MSOUT, MP2_ANTIS, MP2NOP, SRTSHL, CMMNFO
   character(len=72) :: MP2_INDSTR(2,2)
   end module dcbmp2_mod
