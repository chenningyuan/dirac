!-----------------------------------------
! Module with variables previously defined
! at amfi_if.h include file (COMMON /AMFI_RELSCF/, /AMFI_RELSCF_LOG/, /AMFI_RELSCF_CHR/)
!-----------------------------------------
!----------------------------------------------------------------------------
!  Common file for passing information from DIRAC to AMFI/RELSCF
! ---------------------------------------------------------------
!   IPR_AMFI   - print level variable for AMFI part of B.Schimm.
!   IPR_RELSCF - print level for relscf program interfaced to AMFI
!   MXXTRP_INP - maximum number of iterations for relscf
!   ONRELSCF   - stop after relscf run
!   MXRELSCF_INFO - dimension of info array IRELSCF_INFO
!                   dimension is at present LMAX (see para.h) * 2 + 7
!                   but only LMAX*2 + 2 is used; extra space reserved
!                   when LMAX will be extended to at least LMAX = 6
!                   (h and i functions)
!   IRELSCF_INFO    info array for RELSCF run
!   RELSCFLABx      label on file RELSCF_COEF for integer (x = I)
!                   and for REAL*8 (x = X)
! ... other parameters can be added here as well ...
!----------------------------------------------------------------------------
   module amfi_if_mod
   implicit none
   integer, parameter :: MXRELSCF_INFO = 15
   integer :: IPR_AMFI, IPR_RELSCF, MXXTRP_INP, IRELSCF_INFO(MXRELSCF_INFO)
   logical :: ONRELSCF
   character(len=8) :: RELSCFLABI, RELSCFLABX
   end module amfi_if_mod
