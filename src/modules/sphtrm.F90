!-----------------------------------------
! Module with variables previously defined
! at sphtrm.h include file (COMMON /SPHTRM/)
!-----------------------------------------
!  CSP : Cartesian->Spherical transformation matrices
!  CSP( ISPADR(L+1) ) : first element of transformation matrix for ang.mom. = L
!  -- CSP and ISPADR are set in SPHINP(herrdn.F)
   module sphtrm_mod
   use maxaqn_mod, only: MXQN
   implicit none
   integer, parameter :: NCSP = ((MXQN+2)*(MXQN+1)*(3*MXQN**2+6*MXQN+1)*MXQN)/60
   real*8 :: CSP(NCSP)
   integer :: ISPADR(MXQN)
   end module sphtrm_mod
