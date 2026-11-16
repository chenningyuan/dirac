!-----------------------------------------
! Module with variables previously defined
! at ired.h include file (COMMON /IREDUCED2H/)
!-----------------------------------------
!Comdeck ired $Revision$
!bs   irreducible representation of the cartesian functions in D2H
!bs   order taken from Tables for Group theory:
!bs   Atkins, Child and Phillips   Oxford University Press 1970
!bs   1. AG: only even powers (0,0,0)
!bs   2. B1G: (1,1,0)    L_z 
!bs   3. B2G: (1,0,1)    L_y
!bs   4. B3G: (0,1,1)    L_x
!bs   5. AU:  (1,1,1)
!bs   6. B1U: (0,0,1)
!bs   7. B2U: (0,1,0)
!bs   8. B3U: (1,0,0)
   module ired_mod
   use para_mod, only: LMAX, MXCART
   implicit none
   integer :: NUMBOFSYM, IPOW2IRED(0:1,0:1,0:1), IREDORDER(8), IREDORDERINV(8), &
              NFUNCTIONS(8,0:Lmax), NFUNCTPERIRED(8), INCRLM(-Lmax:Lmax,0:Lmax), &
              SHIFTIRED(8), IREDLM(-Lmax:Lmax,0:Lmax), SHIFTIRIR(36), LOFFUNCTION(Mxcart), &
              MOFFUNCTION(Mxcart), IREDOFFUNCTION(Mxcart), IREDOFFUNCTNEW(Mxcart), &
              IRWITHLX(8), IRWITHLY(8), IRWITHLZ(8), MULT(8,8), ITOTALPERIR(8), &
              NMBMPERIRL(8,0:Lmax), NUMBOFCART
   end module ired_mod
