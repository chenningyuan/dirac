!-----------------------------------------
! Module with variables previously defined
! at orgcom.h include file (COMMON /ORGCOM/)
!-----------------------------------------
!   CMXYZ  : center-of-mass
!   NCHBC  : nuclei charge barycenter
!   ORIGIN : general one-electron operator origin
!   DIPORG : origin for multipole operators, including dipole
!   GAGORG : magnetic gauge origin
!     GAGORG_SET true: GAGORG has been set by user
!   CAVORG : cavity origin for spherical cavity model
   module orgcom_mod
   implicit none
   real*8 :: CMXYZ(3), ORIGIN(3), DIPORG(3), GAGORG(3), CAVORG(3), NCHBC(3)
   logical :: GAGORG_SET
   end module orgcom_mod
